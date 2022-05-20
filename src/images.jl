using Plots
using DataFrames
using Images
using Statistics
using ImageBinarization


mutable struct Setup
    cmdir::String
    series::String
    measurement::String
    point::String
    filter::String
    rep::String
    imgid::String 
    texpsec::Float64 
    darktype::String 
    sfilter::String 
    srep::String 
    img::String 
    dark::String 
    exp::String 
    fext::String 
    pext::String 
end

function select_image_name(setup::Setup, imgid="Dark")

    if imgid == "Dark"
        if setup.darktype == "before"
            string(setup.sfilter, "_" , setup.filter, "_" , setup.exp, "_", setup.darktype, 
                                  "_", setup.dark,  setup.fext)
        else
            string(setup.sfilter, "_" , setup.filter, "_" , setup.exp, "_", setup.dark,  setup.fext)
        end
    else
        string(setup.sfilter, "_", setup.filter, "_", setup.srep, "_", setup.rep, "_", setup.exp, 
               "_", setup.img,  setup.fext)
    end

end


function get_image_name(inames::Vector{String}, imgn::String)
	names = [split(f,"/")[end] for f in inames]
	indx = findall([occursin(imgn, name) for name in names])[1]
	inames[indx]
end


function get_image(ffimg::String)
	imgdf = DataFrame(CSV.File(ffimg, header=false,delim="\t"))
	img   = dftomat(imgdf)
    imgn = img ./maximum(img)
    (img = img, imgn = imgn)
end



function select_image(xfiles::Vector{String}, setup::Setup, imgid::String) 
	imgn = select_image_name(setup, imgid)
    println("imgn = ", imgn)
	ximg = get_image_name(xfiles, imgn)
    #println("ximg = ", ximg)
	get_image(ximg)
end


function get_corrected_image(files::Vector{String}, setup::Setup) 
	imgm = select_image(files, setup, "Imag")
	drkm = select_image(files, setup, "Dark")
	
	img = imgm.img .- drkm.img
	imgn = img ./maximum(img)
    (img = img, imgn = imgn, dark=drkm.img)
end


function spectrum_max(xfiles::Vector{String}, setup::Setup,  
                      xfn::Vector{String}, filtnm::NamedTuple, adctopes::Float64; nsigma::Float64=3.0)
	ZMX = Vector{Float64}()
	ZSM = Vector{Float64}()
	ZI = Vector{Float64}()
	ZJ = Vector{Float64}()
	
    fold = setup.filter
	for fltr in xfn
        setup.filter = fltr 
        #println("in spectrum_max: setup point = ", setup.point)
		cimgz = get_corrected_image(xfiles, setup)
		imgmz, imgp = signal_around_maximum(cimgz.img, cimgz.dark; nsigma=nsigma)
		push!(ZMX,imgp.max)
		push!(ZSM,sum(imgmz.img))
		push!(ZI,imgp.imax)
		push!(ZJ,sum(imgp.jmax))

        #println(" for filter ",fltr, " sum = ", sum(imgmz.img))
	end
    setup.filter = fold
	DataFrame("fltn" => xfn, "cflt" => filtnm.center, "lflt" => filtnm.left, "rflt" => filtnm.right, "wflt" => filtnm.width,
		      "sum"=> ZSM, "sumpes"=> adctopes *(ZSM ./filtnm.width), "max"=> ZMX, "imax" => ZI, "jmax" => ZJ)	
end


function spectrum_max_allpoints!(setup::Setup, 
                                xpt::Vector{String}, xfn::Vector{String}, 
                                filtnm::NamedTuple, adctopes::Float64; nsigma::Float64=3.0, odir, etype="csv")
    pold = setup.point 
    for pt in xpt
        setup.point = pt
        xpath = joinpath(setup.cmdir,setup.series,setup.measurement, setup.point)
        dfiles = string("*",setup.fext)
	    xfiles = Glob.glob(dfiles, xpath)
        #println("in spectrum_max_allpoints: xpath = ", xpath)
        sdf = spectrum_max(xfiles, setup, xfn, filtnm, adctopes; nsigma=nsigma)
        sdfnm = get_outpath(setup, etype)
        sdff = joinpath(odir, sdfnm)
	    println("Writing point to  =", sdff)
	    CSV.write(sdff, sdf)
    end
    setup.point = pold 
end


function spectrum_sum(xfiles::Vector{String}, setup::Setup, 
                      xfn::Vector{String}, filtnm::NamedTuple, adctopes::Float64)
    
    ZSM = Vector{Float64}()

    #println("in spectrum_sum: setup = ", setup)
    for fltr in xfn
        setup.filter = fltr
        cimgz = get_corrected_image(xfiles, setup)
        push!(ZSM,sum(cimgz.img))
    end

    #println("in spectrum_sum (after): setup = ", setup)
    DataFrame("fltn" => xfn, "cflt" => filtnm.center, "lflt" => filtnm.left, "rflt" => filtnm.right, "wflt" => filtnm.width,
    "sum"=> ZSM, "sumpes"=> adctopes *(ZSM ./filtnm.width))	
end



function get_outpath(setup::Setup, etype="csv")
    
    if etype == "csv"
	    sdfnm = string(setup.series, "_", setup.measurement, "_", 
                   setup.point,  "_", setup.srep,  "_", setup.rep, setup.fext)
    else
        sdfnm = string(setup.series, "_", setup.measurement, "_", 
                   setup.point,  "_", setup.srep,        "_", setup.rep, setup.pext)
    end
end


# function write_spectrum!(sdf::DataFrame, setup::Setup; odir, etype="csv")
    
#     if etype == "csv"
# 	    sdfnm = string(setup.series, "_", setup.measurement, "_", 
#                    setup.point,  "_", setup.srep,        "_", setup.rep, setup.fext)
#     else
#         sdfnm = string(setup.series, "_", setup.measurement, "_", 
#                    setup.point,  "_", setup.srep,        "_", setup.rep, setup.pext)
#     end
                  
# 	sdff = joinpath(odir, sdfnm)
# 	#println("Writing point to  =", sdff)
# 	CSV.write(sdff, sdf)

# end

function read_spectrum(setup::Setup, csvdir)
    
    
	sdfnm = string(setup.series, "_", setup.measurement, "_", 
                   setup.point,  "_", setup.srep,        "_", setup.rep, setup.fext)
    
                  
	sdff = joinpath(csvdir, sdfnm)
	#println("reading file =", sdff)

    load_df_from_csv(csvdir, sdfnm, enG)
	
end


function spectrum_fromfile_allpoints(setup::Setup, xpt::Vector{String}, csvdir)
    dfdict = Dict()
    for pt in xpt
        setup.point = pt
        df = read_spectrum(setup, csvdir)
        dfdict[pt] = df
    end

    dfdict
end

"""
    dftomat(img::DataFrame)

Copy a dataframe into a 2x2 matrix.    

# Fields
- `img::DataFrame`: The "image dataframe", a sqaure dataframe (typically 512 x 512)

"""
function dftomat(img::DataFrame)
	mimg = Array{Float64}(undef, size(img)...);
	for i in 1:size(img)[1]
		for j in 1:size(img)[2]
			if img[i,j] >= 0.0 && !isnan(img[i,j]) 
    			mimg[i, j] = img[i,j]
			else
				mimg[i, j] = 0.0
			end
		end
	end
	mimg
end
"""
    signal_around_maximum(img::Matrix{Float64}; xrange::Int64=10)

Finds the maximum of the image, then adds signal in a window defined by xwindow    

# Fields
- `img::DataFrame`: The "corrected image matrix"
- `dark::DataFrame`: The "corrected dark matrix"
- `nsigma::Int64`: number of sigmas for noise cutoff

"""
function signal_around_maximum(img::Matrix{Float64}, dimg::Matrix{Float64}; nsigma::Float64=2.0)

    cutoff::Float64 =  nsigma * std(dimg)

    #println("mean dimg =", mean(dimg), " std  = ", std(dimg), " cutoff =", cutoff)

	mxx::Float64, indxt = findmax(img)
    sx::Int = size(img)[1]
    sy::Int = size(img)[2]
    kx::Int = indxt[1]
    ky::Int = indxt[2]

    #println("max =", mxx)
    #println("sx =", sx, " sy = ", sy, " kx = ", kx, " ky = ", ky)

    if mxx < cutoff
        #println("max =", mxx, " less than cutoff! ", cutoff, " taking max")
        img2 = Array{Float64}(undef,1, 1)
        img2[1,1] = mxx 
        imgn = 1.0
        return (img = img2, imgn = imgn), (max=mxx, imax=kx, jmax=ky)
    end

    ik::Int = 0
    for ii in kx:sx
        if img[ii,ky] > cutoff
            ik = ii
        else
            break
        end
    end

    jk::Int = 0
    for jj in ky:sy
        if img[kx,jj] > cutoff
            jk = jj
        else
            break
        end
    end

    ixr::Int = min(sx, ik)
    iyr::Int = min(sy, jk)

    #println("ik = ", ik, " jk = ", jk)
    #println("ixr = ", ixr, " iyr = ", iyr)

	for ii in kx:-1:1
        if img[ii,ky] > cutoff
            ik = ii
        else
            break
        end
    end

    for jj in ky:-1:1
        if img[kx,jj] > cutoff
            jk = jj
        else
            break
        end
    end

	ixl = max(1, ik)
	iyl = max(1, jk)

    #println("ik = ", ik, " jk = ", jk)
    #println("ixl = ", ixl, " iyl = ", iyl)

	img2 = Array{Float64}(undef, ixr - ixl + 1, iyr - iyl + 1)
	
	for (i, ix) in enumerate(ixl:ixr)
		for (j, iy) in enumerate(iyl:iyr)
			img2[i,j]=img[ix,iy]
		end
	end

    #println("size of img2 =", size(img2))

    imgn = img2 ./maximum(img2)
    (img = img2, imgn = imgn), (max=mxx, imax=kx, jmax=ky)

end






"""
    sum_edge(img::Matrix{Float64}, edge::Matrix{Float64})

Adds all the pixels in matrix img that are identified as "edge" (>0) in matrix edge    

# Fields
- `img::Matrix{Float64}`: Image matrix
- `edge::Matrix{Float64}`: edge matrix

"""
function sum_edge(img::Matrix{Float64}, edge::Matrix{Float64})
	sum::Float64 = 0.0
	nedge::Int   = 0
	for i::Int in 1:size(img)[1]
		for j::Int in 1:size(img)[2]
			if edge[i,j] > 0.0 
				sum += img[i,j]
				nedge += 1
			end
		end
	end
	sum, nedge
end

"""
    edges = sujoy(img; four_connectivity=true)

Compute edges of an image using the Sujoy algorithm.

# Parameters

* `img` (Required): any gray image
* `four_connectivity=true`: if true, kernel is based on 4-neighborhood, else, kernel is based on
   8-neighborhood,

# Returns

* `edges` : gray image
"""
function sujoy(img; four_connectivity=true)
    img_channel = Gray.(img)

    min_val = minimum(img_channel)
    img_channel = img_channel .- min_val
    max_val = maximum(img_channel)

    if max_val == 0
        return img
    end

    img_channel = img_channel./max_val

    if four_connectivity
        krnl_h = centered(Gray{Float32}[0 -1 -1 -1 0; 0 -1 -1 -1 0; 0 0 0 0 0; 0 1 1 1 0; 0 1 1 1 0]./12)
        krnl_v = centered(Gray{Float32}[0 0 0 0 0; -1 -1 0 1 1;-1 -1 0 1 1;-1 -1 0 1 1;0 0 0 0 0 ]./12)
    else
        krnl_h = centered(Gray{Float32}[0 0 -1 0 0; 0 -1 -1 -1 0; 0 0 0 0 0; 0 1 1 1 0; 0 0 1 0 0]./8)
        krnl_v = centered(Gray{Float32}[0 0 0 0 0;  0 -1 0 1 0; -1 -1 0 1 1;0 -1 0 1 0; 0 0 0 0 0 ]./8)
    end

    grad_h = imfilter(img_channel, krnl_h')
    grad_v = imfilter(img_channel, krnl_v')

    grad = (grad_h.^2) .+ (grad_v.^2)

    return grad
end

"""
    plot_images(imgnt)

Plots the images with laser on and off, as well as corrected, as heatmaps. 

# Fields
- `imgnt::NamedTuple`: A named tuple containing the three images
- `whichf::String`: filter number (by convention from 2 to 11)
- `whichr::String`: repetition number (1,2,3): files with repetitions of the same measurement.

"""
function plot_images(imgnt::NamedTuple, whichr::String, whichf::String)
	imghm = heatmap(imgnt.image, title=string("rep=", whichr, " filt=", whichf, " Img" ), titlefontsize=10)

	drkhm = heatmap(imgnt.dark, title=string("rep=", whichr, " filt=", whichf, " Dark" ), titlefontsize=10)
	
	simghm = heatmap(imgnt.cimage, title=string("rep=", whichr, " filt=", whichf, " Img - Dark" ))

    sdarkhm = heatmap(imgnt.cdark, title=string("rep=", whichr, " filt=", whichf, " Dark - Dark" ))
	
	plot(imghm, drkhm, simghm, sdarkhm, layout=(2,2), titlefontsize=10)
    xlims!(0,512)
    yticks!([0,250,512])
    xticks!([0,250,512])
    ylims!(0,512)
end