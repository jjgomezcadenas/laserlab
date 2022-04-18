using Plots
using DataFrames
using Images
using Statistics
using ImageBinarization

"""
    select_image(inames::Vector{String}, fnumber::String, repetition::String,
                 itype::String="Imag",  expo::String="_ExpoT_1000000ms_")

Selects an image to be read from a vector holding the names of files with data. 

# Fields
- `inames::Vector{String}`: Names of files (tipically obtained with glob)
- `fnumber::String`: filter number (by convention from 2 to 11)
- `repetition::String`: repetition number (1,2,3): files with repetitions of the same measurement.
- `itype::String`: type of image, either Imag (with laser on) or Dark (with laser off)
- `expo::String`: A string describing the exposure 
"""
function select_image(inames::Vector{String}, fnumber::String, repetition::String,
                      itype::String="Imag", expo::String="_ExpoT_1000000ms_")

    names = [split(f,"/")[end] for f in inames]

    if repetition == "0"
        fhead = string("Filter_", fnumber, expo)
    else
        fhead = string("Filter_", fnumber, "_rep_", repetition, expo)
    end
    
    
    if itype == "Dark"
        ftail = string("Dark_1.dat")
    elseif itype == "Dark2" 
        ftail = string("Dark_2.dat")
    elseif itype == "Imag" 
        ftail = string("Imag_1.dat")
    else
        @error "itype = $itype not implemented"
    end

    ff = string(fhead, ftail)
    indx = findall([occursin(ff, name) for name in names])[1]
    inames[indx]
end

"""
    get_images(files::Vector{String}, whichf::String, whichr::String)

For a filter number and a repetition number returns 3 images: laser on, laser off, and subctracted.  

# Fields
- `files::Vector{String}}`: Names of files (tipically obtained with glob)
- `whichf::String`: filter number (by convention from 2 to 11)
- `whichr::String`: repetition number (1,2,3): files with repetitions of the same measurement.

"""
function get_images(files::Vector{String}, whichf::String, whichr::String)
	ffimg = select_image(files, whichf, whichr, "Imag");
	imgdf = DataFrame(CSV.File(ffimg, header=false,delim="\t"));
	imgm  = dftomat(imgdf)
	
	ffdrk = select_image(files, whichf, whichr, "Dark");
	drkdf = DataFrame(CSV.File(ffdrk, header=false,delim="\t"));
	drkm  = dftomat(drkdf)

    ffdrk2 = select_image(files, whichf, whichr, "Dark2");
	drkdf2 = DataFrame(CSV.File(ffdrk2, header=false,delim="\t"));
	drkm2  = dftomat(drkdf2)

	simgdf = imgdf .- drkdf
	simgm  = dftomat(simgdf)

    sdarkdf = drkdf .- drkdf2
	sdark  = dftomat(sdarkdf)
	
	(image=imgm, dark=drkm, cimage=simgm, cdark=sdark)
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
function signal_around_maximum(img::Matrix{Float64}, dark::Matrix{Float64}; nsigma::Int64=2)

    cutoff::Float64 = mean(dark) + nsigma * std(dark)

    println("mean dark =", mean(dark), " std dark = ", std(dark), " cutoff =", cutoff)

	mxx::Float64, indxt = findmax(img)
    sx::Int = size(img)[1]
    sy::Int = size(img)[2]
    kx::Int = indxt[1]
    ky::Int = indxt[2]

    println("max =", mxx)
    println("sx =", sx, " sy = ", sy, " kx = ", kx, " ky = ", ky)

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

    println("ik = ", ik, " jk = ", jk)
    println("ixr = ", ixr, " iyr = ", iyr)

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

    println("ik = ", ik, " jk = ", jk)
    println("ixl = ", ixl, " iyl = ", iyl)

	img2 = Array{Float64}(undef, ixr - ixl + 1, iyr - iyl + 1)
	
	for (i, ix) in enumerate(ixl:ixr)
		for (j, iy) in enumerate(iyl:iyr)
			img2[i,j]=img[ix,iy]
		end
	end

    println("size of img2 =", size(img2))

	(max=mxx, imax=kx, jmax=ky, img=img2) 
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