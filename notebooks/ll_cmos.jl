### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 981730a6-61fc-484b-ba3c-66920ee7cf83
using Pkg; Pkg.activate("/Users/jj/JuliaProjects/LaserLab/")

# ╔═╡ 8e7ec382-c738-11ec-3aae-b50d60f15c4f
begin
	using PlutoUI
	using CSV
	using DataFrames
	using Images
	using ImageBinarization
	using Colors
	using Plots
	using Printf
	using Interpolations
	using QuadGK
	using Markdown
	using InteractiveUtils
	using LsqFit
	using Statistics
	using StatsBase
	using Unitful 
	using UnitfulEquivalences 
	using PhysicalConstants
	using Peaks
	using FFTW
	using DSP
	import Glob

end

# ╔═╡ 06b8ed45-43bc-464f-89c0-dc0406312b81
import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
    A, N, mol, mmol, V, L, M

# ╔═╡ 8833b198-03e4-4679-8949-0c76546cb847
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ 6163ba69-1237-4b49-988e-9a73cfef67f6
lfi = ingredients("../src/LaserLab.jl")

# ╔═╡ 5edc41bc-b912-44bf-9be5-a013f27a75ab
PlutoUI.TableOfContents(title="Laser Lab CMOS analysis", indent=true)

# ╔═╡ c9aaf1cc-80c4-475b-8a81-e00918d91b1e
md"""
# Analysis for a single point
"""

# ╔═╡ 7ce42aec-b319-4de9-b70c-84046d45a600
md"""
## Filters
"""

# ╔═╡ 58269465-ba8c-4840-bbfc-0a27897f3e2a
md"""
### Filter central value in nm
"""

# ╔═╡ 0b1c5662-ec4f-486e-9ee6-7fa6ba953e45
begin
	xfnm = [420.0,438.0,465.0,503.0,550.0,600.0,650.0,692.0,732.0,810.0]
	wfnm = [10.0, 24.0, 30.0, 40.0, 49.0, 52.0, 60.0, 40.0, 68.0, 10.0]
	filtnm = (center=xfnm,
	width  = wfnm,
	left = xfnm .- 0.5*wfnm,
	right = xfnm .+ 0.5*wfnm)

	println("Filter central values (nm) = ", filtnm.center, " width (nm) =", filtnm.width)
end

# ╔═╡ b98ad447-6055-46e5-bb4f-9e67f9c3176a
md"""
## Conversion factors
"""

# ╔═╡ 88853edb-dc1f-4e7a-a0ba-1868276d1ada
begin
	texp = 10.0
	adctopes = 0.48
	md"""
- Time of exposition (in sec) = $texp
- adc to photoelect = $adctopes
"""
end

# ╔═╡ 88c44435-9685-41c3-a6c8-ffaf7498d60f
md""" Select nsigma cut: $(@bind spsigma NumberField(1.0:10.0, default=5.0))"""

# ╔═╡ f6dcfbd5-8416-42f5-b029-c794f92ee413
nsigma=promsel = Float64(spsigma)

# ╔═╡ dde07319-e6a0-48ac-94d6-92e0fda4df41
md"""
## Cutoff around maximum (number of sigmas in terms or rms of the dark current)
- nsigma = $nsigma
"""

# ╔═╡ 57d96432-4318-4291-8255-bfa5d6d3635c
md"""
## Define Experiment, run and point
"""

# ╔═╡ 0b4d2c08-a677-492c-b5da-982d3d5096fc
begin
	#cmdir   = "/Users/jj/JuliaProjects/LaserLab/labdata/CMOS"
	odir    = "/Users/jj/JuliaProjects/LaserLab/data/CMOS"
	cmdir   = "/Users/jj/LaserLab/Proyectos/data/CMOS2"
	#odir   = "/Users/jj/LaserLab/Proyectos/pdata/CMOS2"
	
	dfiles  = "*.csv"
	dplots  = "*.png"
	dmtype  = ["Imag","Dark"]
	rep     = "1"  # always take repetition number 1 if there is more than 1
	drkpnt  = "Dark"
	imgp    = "ExpoT_10s_Imag_1"
	darkp   = "ExpoT_10s_dark"
	md"""
	CMOS dir = $cmdir
	"""
end

# ╔═╡ 3ae7315b-2056-4e3f-9d41-b6e0292f30e3
#sexp ="IrSL_BOLD_068_A1"

# ╔═╡ e87f48e3-5e5a-44d5-83de-c520e522e33a
let
	readdir(cmdir)
	dirs = lfi.LaserLab.getbolddirs(cmdir)
	md""" Select experiment : $(@bind sexp Select(dirs))"""
end

# ╔═╡ 7e185082-f3bf-4d95-9d5d-57101f47a684
md"""
## Select image

- Image in the left not corrected 
- Image in the right corrected (dark current subtracted)
"""

# ╔═╡ e9e8e58a-e1db-49f1-8429-420271fb1852
md"""
## Comparison between spectrum computed around maximum and using the full camera. 

- Notice the umphysical shoulder when using the full camera
"""

# ╔═╡ 5477b1dc-ef52-4426-aa47-f513854fcdae
md""" Check to compute sum using full cmos (dark current): $(@bind zfull CheckBox())"""

# ╔═╡ a48af8f4-4ed2-45cd-b4e8-9b3106c885f3
md"""
# Analysis for al points
"""

# ╔═╡ 1794afb6-6ef0-46d6-b182-d54362b9a07d
md""" Check to carry analysis for all points: $(@bind zrec CheckBox())"""

# ╔═╡ 82ddd81b-8aea-4711-97d9-a645121786f8
md""" Check to read and plot data for all points: $(@bind zread CheckBox())"""

# ╔═╡ 155d5066-935b-4643-8aad-f4c635aa7eec
if zread
	md""" Select nnumber of plots in x and y: 
	- select nx $(@bind nx NumberField(1:10, default=3))
	- select ny $(@bind ny NumberField(1:10, default=3))
	"""
end

# ╔═╡ c1b2cf36-d72e-4348-911e-8e44f1834ae4
md"""
# Functions
"""

# ╔═╡ f2714261-7bb0-47d7-8aac-c16bb5d1f891
struct Setup
    cmdir::String
    series::String
    measurement::String
    point::String
    filter::String
    rep::String 
end


# ╔═╡ 78de6bcd-4173-40e3-b500-499568289ba1
function fixdark(nxdrk)
	splt = split(nxdrk[1],"_")
	str = string(splt[3],"_")
	for i in 4:length(splt) -1
		str = string(str, splt[i],"_")
	end
	string(str, splt[end])
end

# ╔═╡ 0d5a7021-6072-464e-836e-f05b1e178b80
function create_dir!(dir)
	if isdir(dir) == false
		mkdir(dir)
	end
end

# ╔═╡ a759ecf7-7373-46bf-ab15-df39cbfa6814
function output_dirs!(sexp)
	namex = split(sexp,"_")
	create_dir!(joinpath(odir, sexp))
	csvdir = joinpath(odir, sexp,"csv")
	pngdir = joinpath(odir, sexp, "png")
	create_dir!(csvdir)
	create_dir!(pngdir)
	csvdir, pngdir
end

# ╔═╡ 161b12b6-88a0-4d4d-bfc5-01310534cbdc
begin 
	csvdir, pngdir = output_dirs!(sexp)
	md"""
	### Output dirs
	- csv dir = $csvdir
	- png dir = $pngdir
	"""
end

# ╔═╡ 46b1d54e-4ebf-45a2-b137-9690b7a51d38
function select_files(cmdir,sexp,srun, spoint, dfiles="*.csv")
	path = joinpath(cmdir,sexp,srun, spoint)
	readdir(path)
	xfiles = Glob.glob(dfiles, path)
	nxfiles = string.([split(f,"/")[end] for f in xfiles])
	xfiles, nxfiles
end


# ╔═╡ 2498aa42-2c24-47e3-bf5b-647377af0dbc
function select_run(cmdir, sexp)
	namex = split(sexp,"_")
	path = joinpath(cmdir,sexp)
	readdir(path)
	lfi.LaserLab.getbolddirs(path)
end

# ╔═╡ 50ea2ecc-970f-4630-8c7e-acf5e69cc4c9
let
	dirs = select_run(cmdir, sexp)
	md""" Select run : $(@bind srun Select(dirs))"""
end

# ╔═╡ afeb42ca-b462-4320-b364-98a0b4730e33
function select_point(cmdir, sexp, srun)
	path = joinpath(cmdir,sexp,srun)
	readdir(path)
	pdirs = lfi.LaserLab.getbolddirs(path)
	points = [split(pd, "Point")[2] for pd in pdirs if pd != "Dark"]
	spdirs = sort(parse.(Int64, points))
	[string("Point", string(i)) for i in spdirs]
end

# ╔═╡ d389f99a-14c2-408f-ad7b-838e00225357
begin
	sspdirs = select_point(cmdir, sexp, srun)
	md""" Select point : $(@bind spoint Select(sspdirs))"""
end

# ╔═╡ 479f8f86-372c-4b91-9f73-e57a85d3d194
begin
	xfiles, nxfiles  = select_files(cmdir,sexp,srun, spoint)
	xfdrk, nxdrk     = select_files(cmdir,sexp,srun, drkpnt)
end

# ╔═╡ f13173e1-088c-4f5c-ae6a-f9aebcd6bc57
begin
	xfa = lfi.LaserLab.findpattern(nxfiles, "Filter")
	xfb = sort(parse.(Int64, xfa))
	xfn = string.(xfb)
	md""" Select filter: $(@bind sfn Select(xfn))"""
end

# ╔═╡ 8dbf64ec-5854-44b1-ac73-7cd0363a1c6d
begin 
	
	md"""
	### Setup
	- run = $srun
	- experiment = $sexp
	- point = $spoint
	- filter = $sfn
	"""
end

# ╔═╡ 54bd1f6c-2b10-47a1-838f-b428fe6b7635
xfn

# ╔═╡ b857231c-f06c-426d-8f2c-5f8007134714
fixdark(nxdrk)

# ╔═╡ 43ccd7f3-140f-4cb8-af41-e13534c454f3
xfiles

# ╔═╡ d762313e-1466-4724-9c35-dd89e657a11c
begin
	setup = Setup(cmdir, string(sexp), string(srun), string(spoint),        
		              string(sfn), string(rep))
	
	#pngpath = joinpath(pngdir,lfi.LaserLab.get_outpath(setup, "png"))
	md"""
	Setup
	- root dir = $(setup.cmdir)
	- Series = $(setup.series)
	- Measurement = $(setup.measurement)
	- Point = $(setup.point)
	- Filter = $(setup.filter)
	- Repetition = $(setup.rep)
	

	"""
end

# ╔═╡ a90c8edc-1405-4a8c-aacf-53cd130910ae
setup

# ╔═╡ f9ea012f-9aae-4a8b-89d9-0f86802bb14f
function select_image(xfiles::Vector{String}, xfdrk::Vector{String}, 
	                  nxdrk::Vector{String}, setup::Setup) 
	
	
	imgn =string("Filter_",setup.filter,"_rep_", setup.rep, "_ExpoT_10s_Imag_1.csv")
	idrk =string("Filter_",setup.filter, "_", fixdark(nxdrk) )
    #println("imgn = ", imgn)
	#println("idrk = ", idrk)
	ximg = lfi.LaserLab.get_image_name(xfiles, imgn)
	xdrk = lfi.LaserLab.get_image_name(xfdrk, idrk)
    #println("ximg = ", ximg)
	#println("xdrk = ", xdrk)
	lfi.LaserLab.get_image(ximg), lfi.LaserLab.get_image(xdrk)
end

# ╔═╡ 58f3c3ac-698b-4381-bdcb-ca8a9cadc8d8
function get_corrected_image(imgm::NamedTuple{(:img, :imgn)}, 
	                         drkm::NamedTuple{(:img, :imgn)}) 
	img = imgm.img .- drkm.img
	imgn = img ./maximum(img)
    (img = img, imgn = imgn, dark=drkm.img)
end

# ╔═╡ 4760fdb6-5a0b-4ba2-89b7-0cc7f764d68e
begin
	image, dimage = select_image(xfiles, xfdrk, nxdrk,setup)
	cimg = get_corrected_image(image, dimage)
	mosaicview(Gray.(image.imgn), Gray.(cimg.imgn); nrow = 1)
end

# ╔═╡ 8467180b-dd7a-4e99-931a-cbb87af90fc2
let
	spngn = string(setup.series, "_", setup.measurement, "_",
		           setup.point,  "_Filter_", setup.filter, "_imageCMOS.png")
	pxth = joinpath(pngdir, spngn)	
	save(pxth, colorview(Gray, map(clamp01nan, cimg.imgn) ))
end

# ╔═╡ cdf03376-4450-4b71-a820-33564d1ed71d
begin
	imgmax, imgpos = lfi.LaserLab.signal_around_maximum(cimg.img, cimg.dark; nsigma=nsigma)
	if typeof(imgmax.imgn) == Matrix{Float64}
		Gray.(imgmax.imgn)
	end
end

# ╔═╡ 986b6b9e-3a02-4703-91d0-8088a8066810
let
	spngn = string(setup.series, "_", setup.measurement, "_",
		           setup.point,  "_Filter_", setup.filter, "_imageZoomCMOS.png")
	pxth = joinpath(pngdir, spngn)	
	if typeof(imgmax.imgn) == Matrix{Float64}
		save(pxth, colorview(Gray, map(clamp01nan, imgmax.imgn) ))
	end
end

# ╔═╡ aba1623c-1fd1-4210-b587-a76e21522397
begin
	stot = sum(imgmax.img)
md"""
- Size of the image =$(size(imgmax.img))
- Total sum around the peak = $(round(stot, sigdigits=3))
- value of max = $(imgpos.max)
- position of max: i = $(imgpos.imax)  j = $(imgpos.jmax)
"""
end

# ╔═╡ a62b6fd5-959c-421a-a160-c420dae4ca99
function spectrum_max(setup::Setup,  
                      xfn::Vector{String}, filtnm::NamedTuple, adctopes::Float64;
					  nsigma::Float64=3.0, drkpnt="Dark")
	
	ZMX = Vector{Float64}()
	ZSM = Vector{Float64}()
	ZI = Vector{Float64}()
	ZJ = Vector{Float64}()

	xfiles, _    = select_files(setup.cmdir,setup.series,
			                  setup.measurement, setup.point)
	xfdrk, nxdrk = select_files(setup.cmdir,setup.series,
			                  setup.measurement, drkpnt)
	for fltr in xfn
		setupf = Setup(setup.cmdir, setup.series, setup.measurement,
			           setup.point, fltr, setup.rep)

		imgm, drkm = select_image(xfiles, xfdrk, nxdrk, setupf)
		cimgz = get_corrected_image(imgm, drkm)
		
		imgmz, imgp = lfi.LaserLab.signal_around_maximum(cimgz.img, 
			                                             cimgz.dark; nsigma=nsigma)
		push!(ZMX,imgp.max)
		push!(ZSM,sum(imgmz.img))
		push!(ZI,imgp.imax)
		push!(ZJ,sum(imgp.jmax))

        #println(" for filter ",fltr, " sum = ", sum(imgmz.img))
	end
    
	DataFrame("fltn" => xfn, "cflt" => filtnm.center, "lflt" => filtnm.left, "rflt" => filtnm.right, "wflt" => filtnm.width,
		      "sum"=> ZSM, "sumpes"=> adctopes *(ZSM ./filtnm.width), "max"=> ZMX, "imax" => ZI, "jmax" => ZJ)	
end

# ╔═╡ d5540947-8b91-4cba-9738-c707e9945eab
begin
	sdf = spectrum_max(setup, xfn, filtnm, adctopes; nsigma=nsigma)

	psing = plot(sdf.cflt, sdf.sumpes, lw=2, label=setup.point, title="Spectrum around maximum")
	scatter!(sdf.cflt, sdf.sumpes, label="")
	xlabel!("λ (nm)")
	ylabel!("pes")
end

# ╔═╡ 03c6566e-f4d3-47b6-8a20-4744efc541a0
begin
	spsngn = string(setup.series, "_", setup.measurement, "_", setup.point, ".png")
	pxxth = joinpath(pngdir, spsngn)	
	png(psing, pxxth)
end

# ╔═╡ ac542727-b476-437e-9bc8-8834a0653355
function spectrum_sum(setup::Setup, 
                      xfn::Vector{String}, filtnm::NamedTuple, adctopes::Float64)
    
    ZSM = Vector{Float64}()
	xfiles, _    = select_files(setup.cmdir,setup.series,
			                  setup.measurement, setup.point)
	xfdrk, nxdrk = select_files(setup.cmdir,setup.series,
			                  setup.measurement, drkpnt)
    #println("in spectrum_sum: setup = ", setup)
    for fltr in xfn
        setupf = Setup(setup.cmdir, setup.series, setup.measurement,
			           setup.point, fltr, setup.rep)

		imgm, drkm = select_image(xfiles, xfdrk, nxdrk, setupf)
		cimgz = get_corrected_image(imgm, drkm)
        push!(ZSM,sum(cimgz.img))
    end

    #println("in spectrum_sum (after): setup = ", setup)
    DataFrame("fltn" => xfn, 
		      "cflt" => filtnm.center, 
		      "lflt" => filtnm.left, 
		      "rflt" => filtnm.right, 
		      "wflt" => filtnm.width,
              "sum"=> ZSM, 
		      "sumpes"=> adctopes *(ZSM ./filtnm.width))	
end

# ╔═╡ 6325910e-b376-42ab-962a-28f749bc27b2
if zfull
	sdf2 = spectrum_sum(setup,  xfn, filtnm, adctopes)
	plot(sdf2.cflt, sdf2.sumpes, lw=2, label=setup.point, title="spectrum full CMOS")
	scatter!(sdf2.cflt, sdf2.sumpes, label="")
	xlabel!("λ (nm)")
	ylabel!("pes")
end

# ╔═╡ 18c767aa-1461-4847-ac0d-26ad5a06dd1c
function get_outpath(setup::Setup, ext="*.csv")
	string(setup.series, "_", setup.measurement, "_", 
                   setup.point, ext)

end

# ╔═╡ effb1278-5896-4b19-a6ad-7f12bf5ba9b5
function spectrum_max_allpoints!(setup::Setup, 
                                xpt::Vector{String}, xfn::Vector{String}, 
                                filtnm::NamedTuple, adctopes::Float64; 
								nsigma::Float64=3.0, odir, 
								etype="csv", drkpnt="Dark")
    for pt in xpt
		setupp = Setup(setup.cmdir, setup.series, setup.measurement,
			           pt, setup.filter, setup.rep)
        xpath = joinpath(setupp.cmdir,setupp.series,setupp.measurement, setupp.point)
        dfiles = string("*.",etype)
		
        sdf = spectrum_max(setupp, xfn, filtnm, adctopes;
		                   nsigma=nsigma)
        sdfnm = get_outpath(setupp, ".csv")
        sdff = joinpath(odir, sdfnm)
	    println("Writing point to  =", sdff)
	    CSV.write(sdff, sdf)
    end
end

# ╔═╡ 25219398-6903-4cb0-a336-127eedbfa902
if zrec
	spectrum_max_allpoints!(setup, sspdirs, xfn, 	
		                    filtnm, adctopes; 
                            nsigma=nsigma, odir=csvdir)
end


# ╔═╡ f9608d49-3604-4c8d-913c-6cbf35f7a85f
function read_spectrum(setup::Setup, csvdir::String, ext=".csv")
    
	sdfnm = string(setup.series, "_", setup.measurement, "_", 
                   setup.point, ext)
    
                  
	sdff = joinpath(csvdir, sdfnm)
	println("reading file =", sdff)

    lfi.LaserLab.load_df_from_csv(csvdir, sdfnm, lfi.LaserLab.enG)
	
end

# ╔═╡ 8b09554f-bf5f-4cc8-ab16-9ac34036f111
function spectrum_fromfile_allpoints(setup::Setup, xpt::Vector{String}, csvdir)
    dfdict = Dict()
    for pt in xpt
        setupp = Setup(setup.cmdir, setup.series, setup.measurement,
			           pt, setup.filter, setup.rep)
        df = read_spectrum(setupp, csvdir)
        dfdict[pt] = df
    end

    dfdict
end

# ╔═╡ 9ad17eec-bb5c-4980-9c32-27e01c5b7fcf
function plot_spectrum_for_point(sdfp, pt, fscale="cflt", escale="sumpes")
	plt = plot(sdfp[!, fscale], sdfp[!, escale], lw=2, label="", 
		       xtickfontsize=8,ytickfontsize=8)
	scatter!(sdfp[!, fscale], sdfp[!, escale], label="")
	xlabel!("λ (nm)")
	ylabel!("pes")
	#yticks!([2e+5,4e+5,6e+5])
    xticks!([0,400,600,800])
	
	#xtickfontsize=18,ytickfontsize=18,xlabel="wavelength",xguidefontsize=18,yscale=:log10,ylabel="flux",yguidefontsize=18,legendfontsize=18) here

	plt
end

# ╔═╡ 1f1b334e-941b-4fbd-b964-b4c098ef3231
if zread
	dfdict = spectrum_fromfile_allpoints(setup, sspdirs, csvdir);
	PLT=[]
	for pt in sspdirs
		sdfp = dfdict[pt]
		push!(PLT, plot_spectrum_for_point(sdfp, pt, "cflt"))
	end
	pall = plot(size=(750,750), PLT[1:end]..., layout=(nx,ny), titlefontsize=8)
	
	
end

# ╔═╡ f8b65718-3e1b-454a-817f-1e78feb43225
if zread
	let
		spngn = string(setup.series, "_", setup.measurement, ".png")
		pxth = joinpath(pngdir, spngn)	
		png(pall, pxth)
	end
end

# ╔═╡ Cell order:
# ╠═981730a6-61fc-484b-ba3c-66920ee7cf83
# ╠═8e7ec382-c738-11ec-3aae-b50d60f15c4f
# ╠═06b8ed45-43bc-464f-89c0-dc0406312b81
# ╠═8833b198-03e4-4679-8949-0c76546cb847
# ╠═6163ba69-1237-4b49-988e-9a73cfef67f6
# ╠═5edc41bc-b912-44bf-9be5-a013f27a75ab
# ╠═c9aaf1cc-80c4-475b-8a81-e00918d91b1e
# ╠═7ce42aec-b319-4de9-b70c-84046d45a600
# ╠═58269465-ba8c-4840-bbfc-0a27897f3e2a
# ╠═0b1c5662-ec4f-486e-9ee6-7fa6ba953e45
# ╠═b98ad447-6055-46e5-bb4f-9e67f9c3176a
# ╠═88853edb-dc1f-4e7a-a0ba-1868276d1ada
# ╠═dde07319-e6a0-48ac-94d6-92e0fda4df41
# ╠═f6dcfbd5-8416-42f5-b029-c794f92ee413
# ╠═88c44435-9685-41c3-a6c8-ffaf7498d60f
# ╠═57d96432-4318-4291-8255-bfa5d6d3635c
# ╠═0b4d2c08-a677-492c-b5da-982d3d5096fc
# ╠═3ae7315b-2056-4e3f-9d41-b6e0292f30e3
# ╠═161b12b6-88a0-4d4d-bfc5-01310534cbdc
# ╠═e87f48e3-5e5a-44d5-83de-c520e522e33a
# ╠═50ea2ecc-970f-4630-8c7e-acf5e69cc4c9
# ╠═d389f99a-14c2-408f-ad7b-838e00225357
# ╠═f13173e1-088c-4f5c-ae6a-f9aebcd6bc57
# ╠═8dbf64ec-5854-44b1-ac73-7cd0363a1c6d
# ╠═479f8f86-372c-4b91-9f73-e57a85d3d194
# ╠═b857231c-f06c-426d-8f2c-5f8007134714
# ╠═d762313e-1466-4724-9c35-dd89e657a11c
# ╟─7e185082-f3bf-4d95-9d5d-57101f47a684
# ╠═43ccd7f3-140f-4cb8-af41-e13534c454f3
# ╠═4760fdb6-5a0b-4ba2-89b7-0cc7f764d68e
# ╠═8467180b-dd7a-4e99-931a-cbb87af90fc2
# ╠═cdf03376-4450-4b71-a820-33564d1ed71d
# ╠═986b6b9e-3a02-4703-91d0-8088a8066810
# ╠═aba1623c-1fd1-4210-b587-a76e21522397
# ╠═e9e8e58a-e1db-49f1-8429-420271fb1852
# ╠═54bd1f6c-2b10-47a1-838f-b428fe6b7635
# ╠═d5540947-8b91-4cba-9738-c707e9945eab
# ╠═03c6566e-f4d3-47b6-8a20-4744efc541a0
# ╠═5477b1dc-ef52-4426-aa47-f513854fcdae
# ╠═6325910e-b376-42ab-962a-28f749bc27b2
# ╠═a90c8edc-1405-4a8c-aacf-53cd130910ae
# ╠═a48af8f4-4ed2-45cd-b4e8-9b3106c885f3
# ╠═1794afb6-6ef0-46d6-b182-d54362b9a07d
# ╠═25219398-6903-4cb0-a336-127eedbfa902
# ╟─82ddd81b-8aea-4711-97d9-a645121786f8
# ╟─155d5066-935b-4643-8aad-f4c635aa7eec
# ╠═1f1b334e-941b-4fbd-b964-b4c098ef3231
# ╠═f8b65718-3e1b-454a-817f-1e78feb43225
# ╠═c1b2cf36-d72e-4348-911e-8e44f1834ae4
# ╠═f2714261-7bb0-47d7-8aac-c16bb5d1f891
# ╠═78de6bcd-4173-40e3-b500-499568289ba1
# ╠═0d5a7021-6072-464e-836e-f05b1e178b80
# ╠═a759ecf7-7373-46bf-ab15-df39cbfa6814
# ╠═46b1d54e-4ebf-45a2-b137-9690b7a51d38
# ╠═2498aa42-2c24-47e3-bf5b-647377af0dbc
# ╠═afeb42ca-b462-4320-b364-98a0b4730e33
# ╠═f9ea012f-9aae-4a8b-89d9-0f86802bb14f
# ╠═58f3c3ac-698b-4381-bdcb-ca8a9cadc8d8
# ╠═a62b6fd5-959c-421a-a160-c420dae4ca99
# ╠═ac542727-b476-437e-9bc8-8834a0653355
# ╠═18c767aa-1461-4847-ac0d-26ad5a06dd1c
# ╠═effb1278-5896-4b19-a6ad-7f12bf5ba9b5
# ╠═f9608d49-3604-4c8d-913c-6cbf35f7a85f
# ╠═8b09554f-bf5f-4cc8-ab16-9ac34036f111
# ╠═9ad17eec-bb5c-4980-9c32-27e01c5b7fcf
