### A Pluto.jl notebook ###
# v0.19.3

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

# ╔═╡ f6dcfbd5-8416-42f5-b029-c794f92ee413
nsigma=5.0

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
	cmdir ="/Users/jj/JuliaProjects/LaserLab/labdata/CMOS"
	odir  = "/Users/jj/JuliaProjects/LaserLab/data/CMOS"
	dfiles ="*.csv"
	dplots ="*.png"
	dmtype = ["Imag","Dark"]
	md"""
	CMOS dir = $cmdir
	"""
end

# ╔═╡ e87f48e3-5e5a-44d5-83de-c520e522e33a
let
	readdir(cmdir)
	dirs = lfi.LaserLab.getbolddirs(cmdir)
	md""" Select experiment : $(@bind sexp Select(dirs))"""
end

# ╔═╡ a759ecf7-7373-46bf-ab15-df39cbfa6814
begin
	namex = split(sexp,"_")
	csvdir = joinpath(odir, string(namex[1], namex[2]),"csv")
	pngdir = joinpath(odir, string(namex[1], namex[2]),"png")
	md"""
	- csv dir = $csvdir
	- png dir = $pngdir
	"""
end

# ╔═╡ 1924796e-1cb9-47ad-bfb6-f01e544bdbc1
namex

# ╔═╡ 50ea2ecc-970f-4630-8c7e-acf5e69cc4c9
let
	path = joinpath(cmdir,sexp)
	dirs = lfi.LaserLab.getbolddirs(path)
	md""" Select run : $(@bind srun Select(dirs))"""
end

# ╔═╡ 4a1c449f-2656-4f05-9bea-783a3c7a4d30
begin
	ppath = joinpath(cmdir,sexp,srun)
	pdirs = lfi.LaserLab.getbolddirs(ppath)
	spdirs = sort(parse.(Int64, [split(pd, "Point")[2] for pd in pdirs]))
	sspdirs = [string("Point", string(i)) for i in spdirs]
	md""" Select point : $(@bind spoint Select(sspdirs))"""
end

# ╔═╡ df94ef34-08c9-4711-ae16-b024cf84b82a
md"""
## Define repetition, filter and image type
"""

# ╔═╡ 479f8f86-372c-4b91-9f73-e57a85d3d194
begin
	xpath = joinpath(cmdir,sexp,srun, spoint)
	xfiles = Glob.glob(dfiles, xpath)
	nxfiles = string.([split(f,"/")[end] for f in xfiles])
	rn = lfi.LaserLab.findpattern(nxfiles, "rep")
	md""" Select repetition number: $(@bind srn Select(rn))"""
end

# ╔═╡ 05bfacd2-047e-4491-a168-96120d394a83
typeof(xfiles)

# ╔═╡ 0b162778-31bd-4066-a576-4066acf0c370
typeof(nxfiles)

# ╔═╡ f13173e1-088c-4f5c-ae6a-f9aebcd6bc57
begin
	xfa = lfi.LaserLab.findpattern(nxfiles, "Filter")
	xfb = sort(parse.(Int64, xfa))
	xfn = string.(xfb)
	md""" Select filter: $(@bind sfn Select(xfn))"""
end

# ╔═╡ b19277b6-bdbc-46ce-ab18-d50e08fae24a
md""" Select Img/dark: $(@bind wid Select(dmtype))"""

# ╔═╡ d762313e-1466-4724-9c35-dd89e657a11c
begin
	setup = lfi.LaserLab.Setup(cmdir,
		                       string(sexp), string(srun), string(spoint),        
		                       string(sfn),  string(srn),  string(wid), texp,
							   "before", "Filter", "rep", "Imag_1", "Dark",
                               string("ExpoT_",string(Int64(texp)),"s"),
							   ".csv", ".png")
	
	pngpath = joinpath(pngdir,lfi.LaserLab.get_outpath(setup, "png"))
	md"""
	Setup
	- root dir = $(setup.cmdir)
	- Series = $(setup.series)
	- Measurement = $(setup.measurement)
	- Point = $(setup.point)
	- Filter = $(setup.filter)
	- Repetition = $(setup.rep)
	- Img/Dark = $(setup.imgid)
	- Time exposition = $(setup.texpsec)

	Grammar
	
	- stype = $(setup.darktype)
	- Filter = $(setup.sfilter)
	- Repetition = $(setup.srep)
	- Img = $(setup.img)
	- Dark = $(setup.dark)
	- Time exposition = $(setup.exp)
	- File extensions = $(setup.fext)
	- png extensions = $(setup.pext)

	Paths
	- png = $pngdir
	"""
end

# ╔═╡ 2e7d4d21-69be-45e8-8d30-440c215ffd3e
setup

# ╔═╡ 7e185082-f3bf-4d95-9d5d-57101f47a684
md"""
## Select image

- Image in the left not corrected 
- Image in the right corrected (dark current subtracted)
"""

# ╔═╡ 4760fdb6-5a0b-4ba2-89b7-0cc7f764d68e
begin
	image = lfi.LaserLab.select_image(xfiles, setup, setup.imgid)
	cimg = lfi.LaserLab.get_corrected_image(xfiles, setup)
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
	Gray.(imgmax.imgn)
end

# ╔═╡ 986b6b9e-3a02-4703-91d0-8088a8066810
let
	spngn = string(setup.series, "_", setup.measurement, "_",
		           setup.point,  "_Filter_", setup.filter, "_imageZoomCMOS.png")
	pxth = joinpath(pngdir, spngn)	
	save(pxth, colorview(Gray, map(clamp01nan, imgmax.imgn) ))
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

# ╔═╡ e9e8e58a-e1db-49f1-8429-420271fb1852
md"""
## Comparison between spectrum computed around maximum and using the full camera. 

- Notice the umphysical shoulder when using the full camera
"""

# ╔═╡ d5540947-8b91-4cba-9738-c707e9945eab
begin
	sdf = lfi.LaserLab.spectrum_max(xfiles, setup, xfn, filtnm, adctopes; nsigma=3.0)

	plot(sdf.cflt, sdf.sumpes, lw=2, label=setup.point, title="Spectrum around maximum")
	scatter!(sdf.cflt, sdf.sumpes, label="")
	xlabel!("λ (nm)")
	ylabel!("pes")
end

# ╔═╡ af081fae-01f0-44ce-9f40-9cb5d3f28695
setup

# ╔═╡ 6325910e-b376-42ab-962a-28f749bc27b2
begin
	sdf2 = lfi.LaserLab.spectrum_sum(xfiles, setup,  xfn, filtnm, adctopes)
	plot(sdf2.cflt, sdf2.sumpes, lw=2, label=setup.point, title="spectrum full CMOS")
	scatter!(sdf2.cflt, sdf2.sumpes, label="")
	xlabel!("λ (nm)")
	ylabel!("pes")
end

# ╔═╡ a90c8edc-1405-4a8c-aacf-53cd130910ae
setup

# ╔═╡ 56ccf00b-8797-4347-b4e9-026a9da14a81
setup.point

# ╔═╡ a48af8f4-4ed2-45cd-b4e8-9b3106c885f3
md"""
# Analysis for al points
"""

# ╔═╡ 1794afb6-6ef0-46d6-b182-d54362b9a07d
md""" Check to carry analysis for all points: $(@bind zrec CheckBox())"""

# ╔═╡ 25219398-6903-4cb0-a336-127eedbfa902
if zrec
	lfi.LaserLab.spectrum_max_allpoints!(setup, sspdirs, xfn,  
		                                 filtnm, adctopes; 
                                         nsigma=nsigma, odir=csvdir, etype="csv")
end


# ╔═╡ 82ddd81b-8aea-4711-97d9-a645121786f8
md""" Check to read and plot data for all points: $(@bind zread CheckBox())"""

# ╔═╡ 0d8e9320-ae68-461f-87a5-8563010a930e
spngn = string(setup.series, "_", setup.measurement, ".png")

# ╔═╡ 20b26dae-9bc4-474e-9b2d-c397f6305a96
joinpath(pngdir, spngn)

# ╔═╡ c1b2cf36-d72e-4348-911e-8e44f1834ae4
md"""
# Functions
"""

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
	dfdict = lfi.LaserLab.spectrum_fromfile_allpoints(setup, sspdirs, csvdir);
	PLT=[]
	for pt in sspdirs
		sdfp = dfdict[pt]
		push!(PLT, plot_spectrum_for_point(sdfp, pt, "cflt"))
	end
	pall = plot(size=(750,750), PLT[1:end-1]..., layout=(3,3), titlefontsize=8)
	
	
end

# ╔═╡ 25a42a00-8ec0-46d8-b9f6-a634475cefd3
length(PLT)

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
# ╟─8833b198-03e4-4679-8949-0c76546cb847
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
# ╠═57d96432-4318-4291-8255-bfa5d6d3635c
# ╠═0b4d2c08-a677-492c-b5da-982d3d5096fc
# ╠═e87f48e3-5e5a-44d5-83de-c520e522e33a
# ╠═a759ecf7-7373-46bf-ab15-df39cbfa6814
# ╠═1924796e-1cb9-47ad-bfb6-f01e544bdbc1
# ╠═50ea2ecc-970f-4630-8c7e-acf5e69cc4c9
# ╠═4a1c449f-2656-4f05-9bea-783a3c7a4d30
# ╠═df94ef34-08c9-4711-ae16-b024cf84b82a
# ╠═479f8f86-372c-4b91-9f73-e57a85d3d194
# ╠═05bfacd2-047e-4491-a168-96120d394a83
# ╠═0b162778-31bd-4066-a576-4066acf0c370
# ╠═f13173e1-088c-4f5c-ae6a-f9aebcd6bc57
# ╠═b19277b6-bdbc-46ce-ab18-d50e08fae24a
# ╠═d762313e-1466-4724-9c35-dd89e657a11c
# ╠═2e7d4d21-69be-45e8-8d30-440c215ffd3e
# ╟─7e185082-f3bf-4d95-9d5d-57101f47a684
# ╠═4760fdb6-5a0b-4ba2-89b7-0cc7f764d68e
# ╠═8467180b-dd7a-4e99-931a-cbb87af90fc2
# ╠═cdf03376-4450-4b71-a820-33564d1ed71d
# ╠═986b6b9e-3a02-4703-91d0-8088a8066810
# ╠═aba1623c-1fd1-4210-b587-a76e21522397
# ╠═e9e8e58a-e1db-49f1-8429-420271fb1852
# ╠═d5540947-8b91-4cba-9738-c707e9945eab
# ╠═af081fae-01f0-44ce-9f40-9cb5d3f28695
# ╠═6325910e-b376-42ab-962a-28f749bc27b2
# ╠═a90c8edc-1405-4a8c-aacf-53cd130910ae
# ╠═56ccf00b-8797-4347-b4e9-026a9da14a81
# ╠═a48af8f4-4ed2-45cd-b4e8-9b3106c885f3
# ╠═1794afb6-6ef0-46d6-b182-d54362b9a07d
# ╠═25219398-6903-4cb0-a336-127eedbfa902
# ╠═82ddd81b-8aea-4711-97d9-a645121786f8
# ╠═1f1b334e-941b-4fbd-b964-b4c098ef3231
# ╠═25a42a00-8ec0-46d8-b9f6-a634475cefd3
# ╠═0d8e9320-ae68-461f-87a5-8563010a930e
# ╠═20b26dae-9bc4-474e-9b2d-c397f6305a96
# ╠═f8b65718-3e1b-454a-817f-1e78feb43225
# ╠═c1b2cf36-d72e-4348-911e-8e44f1834ae4
# ╠═9ad17eec-bb5c-4980-9c32-27e01c5b7fcf
