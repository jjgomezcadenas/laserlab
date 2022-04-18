### A Pluto.jl notebook ###
# v0.18.0

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

# ╔═╡ 65beba5e-776d-40d5-96a6-d0a2f3a9a4ea
using Pkg; Pkg.activate("/Users/jj/JuliaProjects/LaserLab/")

# ╔═╡ 8a67d331-35ae-44d6-942c-a61db4e7afa4
begin
	using PlutoUI
	using CSV
	using DataFrames
	using Images
	using Images
	using TestImages
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
	import Glob
end

# ╔═╡ cbac975d-cda5-43fc-bd22-010ce10e8f92


# ╔═╡ d1b57350-bd9d-11ec-21f4-3d7a8656a31c
import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
    A, N, mol, mmol, V, L, M

# ╔═╡ 53943a33-9cf3-4c9f-add5-278fdc3bc0fe
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


# ╔═╡ b829dbf6-486e-4bd1-914c-54fe6d389b9c
lfi = ingredients("../src/LaserLab.jl")


# ╔═╡ 7aafa82d-4fe3-4bdc-a14c-ec48a7178260
PlutoUI.TableOfContents(title="G2 BaCl Bold73 analysis", indent=true)

# ╔═╡ 839325d1-520a-472d-af4f-b378938ba4c9
timg = testimage("cameraman")

# ╔═╡ b9e1b770-a973-409a-a26c-b1b863a57120
img_edge = lfi.LaserLab.sujoy(timg, four_connectivity=true)

# ╔═╡ 9a5ca198-b6a7-460e-a7f9-aab255815ecc
img_edge₀₁ = binarize(img_edge, Otsu()) # or use other binarization methods provided in ImageBinarization


# ╔═╡ 3f3bf3c3-1999-4ab1-892a-00c56b4fe098
mosaicview(timg, img_edge, img_edge₀₁; nrow = 1)

# ╔═╡ 5a5d50e6-32db-4404-b758-3f045c4f3389
img_edgem = convert(Array{Float64}, img_edge₀₁)

# ╔═╡ 9da70d55-28b9-4c12-ba6a-c560832a5194
md"""
# Read Data
"""

# ╔═╡ 6502fdee-cd60-4d32-ac77-c42c845b12ba
begin
	proot ="/Users/jj/JuliaProjects/LaserLab/labdata"
	sroot = "/Users/jj/JuliaProjects/LaserLab/data/G2Sl"
	psample = "G2_BOLD_073_A1"
	prun = "G2SL-BaCl2_20220408"
	path = joinpath(proot,psample, prun)
	pfiles ="*.dat"
	fls = joinpath(path,pfiles)
	nfiles = length(fls) 
	sflts = string.(collect(2:11));
	sreps = string.(collect(5:10));
	snum = ["1", "2", "c"];
	stypes=["Dark", "Imag"];
md"
- Path for glob: $path
- files : $fls
- number of files: $nfiles
"
end

# ╔═╡ dd51775e-a254-49da-9062-c166ec4a0922
files = Glob.glob(pfiles, path);

# ╔═╡ 02e2aa4a-40ff-499b-ab31-1004023b5faa
names = [split(f,"/")[end] for f in files]

# ╔═╡ a2ce79cf-a447-4cab-8843-a2cd691b4004
md"""
# Select Filter and Rep
"""

# ╔═╡ 90bda847-8afe-4968-9a02-50ad7c7cc572
md""" Select Filter number (1-10 if 0 add all): $(@bind whichf Select(sflts))"""

# ╔═╡ 1cd87997-3ff3-40fa-98b7-ee9655848a76
md""" Select repetition number (1-10 if 0 first): $(@bind whichr Select(sreps))"""

# ╔═╡ 6d08fd5c-6757-43dd-a0e3-bb161e78cd26
#pimgx, stot = findpeaks(simgm; xrange=50);

# ╔═╡ f56df91d-53ce-4ca7-a9c8-7740cc3bc430
#maxx, indxt = findmax(pimgx);

# ╔═╡ 3444eb77-b397-4ac3-99a2-972abc97090d
#begin
#	pimgn = pimgx ./maximum(pimgx);
#	gimg2  =Gray.(pimgn);
#	plot(gimg2,layout=(1,1), titlefontsize=10)
#end

# ╔═╡ 64ed6a8c-1bca-4c19-83ff-d507ba4359a5
#sdrkmn = sdrkm ./maximum(sdrkm);

# ╔═╡ 1313a166-29f3-4f5e-8e79-dbfab22a69d9
#begin
#	gdrk = Gray.(sdrkmn)
#	plot(gdrk, layout=(1,1), titlefontsize=10)
#end

# ╔═╡ e6816855-d168-4317-9ddc-8e6b13c1d0a1
#pimgd, stotd = findpeaks(sdrkm; xrange=30);

# ╔═╡ a9d90521-a7e5-4cd2-92ca-45de8c03e514
#maxd, indxd = findmax(pimgd)

# ╔═╡ 55a34a0a-ae88-4b80-b6e1-515da142faa8
#thrz = maxd + sqrt(maxd)

# ╔═╡ 87369d14-7318-417c-bc33-6a1bc807c0ee
#begin
	#pimgdn = pimgd ./maximum(pimgd)
	#gimgd  =Gray.(pimgdn)
	#plot(gimgd, layout=(1,1), titlefontsize=10)
#end

# ╔═╡ 43fd3321-83ce-4a9b-8632-74dd3cd68ddd
#maximum(pimgd)

# ╔═╡ 8feb58fd-3cb7-4860-bc22-52bcee93d40a
scs = [string(c) for c in collect(5:10)]

# ╔═╡ 0eda3bbd-a4c2-43f4-96d8-e8c2f255e56f
@bind zedegr CheckBox()

# ╔═╡ 01aeef0d-12e1-4e32-974f-f5405a661d95
md"""
# Reconstruct spectra
"""

# ╔═╡ 18a0121f-b0a5-46e5-8178-9bc800476b65
@bind zrepx CheckBox()

# ╔═╡ 3aa1788a-131a-402f-9493-32d353020219
fdf = lfi.LaserLab.load_df_from_csv(sroot, "G2Bold073A1Fluo.csv", lfi.LaserLab.enG)

# ╔═╡ 8740f867-109f-47c1-9e90-57ece3b83886
fbadf = lfi.LaserLab.load_df_from_csv(sroot, "G2BaCl2Bold073A1Fluo.csv", lfi.LaserLab.enG)

# ╔═╡ 8c52eac3-388f-49ea-bcc4-6d9a42085c52
begin
	plot(collect(1:10), fdf.r1, lw=2, label="rep5 S")
	plot!(collect(1:10), fdf.r2, lw=2, label="rep6 S")
	plot!(collect(1:10), fdf.r3, lw=2, label="rep7 S")
	plot!(collect(1:10), fdf.r4, lw=2, label="rep8 S")
	plot!(collect(1:10), fdf.r5, lw=2, label="rep9 S")
	plot!(collect(1:10), fdf.r6, lw=2, label="rep10 S")
	plot!(collect(1:10), fbadf.r1, lw=2, label="rep5 Ba2+")
	plot!(collect(1:10), fbadf.r2, lw=2, label="rep6 Ba2+")
	plot!(collect(1:10), fbadf.r3, lw=2, label="rep7 Ba2+")
	plot!(collect(1:10), fbadf.r4, lw=2, label="rep8 Ba2+")
	plot!(collect(1:10), fbadf.r5, lw=2, label="rep9 Ba2+")
	plot!(collect(1:10), fbadf.r6, lw=2, label="rep10 Ba2+")
end

# ╔═╡ cec73687-de04-46e6-a73e-a4d1d3beafc5
md"# Functions"

# ╔═╡ 6fdaa0d2-b6af-40ca-b8e2-5cb8176bbadb
function select_image(inames::Vector{String}, fnumber::String, repetition::String,
	                  itype::String="Imag", inum::String="1")
	
	names = [split(f,"/")[end] for f in inames]
	if repetition == "0"
		fhead = string("Filter_", fnumber, "_ExpoT_1000000ms_")
	else
		fhead = string("Filter_", fnumber, "_rep_",repetition,"_ExpoT_1000000ms_")
	end
	
	if itype == "Dark"
		ftail = string(itype, "_", inum, ".dat")
	elseif itype == "Imag" && inum=="1"
		ftail = string("Imag_1.dat")
	elseif itype == "Imag" && inum=="c"
		ftail = string("Imag_corrected.dat")
	else
		@info "not implemented"
	end

	ff = string(fhead, ftail)
	indx = findall([occursin(ff, name) for name in names])[1]
	inames[indx]
end

# ╔═╡ edabab28-72bb-43f2-b90b-7b6ba939f2eb
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

# ╔═╡ d8caa384-6cc8-4a96-93ee-3504a95c5214
function sum_edge(img::Matrix{Float64}, edge::Matrix{Float64})
	sum=0.0
	nedge = 0
	for i in 1:size(img)[1]
		for j in 1:size(img)[2]
			if edge[i,j] > 0.0 
				sum += img[i,j]
				nedge += 1
			end
		end
	end
	sum, nedge
end

# ╔═╡ b8af0d87-f9e1-4ea7-9fa0-1c675233a8d3
function findpeaks(img::Matrix{Float64}; xrange::Int64=10)
	maxx, indxt = findmax(img)
	ixmin = max(1, indxt[1] - xrange)
	ixmax = min(size(img)[1], indxt[1] + xrange -1)
	iymin = max(1, indxt[2] - xrange)
	iymax = min(size(img)[2], indxt[2] + xrange -1)

	img2 = Array{Float64}(undef, ixmax-ixmin+1, iymax-iymin+1);
	
	stot = 0.0
	for (i, ix) in enumerate(ixmin:ixmax)
		for (j, iy) in enumerate(iymin:iymax)
			img2[i,j]=img[ix,iy]
			stot += img[ix,iy]
		end
	end
	img2, stot 
end

# ╔═╡ 5dc25eb7-7841-4c04-aa5e-7d63acc5a312
function get_images(files::Vector{String}, whichf::String, whichr::String)
	ffimg = select_image(files, whichf, whichr, "Imag", "1");
	imgdf = DataFrame(CSV.File(ffimg, header=false,delim="\t"));
	imgm  = dftomat(imgdf)
	
	ffdrk = select_image(files, whichf, whichr, "Dark", "1");
	drkdf = DataFrame(CSV.File(ffdrk, header=false,delim="\t"));
	drkm  = dftomat(drkdf)
	
	ffdrk2 = select_image(files, whichf, whichr, "Dark", "2");
	drkdf2 = DataFrame(CSV.File(ffdrk2, header=false,delim="\t"));
	sdrk   = drkdf .- drkdf2
	sdrkm  = dftomat(sdrk)
	
	simgdf = imgdf .- drkdf
	simgm  = dftomat(simgdf)
	
	
	imgm, drkm, sdrkm, simgm
end

# ╔═╡ 21f71c85-9e22-4834-a7d7-76e3d1ab681b
imgm, drkm, sdrkm, simgm = get_images(files, whichf, whichr);

# ╔═╡ 4aab4565-280c-432c-8f6d-6774a4dbe8be
begin
	simgmn = simgm ./maximum(simgm);
	gimg = Gray.(simgmn);
	simgmn_edge = lfi.LaserLab.sujoy(gimg, four_connectivity=true)
	simg_edge = binarize(simgmn_edge, Otsu())
	mosaicview(gimg, simg_edge; nrow = 1)
end

# ╔═╡ 74b2c6cc-da87-4de2-bb24-499a0fa794d6

simg_edgem = convert(Array{Float64}, simg_edge);

# ╔═╡ 3eda9d9c-6489-4f13-8f86-191c51fd8ec7
stot, ntot = sum_edge(simgm, simg_edgem)

# ╔═╡ b776cef6-90d6-4897-bd6e-b81a3a58b416
md"""
- Filter = $whichf
- repetition = $whichr
- **signal**
- Total light in spot (DC subtracted) = $(round(stot, sigdigits=3))
- Number of points in edge = $ntot
"""

# ╔═╡ f41cf5c7-d208-40f9-a4ee-5e992f73ca36
function recospec(files::Vector{String}; rep::String="5", xrange::Int64=30)
	STOT = []
	MAX  = []
	IX   = []
	IY   = []
	STOTd = []
	MAXd  = []
	IXd  = []
	IYd   = []
	for flt in string.(collect(2:11))
		imgm, drkm, sdrkm, simgm = get_images(files, flt, rep)
		pimgx, stot = findpeaks(simgm; xrange=xrange)
		maxx, indxt = findmax(pimgx)
		pimgd, stotd = findpeaks(sdrkm; xrange=xrange)
		maxd, indxd = findmax(pimgd)
		push!(STOT, stot)
		push!(MAX, maxx)
		push!(IX, indxt[1])
		push!(IY, indxt[2])
		push!(STOTd, stotd)
		push!(MAXd, maxd)
		push!(IXd, indxd[1])
		push!(IYd, indxd[2])
	end
	Dict("stot" =>STOT, "smax"=>MAX, "six"=>IX,"siy"=>IY,
	     "dtot" =>STOTd, "dmax"=>MAXd, "dix"=>IXd,"diy"=>IYd)
end

# ╔═╡ 2c71142a-34b2-4804-9af7-dc8147ab30a3
if zrepx
	repx = [recospec(files; rep=c, xrange=30) for c in scs]
	
	fdfx    = DataFrame("r1" => repx[1]["stot"],
		               "r2" => repx[2]["stot"],
		               "r3" => repx[3]["stot"],
		               "r4" => repx[4]["stot"],
		               "r5" => repx[5]["stot"],
		               "r6" => repx[6]["stot"])
	
	dfname = joinpath(sroot, "G2BaCl2Bold073A1Fluo.csv")
	CSV.write(dfname, fdfx)

	plot(collect(1:10),  fdfx.r1, lw=2, label="rep5 S")
	plot!(collect(1:10), fdfx.r2, lw=2, label="rep6 S")
	plot!(collect(1:10), fdfx.r3, lw=2, label="rep7 S")
	plot!(collect(1:10), fdfx.r4, lw=2, label="rep8 S")
	plot!(collect(1:10), fdfx.r5, lw=2, label="rep9 S")
	plot!(collect(1:10), fdfx.r6, lw=2, label="rep10 S")
	
end

# ╔═╡ b017b44e-9860-4f21-9b99-5a0cdf1b265c
function get_signal_image(files::Vector{String}, whichf::String, whichr::String)
	ffimg = select_image(files, whichf, whichr, "Imag", "1");
	imgdf = DataFrame(CSV.File(ffimg, header=false,delim="\t"));
	imgm  = dftomat(imgdf)
	
	ffdrk = select_image(files, whichf, whichr, "Dark", "1");
	drkdf = DataFrame(CSV.File(ffdrk, header=false,delim="\t"));
	drkm  = dftomat(drkdf)
	
	ffdrk2 = select_image(files, whichf, whichr, "Dark", "2");
	drkdf2 = DataFrame(CSV.File(ffdrk2, header=false,delim="\t"));
	sdrk   = drkdf .- drkdf2
	sdrkm  = dftomat(sdrk)
	
	simgdf = imgdf .- (drkdf .+ drkdf2) ./2.0
	simgm  = dftomat(simgdf)
	
	simgm
end

# ╔═╡ 02464496-4c60-4e5f-a3af-0e5c2bb0b937
function recoedge(files::Vector{String}; rep::String="5", frange=2:11)
	STOT = Vector{Float64}()
	
	for flt in string.(collect(frange))
		simgm = get_signal_image(files, flt, rep)
		simgmn = simgm ./maximum(simgm)
		gimg = Gray.(simgmn)
		simgmn_edge = lfi.LaserLab.sujoy(gimg, four_connectivity=true)
		simg_edge = binarize(simgmn_edge, Otsu())
		stot, ntot = sum_edge(simgm, simg_edgem)
		push!(STOT, stot)
	end
	STOT
end

# ╔═╡ ff5ef074-e4d2-4ad0-a93a-2c8aae422105
st5 = recoedge(files; rep="5", frange=2:11)

# ╔═╡ b67c0165-6c9e-43f7-b898-e70c533cc03d
plot(collect(1:10), st5, lw=2, label="rep5 S")

# ╔═╡ f5b8ddb4-d1e5-41ff-9c6f-19b9afc8128b
if zedegr
	stx = [recoedge(files; rep=c, frange=2:11) for c in scs]
	
	edf = DataFrame("r1" => stx[1],
		               "r2" => stx[2],
		               "r3" => stx[3],
		               "r4" => stx[4],
		               "r5" => stx[5],
		               "r6" => stx[6])
	
	dfn = joinpath(sroot, "G2BaCl2Bold073A1Edge.csv")
	CSV.write(dfn, edf)

	plot(collect(1:10),  edf.r1, lw=2, label="rep5 S")
	plot!(collect(1:10), edf.r2, lw=2, label="rep6 S")
	plot!(collect(1:10), edf.r3, lw=2, label="rep7 S")
	plot!(collect(1:10), edf.r4, lw=2, label="rep8 S")
	plot!(collect(1:10), edf.r5, lw=2, label="rep9 S")
	plot!(collect(1:10), edf.r6, lw=2, label="rep10 S")
	
end

# ╔═╡ fca14add-ed54-449b-82a3-bffd03f88cdd
function plot_images(imgm, drkm, sdrkm, simgm)
	imghm = heatmap(imgm, fill_z=imgm, title=string("rep=", whichr, " filt=", whichf, " Img" ), titlefontsize=10)

	drkhm = heatmap(drkm, title=string("rep=", whichr, " filt=", whichf, " Dark" ), titlefontsize=10)
	
	sdrkhm = heatmap(sdrkm, title=string("rep=", whichr, " filt=", whichf, "Dark - dark2" ), titlefontsize=10)

	simghm = heatmap(simgm, title=string("rep=", whichr, " filt=", whichf, "Img - Dark" ))
	
	plot(imghm, drkhm, sdrkhm, simghm, layout=(2,2), titlefontsize=10)
end

# ╔═╡ dfc2478c-9834-427d-a976-38068b9e59fc
plot_images(imgm, drkm, sdrkm, simgm)

# ╔═╡ 4480a6a1-3299-43b9-b871-0b9f525577e5
function img_stats(img)
	xmean = describe(img).mean
	xmax = describe(img).max
	xmin = describe(img).min
	hxmean = lfi.LaserLab.hist1d(xmean, 50, minimum(xmean), maximum(xmean)) 
	hxmin = lfi.LaserLab.hist1d(xmin, 50, minimum(xmin), maximum(xmin)) 
	hxmax = lfi.LaserLab.hist1d(xmax, 50, minimum(xmax), maximum(xmax)) 
	(xmean = xmean, xmin = xmin, xmax = xmax, 
	hxmean = hxmean, hxmin = hxmin, hxmax = hxmax)
end

# ╔═╡ 454cab5c-e715-46bd-a799-0e9e2e9a7f67
function img_title(label)
	string(label,"_Filter_", whichf, "_rep_", whichr, "_type_", whicht, "_number_", whichn)
end

# ╔═╡ 4c1c878c-fcfb-4069-a27a-5b1288896955
with_terminal() do
	for (i, ix) in enumerate(10:12)
		for (j, iy) in enumerate(20:22)
			#println("i = ", i, " j = ",j ," ix = ", ix, " iy = ", iy)
		end
	end
end

# ╔═╡ Cell order:
# ╠═8a67d331-35ae-44d6-942c-a61db4e7afa4
# ╠═cbac975d-cda5-43fc-bd22-010ce10e8f92
# ╠═65beba5e-776d-40d5-96a6-d0a2f3a9a4ea
# ╠═d1b57350-bd9d-11ec-21f4-3d7a8656a31c
# ╠═53943a33-9cf3-4c9f-add5-278fdc3bc0fe
# ╠═b829dbf6-486e-4bd1-914c-54fe6d389b9c
# ╠═7aafa82d-4fe3-4bdc-a14c-ec48a7178260
# ╠═839325d1-520a-472d-af4f-b378938ba4c9
# ╠═b9e1b770-a973-409a-a26c-b1b863a57120
# ╠═9a5ca198-b6a7-460e-a7f9-aab255815ecc
# ╠═3f3bf3c3-1999-4ab1-892a-00c56b4fe098
# ╠═5a5d50e6-32db-4404-b758-3f045c4f3389
# ╠═9da70d55-28b9-4c12-ba6a-c560832a5194
# ╠═6502fdee-cd60-4d32-ac77-c42c845b12ba
# ╠═dd51775e-a254-49da-9062-c166ec4a0922
# ╠═02e2aa4a-40ff-499b-ab31-1004023b5faa
# ╠═a2ce79cf-a447-4cab-8843-a2cd691b4004
# ╠═90bda847-8afe-4968-9a02-50ad7c7cc572
# ╠═1cd87997-3ff3-40fa-98b7-ee9655848a76
# ╠═21f71c85-9e22-4834-a7d7-76e3d1ab681b
# ╠═dfc2478c-9834-427d-a976-38068b9e59fc
# ╠═4aab4565-280c-432c-8f6d-6774a4dbe8be
# ╠═74b2c6cc-da87-4de2-bb24-499a0fa794d6
# ╠═3eda9d9c-6489-4f13-8f86-191c51fd8ec7
# ╠═6d08fd5c-6757-43dd-a0e3-bb161e78cd26
# ╠═f56df91d-53ce-4ca7-a9c8-7740cc3bc430
# ╠═3444eb77-b397-4ac3-99a2-972abc97090d
# ╠═64ed6a8c-1bca-4c19-83ff-d507ba4359a5
# ╠═1313a166-29f3-4f5e-8e79-dbfab22a69d9
# ╠═e6816855-d168-4317-9ddc-8e6b13c1d0a1
# ╠═a9d90521-a7e5-4cd2-92ca-45de8c03e514
# ╠═55a34a0a-ae88-4b80-b6e1-515da142faa8
# ╠═87369d14-7318-417c-bc33-6a1bc807c0ee
# ╠═43fd3321-83ce-4a9b-8632-74dd3cd68ddd
# ╠═b776cef6-90d6-4897-bd6e-b81a3a58b416
# ╠═ff5ef074-e4d2-4ad0-a93a-2c8aae422105
# ╠═b67c0165-6c9e-43f7-b898-e70c533cc03d
# ╠═8feb58fd-3cb7-4860-bc22-52bcee93d40a
# ╠═0eda3bbd-a4c2-43f4-96d8-e8c2f255e56f
# ╠═f5b8ddb4-d1e5-41ff-9c6f-19b9afc8128b
# ╠═01aeef0d-12e1-4e32-974f-f5405a661d95
# ╠═18a0121f-b0a5-46e5-8178-9bc800476b65
# ╠═2c71142a-34b2-4804-9af7-dc8147ab30a3
# ╠═3aa1788a-131a-402f-9493-32d353020219
# ╠═8740f867-109f-47c1-9e90-57ece3b83886
# ╠═8c52eac3-388f-49ea-bcc4-6d9a42085c52
# ╠═cec73687-de04-46e6-a73e-a4d1d3beafc5
# ╠═6fdaa0d2-b6af-40ca-b8e2-5cb8176bbadb
# ╠═edabab28-72bb-43f2-b90b-7b6ba939f2eb
# ╠═d8caa384-6cc8-4a96-93ee-3504a95c5214
# ╠═b8af0d87-f9e1-4ea7-9fa0-1c675233a8d3
# ╠═02464496-4c60-4e5f-a3af-0e5c2bb0b937
# ╠═f41cf5c7-d208-40f9-a4ee-5e992f73ca36
# ╠═5dc25eb7-7841-4c04-aa5e-7d63acc5a312
# ╠═b017b44e-9860-4f21-9b99-5a0cdf1b265c
# ╠═fca14add-ed54-449b-82a3-bffd03f88cdd
# ╠═4480a6a1-3299-43b9-b871-0b9f525577e5
# ╠═454cab5c-e715-46bd-a799-0e9e2e9a7f67
# ╠═4c1c878c-fcfb-4069-a27a-5b1288896955
