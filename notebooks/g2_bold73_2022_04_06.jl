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
PlutoUI.TableOfContents(title="G2 Bold73 analysis", indent=true)

# ╔═╡ 6926c7d6-0c1d-4054-9d02-3b12c21ad81c
collect(10:-1:1 )

# ╔═╡ 9da70d55-28b9-4c12-ba6a-c560832a5194
md"""
# Read Data
"""

# ╔═╡ 6502fdee-cd60-4d32-ac77-c42c845b12ba
begin
	proot ="/Users/jj/JuliaProjects/LaserLab/labdata"
	sroot = "/Users/jj/JuliaProjects/LaserLab/data/G2Sl"
	psample = "G2_BOLD_073_A1"
	xrun = ["G20406","G20407", "G2Ba0408"]
	prun = Dict("G20406"=>"G2SL_20220406", 
		        "G20407"=>"G2SL_20220407",
		        "G2Ba0408"=>"G2SL-BaCl2_20220408")
	pfiles ="*.dat"
	sflts = string.(collect(2:11));
	sreps = string.(collect(0:10));
	
	md""" Select run type : $(@bind whichxr Select(xrun))"""
end

# ╔═╡ 21743b1b-3a48-4eb2-9dcd-3a8262c45dec
begin
	path = joinpath(proot,psample, prun[whichxr])
	fls = joinpath(path,pfiles)
	nfiles = length(fls) 
	files = Glob.glob(pfiles, path);
	#names = [split(f,"/")[end] for f in files];
	
md"
- Path for glob: $path
- files : $fls
- number of files: $nfiles
"
end

# ╔═╡ a2ce79cf-a447-4cab-8843-a2cd691b4004
md"""
# Select Filter and Rep
"""

# ╔═╡ 90bda847-8afe-4968-9a02-50ad7c7cc572
md""" Select Filter number (2-11): $(@bind whichf Select(sflts))"""

# ╔═╡ 1cd87997-3ff3-40fa-98b7-ee9655848a76
md""" Select repetition number (0-10): $(@bind whichr Select(sreps))"""

# ╔═╡ 7d756675-333a-4f8d-aa4d-db84abe72183
whichr

# ╔═╡ 21f71c85-9e22-4834-a7d7-76e3d1ab681b
#imgm, drkm, sdrkm, simgm = get_images(files, whichf, whichr);

# ╔═╡ 0387a9f9-4e0b-4568-a698-5fe9e66e5c39
imgnt = lfi.LaserLab.get_images(files, whichf, whichr);

# ╔═╡ 5c283457-fc1d-4249-bd4e-029204e5c57b
lfi.LaserLab.plot_images(imgnt, whichf, whichr)

# ╔═╡ 60b420bd-9bfd-4361-bec9-11a1ccd205a5
begin
	imgmn = imgnt.cimage ./maximum(imgnt.cimage)
	gimg = Gray.(imgmn);
	imgmn_edge = lfi.LaserLab.sujoy(gimg, four_connectivity=true)
	img_edge = binarize(imgmn_edge, Otsu())
	img_edgem = convert(Array{Float64}, img_edge)
	mosaicview(gimg, img_edge; nrow = 1)
end

# ╔═╡ fa2c16c8-b9e6-43ee-96cd-d91e73d528ca
img_edgem

# ╔═╡ ef2fd803-4833-400f-9a40-61fff3d50942
imgmx = lfi.LaserLab.signal_around_maximum(imgnt.cimage, imgnt.cdark; nsigma=5)

# ╔═╡ 1049bc5f-3fcc-4574-8db2-3988227c5fb5
begin
	imgmxn = imgmx ./maximum(imgmx)
	gimgx = Gray.(imgmxn);
	mosaicview(gimg, gimgx; nrow = 1)
end

# ╔═╡ bd3beecc-2154-4c70-a8b4-7334b8cc643d
gimgx

# ╔═╡ c9a93f4d-cd32-410d-bc7e-61fa253bac4b

sedge, nedge = lfi.LaserLab.sum_edge(imgnt.cimage, img_edgem)

# ╔═╡ eb2f4ed6-1b6a-41e0-bc1b-1de01b2a2cd1
mean(imgmx)

# ╔═╡ 9a40236e-fb18-4e2e-849d-51f4292ec15e
std(imgmx)

# ╔═╡ 20902f12-03c0-4a98-937a-a9ff8433fb1b
sum(imgmx)

# ╔═╡ d2ff1ea4-7c02-4b0d-bd85-2fdc8c0495b6
std(imgnt.cdark)

# ╔═╡ dcdf6ccc-bbe4-4133-9211-5e7648736711
sum(imgnt.image)

# ╔═╡ c9d4401d-69a5-4df2-a645-76bf3e6bcb63
sum(imgnt.cimage)

# ╔═╡ 2644c058-591f-4963-a587-3ce58d1b017d


# ╔═╡ dfc2478c-9834-427d-a976-38068b9e59fc
#plot_images(imgm, drkm, sdrkm, simgm)

# ╔═╡ 4aab4565-280c-432c-8f6d-6774a4dbe8be
#begin
#	simgmn = simgm ./maximum(simgm);
#	gimg = Gray.(simgmn);
#	plot(gimg, layout=(1,1), titlefontsize=10)
#end

# ╔═╡ 17eac4d7-ea71-4205-87d1-e239317c2d7e
typeof(simgm)

# ╔═╡ 64ed6a8c-1bca-4c19-83ff-d507ba4359a5
sdrkmn = sdrkm ./maximum(sdrkm);

# ╔═╡ 1313a166-29f3-4f5e-8e79-dbfab22a69d9
begin
	gdrk = Gray.(sdrkmn)
	plot(gdrk, layout=(1,1), titlefontsize=10)
end

# ╔═╡ 01aeef0d-12e1-4e32-974f-f5405a661d95
md"""
# Reconstruct spectra
"""

# ╔═╡ 8feb58fd-3cb7-4860-bc22-52bcee93d40a
scs = [string(c) for c in collect(5:10)]

# ╔═╡ 18a0121f-b0a5-46e5-8178-9bc800476b65
@bind zrepx CheckBox()

# ╔═╡ 3aa1788a-131a-402f-9493-32d353020219
fdf = lfi.LaserLab.load_df_from_csv(sroot, "G2Bold073A1Fluo.csv", lfi.LaserLab.enG)

# ╔═╡ 8740f867-109f-47c1-9e90-57ece3b83886


# ╔═╡ 8c52eac3-388f-49ea-bcc4-6d9a42085c52
begin
	plot(collect(1:10), fdf.r1, lw=2, label="rep5 S")
	plot!(collect(1:10), fdf.r2, lw=2, label="rep6 S")
	plot!(collect(1:10), fdf.r3, lw=2, label="rep7 S")
	plot!(collect(1:10), fdf.r4, lw=2, label="rep8 S")
	plot!(collect(1:10), fdf.r5, lw=2, label="rep9 S")
	plot!(collect(1:10), fdf.r6, lw=2, label="rep10 S")
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

# ╔═╡ 6d08fd5c-6757-43dd-a0e3-bb161e78cd26
pimgx, stot = findpeaks(simgm; xrange=30);

# ╔═╡ f56df91d-53ce-4ca7-a9c8-7740cc3bc430
maxx, indxt = findmax(pimgx);

# ╔═╡ 3444eb77-b397-4ac3-99a2-972abc97090d
begin
	pimgn = pimgx ./maximum(pimgx);
	gimg2  =Gray.(pimgn);
	plot(gimg2,layout=(1,1), titlefontsize=10)
end

# ╔═╡ e6816855-d168-4317-9ddc-8e6b13c1d0a1
pimgd, stotd = findpeaks(sdrkm; xrange=30);

# ╔═╡ a9d90521-a7e5-4cd2-92ca-45de8c03e514
maxd, indxd = findmax(pimgd)

# ╔═╡ 55a34a0a-ae88-4b80-b6e1-515da142faa8
thrz = maxd + sqrt(maxd)

# ╔═╡ 87369d14-7318-417c-bc33-6a1bc807c0ee
begin
	pimgdn = pimgd ./maximum(pimgd)
	gimgd  =Gray.(pimgdn)
	plot(gimgd, layout=(1,1), titlefontsize=10)
end

# ╔═╡ b776cef6-90d6-4897-bd6e-b81a3a58b416
md"""
- Filter = $whichf
- repetition = $whichr
- **signal**
- value of max = $maxx
- position of max = $indxt
- Total light in spot (DC subtracted) = $(round(stot, sigdigits=3))
- **DC**
- value of max = $maxd
- position of max = $indxd
- DC in spot  = $(round(stotd, sigdigits=3))
"""

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
	rp = [repx[i]["stot"] for i in 1:6]
	
	fluodf = DataFrame("r1" => repx[1]["stot"],
		               "r2" => repx[2]["stot"],
		               "r3" => repx[3]["stot"],
		               "r4" => repx[4]["stot"],
		               "r5" => repx[5]["stot"],
		               "r6" => repx[6]["stot"])
	
	dfname = joinpath(sroot, "G2Bold073A1Fluo.csv")
	CSV.write(dfname, fluodf)

	
	plot(collect(1:10), repx[1]["stot"], lw=2, label="rep5 S")
	plot!(collect(1:10), repx[1]["dtot"], label="rep5 D")
	plot!(collect(1:10), repx[2]["stot"], lw=2, label="rep6 S")
	plot!(collect(1:10), repx[3]["stot"], lw=2, label="rep7 S")
	plot!(collect(1:10), repx[4]["stot"], lw=2, label="rep8 S")
	plot!(collect(1:10), repx[5]["stot"], lw=2, label="rep9 S")
	plot!(collect(1:10), repx[6]["stot"], lw=2, label="rep10 S")
end

# ╔═╡ 3794ff0b-87bf-4b08-9b33-4016280ebbbf
function plot_images2(imgnt)
	imghm = heatmap(imgnt.image, title=string("rep=", whichr, " filt=", whichf, " Img" ), titlefontsize=10)

	drkhm = heatmap(imgnt.dark, title=string("rep=", whichr, " filt=", whichf, " Dark" ), titlefontsize=10)
	
	simghm = heatmap(imgnt.corr, title=string("rep=", whichr, " filt=", whichf, "Img - Dark" ))
	
	plot(imghm, drkhm, simghm, layout=(3,1), titlefontsize=10)
end

# ╔═╡ fca14add-ed54-449b-82a3-bffd03f88cdd
function plot_images(imgm, drkm, sdrkm, simgm)
	imghm = heatmap(imgm, fill_z=imgm, title=string("rep=", whichr, " filt=", whichf, " Img" ), titlefontsize=10)

	drkhm = heatmap(drkm, title=string("rep=", whichr, " filt=", whichf, " Dark" ), titlefontsize=10)
	
	sdrkhm = heatmap(sdrkm, title=string("rep=", whichr, " filt=", whichf, "Dark - dark2" ), titlefontsize=10)

	simghm = heatmap(simgm, title=string("rep=", whichr, " filt=", whichf, "Img - Dark" ))
	
	plot(imghm, drkhm, sdrkhm, simghm, layout=(2,2), titlefontsize=10)
end

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
# ╠═6926c7d6-0c1d-4054-9d02-3b12c21ad81c
# ╠═9da70d55-28b9-4c12-ba6a-c560832a5194
# ╠═6502fdee-cd60-4d32-ac77-c42c845b12ba
# ╠═21743b1b-3a48-4eb2-9dcd-3a8262c45dec
# ╠═7d756675-333a-4f8d-aa4d-db84abe72183
# ╠═a2ce79cf-a447-4cab-8843-a2cd691b4004
# ╠═90bda847-8afe-4968-9a02-50ad7c7cc572
# ╠═1cd87997-3ff3-40fa-98b7-ee9655848a76
# ╠═21f71c85-9e22-4834-a7d7-76e3d1ab681b
# ╠═0387a9f9-4e0b-4568-a698-5fe9e66e5c39
# ╠═5c283457-fc1d-4249-bd4e-029204e5c57b
# ╠═60b420bd-9bfd-4361-bec9-11a1ccd205a5
# ╠═fa2c16c8-b9e6-43ee-96cd-d91e73d528ca
# ╠═ef2fd803-4833-400f-9a40-61fff3d50942
# ╠═1049bc5f-3fcc-4574-8db2-3988227c5fb5
# ╠═bd3beecc-2154-4c70-a8b4-7334b8cc643d
# ╠═c9a93f4d-cd32-410d-bc7e-61fa253bac4b
# ╠═eb2f4ed6-1b6a-41e0-bc1b-1de01b2a2cd1
# ╠═9a40236e-fb18-4e2e-849d-51f4292ec15e
# ╠═20902f12-03c0-4a98-937a-a9ff8433fb1b
# ╠═d2ff1ea4-7c02-4b0d-bd85-2fdc8c0495b6
# ╠═dcdf6ccc-bbe4-4133-9211-5e7648736711
# ╠═c9d4401d-69a5-4df2-a645-76bf3e6bcb63
# ╠═2644c058-591f-4963-a587-3ce58d1b017d
# ╠═dfc2478c-9834-427d-a976-38068b9e59fc
# ╠═4aab4565-280c-432c-8f6d-6774a4dbe8be
# ╠═6d08fd5c-6757-43dd-a0e3-bb161e78cd26
# ╠═17eac4d7-ea71-4205-87d1-e239317c2d7e
# ╠═f56df91d-53ce-4ca7-a9c8-7740cc3bc430
# ╠═3444eb77-b397-4ac3-99a2-972abc97090d
# ╠═64ed6a8c-1bca-4c19-83ff-d507ba4359a5
# ╠═1313a166-29f3-4f5e-8e79-dbfab22a69d9
# ╠═e6816855-d168-4317-9ddc-8e6b13c1d0a1
# ╠═a9d90521-a7e5-4cd2-92ca-45de8c03e514
# ╠═55a34a0a-ae88-4b80-b6e1-515da142faa8
# ╠═87369d14-7318-417c-bc33-6a1bc807c0ee
# ╠═b776cef6-90d6-4897-bd6e-b81a3a58b416
# ╠═01aeef0d-12e1-4e32-974f-f5405a661d95
# ╠═8feb58fd-3cb7-4860-bc22-52bcee93d40a
# ╠═18a0121f-b0a5-46e5-8178-9bc800476b65
# ╠═2c71142a-34b2-4804-9af7-dc8147ab30a3
# ╠═3aa1788a-131a-402f-9493-32d353020219
# ╠═8740f867-109f-47c1-9e90-57ece3b83886
# ╠═8c52eac3-388f-49ea-bcc4-6d9a42085c52
# ╠═cec73687-de04-46e6-a73e-a4d1d3beafc5
# ╠═6fdaa0d2-b6af-40ca-b8e2-5cb8176bbadb
# ╠═edabab28-72bb-43f2-b90b-7b6ba939f2eb
# ╠═b8af0d87-f9e1-4ea7-9fa0-1c675233a8d3
# ╠═f41cf5c7-d208-40f9-a4ee-5e992f73ca36
# ╠═5dc25eb7-7841-4c04-aa5e-7d63acc5a312
# ╠═3794ff0b-87bf-4b08-9b33-4016280ebbbf
# ╠═fca14add-ed54-449b-82a3-bffd03f88cdd
# ╠═4480a6a1-3299-43b9-b871-0b9f525577e5
# ╠═454cab5c-e715-46bd-a799-0e9e2e9a7f67
# ╠═4c1c878c-fcfb-4069-a27a-5b1288896955
