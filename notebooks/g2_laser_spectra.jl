### A Pluto.jl notebook ###
# v0.19.2

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
PlutoUI.TableOfContents(title="G2 Laser Analysis", indent=true)

# ╔═╡ 92a0f8ec-9a1a-414b-8d8f-b96d496164a3
md"""
# Data
"""

# ╔═╡ af6ab2c4-5f03-41ae-9672-033b3916ae1e
md"""
Filter width(s) in nm
"""

# ╔═╡ aa4edf06-e3d0-44ad-92ce-49cec750b817
wfnm = [10.0, 24.0, 30.0, 40.0, 49.0, 52.0, 60.0, 40.0, 68.0, 10.0]

# ╔═╡ 595a19db-ef4e-415f-a29e-ee9ce055b55a
md"""
Filter central value in nm
"""

# ╔═╡ 044027cd-f816-4a4a-abfc-0da8682b270a
xfnm=[420.0,438.0,465.0,503.0,550.0,600.0,650.0,692.0,732.0,810.0]

# ╔═╡ b86e6ae7-678b-4564-a97b-375c6414c861
xfnm .- 0.5*wfnm

# ╔═╡ e29a8698-6676-42a3-b8e2-b1f1f72ed650
xfnm .+ 0.5*wfnm

# ╔═╡ 335d9287-90ff-4398-8b26-d9e41bd1e3c4
load("../notebooks/img/filters.png")

# ╔═╡ 13ecb88f-c1d3-46bb-b1ce-0a5f170fd3de
md"""
Exposure time in seconds
"""

# ╔═╡ 0413846f-c988-4663-b37f-c7fde224d8a5
texps = 10.0 

# ╔═╡ 04564c50-c59a-4fe7-8f70-e702dd0b8701
md"""
Laser power in μW
"""

# ╔═╡ e77e7680-74d8-4f4d-ac8c-e55498d1a9cd
pmuw = 134.0

# ╔═╡ 9da70d55-28b9-4c12-ba6a-c560832a5194
md"""
# Select working directory
"""

# ╔═╡ 6502fdee-cd60-4d32-ac77-c42c845b12ba
begin
	proot ="/Users/jj/JuliaProjects/LaserLab/labdata"
	sroot = "/Users/jj/JuliaProjects/LaserLab/data/G2Sl"
	dmtype = ["Imag","Dark"]
	xbold  = ["G2_BOLD_073_A2"]
	pfiles ="*.csv"
	
	md""" Select series : $(@bind wb Select(xbold))"""
	
end

# ╔═╡ 22e6d658-098d-480d-bbb1-d3831eada9c6
begin
	xmeas = joinpath(proot,wb)
	xnmeas = lfi.LaserLab.getbolddirs(xmeas)
end

# ╔═╡ 8adcb337-f6ee-423e-bf7f-0391f943ba32
md""" Select measurement : $(@bind wm Select(xnmeas))"""

# ╔═╡ ba4d56d7-70ee-4032-bc0c-7445f27571a9
begin
	xp = joinpath(xmeas,wm)
	xnp = lfi.LaserLab.getbolddirs(xp)
end

# ╔═╡ a641da7f-1a23-4cdb-b747-1344d1ff0179
md""" Select point : $(@bind wp Select(xnp))"""

# ╔═╡ 325dedd2-7079-4a95-9286-b91306cec8c8
xpath = joinpath(xp,wp);

# ╔═╡ 4f47264c-2cb4-45f1-93ad-793ca74da08e
md"""
working directory = $xpath
"""

# ╔═╡ 21c86a2b-ff51-4431-91f6-d6185e5f00ae
md"""
# Select Image 
"""

# ╔═╡ 8c9fc82f-280d-48a5-8a7b-9af061c4952b
begin
	xfiles = Glob.glob(pfiles, xpath)
	nxfiles = [split(f,"/")[end] for f in xfiles]
	rn = lfi.LaserLab.findpattern(nxfiles, "rep")
	md""" Select repetition number: $(@bind wr Select(rn))"""
end

# ╔═╡ 6cbfbd75-1885-406a-81ba-f99b94f78824
begin
	xfn = lfi.LaserLab.findpattern(nxfiles, "Filter")
	md""" Select filter: $(@bind wfn Select(xfn))"""
end

# ╔═╡ e0a87f8e-f9f4-4863-a591-1947b3ebab94
md""" Select Img/dark: $(@bind wid Select(dmtype))""" 

# ╔═╡ e62942fa-27ff-483b-bc73-b47d7074a365
begin
	if wid == "Dark"
		xfd = findpattern(nxfiles, "Dark.csv", "_", -1)
		md""" Select dark measurement: $(@bind wfd Select(xfd))"""
	else
		wfd = " "
	end
	wbs = string(wb)
	wms = string(wm)
	wps = string(wp)
	wrs = string(wr)
	wfs = string(wfn)
	wids = string(wid)
	wfds = string(wfd)
end

# ╔═╡ 39c39437-1810-4c80-9576-e6522ef26442
begin
	md"""
	- Series = $wbs
	- Measurement = $wms
	- Point = $wps
	- Filter = $wfs
	- Repetition = $wrs
	- Img/Dark = $wids
	- Dark Before/after (nothing if Imag) = $wfds
	"""
end

# ╔═╡ 945b1c82-dc26-4966-902b-c6d4042f69ac
nxfiles;

# ╔═╡ 5c226cb8-0dca-4ebd-bb35-1d51639a7eb0
xfiles;

# ╔═╡ 1a59c1a9-1905-4afe-af42-c85515c421a7
img = lfi.LaserLab.select_image(xfiles, wfs, wrs, wids, wfds);

# ╔═╡ 2ab028ee-1f9f-4493-a3fe-59602a493716
begin
imgmn2 = img ./maximum(img)
gimg2 = Gray.(imgmn2)
end

# ╔═╡ 934e9039-a774-4d10-bd6d-3e793c471d84
md"""
# Select corrected image from filter and repetition
"""

# ╔═╡ 591e3e72-402f-491a-b8d5-5c02a92b1d1d
cimg = lfi.LaserLab.get_corrected_image(xfiles, wfs, wrs);

# ╔═╡ a05cea87-9d3f-4f84-82ba-fdf38b247d46
begin
	imgmn = cimg.cimage ./maximum(cimg.cimage)
	gimg = Gray.(imgmn);
	imgmxn = cimg.image ./maximum(cimg.image)
	gimgx = Gray.(imgmxn)
	mosaicview(gimg, gimgx; nrow = 1)
end

# ╔═╡ 9c264015-b26c-42e8-8e97-e091c4eab8d1
imgmxnt = lfi.LaserLab.signal_around_maximum(cimg.cimage, cimg.cdark; nsigma=3);

# ╔═╡ b4237359-1f1a-4ae1-932d-5d9d7e14bba4
begin
	imgmxntn = imgmxnt.img ./maximum(imgmxnt.img)
	gimgnt = Gray.(imgmxntn)
end

# ╔═╡ 805196f0-2a11-4409-a18d-7620ad0f0190
begin
	stot = sum(imgmxnt.img)
md"""
- Size of the image =$(size(imgmxnt.img))
- Total sum around the peak = $(round(stot, sigdigits=3))
"""
end

# ╔═╡ b776cef6-90d6-4897-bd6e-b81a3a58b416
md"""
- Filter = $wfs
- repetition = $wrs
- **signal**
- value of max = $(imgmxnt.max)
- position of max: i = $(imgmxnt.imax)  j = $(imgmxnt.jmax)
- Total light in spot (DC subtracted) = $(round(stot, sigdigits=3))
"""

# ╔═╡ 01aeef0d-12e1-4e32-974f-f5405a661d95
md"""
# Reconstruct spectra for a given repetition
"""

# ╔═╡ 9463ffd7-ab4f-4b54-a750-a76c670ee9bb
md""" Check to compute sepctrum for this rep: $(@bind zrec CheckBox())"""

# ╔═╡ e12869d8-7607-4573-9b7c-499f2ba069a7
if zrec
	xfi = sort([parse(Int64, string(x)) for x in xfn])
	xfs = string.(xfi)
	ZMX = Vector{Float64}()
	ZSM = Vector{Float64}()
	for fltr in xfs
		cimgz = lfi.LaserLab.get_corrected_image(xfiles, fltr, wrs)
		imgmz = lfi.LaserLab.signal_around_maximum(cimgz.cimage, cimgz.cdark; nsigma=3)
		#println("filter =",fltr, " max val=", imgmz.max, " max i = ", imgmz.imax, " max j = ",imgmz.jmax, " sum = ", sum(imgmz.img))
		push!(ZMX,imgmz.max)
		push!(ZSM,sum(imgmz.img))
	end
end

# ╔═╡ cc85ccd8-692d-400e-b67a-de7164d7ec8f
md""" Check to save sepctrum for this rep: $(@bind zwrite CheckBox())"""

# ╔═╡ b64a4570-f8f0-4370-89ce-1fa1b21d089f
if zwrite && zrec
	sdf = DataFrame("filters" => xfnm, "sum"=> ZSM, "max"=> ZMX)
	sdfnm = string(wbs,"_",wms,"_",wps,"_rep_",wrs,".csv")
	sdff = joinpath(sroot, sdfnm)
	CSV.write(sdff, sdf)

	md"""
	- File to write = $sdfnm
	- Full path = $sdff
	"""
end

# ╔═╡ 634ed7f2-b95b-46f2-b5fd-f681751da423
md""" Check to plot reps and avg for point: $(@bind zrpavg CheckBox())"""

# ╔═╡ 1f76ac8b-3ef7-4dbc-bf84-03ec1453b961
md""" Check to save average for point: $(@bind zwavg CheckBox())"""

# ╔═╡ ab5ca517-7a95-41be-9887-48bc94dd793e
md""" Check to plot All points for selected measurement: $(@bind zpoint CheckBox())"""

# ╔═╡ e8fe5321-56c6-4fca-8fdd-ed7bff5ff428
if zpoint
	sdfmpy = string(wbs,"_",wms,"_avg.png")
	pngny = joinpath(sroot, sdfmpy)
	println(pngny)
	png(pngny)
end

# ╔═╡ 8740f867-109f-47c1-9e90-57ece3b83886
md"""
# Results
"""

# ╔═╡ c59ee355-abf4-488e-8e95-59e21299baf2
md""" Select point : $(@bind sp Select(xnp))"""

# ╔═╡ bdbc198a-8b46-4054-a714-af1e95113e80
md""" Check to plot G2/G2Ba for selected point: $(@bind zcomp CheckBox())"""

# ╔═╡ 599f682e-976c-49cf-b00a-3099136d80fd
md""" Check to load solution fluo spectra: $(@bind zfluo CheckBox())"""

# ╔═╡ 46783dcf-d452-4ae3-a8eb-efe2bfe9c1bd
if zfluo

	dffluo = lfi.LaserLab.load_df_from_csv(sroot, "Fluo_ACN_emi_340nm.csv", lfi.LaserLab.spG)
	plot(dffluo.W, dffluo[!,"5E_5"], lw=2, label="G2 solution")
	xlabel!("wavelength (nm)")
	ylabel!("counts")
end

# ╔═╡ cec73687-de04-46e6-a73e-a4d1d3beafc5
md"# Functions"

# ╔═╡ 454cab5c-e715-46bd-a799-0e9e2e9a7f67
function img_title(label::String, wbs::String, wms::String, wps::String, wrs::String)
	string(label, "Series: ", wbs, "Meas: ", wms, " Point: ", wps, " Rep: ", wrs)
end

# ╔═╡ f7d8c603-92ea-49d1-81bb-538945fac0b3
if zrec
	
	tit = img_title("", wbs, wms, wps, wrs)
	psm = plot(xfnm, ZSM, lw=2, label=string("rep ",wrs), title=tit, titlefontsize=10)
	ssm = scatter!(xfnm, ZSM, legend=false)
	xlabel!("wavelength (nm)")
	ylabel!("counts")
	#pmx = plot(xfnm, ZMX, lw=2, label=string("rep ",wrs))
	#xlabel!("wavelength (nm)")
	#ylabel!("counts")
	#plot(psm, pmx)
end

# ╔═╡ 2b12326f-166a-4439-9e46-239bcd7e55cf
if zrpavg
	sdfnmr1 = string(wbs,"_",wms,"_",wps,"_rep_1.csv")
	avgdf = lfi.LaserLab.load_df_from_csv(sroot, sdfnmr1, lfi.LaserLab.enG)
	plot(avgdf.filters, avgdf.sum, lw=2, label="rep 1")
	#scatter!(avgdf.filters, avgdf.sum, legend=false)
	
	for rep in rn[2:end]
		sdfnmr = string(wbs,"_",wms,"_",wps,"_rep_", rep,".csv")
		dfx = lfi.LaserLab.load_df_from_csv(sroot, sdfnmr, lfi.LaserLab.enG)
		
		plot!(dfx.filters, dfx.sum, lw=2, label=string("rep ", rep))
		#scatter!(dfx.filters, dfx.sum, legend=false)
		
		avgdf .+= dfx
	end
	avgdf ./= length(rn) 
	tita = img_title("", wbs, wms, wps,  "avg")
	plot!(xfnm, avgdf.sum, lw=2, label="average reps ", title=tita, titlefontsize=10)
	scatter!(xfnm, avgdf.sum, label="average reps ")
	xlabel!("wavelength (nm)")
	ylabel!("counts")
end

# ╔═╡ be5a950b-47b2-402c-8d24-57a0826e53e8
if zwavg
	sdfav = DataFrame("filters" => xfnm, "sum"=> avgdf.sum, "max"=> avgdf.max)
	sdfnmav = string(wbs,"_",wms,"_",wps,"_avg.csv")
	sdffav = joinpath(sroot, sdfnmav)
	CSV.write(sdffav, sdfav)

	md"""
	- File to write = $sdfnmav
	- Full path = $sdffav
	"""
end

# ╔═╡ 84a44685-1651-47e9-bfdf-2b21a537fa93
if zpoint
	titb = img_title("", wbs, wms, "all",  "avg")
	for xp in xnp
		sdfnpt = string(wbs,"_",wms,"_",xp,"_avg.csv")
		dfpnt = lfi.LaserLab.load_df_from_csv(sroot, sdfnpt, lfi.LaserLab.enG)
		if xp=="Point1"
			plot(dfpnt.filters, dfpnt.sum, lw=2, label=xp, title=titb, titlefontsize=10)
		else
			plot!(dfpnt.filters, dfpnt.sum, lw=2, label=xp)
		end
		scatter!(dfpnt.filters, dfpnt.sum, label=xp)
		
	end
	xlabel!("wavelength (nm)")
	ylabel!("counts")
end

# ╔═╡ 41cf396d-01ba-4e09-9ac9-bdf0a713de9e
if zcomp
	xxp = string(sp)
	
	sdfmpt = string(wbs,"_",xnmeas[1],"_",xxp,"_avg.csv")
	dfmpt = lfi.LaserLab.load_df_from_csv(sroot, sdfmpt, lfi.LaserLab.enG)
	titd = img_title("", wbs, " G2/G2Ba ", xxp,  "avg")
	
	pxm = plot(dfmpt.filters, dfmpt.sum, lw=2, label=xnmeas[1], title=titd, titlefontsize=10)
	scatter!(dfmpt.filters, dfmpt.sum, label=xnmeas[1])

	sdfmpt2 = string(wbs,"_",xnmeas[2],"_",xxp,"_avg.csv")
	dfmpt2 = lfi.LaserLab.load_df_from_csv(sroot, sdfmpt2, lfi.LaserLab.enG)
	
	plot!(dfmpt2.filters, dfmpt2.sum, lw=2, label=xnmeas[2])
	scatter!(dfmpt2.filters, dfmpt2.sum, label=xnmeas[2])
	xlabel!("wavelength (nm)")
	ylabel!("counts")
	
end

# ╔═╡ 7eaaa7a9-2a5b-4ed6-b3b8-1a03e0fb8081
if zcomp
	sdfmpx = string(wbs,"_","G2-G2Ba","_",xxp,"_avg.png")
	pngn = joinpath(sroot, sdfmpx)
	println(pngn)
	png(pngn)
end

# ╔═╡ Cell order:
# ╠═8a67d331-35ae-44d6-942c-a61db4e7afa4
# ╠═65beba5e-776d-40d5-96a6-d0a2f3a9a4ea
# ╠═d1b57350-bd9d-11ec-21f4-3d7a8656a31c
# ╟─53943a33-9cf3-4c9f-add5-278fdc3bc0fe
# ╠═b829dbf6-486e-4bd1-914c-54fe6d389b9c
# ╠═7aafa82d-4fe3-4bdc-a14c-ec48a7178260
# ╟─92a0f8ec-9a1a-414b-8d8f-b96d496164a3
# ╟─af6ab2c4-5f03-41ae-9672-033b3916ae1e
# ╠═aa4edf06-e3d0-44ad-92ce-49cec750b817
# ╟─595a19db-ef4e-415f-a29e-ee9ce055b55a
# ╠═044027cd-f816-4a4a-abfc-0da8682b270a
# ╠═b86e6ae7-678b-4564-a97b-375c6414c861
# ╠═e29a8698-6676-42a3-b8e2-b1f1f72ed650
# ╠═335d9287-90ff-4398-8b26-d9e41bd1e3c4
# ╟─13ecb88f-c1d3-46bb-b1ce-0a5f170fd3de
# ╠═0413846f-c988-4663-b37f-c7fde224d8a5
# ╟─04564c50-c59a-4fe7-8f70-e702dd0b8701
# ╠═e77e7680-74d8-4f4d-ac8c-e55498d1a9cd
# ╟─9da70d55-28b9-4c12-ba6a-c560832a5194
# ╟─6502fdee-cd60-4d32-ac77-c42c845b12ba
# ╟─22e6d658-098d-480d-bbb1-d3831eada9c6
# ╟─8adcb337-f6ee-423e-bf7f-0391f943ba32
# ╟─ba4d56d7-70ee-4032-bc0c-7445f27571a9
# ╟─a641da7f-1a23-4cdb-b747-1344d1ff0179
# ╟─325dedd2-7079-4a95-9286-b91306cec8c8
# ╟─4f47264c-2cb4-45f1-93ad-793ca74da08e
# ╟─21c86a2b-ff51-4431-91f6-d6185e5f00ae
# ╟─8c9fc82f-280d-48a5-8a7b-9af061c4952b
# ╠═6cbfbd75-1885-406a-81ba-f99b94f78824
# ╟─e0a87f8e-f9f4-4863-a591-1947b3ebab94
# ╟─e62942fa-27ff-483b-bc73-b47d7074a365
# ╟─39c39437-1810-4c80-9576-e6522ef26442
# ╟─945b1c82-dc26-4966-902b-c6d4042f69ac
# ╟─5c226cb8-0dca-4ebd-bb35-1d51639a7eb0
# ╟─1a59c1a9-1905-4afe-af42-c85515c421a7
# ╟─2ab028ee-1f9f-4493-a3fe-59602a493716
# ╟─934e9039-a774-4d10-bd6d-3e793c471d84
# ╟─591e3e72-402f-491a-b8d5-5c02a92b1d1d
# ╟─a05cea87-9d3f-4f84-82ba-fdf38b247d46
# ╟─9c264015-b26c-42e8-8e97-e091c4eab8d1
# ╟─b4237359-1f1a-4ae1-932d-5d9d7e14bba4
# ╟─805196f0-2a11-4409-a18d-7620ad0f0190
# ╟─b776cef6-90d6-4897-bd6e-b81a3a58b416
# ╟─01aeef0d-12e1-4e32-974f-f5405a661d95
# ╟─9463ffd7-ab4f-4b54-a750-a76c670ee9bb
# ╟─e12869d8-7607-4573-9b7c-499f2ba069a7
# ╟─f7d8c603-92ea-49d1-81bb-538945fac0b3
# ╟─cc85ccd8-692d-400e-b67a-de7164d7ec8f
# ╟─b64a4570-f8f0-4370-89ce-1fa1b21d089f
# ╟─634ed7f2-b95b-46f2-b5fd-f681751da423
# ╟─2b12326f-166a-4439-9e46-239bcd7e55cf
# ╟─1f76ac8b-3ef7-4dbc-bf84-03ec1453b961
# ╟─be5a950b-47b2-402c-8d24-57a0826e53e8
# ╟─ab5ca517-7a95-41be-9887-48bc94dd793e
# ╟─84a44685-1651-47e9-bfdf-2b21a537fa93
# ╟─e8fe5321-56c6-4fca-8fdd-ed7bff5ff428
# ╟─8740f867-109f-47c1-9e90-57ece3b83886
# ╟─c59ee355-abf4-488e-8e95-59e21299baf2
# ╟─bdbc198a-8b46-4054-a714-af1e95113e80
# ╟─41cf396d-01ba-4e09-9ac9-bdf0a713de9e
# ╟─7eaaa7a9-2a5b-4ed6-b3b8-1a03e0fb8081
# ╟─599f682e-976c-49cf-b00a-3099136d80fd
# ╟─46783dcf-d452-4ae3-a8eb-efe2bfe9c1bd
# ╠═cec73687-de04-46e6-a73e-a4d1d3beafc5
# ╠═454cab5c-e715-46bd-a799-0e9e2e9a7f67
