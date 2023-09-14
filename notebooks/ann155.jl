### A Pluto.jl notebook ###
# v0.19.27

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

# ╔═╡ 1459f36b-afa4-405e-a9e8-e6ae32ea752e
using Pkg; Pkg.activate(ENV["JLaserLab"])

# ╔═╡ e6550506-f9a3-4999-9fef-4e0cf70cc1ec
begin
	using PlutoUI
	using CSV
	using DataFrames
	using QuadGK
	using Plots
	using Printf
	using Markdown
	using InteractiveUtils
	using LsqFit
	using EasyFit
	using Statistics
	using StatsBase
	using Unitful 
	using UnitfulEquivalences 
	using PhysicalConstants
	using LaTeXStrings
end

# ╔═╡ 74a2ea89-0b20-4d71-9461-82755b6ed6d8
#using PythonStructs

# ╔═╡ 04e7c8ad-d693-49b0-aa30-cbef6aa35b9b
#using EasyFit

# ╔═╡ 9657e6f5-b3b9-448f-8bb5-28d5a64c4c37
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

# ╔═╡ ea788fbe-160c-40fd-b967-d49bc5403398
lfi = ingredients("../src/LaserLab.jl")

# ╔═╡ ce389036-74fa-4715-88fe-914004f4a5c3
md"""
## Filters
"""

# ╔═╡ 138d8b2f-c148-4155-863f-ca400259a458
begin
	xfnm = [400.0, 438.0, 503.0, 549.0, 575.0, 600.0, 630.0, 676.0, 732.0, 810.0]
	wfnm = [40.0, 24.0,   40.0,   17.0,  15.0,  14.0,  38.0,  29.0,  68.0,  10.0]
end

# ╔═╡ 8e7c68c7-ff2e-4895-aaf7-189f1d9eb39d
md"""

## Analysis: ANN155
"""

# ╔═╡ 8438aa92-3c20-417f-8733-757f7e07d37f
md"""
### Read ACN data (background data)
"""

# ╔═╡ 81549a07-5eda-466d-bafb-24f530041d8c
begin
	acnpath = "/Users/jjgomezcadenas/BoldLab/BOLD/APD/ACN/Filter0/default.ptu"
	PhotonDFAcn, tagDictAcn =  lfi.LaserLab.readHH(acnpath, 100, true, false)
	ttrnsacn = tagDictAcn["MeasDesc_Resolution"]* 1e+9
	tnsacn   = tagDictAcn["MeasDesc_GlobalResolution"] * 1e+9

	md"""
	#### Time bins:
	- dt  time bin(units of ns for time difference) = $ttrnsacn
	- tag time bin(units of ns for time ) = $tnsacn
	"""
end


# ╔═╡ 391fe36a-79b2-49f4-a3f1-c3f4ed83dea6
samples=[ "ZFFL51_ANN265.2",  
	"ANN155_10e-5_after_degas_remeasured",
"ANN155", "ANN155.1",  "ANN155_10e-5_after_degas",         
"ANN155_10e-6",    
"ANN155_10e-6_after_degas",    
"ANN155_10e-6_before_degas", "ANN155_10e-6_after_degas.2",
"ANN155_A_Quartz",
"ANN155_A_Quartz_box_before_degas",
"ANN155_A_Quartz_box_after_degas",
"ANN155_B_Quartz_box_after_degas"
]

# ╔═╡ ca0b6daa-ff0c-499a-a64f-fe6dd5f93d74
samplesba=[ "ZFFL51_ANN265_Ba.2",
	"ANN155_10e-5_Ba_after_degas_remeasured",
"ANN155+Ba", 
"ANN155+Ba.1",       "ANN155_10e-5_Ba_after_degas",     
"ANN155_10e-6_Ba_after_degas",  
"ANN155_10e-6_Ba_before_degas",     
"ANN155_Ba_10e-6", "ANN155_10e-6_Ba_after_degas.2",
"ANN155_A_Ba_Quartz",
"ANN155_A_Ba_Quartz_box_before_degas",
"ANN155_A_Ba_Quartz_box_after_degas",
"ANN155_B_Ba_Quartz_box_after_degas"]

# ╔═╡ ab612ed4-56f1-41b1-afc1-581f386fd0a7
md""" Select sample : $(@bind xsample Select(samples))"""

# ╔═╡ 45e75283-429f-440a-8569-56e79cabccda
md""" Select sample Ba : $(@bind sampleba Select(samplesba))"""

# ╔═╡ 711e562b-a17a-4a20-a03c-236118991ea3
begin
	rxdir ="/Users/jjgomezcadenas/BoldLab/BOLD/APD/"
	rdir = string(rxdir, xsample)
	rdirba = string(rxdir, sampleba)
end

# ╔═╡ 0738a3bd-ddcc-4e90-b519-432bf66d6161
#rpath = "/Users/jjgomezcadenas/BoldLab/BOLD/APD/ANN155_10e-6_after_degas"
#	rpathba = "/Users/jjgomezcadenas/BoldLab/BOLD/APD/ANN155_10e-6_Ba_after_degas/"

# ╔═╡ 14192a57-c137-4a61-9074-98d6a0211ef9
dftpar=Dict(
	"Filter0"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^4, dymax=6*10^4, i0=[1,10, 1], l0=[20,100,100], dcut=100.0),
	"Filter1"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^3, dymax=2*10^3, i0=[5,10, 15], l0=[7,15,100], dcut=100.0),
	"Filter2"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^3, dymax=2*10^3, i0=[5,10, 15], l0=[7,15,100], dcut=100.0),
	"Filter3"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^3, dymax=2*10^3, i0=[1,10, 15], l0=[7,100,100], dcut=100.0),
	"Filter4"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^3, dymax=3*10^3, i0=[1,15, 15], l0=[15,100,100], dcut=100.0),
	"Filter5"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^3, dymax=5*10^3, i0=[1,15, 15], l0=[15,30,100], dcut=100.0),
	"Filter6"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^3, dymax=5*10^3, i0=[1, 11, 20], l0=[10,30,100], dcut=100.0),
	"Filter7"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^3, dymax=1*10^4, i0=[1,15, 20], l0=[15,30,100], dcut=100.0),
	"Filter8"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^3, dymax=1*10^4, i0=[1,15, 20], l0=[15,30,100], dcut=100.0),
	"Filter9"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^3, dymax=6*10^3, i0=[1,15, 20], l0=[15,30,100], dcut=100.0),
	"Filter10"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^3, dymax=6*10^3, i0=[5,10, 15], l0=[7,15,100], dcut=150.0),
)

# ╔═╡ 79876fb7-7177-4277-acc6-71a8a5195868
dftparba=Dict(
	"Filter0"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^4, dymax=2*10^5, i0=[1,10, 1], l0=[20,100,100], dcut=100.0),
	"Filter1"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^3, dymax=5*10^3, i0=[5,10, 15], l0=[7,15,100], dcut=100.0),
	"Filter2"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^3, dymax=5*10^3, i0=[5,10, 15], l0=[7,15,100], dcut=100.0),
	"Filter3"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^4, dymax=5*10^3, i0=[1,10, 15], l0=[7,100,100], dcut=100.0),
	"Filter4"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^3, dymax=1*10^4, i0=[1,15, 15], l0=[15,100,100], dcut=100.0),
	"Filter5"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^3, dymax=1.5*10^4, i0=[1,15, 15], l0=[15,30,100], dcut=100.0),
	"Filter6"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=1*10^4, dymax=1.5*10^4, i0=[1, 20, 20], l0=[30,100,100], dcut=100.0),
	"Filter7"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=1*10^4, dymax=3*10^4, i0=[1,10, 10], l0=[40,100,100], dcut=100.0),
	"Filter8"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=2*10^4, dymax=2*10^4, i0=[1,20, 20], l0=[20,30,100], dcut=100.0),
	"Filter9"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=2*10^4, dymax=2*10^4, i0=[1,15, 15], l0=[15,100,100], dcut=100.0),
	"Filter10"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^3, dymax=6*10^3, i0=[5,10, 15], l0=[7,15,100], dcut=150.0),
)

# ╔═╡ 8a57aa29-a9d9-47e5-b13b-34a62533e54d
dftparacn=Dict(
	"Filter0"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=3*10^3, dymax=3*10^3, i0=[1,20, 1], l0=[20,100,100], dcut=100.0))

# ╔═╡ c792b017-9c33-4d03-853a-9672f8eaf481
begin
	xfilters = readdir(rdir)
	filters = sort(filter(x -> x != ".DS_Store" && x != "Info.txt", xfilters))
	md""" Select filter for free molecule : $(@bind xfilter Select(filters))"""
end

# ╔═╡ 1284ccc9-ea32-4fa7-8350-5ba96fbdacf3
begin
	xfiltersba = readdir(rdirba)
	filtersba = sort(filter(x -> x != ".DS_Store" && x != "info.txt", xfiltersba))
	md""" Select filter for chelated molecule : $(@bind xfilterba Select(filters))"""
end

# ╔═╡ 220f711c-957f-4174-b3c7-f4de4d910328
begin
	fdat = dftpar[xfilter]
	badat = dftparba[xfilter]
end

# ╔═╡ f08c25a5-3b2b-4b43-a021-035ba9539e46
begin
	sdflt = joinpath(rdir, xfilter)
	ff = readdir(sdflt)
	xfiles = readdir(sdflt)
	files = filter(x -> x != ".DS_Store", xfiles)
	md""" Select file : $(@bind file Select(files))"""
end

# ╔═╡ 7f2ed729-c09b-4b63-98aa-16b85200911b
begin
	sdfltba = joinpath(rdirba, xfilterba)
	ffba = readdir(sdfltba)
	xfilesba = readdir(sdfltba)
	filesba = filter(x -> x != ".DS_Store", xfilesba)
	md""" Select file Ba : $(@bind fileba Select(filesba))"""
end

# ╔═╡ c07044d5-7bde-42c2-b104-7035487c41e2
flts = ["Filter0","Filter3","Filter4","Filter5","Filter6","Filter7",
	"Filter8","Filter9"]

# ╔═╡ 69acf984-964b-4589-9b65-c15450c5d1a0
#NEVT, NEVTBA = getnorm(rpath,file, flts)

# ╔═╡ 155a4254-bb97-4761-a6bd-9e0fe1c04422
#begin
#	efff = NEVT[2:end] ./ NEVT[1]
#	effb = NEVTBA[2:end] ./ NEVTBA[1]
#	peff = scatter(flts[2:end], efff, label="ANN155")
#	pefba = scatter(flts[2:end], effb, label="ANN155+BA")
#	plot(peff, pefba, layout = (1, 2), size=(850,550))
#end

# ╔═╡ ffd49f31-6a3d-4800-b257-0f1841c50706
begin
	zpath = joinpath(rdir, xfilter, file)
	zpathba = joinpath(rdirba, xfilterba, fileba)

	md"""
	- path to free ANN =$zpath
	- path to ANN+ Ba =$zpathba
	"""
end

# ╔═╡ 843c6950-d629-495b-9d9e-06d920d80229
begin
	PhotonDF, tagDict =  lfi.LaserLab.readHH(zpath, 100, true, false)
	PhotonDFBa, tagDictBa =  lfi.LaserLab.readHH(zpathba, 100, true, false)
end

# ╔═╡ 11754d33-3316-4163-86c0-0ae25a7d679b
fDict = Dict("fctime"  =>tagDict["File_CreatingTime"],
            "acqtime"  =>tagDict["MeasDesc_AcquisitionTime"]*1e-3,
            "nrecords" =>tagDict["TTResult_NumberOfRecords"],
            "tbin"    =>tagDict["MeasDesc_GlobalResolution"] * 1e+9,
            "dtbin"  =>tagDict["MeasDesc_Resolution"]* 1e+9)

# ╔═╡ a888a973-bb71-4577-9024-b5b788891482
fDictBa = Dict("fctime"  =>tagDictBa["File_CreatingTime"],
            "acqtime"  =>tagDictBa["MeasDesc_AcquisitionTime"]*1e-3,
            "nrecords" =>tagDictBa["TTResult_NumberOfRecords"],
            "tbin"    =>tagDictBa["MeasDesc_GlobalResolution"] * 1e+9,
            "dtbin"  =>tagDictBa["MeasDesc_Resolution"]* 1e+9)

# ╔═╡ cd2fb65b-7e9c-4a3f-bc26-51dc87b1470f
begin
	ttrns = tagDict["MeasDesc_Resolution"]* 1e+9
	tns   = tagDict["MeasDesc_GlobalResolution"] * 1e+9
	ttrnsba = tagDictBa["MeasDesc_Resolution"]* 1e+9
	tnsba   = tagDictBa["MeasDesc_GlobalResolution"] * 1e+9
end

# ╔═╡ 0caf7df2-a7e3-4369-a77c-e923d6bd0772
begin
	pdtf = histogram(PhotonDF.dtime * ttrns, label="ANN155")
	pdtba = histogram(PhotonDFBa.dtime * ttrnsba,  label="ANN155+BA")
	plot(pdtf, pdtba)
end

# ╔═╡ a60b9a00-f386-4e28-bbc8-e68c01c82137
begin
	h1dtns0  = lfi.LaserLab.h1d(PhotonDF.dtime *ttrns, fdat.nbins, 50.0, fdat.dtmax)
	h1dtnsba0  = lfi.LaserLab.h1d(PhotonDFBa.dtime *ttrnsba, badat.nbins, 50.0, 	
			                      badat.dtmax)
	pdtfx1 =lfi.LaserLab.plot_h1dx(h1dtns0; xs = "dt (ns)", ys = "a.u", 
                           xlim=true, xl = (0.0, fdat.dtmax), label="C1",
						   ylim=true, yl=(1e+1,2e+6), ylog=true, legend=true,
						   markersize=1, fwl=false)
	ppdtbax1 =lfi.LaserLab.plot_h1dx(h1dtnsba0; xs = "dt (ns)", ys = "a.u", 
                           xlim=true, xl = (1.0, badat.dtmax), label="C1+Ba",
						   ylim=true, yl=(1e+1,2e+6), ylog=true, legend=true,
						   markersize=1, fwl=false)
	plot(pdtfx1, ppdtbax1, size=(600,400))
end

# ╔═╡ b2d6bb5e-96f4-4fae-8717-7026c652b839
begin
	pttf = histogram(PhotonDF.timeTag * tns*1e-9, label="ANN155")
	pttba = histogram(PhotonDFBa.timeTag * tnsba*1e-9,  label="ANN155+BA")
	plot(pttf, pttba)
end

# ╔═╡ 1cedaad8-c2ae-4fba-8326-967d9e37ea9b
begin
h1dtns  = lfi.LaserLab.h1d(PhotonDF.dtime *ttrns, fdat.nbins, fdat.dcut, fdat.dtmax)
h1dtnsba  = lfi.LaserLab.h1d(PhotonDFBa.dtime *ttrnsba, badat.nbins, badat.dcut, badat.dtmax)
end

# ╔═╡ 5f143f8b-b72f-4c21-8bfc-17cf061653c9
begin
	ppdtf =lfi.LaserLab.plot_h1dx(h1dtns; xs = "dt (ns)", ys = "a.u", 
                           xlim=true, xl = (0.0, fdat.dtmax), label="ANN155",
						   ylim=false, yl=(1,fdat.dymax), ylog=false, legend=true,
						   markersize=2, fwl=true)
	ppdtba =lfi.LaserLab.plot_h1dx(h1dtnsba; xs = "dt (ns)", ys = "a.u", 
                           xlim=true, xl = (0.0, badat.dtmax), label="ANN155+Ba",
						   ylim=false, yl=(1,badat.dymax), ylog=false, legend=true,
						   markersize=2, fwl=true)
	plot(ppdtf, ppdtba)
end

# ╔═╡ bcd73cbc-96c4-4cbf-9bc6-9c7bbc56a018
begin
	tdata = h1dtns.centers
	vdataf0 = h1dtns.weights
	nb = fdat.nbins
	tsgn = mean(vdataf0[nb-10:nb])
	
	vdataba0 = h1dtnsba.weights
	nb = badat.nbins
	tsgnba = mean(vdataba0[nb-10:nb])
	md"""
	- tail signal (ANN155) = $(round(tsgn))
	- tail signal (ANN155+Ba) = $(round(tsgnba))
	"""
end

# ╔═╡ 3111871e-03cd-4519-a408-a009a0e1a776
length(tdata)

# ╔═╡ c15de41a-15c5-4034-a0d4-77888b7784b5
begin
	vdataf = vdataf0 .- tsgn 
	vdataba = vdataba0 .- tsgn # notice I substract the tail of the free!
	p3dtfx = scatter(tdata, vdataf, markersize=1, xlabel="t (ns)", ylabel="a.u.",
		yaxis = (:identity, (1.0, 10000)), label="C1")
	p3dtfbax = scatter(tdata, vdataba, markersize=1, xlabel="t (ns)", ylabel="a.u.",
		yaxis = (:identity, (1.0, 3e+4)), label="C1+Ba")
	p3dtfba = scatter(tdata, vdataba, markersize=1, yaxis = (:identity, (1.0, 10^5)), label="ANN155+Ba")
	p3dtf = scatter!(tdata, vdataf, markersize=1, yaxis = (:identity, (1.0, 10^5)),
		xlabel="t (ns)", ylabel="a.u.", label="ANN155")
end

# ╔═╡ f71f9286-0759-4315-8534-5bf2ce66ded6
plot(p3dtfx)

# ╔═╡ 76a7d6fb-e5f6-4ad8-a17b-bfc548d9fdc1
plot(p3dtfbax)

# ╔═╡ 2cddd392-870c-4c8b-8f0c-a7d0c44502cb
begin
	eps = 1.0
	src=1:5:100 
	yy = [sum(vdataba[i:end]) + eps for i in src]
	xx = [sum(vdataf[i:end]) + eps for i in src]
	stnba = abs.(yy)./sqrt.(abs.(xx))
	ster =sqrt.(stnba./abs.(xx))
	
	
	xt =[tdata[i] for i in src]
	
	prf2ba = scatter(xt, stnba, yerr=ster, markersize=2, yaxis=(:log), label="Ratio Ba/free")
	prf2bal = scatter(xt, stnba, yerr=ster, markersize=2, label="Ratio Ba/free")

	yy2 = [sum(vdataba[i:end]) + eps for i in src]
	ytot =sum(vdataba)
	sterx =sqrt.(abs.(yy2)./ytot^2)
	effba = abs.(yy2/ytot)
	
	peffba = scatter(xt, effba, yerr=sterx, markersize=2, yaxis=:log, label="EffBa")
	peffbal = scatter(xt, effba, yerr=sterx, markersize=2, label="EffBa")

	plot(prf2ba, peffba,prf2bal, peffbal)
	
end

# ╔═╡ 84a56ba2-dba2-49e6-918d-715573db8328
zefit = fitexp(tdata, vdataf,n=2);

# ╔═╡ abeeb9dd-70cc-4d6e-88e6-4f71623b4d90
xefit = fitexp(tdata, vdataba,n=2);

# ╔═╡ 1f524d20-d336-4ae5-be00-781c057c1f9c
begin
	#h1ts  = lfi.LaserLab.h1d(PhotonDF.timeTag * tns*1e-9, zdat.nbins, 0.0, zdat.tmax)
	h1tsf  = lfi.LaserLab.h1d(PhotonDF.timeTag * tns*1e-9, fdat.nbins, 0.0, fdat.tmax)
	h1tsba  = lfi.LaserLab.h1d(PhotonDFBa.timeTag * tns*1e-9, badat.nbins, 
		                       0.0, fdat.tmax)
	pztf = lfi.LaserLab.plot_h1dx(h1tsf; xs = "t (s)", ys = "a.u", 
                           xlim=true, xl = (0.0, fdat.tmax), label="ANN155",
	                       legend=true, 
						   ylim=true, yl=(1,fdat.ymax), markersize=2, fwl=true)
	pztba = lfi.LaserLab.plot_h1dx(h1tsba; xs = "t (s)", ys = "a.u", 
                           xlim=true, xl = (0.0, badat.tmax), label="ANN155+Ba",
	                       legend=true, 
						   ylim=true, yl=(1,badat.ymax), markersize=2, fwl=true)
	plot(pztf, pztba)
end

# ╔═╡ f14d4695-e7d4-4b4a-9393-9bd09aefe6e4
length(zefit.ypred)

# ╔═╡ 363898af-1ee7-4fa9-94db-393a777b0255
#pdf = filter(:dtime  => x -> x * ttrns >= dtcut,  PhotonDF)

# ╔═╡ adef1149-cbf5-4403-88b0-e51196030cec
md"""
## Functions
"""

# ╔═╡ fdf023c2-0186-48f0-acec-f2a14f218be6
function qint(f::Function, t0::Number, t1::Number)
	return quadgk(f, t0, t1)[1]
end

# ╔═╡ 1b038953-ca91-4511-aaea-a0d732be35bd
function fsgntn(zefit, xefit)
	expf(t) = zefit.a[1] * exp(-t/zefit.b[1]) 
	#+ zefit.a[2] * exp(-t/zefit.b[2]) + zefit.c
	expb(t) = xefit.a[1] * exp(-t/xefit.b[1]) + xefit.a[2] * exp(-t/xefit.b[2]) + xefit.c
	expf, expb
			
end

# ╔═╡ 5c2c75d9-123e-46dc-a70b-126f7e881589
begin
	expf, expb = fsgntn(zefit, xefit)
	t0 = 10.0
	tl = 5000.0
	ts = 1000.0
	xxt=t0:ts:tl
	y2b = [qint(expb, ti, ti+ts) for ti in xxt]
	y2f = [qint(expf, ti, ti+ts) for ti in xxt]
	xef = [qint(expb, ti, tl) / qint(expb, t0, tl) for ti in xxt]
	rsn = y2b ./ sqrt.(abs.(y2f))
	py1 = plot(xxt, y2b, lw=2, yaxis = (:log,(1e+3, 5e+6)), label="C1+Ba")
	py2 = plot(xxt, y2f, lw=2, yaxis = (:log,(1, 5e+6)), label="C1")
	py3 = plot(xxt,rsn, lw=2, yaxis = (:log,(1e+3, 1e+6)), label="S/N")
	py4 = plot(xxt,xef, lw=2, yaxis = (:log,(1e-3, 1.0)), label="eff")
	plot( py3, py4)
end

# ╔═╡ 86a3c22e-5675-4127-a61a-532f6bd52367
begin
	zxxt=0.0:10.0:5000.0
	zy2b = [qint(expb, t0, tl) for t0 in zxxt]
	zy2f = [qint(expf, t0, tl) for t0 in zxxt]
	zxef = [qint(expb, t0, tl) / qint(expb, 0.0, 5000.0) for t0 in zxxt]
	zrsn = zy2b ./ sqrt.(abs.(zy2f))
	zpy1 = plot(zxxt, zy2b, lw=2, yaxis = (:log,(1e+3, 5e+7)), label="C1+Ba")
	zpy2 = plot(zxxt, zy2f, lw=2, yaxis = (:log,(1, 5e+6)), label="C1")
	zpy3 = plot(zxxt,zrsn, lw=2, yaxis = (:log,(1e+3, 1e+6)), label="S/N")
	zpy4 = plot(zxxt,zxef, lw=2, yaxis = (:log,(1e-3, 1.0)), label="eff")
	plot( zpy1, zpy2, zpy3, zpy4)
	#plot(zpy1, zpy2, zpy4)
end

# ╔═╡ bf6c9f2c-b014-4fef-9e0c-1fcbc4983b4f
begin
	yy2b = [qint(expb, ti, ti+10) for ti in xxt]
	yy2f = [qint(expf, ti, ti+10) for ti in xxt]
	pyy1 = plot(xxt, expb.(xxt), lw=2, label="C1+Ba")
	pyy2 = plot(xxt, expf.(xxt), lw=2, label="C1")
	rrsn = [expb(xx) /sqrt(expf(xx)) for xx in xxt]
	pyy3 = plot(xxt,rrsn, lw=2, label="S/N")
	plot(pyy1, pyy2, pyy3)
end

# ╔═╡ 7a5fbaab-caca-4429-ad74-b50f9620b4ae
begin
	pex1 = plot(tdata, xefit.ypred, lw=2, label="C1+Ba")
	pex11 = plot(pex1,tdata, expb.(tdata), lw=2, label="C1+Ba")
	pex2 = plot(tdata, zefit.ypred, lw=2, label="C1")
	pex22 = plot(pex2,tdata, expf.(tdata), lw=2, label="C1+Ba")
	plot(pex11, pex22)
end

# ╔═╡ 702f4805-8e1f-4fb8-bada-e95777055cdf
function getnorm(rpath,file, flts)
	
	NEVT = zeros(length(flts))
	NEVTBA = zeros(length(flts))
	for (i,xflt) in enumerate(flts)
		zzpath = joinpath(rpath, xflt, file)
		zzpathba = joinpath(rpathba, xflt, file)
		NEVT[i]   =  lfi.LaserLab.readHH(zzpath, 100, true, false, true, true)
		NEVTBA[i] =  lfi.LaserLab.readHH(zzpathba, 100, true, false, true, true)
	end
	NEVT, NEVTBA
end

# ╔═╡ 3204a561-8244-466e-9dbe-e1302ef363d9
function tfit(xdata, ydata, pa0, ffun, mffun)
	
	fit = curve_fit(mffun, xdata, ydata, pa0)
	coef(fit), stderror(fit), ffun.(tdata, coef(fit)...)
end

# ╔═╡ 37a64b1f-93f6-4402-a32a-be1207352d86
function tfit_bkg(xdata, ydata; i0, l0, pa0=[100.0, 0.5])
	
	ffun(t, a,b) = a + b * t
	mffun(t, p) = p[1] .+ t .* p[2]

	tfit(xdata[i0:l0], ydata[i0:l0], pa0, ffun, mffun)
end

# ╔═╡ 235d81bf-5166-4276-b61c-1e7586f4d8bf
function tfit_sgn(xdata, ydata; cbkg, i0, l0, pa0=[100.0, 0.5])
	
	bkg(t) = cbkg[1] + cbkg[2] * t
	
	ffun(t, N1, λ) = N1 * exp(-t/λ)  +  bkg(t)
	mffun(t, p) = p[1] * exp.(-t/p[2]) +  bkg.(t)

	tfit(xdata[i0:l0], ydata[i0:l0], pa0, ffun, mffun)
end

# ╔═╡ c0448a2e-c467-4cd6-aa44-82b488d5cf06
csgnf, stdsgnf, ysgnf = tfit_sgn(tdata, vdataf; cbkg = [0.0,0.0], i0=fdat.i0[1], l0=fdat.l0[1], pa0=[100.0, 20.5])

# ╔═╡ 16580758-0bcd-410e-b26f-db655868ca4a
csgnf2, stdsgnf2, ysgnf2 = tfit_sgn(tdata, vdataf; cbkg = [0.0,0.0], i0=fdat.i0[2], l0=fdat.l0[2], pa0=[100.0, 50.5])

# ╔═╡ 6e481aaa-8fe3-4c37-a16b-1936e81b3576
csgnba, stdsgnba, ysgnba = tfit_sgn(tdata, vdataba; cbkg = [0.0,0.0], i0=badat.i0[1], l0=badat.l0[1], pa0=[100.0, 20.5])

# ╔═╡ 45355282-c022-405a-aa02-8aa6900ec24c
csgnba2, stdsgnba2, ysgnba2 = tfit_sgn(tdata, vdataba; cbkg = [0.0,0.0], i0=badat.i0[2], l0=badat.l0[2], pa0=[100.0, 50.5])

# ╔═╡ 200deb4a-c9da-4460-8f8f-99ed4ddc9820
begin
	zzpf = plot(p3dtfx, tdata, ysgnf, ylim=(1, 1*10^3), label="exp1")
	zzpf2 = plot(zzpf, tdata, ysgnf2, ylim=(1, 1*10^3), label="exp2")
	sexp1 = round(zefit.b[1], digits=2)
	sexpb1 = round(xefit.b[1], digits=2)
	sexpb2 = round(xefit.b[2], digits=2)
	lbl1 = L"\lambda = %$sexp1 \, \mu s"
	lbl2 = L"\lambda1 = %$sexpb1 \, \mu s, \lambda2 = %$sexpb2 \, \mu s"
	zzpf3 = plot(p3dtfx, tdata, zefit.ypred,  ylim=(1, 1*10^4), 
		     label=lbl1)
	
	zzpba = plot(p3dtfbax, tdata, ysgnba, yaxis=:identity, ylim=(1, 1*10^5), label="exp1")
	zzpba2 = plot(zzpba, tdata, ysgnba2, yaxis=:identity, ylim=(1, 1*10^5), label="exp2")
	zzpba3 = plot(p3dtfbax, tdata, xefit.ypred, yaxis=:identity, ylim=(1, 3*10^4),label=lbl2)
	plot(zzpf3, zzpba3)
end

# ╔═╡ 0096d481-0f1c-4e61-95c1-3d90184de56a
md"""
### Fit results
- Free: exp1:  λ = $(round(csgnf[2])) +- $(round(stdsgnf[2])), N = $(round(csgnf[1])) +- $(round(stdsgnf[1])) in range [ $(tdata[fdat.i0[1]]), $(tdata[fdat.l0[1]])] ns

- Free: exp2:  λ = $(round(csgnf2[2])) +- $(round(stdsgnf2[2])), N = $(round(csgnf2[1])) +- $(round(stdsgnf2[1])) in range [ $(tdata[fdat.i0[2]]), $(tdata[fdat.l0[2]])] ns

- Free : exp1+exp2: λ1 = $(round(zefit.b[1])),  λ1 = $(round(zefit.b[2])), N1 = $(round(zefit.a[1])) N2 = $(round(zefit.a[2])) 

- Free:  normalizations (independent fits): N1 = $(round(csgnf[1]/(csgnf[1]+csgnf[1]), digits=2)), N2 = $(round(csgnf2[1]/(csgnf2[1]+csgnf2[1]), digits=2))

- Free: normalizations (combined fit): N1 = $(round(zefit.a[1]/(zefit.a[1]+zefit.a[2]), digits=2)), N2 = $(round(zefit.a[2]/(zefit.a[1]+zefit.a[2]), digits=2))

- Ba2+ : exp1: λ1 = $(round(csgnba[2])) +- $(round(stdsgnba[2])), N1 = $(round(csgnba[1])) +- $(round(stdsgnba[1])) in range [ $(tdata[badat.i0[1]]), $(tdata[badat.l0[1]])] ns

- Ba2+ : exp2: λ2 = $(round(csgnba2[2])) +- $(round(stdsgnba2[2])), N2 = $(round(csgnba2[1])) +- $(round(stdsgnba2[1])) in range [ $(tdata[badat.i0[2]]), $(tdata[badat.l0[2]])] ns

- Ba2+ : exp1+exp2: λ1 = $(round(xefit.b[1])),  λ1 = $(round(xefit.b[2])), N1 = $(round(xefit.a[1])) N2 = $(round(xefit.a[2])) 

- Ba2+:  normalizations (independent fits): N1 = $(round(csgnba[1]/(csgnba[1]+csgnba2[1]), digits=2)), N2 = $(round(csgnba2[1]/(csgnba[1]+csgnba2[1]), digits=2))

- Ba2++: normalizations (combined fit): N1 = $(round(xefit.a[1]/(xefit.a[1]+xefit.a[2]), digits=2)), N2 = $(round(xefit.a[2]/(xefit.a[1]+xefit.a[2]), digits=2))
"""

# ╔═╡ 949d1557-d86f-4471-9b26-eb3bf215cf33
function tfit_exp2(xdata, ydata; i0, l0, pa0=[100.0, 0.5, 100.0, 0.5])
	
	
	ffun(t, N1, N2, λ1, λ2) = N1 * (exp(-t/λ1)  +  N2 * exp(-t/λ1))
	mffun(t, p) = p[1] * (exp.(-t/p[2]) +  p[3] * exp.(-t/p[4]))

	tfit(xdata[i0:l0], ydata[i0:l0], pa0, ffun, mffun)
end

# ╔═╡ 0cf824bc-506b-4b37-881c-97b051735653
function tfit_sgn2(xdata, ydata, csgn1, csgn2; i0, l0, pa0=[100.0, 0.5])
	
	sgn(t) = csgn[1] * exp(-t/csgn[2])
	ffun(t, N1, N2) = N1 * exp(-t/csgn1[2])  +  N2 * exp(-t/csgn2[2])
	
	mffun(t, p) = p[1] * exp.(-t/csgn1[2]) + p[2] * exp.(-t/csgn2[2]) 

	tfit(xdata[i0:l0], ydata[i0:l0], pa0, ffun, mffun)
end

# ╔═╡ b5d2505d-d677-49d5-a8fb-6c96ab3d025f
function tfit_sgn3(xdata, ydata, csgn; i0, l0, pa0=[100.0, 0.5])
	
	ffun(t, λ1, λ2) = csgn[1] * exp(-t/λ1)  +  csgn2[2] * exp(-t/λ2)
	
	mffun(t, p) = csgn1[1] * exp.(-t/p[1]) + csgn2[2] * exp.(-t/p[2]) 

	tfit(xdata[i0:l0], ydata[i0:l0], pa0, ffun, mffun)
end

# ╔═╡ Cell order:
# ╠═1459f36b-afa4-405e-a9e8-e6ae32ea752e
# ╠═e6550506-f9a3-4999-9fef-4e0cf70cc1ec
# ╠═74a2ea89-0b20-4d71-9461-82755b6ed6d8
# ╠═04e7c8ad-d693-49b0-aa30-cbef6aa35b9b
# ╠═9657e6f5-b3b9-448f-8bb5-28d5a64c4c37
# ╠═ea788fbe-160c-40fd-b967-d49bc5403398
# ╠═ce389036-74fa-4715-88fe-914004f4a5c3
# ╠═138d8b2f-c148-4155-863f-ca400259a458
# ╠═8e7c68c7-ff2e-4895-aaf7-189f1d9eb39d
# ╠═8438aa92-3c20-417f-8733-757f7e07d37f
# ╠═81549a07-5eda-466d-bafb-24f530041d8c
# ╠═391fe36a-79b2-49f4-a3f1-c3f4ed83dea6
# ╠═ca0b6daa-ff0c-499a-a64f-fe6dd5f93d74
# ╠═ab612ed4-56f1-41b1-afc1-581f386fd0a7
# ╠═45e75283-429f-440a-8569-56e79cabccda
# ╠═711e562b-a17a-4a20-a03c-236118991ea3
# ╠═0738a3bd-ddcc-4e90-b519-432bf66d6161
# ╠═14192a57-c137-4a61-9074-98d6a0211ef9
# ╠═79876fb7-7177-4277-acc6-71a8a5195868
# ╠═8a57aa29-a9d9-47e5-b13b-34a62533e54d
# ╠═c792b017-9c33-4d03-853a-9672f8eaf481
# ╠═1284ccc9-ea32-4fa7-8350-5ba96fbdacf3
# ╠═220f711c-957f-4174-b3c7-f4de4d910328
# ╠═f08c25a5-3b2b-4b43-a021-035ba9539e46
# ╠═7f2ed729-c09b-4b63-98aa-16b85200911b
# ╠═c07044d5-7bde-42c2-b104-7035487c41e2
# ╠═69acf984-964b-4589-9b65-c15450c5d1a0
# ╠═155a4254-bb97-4761-a6bd-9e0fe1c04422
# ╠═ffd49f31-6a3d-4800-b257-0f1841c50706
# ╠═843c6950-d629-495b-9d9e-06d920d80229
# ╠═11754d33-3316-4163-86c0-0ae25a7d679b
# ╠═a888a973-bb71-4577-9024-b5b788891482
# ╠═cd2fb65b-7e9c-4a3f-bc26-51dc87b1470f
# ╠═0caf7df2-a7e3-4369-a77c-e923d6bd0772
# ╠═a60b9a00-f386-4e28-bbc8-e68c01c82137
# ╠═b2d6bb5e-96f4-4fae-8717-7026c652b839
# ╠═1cedaad8-c2ae-4fba-8326-967d9e37ea9b
# ╠═5f143f8b-b72f-4c21-8bfc-17cf061653c9
# ╠═bcd73cbc-96c4-4cbf-9bc6-9c7bbc56a018
# ╠═3111871e-03cd-4519-a408-a009a0e1a776
# ╠═c15de41a-15c5-4034-a0d4-77888b7784b5
# ╠═f71f9286-0759-4315-8534-5bf2ce66ded6
# ╠═76a7d6fb-e5f6-4ad8-a17b-bfc548d9fdc1
# ╠═2cddd392-870c-4c8b-8f0c-a7d0c44502cb
# ╠═c0448a2e-c467-4cd6-aa44-82b488d5cf06
# ╠═16580758-0bcd-410e-b26f-db655868ca4a
# ╠═84a56ba2-dba2-49e6-918d-715573db8328
# ╠═6e481aaa-8fe3-4c37-a16b-1936e81b3576
# ╠═45355282-c022-405a-aa02-8aa6900ec24c
# ╠═abeeb9dd-70cc-4d6e-88e6-4f71623b4d90
# ╠═200deb4a-c9da-4460-8f8f-99ed4ddc9820
# ╠═0096d481-0f1c-4e61-95c1-3d90184de56a
# ╠═1f524d20-d336-4ae5-be00-781c057c1f9c
# ╠═f14d4695-e7d4-4b4a-9393-9bd09aefe6e4
# ╠═5c2c75d9-123e-46dc-a70b-126f7e881589
# ╠═86a3c22e-5675-4127-a61a-532f6bd52367
# ╠═bf6c9f2c-b014-4fef-9e0c-1fcbc4983b4f
# ╠═7a5fbaab-caca-4429-ad74-b50f9620b4ae
# ╠═363898af-1ee7-4fa9-94db-393a777b0255
# ╠═adef1149-cbf5-4403-88b0-e51196030cec
# ╠═fdf023c2-0186-48f0-acec-f2a14f218be6
# ╠═1b038953-ca91-4511-aaea-a0d732be35bd
# ╠═702f4805-8e1f-4fb8-bada-e95777055cdf
# ╠═3204a561-8244-466e-9dbe-e1302ef363d9
# ╠═37a64b1f-93f6-4402-a32a-be1207352d86
# ╠═235d81bf-5166-4276-b61c-1e7586f4d8bf
# ╠═949d1557-d86f-4471-9b26-eb3bf215cf33
# ╠═0cf824bc-506b-4b37-881c-97b051735653
# ╠═b5d2505d-d677-49d5-a8fb-6c96ab3d025f
