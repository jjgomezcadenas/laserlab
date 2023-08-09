### A Pluto.jl notebook ###
# v0.19.25

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
	using Statistics
	using StatsBase
	using Unitful 
	using UnitfulEquivalences 
	using PhysicalConstants
end

# ╔═╡ 04e7c8ad-d693-49b0-aa30-cbef6aa35b9b
using EasyFit

# ╔═╡ 74a2ea89-0b20-4d71-9461-82755b6ed6d8
#using PythonStructs

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
samples=[ 
"ANN155", "ANN155.1",           
"ANN155_10e-6",    
"ANN155_10e-6_after_degas",    
"ANN155_10e-6_before_degas"]

# ╔═╡ ca0b6daa-ff0c-499a-a64f-fe6dd5f93d74
samplesba=[ 
"ANN155+Ba", 
"ANN155+Ba.1",            
"ANN155_10e-6_Ba_after_degas",  
"ANN155_10e-6_Ba_before_degas",     
"ANN155_Ba_10e-6"]

# ╔═╡ ab612ed4-56f1-41b1-afc1-581f386fd0a7
md""" Select sample : $(@bind sample Select(samples))"""

# ╔═╡ 45e75283-429f-440a-8569-56e79cabccda
md""" Select sample Ba : $(@bind sampleba Select(samplesba))"""

# ╔═╡ e3a9b4c9-778e-4c64-bd83-a9c9d32cee8e


# ╔═╡ 711e562b-a17a-4a20-a03c-236118991ea3
begin
	rxdir ="/Users/jjgomezcadenas/BoldLab/BOLD/APD/"
	rdir = string(rxdir, sample)
	rdirba = string(rxdir, sampleba)
end

# ╔═╡ 0738a3bd-ddcc-4e90-b519-432bf66d6161
#rpath = "/Users/jjgomezcadenas/BoldLab/BOLD/APD/ANN155_10e-6_after_degas"
#	rpathba = "/Users/jjgomezcadenas/BoldLab/BOLD/APD/ANN155_10e-6_Ba_after_degas/"

# ╔═╡ 14192a57-c137-4a61-9074-98d6a0211ef9
dftpar=Dict(
	"Filter0"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=2*10^4, dymax=6*10^4, i0=[1,20, 1], l0=[90,30,30], dcut=100.0),
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
	"Filter0"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^4, dymax=2*10^5, i0=[1,20, 1], l0=[20,100,100], dcut=100.0),
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
	"Filter0"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=3*10^3, dymax=3*10^3, i0=[1,20, 1], l0=[20,100,100], dcut=100.0),)

# ╔═╡ f2c68ebb-fc5f-4dcf-a3b0-d327b7511dbd
begin
	zxdat =dftparacn["Filter0"]
	h1dtnAcn  = lfi.LaserLab.h1d(PhotonDFAcn.dtime *ttrnsacn, zxdat.nbins, zxdat.dcut, zxdat.dtmax)
	h1tsAcn  = lfi.LaserLab.h1d(PhotonDFAcn.timeTag * tnsacn*1e-9, zxdat.nbins, 0.0, zdat.tmax)
	pacndt =lfi.LaserLab.plot_h1dx(h1dtnAcn; xs = "dt (ns)", ys = "a.u", 
                           xlim=true, xl = (0.0, zxdat.dtmax),
						   ylim=true, yl=(1,zxdat.dymax), ylog=false,
						   markersize=2, fwl=true)
	pacntt= lfi.LaserLab.plot_h1dx(h1tsAcn; xs = "t (s)", ys = "a.u", 
                           xlim=true, xl = (0.0, zxdat.tmax),
						   ylim=true, yl=(1,zxdat.ymax), markersize=2, fwl=true)
	plot(pacndt, pacntt)
end

# ╔═╡ 7f195bc8-9ed0-4f80-8310-06a11a6c0e0e
begin
avgbkg = mean(h1tsAcn.weights)
md"""
- Average response for ACN = $(round(avgbkg)) pe
"""
end

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

# ╔═╡ c15de41a-15c5-4034-a0d4-77888b7784b5
begin
	vdataf = vdataf0 .- tsgn 
	vdataba = vdataba0 .- tsgn # notice I substract the tail of the free!
	p3dtfx = scatter(tdata, vdataf, markersize=2, label="ANN155")
	p3dtfbax = scatter(tdata, vdataba, markersize=2, label="ANN155+Ba")
	p3dtfba = scatter(tdata, vdataba, markersize=2, label="ANN155+Ba")
	p3dtf = scatter!(tdata, vdataf, markersize=2, label="ANN155")
end

# ╔═╡ 2cddd392-870c-4c8b-8f0c-a7d0c44502cb
begin
	eps = 1.0
	src=1:5:100 
	yy = [sum(vdataba[i:end]) + eps for i in src]
	xx = [sum(vdataf[i:end]) + eps for i in src]
	stnba = abs.(yy)./sqrt.(abs.(xx))
	ster =sqrt.(stnba./abs.(xx))
	#stnba= abs.([(sum(vdataba[i:end]) + eps) / (sum(vdataf[i:end]) + eps) for i in src])
	
	xt =[tdata[i] for i in src]
	
	#yer = sqrt.([stnba[i:end]^3 ./(sum(vdataf[i:end]) + eps) for i in src])
	prf2ba = scatter(xt, stnba, yerr=ster, markersize=2, yaxis=(:log, (10^3, 10^5)), label="Ratio Ba/free")

	yy2 = [sum(vdataba[i:end]) + eps for i in src]
	ytot =sum(vdataba)
	effba = abs.(yy2/ytot)
	sterx =sqrt.(yy2./ytot^2)
	
	peffba = scatter(xt, effba, yerr=sterx, markersize=2, yaxis=:log, label="EffBa")

	plot(prf2ba, peffba)
	
end

# ╔═╡ f4e75f3d-fcb1-4cb1-a7f9-e0e79ddd232a
stnba

# ╔═╡ 55772fb1-6356-470a-bcf0-8196acb090ed
effba

# ╔═╡ b5479f38-58fd-44e0-9bf8-7914e141bece
function fsgntn(csgnf, csgnba, csgnba2, xt=1000)
	expf(t) = csgnf[1] * exp(-t/csgnf[2])
	function expb(t) 
		if t <= xt
			return csgnba[1] * exp(-t/csgnba[2]) 
		else
			return csgnba2[1] * exp(-t/csgnba2[2]) 
		end
	end
	expf, expb
			
end

# ╔═╡ 67e9922a-0d15-436f-9e22-509ee65b1eb4
function qint(f::Function, t0::Number, t1::Number)
	return quadgk(f, t0, t1)[1]
end

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

# ╔═╡ 5c2c75d9-123e-46dc-a70b-126f7e881589
dtcut = zdat.dcut # in dt units

# ╔═╡ 363898af-1ee7-4fa9-94db-393a777b0255
#pdf = filter(:dtime  => x -> x * ttrns >= dtcut,  PhotonDF)

# ╔═╡ adef1149-cbf5-4403-88b0-e51196030cec
md"""
## Functions
"""

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

# ╔═╡ dff3ed15-5130-4e76-a0a6-6bf686e2ce16
function cut_scan(df, zdat, tns, srange=100:500:5100)
	SGN = zeros(length(collect(srange)))
	i=0
	for dtcut in srange
		i+=1
		pdf = filter(:dtime  => x -> x * ttrns >= dtcut,  df)
		h2ts  = lfi.LaserLab.h1d(pdf.timeTag * tns*1e-9, zdat.nbins, 0.0, zdat.tmax)
		SGN[i] = mean(h2ts.weights)
	end
	SGN
end

# ╔═╡ 9186b77c-f41c-4d80-98b8-74490caa8869
function sgntn(xfilter, file, zdat, tns, srange=100:500:5100)
	rpath = "/Users/jjgomezcadenas/BoldLab/BOLD/APD/ANN155"
	rpathba = "/Users/jjgomezcadenas/BoldLab/BOLD/APD/ANN155+Ba"
	zpath = joinpath(rpath,xfilter, file)
	zpathba = joinpath(rpathba,xfilter, file)
	
	PhotonDF, tagDict =  lfi.LaserLab.readHH(zpath, 100, true, false)
	PhotonDFBa, tagDictBa =  lfi.LaserLab.readHH(zpathba, 100, true, false)

	eff= cut_scan(PhotonDF, zdat, tns, srange)
	effba= cut_scan(PhotonDFBa, zdat, tns, srange)
	eff, effba
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

# ╔═╡ 6e481aaa-8fe3-4c37-a16b-1936e81b3576
csgnba, stdsgnba, ysgnba = tfit_sgn(tdata, vdataba; cbkg = [0.0,0.0], i0=badat.i0[1], l0=badat.l0[1], pa0=[100.0, 20.5])

# ╔═╡ 45355282-c022-405a-aa02-8aa6900ec24c
csgnba2, stdsgnba2, ysgnba2 = tfit_sgn(tdata, vdataba; cbkg = [0.0,0.0], i0=badat.i0[2], l0=badat.l0[2], pa0=[100.0, 50.5])

# ╔═╡ 200deb4a-c9da-4460-8f8f-99ed4ddc9820
begin
	zzpf = plot(p3dtfx, tdata, ysgnf, yaxis=:log, ylim=(1, 10^5))
	zzpba = plot(p3dtfbax, tdata, ysgnba, yaxis=:log, ylim=(1, 10^5))
	zzpba2 = plot(zzpba, tdata, ysgnba2, yaxis=:log, ylim=(1, 10^5))
	plot(zzpf, zzpba2)
end

# ╔═╡ 4ba0ff3e-36c9-4878-8deb-4cabb7f27c39
expf, expb = fsgntn(csgnf, csgnba, csgnba2)

# ╔═╡ 428e28e0-583c-4226-9999-f256144feb62
begin
	xxt = 1:10:5000
	yyb = expb.(xxt) 
	yyf = expf.(xxt)
	plot(xxt,yyb, lw=2, label="ANN155+Ba")
	plot!(xxt,yyf, lw=2, label="ANN155")
end

# ╔═╡ c05365a6-f993-43f5-a2c3-b5c64ebeed01
begin
	y2b = [qint(expb, ti, ti+10) for ti in 1:10:5000]
	y2f = [qint(expf, ti, ti+10) for ti in 1:10:5000]
	rsn = y2b ./ sqrt.(y2f)
	plot(xxt,rsn, lw=2, label="S/N", ylim=(0., 10^4))
end

# ╔═╡ b5d2505d-d677-49d5-a8fb-6c96ab3d025f
function tfit_sgn3(xdata, ydata, csgn; i0, l0, pa0=[100.0, 0.5])
	
	ffun(t, λ1, λ2) = csgn[1] * exp(-t/λ1)  +  csgn2[2] * exp(-t/λ2)
	
	mffun(t, p) = csgn1[1] * exp.(-t/p[1]) + csgn2[2] * exp.(-t/p[2]) 

	tfit(xdata[i0:l0], ydata[i0:l0], pa0, ffun, mffun)
end

# ╔═╡ a349f303-3005-4c97-809b-3dde8af05ccb
function fsgn(rdir, sample, filters, wfnm, zdat, file="default.ptu")
	dtcut = zdat.dcut
	
	dsgn=zeros(10)
	dsgnflt=zeros(10)
	edsgn=zeros(10)
	edsgnflt=zeros(10)
	zsgn = 0.0
	zsgnflt = 0.0

	sdir = string(rdir, sample)
	xfilters = readdir(sdir)
	filters = filter(x -> x != ".DS_Store", xfilters)
	i=0
	for xfilter in filters
		
		println("filter =", xfilter, " index = ", i)
		sdflt = joinpath(sdir, xfilter)
		zpath = joinpath(sdflt,file)
		
		PhotonDF, tagDict =  lfi.LaserLab.readHH(zpath, 100, true, false, false)
		h1ts  = lfi.LaserLab.h1d(PhotonDF.timeTag * tns*1e-9, zdat.nbins, 0.0, 
			                     zdat.tmax)
		pdf = filter(:dtime  => x -> x * ttrns >= dtcut,  PhotonDF)
		h2ts  = lfi.LaserLab.h1d(pdf.timeTag * tns*1e-9, zdat.nbins, 0.0, zdat.tmax)
		
		if xfilter != "Filter0"
			i+=1
			dsgn[i]  = sum(h1ts.weights) /(zdat.tmax * wfnm[1]) 
			dsgnflt[i] = sum(h2ts.weights) /(zdat.tmax * wfnm[1])
			edsgn[i]  = sqrt(sum(h1ts.weights))/(zdat.tmax * wfnm[1]) 
			edsgnflt[i] = sqrt(sum(h2ts.weights))/(zdat.tmax * wfnm[1]) 
		else
			zsgn = sum(h1ts.weights) /(zdat.tmax * wfnm[1]) 
			zsgnflt = sum(h2ts.weights) /(zdat.tmax * wfnm[1])
		end
	end
	(sgntot = zsgn, sgnflttot = zsgnflt, sgn = dsgn, esgn = edsgn, sgnflt=dsgnflt,
	 esgnflt=edsgnflt)
end

# ╔═╡ 0883e1f1-4a94-4762-be8a-c97a99f2627f
function csgn(rdir, samples, zdat, file="default.ptu")
	dtcut = zdat.dcut

	println(samples)
	dsgn=zeros(length(samples))
	dsgnflt=zeros(length(samples))
	edsgn=zeros(length(samples))
	edsgnflt=zeros(length(samples))
	
	i=0
	for sample in samples
		i+=1
		println(sample)
		sdir = string(rdir, sample)
		sdflt = joinpath(sdir, "Filter0")
		zpath = joinpath(sdflt,file)

		println("i = ", i, " zpath = ", zpath)
		
		PhotonDF, tagDict =  lfi.LaserLab.readHH(zpath, 100, true, false, false)
		h1ts  = lfi.LaserLab.h1d(PhotonDF.timeTag * tns*1e-9, zdat.nbins, 0.0, 
			                     zdat.tmax)
		pdf = filter(:dtime  => x -> x * ttrns >= dtcut,  PhotonDF)
		h2ts  = lfi.LaserLab.h1d(pdf.timeTag * tns*1e-9, zdat.nbins, 0.0, zdat.tmax)
		
		dsgn[i]  = sum(h1ts.weights) /(zdat.tmax) 
		dsgnflt[i] = sum(h2ts.weights) /(zdat.tmax )
		edsgn[i]  = sqrt(sum(h1ts.weights))/(zdat.tmax) 
		edsgnflt[i] = sqrt(sum(h2ts.weights))/(zdat.tmax ) 
		
	end
	(sgn = dsgn, esgn = edsgn, sgnflt=dsgnflt,
	 esgnflt=edsgnflt)
end

# ╔═╡ 0cf824bc-506b-4b37-881c-97b051735653
function tfit_sgn2(xdata, ydata, csgn1, csgn2; i0, l0, pa0=[100.0, 0.5])
	
	sgn(t) = csgn[1] * exp(-t/csgn[2])
	ffun(t, N1, N2) = N1 * exp(-t/csgn1[2])  +  N2 * exp(-t/csgn2[2])
	
	mffun(t, p) = p[1] * exp.(-t/csgn1[2]) + p[2] * exp.(-t/csgn2[2]) 

	tfit(xdata[i0:l0], ydata[i0:l0], pa0, ffun, mffun)
end

# ╔═╡ abeeb9dd-70cc-4d6e-88e6-4f71623b4d90
csgnba3, stdsgnba3, ysgnba3 = tfit_sgn2(tdata, vdataba, csgnba, csgnba2; i0=badat.i0[3], l0=badat.l0[3], pa0=[100.0, 50.5])

# ╔═╡ 0096d481-0f1c-4e61-95c1-3d90184de56a
md"""
### Fit results
- Free: exp1:  λ = $(round(csgnf[2])) +- $(round(stdsgnf[2])), N = $(round(csgnf[1])) +- $(round(stdsgnf[1])) in range [ $(tdata[fdat.i0[1]]), $(tdata[fdat.l0[1]])] ns

- Ba22+ : exp1: λ1 = $(round(csgnba[2])) +- $(round(stdsgnba[2])), N1 = $(round(csgnba[1])) +- $(round(stdsgnba[1])) in range [ $(tdata[badat.i0[1]]), $(tdata[badat.l0[1]])] ns

- Ba2+ : exp2: λ2 = $(round(csgnba2[2])) +- $(round(stdsgnba2[2])), N2 = $(round(csgnba2[1])) +- $(round(stdsgnba2[1])) in range [ $(tdata[badat.i0[2]]), $(tdata[badat.l0[2]])] ns

- normalizations: N1 = $(round(csgnba[1]/(csgnba[1]+csgnba2[1]), digits=2)), N2 = $(round(csgnba2[1]/(csgnba[1]+csgnba2[1]), digits=2))

- intensity double fit: N1 = $(round(csgnba3[1])), N2 = $(round(csgnba3[2]))
- normalizations: N1 = $(round(csgnba3[1]/(csgnba3[1]+csgnba3[2]), digits=2)), N2 = $(round(csgnba3[2]/(csgnba3[1]+csgnba3[2]), digits=2))
"""

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
# ╠═f2c68ebb-fc5f-4dcf-a3b0-d327b7511dbd
# ╠═7f195bc8-9ed0-4f80-8310-06a11a6c0e0e
# ╠═391fe36a-79b2-49f4-a3f1-c3f4ed83dea6
# ╠═ca0b6daa-ff0c-499a-a64f-fe6dd5f93d74
# ╠═ab612ed4-56f1-41b1-afc1-581f386fd0a7
# ╠═45e75283-429f-440a-8569-56e79cabccda
# ╠═e3a9b4c9-778e-4c64-bd83-a9c9d32cee8e
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
# ╠═b2d6bb5e-96f4-4fae-8717-7026c652b839
# ╠═1cedaad8-c2ae-4fba-8326-967d9e37ea9b
# ╠═5f143f8b-b72f-4c21-8bfc-17cf061653c9
# ╠═bcd73cbc-96c4-4cbf-9bc6-9c7bbc56a018
# ╠═c15de41a-15c5-4034-a0d4-77888b7784b5
# ╠═2cddd392-870c-4c8b-8f0c-a7d0c44502cb
# ╠═f4e75f3d-fcb1-4cb1-a7f9-e0e79ddd232a
# ╠═55772fb1-6356-470a-bcf0-8196acb090ed
# ╠═c0448a2e-c467-4cd6-aa44-82b488d5cf06
# ╠═6e481aaa-8fe3-4c37-a16b-1936e81b3576
# ╠═45355282-c022-405a-aa02-8aa6900ec24c
# ╠═abeeb9dd-70cc-4d6e-88e6-4f71623b4d90
# ╠═200deb4a-c9da-4460-8f8f-99ed4ddc9820
# ╠═b5479f38-58fd-44e0-9bf8-7914e141bece
# ╠═67e9922a-0d15-436f-9e22-509ee65b1eb4
# ╠═4ba0ff3e-36c9-4878-8deb-4cabb7f27c39
# ╠═428e28e0-583c-4226-9999-f256144feb62
# ╠═c05365a6-f993-43f5-a2c3-b5c64ebeed01
# ╠═0096d481-0f1c-4e61-95c1-3d90184de56a
# ╠═1f524d20-d336-4ae5-be00-781c057c1f9c
# ╠═5c2c75d9-123e-46dc-a70b-126f7e881589
# ╠═363898af-1ee7-4fa9-94db-393a777b0255
# ╠═adef1149-cbf5-4403-88b0-e51196030cec
# ╠═702f4805-8e1f-4fb8-bada-e95777055cdf
# ╠═9186b77c-f41c-4d80-98b8-74490caa8869
# ╠═dff3ed15-5130-4e76-a0a6-6bf686e2ce16
# ╠═3204a561-8244-466e-9dbe-e1302ef363d9
# ╠═37a64b1f-93f6-4402-a32a-be1207352d86
# ╠═235d81bf-5166-4276-b61c-1e7586f4d8bf
# ╠═0cf824bc-506b-4b37-881c-97b051735653
# ╠═b5d2505d-d677-49d5-a8fb-6c96ab3d025f
# ╠═a349f303-3005-4c97-809b-3dde8af05ccb
# ╠═0883e1f1-4a94-4762-be8a-c97a99f2627f
