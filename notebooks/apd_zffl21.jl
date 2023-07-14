### A Pluto.jl notebook ###
# v0.19.26

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
	#using Dates
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

## Analysis: ZFFL21
"""

# ╔═╡ e3a9b4c9-778e-4c64-bd83-a9c9d32cee8e
rdir ="/Users/jjgomezcadenas/BoldLab/BOLD/APD/"

# ╔═╡ f5c90393-faa4-47dc-a55b-9d66dd1d2c25
dsamples=Dict("ZFFL21_C"=>(name="QUINN", dep="spinCoating", 
                                         c=1.0e-6, Ba2=false, 
                                         grass=false),
	"ZFFL21_C.2"=>(name="QUINN", dep="spinCoating", 
                                         c=1.0e-6, Ba2=false, 
                                         grass=false),
	"ZFFL21_D"=>(name="QUINN", dep="spinCoating", 
                                         c=1.05e-7, Ba2=false, 
                                         grass=false),
	"ZFFL21_E"=>(name="QUINN", dep="spinCoating", 
                                         c=1.0e-8, Ba2=false, 
                                         grass=false),
	"ZFFL21_C.3"=>(name="QUINN", dep="spinCoating", 
                                         c=1.0e-6, Ba2=false, 
                                         grass=false),
	"ZFFL21_D.3"=>(name="QUINN", dep="spinCoating", 
                                         c=1.05e-7, Ba2=false, 
                                         grass=false),
	"ZFFL21_E.3"=>(name="QUINN", dep="spinCoating", 
                                         c=1.0e-8, Ba2=false, 
                                         grass=false),
	"ZFFL21_F.3"=>(name="QUINN", dep="spinCoating", 
                                         c=1.0e-9, Ba2=false, 
                                         grass=false),
	"ZFFL21_G.3"=>(name="QUINN", dep="spinCoating", 
                                         c=1.0e-10, Ba2=false, 
                                         grass=false),
	"ZFFL21_H.3"=>(name="QUINN", dep="spinCoating", 
                                         c=1.0e-11, Ba2=false, 
                                         grass=false),
	"CleanQUARTZ_07072023"=>(name="Quartz", dep="clean", 
                                         c=0, Ba2=false, 
                                         grass=false)

	
	
)

# ╔═╡ 14192a57-c137-4a61-9074-98d6a0211ef9
dftpar=Dict(
	"ZFFL21_C_Filter0"=>(nbins=100, dtmax=1000.0, tmax=120.0, ymax=5*10^4, dymax=2*10^5, i0=[5,10, 15], l0=[7,15,100], dcut=100.0),
	"ZFFL21_C.2_Filter0"=>(nbins=100, dtmax=1000.0, tmax=120.0, ymax=5*10^4, dymax=2*10^5, i0=[5,10, 15], l0=[7,15,100], dcut=100.0),
	"ZFFL21_D_Filter0"=>(nbins=100, dtmax=1000.0, tmax=120.0, ymax=1*10^4, dymax=1*10^5, i0=[5,10, 15], l0=[7,15,100], dcut=100.0),
	"ZFFL21_E_Filter0"=>(nbins=100, dtmax=1000.0, tmax=120.0, ymax=5*10^3, dymax=1*10^5, i0=[5,10, 15], l0=[7,15,100], dcut=100.0),
	"ZFFL21_C.3_Filter0"=>(nbins=100, dtmax=1000.0, tmax=120.0, ymax=5*10^4, dymax=2*10^5, i0=[5,10, 15], l0=[7,15,100], dcut=100.0),
	"ZFFL21_D.3_Filter0"=>(nbins=100, dtmax=1000.0, tmax=120.0, ymax=1*10^4, dymax=1*10^5, i0=[5,10, 15], l0=[7,15,100], dcut=100.0),
	"ZFFL21_E.3_Filter0"=>(nbins=100, dtmax=1000.0, tmax=120.0, ymax=5*10^3, dymax=1*10^5, i0=[5,10, 15], l0=[7,15,100], dcut=100.0),
	"ZFFL21_F.3_Filter0"=>(nbins=100, dtmax=1000.0, tmax=120.0, ymax=1*10^3, dymax=2*10^4, i0=[5,10, 15], l0=[7,15,100], dcut=100.0),
	"ZFFL21_G.3_Filter0"=>(nbins=100, dtmax=1000.0, tmax=120.0, ymax=1*10^3, dymax=1*10^5, i0=[5,10, 15], l0=[7,15,100], dcut=100.0),
	"ZFFL21_H.3_Filter0"=>(nbins=100, dtmax=1000.0, tmax=120.0, ymax=1*10^3, dymax=1*10^5, i0=[5,10, 15], l0=[7,15,100], dcut=100.0),
	"CleanQUARTZ_07072023_Filter0"=>(nbins=100, dtmax=1000.0, tmax=120.0, ymax=2*10^3, dymax=1*10^5, i0=[5,10, 15], l0=[7,15,100], dcut=150.0),
)

# ╔═╡ 12a7c46a-ed5a-40fa-bb07-9dd3811d9dcd
md""" Select sample : $(@bind sample Select(sort(collect(keys(dsamples)))))"""

# ╔═╡ 507702eb-f000-48be-99da-b19180dcda6e
md"""
### Sample
- Sample name = $sample
#### Sample params: 
- Molecule name $(dsamples[sample].name)
- Deposition $(dsamples[sample].dep)
- Concentration $(dsamples[sample].c) M
- Grass? $(dsamples[sample].grass)
- Ba2+? $(dsamples[sample].Ba2)
"""

# ╔═╡ c792b017-9c33-4d03-853a-9672f8eaf481
begin
	sdir = string(rdir, sample)
	xfilters = readdir(sdir)
	filters = sort(filter(x -> x != ".DS_Store", xfilters))
	md""" Select filter : $(@bind xfilter Select(filters))"""
end

# ╔═╡ b748cd45-4b70-4383-991a-2dddd5799dfb
xfilter

# ╔═╡ f08c25a5-3b2b-4b43-a021-035ba9539e46
begin
	sdflt = joinpath(sdir, xfilter)
	ff = readdir(sdflt)
	xfiles = readdir(sdflt)
	files = filter(x -> x != ".DS_Store", xfiles)
	md""" Select file : $(@bind file Select(files))"""
end

# ╔═╡ ffd49f31-6a3d-4800-b257-0f1841c50706
begin
	zpath = joinpath(sdflt,file)
	md"""
	- path to file =$zpath
	"""
end

# ╔═╡ 9ae4acf2-a712-4c50-963d-e33754c01f20
begin
	sexp = string(sample, "_", xfilter)
	zdat = dftpar[sexp]

	md"""
	#### Fit pars
	"""
end

# ╔═╡ 547b9bba-592c-423e-be10-0775d60ed677
sexp

# ╔═╡ 04d53ef8-37a5-4d1b-9cdf-88782923e80a
zdat

# ╔═╡ 843c6950-d629-495b-9d9e-06d920d80229
PhotonDF, tagDict =  lfi.LaserLab.readHH(zpath, 100, true, false)

# ╔═╡ 11754d33-3316-4163-86c0-0ae25a7d679b
fDict = Dict("fctime"  =>tagDict["File_CreatingTime"],
            "acqtime"  =>tagDict["MeasDesc_AcquisitionTime"]*1e-3,
            "nrecords" =>tagDict["TTResult_NumberOfRecords"],
            "tbin"    =>tagDict["MeasDesc_GlobalResolution"] * 1e+9,
            "dtbin"  =>tagDict["MeasDesc_Resolution"]* 1e+9)

# ╔═╡ a0341482-339e-4c10-8055-54aa8d483fc8
begin
	ttrns = fDict["dtbin"]
	tns = fDict["tbin"]
	md"""
	#### Time bins:
	- dt  time bin(units of ns for time difference) = $ttrns
	- tag time bin(units of ns for time difference) = $tns
	"""
end

# ╔═╡ 0caf7df2-a7e3-4369-a77c-e923d6bd0772
histogram(PhotonDF.dtime * ttrns)

# ╔═╡ b2d6bb5e-96f4-4fae-8717-7026c652b839
histogram(PhotonDF.timeTag * tns*1e-9)

# ╔═╡ 1cedaad8-c2ae-4fba-8326-967d9e37ea9b
h1dtns  = lfi.LaserLab.h1d(PhotonDF.dtime *ttrns, zdat.nbins, 0.0, zdat.dtmax)

# ╔═╡ 5f143f8b-b72f-4c21-8bfc-17cf061653c9
pp =lfi.LaserLab.plot_h1dx(h1dtns; xs = "dt (ns)", ys = "a.u", 
                           xlim=true, xl = (0.0, zdat.dtmax),
						   ylim=true, yl=(1,zdat.dymax), ylog=false,
						   markersize=2, fwl=true)

# ╔═╡ c15de41a-15c5-4034-a0d4-77888b7784b5
begin
tdata = h1dtns.centers
vdata = h1dtns.weights
end

# ╔═╡ 39a7ba3a-72f4-4e84-8fe7-cbc62545aa9d
imx = argmax(vdata)

# ╔═╡ bfc13b99-e9cf-455f-8269-1e6234bf0e58
vdata[imx-3:imx+5]

# ╔═╡ df977c06-e6b2-47cf-bc93-87fc788473c6
#cdc, stdc, ydc = tfit_bkg(tdata, vdata; i0=zdat.i0[3], l0=zdat.l0[3], pa0=[100.0, 0.5])

# ╔═╡ 1f524d20-d336-4ae5-be00-781c057c1f9c
h1ts  = lfi.LaserLab.h1d(PhotonDF.timeTag * tns*1e-9, zdat.nbins, 0.0, zdat.tmax)

# ╔═╡ 63af7b4a-8738-4c89-965a-bd12655f2cc3
lfi.LaserLab.plot_h1dx(h1ts; xs = "t (s)", ys = "a.u", 
                           xlim=true, xl = (0.0, zdat.tmax),
						   ylim=true, yl=(1,zdat.ymax), markersize=2, fwl=true)

# ╔═╡ d50f9002-d3e6-4175-8351-7973dc5ea51c
md"""
## dt cutoff & bin time
"""

# ╔═╡ 5c2c75d9-123e-46dc-a70b-126f7e881589
dtcut = zdat.dcut # in dt units

# ╔═╡ 363898af-1ee7-4fa9-94db-393a777b0255
pdf = filter(:dtime  => x -> x * ttrns >= dtcut,  PhotonDF)

# ╔═╡ 9cdbf1a3-f510-4a0d-ac02-1a289eea9750
md"""
- cut of dt: efficiency = $(length(pdf.dtime) / length(PhotonDF.dtime))
"""

# ╔═╡ 2ddecf8b-1a57-483f-8891-500bb02081cf
begin
	h1dtfns  = lfi.LaserLab.h1d(pdf.dtime * ttrns, zdat.nbins, 0.0, zdat.dtmax)
	ppf =lfi.LaserLab.plot_h1dx(h1dtfns; xs = "dt (ns)", ys = "a.u", 
                           xlim=false, xl = (0.0, zdat.dtmax),
						   ylim=false, yl=(0,zdat.dymax), markersize=2, fwl=true)
end

# ╔═╡ dfb87e1b-4966-4bd0-bc31-886665ccb4ac
#csgn3, stdsgn3, ysgn3 = tfit_sgn(h1dtfns.centers, h1dtfns.weights; cbkg = [0.0,0.0], i0=50, l0=100, pa0=[100.0, 50.5])

# ╔═╡ 033e6f4f-3f50-474e-b92e-4bd52c1cb653
begin
	h2ts  = lfi.LaserLab.h1d(pdf.timeTag * tns*1e-9, zdat.nbins, 0.0, zdat.tmax)
	pts =lfi.LaserLab.plot_h1dx(h2ts; xs = "t (s)", ys = "a.u", 
                           xlim=true, xl = (0.0, zdat.tmax),
						   ylim=true, yl=(1,zdat.ymax), markersize=2, fwl=true)
end

# ╔═╡ 5f7027d2-521b-47f0-9573-d301e249d585
begin
	h3ts  = lfi.LaserLab.h1d(pdf.timeTag * tns*1e-9, 20, 0.0, zdat.tmax)
	p3ts =lfi.LaserLab.plot_h1dx(h2ts; xs = "t (s)", ys = "a.u", 
                           xlim=true, xl = (0.0, zdat.tmax),
						   ylim=true, yl=(1,100), markersize=2, fwl=true)
end

# ╔═╡ 6d9fb048-a608-4905-aaf5-b8bd90b239ab
md"""
### Quartz combined analysis
"""

# ╔═╡ b00caa31-64bb-4a74-a31c-b7ad98445746
xsamples =["ZFFL21_C.2", "ZFFL21_C.3", "ZFFL21_D.3", "ZFFL21_E.3", "ZFFL21_F.3", "ZFFL21_G.3", "ZFFL21_H.3"]

# ╔═╡ dff3ed15-5130-4e76-a0a6-6bf686e2ce16
#(sgntot = zsgn, sgnflttot = zsgnflt, sgn = dsgn, esgn = edsgn, sgnflt=dsgnflt,
	 #esgnflt=edsgnflt)

# ╔═╡ adef1149-cbf5-4403-88b0-e51196030cec
md"""
## Functions
"""

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
csgn1, stdsgn1, ysgn1 = tfit_sgn(tdata, vdata; cbkg = [0.0,0.0], i0=zdat.i0[1], l0=zdat.l0[1], pa0=[100.0, 20.5])

# ╔═╡ 785d56f4-e681-404a-b9b6-c139f7379745
csgn2, stdsgn2, ysgn2 = tfit_sgn(tdata, vdata; cbkg = [0.0,0.0], i0=zdat.i0[2], l0=zdat.l0[2], pa0=[100.0, 50.5])

# ╔═╡ 7d7ec498-1fe3-47db-b993-b0e699b6d435
csgn3, stdsgn3, ysgn3 = tfit_sgn(tdata, vdata; cbkg = [0.0,0.0], i0=zdat.i0[3], l0=zdat.l0[3], pa0=[100.0, 50.5])

# ╔═╡ 106dfd23-4fea-48b5-9d55-2b4d49186c93
begin
	pp2 = plot(pp, tdata, ysgn1)
	pp3 = plot(pp2, tdata, ysgn2)
	plot(pp3, tdata, ysgn3)
end

# ╔═╡ 0096d481-0f1c-4e61-95c1-3d90184de56a
md"""
### Fit results
- exp1: $(csgn1[2]) +- $(stdsgn1[2])
- exp2: $(csgn2[2]) +- $(stdsgn2[2])
- exp3: $(csgn3[2]) +- $(stdsgn3[2])
"""

# ╔═╡ 8e1d8113-204c-4421-9024-2b3922a6cc7d
ppf2 = plot(ppf, h1dtfns.centers, ysgn3)

# ╔═╡ 0cf824bc-506b-4b37-881c-97b051735653
function tfit_sgn2(xdata, ydata; cbkg, csgn, i0, l0, pa0=[100.0, 0.5])
	
	bkg(t) = cbkg[1] + cbkg[2] * t
	sgn(t) = csgn[1] * exp(-t/csgn[2])  +  bkg(t)
	ffun(t, N1, λ) = N1 * exp(-t/λ)  +  sgn(t)
	
	mffun(t, p) = p[1] * exp.(-t/p[2]) +  sgn.(t)

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

# ╔═╡ 5c729ea7-3249-4978-adc2-7e7ffbab4083
res = csgn(rdir, xsamples, zdat)

# ╔═╡ 58be53d0-1edb-4aa0-97da-d73766b3f663
res

# ╔═╡ 1132f573-90fb-4be7-873d-05a22699ee4c
res.sgn[1]/res.sgn[4]

# ╔═╡ 111362b5-95ff-472c-a047-7ec6106ffa62
begin
	yaxis = ("counts/second", (1,5.0e+4), :log)
	scatter(xsamples, res.sgn, yerr=res.esgn, yaxis=yaxis, label="QUIN ",  xlabel="sample")
	plot!(xsamples, res.sgnflt, yaxis=yaxis, label="")
	scatter!(xsamples, res.sgnflt, yerr=res.esgnflt, yaxis=yaxis, label="QUIN-cut ",  xlabel="sample")
	plot!(xsamples, res.sgnflt, yaxis=yaxis, label="")
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
# ╠═e3a9b4c9-778e-4c64-bd83-a9c9d32cee8e
# ╠═f5c90393-faa4-47dc-a55b-9d66dd1d2c25
# ╠═14192a57-c137-4a61-9074-98d6a0211ef9
# ╠═12a7c46a-ed5a-40fa-bb07-9dd3811d9dcd
# ╠═507702eb-f000-48be-99da-b19180dcda6e
# ╠═c792b017-9c33-4d03-853a-9672f8eaf481
# ╠═b748cd45-4b70-4383-991a-2dddd5799dfb
# ╠═f08c25a5-3b2b-4b43-a021-035ba9539e46
# ╠═ffd49f31-6a3d-4800-b257-0f1841c50706
# ╠═547b9bba-592c-423e-be10-0775d60ed677
# ╠═9ae4acf2-a712-4c50-963d-e33754c01f20
# ╠═04d53ef8-37a5-4d1b-9cdf-88782923e80a
# ╠═843c6950-d629-495b-9d9e-06d920d80229
# ╠═11754d33-3316-4163-86c0-0ae25a7d679b
# ╠═a0341482-339e-4c10-8055-54aa8d483fc8
# ╠═0caf7df2-a7e3-4369-a77c-e923d6bd0772
# ╠═b2d6bb5e-96f4-4fae-8717-7026c652b839
# ╠═1cedaad8-c2ae-4fba-8326-967d9e37ea9b
# ╠═5f143f8b-b72f-4c21-8bfc-17cf061653c9
# ╠═c15de41a-15c5-4034-a0d4-77888b7784b5
# ╠═39a7ba3a-72f4-4e84-8fe7-cbc62545aa9d
# ╠═bfc13b99-e9cf-455f-8269-1e6234bf0e58
# ╠═c0448a2e-c467-4cd6-aa44-82b488d5cf06
# ╠═785d56f4-e681-404a-b9b6-c139f7379745
# ╠═7d7ec498-1fe3-47db-b993-b0e699b6d435
# ╠═df977c06-e6b2-47cf-bc93-87fc788473c6
# ╠═106dfd23-4fea-48b5-9d55-2b4d49186c93
# ╠═0096d481-0f1c-4e61-95c1-3d90184de56a
# ╠═1f524d20-d336-4ae5-be00-781c057c1f9c
# ╠═63af7b4a-8738-4c89-965a-bd12655f2cc3
# ╠═d50f9002-d3e6-4175-8351-7973dc5ea51c
# ╠═5c2c75d9-123e-46dc-a70b-126f7e881589
# ╠═363898af-1ee7-4fa9-94db-393a777b0255
# ╠═9cdbf1a3-f510-4a0d-ac02-1a289eea9750
# ╠═2ddecf8b-1a57-483f-8891-500bb02081cf
# ╠═dfb87e1b-4966-4bd0-bc31-886665ccb4ac
# ╠═8e1d8113-204c-4421-9024-2b3922a6cc7d
# ╠═033e6f4f-3f50-474e-b92e-4bd52c1cb653
# ╠═5f7027d2-521b-47f0-9573-d301e249d585
# ╠═6d9fb048-a608-4905-aaf5-b8bd90b239ab
# ╠═b00caa31-64bb-4a74-a31c-b7ad98445746
# ╠═5c729ea7-3249-4978-adc2-7e7ffbab4083
# ╠═dff3ed15-5130-4e76-a0a6-6bf686e2ce16
# ╠═58be53d0-1edb-4aa0-97da-d73766b3f663
# ╠═1132f573-90fb-4be7-873d-05a22699ee4c
# ╠═111362b5-95ff-472c-a047-7ec6106ffa62
# ╠═adef1149-cbf5-4403-88b0-e51196030cec
# ╠═3204a561-8244-466e-9dbe-e1302ef363d9
# ╠═37a64b1f-93f6-4402-a32a-be1207352d86
# ╠═235d81bf-5166-4276-b61c-1e7586f4d8bf
# ╠═0cf824bc-506b-4b37-881c-97b051735653
# ╠═a349f303-3005-4c97-809b-3dde8af05ccb
# ╠═0883e1f1-4a94-4762-be8a-c97a99f2627f
