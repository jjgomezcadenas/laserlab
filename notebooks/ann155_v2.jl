### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

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


# ╔═╡ 14192a57-c137-4a61-9074-98d6a0211ef9
dftpar=Dict(
	"Filter0"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^4, dymax=6*10^4, i0=
)

# ╔═╡ 79876fb7-7177-4277-acc6-71a8a5195868
dftparba=Dict(
	"Filter0"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=5*10^4, dymax=2*10^5, i0=
)

# ╔═╡ 8a57aa29-a9d9-47e5-b13b-34a62533e54d
dftparacn=Dict(
	"Filter0"=>(nbins=100, dtmax=5000.0, tmax=120.0, ymax=3*10^3, dymax=3*10^3, i0=[1,20, 1], l0=[20,100,100], dcut=100.0))

# ╔═╡ ffd49f31-6a3d-4800-b257-0f1841c50706
begin
	rdir ="/Users/jjgomezcadenas/BoldLab/BOLD/APD/ANN155_concentraciones_6/1e-5"
	fileba2="Ba_f0.ptu"
	filef  ="free1_f0.ptu"
	zpath = joinpath(rdir, filef)
	zpathba = joinpath(rdir,  fileba2)

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
	h1dtns0  = lfi.LaserLab.h1d(PhotonDF.dtime *ttrns, 50, 50.0, 10000.0)
	h1dtnsba0  = lfi.LaserLab.h1d(PhotonDFBa.dtime *ttrnsba, 50, 50.0, 	
			                      10000.0)
	pdtfx1 =lfi.LaserLab.plot_h1dx(h1dtns0; xs = "dt (ns)", ys = "a.u", 
                           xlim=true, xl = (0.0, 10000.0), label="C1",
						   ylim=true, yl=(1e+2,1e+6), ylog=true, legend=true,
						   markersize=1, fwl=false)
	ppdtbax1 =lfi.LaserLab.plot_h1dx(h1dtnsba0; xs = "dt (ns)", ys = "a.u", 
                           xlim=true, xl = (1.0, 10000.0), label="C1+Ba",
						   ylim=true, yl=(1e+2,1e+6), ylog=true, legend=true,
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
h1dtns  = lfi.LaserLab.h1d(PhotonDF.dtime *ttrns, 50, 100.0, 10000.0)
h1dtnsba  = lfi.LaserLab.h1d(PhotonDFBa.dtime *ttrnsba, 50, 100.0, 10000.0)
end

# ╔═╡ 5f143f8b-b72f-4c21-8bfc-17cf061653c9
begin
	ppdtf =lfi.LaserLab.plot_h1dx(h1dtns; xs = "dt (ns)", ys = "a.u", 
                           xlim=true, xl = (0.0, 10000.0), label="ANN155",
						   ylim=false, yl=(1,5.0e+4), ylog=false, legend=true,
						   markersize=2, fwl=true)
	ppdtba =lfi.LaserLab.plot_h1dx(h1dtnsba; xs = "dt (ns)", ys = "a.u", 
                           xlim=true, xl = (0.0, 10000.0), label="ANN155+Ba",
						   ylim=false, yl=(1,5.0e+5), ylog=false, legend=true,
						   markersize=2, fwl=true)
	plot(ppdtf, ppdtba)
end

# ╔═╡ bcd73cbc-96c4-4cbf-9bc6-9c7bbc56a018
begin
	tdata = h1dtns.centers
	vdataf0 = h1dtns.weights
	nb = 50
	tsgn = mean(vdataf0[nb-10:nb])
	
	vdataba0 = h1dtnsba.weights
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
	p3dtfx = scatter(tdata, vdataf, markersize=2, xlabel="t (ns)", ylabel="a.u.",
		yaxis = (:identity, (1.0, 10000)), label="C1")
	p3dtfbax = scatter(tdata, vdataba, markersize=2, xlabel="t (ns)", ylabel="a.u.",
		yaxis = (:identity, (1.0, 3e+4)), label="C1+Ba")
	p3dtfba = scatter(tdata, vdataba, markersize=2, yaxis = (:identity, (1.0, 10^5)), label="ANN155+Ba")
	p3dtf = scatter!(tdata, vdataf, markersize=2, yaxis = (:identity, (1.0, 10^5)),
		xlabel="t (ns)", ylabel="a.u.", label="ANN155")
end

# ╔═╡ f71f9286-0759-4315-8534-5bf2ce66ded6
plot(p3dtfx)

# ╔═╡ 76a7d6fb-e5f6-4ad8-a17b-bfc548d9fdc1
plot(p3dtfbax)

# ╔═╡ 84a56ba2-dba2-49e6-918d-715573db8328
zefit = fitexp(tdata, vdataf,n=2)

# ╔═╡ abeeb9dd-70cc-4d6e-88e6-4f71623b4d90
xefit = fitexp(tdata, vdataba,n=2)

# ╔═╡ 200deb4a-c9da-4460-8f8f-99ed4ddc9820
begin
	sexp1 = round(zefit.b[1], digits=1)
	sexp2 = round(zefit.b[2], digits=1)
	sexpb1 = round(xefit.b[1], digits=1)
	sexpb2 = round(xefit.b[2], digits=1)
	lnd1 = L"\lambda1 = %$sexp1 \, \mu s, \lambda2 = %$sexp2 \, \mu s"
	lnd2 = L"\lambda1 = %$sexpb1 \, \mu s, \lambda2 = %$sexpb2 \, \mu s"
	iexp1 = round(zefit.a[1]/(zefit.a[1]+zefit.a[2]), digits=2) 
	iexp2 = round(zefit.a[2]/(zefit.a[1]+zefit.a[2]), digits=2)
	iexpb1 = round(xefit.a[1]/(xefit.a[1]+xefit.a[2]), digits=2) 
	iexpb2 = round(xefit.a[2]/(xefit.a[1]+xefit.a[2]), digits=2)
	in1 = L"I(\lambda_1) = %$iexp1, I(\lambda_1) = %$iexp2"
	in2 = L"I(\lambda_2) = %$iexpb1, I(\lambda_2) = %$iexpb2"
	lbl1 = string(lnd1,"\n", in1) 
	lbl2 = string(lnd2,"\n", in2)
	
	zzpf3 = plot(p3dtfx, tdata, zefit.ypred,  ylim=(1, 4*10^4), 
		     label=lbl1)

	zzpba3 = plot(p3dtfbax, tdata, xefit.ypred, yaxis=:identity, ylim=(1, 2*10^5),label=lbl2)
	plot(zzpf3, zzpba3)
end

# ╔═╡ 0096d481-0f1c-4e61-95c1-3d90184de56a
md"""
### Fit results


- Free : exp1+exp2: λ1 = $(round(zefit.b[1])),  λ1 = $(round(zefit.b[2])), N1 = $(round(zefit.a[1])) N2 = $(round(zefit.a[2])) 

- Free: normalizations (combined fit): N1 = $(round(zefit.a[1]/(zefit.a[1]+zefit.a[2]), digits=2)), N2 = $(round(zefit.a[2]/(zefit.a[1]+zefit.a[2]), digits=2))

- Ba2+ : exp1+exp2: λ1 = $(round(xefit.b[1])),  λ1 = $(round(xefit.b[2])), N1 = $(round(xefit.a[1])) N2 = $(round(xefit.a[2])) 

- Ba2++: normalizations (combined fit): N1 = $(round(xefit.a[1]/(xefit.a[1]+xefit.a[2]), digits=2)), N2 = $(round(xefit.a[2]/(xefit.a[1]+xefit.a[2]), digits=2))
"""

# ╔═╡ 1f524d20-d336-4ae5-be00-781c057c1f9c
begin
	#h1ts  = lfi.LaserLab.h1d(PhotonDF.timeTag * tns*1e-9, zdat.nbins, 0.0, zdat.tmax)
	h1tsf  = lfi.LaserLab.h1d(PhotonDF.timeTag * tns*1e-9, 50, 0.0, 300.0)
	h1tsba  = lfi.LaserLab.h1d(PhotonDFBa.timeTag * tns*1e-9, 50, 
		                       0.0, 300.0)
	pztf = lfi.LaserLab.plot_h1dx(h1tsf; xs = "t (s)", ys = "a.u", 
                           xlim=true, xl = (0.0, 300.0), label="ANN155",
	                       legend=true, 
						   ylim=true, yl=(1,5e+4), markersize=2, fwl=true)
	pztba = lfi.LaserLab.plot_h1dx(h1tsba; xs = "t (s)", ys = "a.u", 
                           xlim=true, xl = (0.0, 300.0), label="ANN155+Ba",
	                       legend=true, 
						   ylim=true, yl=(1,5e+4), markersize=2, fwl=true)
	plot(pztf, pztba)
end

# ╔═╡ f14d4695-e7d4-4b4a-9393-9bd09aefe6e4
length(zefit.ypred)

# ╔═╡ 363898af-1ee7-4fa9-94db-393a777b0255
begin
	psn1 = plot(tdata, xefit.ypred ./zefit.ypred, lw=2, label="C1+Ba")
	plot(psn1)
end

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
	py1 = plot(xxt, y2b, lw=2, yaxis = (:log,(1e+3, 5e+7)), label="C1+Ba")
	py2 = plot(xxt, y2f, lw=2, yaxis = (:log,(1, 5e+6)), label="C1")
	py3 = plot(xxt,rsn, lw=2, yaxis = (:log,(1e+3, 1e+6)), label="S/N")
	py4 = plot(xxt,xef, lw=2, yaxis = (:log,(1e-3, 1.0)), label="eff")
	plot( py3, py4)
end

# ╔═╡ 86a3c22e-5675-4127-a61a-532f6bd52367
function plotsnsqrt()
	zxxt=0.0:10.0:5000.0
	zy2b = [qint(expb, t0, tl) for t0 in zxxt]
	zy2f = [qint(expf, t0, tl) for t0 in zxxt]
	zxef = [qint(expb, t0, tl) / qint(expb, 0.0, 5000.0) for t0 in zxxt]
	zrsn = zy2b ./ sqrt.(abs.(zy2f))
	zpy1 = plot(zxxt, zy2b, lw=2, xlabel="t (ns)", ylabel = L"\int I(t) dt",
	yaxis = (:log,(1e+3, 5e+8)), label="C1+Ba")
	zpy2 = plot(zxxt, zy2f, lw=2, xlabel="t (ns)", ylabel = L"\int I(t) dt",
		yaxis = (:log,(1, 5e+7)), label="C1")
	zpy3 = plot(zxxt,zrsn, lw=2, yaxis = (:log,(1e+4, 1e+9)), label="S/N")
	zpy4 = plot(zxxt,zxef, lw=2, yaxis = (:log,(1e-3, 1.0)), label="eff")
	plot( zpy1, zpy2, zpy3, zpy4)
	#plot(zpy1, zpy2, zpy4)
end

# ╔═╡ 7650d8d0-36ca-4141-b9b1-67aecd6fbe86
plotsnsqrt()

# ╔═╡ cdf91103-6022-4991-8c00-e4bb4f2287d9
begin
	emax = 3000.0
	zxxt=0.0:10.0:emax
	zy2b = [qint(expb, t0, tl) for t0 in zxxt]
	zy2f = [qint(expf, t0, tl) for t0 in zxxt]
	zxef = [100.0*qint(expb, t0, tl) / qint(expb, 0.0, emax) for t0 in zxxt]
	zrsn = zy2b ./ (abs.(zy2f))
	zpy1 = plot(zxxt, zy2b, lw=2, xlabel="t (ns)", ylabel = L"\int I(t) dt", 
		        yaxis = (:log,(1e+7, 1e+9)), label="C1+Ba")
	zpy2 = plot(zxxt, zy2f, lw=2, xlabel="t (ns)", ylabel = L"\int I(t) dt",
		        yaxis = (:log,(1, 1e+8)), label="C1")
	zpy3 = plot(zxxt,zrsn, lw=2, xlabel="t (ns)", 
		   ylabel = L"\int I_{Ba2^+}(t) dt/\int I_{free}(t) dt",
		   yaxis = (:log,(10, 1e+9)), label="S/N", legend=:left)
	zpy4 = plot(zxxt,zxef, lw=2, xlabel="t (ns)", ylabel="eff(%)",
		   yaxis = (:log,(1, 200.0)), label="eff")
	plot( zpy1, zpy2, zpy3, zpy4)
	#plot(zpy1, zpy2, zpy4)
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
# ╠═14192a57-c137-4a61-9074-98d6a0211ef9
# ╠═79876fb7-7177-4277-acc6-71a8a5195868
# ╠═8a57aa29-a9d9-47e5-b13b-34a62533e54d
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
# ╠═84a56ba2-dba2-49e6-918d-715573db8328
# ╠═abeeb9dd-70cc-4d6e-88e6-4f71623b4d90
# ╠═200deb4a-c9da-4460-8f8f-99ed4ddc9820
# ╠═0096d481-0f1c-4e61-95c1-3d90184de56a
# ╠═1f524d20-d336-4ae5-be00-781c057c1f9c
# ╠═f14d4695-e7d4-4b4a-9393-9bd09aefe6e4
# ╠═5c2c75d9-123e-46dc-a70b-126f7e881589
# ╠═86a3c22e-5675-4127-a61a-532f6bd52367
# ╠═7650d8d0-36ca-4141-b9b1-67aecd6fbe86
# ╠═cdf91103-6022-4991-8c00-e4bb4f2287d9
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
