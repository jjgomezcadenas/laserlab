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

# ╔═╡ 74bc6fb9-8b52-46bf-b9cc-4a63883c5196
using Pkg; Pkg.activate("/Users/jj/JuliaProjects/LaserLab/")

# ╔═╡ f9cb42a4-188e-4929-8c23-b84cfe9df792
begin
	using PlutoUI
	using Statistics
	using StatsBase
	using Plots
	using Unitful
	using Peaks
	using Glob
	using FFTW
	using DSP
	using DataFrames
	using CSV
	using LsqFit
end

# ╔═╡ d87e4feb-a0d7-47ad-b9e2-2c92f5e467aa
import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
    A, N, mol, mmol, V, L, M

# ╔═╡ 948d5fee-cc84-11ec-1981-e71e1585c8bb
PlutoUI.TableOfContents(title="PMT Laser Analysis", indent=true)

# ╔═╡ 1a05bb05-5231-4168-8eab-ef84660deb9e
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

# ╔═╡ fca3e7c4-c160-43a1-ada1-78836dd2fe43
lfi = ingredients("../src/LaserLab.jl")

# ╔═╡ e56888f9-3747-4085-853f-6114956b72d3
md"""
# Analysis
"""

# ╔═╡ abc85752-30e4-48ac-87da-8b1050e05e43
begin
	tl = 1.0μs
	lf = uconvert(MHz, 1.0/tl)
	
md" - Parameters:
 - laser frequency (MHz) = $lf
 - Laser pulses every    = $tl
"
end

# ╔═╡ cc7ea918-066b-4485-80d2-64eba6af76a7
begin
	pmtdir  ="/Users/jj/JuliaProjects/LaserLab/labdata/PMT2"
	odir  = "/Users/jj/JuliaProjects/LaserLab/data/PMT"
	dfiles ="*.trc"
	dplots ="*.png"
	md"""
	CMOS dir = $pmtdir
	"""
	#close(io)
end

# ╔═╡ c907c77c-0291-4c80-bcba-ccbcd8e864cb
let
	readdir(pmtdir)
	dirs = lfi.LaserLab.getbolddirs(pmtdir)
	md""" Select experiment : $(@bind sexp Select(dirs))"""
end

# ╔═╡ a4fdfa7a-57bc-4d9f-ab90-bf99d3b48277
#sexp = "Bold_068_A1"

# ╔═╡ 740681e5-7dc2-4fa1-bfb7-1ee5557ab885
let
	path = joinpath(pmtdir,sexp)
	readdir(path)
	dirs = lfi.LaserLab.getbolddirs(path)
	md""" Select filter : $(@bind sflt Select(dirs))"""
end

# ╔═╡ d3009c78-6318-4b5c-adfe-fb5905414081
begin
	xpath = joinpath(pmtdir,sexp,sflt)
	xfiles = Glob.glob(dfiles, xpath)
	md""" Select a file to analyze (1 to $(length(xfiles))): $(@bind xnf NumberField(0:length(xfiles), default=1))"""
end

# ╔═╡ 513a0b87-409e-4e1e-b872-6c37fd84e80b
xfiles

# ╔═╡ db6868e3-9f96-49ef-b343-309a4b778591
begin
	zpath = joinpath(pmtdir,xfiles[xnf])
	zio = open(zpath, "r")
	WAVEDESC = lfi.LaserLab.wavedesc(zio)
	address = lfi.LaserLab.scope_rawdata(zio, WAVEDESC)
	wvfm = lfi.LaserLab.xydata(zio, address)
	tw = round(wvfm.time[end] - wvfm.time[1], sigdigits=2) * μs
	sr = lfi.LaserLab.sampling_rate(wvfm, tw)
	sp = lfi.LaserLab.sampling_period(wvfm, tw)
	stats = lfi.LaserLab.wstats(wvfm; nsigma=3.0)
	close(zio)
	md"""
	- Waveform length = $(tw) μs
	- sampling rate   = $(round(sr/GHz, sigdigits=2)) GHz
	- sampling period = $(round(sp/ns, sigdigits=4)) ns
	- baseline at     => $(round(stats.mean, sigdigits=3)) mV
	- baseline std    => $(round(stats.std, sigdigits=3)) mV
	"""
end

# ╔═╡ 880dfedd-2df4-4d43-bd66-13962f910055
typeof(wvfm.ampl)

# ╔═╡ b30c1537-465b-4938-97be-17219bef9402
md""" Select max filter frequency (in MHz): $(@bind ffmhz NumberField(0:1500, default=300))"""

# ╔═╡ 92e61e22-5666-43d9-a965-ef1d030a644b
begin
	ffhz = ffmhz*1e+6
    fwvfm = lfi.LaserLab.filter_signal_lp(wvfm, tw; filtertype= "Butterworth", 
	                                      flhz=ffhz)
	fstats = lfi.LaserLab.wstats(fwvfm; nsigma=3.0)
end

# ╔═╡ afeb9110-2058-4f09-9d27-4872f813256e
md""" Click to plot waveforms: $(@bind zpwf CheckBox(default=true))"""

# ╔═╡ 00778f59-2f53-48db-82a0-2adc18fd1aab
if zpwf
	pwvfm = lfi.LaserLab.plot_waveform(wvfm, stats; sct=true, window=false, ws=1, we=50000, trace=false)
	pfft = lfi.LaserLab.plot_fft(wvfm,tw; sct=true, window=false, ws=1, we=50000)
	pfwvfm = lfi.LaserLab.plot_waveform(fwvfm, fstats; sct=true, window=false, ws=1, we=50000, trace=false)
	
	plot(size=(850,550), pwvfm, pfft, pfwvfm, layout=(1,3), titlefontsize=8)
end

# ╔═╡ cb696494-5030-4568-a65b-aad9cbc01cf1
wvfm

# ╔═╡ 99e8afb0-f38e-4b79-88d2-6e43828b6e2b
md""" Select prominence cut: $(@bind spcut NumberField(0:10, default=2))"""

# ╔═╡ 5d2cce97-ad00-46be-b086-8c4a35cdb008
begin
	promsel = Float64(spcut)  #prominence cut
	 
	md"
	- Peak Prominence cut = $sprom8n5k1m
	"
	rfpeaks = lfi.LaserLab.select_filtered_peaks(fwvfm, fstats.thrp; 
                                       promsel=promsel, wsel=0.0)
	
	if rfpeaks != Nothing
		wvflt, fpeaks = rfpeaks
		lfi.LaserLab.plot_filtered_peaks(wvflt, fpeaks)
	end
end

# ╔═╡ bd6faf80-4f68-44a0-9ad6-319a2b4f2316
md"## Optimizing prominence cut"


# ╔═╡ 551cd586-e702-41fb-989c-5139632e29b6
md""" Click to compute proms histos: $(@bind zprom CheckBox(default=false))"""


# ╔═╡ 68bce42f-1417-4a8c-b8e7-05168fe33b9f
if zprom
	md""" Select min prominence cut: $(@bind spcutm NumberField(0.5:10.0, default=0.5))"""
end

# ╔═╡ c135833d-e9ae-41de-84a3-29b7b2e023ba
if zprom
	md""" Select how many files to analyze: $(@bind nnf NumberField(1:length(xfiles), default=ceil(length(xfiles)/2)))"""
end

# ╔═╡ cb8c269e-67f7-4239-bd56-066fce1886d2
if zprom
	promin = Float64(spcutm)
	if promin < 0.5
		promin = 0.5
	end
	proms, widths = lfi.LaserLab.proms(xfiles, 1, nnf;  
	                                   flhz=ffhz, promsel=promin, nsigma=3.0)
	hproms  = lfi.LaserLab.hist1d(proms, 20, 0.0, 10.0)
	hwidths = lfi.LaserLab.hist1d(widths, 20, 5.0, 15.0)

	phproms = lfi.LaserLab.hist1d(hproms, "Prominence (promsel=$promin)")
    phwidths =lfi.LaserLab.hist1d(hwidths, "Width")
	plot(phproms, phwidths, layout=(1,2), titlefontsize=8)
end

# ╔═╡ 05fe45dc-0745-4d21-8648-a1257625fc72
md"## Selection of events"

# ╔═╡ e6d7c71d-7eac-4f4e-af04-89d4d92d54f4
md""" Select  prominence cut: $(@bind zpcut NumberField(0.5:10.0, default=2.0))"""

# ╔═╡ 0e1b6a85-af4b-4493-92e2-33aa4e1374e9

md""" Select how many files to analyze: $(@bind nnf2 NumberField(1:length(xfiles), default=ceil(length(xfiles)/2)))"""

# ╔═╡ 67dac32c-11df-4843-9e52-e4845bc1a2d4
md""" Click to compute time histos: $(@bind ztimes CheckBox(default=false))"""


# ╔═╡ 01c601a6-1531-41df-b011-614c0b70d476
tl/1.0μs

# ╔═╡ 931aad69-ae0f-4168-8245-619590a046a2
md""" Select  prominence cut: $(@bind pcut NumberField(0.5:10.0, default=2.5))"""

# ╔═╡ dc3b9b68-4478-4f53-875d-3ecb6449edb8
coefw

# ╔═╡ 268f59b5-5b40-4ba0-a064-edb534bb3bab
md"""
# Functions
"""

# ╔═╡ 5c11564b-ba5c-49b4-b23e-6bef4fb242c1
function subtract_laser(tv::Vector{Float64}, ltus::Float64)
	tv - ltus * floor.(tv/ltus)
end

# ╔═╡ d557f2b0-7ccf-45c7-a93b-5bf499cfafad
if rfpeaks != Nothing
	subtract_laser(fpeaks.xs, tl/1.0μs)
end

# ╔═╡ 4c222697-e58c-4199-8114-38bc5a457cbb
function tfit(htime, pa0, i0, l0, ffun, mffun)
	
	tdata = htime.centers
	vdata = htime.weights
	
	fitb = curve_fit(mffun, tdata[i0:l0], vdata[i0:l0], pa0)

	fit = curve_fit(mffun, tdata[i0:l0], vdata[i0:l0], pa0)
	coef(fit), stderror(fit), ffun.(tdata, coef(fit)...)
end

# ╔═╡ 2e4d3f1b-bf12-4729-adcb-8ffb95902ec6
function tfit_bkg(htime; pa0=[100.0, 0.5], i0=10, l0=20)
	
	ffun(t, a,b) = a + b * t
	mffun(t, p) = p[1] .+ t .* p[2]

	tfit(htime, pa0, i0, l0, ffun, mffun)
end

# ╔═╡ 0d61f471-d272-4dfb-9cf9-45e08d7b1c82
function tfit_sgn(htime; cbkg, pa0=[100.0, 0.5], i0=10, l0=20)
	
	bkg(t) = cbkg[1] + cbkg[2] * t
	
	ffun(t, N1, λ) = N1 * exp(-t/λ)  +  bkg(t)
	mffun(t, p) = p[1] * exp.(-t/p[2]) +  bkg.(t)

	tfit(htime, pa0, i0, l0, ffun, mffun)
end

# ╔═╡ ae323f97-3bda-4933-9652-68224c939938
function tfit_sgn2(htime; cbkg, csgn, pa0=[100.0, 0.5], i0=1, l0=2)
	
	bkg(t) = cbkg[1] + cbkg[2] * t
	sgn(t) = csgn[1] * exp(-t/csgn[2])  +  bkg(t)
	ffun(t, N1, λ) = N1 * exp(-t/λ)  +  sgn(t)
	
	mffun(t, p) = p[1] * exp.(-t/p[2]) +  sgn.(t)

	tfit(htime, pa0, i0, l0, ffun, mffun)
end

# ╔═╡ cce4cc88-be0b-47a2-b37f-08cc0e40cf70
function tfit_2exp(htime; pa0=[1000.0, 10., 1000.0, 10.0], i0, l0)
	
	ffun(t, N1, λ1, N2, λ2) = N1 * exp(-t/λ1) + N2 * exp(-t/λ2)
	mffun(t, p) = p[1] * exp.(-t/p[2]) + p[3] * exp.(-t/p[4]) 

	tfit(htime, pa0, i0, l0, ffun, mffun)
end

# ╔═╡ 368cd59a-e9e2-436f-872b-92ef5bb8c1d8
function fit_peaks(htime; pa0=[100.0, 0.5], i0=1, il=-1)
	expo(t, N, λ) = N*exp(-t/λ)
	mexp(t, p) = p[1] * exp.(-t/p[2])
	tdata = htime.centers
	vdata = htime.weights
	if il < 0
		il = length(tdata)
	end
	fit = curve_fit(mexp, tdata[i0:il], vdata[i0:il], pa0)
	coef(fit), stderror(fit), expo.(tdata, coef(fit)...)
end

# ╔═╡ 7fa04ed2-113a-4166-9f29-eef98393ab6d
function create_dir!(dir)
	if isdir(dir) == false
		mkdir(dir)
	end
end

# ╔═╡ e44d66e8-bac5-4ec0-9c94-aa0dbbe05d06
function output_dirs!(odir, sexp)
	create_dir!(joinpath(odir, sexp))
	csvdir = joinpath(odir, sexp,"csv")
	pngdir = joinpath(odir, sexp, "png")
	create_dir!(csvdir)
	create_dir!(pngdir)
	csvdir, pngdir
end

# ╔═╡ ce8097c9-9608-49c8-aa5b-7338d090d743
begin
	#namex = split(sexp,"_")
	csvdir, pngdir = output_dirs!(odir, sexp)
	
	md"""
	- csv dir = $csvdir
	- png dir = $pngdir
	"""
end

# ╔═╡ 8436f3b4-43eb-40bd-afa4-ad3dd87310ea
pmtdf = lfi.LaserLab.load_df_from_csv(csvdir, string(sflt,".csv"), lfi.LaserLab.enG);

# ╔═╡ 3e2d1b5f-de3a-4149-a90f-cf56a2c5a244
pmtFpr = filter(row -> row[:pr] > pcut, pmtdf);

# ╔═╡ 12d461dd-1e80-4aab-90f4-8fd58344ee73
begin
	hpr = lfi.LaserLab.h1d(pmtFpr.pr, 20, 0.0, 10.0)
	phpr = lfi.LaserLab.plot_h1d(hpr, "prom")
	htpr = lfi.LaserLab.h1d(pmtFpr.ts, 20, 0.0, tl/1.0μs)
	phtpr = lfi.LaserLab.plot_h1d(htpr, "times")
	hws = lfi.LaserLab.h1d(pmtFpr.ws, 20, 5.0, 10.0)
	phws = lfi.LaserLab.plot_h1d(hws, "width")
	hvs = lfi.LaserLab.h1d(pmtFpr.vs, 20, 0.0, 20.0)
	phvs = lfi.LaserLab.plot_h1d(hvs, "voltage")
	Nothing	
end

# ╔═╡ 39f8b749-3c41-4bff-9140-22c0be92e8dc
begin
	hprx = lfi.LaserLab.h1d(pmtdf.pr, 20, 0.0, 10.0)
	phprx = lfi.LaserLab.plot_h1d(hpr, "prom")
end

# ╔═╡ 8adbeffc-5ac4-43da-b668-9bc17aca8bfc
lfi.LaserLab.plot_h1d(htpr, "times"; i0=2)

# ╔═╡ a3843632-e13c-47df-a7c7-72b4d841dbad
plot(phpr, phtpr, phws, phvs, size=(750,750), layout=(2,2),titlefontsize=8)

# ╔═╡ adafbbc7-9535-4d04-96e3-27cb386770a1
coefx, stderx, tftx = tfit_bkg(htpr; pa0=[10.0, 0.5], i0=13, l0=20)

# ╔═╡ 1146ca40-9d0d-4779-8d62-1ae4f991ada6
coefy, stdery, tfty = tfit_sgn(htpr; cbkg=coefx, pa0=[1000.0, 0.5], i0=2, l0=20)

# ╔═╡ c6e5ed98-aecd-4bba-9911-e3ecbc57bfd2
coefy

# ╔═╡ e90932f1-9a1b-4a77-8750-1e8b1322b788
pngpath = joinpath(pngdir, string(sflt,".png"))

# ╔═╡ fc65ca7f-aac8-49e5-8c13-2c47fb6a4a83
function runsel(csvf::Vector{String}, fi::Integer, fe::Integer;
	            flhz::Float64, promsel::Float64, tlmus::Float64, nsigma::Float64)

	spks = lfi.LaserLab.select_peaks(csvf, fi, fe; flhz=flhz, promsel=promsel, 
	                                 nsigma=nsigma)
	df = DataFrame()
	df.pr = reduce(vcat,[spks[i].proms  for i in 1:length(spks)])
	df.ws    = reduce(vcat,[spks[i].widths  for i in 1:length(spks)])
	df.ts   = reduce(vcat,[subtract_laser(spks[i].xs, tlmus) for i in 1:length(spks)])
	df.vs   = reduce(vcat,[spks[i].ys  for i in 1:length(spks)])
	df
end

# ╔═╡ 799cb553-5093-4efc-91c6-5dc4161e6961
if ztimes
	ppcut = Float64(zpcut)
	tdf = runsel(xfiles, 1, nnf2;flhz=ffhz, promsel=ppcut, tlmus=tl/1.0μs, nsigma=3.0)
	
	dfpath = joinpath(csvdir,string(sflt,".csv"))
	CSV.write(dfpath, tdf)
end


# ╔═╡ e3b0ed46-122e-44d3-9a9a-01a69e032421
md"""
## Histogram functions
"""

# ╔═╡ 70217e57-b672-4f17-9d4e-ab27e9a47057
function add_h1d(h1::Histo1d, h2::Histo1d)
	Histo1d(h1.edges, h1.weights + h2.weights, h1.centers)
end

# ╔═╡ ec2cd649-7293-426c-ad9f-7ace1bf53f45
md"""
## Plotting functions
"""

# ╔═╡ 11ee36e8-2eca-4577-93eb-207d1efa4cbb
function plot_fit(htime, coeff, tft; savefig=false, fn="")
	tdata =htime.centers
	vdata =htime.weights
	ps1 = scatter(tdata, vdata, yerr=sqrt.(abs.(vdata)), markersize=2,
				  color = :black,
	    		  label="data",
				  fmt = :png)
	pp = plot(ps1, tdata, tft, lw=2, label="μ = $(round(coeff[2]*1000, sigdigits=2)) ns", fmt = :png)
	xlabel!("t (μs)")
	ylabel!("frequency")
	if savefig
		png(pp, fn)
	end
	return pp 
end

# ╔═╡ 2c8543f7-878d-4cdf-bc90-d1a0ef989dd3
function plot_fit2(htime, label, tft; savefig=false, fn="")
	tdata =htime.centers
	vdata =htime.weights
	ps1 = scatter(tdata, vdata, yerr=sqrt.(abs.(vdata)), markersize=2,
				  color = :black,
	    		  label="data",
				  fmt = :png)
	pp = plot(ps1, tdata, tft, lw=2, label=label, fmt = :png)
	xlabel!("t (μs)")
	ylabel!("frequency")
	if savefig
		png(pp, fn)
	end
	return pp 
end

# ╔═╡ 3e9b0679-decd-40b0-81a7-94853101ac85
let
	s1 = " cte = $(round(coefx[1]))\n"
   	ptime = plot_fit2(htpr, s1, tftx; savefig=false,fn=pngpath)
	plot(ptime)
end

# ╔═╡ 37ba3d41-8816-470d-9536-8133bcd05cfd
let
	s1 = " cte = $(round(coefx[1]))\n"
	s2 = " μ1 = $(round(coefy[2]*1000)) ns\n"
	label = string(s1,s2)
    ptime2 = plot_fit2(htpr, label, tfty; savefig=true,
	                 fn=pngpath)
	plot(ptime2)
end


# ╔═╡ 76cdb001-b78c-456b-b285-743abc1bf007
begin
	coefz, stderz, tftz = tfit_sgn2(htpr; cbkg=coefx, csgn=coefy, pa0=[1000.0, 0.5], i0=1, l0=2)
	sf1 = " cte = $(round(coefx[1]))\n"
	sf2 = " μ1 = $(round(coefy[2]*1000)) ns \n"
	sf3 = " μ2 = $(round(coefz[2]*1000)) ns \n"
	label = string(sf1,sf2,sf3)
    ptime3 = plot_fit2(htpr, label, tftz; savefig=true,
	                 fn=pngpath)
	plot(ptime3)
end

# ╔═╡ Cell order:
# ╠═74bc6fb9-8b52-46bf-b9cc-4a63883c5196
# ╠═f9cb42a4-188e-4929-8c23-b84cfe9df792
# ╠═d87e4feb-a0d7-47ad-b9e2-2c92f5e467aa
# ╠═948d5fee-cc84-11ec-1981-e71e1585c8bb
# ╠═1a05bb05-5231-4168-8eab-ef84660deb9e
# ╠═fca3e7c4-c160-43a1-ada1-78836dd2fe43
# ╠═e56888f9-3747-4085-853f-6114956b72d3
# ╠═abc85752-30e4-48ac-87da-8b1050e05e43
# ╠═cc7ea918-066b-4485-80d2-64eba6af76a7
# ╠═c907c77c-0291-4c80-bcba-ccbcd8e864cb
# ╠═a4fdfa7a-57bc-4d9f-ab90-bf99d3b48277
# ╠═ce8097c9-9608-49c8-aa5b-7338d090d743
# ╠═740681e5-7dc2-4fa1-bfb7-1ee5557ab885
# ╠═d3009c78-6318-4b5c-adfe-fb5905414081
# ╠═513a0b87-409e-4e1e-b872-6c37fd84e80b
# ╠═db6868e3-9f96-49ef-b343-309a4b778591
# ╠═880dfedd-2df4-4d43-bd66-13962f910055
# ╠═b30c1537-465b-4938-97be-17219bef9402
# ╠═92e61e22-5666-43d9-a965-ef1d030a644b
# ╠═afeb9110-2058-4f09-9d27-4872f813256e
# ╠═00778f59-2f53-48db-82a0-2adc18fd1aab
# ╠═cb696494-5030-4568-a65b-aad9cbc01cf1
# ╠═99e8afb0-f38e-4b79-88d2-6e43828b6e2b
# ╠═5d2cce97-ad00-46be-b086-8c4a35cdb008
# ╠═d557f2b0-7ccf-45c7-a93b-5bf499cfafad
# ╠═bd6faf80-4f68-44a0-9ad6-319a2b4f2316
# ╠═551cd586-e702-41fb-989c-5139632e29b6
# ╠═68bce42f-1417-4a8c-b8e7-05168fe33b9f
# ╠═c135833d-e9ae-41de-84a3-29b7b2e023ba
# ╠═cb8c269e-67f7-4239-bd56-066fce1886d2
# ╠═05fe45dc-0745-4d21-8648-a1257625fc72
# ╠═e6d7c71d-7eac-4f4e-af04-89d4d92d54f4
# ╠═0e1b6a85-af4b-4493-92e2-33aa4e1374e9
# ╠═67dac32c-11df-4843-9e52-e4845bc1a2d4
# ╠═01c601a6-1531-41df-b011-614c0b70d476
# ╠═799cb553-5093-4efc-91c6-5dc4161e6961
# ╠═8436f3b4-43eb-40bd-afa4-ad3dd87310ea
# ╠═39f8b749-3c41-4bff-9140-22c0be92e8dc
# ╠═931aad69-ae0f-4168-8245-619590a046a2
# ╠═3e2d1b5f-de3a-4149-a90f-cf56a2c5a244
# ╠═12d461dd-1e80-4aab-90f4-8fd58344ee73
# ╠═8adbeffc-5ac4-43da-b668-9bc17aca8bfc
# ╠═a3843632-e13c-47df-a7c7-72b4d841dbad
# ╠═e90932f1-9a1b-4a77-8750-1e8b1322b788
# ╠═adafbbc7-9535-4d04-96e3-27cb386770a1
# ╠═3e9b0679-decd-40b0-81a7-94853101ac85
# ╠═1146ca40-9d0d-4779-8d62-1ae4f991ada6
# ╠═37ba3d41-8816-470d-9536-8133bcd05cfd
# ╠═76cdb001-b78c-456b-b285-743abc1bf007
# ╠═dc3b9b68-4478-4f53-875d-3ecb6449edb8
# ╠═c6e5ed98-aecd-4bba-9911-e3ecbc57bfd2
# ╠═268f59b5-5b40-4ba0-a064-edb534bb3bab
# ╠═5c11564b-ba5c-49b4-b23e-6bef4fb242c1
# ╠═2e4d3f1b-bf12-4729-adcb-8ffb95902ec6
# ╠═0d61f471-d272-4dfb-9cf9-45e08d7b1c82
# ╠═ae323f97-3bda-4933-9652-68224c939938
# ╠═cce4cc88-be0b-47a2-b37f-08cc0e40cf70
# ╠═4c222697-e58c-4199-8114-38bc5a457cbb
# ╠═368cd59a-e9e2-436f-872b-92ef5bb8c1d8
# ╠═7fa04ed2-113a-4166-9f29-eef98393ab6d
# ╠═e44d66e8-bac5-4ec0-9c94-aa0dbbe05d06
# ╠═fc65ca7f-aac8-49e5-8c13-2c47fb6a4a83
# ╠═e3b0ed46-122e-44d3-9a9a-01a69e032421
# ╠═70217e57-b672-4f17-9d4e-ab27e9a47057
# ╠═ec2cd649-7293-426c-ad9f-7ace1bf53f45
# ╠═11ee36e8-2eca-4577-93eb-207d1efa4cbb
# ╠═2c8543f7-878d-4cdf-bc90-d1a0ef989dd3
