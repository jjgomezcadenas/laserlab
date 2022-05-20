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
	tl = 2.0μs
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
	md""" Select experiment : $(@bind sexp2 Select(dirs))"""
end

# ╔═╡ a4fdfa7a-57bc-4d9f-ab90-bf99d3b48277
sexp = "Bold_068_A1"

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

# ╔═╡ 31f7fbcb-0361-4e46-a687-5cc1a3b5d4ce
fwvfm.time[1733]

# ╔═╡ 563ca292-24d9-4358-8c31-04e1267e7c37
fpeaks.xs

# ╔═╡ 2592c338-40a1-42da-8a1b-b6d5352b5543
fpeaks.ys

# ╔═╡ a4ed3140-fd53-4cc6-9fcf-88683abbdda2
fpeaks.peaks

# ╔═╡ 9a97a223-5323-40a2-aac5-eee4d46b2aab
fpeaks.leftedge

# ╔═╡ 5c11564b-ba5c-49b4-b23e-6bef4fb242c1
function subtract_laser(tv::Vector{Float64}, ltus::Float64)
	tv - ltus * floor.(tv/ltus)
end

# ╔═╡ d557f2b0-7ccf-45c7-a93b-5bf499cfafad
subtract_laser(fpeaks.xs, tl/1.0μs)

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

# ╔═╡ 55deba27-0dc9-4de0-b020-e782430cc55c
x = [10,100,4,5,3,10,20,500,1,2,3]

# ╔═╡ 5a08ae47-80da-4ed2-a342-fe05f0466e40
x[x .>3]

# ╔═╡ 4594cc43-bb6e-423a-a7b6-9fabe65c02eb
findall(z->z>3, x)

# ╔═╡ 268f59b5-5b40-4ba0-a064-edb534bb3bab
md"""
# Functions
"""

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
pmtdf = lfi.LaserLab.load_df_from_csv(csvdir, string(sflt,".csv"), lfi.LaserLab.enG)

# ╔═╡ 3e2d1b5f-de3a-4149-a90f-cf56a2c5a244
pmtFpr = filter(row -> row[:pr] > 2.5, pmtdf);

# ╔═╡ d6bce551-f730-4d6d-ae14-fdd2730d68bb
begin
	htpr = lfi.LaserLab.hist1d(pmtFpr.ts, 20, 0.0, tl/1.0μs)
	lfi.LaserLab.hist1d(htpr, "times")
end

# ╔═╡ cb3612d3-8376-456c-95df-a57fa492b35e
begin
	hws = lfi.LaserLab.hist1d(pmtFpr.ws, 20, 5.0, 10.0)
	lfi.LaserLab.hist1d(hws, "width")
end

# ╔═╡ 66eea7b2-6ce2-46ab-aad9-c53951215f27
begin
	hvs = lfi.LaserLab.hist1d(pmtFpr.vs, 20, 0.0, 20.0)
	lfi.LaserLab.hist1d(hvs, "voltage")
end

# ╔═╡ 356ed0ec-b08e-4f82-8df2-3363c5297e3b
pmtFpr[!,"txs"] = pmtFpr[!, "ts"] .- 200.0

# ╔═╡ 1ddddc13-1f9c-4eeb-8695-187b3044dffa
begin
	htx = lfi.LaserLab.hist1d(pmtFpr.txs, 20, 0.0, tl/1.0μs)
	lfi.LaserLab.hist1d(htx, "times")
end

# ╔═╡ 6756692a-75c9-4d4f-b99c-1cb2536ddc6e
pmtFpr[!, "ts"]

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

# ╔═╡ 77f1adeb-73c1-40c0-8690-222656939e11
mutable struct Histo1d
	edges::Vector{Float64}
	weights::Vector{Float64}
	centers::Vector{Float64}
end

# ╔═╡ 8b34877d-19af-4a8b-9435-80dccd137173
function get_histo1d(h::Histogram)
	function centers(edges::Vector{<:Real})
    	edges[1:end-1] + .-(@view(edges[2:end]), @view(edges[1:end-1])) / 2
	end	
	function edges(h::Histogram, index::Int64=1)
    	collect(h.edges[index])
	end
	centers(h::Histogram) = centers(edges(h))
	
	Histo1d(edges(h), h.weights, centers(h))
end

# ╔═╡ 22044fd8-7114-4b06-9615-33516e89a2de
h1tpr = get_histo1d(htpr)

# ╔═╡ dfe9c98f-2255-4716-85f6-6a49542129d3
plot_histo1d(h1tpr, "time")

# ╔═╡ 830961ce-098f-405a-abe8-978254a8856d
function h1d(x::Vector{<:Real}, bins::Vector{<:Real}, norm=false)
    h = fit(Histogram, x, bins)
    if norm
        h = StatsBase.normalize(h, mode=:density)
    end
    return get_histo1d(h)
end

# ╔═╡ 63e736f4-ef09-44c5-9ee0-4a0d43eb4196
function h1d(x::Vector{T}, nbins::Integer,
            xmin::T=typemin(T), xmax::T=typemax(T), norm=false) where T

	function in_range(x::Vector{<:Real}, xmin::Real, xmax::Real,
				  interval::Type{<:ValueBound}=OpenBound)
		mask = broadcast(range_bound(xmin, xmax, interval), x)
	    return x[mask]
	end
	
    xx   = in_range(x, xmin, xmax)
    bins = collect(LinRange(xmin, xmax, nbins + 1))
    h    = h1d(xx, bins, norm)
	
    return h
end

# ╔═╡ 0b64c05c-9011-492d-8090-0be9e6cb6b6d
function h1d_unitweight(h::Histo1d)
	w = ones(length(h.weights))
	Histo1d(h.edges, w, h.centers)
end

# ╔═╡ 9b3bc0bc-cf41-464b-8680-e82ffe763f8a
function h1d_constant_weight(h::Histo1d, wt::Float64)
	w = ones(length(h.weights))
	Histo1d(h.edges, wt .* w, h.centers)
end

# ╔═╡ e0875169-06bb-4931-96be-3649a606af5d
h1uw = h1d_constant_weight(h1tpr, -200.0)

# ╔═╡ 51f3eeb6-bdde-46c4-959e-483c9d09b1dc
plot_histo1d(h1uw, "time")

# ╔═╡ 70217e57-b672-4f17-9d4e-ab27e9a47057
function add_h1d(h1::Histo1d, h2::Histo1d)
	Histo1d(h1.edges, h1.weights + h2.weights, h1.centers)
end

# ╔═╡ d9bdc91d-0e1e-4ea5-b348-88ed298927e7
h1tx = add_h1d(h1tpr, h1uw)

# ╔═╡ 5540befa-8636-46f1-baf0-1327503ec399
plot_histo1d(h1tx, "time")

# ╔═╡ ec2cd649-7293-426c-ad9f-7ace1bf53f45
md"""
## Plotting functions
"""

# ╔═╡ f0006eb0-48a4-4285-a674-3958699ad821
function plot_h1d(h::Histo1d, xs::String;
                      markersize::Int64=3, fwl::Bool=false,
                      label::String="", legend::Bool=false)

    xg = h.centers
    yg = h.weights
    p  = scatter(xg, yg, yerr=sqrt.(abs.(yg)),
                 markersize=markersize, label=label, legend=legend)
    if fwl
        p = plot!(p, xg, yg, yerr=sqrt.(abs.(yg)), fmt=:png,
                  linewidth=1, label=label, legend=legend)
    end
    
    xlabel!(xs)
    ylabel!("frequency")
    return p
end

# ╔═╡ 39f8b749-3c41-4bff-9140-22c0be92e8dc
begin
	hpr = h1d(pmtdf.pr, 20, 0.0, 10.0)
	phpr = plot_h1d(hpr, "prom")
end

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

# ╔═╡ 3e9b0679-decd-40b0-81a7-94853101ac85
begin
	pngpath = joinpath(pngdir, string(sflt,".png"))
    coefx, stderx, tft = fit_peaks(h1tx; pa0=[100.0, 0.5], i0=1, il=8)
    ptime = plot_fit(h1tx, coefx, tft; savefig=true,
	              fn=pngpath)
	plot(ptime)
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
# ╠═31f7fbcb-0361-4e46-a687-5cc1a3b5d4ce
# ╠═563ca292-24d9-4358-8c31-04e1267e7c37
# ╠═2592c338-40a1-42da-8a1b-b6d5352b5543
# ╠═a4ed3140-fd53-4cc6-9fcf-88683abbdda2
# ╠═9a97a223-5323-40a2-aac5-eee4d46b2aab
# ╠═5c11564b-ba5c-49b4-b23e-6bef4fb242c1
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
# ╠═3e2d1b5f-de3a-4149-a90f-cf56a2c5a244
# ╠═d6bce551-f730-4d6d-ae14-fdd2730d68bb
# ╠═cb3612d3-8376-456c-95df-a57fa492b35e
# ╠═66eea7b2-6ce2-46ab-aad9-c53951215f27
# ╠═356ed0ec-b08e-4f82-8df2-3363c5297e3b
# ╠═1ddddc13-1f9c-4eeb-8695-187b3044dffa
# ╠═6756692a-75c9-4d4f-b99c-1cb2536ddc6e
# ╠═22044fd8-7114-4b06-9615-33516e89a2de
# ╠═dfe9c98f-2255-4716-85f6-6a49542129d3
# ╠═e0875169-06bb-4931-96be-3649a606af5d
# ╠═55deba27-0dc9-4de0-b020-e782430cc55c
# ╠═5a08ae47-80da-4ed2-a342-fe05f0466e40
# ╠═4594cc43-bb6e-423a-a7b6-9fabe65c02eb
# ╠═51f3eeb6-bdde-46c4-959e-483c9d09b1dc
# ╠═d9bdc91d-0e1e-4ea5-b348-88ed298927e7
# ╠═5540befa-8636-46f1-baf0-1327503ec399
# ╠═3e9b0679-decd-40b0-81a7-94853101ac85
# ╠═268f59b5-5b40-4ba0-a064-edb534bb3bab
# ╠═368cd59a-e9e2-436f-872b-92ef5bb8c1d8
# ╠═7fa04ed2-113a-4166-9f29-eef98393ab6d
# ╠═e44d66e8-bac5-4ec0-9c94-aa0dbbe05d06
# ╠═fc65ca7f-aac8-49e5-8c13-2c47fb6a4a83
# ╠═e3b0ed46-122e-44d3-9a9a-01a69e032421
# ╠═77f1adeb-73c1-40c0-8690-222656939e11
# ╠═8b34877d-19af-4a8b-9435-80dccd137173
# ╠═830961ce-098f-405a-abe8-978254a8856d
# ╠═63e736f4-ef09-44c5-9ee0-4a0d43eb4196
# ╠═0b64c05c-9011-492d-8090-0be9e6cb6b6d
# ╠═9b3bc0bc-cf41-464b-8680-e82ffe763f8a
# ╠═70217e57-b672-4f17-9d4e-ab27e9a47057
# ╠═ec2cd649-7293-426c-ad9f-7ace1bf53f45
# ╠═f0006eb0-48a4-4285-a674-3958699ad821
# ╠═11ee36e8-2eca-4577-93eb-207d1efa4cbb
