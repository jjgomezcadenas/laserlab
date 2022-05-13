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
	dirs = lfi.LaserLab.getbolddirs(pmtdir)
	md""" Select experiment : $(@bind sexp Select(dirs))"""
end

# ╔═╡ ce8097c9-9608-49c8-aa5b-7338d090d743
begin
	namex = split(sexp,"_")
	csvdir = joinpath(odir, string(namex[1], namex[2]),"csv")
	pngdir = joinpath(odir, string(namex[1], namex[2]),"png")
	md"""
	- csv dir = $csvdir
	- png dir = $pngdir
	"""
end

# ╔═╡ 740681e5-7dc2-4fa1-bfb7-1ee5557ab885
let
	path = joinpath(pmtdir,sexp)
	dirs = lfi.LaserLab.getbolddirs(path)
	md""" Select filter : $(@bind sflt Select(dirs))"""
end

# ╔═╡ d3009c78-6318-4b5c-adfe-fb5905414081
begin
	xpath = joinpath(pmtdir,sexp,sflt)
	xfiles = Glob.glob(dfiles, xpath)
	md""" Select a file to analyze (1 to $(length(xfiles))): $(@bind xnf NumberField(0:length(xfiles), default=1))"""
end

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
	
	#r8n5k1m = select_filtered_peaks(wvfm8n5k1m, fthrp8n5k1m; promsel = sprom8n5k1m, wsel=0.0);
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


# ╔═╡ 799cb553-5093-4efc-91c6-5dc4161e6961
if ztimes
	ppcut = Float64(zpcut)
    spks  = lfi.LaserLab.select_peaks(xfiles, 1, nnf2; flhz=ffhz, promsel=ppcut, nsigma=3.0)
    times = reduce(vcat,[lfi.LaserLab.subtract_laser(spks[i].xs, tl/1.0μs) for i in 1:length(spks)])
    htimes = lfi.LaserLab.hist1d(times, 20, 0.0, 1.0)
	lfi.LaserLab.hist1d(htimes, "times")
end


# ╔═╡ 8ffb2af9-3d87-4e39-8340-88ab20054284
function gettimes(pmtdir, sexp; ffhz, pcut, nsigma, dfiles ="*.trc")
	
	path = joinpath(pmtdir,sexp)
	dirs = lfi.LaserLab.getbolddirs(path)
	T = []
	for dir in dirs
		xpath = joinpath(pmtdir,sexp,dir)
		xfiles = Glob.glob(dfiles, xpath)
		spks  = lfi.LaserLab.select_peaks(xfiles, 1, length(xfiles); flhz=ffhz, 
		                                  promsel=pcut, nsigma=nsigma)
    	times = reduce(vcat,
			[lfi.LaserLab.subtract_laser(spks[i].xs, tl/1.0μs) for i in 1:length(spks)])
		push!(T, times)
	end
	reduce(vcat, [t for t in T])
end

# ╔═╡ 39f8b749-3c41-4bff-9140-22c0be92e8dc
#ttimes = gettimes(pmtdir, sexp; ffhz=ffhz, pcut=2.0, nsigma=3.0)

# ╔═╡ 3e2d1b5f-de3a-4149-a90f-cf56a2c5a244
#httimes = lfi.LaserLab.hist1d(ttimes, 20, 0.0, 1.0)

# ╔═╡ d6bce551-f730-4d6d-ae14-fdd2730d68bb
#lfi.LaserLab.hist1d(httimes, "times")

# ╔═╡ 268f59b5-5b40-4ba0-a064-edb534bb3bab
md"""
# Functions
"""

# ╔═╡ e3b0ed46-122e-44d3-9a9a-01a69e032421
md"""
## Analysis functions
"""

# ╔═╡ 0b64c05c-9011-492d-8090-0be9e6cb6b6d


# ╔═╡ ec2cd649-7293-426c-ad9f-7ace1bf53f45
md"""
## Plotting functions
"""

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
# ╠═ce8097c9-9608-49c8-aa5b-7338d090d743
# ╠═740681e5-7dc2-4fa1-bfb7-1ee5557ab885
# ╠═d3009c78-6318-4b5c-adfe-fb5905414081
# ╠═db6868e3-9f96-49ef-b343-309a4b778591
# ╠═880dfedd-2df4-4d43-bd66-13962f910055
# ╠═b30c1537-465b-4938-97be-17219bef9402
# ╠═92e61e22-5666-43d9-a965-ef1d030a644b
# ╠═afeb9110-2058-4f09-9d27-4872f813256e
# ╠═00778f59-2f53-48db-82a0-2adc18fd1aab
# ╠═99e8afb0-f38e-4b79-88d2-6e43828b6e2b
# ╠═5d2cce97-ad00-46be-b086-8c4a35cdb008
# ╠═bd6faf80-4f68-44a0-9ad6-319a2b4f2316
# ╠═551cd586-e702-41fb-989c-5139632e29b6
# ╠═68bce42f-1417-4a8c-b8e7-05168fe33b9f
# ╠═c135833d-e9ae-41de-84a3-29b7b2e023ba
# ╠═cb8c269e-67f7-4239-bd56-066fce1886d2
# ╠═05fe45dc-0745-4d21-8648-a1257625fc72
# ╠═e6d7c71d-7eac-4f4e-af04-89d4d92d54f4
# ╠═0e1b6a85-af4b-4493-92e2-33aa4e1374e9
# ╠═67dac32c-11df-4843-9e52-e4845bc1a2d4
# ╠═799cb553-5093-4efc-91c6-5dc4161e6961
# ╠═8ffb2af9-3d87-4e39-8340-88ab20054284
# ╠═39f8b749-3c41-4bff-9140-22c0be92e8dc
# ╠═3e2d1b5f-de3a-4149-a90f-cf56a2c5a244
# ╠═d6bce551-f730-4d6d-ae14-fdd2730d68bb
# ╠═268f59b5-5b40-4ba0-a064-edb534bb3bab
# ╠═e3b0ed46-122e-44d3-9a9a-01a69e032421
# ╠═0b64c05c-9011-492d-8090-0be9e6cb6b6d
# ╠═ec2cd649-7293-426c-ad9f-7ace1bf53f45
