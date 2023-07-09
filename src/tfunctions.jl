

using Peaks
using Glob
using FFTW
using DSP
using LsqFit


"""
	SPeaks

A struct to store the peaks found by peak finding algorithm

# Fields

- `peaks::Vector{Int64}`          : indexes of the peaks 
- `proms::Vector{Float64}`        : prominence of the peaks 
- `widths::Vector{Float64}`       : width of the peaks
- `leftedge::Vector{Float64}`     : position of the left edge
- `rightedge::Vector{Float64}`    : position of the right edge
- `xs::Vector{Float64}`           : vector of `x` positions (times)
- `ys::Vector{Float64}`           : vector of `y` positions (voltages)
- `xbase::Vector{Float64}`        : pedestal
- `promsel::Float64`              : selection value of prmominence

"""
struct SPeaks
	peaks::Vector{Int64}
	proms::Vector{Float64} 
	widths::Vector{Float64} 
	leftedge::Vector{Float64} 
	rightedge::Vector{Float64}
	xs::Vector{Float64}
	ys::Vector{Float64}
	xbase::Vector{Float64}
	promsel::Float64
end

"""
	read_waveform(wdir::String, wname::String; wmus::Float64=5, symmetric::Bool=false)

reads a weveform specifying directory and file name 
 
# Fields
  
- `wdir`   : name of directory
- `wname`  : name of file
- `wmus`   : width of time window in μs.
- `wsym`   : If true the window is symmetric (centered at zero)
- `sign`   : reverses (the default) or not the waveform sign
- `vf`     : voltage factor, (the default is volts to mlVolts, 1000 )
- `tf`     : time factor, (the default is seconds to microseconds, 10^6 )

# Returns

Waveform as a data frame. Defaults: amplitude in mV (signal positive), time in μs.

"""
# function read_waveform(wdir::String, wname::String; wmus::Float64, wsym::Bool=true, 
# 	                   sign=-1.0, vf= 1e+3, tf=1e+6)
	
# 	wvfm = load_df_from_csv(wdir, wname, enG,5)
# 	shift_waveform!(wvfm, wmus, wsym, sign, vf, tf)
	
# end


# """
# 	read_waveform(csvf::Vector{String}, fnumber::Integer; wmus::Float64, wsym::Bool)

#   Like previous method but takes a vector of strings with the name of the files 
#   and an integer with the index of the file to be read

# """
# function read_waveform(csvf::Vector{String}, fnumber::Integer; wmus::Float64, 
# 	                   wsym::Bool)
	
# 	wvfm = DataFrame(CSV.File(csvf[fnumber], header=5, delim=enG.delim, decimal=enG.decimal))
# 	shift_waveform!(wvfm, wmus, wsym)
# end


"""
	shift_waveform!(wvfm:DataFrame, wmus::Float64, wsym::Bool)

 Shifts waveform to mV in Y (reverses sign) and to μs in x.
 If the window is centered in zero, shifts it, so that first time is zero.

 # Fields
  
- `wvfm`   : Waveform (a Data Frame with time/amplitude columns)
- `wmus`   : width of time window in μs.
- `wsym`   : If true the window is symmetric (centered at zero)
- `sign`   : reverses (the default) or not the waveform sign
- `vf`     : voltage factor, (the default is volts to mlVolts, 1000 )
- `tf`     : time factor, (the default is seconds to microseconds, 10^6 )

# Returns

Shifted waveform 

"""
function shift_waveform!(wvfm::DataFrame, wmus::Float64, wsym::Bool, sign=-1.0, vf= 1e+3, tf=1e+6)
	wvfm[!, :Ampl] = sign * vf * wvfm[!, :Ampl] 
	wvfm[!, :Time] = tf  .* wvfm[!, :Time]
	if wsym
		wvfm[!, :Time] = wmus/2.0  .+ wvfm[!, :Time]
	end
	wvfm
end


"""
		sampling_rate(wvfm::DataFrame, wmus::Float64)

# Fields
  
- `wvfm`  : Waveform (a Data Frame with time/amplitude columns)
- `wmus`  : width of time window in μs.

# Returns 

The sampling rate

"""
function sampling_rate(wvfm::DataFrame, wmus::Float64) 
	uconvert(GHz, (length(wvfm.Time) / (wmus*1.0e-6)) * Hz)
end


"""
sampling_period(wvfm::DataFrame, wmus::Float64)

# Fields
  
- `wvfm`  : Waveform (a Data Frame with time/amplitude columns)
- `wmus`  : width of time window in μs.

# Returns 

The sampling period

"""
function sampling_period(wvfm::DataFrame, wmus::Float64) 
	uconvert(ns, 1.0 / sampling_rate(wvfm, wmus))
end


"""
	wstats(wvfm::DataFrame; nsigma=1, waveform='filtered')

Return the mean, std, and (+-) threshold at nsigma, for raw or filtered waveform

"""
function wstats(wvfm::DataFrame; nsigma=1, waveform="filtered")
	if waveform=="filtered"
		ampl = "fAmpl"
	else
		ampl = "Ampl"
	end
	
	meanfs = mean(wvfm[!, ampl])
	stdfs = std(wvfm[!, ampl])
	thrp = meanfs + stdfs * nsigma
	thrn = meanfs - stdfs * nsigma
	meanfs, stdfs, thrp, thrn
end


"""
	fourier_transform(wvfm1::DataFrame, ws::Float64)
	
	Computes the fourier transform of wvfm1

	# Fields
  
- `wvfm1`   : Waveform
- `ws`      : width of the data window in μs 
  	
"""
function fourier_transform(wvfm1::DataFrame, ws::Float64)
	ts = sampling_period(wvfm1, ws)
	tss = uconvert(s, ts)/s           # tss is the sampling period in seconds
	F = fft(wvfm1.Ampl) |> fftshift   # fourier transform
	freqs = fftfreq(length(wvfm1.Time), 1.0/tss) |> fftshift  # frequencies in Hz
	ts, F, freqs
end


"""
	filter_signal_lowpass(wvfm::DataFrame, 
	                      flh::Float64, fs::Float64, filtertype::String, 
	                      n=4, ϵ=1.0)
  

`flh`         : higer frequency accepted by the filter in Hz
`fs`          : sampling frequency
`filtertype`  : two possible filters: Butterworth or Chebyshev1
`n`           : filter order
`ϵ`           : ripple factor for Chebyshev1 filter

Returns 

frequency response, filter (Butterworth, or Chebyshev) and filtered signal

"""
function filter_signal_lowpass(wvfm::DataFrame, 
	                           flh::Float64, fs::Float64, filtertype::String, 
	                           n=4, ϵ=1.0)

	responsetype = Lowpass(flh; fs=fs)
	
	if filtertype == "Butterworth"
		designmethod = Butterworth(n)
	elseif filtertype == "Chebyshev1"
		ripple = 10.0*log10(ϵ^2 + 1.0)
		designmethod = Chebyshev1(n, ripple)
	end
	filter = digitalfilter(responsetype, designmethod)
	H, _ = freqresp(filter)	
	H, filter, filt(filter, wvfm.Ampl)
end


"""
	filter_signal_lp!(wvfm::DataFrame, wmus::Float64; 
                           filtertype::String, flhz::Float64, 
                           n=4, ϵ=1.0)
Fields
  
`wvfm`        : Waveform (a Data Frame with time/amplitude columns)
`wmus`        : width of time window in μs.
`filtertype`  : two possible filters: Butterworth or Chebyshev1
`flhz`        : higer frequency accepted by the filter in Hz
`n`           : filter order
`ϵ`           : ripple factor for Chebyshev1 filter

Returns

A column with the filtered amplitude is added to the waveform
"""
function filter_signal_lp!(wvfm::DataFrame, wmus::Float64; 
                           filtertype::String, flhz::Float64, 
                           n=4, ϵ=1.0)
	fs = uconvert(Hz, sampling_rate(wvfm, wmus))/Hz
	_, _, fsgncb = filter_signal_lowpass(wvfm, flhz, fs, filtertype, n, ϵ)
	wvfm[!, "fAmpl"] = fsgncb
end


"""
	select_filtered_peaks(wvfm::DataFrame, thrp::Float64; 
                          promsel::Float64, wsel=0.0)


	Fields
  
`wvfm`        : Waveform (a Data Frame with time/amplitude columns)
`thrp`        : threshold for peak search.
`promsel`     : cut on peak prominence
`wsel`        : cut on peak width

Returns

A filtered waveform and a SPeaks struct 

"""
function select_filtered_peaks(wvfm::DataFrame, thrp::Float64; 
                               promsel::Number, wsel=0.0)

	wvfmFlt = filter(row -> row[:fAmpl] > thrp, wvfm)
	pks, _ = findmaxima(wvfmFlt.fAmpl)
	peaks2, proms = peakproms(pks, wvfmFlt.fAmpl; minprom=promsel)
	
	if length(peaks2) == 0
		return Nothing 
	end
	
	peaks2, widths, leftedge, rightedge = peakwidths(peaks2, wvfmFlt.fAmpl, proms; minwidth=wsel)

	ys = [wvfmFlt.fAmpl[i] for i in peaks2]
	xs = [wvfmFlt.Time[i] for i in peaks2]
	npeaks = length(peaks2)
	xbase = ones(npeaks,1)*[promsel]
	
	return wvfmFlt, SPeaks(peaks2, proms, widths, leftedge, rightedge, xs,ys, xbase, promsel)
end


"""
	select_peaks(csvf::Vector{String}, i0::Integer, il::Integer; 
	             wmus::Float64, wsym::Bool, flhz::Float64, promsel::Number, 
	             nsigma=2.0, filtertype="Butterworth")


Fields
  
`csvf`        : A vector of strings with the names of the files to be readout 
`i0`          : index of first file to read 
`il`          : index of last file to read 
`wmus`        : width of time window in μs.
`wsym`        : If true the window is symmetric (centered at zero)
`flhz`        : higer frequency accepted by the filter in Hz
`promsel`     : cut on peak prominence
`nsigma`      : number of sigmas to compute stats 
`filtertype`  : filter type  

Returns

A vector of SPeaks 

"""
# function select_peaks(csvf::Vector{String}, i0::Integer, il::Integer; 
# 	                  wmus::Float64, wsym::Bool, flhz::Float64, promsel::Number, 
# 	                  nsigma=2.0, filtertype="Butterworth")

# 	SPK = []
# 	for fileNumber in i0:il
# 		wvfm = read_waveform(csvf,fileNumber; wmus=wmus, wsym=wsym)
# 		filter_signal_lp!(wvfm, wmus; filtertype=filtertype, flhz=flhz, n=4)
# 		_, _, thrp, _ = wstats(wvfm; nsigma=nsigma, waveform="filtered")
# 		result = select_filtered_peaks(wvfm, thrp; promsel = promsel, wsel=0.0)

# 		if result != Nothing
# 			_, speaks = result
# 			push!(SPK,speaks)
# 		end
# 	end
# 	SPK
# end

"""
	proms(csvf::Vector{String}, fi::Integer, fe::Integer; 
	      wmus::Float64, wsym::Bool, flhz::Float64, promsel::Number, 
	      nsigma=2.0, filtertype="Butterworth")

Fields
  
`csvf`        : A vector of strings with the names of the files to be readout 
`i0`          : index of first file to read 
`il`          : index of last file to read 
`wmus`        : width of time window in μs.
`wsym`        : If true the window is symmetric (centered at zero)
`flhz`        : higer frequency accepted by the filter in Hz
`promsel`     : cut on peak prominence
`nsigma`      : number of sigmas to compute stats 
`filtertype`  : filter type  

Returns

A vector of proms and a vector of widths  

"""
# function proms(csvf::Vector{String}, fi::Integer, fe::Integer; 
# 	wmus::Float64, wsym::Bool, flhz::Float64, promsel::Number, 
# 	nsigma=2.0, filtertype="Butterworth")

# 	spks0 = select_peaks(csvf, fi, fe; wmus=wmus, wsym=wsym, flhz=flhz,
# 				promsel=promsel, nsigma=nsigma, filtertype=filtertype)

# 	sP = reduce(vcat,[spks0[i].proms for i in 1:length(spks0)])
# 	sW = reduce(vcat,[spks0[i].widths for i in 1:length(spks0)])
# 	return sP, sW
# end


"""
	runproms(csvf::Vector{String}, fi::Integer, fe::Integer;
             wmus::Float64, wsym::Bool=true, flhz::Float64, promsel::Number)

Fields
  
`csvf`        : A vector of strings with the names of the files to be readout 
`fi`          : index of first file to read 
`fe`          : index of last file to read 
`wmus`        : width of time window in μs.
`wsym`        : If true the window is symmetric (centered at zero)
`flhz`        : higer frequency accepted by the filter in Hz
`promsel`     : cut on peak prominence


Returns

A vector of proms and a vector of widths  

"""
function runsel(csvf::Vector{String}, fi::Integer, fe::Integer;
	            wmus::Float64, wsym::Bool, flhz::Float64, promsel::Number, tlmus::Float64)


	spks = select_peaks(csvf, fi, fe; wmus=wmus, wsym=wsym, flhz=flhz, promsel=promsel)
	df = DataFrame()
	df.pr = reduce(vcat,[spks[i].proms  for i in 1:length(spks)])
	df.ws    = reduce(vcat,[spks[i].widths  for i in 1:length(spks)])
	df.ts   = reduce(vcat,[subtract_laser(spks[i].xs, tlmus) for i in 1:length(spks)])
	df.vs   = reduce(vcat,[spks[i].ys  for i in 1:length(spks)])
	df
end


"""
	function subtract_laser(tv::Vector{Float64}, ltus::Float64)

Fields
  
`tv`   : A vector measured times 
`ltus` : time (period) of laser, in μs

Returns

A vector of times where the time of the laser has been subtracted. 

"""
function subtract_laser(tv::Vector{Float64}, ltus::Float64)
	tv - ltus * floor.(tv/ltus)
end


"""
	fit_peaks(htime; pa0=[100.0, 0.5], i0=1)
Fields

	`htime` : An histogram of times 
	`pa0`   : initial value of fit parameters 
	`i0`    : first bin to fit  

Returns

Fits an exponential to the histogramm, returns the coefficients
of the fit, the errors of the fit and a vector of predictions from the fit 
"""
function fit_peaks(htime; pa0=[100.0, 0.5], i0=1)
	expo(t, N, λ) = N*exp(-t/λ)
	mexp(t, p) = p[1] * exp.(-t/p[2])
	tdata = centers(htime)
	vdata = htime.weights
	il = length(tdata)
	fit = curve_fit(mexp, tdata[i0:il], vdata[i0:il], pa0)
	coef(fit), stderror(fit), expo.(tdata, coef(fit)...)
end


"""
	plot_waveform(wvfm; window=false, waveform='filtered', wstart=-1, wend=-1, thrp=0, thrn=0)

"""
function plot_waveform(wvfm::DataFrame; 
	                   window=false, waveform="filtered", wstart=-1, wend=-1, thrp=0, thrn=0, scatter=true)
	
	if waveform=="filtered"
		ampl = "fAmpl"
	else
		ampl = "Ampl"
	end
	
	if window
		p1 = plot(wvfm.Time[wstart:wend], wvfm[wstart:wend, ampl], lw=2, label=waveform, fmt = :png)
		if scatter
			scatter!(wvfm.Time[wstart:wend], wvfm[wstart:wend, ampl],  markersize=2,
					  color = :black, label=false, fmt = :png)
		end
	else
		p1 = plot(wvfm.Time, wvfm[!, ampl], lw=2, label=waveform, fmt = :png)
		if scatter
			scatter!(wvfm.Time, wvfm[!, ampl], markersize=2,
					  color = :black, label=false, fmt = :png)
		end
	end
	
	hline!([thrp], label=false, fmt = :png)
	hline!([thrn], label=false, fmt = :png)
	xlabel!("Time (μs or index)")
	ylabel!("Amplitude")
	p1
end

"""
	plot_fft(wvfm1, ws; window=false, is=-1, il=-1)

"""
function plot_fft(wvfm1::DataFrame, ws::Float64; window=false, is=-1, il=-1)
	ts, F, freqs = fourier_transform(wvfm1, ws)
	xl = Integer(length(freqs))
	xi = Integer(2 + xl/2)
	fmhz = freqs[xi:xl] * 1e-6      # frequencies in MHz
	ampl = abs.(F)[xi:xl]
	if window
		scatter(fmhz[is:il], ampl[is:il],  markersize=2,
				  color = :black,
	    		  label=false,
				  fmt = :png)
	else
		plot(fmhz, ampl)
	end
	
	xlabel!("frequency (MHz)")
	ylabel!("Amplitude")
end


"""
	plot_filtered_peaks(wvfmFlt, spks)

"""
function plot_filtered_peaks(wvfmFlt::DataFrame, spks::SPeaks)
	plot(wvfmFlt.fAmpl, lw=2, label=false)
	scatter!(spks.peaks, spks.ys, label="peak")
	scatter!(spks.leftedge, spks.xbase, label="leftedge")
	scatter!(spks.rightedge, spks.xbase, label="rightedge")
	hline!([spks.promsel], label=false)
	xlabel!("peak index")
	ylabel!("peak amplitude")
end


"""
	plot_fit(htime, coeff, tft)

"""
function plot_fit(htime, coeff, tft; savefig=false, fn="")
	tdata = centers(htime)
	vdata = htime.weights
	ps1 = scatter(tdata, vdata, yerr=sqrt.(vdata), markersize=2,
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

