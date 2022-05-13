
"""
    scope_rawdata(WAVEDESC::Int)

Returns data needed to read PMT rawdata format

# Fields

- `WAVEDESC::String`   : Base address that serves as reference to all other addresses.    
"""
scope_rawdata(io::IOStream, WAVEDESC::Int) = Dict("WAVEDESC"         => WAVEDESC,
                                    "COMM_TYPE"	          => readword(io, WAVEDESC + 32), 
                                    "COMM_ORDER"          => readlong(io, WAVEDESC + 34), 
                                    "WAVE_DESCRIPTOR"     => readlong(io, WAVEDESC + 36), 
                                    "USER_TEXT"           => readlong(io, WAVEDESC + 40), 
                                    "TRIGTIME_ARRAY"      => readlong(io, WAVEDESC + 48), 
                                    "WAVE_ARRAY_1"	      => readlong(io, WAVEDESC+ 60), 
                                    "INSTRUMENT_NAME"     => readstring(io, WAVEDESC+ 76, 14), 
                                    "INSTRUMENT_NUMBER"	  => readlong(io, WAVEDESC+ 92), 
                                    "WAVE_ARRAY_COUNT"	  => readlong(io, WAVEDESC+ 116), 
                                    "VERTICAL_GAIN"	      => readfloat(io, WAVEDESC+ 156), 
                                    "VERTICAL_OFFSET"	  => readfloat(io, WAVEDESC+ 160), 
                                    "NOMINAL_BITS"	      => readword(io, WAVEDESC+ 172), 
                                    "HORIZ_INTERVAL"	  => readfloat(io, WAVEDESC+ 176), 
                                    "HORIZ_OFFSET"	      => readdouble(io, WAVEDESC+ 180), 
                                    "TRIGGER_TIME"	      => readtimestamp(io, WAVEDESC+ 296), 
                                    "RECORD_TYPE"	      => readword(io, WAVEDESC+ 316), 
                                    "PROCESSING_DONE"	  => readword(io, WAVEDESC+ 318), 
                                    "TIMEBASE"	          => readword(io, WAVEDESC+ 324), #rg
                                    "VERT_COUPLING"	      => readword(io, WAVEDESC+ 326), 
                                    "PROBE_ATT"	          => readfloat(io, WAVEDESC+ 328), 
                                    "FIXED_VERT_GAIN"	  => readword(io, WAVEDESC+ 332), #rg
                                    "BANDWIDTH_LIMIT"	  => readword(io, WAVEDESC+ 334), 
                                    "WAVE_SOURCE"	      => readword(io, WAVEDESC+ 344))


"""Returns the adress of the reference mark in the stream"""
function wavedesc(io::IOStream, strn="WAVEDESC", strl=20)
	seek(io, 0)
	header = Vector{UInt8}(undef, strl)
	readbytes!(io, header)
	header
	sh = String(header)
	result = findfirst(strn, sh)
	first(result) - 1 
end


""" reads two bytes at location s and unpack them as an int16"""
function readword(s::IOStream, pos::Int)
	seek(s, pos)
	bytes = Vector{UInt8}(undef, 2)
	readbytes!(s, bytes)
    collect(reinterpret(Int16, bytes))[1]
end


""" reads four bytes at location s and unpack them as an Int32"""
function readlong(s::IOStream, pos::Int)
	seek(s, pos)
	bytes = Vector{UInt8}(undef, 4)
	readbytes!(s, bytes)
    collect(reinterpret(Int32, bytes))[1]
end


""" reads four bytes at location s and unpack them as an Float32"""
function readfloat(s::IOStream, pos::Int)
	seek(s, pos)
	bytes = Vector{UInt8}(undef, 4)
	readbytes!(s, bytes)
    collect(reinterpret(Float32, bytes))[1]
end


""" reads eight bytes at location s and unpack them as an Float64"""
function readdouble(s::IOStream, pos::Int)
	seek(s, pos)
	bytes = Vector{UInt8}(undef, 8)
	readbytes!(s, bytes)
    collect(reinterpret(Float64, bytes))[1]
end


""" reads eight bytes at location s and unpack them as an Float64"""
function readstring(s::IOStream, pos::Int, strlen::Int)
	seek(s, pos)
	bytes = Vector{UInt8}(undef, strlen)
	readbytes!(s, bytes)
    String(bytes)
end


""" reads the time stamp"""
function readtimestamp(s::IOStream, pos::Int)
	seek(s,pos)
	bytes = Vector{UInt8}(undef, 8)
	readbytes!(s, bytes)
	seconds = collect(reinterpret(Float64, bytes))[1] 
	
	bytes = Vector{UInt8}(undef, 1)
	readbytes!(s, bytes)
	minutes = collect(reinterpret(Int8, bytes))[1]
	
	readbytes!(s, bytes)
	hours = collect(reinterpret(Int8, bytes))[1]
	
	readbytes!(s, bytes)
	days = collect(reinterpret(Int8, bytes))[1]
	
	readbytes!(s, bytes)
	months = collect(reinterpret(Int8, bytes))[1]
	
	bytes = Vector{UInt8}(undef, 2)
	readbytes!(s, bytes)
	years = collect(reinterpret(Int16, bytes))[1]
	
	(year=years,month=months,day=days,hour=hours,minute=minutes,second=seconds)
end


""" Returns the (time, ampl) data (t in mus ampl in mV)"""
function xydata(io::IOStream, addr::Dict{String, Any})
	yw = read_ywave(io, addr)
	wy = addr["VERTICAL_GAIN"] * yw .- addr["VERTICAL_OFFSET"]
	xw = collect(1:addr["WAVE_ARRAY_COUNT"]) * addr["HORIZ_INTERVAL"] .+ addr["HORIZ_OFFSET"]
	wx = -xw[1] .+ xw
	(time=wx *1e+6, ampl = -wy*1e+3 )
end


""" reads the waveform in y"""
function read_ywave(s::IOStream, addr::Dict{String, Any})
	seek(s, addr["WAVEDESC"] + addr["WAVE_DESCRIPTOR"] + addr["USER_TEXT"] + addr["TRIGTIME_ARRAY"])

	bytes = Vector{UInt8}(undef, addr["WAVE_ARRAY_1"])
	readbytes!(s, bytes)
	if addr["COMM_TYPE"] == 0  # file is type byte (Int8)
		collect(reinterpret(Int8, bytes))
	else
		collect(reinterpret(Int16, bytes))
	end
end


### Analysis functions 

"""
Read  weveform specifying an iostream  
"""
function read_waveform(zio::IOStream)
	WAVEDESC = wavedesc(zio)
	address = scope_rawdata(zio, WAVEDESC)
	xydata(zio, address)
end


"""
Reads  weveform specifying directory and file name 
"""
function read_waveform(pmtdir::String, fname::String)
    zpath = joinpath(pmtdir,fname)
	zio = open(zpath, "r")
	wvfm = read_waveform(zio)
	close(zio)
	wvfm 
end


"""
Read waveform specified by fnumber taking the full path from a vector of strings 
"""
function read_waveform(csvf::Vector{String}, fnumber::Integer)
	zio = open(csvf[fnumber], "r")
	wvfm = read_waveform(zio)
	close(zio)
	wvfm 
end


"""Return the mean, std, and (+-) threshold at nsigma, for raw or filtered waveform"""
function wstats(wvfm::NamedTuple{(:time, :ampl)}; nsigma::Float64=1.0)
	meanfs = mean(wvfm.ampl)
	stdfs = std(wvfm.ampl)
	thrp = meanfs + stdfs * nsigma
	thrn = meanfs - stdfs * nsigma
	(mean=meanfs, std=stdfs, thrp=thrp, thrn=thrn)
end


"""Return the sampling rate"""
function sampling_rate(wvfm::NamedTuple{(:time, :ampl)}, wmus::typeof(1.0μs)) 
	uconvert(GHz, (length(wvfm.time) / (wmus/1.0s)) * Hz)
end


"""Return the sampling period"""
function sampling_period(wvfm::NamedTuple{(:time, :ampl)}, wmus::typeof(1.0μs)) 
	uconvert(ns, 1.0 / sampling_rate(wvfm, wmus))
end


"""Return the fourier transform of the waveform wvfm"""
function fourier_transform(wvfm::NamedTuple{(:time, :ampl)}, wmus::typeof(1.0μs))
	#ts = sampling_period(wvfm, wmus)
	#tss = uconvert(s, ts)/s           # tss is the sampling period in seconds
    sr = uconvert(Hz, sampling_rate(wvfm, wmus)) /Hz
	F = fft(wvfm.ampl) |> fftshift   # fourier transform
	#freqs = fftfreq(length(wvfm.time), 1.0/tss) |> fftshift  # frequencies in Hz
    freqs = fftfreq(length(wvfm.time), sr) |> fftshift  # frequencies in Hz
	F, freqs
end


"""
	Returns the amplitude filtered by a low pass filter 

*Fields*
  
`wvfm`        : Waveform 
`wmus`        : width of time window in μs.
`filtertype`  : two possible filters: Butterworth or Chebyshev1
`flhz`        : higer frequency accepted by the filter in Hz
`n`           : filter order
`ϵ`           : ripple factor for Chebyshev1 filter
"""
function filter_signal_lp(wvfm::NamedTuple{(:time, :ampl)}, wmus::typeof(1.0μs); 
                          filtertype::String, flhz::Float64, n=4, ϵ=1.0)

	fs = uconvert(Hz, sampling_rate(wvfm, wmus))/Hz
	_, _, fampl = filter_signal_lowpass(wvfm, flhz, fs, filtertype, n, ϵ)
    
    (time=wvfm.time, ampl=fampl)
end


"""
Returns frequency response, filter (Butterworth, or Chebyshev) and filtered signal

*Fields*

`flh`         : higer frequency accepted by the filter in Hz
`fs`          : sampling rate in Hz 
`filtertype`  : two possible filters: Butterworth or Chebyshev1
`n`           : filter order
`ϵ`           : ripple factor for Chebyshev1 filter

"""
function filter_signal_lowpass(wvfm::NamedTuple{(:time, :ampl)}, 
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
	H, filter, filt(filter, wvfm.ampl)
end


"""
Selects the peaks in the waveform above prominence cut 

	Fields
  
`wvfm`        : Waveform 
`thrp`        : threshold for peak search.
`promsel`     : cut on peak prominence
`wsel`        : cut on peak width

Returns

A filtered waveform and a SPeaks struct 

"""
function select_filtered_peaks(wvfm::NamedTuple{(:time, :ampl)}, thrp::Float64; 
                               promsel::Float64, wsel=0.0)

    amplth = wvfm.ampl[wvfm.ampl .> thrp]

	#wvfmFlt = filter(row -> row[:fAmpl] > thrp, wvfm)
	pks, _ = findmaxima(amplth)
	peaks2, proms = peakproms(pks,amplth; minprom=promsel)
	
	if length(peaks2) == 0
		return Nothing 
	end
	
	peaks2, widths, leftedge, rightedge = peakwidths(peaks2, amplth, proms; minwidth=wsel)

	ys = [amplth[i] for i in peaks2]
	xs = [wvfm.time[i] for i in peaks2]

	npeaks = length(peaks2)
	xbase = ones(npeaks,1)*[promsel]
	
	return (ampl = amplth, speaks = SPeaks(peaks2, proms, widths, leftedge, rightedge, xs,ys, xbase, promsel))
end


"""
Computes the proms and widths of the events in the input vector 
"""
function proms(csvf::Vector{String}, fi::Integer, fe::Integer; 
	           flhz::Float64, promsel::Float64, 
	           nsigma=2.0, filtertype="Butterworth")

	spks0 = select_peaks(csvf, fi, fe;  flhz=flhz,
				         promsel=promsel, nsigma=nsigma, filtertype=filtertype)

	sP = reduce(vcat,[spks0[i].proms for i in 1:length(spks0) if spks0[i] !=Nothing])
	sW = reduce(vcat,[spks0[i].widths for i in 1:length(spks0) if spks0[i] !=Nothing])
	return sP, sW
end

"""
Returns a vector o SPeaks, one per event in the input vector
"""
function select_peaks(csvf::Vector{String}, i0::Integer, il::Integer; 
                      flhz::Float64, promsel::Float64, 
	                  nsigma=2.0, filtertype="Butterworth")

	ndim = il - i0 + 1
    SPK = Vector{SPeaks}(undef, ndim)
	for fileNumber in i0:il
		wvfm   = read_waveform(csvf,fileNumber)
        tw     = (wvfm.time[end] - wvfm.time[1]) * μs
        fwvfm  = filter_signal_lp(wvfm, tw; filtertype= filtertype, flhz=flhz)
        fstats = wstats(fwvfm; nsigma=nsigma)

		result = select_filtered_peaks(fwvfm, fstats.thrp; promsel=promsel, wsel=0.0)
        if result == Nothing 
            SPK[fileNumber] = SPeaks([1], [0.0], [0.0], [0.0], [0.0], [0.0],[0.0], [0.0], 0.0)
        else 
            SPK[fileNumber] = result.speaks 
        end
	end
	SPK
end


### Plots

"""
Plot the  waveform wvfm
"""
function plot_waveform(wvfm::NamedTuple{(:time, :ampl)}, stats::NamedTuple{(:mean, :std, :thrp, :thrn)}; 
                       window=false, sct=false, trace=false, ws=-1, we=1)

    xend = length(wvfm.time)
    if window
        range  = ws:we
    else
        range = 1:xend 
    end

    if sct
        p1 = scatter(wvfm.time[range], wvfm.ampl[range],  markersize=2,
        color = :black,
        label=false,
        fmt = :png)
    elseif trace
        p1 = plot(wvfm.ampl, lw=2, label=false, fmt = :png)
    else
        p1 = plot(wvfm.time[range], wvfm.ampl[range], lw=2, label=false, fmt = :png)
    end

    hline!([stats.thrp], label=false, fmt = :png)
    hline!([stats.thrn], label=false, fmt = :png)
    xlabel!("t (μs)")
    ylabel!("A (mV)")
    p1
end


"""
Plot the fourier transform of the waveform wvfm
"""
function plot_fft(wvfm::NamedTuple{(:time, :ampl)}, wmus::typeof(1.0μs); sct=true, window=false, ws=-1, we=1)
	F, freqs = fourier_transform(wvfm, wmus)
	xl = Integer(length(freqs))  #positive freqs
	xi = Integer(xl/2)
	fmhz = freqs[xi:xl] * 1.0e-6      # frequencies in MHz
	ampl = abs.(F)[xi:xl]  # abs amplitude 


    if window
        xi = xi + ws
        xl = min(xl, we)
    end
    fmhz = freqs[xi:xl] * 1.0e-6      # frequencies in MHz
	ampl = abs.(F)[xi:xl]  # abs amplitude 
    

    if sct
        p1 = scatter(fmhz, ampl,  markersize=2,
                     color = :black,
                     label=false,
                     fmt = :png)
    else
        p1 = plot(fmhz, ampl, lw=2, label=false, fmt = :png)
    end
	
	xlabel!("frequency (MHz)")
	ylabel!("A")
    p1
end

"""
	Plots filtered peaks 
"""
function plot_filtered_peaks(ampl::Vector{Float64}, spks::SPeaks)
	plot(ampl, lw=2, label=false)
	scatter!(spks.peaks, spks.ys, label="peak")
	scatter!(spks.leftedge, spks.xbase, label="leftedge")
	scatter!(spks.rightedge, spks.xbase, label="rightedge")
	hline!([spks.promsel], label=false)
	xlabel!("peak index")
	ylabel!("peak amplitude")
end


