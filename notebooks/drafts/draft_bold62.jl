begin
	f500kHz = 500.0kHz
	w100mus = 100.0μs
	wmus100mus=w100mus/μs # window in μs
	tlaser500kHz = uconvert(μs, 1.0/f500kHz)
	tlus2 = tlaser500kHz/μs
end

md" - Parameters:
 - laser frequency= $f500kHz
 - DAQ window= $w100mus
 - Laser pulses every $tlaser500kHz
"

md"# BOLD062\_A2\_500kHz\_730nm\_100us"

md"## Read waveform"

fn8n5k1m = 1  # test file 

begin
	fls8n5k1m = glob("../labdata/BOLD062_A2_500kHz_730nm_100us/C1*.csv")
	n8n5k1m = length(fls8n5k1m) 
	
	ffl8n5k1m = 200.0e+6
	wvfm8n5k1m = read_waveform(fls8n5k1m,fn8n5k1m; wmus=wmus100mus, wsym=true)
	first(wvfm8n5k1m,3)
end


md"
- Files in directory = $n8n5k1m
- sampling rate= $(round(sr8n5k1m, sigdigits=2)) GHz
- sampling period = $(round(sp8n5k1m, sigdigits=2) )ns
"

begin
	sr8n5k1m = sampling_rate(wvfm8n5k1m,wmus100mus)/GHz
	sp8n5k1m = sampling_period(wvfm8n5k1m,wmus100mus)/ns
end

md"
Unfiltered wvfm:
- mean = $(round(mean8n5k1m, sigdigits=4))
- std = $(round(std8n5k1m, sigdigits=4))
"
begin
	mean8n5k1m, std8n5k1m, thrp8n5k1m, thrn8n5k1m = thrx(wvfm8n5k1m; nsigma=3.0, waveform="raw");
	pfw8n5k1m = plot_filtered_waveform(wvfm8n5k1m; window=false, waveform="raw", thrp=thrp8n5k1m, thrn=thrn8n5k1m, scatter=false)
	pfft8n5k1m =plot_fft(wvfm8n5k1m,wmus100mus)
	plot(pfw8n5k1m, pfft8n5k1m)
end

md"### Filter signal
- Filter frequency (in Hz) = $ffl8n5k1m

Filtered wvfm:
- mean = $fmean8n5k1m
- std = $fstd8n5k1m
"
begin
    filter_signal_lp!(wvfm8n5k1m, wmus100mus; filtertype="Butterworth", flhz=ffl8n5k1m, n=4)
    first(wvfm8n5k1m,3)
end

begin
    fmean8n5k1m, fstd8n5k1m, fthrp8n5k1m, fthr8n5k1m = thrx(wvfm8n5k1m; nsigma=2.0, waveform="filtered")
    plot_filtered_waveform(wvfm8n5k1m; window=false, thrp=fthrp8n5k1m, thrn=fthr8n5k1m)
end

md"### Select peaks"

sprom8n5k1m = 2.0  #prominence cut
 
md"
- Peak Prominence cut = $sprom8n5k1m
"

begin
	r8n5k1m = select_filtered_peaks(wvfm8n5k1m, fthrp8n5k1m; promsel = sprom8n5k1m, wsel=0.0);
	if r8n5k1m != Nothing
		wvflt8n5k1m, spk8n5k1m = r8n5k1m
		plot_filtered_peaks(wvflt8n5k1m, spk8n5k1m)
	end
end

md"### Optimizing prominence cut"

begin
	P8n5k1m, W8n5k1m = proms(fls8n5k1m, 1, n8n5k1m;  
	                 wmus=wmus100mus, wsym=true, flhz=ffl8n5k1m, promsel=0.5)
	hP8n5k1m = lfi.LaserLab.hist1d(P8n5k1m, 20, 0.0, 10.0)
	hW8n5k1m = lfi.LaserLab.hist1d(W8n5k1m, 20, 5.0, 15.0)
end

lfi.LaserLab.hist1d(hP8n5k1m, "Prominence (promsel=0.5)")

lfi.LaserLab.hist1d(hW8n5k1m, "Width")

md"### Selection"

begin
    spks8n5k1m = select_peaks(fls8n5k1m, 1, n8n5k1m; wmus=wmus100mus, wsym=true, flhz=ffl8n5k1m, promsel=sprom8n5k1m)
    t8n5k1m = reduce(vcat,[subtract_laser(spks8n5k1m[i].xs, tlus2) for i in 1:length(spks8n5k1m)])
    hp8n5k1m = lfi.LaserLab.hist1d(t8n5k1m, 20, 0.0, 2.0)
end
lfi.LaserLab.hist1d(hp8n5k1m, "times")

md"### Fits"

begin
    cofe8n5k1m, stder8n5k1m, tft8n5k1m = fit_peaks(hp8n5k1m; pa0=[100.0, 0.5], i0=1)
    pp8n5k1m = plot_fit(hp8n5k1m, cofe8n5k1m, tft8n5k1m; savefig=true,
	              fn="b62A2_500kHz_730nm_100mus_fit.png")
	plot(pp8n5k1m)
end