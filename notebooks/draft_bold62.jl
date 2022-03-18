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

md"# BOLD\_062\_A2\_500kHz\_700nm\_100us"


md"## Read waveform"

fn7n5k1m = 1  # test file 

begin
	fls7n5k1m = glob("../labdata/BOLD_062_A2_500kHz_700nm_100us/C1*.csv")
	n7n5k1m = length(fls7n5k1m) 
	
	ffl7n5k1m = 200.0e+6
	wvfm7n5k1m = read_waveform(fls7n5k1m,fn7n5k1m; wmus=wmus100mus, wsym=true)
	first(wvfm7n5k1m,3)
end


md"
- Files in directory = $n7n5k1m
- sampling rate= $(round(sr7n5k1m, sigdigits=2)) GHz
- sampling period = $(round(sp7n5k1m, sigdigits=2) )ns
"

begin
	sr7n5k1m = sampling_rate(wvfm7n5k1m,wmus100mus)/GHz
	sp7n5k1m = sampling_period(wvfm7n5k1m,wmus100mus)/ns
end

md"
Unfiltered wvfm:
- mean = $(round(mean7n5k1m, sigdigits=4))
- std = $(round(std7n5k1m, sigdigits=4))
"
begin
	mean7n5k1m, std7n5k1m, thrp7n5k1m, thrn7n5k1m = thrx(wvfm7n5k1m; nsigma=3.0, waveform="raw");
	pfw7n5k1m = plot_filtered_waveform(wvfm7n5k1m; window=false, waveform="raw", thrp=thrp7n5k1m, thrn=thrn7n5k1m, scatter=false)
	pfft7n5k1m =plot_fft(wvfm7n5k1m,wmus100mus)
	plot(pfw7n5k1m, pfft7n5k1m)
end

md"### Filter signal
- Filter frequency (in Hz) = $ffl7n5k1m

Filtered wvfm:
- mean = $fmean7n5k1m
- std = $fstd7n5k1m
"
begin
    filter_signal_lp!(wvfm7n5k1m, wmus100mus; filtertype="Butterworth", flhz=ffl7n5k1m, n=4)
    first(wvfm7n5k1m,3)
end

begin
    fmean7n5k1m, fstd7n5k1m, fthrp7n5k1m, fthr7n5k1m = thrx(wvfm7n5k1m; nsigma=2.0, waveform="filtered")
    plot_filtered_waveform(wvfm7n5k1m; window=false, thrp=fthrp7n5k1m, thrn=fthr7n5k1m)
end

md"### Select peaks"

sprom7n5k1m = 2.0  #prominence cut
 
md"
- Peak Prominence cut = $sprom7n5k1m
"

begin
	r7n5k1m = select_filtered_peaks(wvfm7n5k1m, fthrp7n5k1m; promsel = sprom7n5k1m, wsel=0.0);
	if r7n5k1m != Nothing
		wvflt7n5k1m, spk7n5k1m = r7n5k1m
		plot_filtered_peaks(wvflt7n5k1m, spk7n5k1m)
	end
end

md"### Optimizing prominence cut"

begin
	P7n5k1m, W7n5k1m = proms(fls7n5k1m, 1, n7n5k1m;  
	                 wmus=wmus100mus, wsym=true, flhz=ffl7n5k1m, promsel=0.5)
	hP7n5k1m = lfi.LaserLab.hist1d(P7n5k1m, 20, 0.0, 10.0)
	hW7n5k1m = lfi.LaserLab.hist1d(W7n5k1m, 20, 5.0, 15.0)
end

lfi.LaserLab.hist1d(hP7n5k1m, "Prominence (promsel=0.5)")

lfi.LaserLab.hist1d(hW7n5k1m, "Width")

md"### Selection"

begin
    spks7n5k1m = select_peaks(fls7n5k1m, 1, n7n5k1m; wmus=wmus100mus, wsym=true, flhz=ffl7n5k1m, promsel=sprom7n5k1m)
    t7n5k1m = reduce(vcat,[subtract_laser(spks7n5k1m[i].xs, tlus2) for i in 1:length(spks7n5k1m)])
    hp7n5k1m = lfi.LaserLab.hist1d(t7n5k1m, 20, 0.0, 2.0)
end
lfi.LaserLab.hist1d(hp7n5k1m, "times")

md"### Fits"

begin
    cofe7n5k1m, stder7n5k1m, tft7n5k1m = fit_peaks(hp7n5k1m; pa0=[100.0, 0.5], i0=1)
    plot_fit(hp7n5k1m, cofe7n5k1m, tft7n5k1m)
end