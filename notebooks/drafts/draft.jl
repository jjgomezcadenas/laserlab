md"# BOLD\_055\_500kHZ\_732nm\_100us"

md"## Read waveform"

md" - Parameters:
 - laser frequency= $f500kHz
 - DAQ window= $w100mus
 - Laser pulses every $tlaser500kHz
 - sampling rate= $(round(sr73n5k1m, sigdigits=2)) GHz
 - sampling period = $(round(sp73n5k1m, sigdigits=2) )ns
"

fn73n5k1m = 1  # test file 

begin
	fls73n5k1m = glob("../labdata/BOLD_055_500kHZ_750nm_100us/C1*.csv")
	n73n5k1m = length(fls73n5k1m) 
	
	ffl73n5k1m = 200.0e+6
	wvfm73n5k1m = read_waveform(fls73n5k1m,fn73n5k1m; wmus=wmus100mus, wsym=true)
	first(wvfm73n5k1m,3)
end
begin
	sr73n5k1m = sampling_rate(wvfm73n5k1m,wmus100mus)/GHz
	sp73n5k1m = sampling_period(wvfm73n5k1m,wmus100mus)/ns
end
md"
Unfiltered wvfm:
- mean = $(round(mean73n5k1m, sigdigits=4))
- std = $(round(std73n5k1m, sigdigits=4))
"
begin
	mean73n5k1m, std73n5k1m, thrp73n5k1m, thrn73n5k1m = thrx(wvfm73n5k1m; nsigma=3.0, waveform="raw");
	pfw73n5k1m = plot_filtered_waveform(wvfm73n5k1m; window=false, waveform="raw", thrp=thrp73n5k1m, thrn=thrn73n5k1m, scatter=false)
	pfft73n5k1m =plot_fft(wvfm73n5k1m,wmus100mus)
	plot(pfw73n5k1m, pfft73n5k1m)
end

md"### Filter signal
- Filter frequency (in Hz) = $ffl73n5k1m

Filtered wvfm:
- mean = $fmean73n5k1m
- std = $fstd73n5k1m
"
begin
    filter_signal_lp!(wvfm73n5k1m, wmus100mus; filtertype="Butterworth", flhz=ffl73n5k1m, n=4)
    first(wvfm73n5k1m,3)
end

begin
    fmean73n5k1m, fstd73n5k1m, fthrp73n5k1m, fthr73n5k1m = thrx(wvfm73n5k1m; nsigma=2.0, waveform="filtered")
    plot_filtered_waveform(wvfm73n5k1m; window=false, thrp=fthrp73n5k1m, thrn=fthr73n5k1m)
end

md"### Select peaks"

sprom73n5k1m = 2.0  #prominence cut
 
md"
- Peak Prominence cut = $sprom73n5k1m
"
begin
	r73n5k1m = select_filtered_peaks(wvfm73n5k1m, fthrp73n5k1m; promsel = sprom73n5k1m, wsel=0.0);
	if r73n5k1m != Nothing
		wvflt73n5k1m, spk73n5k1m = r73n5k1m
		plot_filtered_peaks(wvflt73n5k1m, spk73n5k1m)
	end
end

md"### Optimizing prominence cut"

begin
	P73n5k1m, W73n5k1m = proms(fls73n5k1m, 1, n73n5k1m;  
	                 wmus=wmus100mus, wsym=true, flhz=ffl73n5k1m, promsel=0.5)
	hP73n5k1m = lfi.LaserLab.hist1d(P73n5k1m, 20, 0.0, 10.0)
	hW73n5k1m = lfi.LaserLab.hist1d(W73n5k1m, 20, 5.0, 15.0)
end
lfi.LaserLab.hist1d(hP73n5k1m, "Prominence (promsel=0.5)")
lfi.LaserLab.hist1d(hW73n5k1m, "Width")

md"### Selection"
spks73n5k1m = select_peaks(fls73n5k1m, 1, n73n5k1m; wmus=wmus100mus, wsym=true, flhz=ffl73n5k1m, promsel=sprom73n5k1m)
t73n5k1m = reduce(vcat,[subtract_laser(spks73n5k1m[i].xs, tlus2) for i in 1:length(spks73n5k1m)])
hp73n5k1m = lfi.LaserLab.hist1d(t73n5k1m, 20, 0.0, 2.0)
lfi.LaserLab.hist1d(hp73n5k1m, "times")

md"### Fits"
begin
	tfit73n5k1m = lfi.LaserLab.centers(hp73n5k1m)
	vfit73n5k1m = hp73n5k1m.weights
	fit73n5k1m = curve_fit(mexp, tfit73n5k1m, vfit73n5k1m, pa0)
	cofe73n5k1m = coef(fit73n5k1m)
	stder73n5k1m = stderror(fit73n5k1m)
	tft73n5k1m = expo.(tfit73n5k1m, coef(fit73n5k1m)...)
	println("")
end
begin
    ps73n5k1m = scatter(tfit73n5k1m, vfit73n5k1m, yerr=sqrt.(vfit73n5k1m), markersize=2,
            color = :black,
            label="data",
            fmt = :png)
    pp73n5k1m = plot(ps73n5k1m, tfit73n5k1m, tft73n5k1m, lw=2, label="μ = $(round(cofe73n5k1m[2]*1000, sigdigits=2)) ns", fmt = :png)
    xlabel!("t (μs)")
    ylabel!("frequency")
    end