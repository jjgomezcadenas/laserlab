### A Pluto.jl notebook ###
# v0.18.0

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

# ╔═╡ 52d4f953-98d1-4016-b2b3-7c234b012189
begin
	using Plots
	using Printf
	using Unitful
	using DataFrames
	using CSV
	using Images
	using Interpolations
	using QuadGK
	using UnitfulEquivalences
	using Markdown
	using InteractiveUtils
	using PlutoUI
	using LsqFit
	using Distributions
	using Statistics
	using StatsBase
end

# ╔═╡ 91726783-795a-4da6-9ca7-2d176ba328ca
begin
	using Peaks
	using Glob
	using FFTW
	using DSP
end

# ╔═╡ b4f98dfc-9e44-11ec-3621-4526f6f187bd
md"# RuSL Lifetime: Bold55 200/500 kHz 

- JJGC, 18-03-2022
"

# ╔═╡ 1e429641-5615-4c2f-b466-a5e1d5cbc297
import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
    A, N, mol, mmol, V, L, M

# ╔═╡ 860a66e0-9fdb-4537-865e-bf70e7fa289a
import PhysicalConstants.CODATA2018: N_A

# ╔═╡ 1665009d-f4a8-4d3e-bc9a-3cca8c26c213
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

# ╔═╡ c9e85963-fe90-48c4-8e20-fcc92c31ce99
lfi = ingredients("../src/LaserLab.jl")

# ╔═╡ f3fea484-8430-47ab-986b-7b170bef4ba6
PlutoUI.TableOfContents(title="Table of Contents", indent=true)

# ╔═╡ 9e58009d-0244-49db-801e-362bcb895eb0
md" ## Introduction
Data was taken with the TOPATU setup, with the following elements:

1. EPL480 nm laser operating at 500 kHz (2 μs pulses).
2. A PMT H13543 from Hamamatsu, with 20 % QE at 500 nm, 10.6 % QE at 700 nm and 8 % QE at 800 nm.
3. Data was taken with three different band filters, centered at 660, 700 and 730 nm.
4. Coincidence were triggered in a scope, at 10 MHz rate, the trigger signal was given by the laser.

   **Features of the data**

- The PMT data is somewhat noisy. Thus, it is necessary to characterise single electron pulses by a peak search analysis. 
"

# ╔═╡ aedcf3c2-7d51-466c-a88e-2c681679b8b0
load("../notebooks/img/setupPMTScope.png") 

# ╔═╡ a6cd504a-1afc-4be0-81d2-0de827650526
begin
	spG = lfi.LaserLab.CsvG(';',',')  # spanish: decimals represented with ',' delimited with ';'
	enG = lfi.LaserLab.CsvG(',','.')  # english
	enG2 = lfi.LaserLab.CsvG('\t','.')  # english
	println("")
end

# ╔═╡ 01355b47-de67-4104-a6bc-3d773b156bae
md"# BOLD\_055\_200kHz\_650nm\_50us
"

# ╔═╡ 2e0a2eb7-1b0e-486a-bd0c-a6937a93c350
csvf = glob("../labdata/BOLD_055_200kHz_650nm_50us/C1*.csv") 


# ╔═╡ a097a64b-dcb6-4d52-b88e-c2d2669222ff
w50mus = 50.0e-6

# ╔═╡ 6e35abaf-4392-4882-9531-f6d629373299
n650nm200khz = length(csvf)

# ╔═╡ ceac7188-e300-47b7-b2b4-87c52d3e48d4
tlaser = uconvert(μs, 1.0/200.0kHz)

# ╔═╡ 992eb8da-2d8d-4ec9-b071-d30d83dbea42
tlus = tlaser/μs

# ╔═╡ 2e55dc03-3ca2-4c55-ad36-290b20db5d66
wmus=50.0  # window in μs

# ╔═╡ ca391e17-2f39-4cc3-9487-15449b2d7d69


# ╔═╡ ed36670d-74fa-4ca7-b815-63022739ec48
md"### read waveform, compute stats
- We can also read waveforms specifying a file number in csvg vector (obtained by glob)
"

# ╔═╡ 0bf5a3d6-a66e-458a-9852-93bb4f00d54c
@bind fileNumber NumberField(1:900; default=1)

# ╔═╡ 3e3567bf-70bc-4971-83f6-3614981af6b4
fileNumber

# ╔═╡ b35c0b28-6fd2-4490-9c22-8692a798575a
md"### Filter signal
"

# ╔═╡ 8e17f659-eb82-416c-bd8c-21060a957bc3
@bind flf NumberField(1.0:300.0; default=200.0)

# ╔═╡ 28b3242f-2ddc-4d36-9869-70f57ee98952
ffl = flf*1.0e+6  # in Hz

# ╔═╡ a1a50633-bd5d-4148-8825-c0b0855087a9
md"### Subtract baseline
"

# ╔═╡ 2673dfd5-b088-4f35-8669-29e8f9806fe4
md"### Select peaks
"

# ╔═╡ bf46ed11-6f1c-4f1d-b010-82f6e9988f7c
@bind sprom NumberField(0.0:0.1:10.0; default=0.5)

# ╔═╡ ad1aa93f-ca21-454f-b9ac-10750495113e
sprom

# ╔═╡ 0aacd779-d9c5-4850-8d20-b766496a19e7
typeof(sprom)

# ╔═╡ d9c433bf-4d96-4ee9-a98c-056454db37ea
@bind sw NumberField(0.0:10.0; default=0.0)

# ╔═╡ ecd4080b-38f0-4191-8cd0-d64c6d3674be
md"- Adjust sliders to select on:
- prominence => : $sprom
- Peak width => : $sw

- Set for example the prominence cut at 0.1. One observes that the reflection peaks are not eliminated. This is easily done setting prominence at 0.5
"

# ╔═╡ fa13e453-124b-43c1-84c0-4d0e7ce1ea94
md"### Statistics
- The next step is to study the cuts on a statistical basis. 
"

# ╔═╡ b0a7208d-26c1-4ea5-9072-f43a70d9e1ca
md"- Cut on prominence = $sprom"

# ╔═╡ ab03fda5-fbf7-4ae2-a743-b57e95211dad
md"### Selection"

# ╔═╡ b66b58c1-99f7-40c9-824e-74fdf3da4220
md"### Fits"

# ╔═╡ b3a05a67-3844-4ef8-9bdf-e631477e2b15
md"# BOLD\_055\_500kHz\_650nm\_100us
"

# ╔═╡ e9b35fca-b72a-4b2e-8167-92e070017e25
md"### Read waveform
"

# ╔═╡ 92562f1f-3513-4450-af37-a581a8def5f2
begin
	f500kHz = 500.0kHz
	w100mus = 100.0μs
end

# ╔═╡ f00f105c-3c1f-41f0-a3bc-681af74a1f3d
wmus100mus=w100mus/μs  # window in μs

# ╔═╡ bf287007-369e-4995-a3ad-8d5b96aa1b90
tlaser500kHz = uconvert(μs, 1.0/f500kHz)

# ╔═╡ 569f7f0d-c423-4fce-b8f6-33fc2af07528
tlus2 = tlaser500kHz/μs

# ╔═╡ b4c777b1-1245-4ce9-89f5-d0110c453d60
fls6n5k1m = glob("../labdata/BOLD_055_500kHZ_650nm_100us/C1*.csv") 

# ╔═╡ 938ff2ea-13af-4bce-98d9-8ef5349a7f3d
n6n5k1m = length(fls6n5k1m)

# ╔═╡ d8f37622-4d4c-4983-9184-b7892559cca3
md"
- Files in directory = $n6n5k1m
"

# ╔═╡ 81ac944f-7678-481e-a102-d10cef338ea5
md"- File number to read"

# ╔═╡ 94ea1009-a504-44f0-87a6-6449958bd180
@bind fn6n5k1m NumberField(1:n6n5k1m; default=1)

# ╔═╡ 13fb12a9-5874-4d1d-9d39-f5903be599b6
@bind fl6n5k1m NumberField(1.0:300.0; default=200.0)

# ╔═╡ dc7967a5-9713-4610-ba43-82d24b2481f7
md"- Filter frequency (in MHz) = $fl6n5k1m"

# ╔═╡ ef31cd78-5e84-4f6e-bd17-4b575bb17775
ffl6n5k1m = fl6n5k1m*1.0e+6  # in Hz

# ╔═╡ bd9d7873-544c-4f4c-acaa-22eac4ae7bc1
md"### Filter signal
- Filter frequency (in Hz) = $ffl6n5k1m
"

# ╔═╡ 8d475de0-2018-41a2-9a4c-d7378e9375fd
md"### Select peaks"

# ╔═╡ 34f29042-70f2-4a28-8d61-ded8ee6d6d5a
@bind sprom6n5k1m NumberField(0.0:0.1:10.0; default=0.5)


# ╔═╡ 025a1b05-7831-485a-9cb7-d580f84f7115
md"
- Peak Prominence cut = $sprom6n5k1m
"

# ╔═╡ 4ff4855f-ccbc-4b55-884c-2f35af2f43d5
md"### Optimizing prominence cut"

# ╔═╡ 97256e97-6df1-49cc-a8d3-ebf42041919b
md"### Selection"

# ╔═╡ 50502497-7576-40a8-917c-840e7c4c52f9
md"### Fits"

# ╔═╡ 4a6b9446-8724-4d03-a2d0-48c3d26057b4
md"# BOLD\_055\_500kHZ\_700nm\_100us
"

# ╔═╡ fa381020-606f-4421-9b56-cb5a0360fada
md"### Read waveform"

# ╔═╡ 087d894a-4ab2-4ef9-8824-da29d5755a44
sprom7n5k1m = 2.0

# ╔═╡ cdf477d2-aa0f-464f-9dc1-9aec9aa583ef
md"### Select peaks
- Peak Prominence cut = $sprom7n5k1m
"


# ╔═╡ a7226039-1bb3-4d2e-9a43-4a8588cd7b2d
md"### Optimizing prominence cut"

# ╔═╡ 673b7b87-3e74-40f8-978a-efac0f265b2f
md"### Selection"

# ╔═╡ 133a7d51-791c-4952-8498-b8c938b2fa0d
md"### Fits"

# ╔═╡ e424d154-7090-4311-bb03-ecd06bc10f21
md"# BOLD\_055\_500kHZ\_732nm\_100us"

# ╔═╡ 9c970c92-9cd2-4440-ae40-b41139340c40
md"## Read waveform"

# ╔═╡ 12b0d4b3-75f1-41a5-8b22-5c307105f3e5
fn73n5k1m = 10    # test file to read

# ╔═╡ 1f570b4c-a5e2-40e3-93b0-6be463215d27
md"### Select peaks"

# ╔═╡ e4ba9f6f-e742-4600-864e-50e6794f74d3
sprom73n5k1m = 2.0   # prom cut

# ╔═╡ 88e76067-0b95-4e97-8f3b-bd68b6e0409d
md"
- Peak Prominence cut = $sprom73n5k1m
"

# ╔═╡ 690ada82-db76-4e9b-9795-8f60857f6d94
md"### Optimizing prominence cut"

# ╔═╡ 7a068208-2caf-468a-ae88-df1809cb94d7
md"### Selection"

# ╔═╡ 7e37e56a-f536-49bd-8839-a4af59ff8bf3
md"### Fits"

# ╔═╡ fde49db9-165b-4b55-928b-b0d3f2610a31
md"## Functions and data structures"

# ╔═╡ 500c0b3b-fc21-45a7-b03d-07c413a947c9
md"
	struct SPeaks
"

# ╔═╡ 79c95b20-7b71-414f-b5d6-41e4cfe373ca
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

# ╔═╡ 50e2ec32-be70-444a-9c80-38770ce5da7e
md"
	shift_waveform!(wvfm:DataFrame, wmus::Float64, wsym::Bool)

 Shifts waveform to mV in Y (reverses sign) and to μs in x.
 If the window is centered in zero, shifts it, so that first time is zero.

 **Fields**:
  
*wvfm*  : Waveform (a Data Frame with time/amplitude columns)

*wmus*  : width of time window in μs.

*wsym* : true if the window is symmetric (centered at zero)

**Returns**

Shifts waveform to mV in Y (reverses sign) and to μs in x.
 If the window is centered in zero, shifts it, so that first time is zero.

"

# ╔═╡ 6be64cc8-a357-477c-ae0e-cb0405de4f4b
function shift_waveform!(wvfm::DataFrame, wmus::Float64, wsym::Bool)
	wvfm[!, :Ampl] = -1000.0 * wvfm[!, :Ampl] 
	wvfm[!, :Time] = 1e+6  .* wvfm[!, :Time]
	if wsym
		wvfm[!, :Time] = wmus/2.0  .+ wvfm[!, :Time]
	end
	wvfm
end

# ╔═╡ fe91ac2c-b95a-4d5e-a2fe-f09f43297d4d
md"
	read_waveform(wdir::String, wname::String; wmus::Float64=5, symmetric::Bool=false)
 
**Fields**:
  
*wdir*   : name of directory

*wname*  : name of file
*wmus*  : width of time window in μs.

*wsym* : If true the window is symmetric (centered at zero)

**Returns**

Waveform as a data frame, amplitude in mV (signal positive), time in μs.
"

# ╔═╡ e5bf7e93-f089-4b77-a2e0-d014bdc309e1
function read_waveform(wdir::String, wname::String; wmus::Float64, 
	                   wsym::Bool=false)
	
	wvfm = lfi.LaserLab.load_df_from_csv(wdir, wname, enG,5)
	shift_waveform!(wvfm, wmus, wsym)
	
end

# ╔═╡ 0eded161-9ad3-475e-a6d2-d8068f15bad5
md"
	read_waveform(csvf::Vector{String}, fnumber::Integer; wmus::Float64, 
	              wsym::Bool)

  Like previous method but takes a vector of strings with the name of the files and an integer with the index of the file to be read
"

# ╔═╡ e673cdda-7032-4600-bbd0-01a291aa3eb5
function read_waveform(csvf::Vector{String}, fnumber::Integer; wmus::Float64, 
	                   wsym::Bool)
	
	wvfm = DataFrame(CSV.File(csvf[fnumber], header=5, delim=enG.delim, decimal=enG.decimal))
	shift_waveform!(wvfm, wmus, wsym)
end

# ╔═╡ 8fcae1eb-ffd9-4353-953d-7245559bdd19
wvfm = read_waveform(csvf,fileNumber; wmus=wmus, wsym=true);

# ╔═╡ 8311c814-c9e7-4845-926e-806dc2875dbe
first(wvfm,3)

# ╔═╡ 159234a8-6db6-4006-ab1a-ebea1b39be14
first(wvfm,3)

# ╔═╡ e5033393-655f-405e-978a-e3b904d35b03
wvfm2 = read_waveform(csvf,fileNumber; wmus=wmus, wsym=true);

# ╔═╡ e523ec46-d420-42c2-979a-b7ec4ad2010d
wvfm3 = read_waveform(csvf,fileNumber; wmus=wmus, wsym=true);

# ╔═╡ 9d7c75f9-ffe1-47fe-9bab-784f2a7945c1
wvfm6n5k1m = read_waveform(fls6n5k1m,fn6n5k1m; wmus=wmus100mus, wsym=true);

# ╔═╡ d2a192c8-7039-493b-8c67-d40ac9c366b6
first(wvfm6n5k1m,3)

# ╔═╡ e3413e4a-88ce-47b0-9514-70b1e9245508
begin
	fls7n5k1m = glob("../labdata/BOLD_055_500kHZ_700nm_100us/C1*.csv")
	n7n5k1m = length(fls7n5k1m) 
	fn7n5k1m = 10
	ffl7n5k1m = 200.0e+6
	wvfm7n5k1m = read_waveform(fls7n5k1m,fn7n5k1m; wmus=wmus100mus, wsym=true)
	first(wvfm7n5k1m,3)
end

# ╔═╡ 19257a48-69b7-40a5-859e-483d08928cbf
begin
	# name of file lies, it corresponds to 732 nm
	fls73n5k1m = glob("../labdata/BOLD_055_500kHZ_750nm_100us/C1*.csv")
	n73n5k1m = length(fls73n5k1m) 
	ffl73n5k1m = 200.0e+6
	wvfm73n5k1m = read_waveform(fls73n5k1m,fn73n5k1m; wmus=wmus100mus, wsym=true)
	first(wvfm73n5k1m,3)
end

# ╔═╡ 96dad766-398d-414c-ac66-dfaa76947ffc
md"
		sampling_rate(wvfm::DataFrame, wmus::Float64)
**Fields**:
  
*wvfm*  : Waveform (a Data Frame with time/amplitude columns)

*wmus*  : width of time window in μs.

**Returns**

The sampling rate
"

# ╔═╡ b43916ba-0499-4d82-b81c-49aba21ee145
function sampling_rate(wvfm::DataFrame, wmus::Float64) 
	uconvert(GHz, (length(wvfm.Time) / (wmus*1.0e-6)) * Hz)
end


# ╔═╡ 2432ba2b-fe47-4fb8-bbbb-1b5b971cd96f
sampling_rate(wvfm, wmus)

# ╔═╡ 5dfbb58a-1411-4264-9540-8cc26263910f
md"
		sampling_period(wvfm::DataFrame, wmus::Float64)
  **Fields**:
  
*wvfm*  : Waveform (a Data Frame with time/amplitude columns)

*wmus*  : width of time window in μs.

**Returns**

The sampling period
"

# ╔═╡ ff872f5f-b10b-4065-8850-d81f9d6aece3
function sampling_period(wvfm::DataFrame, wmus::Float64) 
	uconvert(ns, 1.0 / sampling_rate(wvfm, wmus))
end

# ╔═╡ d7e37bdd-a3c3-485b-accd-46c2bd6f3f62
sampling_period(wvfm, wmus)

# ╔═╡ 53849542-08be-4c55-9675-f70bded3546b
begin
	sr6n5k1m = sampling_rate(wvfm6n5k1m,wmus100mus)/GHz
	sp6n5k1m = sampling_period(wvfm6n5k1m,wmus100mus)/ns
end

# ╔═╡ 01333780-fbab-41cf-8cc5-45677ef971e2
md" - Parameters:
 - laser frequency= $f500kHz
 - DAQ window= $w100mus
 - Laser pulses every $tlaser500kHz
 - sampling rate= $(round(sr6n5k1m, sigdigits=2)) GHz
 - sampling period = $(round(sp6n5k1m, sigdigits=2) )ns
   "

# ╔═╡ 0b19c137-824f-4d56-a1ee-a08c4133f602
begin
	sr7n5k1m = sampling_rate(wvfm7n5k1m,wmus100mus)/GHz
	sp7n5k1m = sampling_period(wvfm7n5k1m,wmus100mus)/ns
end

# ╔═╡ b1707796-6aaf-408f-a69d-1414f35d5093
md" - Parameters:
 - laser frequency= $f500kHz
 - DAQ window= $w100mus
 - Laser pulses every $tlaser500kHz
 - sampling rate= $(round(sr7n5k1m, sigdigits=2)) GHz
 - sampling period = $(round(sp7n5k1m, sigdigits=2) )ns
"

# ╔═╡ 3999479a-b980-4e81-a597-f9775896fe56
begin
	sr73n5k1m = sampling_rate(wvfm73n5k1m,wmus100mus)/GHz
	sp73n5k1m = sampling_period(wvfm73n5k1m,wmus100mus)/ns
end

# ╔═╡ a9f1938c-93c8-4e6f-91ad-8b38d0b2f6ff
md" - Parameters:
 - laser frequency= $f500kHz
 - DAQ window= $w100mus
 - Laser pulses every $tlaser500kHz
 - sampling rate= $(round(sr73n5k1m, sigdigits=2)) GHz
 - sampling period = $(round(sp73n5k1m, sigdigits=2) )ns
"

# ╔═╡ 5ff2bb27-3db8-4b05-bd90-c73b493a9c76
md"
	thrx(wvfm::DataFrame; nsigma=1, waveform='filtered')
Return the mean, std, and (+-) threshold at nsigma, for raw or filtered waveform
"

# ╔═╡ f20b33b0-c0ec-4ee7-8273-6b738eb44dd7
function thrx(wvfm::DataFrame; nsigma=1, waveform="filtered")
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

# ╔═╡ 9fa21941-3586-43bc-adeb-92cc1ae805a7
meanAmp, stdAmp, thrp3s, thrn3s = thrx(wvfm; nsigma=3.0, waveform="raw")

# ╔═╡ c31442a6-7480-422f-926d-f97c44996e22
fmeanAmp, fstdAmp, thrp, thrn = thrx(wvfm; nsigma=2.0, waveform="filtered")

# ╔═╡ 3979e3b9-dddc-433f-ab69-6fc1265bdf4f
wvfmFlt = filter(row -> row[:fAmpl] > thrp, wvfm);

# ╔═╡ 99c9d209-18b4-4796-8eed-9786cb13a780
f2meanAmp, f2stdAmp, t2hrp, t2hrn = thrx(wvfm2; nsigma=2.0, waveform="filtered")

# ╔═╡ 20723761-265b-4240-a4a4-132705481f55
f3meanAmp, f3stdAmp, t3hrp, t3hrn = thrx(wvfm3; nsigma=2.0, waveform="filtered")

# ╔═╡ a6aa5ce8-e38f-466b-843e-c98cd3a0b556
mean6n5k1m, std6n5k1m, thrp6n5k1m, thrn6n5k1m = thrx(wvfm6n5k1m; nsigma=3.0, waveform="raw");

# ╔═╡ 017615a3-98ee-4115-9dad-567121079952
md"
Unfiltered wvfm:
- mean = $(round(mean6n5k1m, sigdigits=4))
- std = $(round(std6n5k1m, sigdigits=4))
"

# ╔═╡ 80d9ed58-f4a9-43cd-ab94-251a190f2b3e
fmean6n5k1m, fstd6n5k1m, fthrp6n5k1m, fthr6n5k1m = thrx(wvfm6n5k1m; nsigma=2.0, waveform="filtered")


# ╔═╡ 8c1e081f-c5d3-42e6-8d8b-8f65dddea440
md"
Filtered wvfm:
- mean = $fmean6n5k1m
- std = $fstd6n5k1m
"

# ╔═╡ a564e0f3-f876-4b37-bb33-83fbdafdd513
mean7n5k1m, std7n5k1m, thrp7n5k1m, thrn7n5k1m = thrx(wvfm7n5k1m; nsigma=3.0, waveform="raw");

# ╔═╡ d64991e3-bb7c-4702-a119-f02293a3f686
md"
Unfiltered wvfm:
- mean = $(round(mean7n5k1m, sigdigits=4))
- std = $(round(std7n5k1m, sigdigits=4))
"


# ╔═╡ c548e453-f5bc-40cc-a3d3-b3c0c72d6f9d
md"
	filter_signal_lowpass(wvfm::DataFrame, 
	                      flh::Float64, fs::Float64, filtertype::String, 
	                      n=4, ϵ=1.0)
  **Fields**:
  
*wvfm*  : Waveform (a Data Frame with time/amplitude columns)

*flh*  : higher frequency accepted by the filter

*fs*  : sampling frequency

*filtertype*  : two possible filters: Butterworth or Chebyshev1

*n* : filter order

*ϵ* : ripple factor for Chebyshev1 filter

**Returns**

frequency response, filter (Butterworth, or Chebyshev) and filtered signal
"

# ╔═╡ b5e2dfb7-c325-4898-b45d-8f071f1958b6

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
	H, w = freqresp(filter)	
	H, filter, filt(filter, wvfm.Ampl)
end

# ╔═╡ c5c7a723-6975-4903-ae84-aee4a1bb5441
md"
	filter_signal_bandpass(wvfm::DataFrame, 
	                       fll::Float64, flh::Float64, fs::Float64, 
	                       filtertype::String, n=4, ϵ=1.0)

  **Fields**:
  
*wvfm*  : Waveform (a Data Frame with time/amplitude columns)

*flh*  : lower frequency accepted by the filter

*flh*  : higher frequency accepted by the filter

*fs*  : sampling frequency

*filtertype*  : two possible filters: Butterworth or Chebyshev1

*n* : filter order

*ϵ* : ripple factor for Chebyshev1 filter

**Returns**

frequency response, filter (Butterworth, or Chebyshev) and filtered signal
"

# ╔═╡ 1f950122-3c05-4d4c-85c2-1ed2094b1448

function filter_signal_bandpass(wvfm::DataFrame, 
	                            fll::Float64, flh::Float64, fs::Float64, 
	                            filtertype::String, n=4, ϵ=1.0)

	responsetype = Bandpass(fll, flh; fs=fs)
	
	if filtertype == "Butterworth"
		designmethod = Butterworth(n)
	elseif filtertype == "Chebyshev1"
		ripple = 10.0*log10(ϵ^2 + 1.0)
		designmethod = Chebyshev1(n, ripple)
	end
	filter = digitalfilter(responsetype, designmethod)
	H, w = freqresp(filter)	
	H, filter, filt(filter, wvfm.Ampl)
end

# ╔═╡ d686b4a3-03fd-4213-9b82-45fd395f96ee
md"
	filter_signal_lp!(wvfm::DataFrame, wmus::Float64; 
                           filtertype::String, flhz::Float64, 
                           n=4, ϵ=1.0)
   **Fields**:
  
*wvfm*  : Waveform (a Data Frame with time/amplitude columns)

*wmus*  : width of time window in μs.

*filtertype*  : two possible filters: Butterworth or Chebyshev1

*flhz*  : higer frequency accepted by the filter in Hz

*n* : filter order

*ϵ* : ripple factor for Chebyshev1 filter

**Returns**

A column with the filtered amplitude is added to the waveform
"

# ╔═╡ 35f0fda8-7e9a-4f5c-bbd3-979e8d04ebfa
function filter_signal_lp!(wvfm::DataFrame, wmus::Float64; 
                           filtertype::String, flhz::Float64, 
                           n=4, ϵ=1.0)
	fs = uconvert(Hz, sampling_rate(wvfm, wmus))/Hz
	hcb, fltcb, fsgncb = filter_signal_lowpass(wvfm, flhz, fs, filtertype, n, ϵ)
	wvfm[!, "fAmpl"] = fsgncb
end


# ╔═╡ f184a4be-2a7f-49ef-9d9e-bb8faa087fec
filter_signal_lp!(wvfm, wmus; filtertype="Butterworth", flhz=ffl, n=4)

# ╔═╡ 71c1bbba-8c08-41f5-a1d2-a4de3b6f8890
filter_signal_lp!(wvfm2, wmus; filtertype="Chebyshev1", flhz=ffl, n=4, ϵ=1)

# ╔═╡ 3b567b7b-6fc9-48e8-81ef-90dc6ff4b926
filter_signal_lp!(wvfm6n5k1m, wmus100mus; filtertype="Butterworth", flhz=ffl6n5k1m, n=4)

# ╔═╡ 612f3f29-a64a-4adc-802f-e5492767f933
begin
	filter_signal_lp!(wvfm7n5k1m, wmus100mus; filtertype="Butterworth", flhz=ffl7n5k1m, n=4)
	first(wvfm7n5k1m,3)
end



# ╔═╡ fc9bee5c-26aa-4e91-aaa0-c071e8a63737
begin
    filter_signal_lp!(wvfm73n5k1m, wmus100mus; filtertype="Butterworth", flhz=ffl73n5k1m, n=4)
    first(wvfm73n5k1m,3)
end

# ╔═╡ d8a84dea-adbb-48c5-809b-309236f3dfc5
md"
	filter_signal_bp!(wvfm::DataFrame, wmus::Float64; 
                           filtertype:String, fihz::Float64, flhz::Float64, 
                           n=4, ϵ=1)
"

# ╔═╡ 3d5ef2d8-1d1c-4f90-9178-f72ce4d67f51
function filter_signal_bp!(wvfm::DataFrame, wmus::Float64; 
                           filtertype::String, fihz::Float64, flhz::Float64, 
                           n=4, ϵ=1)
	fs = uconvert(Hz, sampling_rate(wvfm, wmus))/Hz
	hcb, fltcb, fsgncb = filter_signal_bandpass(wvfm, fihz, flhz, fs, 
		                                        filtertype, n, ϵ)
	wvfm[!, "fAmpl"] = fsgncb
end

# ╔═╡ b55852f8-9cbc-4346-a938-b820f04b7ce0
filter_signal_bp!(wvfm3, wmus; 
                  filtertype="Butterworth", fihz=1.0e+6, flhz=290.0e+6, n=4)

# ╔═╡ 647e0043-3389-43f2-bfd6-329d88fce931
md"
	plot_filtered_waveform(wvfm; window=false, waveform='filtered', wstart=-1, wend=-1, thrp=0, thrn=0)
"

# ╔═╡ 843d888f-1e8d-4875-b175-1bc0a98f3320
function plot_filtered_waveform(wvfm; window=false, waveform="filtered", wstart=-1, wend=-1, thrp=0, thrn=0, scatter=true)
	
	if waveform=="filtered"
		ampl = "fAmpl"
	else
		ampl = "Ampl"
	end
	
	if window
		plot(wvfm.Time[wstart:wend], wvfm[wstart:wend, ampl], lw=2, label=waveform, fmt = :png)
		if scatter
			scatter!(wvfm.Time[wstart:wend], wvfm[wstart:wend, ampl],  markersize=2,
					  color = :black, label=false, fmt = :png)
		end
	else
		plot(wvfm.Time, wvfm[!, ampl], lw=2, label=waveform, fmt = :png)
		if scatter
			scatter!(wvfm.Time, wvfm[!, ampl], markersize=2,
					  color = :black, label=false, fmt = :png)
		end
	end
	
	hline!([thrp], label=false, fmt = :png)
	hline!([thrn], label=false, fmt = :png)
	xlabel!("Time (μs or index)")
	ylabel!("Amplitude")
end

# ╔═╡ b3ddc6fe-6752-4f53-8eee-da30555fd95f
plot_filtered_waveform(wvfm; window=false, waveform="raw", thrp=thrp3s, thrn=thrn3s)

# ╔═╡ 218b7cbc-641d-4b55-bcb3-f5572fffaa32
plot_filtered_waveform(wvfm; window=false, thrp=thrp, thrn=thrn)

# ╔═╡ f80a6c49-ad44-420a-bb53-4f196c26e260
plot_filtered_waveform(wvfm2; window=false, thrp=t2hrp, thrn=t2hrn)

# ╔═╡ d5bd5954-4620-4855-87a6-408ce8479552
plot_filtered_waveform(wvfm3; window=false, thrp=t3hrp, thrn=t3hrn)

# ╔═╡ a2f24a67-c112-4b10-94e8-10d078c29181
plot_filtered_waveform(wvfmFlt; window=false, thrp=thrp, thrn=thrn)

# ╔═╡ 9cb01ba4-2a1b-47f4-83db-859bdd6e8178
plot_filtered_waveform(wvfm6n5k1m; window=false, waveform="raw", thrp=thrp6n5k1m, thrn=thrn6n5k1m, scatter=false)

# ╔═╡ 616f1a1c-86dc-41a4-909b-d78ef0663977
plot_filtered_waveform(wvfm6n5k1m; window=false, thrp=fthrp6n5k1m, thrn=fthr6n5k1m)

# ╔═╡ 61dbb1ea-e902-493f-9273-27d7d381ed99
plot_filtered_waveform(wvfm7n5k1m; window=false, waveform="raw", thrp=thrp7n5k1m, thrn=thrn7n5k1m, scatter=false)

# ╔═╡ 16de3662-eac1-4dd8-ae77-de224e935add
begin
	fmean7n5k1m, fstd7n5k1m, fthrp7n5k1m, fthr7n5k1m = thrx(wvfm7n5k1m; nsigma=2.0, waveform="filtered")
	plot_filtered_waveform(wvfm7n5k1m; window=false, thrp=fthrp7n5k1m, thrn=fthr7n5k1m)
end

# ╔═╡ b4aa6c2e-f682-4ae4-873a-dcabd44fc528
md"### Filter signal
- Filter frequency (in Hz) = $ffl7n5k1m

Filtered wvfm:
- mean = $fmean7n5k1m
- std = $fstd7n5k1m
"

# ╔═╡ 3389bf09-0976-4dde-99e2-7500ea80d695
begin
    fmean73n5k1m, fstd73n5k1m, fthrp73n5k1m, fthr73n5k1m = thrx(wvfm73n5k1m; nsigma=2.0, waveform="filtered")
    plot_filtered_waveform(wvfm73n5k1m; window=false, thrp=fthrp73n5k1m, thrn=fthr73n5k1m)
end

# ╔═╡ 78fb4680-4b8a-42f0-a69d-c68cd2917168
md"### Filter signal
- Filter frequency (in Hz) = $ffl73n5k1m

Filtered wvfm:
- mean = $fmean73n5k1m
- std = $fstd73n5k1m
"

# ╔═╡ 99abe607-0a0e-411a-8784-df5318be1db9
md"
	select_peaks(csvf, i0, il; fl=200.0MHz, nsigma=2.0, promsel=0.5, wsel=6.0)
"

# ╔═╡ c126abdc-a067-45fe-be68-817399f31e65
md"
	proms(csvf::Vector{String}, fi::Integer, fe::Integer; 
               wmus::Float64, wsym::Bool, flhz::Float64, promsel::Float64, 
               nsigma=2.0, filtertype='Butterworth')
"

# ╔═╡ 9d34622a-d7c4-4621-96e8-51e9b801a829
md"
	select_filtered_peaks(wvfm::DataFrame, thrp::Float64; 
                               promsel::Number, wsel=0.0)
"

# ╔═╡ 58d706c0-44f5-4514-b87b-eaf1f1940448
function select_filtered_peaks(wvfm::DataFrame, thrp::Float64; 
                               promsel::Number, wsel=0.0)

	wvfmFlt = filter(row -> row[:fAmpl] > thrp, wvfm)
	pks, vals = findmaxima(wvfmFlt.fAmpl)
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

# ╔═╡ 853787cc-229c-4b0b-afe6-4bc36666dec6
result = select_filtered_peaks(wvfm, thrp; promsel = sprom, wsel=sw);

# ╔═╡ 9428c207-98e0-4fcd-8477-d2d92a15d3c9
if result != Nothing
	wvflt, speaks = result
	println("")
end

# ╔═╡ 705652f8-be2b-42b6-b76f-672015f78f19
speaks.proms

# ╔═╡ 6a07f88c-5012-4bd4-b4f9-0cd28e068f12
begin
	r6n5k1m = select_filtered_peaks(wvfm6n5k1m, fthrp6n5k1m; promsel = sprom6n5k1m, wsel=0.0);
	if r6n5k1m != Nothing
		wvflt6n5k1m, spk6n5k1m = r6n5k1m
		println("")
	end
end

# ╔═╡ 0b2d1160-9c7c-4e9d-9545-4c183fb193d2
begin
	r7n5k1m = select_filtered_peaks(wvfm7n5k1m, fthrp7n5k1m; promsel = sprom7n5k1m, wsel=0.0);
	if r7n5k1m != Nothing
		wvflt7n5k1m, spk7n5k1m = r7n5k1m
		println("")
	end
end

# ╔═╡ f6902608-37a4-4c63-98d4-6ad7a704907c
function select_peaks(csvf::Vector{String}, i0::Integer, il::Integer; 
					  wmus::Float64, wsym::Bool, flhz::Float64, promsel::Number, 
               		  nsigma=2.0, filtertype="Butterworth")
	SPK = []
	for fileNumber in i0:il
		wvfm = read_waveform(csvf,fileNumber; wmus=wmus, wsym=wsym)
		filter_signal_lp!(wvfm, wmus; filtertype=filtertype, flhz=flhz, n=4)
		_, _, thrp, _ = thrx(wvfm; nsigma=nsigma, waveform="filtered")
		result = select_filtered_peaks(wvfm, thrp; promsel = promsel, wsel=0.0)
		if result != Nothing
			_, speaks = result
			push!(SPK,speaks)
		end
	end
	SPK
end

# ╔═╡ cea977d9-69a5-4b1a-8e54-c8389d95e2d4
spks = select_peaks(csvf, 1, n650nm200khz; 
					  wmus=wmus, wsym=true, flhz=ffl, promsel=sprom)

# ╔═╡ f8b6d709-1f14-4691-b7a2-93932a972ce0
spks[1].xs

# ╔═╡ d1409c38-01ae-4de2-b4b1-de3ead44c98b
spks6n5k1m = select_peaks(fls6n5k1m, 1, n6n5k1m; 
					  wmus=wmus100mus, wsym=true, flhz=ffl6n5k1m, promsel=sprom6n5k1m)


# ╔═╡ c1c210fe-03ca-4139-91ac-c9687ebc4ca1
spks7n5k1m = select_peaks(fls7n5k1m, 1, n7n5k1m; wmus=wmus100mus, wsym=true, flhz=ffl7n5k1m, promsel=sprom7n5k1m)

# ╔═╡ 00c178f8-d4b9-4822-937a-f4e82da84f11
function proms(csvf::Vector{String}, fi::Integer, fe::Integer; 
               wmus::Float64, wsym::Bool, flhz::Float64, promsel::Number, 
               nsigma=2.0, filtertype="Butterworth")
	
	spks0 = select_peaks(csvf, fi, fe; wmus=wmus, wsym=wsym, flhz=flhz,
	                     promsel=promsel)
	sP = reduce(vcat,[spks0[i].proms for i in 1:length(spks0)])
	sW = reduce(vcat,[spks0[i].widths for i in 1:length(spks0)])
	return sP, sW
end

# ╔═╡ 66bf47c3-d538-40de-9dca-aeb560dfea15
P2P, W2P = proms(csvf, 1, n650nm200khz;  
                 wmus=wmus, wsym=true, flhz=ffl, promsel=0.5)

# ╔═╡ 67cd1b0d-89cc-46ca-ae4b-3900d072ea68
hP = lfi.LaserLab.hist1d(P2P, 20, 0.0, 10.0)

# ╔═╡ 321b373a-5fda-4c9c-b1ad-34d1f82d086b
lfi.LaserLab.hist1d(hP, "Prominence (+)")

# ╔═╡ aad39fc9-433c-4155-92a7-6462971dbf22
hW = lfi.LaserLab.hist1d(W2P, 20, 5.0, 15.0)

# ╔═╡ 015712c4-b2b9-48d4-a076-427a8d0a0baa
lfi.LaserLab.hist1d(hW, "Width (+)")

# ╔═╡ 681425bc-3bc4-4f34-802a-af735d745baf
PP2, WP2 = proms(csvf, 1, n650nm200khz;  
                 wmus=wmus, wsym=true, flhz=ffl, promsel=sprom)

# ╔═╡ 07c0a448-8acd-4387-a755-4bdfa06254df
hP2 = lfi.LaserLab.hist1d(PP2, 20, 0.0, 10.0)

# ╔═╡ f93d9148-1fca-4161-8d27-03d6635da764
lfi.LaserLab.hist1d(hP2, "Prominence (+)")

# ╔═╡ e9144680-0127-4cbd-8408-6fec1f412fc8
hW2 = lfi.LaserLab.hist1d(WP2, 20, 0.0, 20.0)

# ╔═╡ d87eaa68-8c56-48c5-b19f-2f95c70f0e6b
lfi.LaserLab.hist1d(hW2, "Width (+)")

# ╔═╡ 1fe12f82-acce-4a0f-a767-6f33945db454
begin
	P6n5k1m, W6n5k1m = proms(fls6n5k1m, 1, n6n5k1m;  
	                 wmus=wmus100mus, wsym=true, flhz=ffl6n5k1m, promsel=0.5)
	hP6n5k1m = lfi.LaserLab.hist1d(P6n5k1m, 20, 0.0, 10.0)
	hW6n5k1m = lfi.LaserLab.hist1d(W6n5k1m, 20, 5.0, 15.0)
end

# ╔═╡ d8801d90-f74e-4913-8963-75484743303f
lfi.LaserLab.hist1d(hP6n5k1m, "Prominence (promsel=0.5)")

# ╔═╡ 279cc721-78cf-41be-9120-0e2ceea7300d
lfi.LaserLab.hist1d(hW6n5k1m, "Width")

# ╔═╡ 0d49b4d3-6e86-414a-9db4-76d059c2b67c
begin
	P7n5k1m, W7n5k1m = proms(fls7n5k1m, 1, n7n5k1m;  
	                 wmus=wmus100mus, wsym=true, flhz=ffl7n5k1m, promsel=0.5)
	hP7n5k1m = lfi.LaserLab.hist1d(P7n5k1m, 20, 0.0, 10.0)
	hW7n5k1m = lfi.LaserLab.hist1d(W7n5k1m, 20, 5.0, 15.0)
end

# ╔═╡ 7a83f2f4-df2f-4b81-b61c-62095e18d47f
lfi.LaserLab.hist1d(hP7n5k1m, "Prominence (promsel=0.5)")

# ╔═╡ 75a6ea58-dd9b-4e61-b483-de37856e9229
lfi.LaserLab.hist1d(hW7n5k1m, "Width")

# ╔═╡ 0956e274-9339-4314-9baa-7174bdbcde4e
begin
	P73n5k1m, W73n5k1m = proms(fls73n5k1m, 1, n73n5k1m;  
	                 wmus=wmus100mus, wsym=true, flhz=ffl73n5k1m, promsel=0.5)
	hP73n5k1m = lfi.LaserLab.hist1d(P73n5k1m, 20, 0.0, 10.0)
	hW73n5k1m = lfi.LaserLab.hist1d(W73n5k1m, 20, 5.0, 15.0)
end

# ╔═╡ d8528f90-0e76-403e-8131-762be3fef133
lfi.LaserLab.hist1d(hP73n5k1m, "Prominence (promsel=0.5)")

# ╔═╡ 3ead5520-6e2e-4189-806d-a8ad6c8b18e2
lfi.LaserLab.hist1d(hW73n5k1m, "Width")

# ╔═╡ a840300f-0724-44c8-ba08-d7b59491e0e4
md"
	plot_filtered_peaks(wvfmFlt, spks)
"

# ╔═╡ c1894849-5d29-4857-a715-9b91ffe8163d
function plot_filtered_peaks(wvfmFlt, spks)
	plot(wvfmFlt.fAmpl, lw=2, label=false)
	scatter!(spks.peaks, spks.ys, label="peak")
	scatter!(spks.leftedge, spks.xbase, label="leftedge")
	scatter!(spks.rightedge, spks.xbase, label="rightedge")
	hline!([spks.promsel], label=false)
	xlabel!("peak index")
	ylabel!("peak amplitude")
end

# ╔═╡ b49aaeaf-63ea-4bf1-a61b-beca60c8df98
plot_filtered_peaks(wvflt, speaks)

# ╔═╡ 608ba437-7c4e-40cc-a916-eca292a637f2
plot_filtered_peaks(wvflt6n5k1m, spk6n5k1m)

# ╔═╡ 4a435210-a6de-4c61-8ef8-46fba6931b0f
plot_filtered_peaks(wvflt7n5k1m, spk7n5k1m)

# ╔═╡ 00b257ad-c634-4f01-b195-4bddb0c074fe
begin
	r73n5k1m = select_filtered_peaks(wvfm73n5k1m, fthrp73n5k1m; promsel = sprom73n5k1m, wsel=0.0);
	if r73n5k1m != Nothing
		wvflt73n5k1m, spk73n5k1m = r73n5k1m
		plot_filtered_peaks(wvflt73n5k1m, spk73n5k1m)
	end
end

# ╔═╡ 083a5ba2-7de3-4e51-962f-08289fe987ce
md"
	function subtract_laser(tv::Vector{Float64})
"

# ╔═╡ 773563de-8eb2-4ea4-87fc-86609e35370f
function subtract_laser(tv::Vector{Float64}, ltus::Float64)
	tv - ltus * floor.(tv/ltus)
end

# ╔═╡ a2d9b266-a9ca-4c88-9dbc-26815fafac21
subtract_laser(speaks.xs, tlus)

# ╔═╡ ab7c483e-f8db-4556-9d78-8a87a7c2cf77
times = reduce(vcat,[subtract_laser(spks[i].xs, tlus) for i in 1:length(spks)])

# ╔═╡ fee7bf69-da15-484b-82ce-de0246a49b75
hpkp = lfi.LaserLab.hist1d(times, 20, 0.0, 5.0)

# ╔═╡ 6d5bc9d2-3c7a-497b-9372-e4d43eabf3c4
lfi.LaserLab.hist1d(hpkp, "Peaks (+)")

# ╔═╡ 8065bf42-8d71-4c99-85c6-a2b6fd2e2087
begin
	expo(t, N, λ) = N*exp(-t/λ)
	mexp(t, p) = p[1] * exp.(-t/p[2])
	pa0 = [100.0, 0.5]
	tdata = lfi.LaserLab.centers(hpkp)
	vdata = hpkp.weights
	fit = curve_fit(mexp, tdata, vdata, pa0)
	cofe = coef(fit)
	stder = stderror(fit)
	tft = expo.(tdata, coef(fit)...)
	println("")
end

# ╔═╡ 73878b74-389a-4306-834b-9d9a126ba95d
cofe

# ╔═╡ d73cae53-b429-4c06-a63b-9bc876387256
stder

# ╔═╡ baba5028-d85c-415f-95f6-0393ce686249
begin
ps1 = scatter(tdata, vdata, yerr=sqrt.(vdata), markersize=2,
		color = :black,
	    label="data",
		fmt = :png)
pp = plot(ps1, tdata, tft, lw=2, label="μ = $(round(cofe[2]*1000, sigdigits=2)) ns", fmt = :png)
xlabel!("t (μs)")
ylabel!("frequency")
end

# ╔═╡ fd31ac5a-7248-41fc-b6cc-0563790f07b6
t6n5k1m = reduce(vcat,[subtract_laser(spks6n5k1m[i].xs, tlus2) for i in 1:length(spks6n5k1m)])


# ╔═╡ a09a9d75-79e7-4f6a-a058-0309901838c4
hp6n5k1m = lfi.LaserLab.hist1d(t6n5k1m, 20, 0.0, 2.0)

# ╔═╡ 965439d1-698f-47ad-be87-04d7d27e42fc
lfi.LaserLab.hist1d(hp6n5k1m, "times")

# ╔═╡ 9637e720-24f1-4321-bb77-72316dfb7bbd
begin
	tfit6n5k1m = lfi.LaserLab.centers(hp6n5k1m)
	vfit6n5k1m = hp6n5k1m.weights
	fit6n5k1m = curve_fit(mexp, tfit6n5k1m, vfit6n5k1m, pa0)
	cofe6n5k1m = coef(fit6n5k1m)
	stder6n5k1m = stderror(fit6n5k1m)
	tft6n5k1m = expo.(tfit6n5k1m, coef(fit6n5k1m)...)
	println("")
end

# ╔═╡ 57008d1a-1b6d-42a5-86bb-62a9bd8c2f63
begin
ps6n5k1m = scatter(tfit6n5k1m, vfit6n5k1m, yerr=sqrt.(vfit6n5k1m), markersize=2,
		color = :black,
	    label="data",
		fmt = :png)
pp6n5k1m = plot(ps6n5k1m, tfit6n5k1m, tft6n5k1m, lw=2, label="μ = $(round(cofe6n5k1m[2]*1000, sigdigits=2)) ns", fmt = :png)
xlabel!("t (μs)")
ylabel!("frequency")
end

# ╔═╡ 286c10f9-b1e5-4088-926c-1b86c11107a3
begin
	t7n5k1m = reduce(vcat,[subtract_laser(spks7n5k1m[i].xs, tlus2) for i in 1:length(spks7n5k1m)])
	hp7n5k1m = lfi.LaserLab.hist1d(t7n5k1m, 20, 0.0, 2.0)
	lfi.LaserLab.hist1d(hp7n5k1m, "times")
end

# ╔═╡ e2914357-c8c2-4fbf-a0f9-308477242b73
begin
	tfit7n5k1m = lfi.LaserLab.centers(hp7n5k1m)
	vfit7n5k1m = hp7n5k1m.weights
	fit7n5k1m = curve_fit(mexp, tfit7n5k1m, vfit7n5k1m, pa0)
	cofe7n5k1m = coef(fit7n5k1m)
	stder7n5k1m = stderror(fit7n5k1m)
	tft7n5k1m = expo.(tfit7n5k1m, coef(fit7n5k1m)...)
	println("")
end

# ╔═╡ 7e9fa3f9-a1a9-45c3-a0cc-1591d320f8df
begin
    ps7n5k1m = scatter(tfit7n5k1m, vfit7n5k1m, yerr=sqrt.(vfit7n5k1m), markersize=2,
            color = :black,
            label="data",
            fmt = :png)
    pp7n5k1m = plot(ps7n5k1m, tfit7n5k1m, tft7n5k1m, lw=2, label="μ = $(round(cofe7n5k1m[2]*1000, sigdigits=2)) ns", fmt = :png)
    xlabel!("t (μs)")
    ylabel!("frequency")
end

# ╔═╡ 718a44cf-d893-4574-95d8-7e9c23f47fbc
begin
    spks73n5k1m = select_peaks(fls73n5k1m, 1, n73n5k1m; wmus=wmus100mus, wsym=true, flhz=ffl73n5k1m, promsel=sprom73n5k1m)
    t73n5k1m = reduce(vcat,[subtract_laser(spks73n5k1m[i].xs, tlus2) for i in 1:length(spks73n5k1m)])
    hp73n5k1m = lfi.LaserLab.hist1d(t73n5k1m, 20, 0.0, 2.0)
end

# ╔═╡ 93e82b00-94d9-4691-8e36-653b4f2c3526
lfi.LaserLab.hist1d(hp73n5k1m, "times")

# ╔═╡ 3744cb91-f3f4-481f-964f-40574066fb86
begin
	tfit73n5k1m = lfi.LaserLab.centers(hp73n5k1m)
	vfit73n5k1m = hp73n5k1m.weights
	fit73n5k1m = curve_fit(mexp, tfit73n5k1m, vfit73n5k1m, pa0)
	cofe73n5k1m = coef(fit73n5k1m)
	stder73n5k1m = stderror(fit73n5k1m)
	tft73n5k1m = expo.(tfit73n5k1m, coef(fit73n5k1m)...)
	println("")
end

# ╔═╡ bfb83129-5f50-40fe-87d5-793cfd292bc7
begin
    ps73n5k1m = scatter(tfit73n5k1m, vfit73n5k1m, yerr=sqrt.(vfit73n5k1m), markersize=2,
            color = :black,
            label="data",
            fmt = :png)
    pp73n5k1m = plot(ps73n5k1m, tfit73n5k1m, tft73n5k1m, lw=2, label="μ = $(round(cofe73n5k1m[2]*1000, sigdigits=2)) ns", fmt = :png)
    xlabel!("t (μs)")
    ylabel!("frequency")
end

# ╔═╡ cfaf2b26-e39d-47b5-8e99-8d35dbd7b3ac
md"
		fourier_transform(wvfm1, ws)
  	
   Given waveform **wvfm1** and the data window **ws** in seconds, returns the sampling period, fourier transform amplitude and frequency
"

# ╔═╡ 13e729c1-15ff-4fcf-b3df-9247fda54109
function fourier_transform(wvfm1, ws)
	ts = sampling_period(wvfm1, ws)
	tss = uconvert(s, ts)/s           # tss is the sampling period in seconds
	F = fft(wvfm1.Ampl) |> fftshift   # fourier transform
	freqs = fftfreq(length(wvfm1.Time), 1.0/tss) |> fftshift  # frequencies in Hz
	ts, F, freqs
end

# ╔═╡ 52e6ee88-f49d-420e-9a99-9c8c54910362
md"
	plot_fft(wvfm1, ws; window=false, is=-1, il=-1)
"

# ╔═╡ 604081ba-c975-414a-a8b6-68a6c401172f
function plot_fft(wvfm1, ws; window=false, is=-1, il=-1)
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

# ╔═╡ 34b066df-4f23-4bc8-af26-1a353bfedaf6
plot_fft(wvfm, wmus)

# ╔═╡ fbeb7ba2-00c4-4d52-a0a3-cd39c38de58b
plot_fft(wvfm, wmus; window=true, is=1, il=18000)

# ╔═╡ 7bb296f1-6f93-4007-8503-4e5e85d26751
plot_fft(wvfm6n5k1m,wmus100mus)

# ╔═╡ bc862bfe-6b3b-444f-8951-653b204ebbf6
plot_fft(wvfm7n5k1m,wmus100mus)

# ╔═╡ b31b5d94-da48-4e84-8c8f-10a5185b32e9
begin
	mean73n5k1m, std73n5k1m, thrp73n5k1m, thrn73n5k1m = thrx(wvfm73n5k1m; nsigma=3.0, waveform="raw");
	pfw73n5k1m = plot_filtered_waveform(wvfm73n5k1m; window=false, waveform="raw", thrp=thrp73n5k1m, thrn=thrn73n5k1m, scatter=false)
	pfft73n5k1m =plot_fft(wvfm73n5k1m,wmus100mus)
	plot(pfw73n5k1m, pfft73n5k1m)
end

# ╔═╡ 0b083e48-8c7a-4f16-bf28-b5470f23acdb
md"
Unfiltered wvfm:
- mean = $(round(mean73n5k1m, sigdigits=4))
- std = $(round(std73n5k1m, sigdigits=4))
"

# ╔═╡ 1928ef8d-d7df-4b2b-9b9a-f958cacecaa3
md"
	plot_waveform(wvfm, mean, std, nsigma, nsigma2; sctter=false)
"

# ╔═╡ 17fffef0-699e-4aa2-8931-41afa2ef4da1
function plot_waveform(wvfm, mean, std, nsigma, nsigma2; sct=false, trace=false)
	if sct
		p1 = scatter(wvfm.Time, wvfm.Ampl,  markersize=2,
				  color = :black,
	    		  label=false,
				  fmt = :png)
	elseif trace
		p1 = plot(wvfm.Ampl, lw=2, label=false, fmt = :png)
	else
		p1 = plot(wvfm.Time, wvfm.Ampl, lw=2, label=false, fmt = :png)
	end
	hline!([mean - nsigma * std], label=false, fmt = :png)
	hline!([mean + nsigma * std], label=false, fmt = :png)
	hline!([mean - nsigma2 * std], label=false, fmt = :png)
	hline!([mean + nsigma2 * std], label=false, fmt = :png)
	xlabel!("t (μs)")
	ylabel!("I (mV)")
	p1
end

# ╔═╡ 0bb31d5b-6d5c-49cc-9ddc-1675dfe4c98f
md"
	function mean_and_std(wvfm)
"

# ╔═╡ 8fe40af8-5374-48ed-831e-69b52a9e3ffc
function mean_and_std(wvfm)
	mean(wvfm.Ampl), std(wvfm.Ampl)
end

# ╔═╡ f3b5c2b0-1039-4edf-9efe-6da3665e5b39



# ╔═╡ 84c75a93-92ce-4df4-be33-9fd71bffb2ca


# ╔═╡ b0738b82-5ce0-46de-828b-5c25477d2449


# ╔═╡ 70f3237c-c051-494a-91ef-ad97bcc07e94


# ╔═╡ c056af02-693b-4837-ab7b-1c4c7b7a9ca5
md"### glob_files
	function glob_files(path)
"

# ╔═╡ 8d02c42c-452e-4ad3-9b27-95ac7ea5e210
function glob_files(path)
	sdr = string(path, "/rep*")
	drs = glob(sdr) 
	FLS = String[]
	for dr in drs
		fls   = string(dr, "/C1*.csv")
		csf = glob(fls) 
		append!(FLS, csf)
	end
	FLS
end

# ╔═╡ 380a7a82-0924-42b6-8698-7573b389e9ad
function glob_files2(path; wvl="650nm")
	sdr = string(path, "/", wvl,"*")
	drs = glob(sdr) 
	FLS = String[]
	for dr in drs
		fls   = string(dr, "/C1*.csv")
		csf = glob(fls) 
		append!(FLS, csf)
	end
	FLS
end

# ╔═╡ f62e8d24-10c6-42e9-a68e-f2d68f950ed1
md"### fit_peaks
	fit_peaks(htime; pa0=[100.0, 0.5], i0=1)
"

# ╔═╡ e45bd546-c5fc-4914-b237-3df5f473363d
function fit_peaks(htime; pa0=[100.0, 0.5], i0=1)
	expo(t, N, λ) = N*exp(-t/λ)
	mexp(t, p) = p[1] * exp.(-t/p[2])
	tdata = lfi.LaserLab.centers(htime)
	vdata = htime.weights
	il = length(tdata)
	fit = curve_fit(mexp, tdata[i0:il], vdata[i0:il], pa0)
	coef(fit), stderror(fit), expo.(tdata, coef(fit)...)
end

# ╔═╡ 5e6ea853-c7a2-44ce-acb5-628ed36f8180
md"### plot_fit
	plot_fit(htime, coeff, tft)
"

# ╔═╡ 8953d6fb-f9ad-4f95-b022-5ac95058dbf7
function plot_fit(htime, coeff, tft; savefig=false, fn="")
	tdata =lfi.LaserLab.centers(htime)
	vdata =htime.weights
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

# ╔═╡ 2174c692-ccc7-490b-b2c4-940876c579dd
begin
    ppx2 = plot_fit(hpkp, cofe, tft; savefig=true,
	              fn="b55A_200kHz_650nm_50mus_fit.png")
	plot(ppx2)
end

# ╔═╡ 5dd92589-d1a9-4eb0-9e58-d173b8ff3a04
begin
    ppx3 = plot_fit(hp6n5k1m, cofe6n5k1m, tft6n5k1m; savefig=true,
	              fn="b55A_500kHz_650nm_100mus_fit.png")
	plot(ppx3)
end

# ╔═╡ b375e4a5-5e99-47da-841b-2475936c437d
begin
    ppx4 = plot_fit(hp7n5k1m, cofe7n5k1m, tft7n5k1m; savefig=true,
	              fn="b55A_500kHz_700nm_100mus_fit.png")
	plot(ppx4)
end

# ╔═╡ c8acd563-df8e-4404-904f-6e0b97a13f55
begin
    ppx5 = plot_fit(hp73n5k1m, cofe73n5k1m, tft73n5k1m; savefig=true,
	              fn="b55A_500kHz_732nm_100mus_fit.png")
	plot(ppx5)
end

# ╔═╡ 6c95fd00-0023-4f00-b495-050ca098cfb6


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
Glob = "c27321d9-0574-5035-807b-f59d2c89b15c"
Images = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
InteractiveUtils = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
LsqFit = "2fda8390-95c7-5789-9bda-21331edee243"
Markdown = "d6f4376e-aef5-505a-96c1-9c027394607a"
Peaks = "18e31ff7-3703-566c-8e60-38913d67486b"
PhysicalConstants = "5ad8b20f-a522-5ce9-bfc9-ddf1d5bda6ab"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
QuadGK = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
UnitfulEquivalences = "da9c4bc3-91c8-4f02-8a40-6b990d2a7e0c"

[compat]
CSV = "~0.10.2"
DSP = "~0.7.5"
DataFrames = "~1.3.2"
Distributions = "~0.25.49"
FFTW = "~1.4.6"
Glob = "~1.3.0"
Images = "~0.25.1"
Interpolations = "~0.13.5"
LsqFit = "~0.12.1"
Peaks = "~0.4.0"
PhysicalConstants = "~0.2.1"
Plots = "~1.26.0"
PlutoUI = "~0.7.35"
QuadGK = "~2.4.2"
StatsBase = "~0.33.16"
Unitful = "~1.11.0"
UnitfulEquivalences = "~0.2.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "9f8186bc19cd1c129d367cb667215517cc03e144"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "5.0.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "d127d5e4d86c7680b20c35d40b503c74b9a39b5e"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.4"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "9519274b50500b8029973d241d32cfbf0b127d97"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.2"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c9a6160317d1abe9c44b3beb367fd448117679ca"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.13.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "75479b7df4167267d75294d14b58244695beb2ac"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.14.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "3f1f500312161f1ae067abe07d13b40f78f32e07"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.8"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.CommonSolve]]
git-tree-sha1 = "68a0743f578349ada8bc911a5cbd5a2ef6ed6d1f"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.0"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "681ea870b918e7cff7111da58791d7f718067a19"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

[[deps.DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "3e03979d16275ed5d9078d50327332c546e24e68"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.5"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "ae02104e835f219b8930c7664b8012c93475c340"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.2"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "dd933c4ef7b4c270aacd4eb88fa64c147492acf0"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.10.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "9d3c0c762d4666db9187f363a76b47f7346e673b"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.49"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "84f04fe68a3176a583b864e492578b9466d87f1e"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.EllipsisNotation]]
git-tree-sha1 = "18ee049accec8763be17a933737c1dd0fdf8673a"
uuid = "da5c29d0-fa7d-589e-88eb-ea29b0a81949"
version = "1.0.0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ae13fcbc7ab8f16b0856729b050ef0c446aa3492"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.4+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "cbdf14d1e8c7c8aacbe8b19862e0179fd08321c2"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "505876577b5481e50d089c1c68899dfb6faebc62"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.6"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "80ced645013a5dbdc52cf70329399c35ce007fae"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.13.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "04d13bfa8ef11720c24e4d840c0033d145537df7"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.17"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "4c7d3757f3ecbcb9055870351078552b7d1dbd2d"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.0"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "56956d1e4c1221000b7781104c58c34019792951"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "1bd6fc0c344fc0cbee1f42f8d2e7ec8253dda2d2"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.25"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "9f836fb62492f4b0f0d3b06f55983f2704ed0883"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.0"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a6c850d77ad5118ad3be4bd188919ce97fffac47"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.0+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78e2c69783c9753a91cdae88a8d432be85a2ab5e"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Glob]]
git-tree-sha1 = "4df9f7e06108728ebf00a0a11edee4b29a482bb2"
uuid = "c27321d9-0574-5035-807b-f59d2c89b15c"
version = "1.3.0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "1c5a84319923bea76fa145d49e93aa4394c73fc2"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.1"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "57c021de207e234108a6f1454003120a1bf350c4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.6.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "SpecialFunctions", "Test"]
git-tree-sha1 = "65e4589030ef3c44d3b90bdc5aac462b4bb05567"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.8"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "c54b581a83008dc7f292e205f4c409ab5caa0f04"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.10"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageContrastAdjustment]]
deps = ["ImageCore", "ImageTransformations", "Parameters"]
git-tree-sha1 = "0d75cafa80cf22026cea21a8e6cf965295003edc"
uuid = "f332f351-ec65-5f6a-b3d1-319c6670881a"
version = "0.3.10"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "9a5c62f231e5bba35695a20988fc7cd6de7eeb5a"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.3"

[[deps.ImageDistances]]
deps = ["Distances", "ImageCore", "ImageMorphology", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "7a20463713d239a19cbad3f6991e404aca876bda"
uuid = "51556ac3-7006-55f5-8cb3-34580c88182d"
version = "0.2.15"

[[deps.ImageFiltering]]
deps = ["CatIndices", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageBase", "ImageCore", "LinearAlgebra", "OffsetArrays", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "15bd05c1c0d5dbb32a9a3d7e0ad2d50dd6167189"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.7.1"

[[deps.ImageIO]]
deps = ["FileIO", "JpegTurbo", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "464bdef044df52e6436f8c018bea2d48c40bb27b"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.1"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils", "Libdl", "Pkg", "Random"]
git-tree-sha1 = "5bc1cb62e0c5f1005868358db0692c994c3a13c6"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.2.1"

[[deps.ImageMagick_jll]]
deps = ["Artifacts", "Ghostscript_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f025b79883f361fa1bd80ad132773161d231fd9f"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.12+2"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "36cbaebed194b292590cba2593da27b34763804a"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.8"

[[deps.ImageMorphology]]
deps = ["ImageCore", "LinearAlgebra", "Requires", "TiledIteration"]
git-tree-sha1 = "7668b123ecfd39a6ae3fc31c532b588999bdc166"
uuid = "787d08f9-d448-5407-9aad-5290dd7ab264"
version = "0.3.1"

[[deps.ImageQualityIndexes]]
deps = ["ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "OffsetArrays", "Statistics"]
git-tree-sha1 = "1d2d73b14198d10f7f12bf7f8481fd4b3ff5cd61"
uuid = "2996bd0c-7a13-11e9-2da2-2f5ce47296a9"
version = "0.3.0"

[[deps.ImageSegmentation]]
deps = ["Clustering", "DataStructures", "Distances", "Graphs", "ImageCore", "ImageFiltering", "ImageMorphology", "LinearAlgebra", "MetaGraphs", "RegionTrees", "SimpleWeightedGraphs", "StaticArrays", "Statistics"]
git-tree-sha1 = "36832067ea220818d105d718527d6ed02385bf22"
uuid = "80713f31-8817-5129-9cf8-209ff8fb23e1"
version = "1.7.0"

[[deps.ImageShow]]
deps = ["Base64", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "d0ac64c9bee0aed6fdbb2bc0e5dfa9a3a78e3acc"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.3"

[[deps.ImageTransformations]]
deps = ["AxisAlgorithms", "ColorVectorSpace", "CoordinateTransformations", "ImageBase", "ImageCore", "Interpolations", "OffsetArrays", "Rotations", "StaticArrays"]
git-tree-sha1 = "42fe8de1fe1f80dab37a39d391b6301f7aeaa7b8"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.9.4"

[[deps.Images]]
deps = ["Base64", "FileIO", "Graphics", "ImageAxes", "ImageBase", "ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "ImageIO", "ImageMagick", "ImageMetadata", "ImageMorphology", "ImageQualityIndexes", "ImageSegmentation", "ImageShow", "ImageTransformations", "IndirectArrays", "IntegralArrays", "Random", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "TiledIteration"]
git-tree-sha1 = "11d268adba1869067620659e7cdf07f5e54b6c76"
uuid = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
version = "0.25.1"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "61feba885fac3a407465726d0c330b3055df897f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.2"

[[deps.IntegralArrays]]
deps = ["ColorTypes", "FixedPointNumbers", "IntervalSets"]
git-tree-sha1 = "cf737764159c66b95cdbf5c10484929b247fecfe"
uuid = "1d092043-8f09-5a30-832f-7509e371ab51"
version = "0.1.3"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b15fc0a95c564ca2e0a7ae12c1f095ca848ceb31"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.5"

[[deps.IntervalSets]]
deps = ["Dates", "EllipsisNotation", "Statistics"]
git-tree-sha1 = "3cc368af3f110a767ac786560045dceddfc16758"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.5.3"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "28b114b3279cdbac9a61c57b3e6548a572142b34"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.21"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "a77b273f1ddec645d1b7c4fd5fb98c8f90ad10a5"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a6552bfeab40de157a297d84e03ade4b8177677f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.12"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LsqFit]]
deps = ["Distributions", "ForwardDiff", "LinearAlgebra", "NLSolversBase", "OptimBase", "Random", "StatsBase"]
git-tree-sha1 = "91aa1442e63a77f101aff01dec5a821a17f43922"
uuid = "2fda8390-95c7-5789-9bda-21331edee243"
version = "0.12.1"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "e595b205efd49508358f7dc670a940c790204629"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.0.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf", "RecipesBase", "Requires"]
git-tree-sha1 = "88cd033eb781c698e75ae0b680e5cef1553f0856"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.7.1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.MetaGraphs]]
deps = ["Graphs", "JLD2", "Random"]
git-tree-sha1 = "2af69ff3c024d13bde52b34a2a7d6887d4e7b438"
uuid = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
version = "0.7.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "b34e3bc3ca7c94914418637cb10cc4d1d80d877d"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.3"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "ba8c0f8732a24facba709388c74ba99dcbfdda1e"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.0.0"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "16baacfdc8758bc374882566c9187e785e85c2f0"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.9"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore"]
git-tree-sha1 = "18efc06f6ec36a8b801b23f076e3c6ac7c3bf153"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "648107615c15d4e09f7eca16307bc821c1f718d8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.13+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OptimBase]]
deps = ["NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "9cb1fee807b599b5f803809e85c81b582d2009d6"
uuid = "87e2bd06-a317-5318-96d9-3ecbac512eee"
version = "2.0.2"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "7e2166042d1698b6072352c74cfd1fca2a968253"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.6"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "eb4dbb8139f6125471aa3da98fb70f02dc58e49c"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.14"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "13468f237353112a01b2d6b32f3d0f80219944aa"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.2"

[[deps.Peaks]]
deps = ["Compat"]
git-tree-sha1 = "79e1f108ef46e9393bc670440c5f3ec78d23eb78"
uuid = "18e31ff7-3703-566c-8e60-38913d67486b"
version = "0.4.0"

[[deps.PhysicalConstants]]
deps = ["Measurements", "Roots", "Unitful"]
git-tree-sha1 = "2bc26b693b5cbc823c54b33ea88a9209d27e2db7"
uuid = "5ad8b20f-a522-5ce9-bfc9-ddf1d5bda6ab"
version = "0.2.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "6f1b25e8ea06279b5689263cc538f51331d7ca17"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "23d109aad5d225e945c813c6ebef79104beda955"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.26.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "85bf3e4bd279e405f91489ce518dedb1e32119cb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.35"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "MutableArithmetics", "RecipesBase"]
git-tree-sha1 = "0107e2f7f90cc7f756fee8a304987c574bbd7583"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "3.0.0"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "db3a23166af8aebf4db5ef87ac5b00d36eb771e2"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "de893592a221142f3db370f48290e3a2ef39998f"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.4"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "afadeba63d90ff223a6a48d2009434ecee2ec9e8"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.1"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.Quaternions]]
deps = ["DualNumbers", "LinearAlgebra", "Random"]
git-tree-sha1 = "d0baaa6bcbac4369f1ecfb4a8c44b96ef3e5acb9"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.5.1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "995a812c6f7edea7527bb570f0ac39d0fb15663c"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.1"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegionTrees]]
deps = ["IterTools", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "4618ed0da7a251c7f92e869ae1a19c74a7d2a7f9"
uuid = "dee08c22-ab7f-5625-9660-a9af2021b33f"
version = "0.3.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.Roots]]
deps = ["CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "0abe7fc220977da88ad86d339335a4517944fea2"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "1.3.14"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "a167638e2cbd8ac41f9cd57282cab9b042fa26e6"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "6a2f7d70512d205ca8c7ee31bfa9f142fe74310c"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.12"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "38d88503f695eb0301479bc9b0d4320b378bafe5"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays", "Test"]
git-tree-sha1 = "a6f404cc44d3d3b28c793ec0eb59af709d827e4e"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.2.1"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "87e9954dfa33fd145694e42337bdd3d5b07021a6"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.6.0"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "74fb527333e72ada2dd9ef77d98e4991fb185f04"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.1"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c3d8ba7f3fa0625b062b82853a7d5229cb728b6b"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.1"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "25405d7016a47cf2bd6cd91e66f4de437fd54a07"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.16"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "991d34bbff0d9125d93ba15887d6594e8e84b305"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.5.3"

[[deps.TiledIteration]]
deps = ["OffsetArrays"]
git-tree-sha1 = "5683455224ba92ef59db72d10690690f4a8dc297"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.3.1"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "b649200e887a487468b71821e2644382699f1b0f"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.11.0"

[[deps.UnitfulEquivalences]]
deps = ["Unitful"]
git-tree-sha1 = "76fc2f7fdc87531a1018eb7d647df7c29daf36b7"
uuid = "da9c4bc3-91c8-4f02-8a40-6b990d2a7e0c"
version = "0.2.0"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "c69f9da3ff2f4f02e811c3323c22e5dfcb584cfa"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.1"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78736dab31ae7a53540a6b752efc61f77b304c5b"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.8.6+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╠═b4f98dfc-9e44-11ec-3621-4526f6f187bd
# ╠═52d4f953-98d1-4016-b2b3-7c234b012189
# ╠═91726783-795a-4da6-9ca7-2d176ba328ca
# ╠═1e429641-5615-4c2f-b466-a5e1d5cbc297
# ╠═860a66e0-9fdb-4537-865e-bf70e7fa289a
# ╟─1665009d-f4a8-4d3e-bc9a-3cca8c26c213
# ╠═c9e85963-fe90-48c4-8e20-fcc92c31ce99
# ╠═f3fea484-8430-47ab-986b-7b170bef4ba6
# ╠═9e58009d-0244-49db-801e-362bcb895eb0
# ╟─aedcf3c2-7d51-466c-a88e-2c681679b8b0
# ╠═a6cd504a-1afc-4be0-81d2-0de827650526
# ╠═01355b47-de67-4104-a6bc-3d773b156bae
# ╠═2e0a2eb7-1b0e-486a-bd0c-a6937a93c350
# ╠═a097a64b-dcb6-4d52-b88e-c2d2669222ff
# ╠═6e35abaf-4392-4882-9531-f6d629373299
# ╠═ceac7188-e300-47b7-b2b4-87c52d3e48d4
# ╠═992eb8da-2d8d-4ec9-b071-d30d83dbea42
# ╠═2e55dc03-3ca2-4c55-ad36-290b20db5d66
# ╠═ca391e17-2f39-4cc3-9487-15449b2d7d69
# ╠═ed36670d-74fa-4ca7-b815-63022739ec48
# ╠═0bf5a3d6-a66e-458a-9852-93bb4f00d54c
# ╠═3e3567bf-70bc-4971-83f6-3614981af6b4
# ╠═8fcae1eb-ffd9-4353-953d-7245559bdd19
# ╠═8311c814-c9e7-4845-926e-806dc2875dbe
# ╠═9fa21941-3586-43bc-adeb-92cc1ae805a7
# ╠═b3ddc6fe-6752-4f53-8eee-da30555fd95f
# ╠═2432ba2b-fe47-4fb8-bbbb-1b5b971cd96f
# ╠═d7e37bdd-a3c3-485b-accd-46c2bd6f3f62
# ╠═34b066df-4f23-4bc8-af26-1a353bfedaf6
# ╠═fbeb7ba2-00c4-4d52-a0a3-cd39c38de58b
# ╠═b35c0b28-6fd2-4490-9c22-8692a798575a
# ╠═8e17f659-eb82-416c-bd8c-21060a957bc3
# ╠═28b3242f-2ddc-4d36-9869-70f57ee98952
# ╠═f184a4be-2a7f-49ef-9d9e-bb8faa087fec
# ╠═159234a8-6db6-4006-ab1a-ebea1b39be14
# ╠═c31442a6-7480-422f-926d-f97c44996e22
# ╠═218b7cbc-641d-4b55-bcb3-f5572fffaa32
# ╠═e5033393-655f-405e-978a-e3b904d35b03
# ╠═71c1bbba-8c08-41f5-a1d2-a4de3b6f8890
# ╠═99c9d209-18b4-4796-8eed-9786cb13a780
# ╠═f80a6c49-ad44-420a-bb53-4f196c26e260
# ╠═e523ec46-d420-42c2-979a-b7ec4ad2010d
# ╠═b55852f8-9cbc-4346-a938-b820f04b7ce0
# ╠═20723761-265b-4240-a4a4-132705481f55
# ╠═d5bd5954-4620-4855-87a6-408ce8479552
# ╠═a1a50633-bd5d-4148-8825-c0b0855087a9
# ╠═3979e3b9-dddc-433f-ab69-6fc1265bdf4f
# ╠═a2f24a67-c112-4b10-94e8-10d078c29181
# ╠═2673dfd5-b088-4f35-8669-29e8f9806fe4
# ╠═bf46ed11-6f1c-4f1d-b010-82f6e9988f7c
# ╠═ad1aa93f-ca21-454f-b9ac-10750495113e
# ╠═0aacd779-d9c5-4850-8d20-b766496a19e7
# ╠═d9c433bf-4d96-4ee9-a98c-056454db37ea
# ╠═ecd4080b-38f0-4191-8cd0-d64c6d3674be
# ╠═853787cc-229c-4b0b-afe6-4bc36666dec6
# ╠═9428c207-98e0-4fcd-8477-d2d92a15d3c9
# ╠═b49aaeaf-63ea-4bf1-a61b-beca60c8df98
# ╠═705652f8-be2b-42b6-b76f-672015f78f19
# ╠═a2d9b266-a9ca-4c88-9dbc-26815fafac21
# ╠═fa13e453-124b-43c1-84c0-4d0e7ce1ea94
# ╠═66bf47c3-d538-40de-9dca-aeb560dfea15
# ╠═67cd1b0d-89cc-46ca-ae4b-3900d072ea68
# ╠═321b373a-5fda-4c9c-b1ad-34d1f82d086b
# ╠═aad39fc9-433c-4155-92a7-6462971dbf22
# ╠═015712c4-b2b9-48d4-a076-427a8d0a0baa
# ╠═b0a7208d-26c1-4ea5-9072-f43a70d9e1ca
# ╠═681425bc-3bc4-4f34-802a-af735d745baf
# ╠═07c0a448-8acd-4387-a755-4bdfa06254df
# ╠═f93d9148-1fca-4161-8d27-03d6635da764
# ╠═e9144680-0127-4cbd-8408-6fec1f412fc8
# ╠═d87eaa68-8c56-48c5-b19f-2f95c70f0e6b
# ╠═ab03fda5-fbf7-4ae2-a743-b57e95211dad
# ╠═cea977d9-69a5-4b1a-8e54-c8389d95e2d4
# ╠═f8b6d709-1f14-4691-b7a2-93932a972ce0
# ╠═ab7c483e-f8db-4556-9d78-8a87a7c2cf77
# ╠═fee7bf69-da15-484b-82ce-de0246a49b75
# ╠═6d5bc9d2-3c7a-497b-9372-e4d43eabf3c4
# ╠═b66b58c1-99f7-40c9-824e-74fdf3da4220
# ╠═8065bf42-8d71-4c99-85c6-a2b6fd2e2087
# ╠═73878b74-389a-4306-834b-9d9a126ba95d
# ╠═d73cae53-b429-4c06-a63b-9bc876387256
# ╠═baba5028-d85c-415f-95f6-0393ce686249
# ╠═2174c692-ccc7-490b-b2c4-940876c579dd
# ╠═b3a05a67-3844-4ef8-9bdf-e631477e2b15
# ╠═e9b35fca-b72a-4b2e-8167-92e070017e25
# ╠═01333780-fbab-41cf-8cc5-45677ef971e2
# ╠═92562f1f-3513-4450-af37-a581a8def5f2
# ╠═f00f105c-3c1f-41f0-a3bc-681af74a1f3d
# ╠═bf287007-369e-4995-a3ad-8d5b96aa1b90
# ╠═569f7f0d-c423-4fce-b8f6-33fc2af07528
# ╟─d8f37622-4d4c-4983-9184-b7892559cca3
# ╠═b4c777b1-1245-4ce9-89f5-d0110c453d60
# ╠═938ff2ea-13af-4bce-98d9-8ef5349a7f3d
# ╟─81ac944f-7678-481e-a102-d10cef338ea5
# ╠═94ea1009-a504-44f0-87a6-6449958bd180
# ╠═dc7967a5-9713-4610-ba43-82d24b2481f7
# ╠═13fb12a9-5874-4d1d-9d39-f5903be599b6
# ╠═ef31cd78-5e84-4f6e-bd17-4b575bb17775
# ╠═9d7c75f9-ffe1-47fe-9bab-784f2a7945c1
# ╠═53849542-08be-4c55-9675-f70bded3546b
# ╠═017615a3-98ee-4115-9dad-567121079952
# ╠═a6aa5ce8-e38f-466b-843e-c98cd3a0b556
# ╠═9cb01ba4-2a1b-47f4-83db-859bdd6e8178
# ╠═7bb296f1-6f93-4007-8503-4e5e85d26751
# ╠═bd9d7873-544c-4f4c-acaa-22eac4ae7bc1
# ╠═8c1e081f-c5d3-42e6-8d8b-8f65dddea440
# ╠═3b567b7b-6fc9-48e8-81ef-90dc6ff4b926
# ╠═d2a192c8-7039-493b-8c67-d40ac9c366b6
# ╠═80d9ed58-f4a9-43cd-ab94-251a190f2b3e
# ╠═616f1a1c-86dc-41a4-909b-d78ef0663977
# ╠═8d475de0-2018-41a2-9a4c-d7378e9375fd
# ╠═34f29042-70f2-4a28-8d61-ded8ee6d6d5a
# ╠═025a1b05-7831-485a-9cb7-d580f84f7115
# ╠═6a07f88c-5012-4bd4-b4f9-0cd28e068f12
# ╠═608ba437-7c4e-40cc-a916-eca292a637f2
# ╠═4ff4855f-ccbc-4b55-884c-2f35af2f43d5
# ╠═1fe12f82-acce-4a0f-a767-6f33945db454
# ╠═d8801d90-f74e-4913-8963-75484743303f
# ╠═279cc721-78cf-41be-9120-0e2ceea7300d
# ╠═97256e97-6df1-49cc-a8d3-ebf42041919b
# ╠═d1409c38-01ae-4de2-b4b1-de3ead44c98b
# ╠═fd31ac5a-7248-41fc-b6cc-0563790f07b6
# ╠═a09a9d75-79e7-4f6a-a058-0309901838c4
# ╠═965439d1-698f-47ad-be87-04d7d27e42fc
# ╠═50502497-7576-40a8-917c-840e7c4c52f9
# ╠═9637e720-24f1-4321-bb77-72316dfb7bbd
# ╠═57008d1a-1b6d-42a5-86bb-62a9bd8c2f63
# ╠═5dd92589-d1a9-4eb0-9e58-d173b8ff3a04
# ╠═4a6b9446-8724-4d03-a2d0-48c3d26057b4
# ╠═fa381020-606f-4421-9b56-cb5a0360fada
# ╠═e3413e4a-88ce-47b0-9514-70b1e9245508
# ╠═0b19c137-824f-4d56-a1ee-a08c4133f602
# ╠═b1707796-6aaf-408f-a69d-1414f35d5093
# ╠═d64991e3-bb7c-4702-a119-f02293a3f686
# ╠═a564e0f3-f876-4b37-bb33-83fbdafdd513
# ╠═61dbb1ea-e902-493f-9273-27d7d381ed99
# ╠═bc862bfe-6b3b-444f-8951-653b204ebbf6
# ╠═b4aa6c2e-f682-4ae4-873a-dcabd44fc528
# ╠═612f3f29-a64a-4adc-802f-e5492767f933
# ╠═16de3662-eac1-4dd8-ae77-de224e935add
# ╠═cdf477d2-aa0f-464f-9dc1-9aec9aa583ef
# ╠═087d894a-4ab2-4ef9-8824-da29d5755a44
# ╠═0b2d1160-9c7c-4e9d-9545-4c183fb193d2
# ╠═4a435210-a6de-4c61-8ef8-46fba6931b0f
# ╠═a7226039-1bb3-4d2e-9a43-4a8588cd7b2d
# ╠═0d49b4d3-6e86-414a-9db4-76d059c2b67c
# ╠═7a83f2f4-df2f-4b81-b61c-62095e18d47f
# ╠═75a6ea58-dd9b-4e61-b483-de37856e9229
# ╠═673b7b87-3e74-40f8-978a-efac0f265b2f
# ╠═c1c210fe-03ca-4139-91ac-c9687ebc4ca1
# ╠═286c10f9-b1e5-4088-926c-1b86c11107a3
# ╠═133a7d51-791c-4952-8498-b8c938b2fa0d
# ╠═e2914357-c8c2-4fbf-a0f9-308477242b73
# ╠═7e9fa3f9-a1a9-45c3-a0cc-1591d320f8df
# ╠═b375e4a5-5e99-47da-841b-2475936c437d
# ╠═e424d154-7090-4311-bb03-ecd06bc10f21
# ╠═9c970c92-9cd2-4440-ae40-b41139340c40
# ╠═a9f1938c-93c8-4e6f-91ad-8b38d0b2f6ff
# ╠═12b0d4b3-75f1-41a5-8b22-5c307105f3e5
# ╠═19257a48-69b7-40a5-859e-483d08928cbf
# ╠═3999479a-b980-4e81-a597-f9775896fe56
# ╠═0b083e48-8c7a-4f16-bf28-b5470f23acdb
# ╠═b31b5d94-da48-4e84-8c8f-10a5185b32e9
# ╠═78fb4680-4b8a-42f0-a69d-c68cd2917168
# ╠═fc9bee5c-26aa-4e91-aaa0-c071e8a63737
# ╠═3389bf09-0976-4dde-99e2-7500ea80d695
# ╠═1f570b4c-a5e2-40e3-93b0-6be463215d27
# ╠═e4ba9f6f-e742-4600-864e-50e6794f74d3
# ╠═88e76067-0b95-4e97-8f3b-bd68b6e0409d
# ╠═00b257ad-c634-4f01-b195-4bddb0c074fe
# ╠═690ada82-db76-4e9b-9795-8f60857f6d94
# ╠═0956e274-9339-4314-9baa-7174bdbcde4e
# ╠═d8528f90-0e76-403e-8131-762be3fef133
# ╠═3ead5520-6e2e-4189-806d-a8ad6c8b18e2
# ╠═7a068208-2caf-468a-ae88-df1809cb94d7
# ╠═718a44cf-d893-4574-95d8-7e9c23f47fbc
# ╠═93e82b00-94d9-4691-8e36-653b4f2c3526
# ╠═7e37e56a-f536-49bd-8839-a4af59ff8bf3
# ╠═3744cb91-f3f4-481f-964f-40574066fb86
# ╠═bfb83129-5f50-40fe-87d5-793cfd292bc7
# ╠═c8acd563-df8e-4404-904f-6e0b97a13f55
# ╠═fde49db9-165b-4b55-928b-b0d3f2610a31
# ╠═500c0b3b-fc21-45a7-b03d-07c413a947c9
# ╠═79c95b20-7b71-414f-b5d6-41e4cfe373ca
# ╟─50e2ec32-be70-444a-9c80-38770ce5da7e
# ╠═6be64cc8-a357-477c-ae0e-cb0405de4f4b
# ╟─fe91ac2c-b95a-4d5e-a2fe-f09f43297d4d
# ╠═e5bf7e93-f089-4b77-a2e0-d014bdc309e1
# ╟─0eded161-9ad3-475e-a6d2-d8068f15bad5
# ╠═e673cdda-7032-4600-bbd0-01a291aa3eb5
# ╟─96dad766-398d-414c-ac66-dfaa76947ffc
# ╠═b43916ba-0499-4d82-b81c-49aba21ee145
# ╟─5dfbb58a-1411-4264-9540-8cc26263910f
# ╠═ff872f5f-b10b-4065-8850-d81f9d6aece3
# ╟─5ff2bb27-3db8-4b05-bd90-c73b493a9c76
# ╠═f20b33b0-c0ec-4ee7-8273-6b738eb44dd7
# ╟─c548e453-f5bc-40cc-a3d3-b3c0c72d6f9d
# ╠═b5e2dfb7-c325-4898-b45d-8f071f1958b6
# ╟─c5c7a723-6975-4903-ae84-aee4a1bb5441
# ╠═1f950122-3c05-4d4c-85c2-1ed2094b1448
# ╟─d686b4a3-03fd-4213-9b82-45fd395f96ee
# ╠═35f0fda8-7e9a-4f5c-bbd3-979e8d04ebfa
# ╠═d8a84dea-adbb-48c5-809b-309236f3dfc5
# ╠═3d5ef2d8-1d1c-4f90-9178-f72ce4d67f51
# ╠═647e0043-3389-43f2-bfd6-329d88fce931
# ╠═843d888f-1e8d-4875-b175-1bc0a98f3320
# ╠═99abe607-0a0e-411a-8784-df5318be1db9
# ╠═f6902608-37a4-4c63-98d4-6ad7a704907c
# ╠═c126abdc-a067-45fe-be68-817399f31e65
# ╠═00c178f8-d4b9-4822-937a-f4e82da84f11
# ╠═9d34622a-d7c4-4621-96e8-51e9b801a829
# ╠═58d706c0-44f5-4514-b87b-eaf1f1940448
# ╟─a840300f-0724-44c8-ba08-d7b59491e0e4
# ╠═c1894849-5d29-4857-a715-9b91ffe8163d
# ╟─083a5ba2-7de3-4e51-962f-08289fe987ce
# ╠═773563de-8eb2-4ea4-87fc-86609e35370f
# ╠═cfaf2b26-e39d-47b5-8e99-8d35dbd7b3ac
# ╠═13e729c1-15ff-4fcf-b3df-9247fda54109
# ╠═52e6ee88-f49d-420e-9a99-9c8c54910362
# ╠═604081ba-c975-414a-a8b6-68a6c401172f
# ╠═1928ef8d-d7df-4b2b-9b9a-f958cacecaa3
# ╠═17fffef0-699e-4aa2-8931-41afa2ef4da1
# ╠═0bb31d5b-6d5c-49cc-9ddc-1675dfe4c98f
# ╠═8fe40af8-5374-48ed-831e-69b52a9e3ffc
# ╠═f3b5c2b0-1039-4edf-9efe-6da3665e5b39
# ╠═84c75a93-92ce-4df4-be33-9fd71bffb2ca
# ╠═b0738b82-5ce0-46de-828b-5c25477d2449
# ╠═70f3237c-c051-494a-91ef-ad97bcc07e94
# ╠═c056af02-693b-4837-ab7b-1c4c7b7a9ca5
# ╠═8d02c42c-452e-4ad3-9b27-95ac7ea5e210
# ╠═380a7a82-0924-42b6-8698-7573b389e9ad
# ╠═f62e8d24-10c6-42e9-a68e-f2d68f950ed1
# ╠═e45bd546-c5fc-4914-b237-3df5f473363d
# ╠═5e6ea853-c7a2-44ce-acb5-628ed36f8180
# ╠═8953d6fb-f9ad-4f95-b022-5ac95058dbf7
# ╠═6c95fd00-0023-4f00-b495-050ca098cfb6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
