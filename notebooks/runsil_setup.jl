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

# ╔═╡ 0b4c1052-9180-11ec-2cdd-1b61bde8bdc6
begin
using Plots
using Printf
using Unitful
using DataFrames
using CSV
using Images
end

# ╔═╡ 2d15cff8-622e-4f60-a16f-77a4d16ccc0d
begin
using Interpolations
using QuadGK
using UnitfulEquivalences
end

# ╔═╡ ffc498a9-3522-4309-9ced-48955cab1759
begin
using Markdown
using InteractiveUtils
using PlutoUI
end

# ╔═╡ 91aa2620-9728-4de7-9c71-f01a6c765350
begin
using LsqFit
using Distributions
end

# ╔═╡ 9d66c472-b742-49b3-9495-1c3be998108b
import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
    A, N, mol, mmol, V, L, M

# ╔═╡ c3b89fa7-67f8-4cc3-b3ba-fcb304f44669
import PhysicalConstants.CODATA2018: N_A

# ╔═╡ bd9eed27-fc29-41ff-bf12-ec0cd3881ea2
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

# ╔═╡ ca6f7233-a1d1-403d-9416-835a51d360db
lfi = ingredients("../src/LaserLab.jl")

# ╔═╡ 4c6c9627-d80d-41ac-8141-f935df893b73
md"# Observing fluorescence in a ML of Ru-sil molecules (RUSL)"

# ╔═╡ 34b3d66e-3b39-43dc-922b-9ba6c15c1bac
md"## The RuSL molecules

RuSL molecules are a type of phosphorescente molecules which can be used to calibrate the TOPATU laser setup.

The molecules emit ligh peaked in the red, with a lifetime in the range of 500 ns.
"

# ╔═╡ 7e150a69-0070-4371-b5c0-02b7ad70d813
md"### Fluorescence cross section and quantum yield"

# ╔═╡ 16cdbadf-4deb-4a66-8ab7-84437a4fe3d4
load("../notebooks/img/RuSlAbs.png")

# ╔═╡ 9cc3dfcb-b95d-4086-a228-ed4753f6ca0d
begin
	ϵabs = 8801.5/(M*cm)
	Q    = 0.9
	λexc = 485.0nm
	λEM  = 690.0nm
	println("")
end

# ╔═╡ 4a0e30fe-398e-4d80-86a2-bfc384cbc6e6
md"
- The absorption cross section at 469 nm is ϵ = $ϵabs
- Measured quantum yield is Q= $Q
- The excitation wavelength of the laser is $λexc
"

# ╔═╡ 591f9fcb-7ad1-400e-af0a-29686a4914ec
md"### Emission spectrum on ML

- ML of RuSL, at a nominal packing of 1 molecule per nm
- Peak emission around $λEM
"

# ╔═╡ 768feb2a-855f-4db6-9ae0-6abec23707c1
begin
	spG = lfi.LaserLab.CsvG(';',',')  # spanish: decimals represented with ',' delimited with ';'
	enG = lfi.LaserLab.CsvG(',','.')  # english
	enG2 = lfi.LaserLab.CsvG('\t','.')  # english
	println("")
end

# ╔═╡ 9cc092e9-40d2-4c5a-9b55-4e2adccb3382
begin
	dfrusl = lfi.LaserLab.load_df_from_csv("/Users/jj/JuliaProjects/LaserLab/data/fluorimeter", "RuSL_quartz.csv", spG)
	dfrusl[!, :λ] = convert.(Float64, dfrusl[:, :λ]) 
	first(dfrusl,10)
	println("")
end

# ╔═╡ d75a8426-c97b-422e-8c69-3f3ce94c5370
begin
plot(dfrusl.λ, dfrusl.QUARTZ_Rusilatrane_A, lw=2, label="RuSL on quartz")
xlabel!("λ (nm)")
ylabel!("I (a.u.)")
end

# ╔═╡ a7d330cf-5093-490f-adc8-e482e3084806
md"## Temporal dependence of the phosoprescence"

# ╔═╡ 280a17c4-6d5e-4ff3-a5ab-dbbcb7199bb2
begin
	dfruslt = lfi.LaserLab.load_df_from_csv("/Users/jj/JuliaProjects/LaserLab/data/fluorimeter", "Ru_SL_time.csv", spG)
	dfruslt[!, :λ] = convert.(Float64, dfruslt[:, :λ]) 
	first(dfruslt,10)
	println("")
end

# ╔═╡ a230e061-b6ab-49d4-bdad-730d75e20e9c
md"### The data fits (but not too well) to a single exponential

- At large times the fluorimeter may be measuring a constant pedestal
  "

# ╔═╡ 6d5e295c-1c57-4a92-aee2-3ffcbb3306df
begin
	mexp(t, p) = p[1] * exp.(-t/p[2])
	pa0 = [3500.0, 500.0]
	tdata = dfruslt.λ
	ydata= dfruslt.Ru_Quartz
	fit = curve_fit(mexp, tdata, ydata, pa0)
	cofe = coef(fit)
	stder = stderror(fit)
	println("")
end

# ╔═╡ 80447dc7-1103-4019-bc64-283d4a9368a0
@info "fit coefficients" cofe

# ╔═╡ 08ec5b69-946f-407e-ac8d-bf8814cc0121
@info "coefficient errors (std)" stder

# ╔═╡ 58b227ac-627f-49ba-bc58-75039c65733b
md"## The RUSL experiment

The goal of the RUSL experiment is to use the (modified) TOPATU setup to measure the spectrum (color) and temporal response of the RuSL ML on quartz. This experiment has already been performed in the fluorimeter (see results below), and the main purpose of RUSL is to repeat it with the laser setup to calibrate the system (including temporal response)

The experiment ingredients are:

- A ML of RuSL on quartz. 
- A pulsed laser of 480 nm (EPL 480 from Edimburgh).
- The topatu setup, modified to allow measurement of spectrum (with filters) and time response. 
- The measurement must include spectral response (number of photons observed as a function of the wavelength) and time response (which requires time-stamps for the observed photons)
"

# ╔═╡ 4223f9f8-8ee6-4ea9-a7bc-6379c81c48c0
load("../notebooks/img/Nuevo_SET_UP_RuSl.png") 

# ╔═╡ fabba61a-f28d-4535-b0e3-e94f5098bb0a
md"## The Monolayers

A Monolayer (ML) of fluorophores is characterised by the packing of the fluorophores, the absorption cross section of the individual fluorophores and their emission spectra. These properties may also be correlated or can be modified when moving from solution to ML. For example colective effects of the interference with substrate can modify the naive assumption that the total fluorescence is the product of the number of fluorophores and the fluorescence per fluorophore. However, it is useful to start with the simplest assumption that all cross sections measured in solution hold in solid/gas interface

"

# ╔═╡ 7a1fcff4-6948-4333-b971-8fbf31c7d3ee
md"#### Define the fluorophores at nominal excitation ($λexc)"

# ╔═╡ a917927d-d471-4465-abd1-30d959af0b45
frusl = lfi.LaserLab.Fluorophore(λexc, λEM, ϵabs, Q)

# ╔═╡ aa4283d7-37d4-4aaf-a5b7-b515466e545d
begin
	λepl    = 485.0nm
	Pkepl   = 35.0mW
	P200khz = 0.6μW
	fepl    = 200.0kHz
	wepl    = 140.0ps
	dc      = uconvert(μs, 1.0/fepl)
end

# ╔═╡ 93b78d1a-ee07-4d66-b5a4-c7a109a24f81
md"## Laser
- The experiment requires a VUV pulsed laser. For the experiment we will use a repetition rate of $fepl (pulsed each 1 μs)
- 
- The nominal laser for the experiment is the EPL485
  
| λepl   | Pkepl      | fepl  | wepl | dc| P|
|:-------| ---------- |:-----:|:-----:|:-----:|:-----:|
| $λepl  | $Pkepl| $fepl | $wepl |$dc| $P200khz

"

# ╔═╡ 19b7d7fb-adab-4f29-82af-d059525921ee
md"#### Define laser"

# ╔═╡ 102a2054-2fd5-406e-8bb4-20ecb47c278f
epl485 = lfi.LaserLab.PulsedLaser(λepl, Pkepl, fepl, P200khz, wepl)

# ╔═╡ e97474ef-145f-45c8-96f6-213acf4a3b41
md"### Beam shape
- The beam has an oval shape, with a long axis of 3.5 mm and short axis of 1.5 mm
- The entrance iris diameter of the objective is 2.5 mm. The beam overfills one of the dimensions and not quite the other.
- We will make the approch that the beam fills the entrance iris
"

# ╔═╡ 106f7063-719e-4eff-857d-3566d6e6d4c8
load("../notebooks/img/beamspot.png") 

# ╔═╡ 177762d5-b9f2-449b-abba-256f3d4318ad
md"## Objective

The objective directs the laser light into the ML. For the first round of experiments, it is convenient to focus the laser to the smaller possible spot (diffraction limit). This is done by filling the entrance pupil of the objective with the laser. 

Let's assume a setup in which the back lens of the objective is filled up with a laser beam (assumed to be gaussian). The waist of the beam, assuming $z_r >> f$ (where $f$ is the focal distance of the objective and $z_r$ is the depth of focus of the gaussian beam) is then $w_0 = d/2$, where $d$ is the diameter of the back lens of the objective

The experiment must be conducted in a dry atmosphere, to avoid quenching the phosporescence. Thus the sample must be in a box at vacuum or filled with an inert gas (e.g, argon, N2). The objective may be inside the box (if working distance is small) or outside (if working distance is large)

We will use a reflection objective, the MM40XF-VUV, characterized by a large working distance and large (for an air coupled) NA. 
"

# ╔═╡ 4e4ab4e5-bacc-4c58-8000-10602a1f8465
#load("../notebooks/img/LMM40XVUV.png")  

# ╔═╡ e35f8650-dc2e-46f3-9220-26754e8860e0
#md"### Characteristics of the MM40xVUV

#| Feature   | Value     
#|:-------| ---------- |
#| Entrance pupil diameter  | $dd| 
#| Focal length  | $fl| 
#|NA  | $NA|
#|M (magnification)  | $MM|
#|working distance  | $wd|
#|Tranmission (250-1000 nm)  | $T|
#|Damage threshold  | $dth|
#"

# ╔═╡ 13827ddc-bebd-44d5-929f-f4ac6f43b093
#begin
#	dd    = 5.1mm
#	fl    = 5.0mm
#	NA    = 0.5
#	dth   = 0.3J/cm^2
#	wd    = 7.8mm
#	T     = 0.85
#	MM     = 40.0
#end

# ╔═╡ 2a07dd38-bac4-410c-8135-5c4e6f851df6
#lmm40xf_uvv  = lfi.LaserLab.Objective("LMM40XF-UVVV", fl, dd, MM)

# ╔═╡ ebe4e491-1c61-45cb-a559-64fa3f5fbb9d
load("../notebooks/img/NikonMUE31900.png") 

# ╔═╡ d6d6dc94-5b99-4de6-a5c3-ab2be2b49d32
begin
	dd    = 2.4mm
	fl    = 2.0mm
	NA    = 0.6
	wd    = 10.0mm
	MM    = 100.0
end

# ╔═╡ a3703591-17d4-4049-8c1b-21a4c8329ffd
md"### Characteristics of the NikonMUE31900

| Feature   | Value     
|:-------| ---------- |
| Entrance pupil diameter  | $dd| 
| Focal length  | $fl| 
|NA  | $NA|
|M (magnification)  | $MM|
|working distance  | $wd|
"

# ╔═╡ b0112287-c957-4ce3-ba24-9c47c8128b2d
tobj  = lfi.LaserLab.Objective("Nikon", fl, dd, MM)

# ╔═╡ 417fecfd-b316-4e63-8e68-303d4c568434
md"## The laser beam as a Gaussian Laser"

# ╔═╡ bbd09112-d6cd-4f0f-9a87-7b85be983e09
md"### Focusing the beam"

# ╔═╡ 8c405cc2-e8b5-4125-b21d-f7a22f94fb59
md"The beam is now focused in a narrow spot by the objective. "

# ╔═╡ ebe77978-7d7f-407d-9a73-4be3568265ef
gepl485 = lfi.LaserLab.propagate_paralell_beam(epl485, tobj)

# ╔═╡ e7e7b1ca-b27f-406c-9df0-5ffea702719f
md"### Spot size and Depth of focus

The spot size is $(round(lfi.LaserLab.spot_size(gepl485)/nm,digits=1)) nm, while the depth of focus is $(round(lfi.LaserLab.depth_of_focus(gepl485)/nm, digits=1)) nm. "

# ╔═╡ f48eb494-f9f5-4e4f-9fc0-ffa20312155a
i0mWcm2f = round(uconvert(W/cm^2, gepl485.I0)/(W*cm^-2), sigdigits=2);

# ╔═╡ fc5b2f60-2c04-4835-a222-9971189ac54e
ng0 = round(gepl485.γ0/(Hz*cm^-2), sigdigits=2);

# ╔═╡ 4eeabd93-cd82-4051-b2d0-bee3db1723e4
md"### Power density 

- The power density in the spot is now much larger: $(i0mWcm2f ) W/cm2
- Or in term of photon density: $(ng0) Hz/cm2
"

# ╔═╡ e40a68c6-5885-49d1-9259-fdd3be09fee4
md"In a gaussian laser, the intensity as a function of the radial direction ($\rho$) and the direction of propagation (z) is:

$I(\rho, z) = I_0 ( W_0 / W(z))^2 \exp{-2 \rho^2/W^2(z)}$
"

# ╔═╡ f79fb47b-b4f7-442a-927f-0c186e673b80
md"Since both the spot size and the depth of focus are small, we can approximate the intensity that will illuminate the molecules of the mono layer with I(0,0)"

# ╔═╡ 8da11b18-5934-4207-941b-795969173f21
fI = lfi.LaserLab.I(gepl485) ;

# ╔═╡ ce47be17-7773-4797-beb7-202d3c02555a
begin()
zl=-5.0:0.01:5.0
p1 = plot(zl, fI.(0.0*μm, zl*μm)/(mW * cm^-2), label="I(0,z)")
xlabel!("z (μm)")
ylabel!("I(0,z)")

rl=-1.0:0.01:1.0
p2 = plot(rl, fI.(rl*μm, 0*μm)/(mW * cm^-2), label="I(ρ,0)")
xlabel!("ρ (μm)")
ylabel!("I(0,z)")

plot(p1,p2, layout = (1, 2), legend=false, fmt = :png)
end

# ╔═╡ 498d145b-ce5f-46a5-b77e-31f412e04eb9
md"### Packing of the monolayer

- One can define the packing of the ML in terms of the 'pitch' or distance separating two molecules. This is a crucial parameter of the experiment.

"

# ╔═╡ 0d4a0651-4ce4-497e-8764-2bcbbef83cf1
md"##### molecular pitch (in nm)"

# ╔═╡ d437cc55-8c19-4d0b-86d8-760da3956895
@bind mp NumberField(1.0:10.0^3; default=1.0)

# ╔═╡ 7dedd16b-4279-4f50-8120-d9ef394f3e13
pitch = mp*nm;

# ╔═╡ 649fa1b3-0dd2-487d-9bdf-7db97a0ec178
ml = lfi.LaserLab.Monolayer(pitch);

# ╔═╡ 92a5b436-1d8e-4436-8dda-8b1d3518bdea
md"This corresponds to $(uconvert(cm^-2, ml.σ)) molecules"

# ╔═╡ 8d4516fd-8838-4e5f-a61d-9dc1f65b31ad
aspot = π * gepl485.w0^2;

# ╔═╡ c5b0405a-a991-45a5-aa04-d09020b0c7f0
md"### Spot area
- The area of the spot iluminated by the beam is $(round(uconvert(μm^2,aspot)/μm^2, digits=1)) μm2"

# ╔═╡ 03e15d19-ea63-413e-8ec5-4d50e28558ac
nmol = uconvert(μm^2, aspot) * ml.σ;

# ╔═╡ 4a1cf5fc-b44f-4b19-ac8e-d610da0a17cb
md"### Number of molecules in the spot
- The number of molecules in the spot illuminated by the beam is: $(round(nmol))"

# ╔═╡ 39a7f7b4-058b-4283-b35c-4409ea9e478a
md"### Fluorescence per molecule
- The fluorescence of each molecule is the product of the beam density, the fluorescence cross section and the quantum yield:

$f = I_0 \cdot \sigma \cdot Q$
"

# ╔═╡ abe34c91-9824-45e4-8865-7b3f73ff8758
fmrs = lfi.LaserLab.fluorescence(frusl, gepl485.γ0);

# ╔═╡ 52aeae45-0dfc-48f2-a362-20d1add0ff7f
md"### Fluorescence per molecule for free and chelated species
- The fluorescence per molecule for RuSL is: $(round(fmrs/Hz)) Hz, 
"

# ╔═╡ 9d97d16a-8f90-4ff7-ac0b-ee610c20ee32
sfrs = fmrs * nmol;

# ╔═╡ 960ba1ef-495d-4bb2-8aff-73cb68ae440e
md"### Total fluorescence in the spot

The total fluorescence in the spot is the product of the number of molecules and the fluorescence per molecule:

- Fluorescence in spot for RuSL = $(round(sfrs/Hz, sigdigits=1)) Hz

"

# ╔═╡ 8815a3e2-e559-4857-b20e-14aa5f91d341
md"## Filters and dichroics"

# ╔═╡ 89aa75b3-2263-450e-9f2e-cd7c875a7919
md"### Dichroic DMLP 567
- Transmits > 98% of the light above 600 nm.
"

# ╔═╡ b5045344-fb52-4b2e-88b9-1c1d4d1d50f6
begin
	dmlp567 = lfi.LaserLab.load_df_from_csv("/Users/jj/JuliaProjects/LaserLab/data/Filters", "dmlp567.csv", spG)
	dmlp567[!, :λ] = convert.(Float64, dmlp567[:, :λ]) 
	first(dmlp567,10)
	println("")
end

# ╔═╡ f6fefa8e-8490-4c0d-b707-c9a95011cbb1
begin
	dp1 = plot(dmlp567.λ, dmlp567.T, lw=2, label = "T", fmt = :png)
	plot(dp1, dmlp567.λ, dmlp567.R, lw=2, label = "R", fmt = :png)
	xlabel!("λ (nm)")
	ylabel!("T (R) (%)")
end

# ╔═╡ 3555a982-9879-460f-babf-6f31c11c6f7f
md"### High band pass filter FGL550"

# ╔═╡ 954d503d-19f2-4575-84d7-2d1722a28a7f
begin
	fgl550 = lfi.LaserLab.load_df_from_csv("/Users/jj/JuliaProjects/LaserLab/data/Filters", "fgl550.csv", spG)
	fgl550[!, :λ] = convert.(Float64, fgl550[:, :λ]) 
	first(fgl550,10)
	println("")
end

# ╔═╡ 41eb41f7-bc6e-4dce-8d39-68142fd329db
begin
	plot(fgl550.λ, fgl550.T, lw=2, label = "T", fmt = :png)
	xlabel!("λ (nm)")
	ylabel!("T  (%)")
end

# ╔═╡ bc6ba394-d1b8-4c1b-8a01-c23bcd29178c
md"### Color filters FF01"

# ╔═╡ 7140a27e-4cde-4c7d-a794-9a624e540677
begin
	ff550 = lfi.LaserLab.load_df_from_csv("/Users/jj/JuliaProjects/LaserLab/data/Filters", "FF01-550_49_Spectrum.csv", enG2)
	ff600 = lfi.LaserLab.load_df_from_csv("/Users/jj/JuliaProjects/LaserLab/data/Filters", "FF01-600_52-25.csv", enG2)
	ff650 = lfi.LaserLab.load_df_from_csv("/Users/jj/JuliaProjects/LaserLab/data/Filters", "FF01-650_54_Spectrum.csv", enG2)
	ff692 = lfi.LaserLab.load_df_from_csv("/Users/jj/JuliaProjects/LaserLab/data/Filters", "FF01-692_40_Spectrum.csv", enG2)
	ff732 = lfi.LaserLab.load_df_from_csv("/Users/jj/JuliaProjects/LaserLab/data/Filters", "FF01-732_68_Spectrum.csv", enG2)
	ff810 = lfi.LaserLab.load_df_from_csv("/Users/jj/JuliaProjects/LaserLab/data/Filters", "FF01-810_Spectrum.csv", spG)
	ff810[!, :T] = ff810[!, :T] * 0.01 
	first(ff550,10)
	println("")
end

# ╔═╡ b0bf01be-6a53-4c42-9d9a-752bf4652986
begin
	ff550s = filter(df -> df.T >= 0.01, ff550)
	ff600s = filter(df -> df.T >= 0.01, ff600)
	ff650s = filter(df -> df.T >= 0.01, ff650)
	ff692s = filter(df -> df.T >= 0.01, ff692)
	ff732s = filter(df -> df.T >= 0.01, ff732)
	ff810s = filter(df -> df.T >= 0.01, ff810)
	println("")
end

# ╔═╡ 0ee9bd36-7192-4922-abaf-7d7d08605915
begin
	p550 = plot(ff550s.λ, ff550s.T, lw=2, label = "FF01-550", fmt = :png)
	p600 = plot(p550, ff600s.λ, ff600s.T, lw=2, label = "FF01-600", fmt = :png)
	p650 = plot(p600, ff650s.λ, ff650s.T, lw=2, label = "FF01-650", fmt = :png)
	p692 = plot(p650, ff692s.λ, ff692s.T, lw=2, label = "FF01-692", fmt = :png)
	p732 = plot(p692, ff732s.λ, ff732s.T, lw=2, label = "FF01-732", fmt = :png)
	p810 = plot(p732, ff810s.λ, ff810s.T, lw=2, label = "FF01-810", fmt = :png)
	xlabel!("λ (nm)")
	ylabel!("T  ")
end

# ╔═╡ df1ff4c2-8c84-455a-9cc4-7fdc34e2cb83
begin
	λm = [550.0, 600.0, 650.0, 692.0, 732.0] * nm
	Im = [0.23, 0.43, 3.44, 1.7, 1.22] * MHz
	wm = [size(ff550s)[1], size(ff600s)[1], size(ff650s)[1], size(ff692s)[1], size(ff732s)[1]] * nm 
	Ilm = uconvert.(MHz*μm^-1, Im ./λm) 
	Ilmx = round.(Ilm ./ (MHz*μm^-1),  sigdigits=2)
	It = sum(Im)
end

# ╔═╡ d5ef9a03-f716-454b-a5de-0c67ee935679
md"## Summary of measurements with TOPATU setup and RuSL ML
- P = $(epl485.P)
- f = $(epl485.f) 
  
| λ filter (nm)   | FF550      | FF600  | FF650 | FF692| FF732|
|:-------| ---------- |:-----:|:-----:|:-----:|:-----:|
| λ0 (nm)  | $(λm[1])| $(λm[2]) | $(λm[3]) |$(λm[4])| $(λm[5])
| w (nm)  | $(wm[1])| $(wm[2]) | $(wm[3]) |$(wm[4])| $(wm[5])
| I (MHz)  | $(Im[1])| $(Im[2]) | $(Im[3]) |$(Im[4])| $(Im[5])
| Ilm (MHz/nm)  | $(Ilmx[1])| $(Ilmx[2]) | $(Ilmx[3]) |$(Ilmx[4])| $(Ilmx[5])

Total observed light = $It

"

# ╔═╡ 6396d4af-e668-4a73-b8be-b73b8b41267c
begin
ps2 = scatter(λm/nm, Ilmx, markersize=3,
		color = :black,
	    legend=false,
		fmt = :png)
plot(ps2, λm/nm, Ilmx, lw=2)
xlabel!("λ (nm)")
ylabel!("Il (MHz/μm)")
end

# ╔═╡ 6445e2c6-4dcc-44dd-bdda-564e4c9b3911
effCCD = lfi.LaserLab.ccd()

# ╔═╡ 7e0b0616-5dd0-44d3-beab-4fa32521d3ff
begin
	wl = 350.0:10.0:1000.0
	eccd = effCCD.(collect(wl))
	plot(wl, eccd, lw=2)
	xlabel!("λ (nm)")
	ylabel!("ϵ")
end

# ╔═╡ b49c15a7-d9de-4942-b09c-7f2ed9b4550e
function ccd_eff(lmn,wmn)
	effx = zeros(1,5)
	for i in 1:5
		effx[i] = (effCCD(lmn[i] - wmn[i]/2.0) + effCCD(lmn[i] + wmn[i]/2.0))/2.0
	end
	#@info "effx" effx
	sum(effx)/length(effx)
end

# ╔═╡ 356b7b35-6794-40d0-8c88-b8e066f086a6
begin
	lmn = λm/nm
	wmn = wm/nm
	ϵobj = 0.95
	ϵd = 0.95
	ϵf = 0.85
	ϵPMT = 0.1
	ϵNA      = lfi.LaserLab.transmission(tobj)
	efccd = ccd_eff(lmn,wmn)
	ϵT =  ϵobj^2 * ϵd^2 * ϵf * ϵPMT * ϵNA
	ϵTc =  ϵobj^2 * ϵd^2 * ϵf * efccd * ϵNA
	qf = 0.1
	println("")
end

# ╔═╡ f38318dd-f7bb-4bf1-9bc9-d1f7ed8a8397
md"## Detected Light

The detected light is the product of the emitted fluorescence and the detection efficiency, which in turns includes:

- Transmission efficiency of the objective : ϵ_obj_vuv =$(ϵobj^2)
- Tranmission due to the NA ϵ_NA = $(round(ϵNA, digits=2))
- Transmission due to the dichroic ϵd =$(ϵd^2)
- Tranmission due to the filters ϵf = $ϵf
- Transmission due to PMT ϵPMT = $ϵPMT
- Transmission due to CCD ϵCCD = $(round(efccd, digits=2))
- Quenching factor of phosphorescence due to oxygen = $qf
- The total tranmission (with PMT) is $(round(ϵT, digits=3))
- The total tranmission (with CCD) is $(round(ϵTc, digits=3))
"

# ╔═╡ f615fc39-bfd9-45ce-84fe-a28921bde525


# ╔═╡ 43e8a5f1-512c-46b7-91f6-89d3c7e81368
begin
	osf  = sfrs * ϵT * qf
	osfc = sfrs * ϵTc * qf
end

# ╔═╡ d49e7e9b-6487-407b-aef7-2884461879a0
md" ### Expected observed light in the PMT

Thus the total expected light in the PMT for a ML of ~ $(uconvert(cm^-2, ml.σ)) molecules is:

- Expected observed light:
  -  with PMT= $(round(osf/Hz, sigdigits=1)) Hz
  -  with CCD= $(round(osfc/Hz, sigdigits=1)) Hz

"

# ╔═╡ ae5f17d2-3c31-4a9c-85e4-f2f22625d86b
begin
	P0 = 1.5μW
	f0  = 500kHz
	println("")
end

# ╔═╡ c3952378-79bb-4605-8ffb-828c1e3e3321
#load("../notebooks/img/RuSlCCD250222.png")   

# ╔═╡ 0f7ff8a4-43e4-4829-a4b4-78594664cee2
md"## TCSP experiments"

# ╔═╡ 131f4e45-85f5-43bc-8c32-d59f8803bfd6
md"### Typical setup"

# ╔═╡ a4066207-4c50-4360-b43f-218d5355ff3e
load("../notebooks/img/tcsp_setup.png")   

# ╔═╡ 218fc76e-1060-40ff-a52f-d884441380e2
md"### Photon counting technique"

# ╔═╡ 58b61eac-e1de-48f3-a37b-ecc1f8d5f8ab
load("../notebooks/img/tcspc.png")  

# ╔═╡ 320f9c75-5047-4a4c-890b-ffdc70767634
md" The phothon counting technique requires that one photon is recorded on average per pulse (in this case in the period 1-5 mus). This is to avoid inefficiencies which can bias the measurement as ilustrated below."

# ╔═╡ b97cd0a9-ee46-4520-a408-f60150bbea74
load("../notebooks/img/tcspc_deadtime.png")  

# ╔═╡ b96961fd-6a9a-40bb-b582-2b3586f56edb
md"However, most TCSPC techniques are designed for short interval times. The dead time is typically a few nanoseconds, thus no loss for dead time is expected here. Nevertheless the system is designed to record just one photon per pulse. This implies that the average number of photons must be reduced to 1 per 1 (5) μs. This can be done by:
- attenuating the laser light
- reducing the laser power
- spacing the molecules in the ML
"

# ╔═╡ eba3554f-ae6f-4a3d-8131-43ffa2743977
md"# Appendix"

# ╔═╡ f57ca807-9b3a-4f77-9607-74ee3a411990
md"## Fitting"

# ╔═╡ 7a6a3b8a-8b1c-41d5-ab25-7d294f1bee3d
md"### *func1dfit* is a light wrapper to curve_fit"

# ╔═╡ 92698c05-edf3-4b9e-a10a-1da9ed0dc82a
"""
    func1dfit(ffit::Function, x::Vector{<:Real},
              y::Vector{<:Real}, p0::Vector{<:Real},
              lb::Vector{Float64}, ub::Vector{Float64})

Fit a function to the data x, y with start prediction p0
and return coefficients and errors.
"""
func1dfit(ffit::Function, x::Vector{<:Real},
          y::Vector{<:Real}, p0::Vector{<:Real},
          lb::Vector{<:Real}, ub::Vector{<:Real}) = curve_fit(ffit, x, y, p0, 
			                                                 lower=lb, upper=ub)
    

# ╔═╡ 0dd448d1-9952-438c-9284-e0c320c955aa
md"### Example: fit to a polynomial"

# ╔═╡ f8b61a8e-2679-48cd-842c-17d6e6ee760e

pol3(x, a, b, c, d) = a + b*x + c*x^2 + d*x^3


# ╔═╡ bfc6e346-3525-47d7-a436-ffa02dab11b9
begin
	err_sigma = 0.04
	x=collect(LinRange(0., 10., 100))
	p0 = [10.0, 1.0, 0.7, 0.5]
	y = pol3.(x, p0...)
    y += rand(Normal(0, err_sigma), length(y))
    lb = fill( 0.0, length(p0))
    ub = fill(20.0, length(p0))
    pol3_fit = @. pol(x, p) = p[1] + p[2] * x + p[3] * x^2 + p[4] * x^3
    fq = func1dfit(pol3_fit, x, y, p0, lb, ub)
	cfq = coef(fq)
	sfq = stderror(fq)
end

# ╔═╡ 2ff6f551-2662-43a5-9776-70951ff40364
@info "fit coefficients" cfq

# ╔═╡ ec3cd40c-6a52-48fd-a3ce-9e091250e981
yf = pol3.(x, cfq...);

# ╔═╡ 525b0316-8374-4882-aa5b-84f2bf33c5c3
@info "coefficient errors (std)" sfq

# ╔═╡ 840fa50e-bf0f-42f9-8e7f-ec3cc4bdb9af
sfq

# ╔═╡ 5738e27c-5307-4623-bd43-67f66b7b97d2
@info "margin_of_error (90%)" margin_error(fq, 0.1)

# ╔═╡ 14360d77-340e-488d-bed0-00f408ef1dd4
all(isapprox.(cfq, p0; atol=err_sigma))

# ╔═╡ 685fd840-64b1-427d-83f3-217664ea9798
begin
pp1 = scatter(x, y,
	          label="p3",
			  markersize=2,
			  color = :black,
	          legend=false,
		   	  fmt = :png)
pp2 = plot(pp1, x, yf, lw=2)
	
xlabel!("x")
ylabel!("p3(x)")
end

# ╔═╡ 8ecc8206-275b-4d55-9b01-e82e7f2df5dc
md"### Fit an exponential"

# ╔═╡ aac7a640-97fe-46c0-89b3-7e76978baf0d
exp(1.0)

# ╔═╡ 89ef5f08-7a1e-42a5-8fcf-b7efa63a8a68
expo(t, N, λ) = N*exp(-t/λ)

# ╔═╡ 426f8190-7865-4bd1-bc23-0adb5fc1892c
begin
	err_sigma2 = 0.05
	t=collect(LinRange(0.0, 5000.0, 1000))
	p0t = [3500.0, 500.0]
	yt = expo.(t, p0t...)
	yts = yt + rand(Normal(0, err_sigma2), length(yt)) .* yt
    lbt = [0.0, 0.0]
    ubt = [50000.0, 50000.0]
    expo_fit = @. expo(t, p) = p[1]*exp(-t/p[2])
    fqe = func1dfit(expo_fit, t, yts, p0t, lbt, ubt)
	cfqe = coef(fqe)
	sfqe = stderror(fqe)
end

# ╔═╡ 2fdad658-0b71-4d42-8591-9284fee3aeb6

begin
	tft = expo.(tdata, coef(fit)...);
	println("")
end

# ╔═╡ 5a0cc7d8-a3c3-45b9-8b91-9528e4938fb0
begin
ps1 = scatter(tdata, ydata, markersize=1,
		color = :black,
	    legend=false,
		fmt = :png)
pp = plot(ps1, tdata, tft, lw=2,fmt = :png)
xlabel!("t (ns)")
ylabel!("I (a.u.)")
end

# ╔═╡ f1e6c6b6-49fa-44bb-bb9a-07986f6a1d13
expo(0,3000.,500.)

# ╔═╡ c54205d5-ea1a-4416-b4b0-3bdb168dae61
fqe.converged

# ╔═╡ b83d4ff4-cc73-4cf7-abad-78da291eb404
@info "fit coefficients" cfqe

# ╔═╡ 13ad7221-0d2d-4a77-9258-9edace85fde0
@info "coefficient errors (std)" sfqe

# ╔═╡ f2b2cc4d-54ed-4f2f-80bc-bc3bf82bb2e8
begin
pp3 = scatter(t, yt,
	          label="expo",
			  markersize=2,
			  color = :black,
	          legend=false,
		   	  fmt = :png)
#pp2 = plot(pp1, x, yf, lw=2)
	
xlabel!("t")
ylabel!("expo(t)")
end

# ╔═╡ Cell order:
# ╠═0b4c1052-9180-11ec-2cdd-1b61bde8bdc6
# ╠═2d15cff8-622e-4f60-a16f-77a4d16ccc0d
# ╠═ffc498a9-3522-4309-9ced-48955cab1759
# ╠═91aa2620-9728-4de7-9c71-f01a6c765350
# ╠═9d66c472-b742-49b3-9495-1c3be998108b
# ╠═c3b89fa7-67f8-4cc3-b3ba-fcb304f44669
# ╠═bd9eed27-fc29-41ff-bf12-ec0cd3881ea2
# ╠═ca6f7233-a1d1-403d-9416-835a51d360db
# ╟─4c6c9627-d80d-41ac-8141-f935df893b73
# ╟─34b3d66e-3b39-43dc-922b-9ba6c15c1bac
# ╟─7e150a69-0070-4371-b5c0-02b7ad70d813
# ╠═16cdbadf-4deb-4a66-8ab7-84437a4fe3d4
# ╟─4a0e30fe-398e-4d80-86a2-bfc384cbc6e6
# ╠═9cc3dfcb-b95d-4086-a228-ed4753f6ca0d
# ╟─591f9fcb-7ad1-400e-af0a-29686a4914ec
# ╠═768feb2a-855f-4db6-9ae0-6abec23707c1
# ╠═9cc092e9-40d2-4c5a-9b55-4e2adccb3382
# ╠═d75a8426-c97b-422e-8c69-3f3ce94c5370
# ╟─a7d330cf-5093-490f-adc8-e482e3084806
# ╟─280a17c4-6d5e-4ff3-a5ab-dbbcb7199bb2
# ╟─a230e061-b6ab-49d4-bdad-730d75e20e9c
# ╠═6d5e295c-1c57-4a92-aee2-3ffcbb3306df
# ╠═80447dc7-1103-4019-bc64-283d4a9368a0
# ╠═08ec5b69-946f-407e-ac8d-bf8814cc0121
# ╠═2fdad658-0b71-4d42-8591-9284fee3aeb6
# ╠═5a0cc7d8-a3c3-45b9-8b91-9528e4938fb0
# ╟─58b227ac-627f-49ba-bc58-75039c65733b
# ╠═4223f9f8-8ee6-4ea9-a7bc-6379c81c48c0
# ╟─fabba61a-f28d-4535-b0e3-e94f5098bb0a
# ╟─7a1fcff4-6948-4333-b971-8fbf31c7d3ee
# ╠═a917927d-d471-4465-abd1-30d959af0b45
# ╠═93b78d1a-ee07-4d66-b5a4-c7a109a24f81
# ╟─aa4283d7-37d4-4aaf-a5b7-b515466e545d
# ╟─19b7d7fb-adab-4f29-82af-d059525921ee
# ╠═102a2054-2fd5-406e-8bb4-20ecb47c278f
# ╟─e97474ef-145f-45c8-96f6-213acf4a3b41
# ╟─106f7063-719e-4eff-857d-3566d6e6d4c8
# ╟─177762d5-b9f2-449b-abba-256f3d4318ad
# ╟─4e4ab4e5-bacc-4c58-8000-10602a1f8465
# ╟─e35f8650-dc2e-46f3-9220-26754e8860e0
# ╟─13827ddc-bebd-44d5-929f-f4ac6f43b093
# ╟─2a07dd38-bac4-410c-8135-5c4e6f851df6
# ╟─a3703591-17d4-4049-8c1b-21a4c8329ffd
# ╟─ebe4e491-1c61-45cb-a559-64fa3f5fbb9d
# ╟─d6d6dc94-5b99-4de6-a5c3-ab2be2b49d32
# ╠═b0112287-c957-4ce3-ba24-9c47c8128b2d
# ╟─417fecfd-b316-4e63-8e68-303d4c568434
# ╟─bbd09112-d6cd-4f0f-9a87-7b85be983e09
# ╟─8c405cc2-e8b5-4125-b21d-f7a22f94fb59
# ╟─ebe77978-7d7f-407d-9a73-4be3568265ef
# ╟─e7e7b1ca-b27f-406c-9df0-5ffea702719f
# ╟─f48eb494-f9f5-4e4f-9fc0-ffa20312155a
# ╟─4eeabd93-cd82-4051-b2d0-bee3db1723e4
# ╟─fc5b2f60-2c04-4835-a222-9971189ac54e
# ╟─e40a68c6-5885-49d1-9259-fdd3be09fee4
# ╟─f79fb47b-b4f7-442a-927f-0c186e673b80
# ╠═8da11b18-5934-4207-941b-795969173f21
# ╟─ce47be17-7773-4797-beb7-202d3c02555a
# ╟─498d145b-ce5f-46a5-b77e-31f412e04eb9
# ╟─0d4a0651-4ce4-497e-8764-2bcbbef83cf1
# ╟─d437cc55-8c19-4d0b-86d8-760da3956895
# ╟─7dedd16b-4279-4f50-8120-d9ef394f3e13
# ╟─649fa1b3-0dd2-487d-9bdf-7db97a0ec178
# ╟─92a5b436-1d8e-4436-8dda-8b1d3518bdea
# ╟─c5b0405a-a991-45a5-aa04-d09020b0c7f0
# ╟─8d4516fd-8838-4e5f-a61d-9dc1f65b31ad
# ╟─4a1cf5fc-b44f-4b19-ac8e-d610da0a17cb
# ╟─03e15d19-ea63-413e-8ec5-4d50e28558ac
# ╟─39a7f7b4-058b-4283-b35c-4409ea9e478a
# ╟─52aeae45-0dfc-48f2-a362-20d1add0ff7f
# ╟─abe34c91-9824-45e4-8865-7b3f73ff8758
# ╟─960ba1ef-495d-4bb2-8aff-73cb68ae440e
# ╟─9d97d16a-8f90-4ff7-ac0b-ee610c20ee32
# ╟─8815a3e2-e559-4857-b20e-14aa5f91d341
# ╟─89aa75b3-2263-450e-9f2e-cd7c875a7919
# ╟─b5045344-fb52-4b2e-88b9-1c1d4d1d50f6
# ╟─f6fefa8e-8490-4c0d-b707-c9a95011cbb1
# ╟─3555a982-9879-460f-babf-6f31c11c6f7f
# ╟─954d503d-19f2-4575-84d7-2d1722a28a7f
# ╟─41eb41f7-bc6e-4dce-8d39-68142fd329db
# ╟─bc6ba394-d1b8-4c1b-8a01-c23bcd29178c
# ╟─7140a27e-4cde-4c7d-a794-9a624e540677
# ╠═b0bf01be-6a53-4c42-9d9a-752bf4652986
# ╟─0ee9bd36-7192-4922-abaf-7d7d08605915
# ╠═d5ef9a03-f716-454b-a5de-0c67ee935679
# ╠═df1ff4c2-8c84-455a-9cc4-7fdc34e2cb83
# ╠═6396d4af-e668-4a73-b8be-b73b8b41267c
# ╠═f38318dd-f7bb-4bf1-9bc9-d1f7ed8a8397
# ╠═356b7b35-6794-40d0-8c88-b8e066f086a6
# ╠═6445e2c6-4dcc-44dd-bdda-564e4c9b3911
# ╠═7e0b0616-5dd0-44d3-beab-4fa32521d3ff
# ╠═b49c15a7-d9de-4942-b09c-7f2ed9b4550e
# ╠═f615fc39-bfd9-45ce-84fe-a28921bde525
# ╠═d49e7e9b-6487-407b-aef7-2884461879a0
# ╠═43e8a5f1-512c-46b7-91f6-89d3c7e81368
# ╟─ae5f17d2-3c31-4a9c-85e4-f2f22625d86b
# ╠═c3952378-79bb-4605-8ffb-828c1e3e3321
# ╟─0f7ff8a4-43e4-4829-a4b4-78594664cee2
# ╟─131f4e45-85f5-43bc-8c32-d59f8803bfd6
# ╠═a4066207-4c50-4360-b43f-218d5355ff3e
# ╟─218fc76e-1060-40ff-a52f-d884441380e2
# ╟─58b61eac-e1de-48f3-a37b-ecc1f8d5f8ab
# ╠═320f9c75-5047-4a4c-890b-ffdc70767634
# ╠═b97cd0a9-ee46-4520-a408-f60150bbea74
# ╠═b96961fd-6a9a-40bb-b582-2b3586f56edb
# ╠═eba3554f-ae6f-4a3d-8131-43ffa2743977
# ╠═f57ca807-9b3a-4f77-9607-74ee3a411990
# ╠═7a6a3b8a-8b1c-41d5-ab25-7d294f1bee3d
# ╠═92698c05-edf3-4b9e-a10a-1da9ed0dc82a
# ╠═0dd448d1-9952-438c-9284-e0c320c955aa
# ╠═f8b61a8e-2679-48cd-842c-17d6e6ee760e
# ╠═bfc6e346-3525-47d7-a436-ffa02dab11b9
# ╠═2ff6f551-2662-43a5-9776-70951ff40364
# ╠═ec3cd40c-6a52-48fd-a3ce-9e091250e981
# ╠═525b0316-8374-4882-aa5b-84f2bf33c5c3
# ╠═840fa50e-bf0f-42f9-8e7f-ec3cc4bdb9af
# ╠═5738e27c-5307-4623-bd43-67f66b7b97d2
# ╠═14360d77-340e-488d-bed0-00f408ef1dd4
# ╠═685fd840-64b1-427d-83f3-217664ea9798
# ╠═8ecc8206-275b-4d55-9b01-e82e7f2df5dc
# ╠═aac7a640-97fe-46c0-89b3-7e76978baf0d
# ╠═89ef5f08-7a1e-42a5-8fcf-b7efa63a8a68
# ╠═f1e6c6b6-49fa-44bb-bb9a-07986f6a1d13
# ╠═426f8190-7865-4bd1-bc23-0adb5fc1892c
# ╠═c54205d5-ea1a-4416-b4b0-3bdb168dae61
# ╠═b83d4ff4-cc73-4cf7-abad-78da291eb404
# ╠═13ad7221-0d2d-4a77-9258-9edace85fde0
# ╠═f2b2cc4d-54ed-4f2f-80bc-bc3bf82bb2e8
