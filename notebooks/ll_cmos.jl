### A Pluto.jl notebook ###
# v0.19.11

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

# ╔═╡ 981730a6-61fc-484b-ba3c-66920ee7cf83
using Pkg; Pkg.activate(ENV["JLaserLab"])

# ╔═╡ 8e7ec382-c738-11ec-3aae-b50d60f15c4f
begin
	using PlutoUI
	using CSV
	using DataFrames
	using Images
	using ImageBinarization
	using Colors
	using Plots
	using Printf
	using Interpolations
	using QuadGK
	using Markdown
	using InteractiveUtils
	using LsqFit
	using Statistics
	using StatsBase
	using Unitful 
	using UnitfulEquivalences 
	using PhysicalConstants
	using Peaks
	using FFTW
	using DSP
	using Clustering
	import Glob
end

# ╔═╡ f7bb5111-2fc9-49df-972a-0737182da98c
ENV["JLaserLab"]

# ╔═╡ 06b8ed45-43bc-464f-89c0-dc0406312b81
import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
    A, N, mol, mmol, V, L, M

# ╔═╡ 8833b198-03e4-4679-8949-0c76546cb847
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

# ╔═╡ 6163ba69-1237-4b49-988e-9a73cfef67f6
lfi = ingredients("../src/LaserLab.jl")

# ╔═╡ 5edc41bc-b912-44bf-9be5-a013f27a75ab
PlutoUI.TableOfContents(title="Laser Lab CMOS analysis", indent=true)

# ╔═╡ 90c97a39-b35f-44ba-9646-f0bb9eead338
md"""
# G2SL Characteristics
"""

# ╔═╡ d30e1ceb-2e90-438e-8554-228aa5dc2a59
begin
g2df = lfi.LaserLab.load_df_from_csv("/Users/jjgomezcadenas/LaserLab/Proyectos/FLUORI/G2/G2SL/BOLD_104_SIL_GB_onquartz", 
"Emisi_espec_G2SIL_quartz.csv", lfi.LaserLab.enG)
	g2df = select!(g2df, [:L, :I])
end

# ╔═╡ 8bda9d8a-c928-431d-a441-ca64bacf7651
begin
	plot(g2df.L,g2df.I, lw=2, label="G2SL")
	xlabel!("λ (nm)")
	ylabel!("arbitrary units")
end

# ╔═╡ 2be5dcc0-e7c4-412b-990a-d7edb9967186
md"""
- In the spectrum shown above, signal below 500 nm is most likely an artifact.
"""

# ╔═╡ 7ce42aec-b319-4de9-b70c-84046d45a600
md"""
## Filters
"""

# ╔═╡ 58269465-ba8c-4840-bbfc-0a27897f3e2a
md"""
### Filter central values in nm
"""

# ╔═╡ 0b1c5662-ec4f-486e-9ee6-7fa6ba953e45
begin
	xfnm = [420.0,438.0,465.0,503.0,550.0,600.0,650.0,692.0,732.0,810.0]
	wfnm = [10.0, 24.0, 30.0, 40.0, 49.0, 52.0, 60.0, 40.0, 68.0, 10.0]
	filtnm = (center=xfnm,
	width  = wfnm,
	left = xfnm .- 0.5*wfnm,
	right = xfnm .+ 0.5*wfnm)

	println("Filter central values (nm) = ", filtnm.center, " width (nm) =", filtnm.width)
end

# ╔═╡ b07466c0-dfcd-4c10-ae86-45e71a832476
md"""
### Effect of filters in spectrum
"""

# ╔═╡ 2c75e750-854e-459f-91a6-ba135ae263cf
begin
	wr=396.0:1.0:850.0
	fg2 = lfi.LaserLab.dftof(wr, g2df, "I")
	qs = [lfi.LaserLab.qpdf(fg2, filtnm.left[l], filtnm.right[l])/filtnm.width[l] for l in 1:length(xfnm)]
	pqyd = plot(g2df.L,g2df.I, lw=2, label="G2SL")	
	pfy = plot!(collect(wr), fg2.(wr), label="")
	pqs = scatter!(xfnm, qs, label="Filters")
	plot!(xfnm, qs, lw=2, label="")
end


# ╔═╡ c892d4f2-2678-41eb-8724-6d366178f491
md"""
- The plot shows the expected discretized distribution of G2SL passed by the filters of the laser setup.  
"""

# ╔═╡ eab79cba-ca3d-40d8-9961-257e711bb9ae
begin
	qwl = lfi.LaserLab.qpdf(fg2, 0.0, 850.0)
	qflt = [lfi.LaserLab.qpdf(fg2, filtnm.left[l], filtnm.right[l]) for l in 1:length(xfnm)]
	qx = qflt ./qwl
	scatter(xfnm, qx, label="Fraction of total charge per filter", legend=:topleft)
end

# ╔═╡ e65eea70-46b3-4852-85d1-5edef9b21b37
md"""
- Plot shows the fraction of charge expected in each filter bin.
"""

# ╔═╡ c9aaf1cc-80c4-475b-8a81-e00918d91b1e
md"""
# Analysis for a single point
"""

# ╔═╡ b98ad447-6055-46e5-bb4f-9e67f9c3176a
md"""
## Parameters of the CCD
"""

# ╔═╡ 88853edb-dc1f-4e7a-a0ba-1868276d1ada
begin
	texp = 10.0
	adctopes = 0.48
	pixelsize = 16.0
	pedestal = 100.0 * pixelsize
	dcperpixel = 0.06 * pixelsize * texp
	dcperpixelc = 0.06 * pixelsize * texp / adctopes
	enoise = 1.0 * pixelsize
	enoisec = enoise / adctopes
	
	md"""
- Time of exposition (in sec) = $texp
- adc to photoelect = $adctopes
- pixelsize = $pixelsize
- pedestal = $pedestal
- dark current per pixel in e- = $dcperpixel
- dark current + noise per pixel in e- = $(dcperpixel + enoise)
- dark current in counts $dcperpixelc 
- dark current + noise per pixel in e- = $(dcperpixelc + enoisec)
	
"""
end

# ╔═╡ 57d96432-4318-4291-8255-bfa5d6d3635c
md"""
## Define input and output directories
"""

# ╔═╡ 0b4d2c08-a677-492c-b5da-982d3d5096fc
begin
	odir    = "/Users/jjgomezcadenas/LaserLab/Proyectos/pdata/CMOS2/ITO"
    cmdir   ="/Users/jjgomezcadenas/LaserLab/Proyectos/data/CMOS2/ITO"
	
	dfiles  = "*.csv"
	dplots  = "*.png"
	dmtype  = ["Imag","Dark"]
	rep     = "1"  # always take repetition number 1 if there is more than 1
	drkpnt  = "Dark"
	imgp    = "ExpoT_10s_Imag_1"
	darkp   = "ExpoT_10s_dark"
	md"""
	CMOS dir = $cmdir
	"""
end

# ╔═╡ 87519878-fb5d-405c-9222-a71872216ce7
md"""
	#### Experiment and run 
	"""

# ╔═╡ 308649b5-5c65-40dd-bc66-5b0273648341
md"""
## Algorithm
- Measure the dark current (filter and point independent)
- Select the spot (filter independent)
- For each point, measure the response for all the filters
"""

# ╔═╡ 4218a405-5cb5-464f-9ce1-5d23daeabbef
md"""
### Dark current 
- Dark folder contains measurements of the dark current for each filter. 
"""

# ╔═╡ 347a0f01-fbee-4195-b3a3-55a29285298d
md"""
The measurements above:
- Show a peak ~1600 (this is the pedestal that corresponds to the logical zero).
- The peak has an rms (std) which is defined by the DC.
- The distribution has long tails (noisy pixels).
"""

# ╔═╡ a0b9922e-189f-445e-8508-130e8e561aba
md"""
In order to supress the DK, a cutoff is needed. It can be set in terms of 
nsigma x std, where std refers to the standard deviation of the dc and nsigma is set by hand
"""

# ╔═╡ 0d846a7e-cd53-40a9-8fd5-c5be630790bb
md""" Select nsigma cut: $(@bind spsigma NumberField(1.0:100.0, default=5.0))"""

# ╔═╡ b270c34c-177b-41ba-8024-56576770b45c
md"""
### Select ROI
- Folder Filter1 contains the measurements carried out with unfiltered light. Those measurements are used to focus the beam and define the spot. 

- The algorithm to find the spot (for each point) is as follows:
	1. Compute the edge of the image (suing a sujoy algorithm, which finds the gradient and binarizes it)
	2. Find the spot delimited by the edge. The simple version implemented here finds a square spot defined between points (xmin, ymin) and (xmax, ymax) in the edge.
	
"""

# ╔═╡ ecbc5d11-3397-4495-a21f-fa7151dabcd1
md"""
#### Find the edge of the image using sujoy algorithm.
- **iedge** is a vector of tuples which contains the coordinates (indexes of the matrix) of the edge. The first coordinate runs through the y axis and the second through the x axys. 
- **yedge** and **xedge** are vectors of coordinates corresponding to the first and second coordinates of the edge
- **medge** is a matrix in which each column correspond to a data point (y, x). The matrix is used by the clusterization algorithm

The image show the edge.
"""

# ╔═╡ 19cd1ee2-f0a4-4620-89c3-e916c4551246
md"""
The image shows the edge as a scatter plot
"""

# ╔═╡ 4615ac40-0164-4297-9aad-66d39d289d15
md"""
#### Clusterize the edge

##### Customize DBSCAN algorithm 
"""

# ╔═╡ c156b901-e9e2-4576-858f-dc3aa9ae65ee
md""" Select clustering radius: $(@bind crad NumberField(0:100, default=10))"""

# ╔═╡ 8774dd9a-98d3-432e-9ff6-51190e6d4326
md""" Select min number of neighbors: $(@bind nmin NumberField(1:10, default=5))"""

# ╔═╡ 455952f8-9cf4-484a-beef-1fc2810e3b89
md""" Select min cluster size: $(@bind csize NumberField(1:100, default=10))"""

# ╔═╡ 607b0b63-f80e-4d57-bd2f-ccdf74a9af3b
md"""
#### Cluster selection 
"""

# ╔═╡ ca926619-b396-416e-9f72-4c3254f94f80
md"""
#### Find the ROI
"""

# ╔═╡ ff856515-001c-4251-ae41-a763171949e5
#indroi = indxroi(xc1, yc1, ecorner)

# ╔═╡ ba8facfa-970e-4265-95e5-58bbe27f273d
#xxe, yye = edge_to_vct(iedge, indroi)

# ╔═╡ cecdf185-4a1b-481e-8bc7-a8c4bb7d4990
md"""
##### Plot the signal in the ROI, subtracting the average value of dk
"""

# ╔═╡ a5dc8f3a-420b-4676-93e2-b6d947f26d4c
md"""
### Analysis for Filters 2-11
"""

# ╔═╡ e9e8e58a-e1db-49f1-8429-420271fb1852
md"""
## Spectrum using ROI
"""

# ╔═╡ 54bd1f6c-2b10-47a1-838f-b428fe6b7635
md""" Check to compute sum using full ROI: $(@bind zroi CheckBox())"""

# ╔═╡ a48af8f4-4ed2-45cd-b4e8-9b3106c885f3
md"""
# Analysis for al points
"""

# ╔═╡ 1794afb6-6ef0-46d6-b182-d54362b9a07d
md""" Check to carry analysis for all points: $(@bind zrec CheckBox())"""

# ╔═╡ be87f23e-18b0-4927-ba6d-902180f05489
function signalnorm(sumf::Vector{Float64}, sumtot::Float64, filtnm; scale=1e+4)
	scale*(sumf/sumtot) ./filtnm.width
end

# ╔═╡ b2d6792d-9de4-4d7e-beb3-879459d3709c
if zrec
	#PFLT = sum_allpoints(cmdir, string(sexp), string(srun), 
		#                 string.(fpoints), string.(xfn),
		#	             dkavg, ctx, crad, nmin, csize, scl)
end

# ╔═╡ 50e69dc7-0027-4d0c-bfa8-f6778a5754f3
if zrec
	#PXFLT = []
	#for (i, plft) in enumerate(PFLT)
	#	lbl = string("Point",i)
	#	pfxx = plot(filtnm.center, plft, lw=2, label=lbl, legend=:topright, title="")
	#	scatter!(filtnm.center, plft, label="")
	#	xlabel!("λ (nm)")
	#	ylabel!("counts")
	#	push!(PXFLT, pfxx)
	#end
	#plot(size=(1050,1050), PXFLT..., layout=(2,2), titlefontsize=8)
end

# ╔═╡ 82ddd81b-8aea-4711-97d9-a645121786f8
md""" Check to read and plot data for all points: $(@bind zread CheckBox())"""

# ╔═╡ 155d5066-935b-4643-8aad-f4c635aa7eec
if zread
	md""" Select nnumber of plots in x and y: 
	- select nx $(@bind nx NumberField(1:10, default=3))
	- select ny $(@bind ny NumberField(1:10, default=3))
	"""
end

# ╔═╡ 5a88cb1e-47b2-45cc-965c-2af9a45e72f5
md"""
# On image representation

Images are stored in Matrices (e.g, a 512 x 512 matrix in our case). One has to be careful with the index convention to handle them:

1. First index runs over rows (thus it corresponds to y), second over columns (x)
2. In the image, pixel (1,1) corresponds to the bottom-left of the image, while in an ordinary array arrangement, one expects pixel (1,1) in the upper top. One can visualize both reference systems are rotated 90 degrees with respect each other

One has to be careful with this shift when handly images, as illustrated below
"""

# ╔═╡ 7d3bd063-e821-4bd2-b375-5b0989e49270
function test_matrix(xs, ys, irng, jrng)
    tmx = zeros(xs, ys)
	indx = []
    for i in irng
        for j in jrng
            tmx[i,j] = 1.0
			push!(indx, (i,j))
        end
    end
    tmx, indx
end


# ╔═╡ 7083fcc2-d2f0-44fa-b85e-0000bb100c0a
function transform_indx(indx, size)
	[(size-I[1] + 1, I[2]) for I in indx]
end

# ╔═╡ ed902bce-fb55-4b96-a0ea-cf335e529531
function transform_mtrx(mtrx, xsz)
	tedge = lfi.LaserLab.indx_from_edge(mtrx, xsz)
	tte = transform_indx(tedge, xsz)
	
	tmx = zeros(xsz, xsz)

    for I in tte
    	tmx[I[1],I[2]] = 1.0
    end
    tmx
end

# ╔═╡ 46f627cd-4166-40c7-8330-d72ac586d3c0
begin
	xsz=6
	irng = 2:4
	jrng = 3:4
end

# ╔═╡ d31856a8-4f2b-4f8b-ab9c-20b4cbb643ea
md"""
Create an image of sz =$xsz filling it with ranges $irng and $jrng
"""

# ╔═╡ 61198991-dad0-44e6-9715-a599d4dac0c9
tmrx, indx = test_matrix(xsz, xsz, irng, jrng);

# ╔═╡ 7ea8537e-3951-4453-8140-7e2f31f5d900
md"""
This is the matrix
"""

# ╔═╡ 6d2f6df0-c7b3-4f05-9c07-4f9690372c19
tmrx

# ╔═╡ 2ff4a48e-9621-4375-9b05-ab7424ba98fa
md"""
And these are the indexes in which the value of the pixel equals one in "upper top" reference system 
"""

# ╔═╡ 9a153985-5b8d-4686-99e6-a8038965dddd
indx

# ╔═╡ 087751d5-087a-4a88-9dc1-a599fbbfcae3
md"""
Function **indx\_from\_edge** find those indexes
"""

# ╔═╡ 9752523c-7f50-45cb-9341-c0d59e35f772
tedge = lfi.LaserLab.indx_from_edge(tmrx, xsz)

# ╔═╡ 436dbab0-5078-46dc-be07-9f04cdf4c46a
function edge_corners2(iedge)
	indxmin = minimum([ii[1] for ii in tedge])
	lindxmin = [ii for ii in tedge if ii[1] == indxmin ]
	indymin = minimum([ii[2] for ii in lindxmin])
	indxmax = maximum([ii[1] for ii in tedge])
	lindxmax = [ii for ii in tedge if ii[1] == indxmax ]
	indymax = maximum([ii[2] for ii in lindxmax])

	(minvx=(indxmin,indymin ), maxvx = (indxmax, indymax))
end

# ╔═╡ 867d2595-632c-477b-89b7-85a0dd8a8941
md"""
The edge corners (upper left and bottom right) in the top-left convention
"""

# ╔═╡ 47161d36-4c22-4ca0-a580-24902fc4e1c4
edge_corners2(tedge)

# ╔═╡ dd1a6f48-cba1-4896-91c9-dfa0ee51b765
md"""
Notice that in this convention minvx is upper left, maxvx is bottom right and the ROI is defined by indexes smaller than minvx and larger than maxvx
"""

# ╔═╡ 23c3ee67-80e1-48d1-8296-05c814d30c76
md"""
Function **transform\_indx** transforms to bottom-left reference system
"""

# ╔═╡ b1c38ab5-0018-4200-a5b2-f8a7b24bc129
tte = transform_indx(indx, xsz)

# ╔═╡ 6541aa1a-cbcb-47ec-baed-62c58f4f8ae3
md"""
This is the matrix in the bottom-left reference system
"""

# ╔═╡ 07e6d6a8-4556-423d-8600-281750f04707
t2mrx = transform_mtrx(tmrx, xsz)

# ╔═╡ 72f975a9-f9ad-414b-8aee-4f7820fcf3de
md"""
Notice that in this convention minvx is bottom left, maxvx is upper right and the ROI is defined by indexes larger than minvx and smaller than maxvx
"""

# ╔═╡ e82ede75-0b64-4eb0-8130-748cfdf69945
#imgbox(tmrx, (3,3), 1; isize=xsz)

# ╔═╡ c1b2cf36-d72e-4348-911e-8e44f1834ae4
md"""
# Functions
"""

# ╔═╡ 56771073-b0b3-47bf-8578-2dd11a59a9b2
md"""
## Data structures
"""

# ╔═╡ f2714261-7bb0-47d7-8aac-c16bb5d1f891
struct Setup
    cmdir::String
    series::String
    measurement::String
    point::String
    filter::String
    rep::String 
end


# ╔═╡ 81a2bc2a-2821-4aa6-86e9-97eedc0bc51c
md"""
## Management and intendence functions
"""

# ╔═╡ 0d5a7021-6072-464e-836e-f05b1e178b80
"""
Create a dir if it doesn't exist yet
"""
function create_dir!(dir)
	if isdir(dir) == false
		mkdir(dir)
	end
end

# ╔═╡ a759ecf7-7373-46bf-ab15-df39cbfa6814
"""
Create output directories i they don't exist
"""
function output_dirs!(sexp)
	namex = split(sexp,"_")
	create_dir!(joinpath(odir, sexp))
	csvdir = joinpath(odir, sexp,"csv")
	pngdir = joinpath(odir, sexp, "png")
	create_dir!(csvdir)
	create_dir!(pngdir)
	csvdir, pngdir
end

# ╔═╡ bc421f3e-e091-4f83-bc47-ab8572570e1b
"""
Returns dirs below cmdir (experiments by convention)
"""
function select_exp(cmdir)
	readdir(cmdir)
	lfi.LaserLab.getbolddirs(cmdir)
end

# ╔═╡ e87f48e3-5e5a-44d5-83de-c520e522e33a
let
	dirs = select_exp(cmdir)
	md""" Select experiment : $(@bind sexp Select(dirs))"""
end

# ╔═╡ 161b12b6-88a0-4d4d-bfc5-01310534cbdc
begin 
	csvdir, pngdir = output_dirs!(sexp)
	md"""
	#### Output dirs
	- csv dir = $csvdir
	- png dir = $pngdir
	"""
end

# ╔═╡ 2498aa42-2c24-47e3-bf5b-647377af0dbc
"""
Returns dirs defined by cmdir and sexp (run by convention)
"""
function select_run(cmdir, sexp)
	namex = split(sexp,"_")
	path = joinpath(cmdir,sexp)
	readdir(path)
	lfi.LaserLab.getbolddirs(path)
end

# ╔═╡ 50ea2ecc-970f-4630-8c7e-acf5e69cc4c9
let
	dirs = select_run(cmdir, sexp)
	md""" Select run : $(@bind srun Select(dirs))"""
end

# ╔═╡ f26bb6e0-45ac-4419-bcb2-46e2cac1f75b
begin
	md"""
	
	- experiment = $sexp
	- run = $srun
	"""
	#sexp = "AYN_05_Ba"
	#srun = "AYN_05_Ba2_45x_G2SL_rep2_P9_220901"
end

# ╔═╡ afeb42ca-b462-4320-b364-98a0b4730e33
"""
Returns dirs defined by cmdir, sexp and srun (point by convention)
"""
function select_point(cmdir, sexp, srun)
	path = joinpath(cmdir,sexp,srun)
	pdirs = readdir(path)
	pdirs = lfi.LaserLab.getbolddirs(path)
	points = [split(pd, "Point")[2] for pd in pdirs if findfirst("Point", pd) != nothing]
	spdirs = sort(parse.(Int64, points))
	[string("Point", string(i)) for i in spdirs]
end

# ╔═╡ d389f99a-14c2-408f-ad7b-838e00225357
begin
	sspdirs = select_point(cmdir, sexp, srun)
	#md""" Select point : $(@bind spoint Select(sspdirs))"""
end

# ╔═╡ 46b1d54e-4ebf-45a2-b137-9690b7a51d38
"""
Given a directory tree defined by root directory (cmdir), experiment (sexp),
run (srun) and point (spoint) returns a list of files (of type dfiles) found
in the directory

"""
function select_files(cmdir::String, sexp::String,srun::String, spoint::String, 
	                  dfiles="*.csv")
	path = joinpath(cmdir,sexp,srun, spoint)
	readdir(path)
	xfiles = Glob.glob(dfiles, path)
	nxfiles = string.([split(f,"/")[end] for f in xfiles])
	xfiles, nxfiles
end


# ╔═╡ 077d7e4e-1b94-4b30-a70a-f3b5d3a6fc46
ff1, nff1 = select_files(cmdir,string(sexp),string(srun), "Filter1")

# ╔═╡ c69c8b9d-50bf-46ce-8614-1fee1661e424
fpoints = [split(pd, "_")[1] for pd in nff1] 

# ╔═╡ aa5bdf46-a4db-4ff7-9e51-8d55bc6c203d
begin
	md""" #### Select point  : $(@bind spointf1 Select(fpoints))"""
end

# ╔═╡ 95100b23-f017-4861-93c1-4adc571e467e
md"""
#### Image for $spointf1
"""

# ╔═╡ e38d2a94-008e-46db-a35e-176bb5ae297a
fpoints

# ╔═╡ 24d380e0-18d2-48ae-a2ea-3cb42d5d75e0
fpoints

# ╔═╡ e8b6d609-7357-4996-bda5-f1119b4b0263
f1img = lfi.LaserLab.select_image(ff1, string(spointf1));

# ╔═╡ 8c911d89-107b-4445-9bda-e7e9b50ff051
heatmap(f1img.img)

# ╔═╡ 479f8f86-372c-4b91-9f73-e57a85d3d194
begin
	spoint = spointf1
	xfiles, nxfiles  = select_files(cmdir,string(sexp),string(srun), string(spoint))
	#xfdrk, nxdrk     = select_files(cmdir,sexp,srun, drkpnt)
end

# ╔═╡ f13173e1-088c-4f5c-ae6a-f9aebcd6bc57
begin
	xfa = lfi.LaserLab.findpattern(nxfiles, "Filter")
	xfb = sort(parse.(Int64, xfa))
	xfn = string.(xfb)
	md""" ##### Select filter: $(@bind sfn Select(xfn))"""
end

# ╔═╡ ec95a04a-bda9-429b-92f2-c79f28a322f0
"""
Given a list **xnames** (full path to files) returns the filter names
"""
function flt_names(xnames::Vector{String})
		fxnm = [split(fn, "/")[end] for fn in xnames]
		fxnb = [split(fn, "_")[2] for fn in fxnm]
		fxint = sort([parse(Int64, ff) for ff in fxnb])
		[string("Filter_", string(i), "_") for i in fxint]
	end

# ╔═╡ d71e7c72-d7ae-4915-bc36-aed347d09450
"""
Return the mean and std of a matrix if the elements are in an interval
"""
function meanstd_interval(img::Matrix{Float64}, int::Tuple{Float64, Float64})
	
	sgn = Vector{Float64}(undef,0)

	for i in 1:size(img)[1]
		for j in 1:size(img)[2]
			if img[i,j] > int[1] && img[i,j] < int[2]
    			push!(sgn,img[i,j])
			end
		end
	end
	mean(sgn), std(sgn)
end

# ╔═╡ 5f89baa6-d8f4-4d9c-a625-81cf9375f89c
"""
Return the sum of a matrix above threshold
"""
function sum_thr(img::Matrix{Float64}, thr::Float64)
	sumx = 0.0
	for i in 1:size(img)[1]
		for j in 1:size(img)[2]
			if img[i,j] > thr
    			sumx += img[i,j] 
			end
		end
	end
	sumx
end

# ╔═╡ 85d6c31f-2fc8-486b-98e9-e5b9e134df6a
"""
Return the sum of a matrix above threshold
"""
function sum_interval(img::Matrix{Float64}, int::Tuple{Float64, Float64})
	sumx = 0.0
	for i in 1:size(img)[1]
		for j in 1:size(img)[2]
			if img[i,j] > int[1] && img[i,j] < int[2]
    			sumx += img[i,j] 
			end
		end
	end
	sumx
end

# ╔═╡ 03b663a3-8e01-4720-ad55-b9085ee5115d
"""
Returns the regularised sum in the ROI. 

- roidks is a matrix with active pixels
- dkavg is the average of the dc which must be substracte
- ctx is the effective cut to suppress dark current
- nsigmaT is the interval to recompute mean and std of the roi suppressing tails
- nsigmaS is the right-wig of the distribution to sum (mean + std* nsigmaS)
"""
function roisum(roidks::Matrix{Float64}, 
				dkavg::Float64, ctx::Float64, 
                nsigmaT::Float64 = 5.0, nsigmaS::Float64 = 2.5)

	roi = roidks .- dkavg
	meanf = mean(roi)
	stdf  = std(roi)
	meanflt, stdflt = meanstd_interval(roi, (-nsigmaT*stdf,nsigmaT*stdf))
	sum_interval(roi,(ctx, meanflt + nsigmaS*stdflt))
end

# ╔═╡ 84170688-fbc6-4676-84f1-126ecce4f5f2
"""
Indexes are colum wise in julia, so (x,y) coordinates corresponde to (second, first) index 
"""
function get_coord_from_indx(xyindx)
	xyindx[2], xyindx[1]
end

# ╔═╡ 0ad139e9-85bc-421b-bdd7-e61711d47454
"""
Return a tuple ((x1,x1)...(xn,yn)) with the coordinates of the points in edge
"""
function indx_from_edge(iedge::Matrix{Float64}, isize=512)
	indx = Vector{Tuple{Int, Int}}(undef, 0)
	for i in 1:isize
		for j in 1:isize
			if iedge[i,j] == 1
				push!(indx,(i,j))
			end
		end
	end
	indx
end

# ╔═╡ 65402da4-601a-4766-ba2c-6735f201ef6a
"""
Take the vector of indexes defining the edge (iedge) and return a matrix m(2, n),
where n is the size of iedge and two vectors x(n), y(n). These objects are needed
for DBDSCAN clustering 
"""
function edge_to_mtrx(iedge::Vector{Tuple{Int, Int}})
	nedge = length(iedge)
	medge = zeros(2, nedge )
	xedge = zeros(nedge)
	yedge = zeros(nedge)
	for i in 1:nedge
		medge[1, i] = iedge[i][1]
		medge[2, i] = iedge[i][2]
		xedge[i] = iedge[i][1]
		yedge[i] = iedge[i][2]
	end
	medge, xedge, yedge
end

# ╔═╡ 9b118763-2739-4535-99c3-da6245ba1eae
begin
	img_edge = Float64.(lfi.LaserLab.sujoy(f1img.imgn, four_connectivity=true))
	img_edgeb = Float64.(lfi.LaserLab.binarize(img_edge, Otsu()))
	mxedge, imxedge =findmax(img_edge)
	vxmaxed, vymaxed = get_coord_from_indx(imxedge)
	iedge = indx_from_edge(img_edgeb)  #indexed of the edge
	medge, xedge, yedge = edge_to_mtrx(iedge)  # edge expressed as matrix, x,y
	heatmap(img_edgeb)
	#mosaicview(Gray.(f1img.imgn), Gray.(img_edgeb); nrow = 1)
	#scatter!([vxmaxed], [vymaxed], label="edge max",markersize=3)
end

# ╔═╡ 1154dba4-fc68-48c4-ab55-c4a5687d0547
scatter([yedge], [xedge], label="edge",markersize=2)

# ╔═╡ ddf630eb-6be5-422e-9db5-b1a4348adf01
clusters = dbscan(medge, crad, min_neighbors = nmin, min_cluster_size = csize)

# ╔═╡ 8414f6ab-c4f4-4f01-8082-80468ac4e04b
md"""
- DBSCAN found $(length(clusters)) clusters
"""

# ╔═╡ 0ca0299c-6c33-4e35-9642-66e77b1cedba
begin
	fcl = 1:length(clusters)
	md""" ##### Select cluster  : $(@bind scl Select(fcl))"""
end

# ╔═╡ e956095a-e468-4456-9659-b1acd6bd507d
"""
Given vector xedge and yedge containing (x,y) coordinates of the edge and the
clusters object returned by DBDSCAN, returns vectors of cluster points, cx, yc (for cluster number nc)

"""
function getclusters(iedge::Vector{Tuple{Int, Int}}, 
	                 clusters::Vector{DbscanCluster}, nc::Int64)
	cs = clusters[nc].size
	ci = clusters[nc].core_indices
	cb = clusters[nc].boundary_indices
	
	xc = Vector{Int64}(undef, cs)
	yc = Vector{Int64}(undef, cs)
	ii=1
	for (i, indx) in enumerate(ci)
		xc[ii] = iedge[indx][1]
		yc[ii] = iedge[indx][2]
		ii+=1
	end
	for (i, indx) in enumerate(cb)
		xc[ii] = iedge[indx][1]
		yc[ii] = iedge[indx][2]
		ii+=1
	end
	xc,yc
end

# ╔═╡ cf1adbfc-e494-44e2-a6d6-7b89fb9a6be6
xc1, yc1 = getclusters(iedge, clusters, scl);

# ╔═╡ 8946f9d0-7108-4af5-90a7-8f338f59009b
md"""
- Cluster $scl containd $(length(xc1)) points
"""

# ╔═╡ 58bdc137-fc29-4f5e-9ab0-d6f959cfa71c
"""
Given a vector defining an edge, of two vectors of coordinates, return the topleft and bottomright corners
"""
function find_edge_corners(xc1::Vector{Int}, yc1::Vector{Int})
	zcorner = zip(xc1, yc1)
	topy = maximum([ii[1] for ii in zcorner])
	lfty = minimum([ii[2] for ii in zcorner])
	boty = minimum([ii[1] for ii in zcorner])
	rgty = maximum([ii[2] for ii in zcorner])

	topy, lfty, boty, rgty
	(topleft=(topy, lfty), botright=(boty, rgty))
end


# ╔═╡ e3b52355-a912-42ff-ac39-08a79a1ccee5
ecorner = find_edge_corners(xc1, yc1);

# ╔═╡ ca7eb636-f817-4a1a-8b25-7e5ff2d3d0e5
md"""
##### corners of selected cluster:
- top-left =$(ecorner.topleft)
- bottom-right =$(ecorner.botright)
"""

# ╔═╡ 1658ecd9-8949-4052-9874-a31248c45821
"""
Given point pt, return the edge corner of cluster scl
crad, nmin and csize are parameters of the DBDSCAN algorithm 
"""
function select_edge_corners(cmdir::String, sexp::String, srun::String, pt::String,
                             crad::Int64, nmin::Int64, csize::Int64, scl::Int64)

		# get image
		flt1, _   = select_files(cmdir, sexp, srun, "Filter1")
		f1img     = lfi.LaserLab.select_image(flt1, pt)

		#image edge
		img_edge  = Float64.(lfi.LaserLab.sujoy(f1img.imgn, four_connectivity=true))
		img_edgeb = Float64.(lfi.LaserLab.binarize(img_edge, Otsu()))
		iedge     = indx_from_edge(img_edgeb)  
		medge, xedge, yedge = edge_to_mtrx(iedge)

		#find clusters and return edge
		clusters = dbscan(medge, crad, min_neighbors = nmin, min_cluster_size = csize)
		xc1, yc1 = getclusters(iedge, clusters, scl)
		find_edge_corners(xc1, yc1)
	end

# ╔═╡ fa052909-2ed8-4d90-8144-73147682df0e
"""
Returns true if coordinates x,y are inside the box defined by the corners of the edge
"""
function inedgex(vxe::NamedTuple{(:topleft, :botright)}, i::Int64,j::Int64)
	if i <= vxe.topleft[1] &&  j >= vxe.topleft[2]   
		if i >= vxe.botright[1] && j <= vxe.botright[2] 
			return true
		end
	end
	return false
end

# ╔═╡ 369dcaf2-ce5d-4641-a8cb-274161749c19
"""
Given an image (img) and the vertex of the image edge (ecorn), returns the image
in the roi defined by a square box around the edge. The returned matrix is sparse
"""
function imgroix(img::Matrix{Float64}, 
	             ecorn::NamedTuple{(:topleft, :botright)}; isize=512, prtlv=0)

	#vxe = lfi.LaserLab.indx_from_edge(img, isize)
	#ecorn = lfi.LaserLab.edge_corners(vxe)
	
	sx = abs(ecorn.botright[2] - ecorn.topleft[2]) + 1
	sy = abs(ecorn.topleft[1] - ecorn.botright[1]) + 1

	if prtlv == 1
		#println("vxe =", vxe)
		println("ecorn = ", ecorn)
		println("sx = ", sx, " sy = ", sy)
	end 

	roi = zeros(sx, sy)
	iroi = Vector{Tuple{Int, Int}}(undef, 0)
	
	ix = 0
	jx = 1
	for il in 1:isize
		for jl in 1:isize
			if prtlv == 3
				println("i = ", il, " j = ", jl, " ic = ", " in edge? =", 
				     inedgex(ecorn, il,jl))
			end
			if inedgex(ecorn, il,jl) 
				if ix < sy
					ix+=1
				else
					ix = 1
					jx+=1
				end
				if prtlv == 2
					println("il = ", il, " jl= ", jl, " ix = ", ix, " jx = ", jx)
				end
				roi[jx,ix] = img[il,jl]
				push!(iroi, (il,jl))
			end
		end
	end
	roi, iroi
end

# ╔═╡ 1968c0dd-8cf9-476a-84a8-b92b37ac02ad
roixx, iroixx = imgroix(f1img.img, ecorner, prtlv=0);

# ╔═╡ a1817474-d089-4245-9263-70314969b43e
md"""
- ROI has dimensions: $(size(roixx)) 
"""

# ╔═╡ 205fd30f-6a4e-473f-9eb5-e826b8c480b1
"""
Returns the sum of the signal in the ROI

The sum is computed between the threshold to suppress dark current (ctx) and
a maximum value of nsigmaS * std, where std is computed suppressing previously tails. 
This guarantees that the impact of hot pixels is mimimised. 

"""
function signal_roi(xfiles::Vector{String}, xfn::Vector{String}, 
					ecorn::NamedTuple{(:topleft, :botright)},
	                dkavg::Float64, ctx::Float64, 
	                nsigmaT::Float64 = 5.0, nsigmaS::Float64 = 2.5)

	function signal_flt(flt::String)
		fltimg = lfi.LaserLab.select_image(xfiles, flt)
		#roisum(fltimg.img, dkavg, ctx)
		fltcimg = fltimg.img .- dkavg
		fltroi, iroiflt = imgroix(fltcimg, ecorn)
		meanfltt = mean(fltroi)
		stdfltt = std(fltroi)
		meanflt, stdflt = meanstd_interval(fltroi, 
			                               (-nsigmaT*stdfltt,nsigmaT*stdfltt))
		sum_interval(fltroi,(ctx, meanflt + nsigmaS*stdflt))
	end
	
	map(flt->signal_flt(flt), xfn)

end

# ╔═╡ f649a197-c7a8-407d-9e11-6fd811c9d375
"""
Returns the sum of the signal in the ROI for all points

"""
function signal_roi_allpoints(cmdir::String, sexp::String, srun::String,
	                          xpt::Vector{String}, xfn::Vector{String}, 
	                          dkavg::Float64, ctx::Float64, 
                              crad::Int64, nmin::Int64, csize::Int64, scl::Int64,
	                          nsigmaT::Float64 = 5.0, nsigmaS::Float64 = 2.5)

	function signal_pt(pt::String)
		ecorn  = select_edge_corners(cmdir, sexp, srun, pt,
                                      crad, nmin, csize, scl)
		xfiles, _ = select_files(cmdir,sexp,srun, pt)
		signal_roi(xfiles, xfn, ecorn, dkavg, ctx, nsigmaT, nsigmaS)
	end
	
	map(pt->signal_pt(pt), xpt)
end

# ╔═╡ dfda19e3-e772-4458-b29d-2885a020f77b
"""
Given an idege vector and a set of indexes return coordiantes
"""
function xy_from_tuplelist(iedge::Vector{Tuple{Int, Int}})
	xedge = Vector{Int}(undef, length(iedge))
	yedge = Vector{Int}(undef, length(iedge))
	for ii in 1:length(iedge)
		xedge[ii] = iedge[ii][1]
		yedge[ii] = iedge[ii][2]
	end
	xedge, yedge
end

# ╔═╡ 939c1d23-35d3-41bf-a854-395457e9ad59
begin
	xxe, yye = xy_from_tuplelist(iroixx)
	
	exminc, eyminc = get_coord_from_indx(ecorner.topleft)
	exmaxc, eymaxc = get_coord_from_indx(ecorner.botright)
	
	scatter([yedge], [xedge], label="edge",markersize=2)
	scatter!([yc1], [xc1], label="cluster $scl",markersize=2)
	
	scatter!([ecorner.topleft[2]], [ecorner.topleft[1]], 
		     label="top left clust $scl",markersize=4)
	
	scatter!([ecorner.botright[2]], [ecorner.botright[1]], 
	         label="bottom right clust $scl",markersize=4)
	hline!([ecorner.topleft[1]], label="")
	vline!([ecorner.topleft[2]], label="")
	hline!([ecorner.botright[1]], label="")
	vline!([ecorner.botright[2]], label="")
	scatter!([yc1], [xc1], label="ROI edge for clust $scl",markersize=4)
	scatter!([yye], [xxe], label="ROI points for clust $scl",markersize=1)
	
	
end

# ╔═╡ 92a5b0b7-5dbb-4e98-a33e-bef1f6992b40
md"""
## Manipulation of images
"""

# ╔═╡ aa088449-b0e0-414f-9a96-1a2c5d6656ff
"""
Return the mean and std of a matrix above threshold
"""
function meanstd_thr(img::Matrix{Float64}, thr::Float64)
	sgn = Vector{Float64}[]
	for i in 1:size(img)[1]
		for j in 1:size(img)[2]
			if img[i,j] > thr
    			push!(sgn,img[i,j])
			end
		end
	end
	mean(sgn), std(sgn)
end

# ╔═╡ 8496f1e5-ad5a-4469-851e-0795c8cf6c6b
"""
Returns the indexes of the roi contained between the corners topleft and bottomright
"""
function indxroi(xc1::Vector{Int}, yc1::Vector{Int}, 
	             ecorn::NamedTuple{(:topleft, :botright)})
	
	zcorner = zip(xc1, yc1)
	findall(i->(i==true), [II[1] <= ecorn.topleft[1] && II[2] >= ecorn.topleft[2] && II[1] >= ecorn.botright[1] && II[2] <= ecorn.botright[2] for II in zcorner ])
end


# ╔═╡ 0bf2d20e-dff6-48c2-b059-38b87ed46bb6
"""
Given an idege vector and a set of indexes return coordiantes
"""
function edge_to_vct(iedge::Vector{Tuple{Int, Int}}, indx::Vector{Int})
	xedge = Vector{Int64}(undef, length(indx))
	yedge = Vector{Int64}(undef, length(indx))
	for ii in 1:length(indx)
		xedge[ii] = iedge[ii][1]
		yedge[ii] = iedge[ii][2]
	end
	xedge, yedge
end



# ╔═╡ c3cb6332-0b65-4a49-aa63-a846df3a5277
"""
Given a roi matrix, return the indexes
"""
function indx_from_roimtrx(imgroi)
	INDX = []
	sz= size(imgroi)
	for i in 1:sz[1]
		for j in 1:sz[2]
			if imgroi[i,j] != 0.0 
				push!(INDX, (i,j))
			end
		end
	end
	INDX
end

# ╔═╡ ae5cf21c-7656-4f9f-bd04-73e0e0d8fbee
"""
Dumps a matrix into a vector
"""
function mtrxtovct(mtrx::Matrix{<:Real})
	sz = size(mtrx)
	vx = zeros(sz[1]*sz[2])
	ii = 1
	for i in 1:sz[1]
		for j in 1:sz[2]
			vx[ii] = mtrx[i,j]
			ii+=1
		end
	end
	vx	
end

# ╔═╡ 54664106-32d8-4ba9-b3e4-0a926a02309c
"""
Given a ROI (eg, a matrix representing a section of an image) it returns an histogram of the signal in the ROI
"""
function histo_signal(iroi::Matrix{Float64}, nbin::Int=100)
	vroi = mtrxtovct(iroi)
	mxvroi = maximum(vroi)
	mnvroi = minimum(vroi)
	lfi.LaserLab.hist1d(vroi, "signal in roi", nbin, mnvroi,mxvroi)
end



# ╔═╡ 6b6a3651-62e8-4449-9f9d-d8cf3e9ee04a
function histo_signal(iroi::Matrix{Float64}, nbin::Int, min::Float64, max::Float64)
	vroi = mtrxtovct(iroi)
	lfi.LaserLab.hist1d(vroi, "signal in roi", nbin, min,max)
end

# ╔═╡ 0ee6da42-aa03-4923-bee9-15b92b6583c5
begin
	fd1, nfd1 = select_files(cmdir,string(sexp),string(srun), "Dark")
	fltn = flt_names(fd1)
	DRK = [lfi.LaserLab.select_image(fd1, flt).img for flt in fltn];
	HDRK = [histo_signal(drk) for drk in DRK]
	PDRK = [HDRK[i][2] for i in 1:length(HDRK)]
	AVG = [mean(drk) for drk in DRK]
	STD = [std(drk) for drk in DRK]
	dkavg = mean(AVG)
	dkstd = mean(STD)
	phdrk = plot(size=(750,750), PDRK..., layout=(5,2), titlefontsize=8)
end

# ╔═╡ e39d2f51-717a-43d0-90b2-25579e625844
if zroi
	pdavg = plot(filtnm.center, AVG, lw=2, label=spointf1, title="mean")
	scatter!(filtnm.center, AVG,label="")
	ylims!(1600., 1700.0)
	xlabel!("counts")
	ylabel!("frequency")
	pdstd = plot(filtnm.center, STD, lw=2, label=spointf1, title="std")
	scatter!(filtnm.center, STD,label="")
	xlabel!("counts")
	ylabel!("frequency")
	ylims!(10., 40.0)
	plot( pdavg, pdstd, layout=(1,2), titlefontsize=8)
end

# ╔═╡ 7472f7b8-b50b-4832-b7f7-0134d1c5ed8f
begin
	dkstat = [meanstd_interval(drk .- dkavg, (-3.0* dkstd, 3.0* dkstd)) for drk in DRK]
	STDX = [dk[2] for dk in dkstat]
	dkstdx = mean(STDX)
	plot(filtnm.center, STDX, lw=2, label=spointf1, title="std DC no tail")
	scatter!(filtnm.center, STDX,label="")
	ylims!(0., 25.0)
	xlabel!("Filter (nm)")
	ylabel!("STD DC counts")
end

# ╔═╡ dbc8bed9-ba70-443c-9bd3-eaf91336a603

begin
nsigma=promsel = Float64(spsigma)
ctoff = dkstdx * nsigma	
md"""
- With nsigma = $nsigma
- Cutoff value (nsigma x dkstd) = $ctoff
"""
end

# ╔═╡ 9957adcf-48e2-4930-849e-b7601f1d8346
ctoff

# ╔═╡ 952ecc28-e2d7-4833-b820-73689a950ea6
nsigma* dkstdx

# ╔═╡ 9e4e99fe-351d-4a4c-afb2-3a9ed749209c
md""" 
STD of dark current:

- Not suppressing tails = $dkstd
- Supressing tails = $dkstdx
"""

# ╔═╡ 406cc319-e7a9-4c68-b732-774b7d1a7e59
begin
	zimg = lfi.LaserLab.select_image(xfiles, string(sfn));
	zcimg = zimg.img .- dkavg;
	roiflt, iroiflt = imgroix(zcimg, ecorner, prtlv=0)
	heatmap(zcimg)
end

# ╔═╡ 904d10aa-aea9-46e7-a558-73bcd061a0ec
begin
	HDRK2 = [histo_signal(drk .- dkavg, 100, -3.0* dkstd, 3.0* dkstd) for drk in DRK]
	PDRK2 = [HDRK2[i][2] for i in 1:length(HDRK)]
	phdrk2 = plot(size=(750,750), PDRK2..., layout=(5,2), titlefontsize=8)
end

# ╔═╡ b81df45b-e64e-4a07-982f-368ae03353c2
begin
	roidks = roixx .- dkavg
	hroi, proi = histo_signal(roidks, 50, 0.0, 5e+4)
	plot(proi)
end

# ╔═╡ 918126c6-b25f-49b5-b5a2-3c15ba1d9cfd
begin
	ctx = nsigma * dkstdx
	ntsigma = 5.0
	nssigma = 3.0
	#sumnf = roisum(roixx, dkavg, ctx)
	#sumtl = sum(roidks)

	meanwltt = mean(roidks)
	stdwltt = std(roidks)
	#meanwlt, stdwlt = meanstd_interval(roidks, (-ntsigma*stdwltt, ntsigma*stdwltt))
	meanwlt, stdwlt = meanstd_interval(roidks, (0.0, ntsigma* stdwltt))
	sumtl = sum_thr(roidks, 0.0)
	sumnf = sum_interval(roidks,(ctx, meanwlt +nssigma*stdwlt))
md"""
- Number of sigmas for dc suppression = $nsigma
- Number of sigmas to recompute mean/std no tails = $ntsigma
- Number of sigmas to sum signal = $nssigma
- std of dark current =$(round(dkstdx, sigdigits=3))

- DC cutoff = $(round(ctoff, sigdigits=4))
- signal mean (tails) =$(round(meanwltt, sigdigits=3))
- signal std (tails) =$(round(stdwltt, sigdigits=3))
- signal mean (no tails) =$(round(meanwlt, sigdigits=3))
- signal std (no tails) =$(round(stdwlt, sigdigits=3))
- Total signal (no cuts) = $(round(sumtl, sigdigits=2))
- Total signal (no dc, no tails) = $(round(sumnf, sigdigits=2))
"""
end

# ╔═╡ 3549eafe-70c4-49d5-a121-11b8d2860804
begin
	DKSUMZ = [sum_interval(drk .- dkavg, (ctx, nsigma* dkstdx)) for drk in DRK]
	#dksumz = mean(DKSUMZ)
	DKSUMT = [sum_thr(drk .- dkavg, 0.0) for drk in DRK]
	#dksumt = mean(DKSUMT)
	
	pdksumz = plot(filtnm.center, DKSUMZ, lw=2, label=spointf1, title="sum DC no tail")
	scatter!(filtnm.center, DKSUMZ,label="")
	ylims!(-1.0e3, 1.0e3)
	xlabel!("value in filter (nm)")
	ylabel!(" Sum in filter DC/tail suppressed (counts)")
	
	pdksumt = plot(filtnm.center, DKSUMT, lw=2, label=spointf1, title="sum DC int tail")
	scatter!(filtnm.center, DKSUMT,label="")
	ylims!(0.0, 1.0e7)
	xlabel!("Value in filter (nm)")
	ylabel!("SUM in filter DC (counts)")
	plot(size=(750,750), pdksumz, pdksumt, layout=(1,2), titlefontsize=8)
end

# ╔═╡ 2ad43913-db20-4c21-8a9d-cae776aa9b92
ctx

# ╔═╡ da7e09c2-cf28-414d-bca2-3b9a35275857
begin
	nflt = parse(Int64, sfn)
	meanfltt = mean(roiflt)
	stdfltt = std(roiflt)
	meanflt, stdflt = meanstd_interval(roiflt, (-5.0*stdfltt, 5.0*stdfltt))
	sumtt = sum_thr(roiflt, 0.0)
	sumt = sum_interval(roiflt,(ctx, meanflt + 2.5*stdflt))

	sumtf = sumt/filtnm.width[nflt]
md"""
##### Image for filter $sfn in ROI
- average: $(round(meanfltt, sigdigits=2))    
- std = $(round(stdfltt,sigdigits=3))
- average (tail suppressed): $(round(meanflt, sigdigits=2))    
- std (tail suppressed) = $(round(stdflt,sigdigits=3))

- Total charge (no interval) : $(round(sumtt, sigdigits=3))

- lower cutoff for sum = $(round(ctx, sigdigits=3))
- upper cutoff for sum= $(round(meanflt + 2.5*stdflt, sigdigits=3))
- Total charge (interval) = $(round(sumt, sigdigits=3))

- Filter width = $(filtnm.width[nflt])
- Charge/nm = $(round(sumtf, sigdigits=3))
"""
end

# ╔═╡ 1b6bd414-ddc4-4e6c-a4be-2e1d8fe623ef
if zrec
	SPFLT = signal_roi_allpoints(cmdir, string(sexp), string(srun), 
		                 string.(fpoints), string.(xfn),
			             dkavg, ctx, crad, nmin, csize, scl)
end

# ╔═╡ 6461b6c0-4741-4da6-8db2-5c0a809d3eea
if zrec
	SPXFLT = []
	for (i, sf) in enumerate(SPFLT)
		lbl = string("Point",i)
		sc = 1e+4
		sgnorm = signalnorm(sf, sumnf, filtnm; scale=sc)
		scstr = @sprintf "%.1E" sc
		xlb = string(scstr, " x sgn/nm")
		pfxx = plot(filtnm.center, sgnorm, lw=2, label=lbl, legend=:topleft, title="")
		scatter!(filtnm.center, sgnorm, label="")
		xlabel!("λ (nm)")
		ylabel!(xlb)
		push!(SPXFLT, pfxx)
	end
	plot(size=(1050,1050), SPXFLT..., layout=(4,4), titlefontsize=8)
end

# ╔═╡ 2d6fba42-0e9e-4714-bc79-d47b9bab6c5c
begin	
	hroiz, proiz = histo_signal(roiflt,50)
	#prz1 = plot(proiz)
	hroizx, proizx = histo_signal(roiflt, 50, -5.0*stdfltt, 5.0*stdfltt)
	plot(proiz, proizx, layout=(1,2), titlefontsize=8)
end

# ╔═╡ ed8ab0b7-9fc5-4cb3-8436-79a056a5667f
"""
Returns the histogram of the signal per pixel in the ROI


"""
function signal_histos(xfiles::Vector{String}, xfn::Vector{String}, 
					ecorn::NamedTuple{(:topleft, :botright)},
	                dkavg::Float64, nbins::Int=50)

	function histos_flt(flt::String)
		fltimg = lfi.LaserLab.select_image(xfiles, flt)
		fltcimg = fltimg.img .- dkavg
		fltroi, iroiflt = imgroix(fltcimg, ecorn)
		hfltroi, pfltroi = histo_signal(fltroi, nbins)
		pfltroi
	end

	map(flt->histos_flt(flt), xfn)
	
end

# ╔═╡ aa820b03-361a-4fc1-beb2-9e8b9d1419ee
if zroi
	sfroi = signal_roi(xfiles, xfn, ecorner, dkavg, ctx)
	pfhroi = signal_histos(xfiles, xfn, ecorner, dkavg)
	psfroi = plot(filtnm.center, sfroi, lw=2, label=spointf1, title="Sum (no tails)")
	scatter!(filtnm.center, sfroi, label="")
	xlabel!("λ (nm)")
	ylabel!("counts")
	
	psfroin = plot(filtnm.center, sfroi ./filtnm.width, lw=2, label=spointf1, title="Sum (per nm)")
	scatter!(filtnm.center, sfroi ./filtnm.width, label="")
	xlabel!("λ (nm)")
	ylabel!("counts/nm")

	psfroix = plot(filtnm.center, sfroi /sumnf, lw=2, label=spointf1, title="Sum (fraction of white light)")
	scatter!(filtnm.center, sfroi /sumnf, label="")
	xlabel!("λ (nm)")
	ylabel!("fraction of WL")

	psfroiy = plot(filtnm.center, (sfroi/sumnf) ./filtnm.width, lw=2, label=spointf1, title="Sum (fraction of white light per nm)")
	scatter!(filtnm.center, (sfroi/sumnf) ./filtnm.width, label="")
	xlabel!("λ (nm)")
	ylabel!("fraction of WL/nm")

	plot(size=(750,750), psfroi, psfroin, psfroix, psfroiy,
		 layout=(2,2), titlefontsize=8)
end

# ╔═╡ b842a9b8-1642-480c-a5b0-0b5228832c16
if zroi
	plot(size=(750,750), pfhroi..., layout=(5,2), titlefontsize=10)
end

# ╔═╡ 92900aa3-c295-4def-8700-384ddbea43d9
"""
Returns the indexes of the (xmin, ymin), (xmax,ymax) corners from the edge 
"""
function edge_corners(iedge::Vector{Tuple{Int, Int}})
	indxmin = maximum([ii[1] for ii in iedge])
	lindxmin = [ii for ii in iedge if ii[1] == indxmin ]
	zmindy = minimum([ii[2] for ii in lindxmin])
	indxmax = minimum([ii[1] for ii in iedge])
	lindxmax = [ii for ii in iedge if ii[1] == indxmax ]
	zmaxy = maximum([ii[2] for ii in lindxmax])
	(minvx=(indxmin,zmindy ), maxvx = (indxmax, zmaxy))
end

# ╔═╡ f52ea16e-2aaa-41c4-802d-1481ad1f1fb2
function edge_corners(xc1, yc1)
	iedge = zip(xc1, yc1)
	edge_corners(iedge)
end

# ╔═╡ 376f9591-a2e7-464f-bedf-9e4e4e9ed600
ec2 = edge_corners(tte)

# ╔═╡ 11cc6f74-6e79-4e73-8f75-cd8312298cf5
sy = ec2.minvx[1] - ec2.maxvx[1] +1

# ╔═╡ 99fa3559-3df3-4172-9c5b-2f64a87d9603
sx = -ec2.minvx[2] + ec2.maxvx[2] +1

# ╔═╡ aeb87084-72fe-4e82-acfc-0a92bd534bcd
"""
Returns true if coordinates x,y are inside the box defined by the corners of the edge
"""
function inedge(vxe, i::Int64,j::Int64)
	if i <= vxe.minvx[1] &&  j >= vxe.minvx[2]   
		if i >= vxe.maxvx[1] && j <= vxe.maxvx[2] 
			return true
		end
	end
	return false
end

# ╔═╡ 44107b67-1a43-4b78-a928-71f96b5d6ab8


# ╔═╡ 89b2979a-6e18-40df-8ac4-866e3a0639d6
"""
Given an image (img) and the vertex of the image edge (ecorn), returns the image
in the roi defined by a square box around the edge. The returned matrix is sparse
"""
function imgroi2(img::Matrix{Float64}, ecorn; isize=512, prtlv=0)

	#vxe = lfi.LaserLab.indx_from_edge(img, isize)
	#ecorn = lfi.LaserLab.edge_corners(vxe)
	
	sx = abs(ecorn.maxvx[2] - ecorn.minvx[2]) + 1
	sy = abs(ecorn.minvx[1] - ecorn.maxvx[1]) + 1

	if prtlv == 1
		#println("vxe =", vxe)
		println("ecorn = ", ecorn)
		println("sx = ", sx, " sy = ", sy)
	end 

	roi = zeros(sx, sy)
	
	ix = 0
	jx = 1
	for il in 1:isize
		for jl in 1:isize
			if prtlv == 3
				println("i = ", il, " j = ", jl, " ic = ", " in edge? =", 
				     inedge(ecorn, il,jl))
			end
			if inedge(ecorn, il,jl) 
				if ix < sy
					ix+=1
				else
					ix = 1
					jx+=1
				end
				if prtlv == 2
					println("il = ", il, " jl= ", jl, " ix = ", ix, " jx = ", jx)
				end
				roi[jx,ix] = img[il,jl]
			end
		end
	end
	roi
end

# ╔═╡ 07084d5b-ed96-4f43-863d-c5f67b2bc05f
imgroi2(t2mrx, ec2; isize=xsz, prtlv=0)

# ╔═╡ fda4ec0c-ca4f-48e3-ac98-85644cc2ba67
"""
Returns true if coordinates x,y are contained in a square box of size rbox defined 
from point xybox
"""
function in_box(xybox::Tuple{Int64, Int64}, rbox::Int64, i::Int64,j::Int64)
	if i >= xybox[1] - rbox  &&  i <= xybox[1] + rbox   
		if j >= xybox[2] - rbox  && j <= xybox[2] + rbox 
			return true
		end
	end
	return false
end

# ╔═╡ 50c2178f-eec0-4d94-b621-db383a509789
"""
Given an image (img) and the vertex of the image edge (vxe), returns the image
in the roi defined by a square box around the edge
"""
function imgbox(img::Matrix{Float64}, xybox::Tuple{Int64, Int64}, 
	            rbox::Int64; isize=512)

	roi = zeros(isize, isize)
	
	for il in 1:isize
		for jl in 1:isize
			if in_box(xybox, rbox, il,jl) 
				roi[il,jl] = img[il,jl]
			end
		end
	end
	roi
end

# ╔═╡ dfc058a5-a988-4244-801d-7685dbe47e1b
md"""
## Computation of spectra
"""

# ╔═╡ 136289e4-f0ca-48f4-8ea0-463f7c49dbbf
"""
Gives the sum of the image, correcting with dark current and adding only
pixels above threshold
"""
function sum_dkth(img::Matrix{Float64}, dimg::Matrix{Float64}, thr::Float64)
	sumx = 0.0
	for i in 1:size(img)[1]
		for j in 1:size(img)[2]
			if img[i,j] - dimg[i,j] > thr
    			sumx += (img[i,j] - dimg[i,j])
			end
		end
	end
	sumx
end

# ╔═╡ 0d02c70b-af32-442e-adb8-fd5a666ba574
"""
Gives the sum of the image, adding only
pixels above threshold
"""
function sum_ovth(img::Matrix{Float64}, thr::Float64)
	sumx = 0.0
	for i in 1:size(img)[1]
		for j in 1:size(img)[2]
			if img[i,j]  > thr
    			sumx += img[i,j]
			end
		end
	end
	sumx
end

# ╔═╡ eac4c715-03ec-4b52-92e5-4dfc4de6b7be
"""
Returns the sum of the signal in the ROI and the histograms of the signal 
"""
function roi_sum_and_signal_histos(xfiles::Vector{String}, xfn::Vector{String}, 
								   ecorn::NamedTuple{(:topleft, :botright)},
	                               dkavg::Float64, ctx::Float64, 
	                               nsigmaT::Float64 = 5.0, nsigmaS::Float64 = 2.5)
	SUMTT= Vector{Float64}(undef,0)
	SUMT= Vector{Float64}(undef,0)
	PLTF = []
	for flt in xfn
		fltimg = lfi.LaserLab.select_image(xfiles, string(flt))
		fltcimg = fltimg.img .- dkavg
		fltroi, iroiflt = imgroix(fltcimg, ecorn)
		hfltroi, pfltroi = histo_signal(fltroi,50)

		meanfltt = mean(fltroi)
		stdfltt = std(fltroi)
		#println("Filter ", flt, "tail: mean =", meanfltt, " std =", stdfltt)
		
		meanflt, stdflt = meanstd_interval(fltroi, 
			                               (-nsigmaT*stdfltt,nsigmaT*stdfltt))
		#println("no tail: mean =", meanflt, " std =", stdflt)
		
		sumtt = sum_thr(fltroi, 0.0)
		sumt = sum_interval(fltroi,(ctx, meanflt + nsigmaS*stdflt))
		
		#println("tail: sum =", sumtt, " no tail =", sumt)
		
		push!(PLTF,pfltroi)
		push!(SUMTT, sumtt)
		push!(SUMT, sumt)
	end
	SUMT, SUMTT, PLTF
end

# ╔═╡ 902e4aaf-346c-4698-81a9-6872713fb18e
function sum_allpoints(cmdir::String, sexp::String, srun::String,
	          xpt::Vector{String}, xfn::Vector{String}, 
	          dkavg::Float64, ctx::Float64, 
              crad::Int64, nmin::Int64, csize::Int64, scl::Int64,
	          nsigmaT::Float64 = 5.0, nsigmaS::Float64 = 2.5)

	function select_roi(pt::String)
		flt1, _   = select_files(cmdir, sexp, srun, "Filter1")
		f1img     = lfi.LaserLab.select_image(flt1, pt)
		img_edge  = Float64.(lfi.LaserLab.sujoy(f1img.imgn, four_connectivity=true))
		img_edgeb = Float64.(lfi.LaserLab.binarize(img_edge, Otsu()))
		iedge     = indx_from_edge(img_edgeb)  
		medge, xedge, yedge = edge_to_mtrx(iedge)
		clusters = dbscan(medge, crad, min_neighbors = nmin, min_cluster_size = csize)
		xc1, yc1 = getclusters(iedge, clusters, scl)
		find_edge_corners(xc1, yc1)
	end
	
	PFLT = []
    for (i,pt) in enumerate(xpt)
		ecorner = select_roi(pt)
		#println("ROI defined by topleft, bottomright: =", ecorner)

		xfiles, _  = select_files(cmdir,sexp,srun, pt)

		fsum, _, _ = roi_sum_and_signal_histos(xfiles, xfn, ecorner, dkavg, ctx,
		                                       nsigmaT, nsigmaS)
		push!(PFLT, fsum)
    end
	PFLT
end

# ╔═╡ 78de6bcd-4173-40e3-b500-499568289ba1
function fixdark(nxdrk)
	splt = split(nxdrk[1],"_")
	str = string(splt[3],"_")
	for i in 4:length(splt) -1
		str = string(str, splt[i],"_")
	end
	string(str, splt[end])
end

# ╔═╡ f9ea012f-9aae-4a8b-89d9-0f86802bb14f
function select_image(xfiles::Vector{String}, xfdrk::Vector{String}, 
	                  nxdrk::Vector{String}, setup::Setup) 
	
	
	imgn =string("Filter_",setup.filter,"_rep_", setup.rep, "_ExpoT_10s_Imag_1.csv")
	idrk =string("Filter_",setup.filter, "_", fixdark(nxdrk) )
    #println("imgn = ", imgn)
	#println("idrk = ", idrk)
	ximg = lfi.LaserLab.get_image_name(xfiles, imgn)
	xdrk = lfi.LaserLab.get_image_name(xfdrk, idrk)
    #println("ximg = ", ximg)
	#println("xdrk = ", xdrk)
	lfi.LaserLab.get_image(ximg), lfi.LaserLab.get_image(xdrk)
end

# ╔═╡ 1be45690-2ce1-45fb-bbc6-f9f4f1bf11e4
function select_f1point_image(xfiles::Vector{String}, point::String) 
	function getxfile(files::Vector{String})
		for (i, xf) in enumerate(files)
			if findfirst(point, xf) != nothing
				return files[i]
			end
		end
	end
	lfi.LaserLab.get_image(getxfile(xfiles))
end

# ╔═╡ 58f3c3ac-698b-4381-bdcb-ca8a9cadc8d8
function get_corrected_image(imgm::NamedTuple{(:img, :imgn)}, 
	                         drkm::NamedTuple{(:img, :imgn)}) 
	img = imgm.img .- drkm.img
	imgn = img ./maximum(img)
    (img = img, imgn = imgn, dark=drkm.img)
end

# ╔═╡ e0f20bf3-cf17-4718-85dc-999b71d391d0
"""
Returns the sum of the signal in the ROI or each filter
"""
function spectrum_roi(xfiles::Vector{String}, xfn::Vector{String}, 
	                  darkimg::Matrix{Float64}, ecorn::NamedTuple, nsigma::Float64)
	SUM= []
	for flt in xfn
		fltimg = lfi.LaserLab.select_image(xfiles, string(flt));
		fltcimg = fltimg.img .- darkimg
		fltroi = imgroi2(fltcimg, ecorn)
		push!(SUM, sum_ovth(fltroi, nsigma * std(darkimg)))
	end
	SUM
end

# ╔═╡ 57fc1b69-1213-4ced-967f-4e54c8b6624a


# ╔═╡ a62b6fd5-959c-421a-a160-c420dae4ca99
"""
Returns the signal in a ROI around the maximum
"""
function spectrum_max(setup::Setup,  
                      xfn::Vector{String}, filtnm::NamedTuple, adctopes::Float64;
					  nsigma::Float64=3.0, drkpnt="Dark")
	
	ZMX = Vector{Float64}()
	ZSM = Vector{Float64}()
	ZI = Vector{Float64}()
	ZJ = Vector{Float64}()

	xfiles, _    = select_files(setup.cmdir,setup.series,
			                  setup.measurement, setup.point)
	xfdrk, nxdrk = select_files(setup.cmdir,setup.series,
			                  setup.measurement, drkpnt)
	for fltr in xfn
		setupf = Setup(setup.cmdir, setup.series, setup.measurement,
			           setup.point, fltr, setup.rep)

		imgm, drkm = select_image(xfiles, xfdrk, nxdrk, setupf)
		cimgz = get_corrected_image(imgm, drkm)
		
		imgmz, imgp = lfi.LaserLab.signal_around_maximum(cimgz.img, 
			                                             cimgz.dark; nsigma=nsigma)
		push!(ZMX,imgp.max)
		push!(ZSM,sum(imgmz.img))
		push!(ZI,imgp.imax)
		push!(ZJ,sum(imgp.jmax))

        #println(" for filter ",fltr, " sum = ", sum(imgmz.img))
	end
    
	DataFrame("fltn" => xfn, "cflt" => filtnm.center, "lflt" => filtnm.left, "rflt" => filtnm.right, "wflt" => filtnm.width,
		      "sum"=> ZSM, "sumpes"=> adctopes *(ZSM ./filtnm.width), "max"=> ZMX, "imax" => ZI, "jmax" => ZJ)	
end

# ╔═╡ ac542727-b476-437e-9bc8-8834a0653355
function spectrum_sum(setup::Setup, 
                      xfn::Vector{String}, filtnm::NamedTuple, adctopes::Float64)
    
    ZSM = Vector{Float64}()
	xfiles, _    = select_files(setup.cmdir,setup.series,
			                  setup.measurement, setup.point)
	xfdrk, nxdrk = select_files(setup.cmdir,setup.series,
			                  setup.measurement, drkpnt)
    #println("in spectrum_sum: setup = ", setup)
    for fltr in xfn
        setupf = Setup(setup.cmdir, setup.series, setup.measurement,
			           setup.point, fltr, setup.rep)

		imgm, drkm = select_image(xfiles, xfdrk, nxdrk, setupf)
		cimgz = get_corrected_image(imgm, drkm)
        push!(ZSM,sum(cimgz.img))
    end

    #println("in spectrum_sum (after): setup = ", setup)
    DataFrame("fltn" => xfn, 
		      "cflt" => filtnm.center, 
		      "lflt" => filtnm.left, 
		      "rflt" => filtnm.right, 
		      "wflt" => filtnm.width,
              "sum"=> ZSM, 
		      "sumpes"=> adctopes *(ZSM ./filtnm.width))	
end

# ╔═╡ 18c767aa-1461-4847-ac0d-26ad5a06dd1c
function get_outpath(setup::Setup, ext="*.csv")
	string(setup.series, "_", setup.measurement, "_", 
                   setup.point, ext)

end

# ╔═╡ effb1278-5896-4b19-a6ad-7f12bf5ba9b5
function spectrum_max_allpoints!(setup::Setup, 
                                xpt::Vector{String}, xfn::Vector{String}, 
                                filtnm::NamedTuple, adctopes::Float64; 
								nsigma::Float64=3.0, odir, 
								etype="csv", drkpnt="Dark")
    for pt in xpt
		setupp = Setup(setup.cmdir, setup.series, setup.measurement,
			           pt, setup.filter, setup.rep)
        xpath = joinpath(setupp.cmdir,setupp.series,setupp.measurement, setupp.point)
        dfiles = string("*.",etype)
		
        sdf = spectrum_max(setupp, xfn, filtnm, adctopes;
		                   nsigma=nsigma)
        sdfnm = get_outpath(setupp, ".csv")
        sdff = joinpath(odir, sdfnm)
	    println("Writing point to  =", sdff)
	    CSV.write(sdff, sdf)
    end
end

# ╔═╡ f9608d49-3604-4c8d-913c-6cbf35f7a85f
function read_spectrum(setup::Setup, csvdir::String, ext=".csv")
    
	sdfnm = string(setup.series, "_", setup.measurement, "_", 
                   setup.point, ext)
    
                  
	sdff = joinpath(csvdir, sdfnm)
	println("reading file =", sdff)

    lfi.LaserLab.load_df_from_csv(csvdir, sdfnm, lfi.LaserLab.enG)
	
end

# ╔═╡ 8b09554f-bf5f-4cc8-ab16-9ac34036f111
function spectrum_fromfile_allpoints(setup::Setup, xpt::Vector{String}, csvdir)
    dfdict = Dict()
    for pt in xpt
        setupp = Setup(setup.cmdir, setup.series, setup.measurement,
			           pt, setup.filter, setup.rep)
        df = read_spectrum(setupp, csvdir)
        dfdict[pt] = df
    end

    dfdict
end

# ╔═╡ 9ad17eec-bb5c-4980-9c32-27e01c5b7fcf
function plot_spectrum_for_point(sdfp, pt, fscale="cflt", escale="sumpes")
	plt = plot(sdfp[!, fscale], sdfp[!, escale], lw=2, label="", 
		       xtickfontsize=8,ytickfontsize=8)
	scatter!(sdfp[!, fscale], sdfp[!, escale], label="")
	xlabel!("λ (nm)")
	ylabel!("pes")
	#yticks!([2e+5,4e+5,6e+5])
    xticks!([0,400,600,800])
	
	#xtickfontsize=18,ytickfontsize=18,xlabel="wavelength",xguidefontsize=18,yscale=:log10,ylabel="flux",yguidefontsize=18,legendfontsize=18) here

	plt
end

# ╔═╡ Cell order:
# ╠═f7bb5111-2fc9-49df-972a-0737182da98c
# ╠═981730a6-61fc-484b-ba3c-66920ee7cf83
# ╠═8e7ec382-c738-11ec-3aae-b50d60f15c4f
# ╠═06b8ed45-43bc-464f-89c0-dc0406312b81
# ╠═8833b198-03e4-4679-8949-0c76546cb847
# ╠═6163ba69-1237-4b49-988e-9a73cfef67f6
# ╠═5edc41bc-b912-44bf-9be5-a013f27a75ab
# ╠═90c97a39-b35f-44ba-9646-f0bb9eead338
# ╠═d30e1ceb-2e90-438e-8554-228aa5dc2a59
# ╠═8bda9d8a-c928-431d-a441-ca64bacf7651
# ╠═2be5dcc0-e7c4-412b-990a-d7edb9967186
# ╠═7ce42aec-b319-4de9-b70c-84046d45a600
# ╠═58269465-ba8c-4840-bbfc-0a27897f3e2a
# ╠═0b1c5662-ec4f-486e-9ee6-7fa6ba953e45
# ╠═b07466c0-dfcd-4c10-ae86-45e71a832476
# ╠═2c75e750-854e-459f-91a6-ba135ae263cf
# ╠═c892d4f2-2678-41eb-8724-6d366178f491
# ╠═eab79cba-ca3d-40d8-9961-257e711bb9ae
# ╠═e65eea70-46b3-4852-85d1-5edef9b21b37
# ╠═c9aaf1cc-80c4-475b-8a81-e00918d91b1e
# ╠═b98ad447-6055-46e5-bb4f-9e67f9c3176a
# ╠═88853edb-dc1f-4e7a-a0ba-1868276d1ada
# ╠═57d96432-4318-4291-8255-bfa5d6d3635c
# ╠═0b4d2c08-a677-492c-b5da-982d3d5096fc
# ╠═161b12b6-88a0-4d4d-bfc5-01310534cbdc
# ╠═87519878-fb5d-405c-9222-a71872216ce7
# ╠═e87f48e3-5e5a-44d5-83de-c520e522e33a
# ╠═50ea2ecc-970f-4630-8c7e-acf5e69cc4c9
# ╠═f26bb6e0-45ac-4419-bcb2-46e2cac1f75b
# ╠═308649b5-5c65-40dd-bc66-5b0273648341
# ╠═4218a405-5cb5-464f-9ce1-5d23daeabbef
# ╠═0ee6da42-aa03-4923-bee9-15b92b6583c5
# ╠═347a0f01-fbee-4195-b3a3-55a29285298d
# ╠═e39d2f51-717a-43d0-90b2-25579e625844
# ╠═904d10aa-aea9-46e7-a558-73bcd061a0ec
# ╠═7472f7b8-b50b-4832-b7f7-0134d1c5ed8f
# ╠═9e4e99fe-351d-4a4c-afb2-3a9ed749209c
# ╠═a0b9922e-189f-445e-8508-130e8e561aba
# ╠═0d846a7e-cd53-40a9-8fd5-c5be630790bb
# ╠═dbc8bed9-ba70-443c-9bd3-eaf91336a603
# ╠═3549eafe-70c4-49d5-a121-11b8d2860804
# ╠═2ad43913-db20-4c21-8a9d-cae776aa9b92
# ╠═952ecc28-e2d7-4833-b820-73689a950ea6
# ╠═9957adcf-48e2-4930-849e-b7601f1d8346
# ╠═b270c34c-177b-41ba-8024-56576770b45c
# ╠═077d7e4e-1b94-4b30-a70a-f3b5d3a6fc46
# ╠═c69c8b9d-50bf-46ce-8614-1fee1661e424
# ╠═aa5bdf46-a4db-4ff7-9e51-8d55bc6c203d
# ╠═95100b23-f017-4861-93c1-4adc571e467e
# ╠═e8b6d609-7357-4996-bda5-f1119b4b0263
# ╠═8c911d89-107b-4445-9bda-e7e9b50ff051
# ╠═ecbc5d11-3397-4495-a21f-fa7151dabcd1
# ╠═9b118763-2739-4535-99c3-da6245ba1eae
# ╠═19cd1ee2-f0a4-4620-89c3-e916c4551246
# ╠═1154dba4-fc68-48c4-ab55-c4a5687d0547
# ╠═4615ac40-0164-4297-9aad-66d39d289d15
# ╠═c156b901-e9e2-4576-858f-dc3aa9ae65ee
# ╠═8774dd9a-98d3-432e-9ff6-51190e6d4326
# ╠═455952f8-9cf4-484a-beef-1fc2810e3b89
# ╠═ddf630eb-6be5-422e-9db5-b1a4348adf01
# ╠═8414f6ab-c4f4-4f01-8082-80468ac4e04b
# ╠═607b0b63-f80e-4d57-bd2f-ccdf74a9af3b
# ╠═0ca0299c-6c33-4e35-9642-66e77b1cedba
# ╠═cf1adbfc-e494-44e2-a6d6-7b89fb9a6be6
# ╠═8946f9d0-7108-4af5-90a7-8f338f59009b
# ╠═e3b52355-a912-42ff-ac39-08a79a1ccee5
# ╠═ca7eb636-f817-4a1a-8b25-7e5ff2d3d0e5
# ╠═ca926619-b396-416e-9f72-4c3254f94f80
# ╠═ff856515-001c-4251-ae41-a763171949e5
# ╠═1968c0dd-8cf9-476a-84a8-b92b37ac02ad
# ╠═a1817474-d089-4245-9263-70314969b43e
# ╠═ba8facfa-970e-4265-95e5-58bbe27f273d
# ╠═939c1d23-35d3-41bf-a854-395457e9ad59
# ╠═cecdf185-4a1b-481e-8bc7-a8c4bb7d4990
# ╠═b81df45b-e64e-4a07-982f-368ae03353c2
# ╠═918126c6-b25f-49b5-b5a2-3c15ba1d9cfd
# ╠═a5dc8f3a-420b-4676-93e2-b6d947f26d4c
# ╠═d389f99a-14c2-408f-ad7b-838e00225357
# ╠═f13173e1-088c-4f5c-ae6a-f9aebcd6bc57
# ╠═479f8f86-372c-4b91-9f73-e57a85d3d194
# ╠═406cc319-e7a9-4c68-b732-774b7d1a7e59
# ╠═2d6fba42-0e9e-4714-bc79-d47b9bab6c5c
# ╠═da7e09c2-cf28-414d-bca2-3b9a35275857
# ╠═e9e8e58a-e1db-49f1-8429-420271fb1852
# ╠═54bd1f6c-2b10-47a1-838f-b428fe6b7635
# ╠═aa820b03-361a-4fc1-beb2-9e8b9d1419ee
# ╠═b842a9b8-1642-480c-a5b0-0b5228832c16
# ╠═a48af8f4-4ed2-45cd-b4e8-9b3106c885f3
# ╠═e38d2a94-008e-46db-a35e-176bb5ae297a
# ╠═1794afb6-6ef0-46d6-b182-d54362b9a07d
# ╠═1b6bd414-ddc4-4e6c-a4be-2e1d8fe623ef
# ╠═6461b6c0-4741-4da6-8db2-5c0a809d3eea
# ╠═ed8ab0b7-9fc5-4cb3-8436-79a056a5667f
# ╠═205fd30f-6a4e-473f-9eb5-e826b8c480b1
# ╠═03b663a3-8e01-4720-ad55-b9085ee5115d
# ╠═1658ecd9-8949-4052-9874-a31248c45821
# ╠═f649a197-c7a8-407d-9e11-6fd811c9d375
# ╠═be87f23e-18b0-4927-ba6d-902180f05489
# ╠═902e4aaf-346c-4698-81a9-6872713fb18e
# ╠═b2d6792d-9de4-4d7e-beb3-879459d3709c
# ╠═50e69dc7-0027-4d0c-bfa8-f6778a5754f3
# ╠═82ddd81b-8aea-4711-97d9-a645121786f8
# ╟─155d5066-935b-4643-8aad-f4c635aa7eec
# ╠═24d380e0-18d2-48ae-a2ea-3cb42d5d75e0
# ╠═5a88cb1e-47b2-45cc-965c-2af9a45e72f5
# ╠═7d3bd063-e821-4bd2-b375-5b0989e49270
# ╠═7083fcc2-d2f0-44fa-b85e-0000bb100c0a
# ╠═ed902bce-fb55-4b96-a0ea-cf335e529531
# ╠═436dbab0-5078-46dc-be07-9f04cdf4c46a
# ╠═46f627cd-4166-40c7-8330-d72ac586d3c0
# ╟─d31856a8-4f2b-4f8b-ab9c-20b4cbb643ea
# ╠═61198991-dad0-44e6-9715-a599d4dac0c9
# ╟─7ea8537e-3951-4453-8140-7e2f31f5d900
# ╠═6d2f6df0-c7b3-4f05-9c07-4f9690372c19
# ╟─2ff4a48e-9621-4375-9b05-ab7424ba98fa
# ╠═9a153985-5b8d-4686-99e6-a8038965dddd
# ╟─087751d5-087a-4a88-9dc1-a599fbbfcae3
# ╠═9752523c-7f50-45cb-9341-c0d59e35f772
# ╟─867d2595-632c-477b-89b7-85a0dd8a8941
# ╠═47161d36-4c22-4ca0-a580-24902fc4e1c4
# ╟─dd1a6f48-cba1-4896-91c9-dfa0ee51b765
# ╟─23c3ee67-80e1-48d1-8296-05c814d30c76
# ╠═b1c38ab5-0018-4200-a5b2-f8a7b24bc129
# ╟─6541aa1a-cbcb-47ec-baed-62c58f4f8ae3
# ╠═07e6d6a8-4556-423d-8600-281750f04707
# ╠═376f9591-a2e7-464f-bedf-9e4e4e9ed600
# ╠═11cc6f74-6e79-4e73-8f75-cd8312298cf5
# ╠═99fa3559-3df3-4172-9c5b-2f64a87d9603
# ╠═72f975a9-f9ad-414b-8aee-4f7820fcf3de
# ╠═07084d5b-ed96-4f43-863d-c5f67b2bc05f
# ╠═e82ede75-0b64-4eb0-8130-748cfdf69945
# ╠═c1b2cf36-d72e-4348-911e-8e44f1834ae4
# ╠═56771073-b0b3-47bf-8578-2dd11a59a9b2
# ╠═f2714261-7bb0-47d7-8aac-c16bb5d1f891
# ╠═81a2bc2a-2821-4aa6-86e9-97eedc0bc51c
# ╠═0d5a7021-6072-464e-836e-f05b1e178b80
# ╠═a759ecf7-7373-46bf-ab15-df39cbfa6814
# ╠═bc421f3e-e091-4f83-bc47-ab8572570e1b
# ╠═2498aa42-2c24-47e3-bf5b-647377af0dbc
# ╠═afeb42ca-b462-4320-b364-98a0b4730e33
# ╠═46b1d54e-4ebf-45a2-b137-9690b7a51d38
# ╠═ec95a04a-bda9-429b-92f2-c79f28a322f0
# ╠═54664106-32d8-4ba9-b3e4-0a926a02309c
# ╠═6b6a3651-62e8-4449-9f9d-d8cf3e9ee04a
# ╠═d71e7c72-d7ae-4915-bc36-aed347d09450
# ╠═5f89baa6-d8f4-4d9c-a625-81cf9375f89c
# ╠═85d6c31f-2fc8-486b-98e9-e5b9e134df6a
# ╠═84170688-fbc6-4676-84f1-126ecce4f5f2
# ╠═0ad139e9-85bc-421b-bdd7-e61711d47454
# ╠═65402da4-601a-4766-ba2c-6735f201ef6a
# ╠═e956095a-e468-4456-9659-b1acd6bd507d
# ╠═58bdc137-fc29-4f5e-9ab0-d6f959cfa71c
# ╠═369dcaf2-ce5d-4641-a8cb-274161749c19
# ╠═fa052909-2ed8-4d90-8144-73147682df0e
# ╠═dfda19e3-e772-4458-b29d-2885a020f77b
# ╠═92a5b0b7-5dbb-4e98-a33e-bef1f6992b40
# ╠═aa088449-b0e0-414f-9a96-1a2c5d6656ff
# ╠═8496f1e5-ad5a-4469-851e-0795c8cf6c6b
# ╠═0bf2d20e-dff6-48c2-b059-38b87ed46bb6
# ╠═c3cb6332-0b65-4a49-aa63-a846df3a5277
# ╠═ae5cf21c-7656-4f9f-bd04-73e0e0d8fbee
# ╠═92900aa3-c295-4def-8700-384ddbea43d9
# ╠═f52ea16e-2aaa-41c4-802d-1481ad1f1fb2
# ╠═aeb87084-72fe-4e82-acfc-0a92bd534bcd
# ╠═44107b67-1a43-4b78-a928-71f96b5d6ab8
# ╠═89b2979a-6e18-40df-8ac4-866e3a0639d6
# ╠═fda4ec0c-ca4f-48e3-ac98-85644cc2ba67
# ╠═50c2178f-eec0-4d94-b621-db383a509789
# ╠═dfc058a5-a988-4244-801d-7685dbe47e1b
# ╠═136289e4-f0ca-48f4-8ea0-463f7c49dbbf
# ╠═0d02c70b-af32-442e-adb8-fd5a666ba574
# ╠═eac4c715-03ec-4b52-92e5-4dfc4de6b7be
# ╠═78de6bcd-4173-40e3-b500-499568289ba1
# ╠═f9ea012f-9aae-4a8b-89d9-0f86802bb14f
# ╠═1be45690-2ce1-45fb-bbc6-f9f4f1bf11e4
# ╠═58f3c3ac-698b-4381-bdcb-ca8a9cadc8d8
# ╠═e0f20bf3-cf17-4718-85dc-999b71d391d0
# ╠═57fc1b69-1213-4ced-967f-4e54c8b6624a
# ╠═a62b6fd5-959c-421a-a160-c420dae4ca99
# ╠═ac542727-b476-437e-9bc8-8834a0653355
# ╠═18c767aa-1461-4847-ac0d-26ad5a06dd1c
# ╠═effb1278-5896-4b19-a6ad-7f12bf5ba9b5
# ╠═f9608d49-3604-4c8d-913c-6cbf35f7a85f
# ╠═8b09554f-bf5f-4cc8-ab16-9ac34036f111
# ╠═9ad17eec-bb5c-4980-9c32-27e01c5b7fcf
