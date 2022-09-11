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

# ╔═╡ c9aaf1cc-80c4-475b-8a81-e00918d91b1e
md"""
# Analysis for a single point
"""

# ╔═╡ 7ce42aec-b319-4de9-b70c-84046d45a600
md"""
## Filters
"""

# ╔═╡ 58269465-ba8c-4840-bbfc-0a27897f3e2a
md"""
### Filter central value in nm
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

# ╔═╡ b98ad447-6055-46e5-bb4f-9e67f9c3176a
md"""
## Conversion factors
"""

# ╔═╡ 88853edb-dc1f-4e7a-a0ba-1868276d1ada
begin
	texp = 10.0
	adctopes = 0.48
	md"""
- Time of exposition (in sec) = $texp
- adc to photoelect = $adctopes
"""
end

# ╔═╡ 88c44435-9685-41c3-a6c8-ffaf7498d60f
md""" Select nsigma cut: $(@bind spsigma NumberField(1.0:10.0, default=5.0))"""

# ╔═╡ f6dcfbd5-8416-42f5-b029-c794f92ee413
nsigma=promsel = Float64(spsigma)

# ╔═╡ dde07319-e6a0-48ac-94d6-92e0fda4df41
md"""
## Cutoff to accept a point (number of sigmas in terms or rms of the dark current)
- nsigma = $nsigma
"""

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

# ╔═╡ f26bb6e0-45ac-4419-bcb2-46e2cac1f75b
begin
	sexp = "AYN_05_Ba"
	srun = "AYN_05_Ba2_45x_G2SL_rep2_P9_220901"
end

# ╔═╡ 308649b5-5c65-40dd-bc66-5b0273648341
md"""
## Algorithm
- Measure the dark current (filter and point independent)
- Select the spot (filter independent)
- For each point, measure the response for all the filters
"""

# ╔═╡ 4218a405-5cb5-464f-9ce1-5d23daeabbef
md"""
### Compute (average) dark current 
- Dark folder contains measurements of the dark current for one or more filter
- All those measurements should be equivalent, since the dark current does not depend of illumination. When more than one measurement exists, they can be averaged
"""

# ╔═╡ b270c34c-177b-41ba-8024-56576770b45c
md"""
### Select spot
- Folder Filter1 contains the measurements carried out with unfiltered light. Those measurements are used to focus the beam and define the spot. 

- The algorithm to find the spot (for each point) is as follows:
	1. Compute the edge of the image (suing a sujoy algorithm, which finds the gradient and binarizes it)
	2. Find the spot delimited by the edge. The simple version implemented here finds a square spot defined between points (xmin, ymin) and (xmax, ymax) in the edge.
	
"""

# ╔═╡ ecbc5d11-3397-4495-a21f-fa7151dabcd1
md"""
#### Find the edge of the image
"""

# ╔═╡ ff46c65c-5b5c-4db5-9265-5f335eca84aa
md"""
#### Find the indexes of the edge
"""

# ╔═╡ 1c2f1eb7-cefa-4df6-8183-edbe99c1e5a7
md"""
#### Find the corners of the edge 
"""

# ╔═╡ ae87b961-58ec-4ca2-b6e9-e4f14b120b87
md"""
#### Find the ROI defined by the corners 
"""

# ╔═╡ 6684ea02-5feb-4328-a174-9ea8e4f2e505
md"""
#### Find the ROI defined by the corners and store it in a sparse matrix
"""

# ╔═╡ cecdf185-4a1b-481e-8bc7-a8c4bb7d4990
md"""
Plot the signal in the ROI
"""

# ╔═╡ ff242019-65d1-449b-bb52-c429966c33e7
md"""
#### Find the ROI defined by a box o radius R around the maximum 
"""

# ╔═╡ a5dc8f3a-420b-4676-93e2-b6d947f26d4c
md"""
### Analysis for Filters 2-11
"""

# ╔═╡ 95e7ab4d-68d5-483b-961b-b20eb9b37d17
md"""
## Select Filter
"""

# ╔═╡ e9e8e58a-e1db-49f1-8429-420271fb1852
md"""
## Comparison between spectrum computed around maximum and using the full camera. 

- Notice the umphysical shoulder when using the full camera
"""

# ╔═╡ 54bd1f6c-2b10-47a1-838f-b428fe6b7635
md""" Check to compute sum using full ROI: $(@bind zroi CheckBox())"""

# ╔═╡ 7e185082-f3bf-4d95-9d5d-57101f47a684
md"""
## Select image

- Image in the left not corrected 
- Image in the right corrected (dark current subtracted)
"""

# ╔═╡ 5477b1dc-ef52-4426-aa47-f513854fcdae
md""" Check to compute sum using full cmos (dark current): $(@bind zfull CheckBox())"""

# ╔═╡ a48af8f4-4ed2-45cd-b4e8-9b3106c885f3
md"""
# Analysis for al points
"""

# ╔═╡ 1794afb6-6ef0-46d6-b182-d54362b9a07d
md""" Check to carry analysis for all points: $(@bind zrec CheckBox())"""

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

# ╔═╡ ae5cf21c-7656-4f9f-bd04-73e0e0d8fbee
"""
Dumps a matrix into a vector
"""
function mtrxtovct(mtrx)
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
function histo_signal(iroi)
	vroi = mtrxtovct(iroi)
	mxvroi = maximum(vroi)
	mnvroi = minimum(vroi)
	lfi.LaserLab.hist1d(vroi, "signal in roi", 100,mnvroi,mxvroi)
end

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

# ╔═╡ 161b12b6-88a0-4d4d-bfc5-01310534cbdc
begin 
	csvdir, pngdir = output_dirs!(sexp)
	md"""
	#### Output dirs
	- csv dir = $csvdir
	- png dir = $pngdir
	"""
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
	#md""" Select experiment : $(@bind sexp Select(dirs))"""
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
	#md""" Select run : $(@bind srun Select(dirs))"""
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
function select_files(cmdir,sexp,srun, spoint, dfiles="*.csv")
	path = joinpath(cmdir,sexp,srun, spoint)
	readdir(path)
	xfiles = Glob.glob(dfiles, path)
	nxfiles = string.([split(f,"/")[end] for f in xfiles])
	xfiles, nxfiles
end


# ╔═╡ def14fbe-f0cf-477f-910a-e8d2ede5eeaf
fd1, nfd1 = select_files(cmdir,sexp,srun, "Dark")

# ╔═╡ 9f898e05-51da-47a9-8654-c592ff5bde01
davgimg, fltnm = lfi.LaserLab.dark_avg(fd1; prnt=false);

# ╔═╡ dfa2081b-dcc2-444c-8e6c-885affb89426
begin
	mxdk, imdk =findmax(davgimg.img)
	heatmap(davgimg.img)
end

# ╔═╡ 842ba4d1-8511-4cef-aca4-f2afffeaa298
md"""
##### Average dark current
- mean = $(mean(davgimg.img))
- std = $(std(davgimg.img))
- max velue $(mxdk) at $imdk
"""

# ╔═╡ dd30581f-7968-4a6e-bdbf-3f69345e372b
nsigma  * std(davgimg.img)

# ╔═╡ 077d7e4e-1b94-4b30-a70a-f3b5d3a6fc46
ff1, nff1 = select_files(cmdir,sexp,srun, "Filter1")

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

# ╔═╡ ba301d3e-3852-4449-bec0-16f53aa7a002
spoint = spointf1

# ╔═╡ e8b6d609-7357-4996-bda5-f1119b4b0263
f1img = lfi.LaserLab.select_image(ff1, string(spointf1));

# ╔═╡ 8c911d89-107b-4445-9bda-e7e9b50ff051
heatmap(f1img.img)

# ╔═╡ 81fb5293-c82f-41ef-aa55-50055f74fcf4
mximg, imximg =findmax(f1img.img)

# ╔═╡ d182ae16-b964-48e2-ab7f-9016e44c5d32
md"""
#### Maximum of the image
- value  = $mximg
- coordinates = $imximg
"""

# ╔═╡ 9b118763-2739-4535-99c3-da6245ba1eae
begin
img_edge = Float64.(lfi.LaserLab.sujoy(f1img.imgn, four_connectivity=true))
img_edgeb = Float64.(lfi.LaserLab.binarize(img_edge, Otsu()))

mosaicview(Gray.(f1img.imgn), Gray.(img_edgeb); nrow = 1)
end

# ╔═╡ 5f8dbfec-6354-49ac-9f72-449c7280ff37
heatmap(img_edge)

# ╔═╡ 7c005ef7-c9f4-4eab-8f98-a002e42138b6
begin
	iedge = lfi.LaserLab.indx_from_edge(img_edgeb)
	xsed = [ii[2] for ii in iedge]
	ysed = [ii[1] for ii in iedge]
end

# ╔═╡ 2864ed9a-143d-4c34-a5ba-d9e5b0f23d27
ecorn = lfi.LaserLab.edge_corners(iedge)

# ╔═╡ 57ca9236-cb82-4429-923e-51ca5feb28bd
begin
	heatmap(img_edge)
	scatter!(xsed, ysed, label="edge",markersize=2)
end

# ╔═╡ e1c80222-4e9f-4dad-9bf9-d799e02b1527
begin
	iroix = lfi.LaserLab.imgroi(f1img.img, ecorn)
	mxxroi, imxroi =findmax(iroix)
end

# ╔═╡ 9ce16f28-f56b-4f74-a7e8-23f909cf2fb6
typeof(imxroi)

# ╔═╡ d2af1b2b-9c86-408c-b008-9e5a8ddd6c0c
iboxx = lfi.LaserLab.imgbox(f1img.img, (imxroi[1],imxroi[2]), 50);

# ╔═╡ 2113ea1a-b784-48fd-bcf1-7b45a85234bd
md"""
Signal in ROIS compared

- sum(roi) from edge --> $(sum(iroix))
- sum(roi) from box  --> $(sum(iboxx))
"""

# ╔═╡ 479f8f86-372c-4b91-9f73-e57a85d3d194
begin
	xfiles, nxfiles  = select_files(cmdir,sexp,srun, spoint)
	xfdrk, nxdrk     = select_files(cmdir,sexp,srun, drkpnt)
end

# ╔═╡ f13173e1-088c-4f5c-ae6a-f9aebcd6bc57
begin
	xfa = lfi.LaserLab.findpattern(nxfiles, "Filter")
	xfb = sort(parse.(Int64, xfa))
	xfn = string.(xfb)
	md""" Select filter: $(@bind sfn Select(xfn))"""
end

# ╔═╡ 8dbf64ec-5854-44b1-ac73-7cd0363a1c6d
begin 
	
	md"""
	### Setup
	- run = $srun
	- experiment = $sexp
	- point = $spoint
	- filter = $sfn
	"""
end

# ╔═╡ d762313e-1466-4724-9c35-dd89e657a11c
begin
	setup = Setup(cmdir, string(sexp), string(srun), string(spoint),        
		              string(sfn), string(rep))
	
	#pngpath = joinpath(pngdir,lfi.LaserLab.get_outpath(setup, "png"))
	md"""
	Setup
	- root dir = $(setup.cmdir)
	- Series = $(setup.series)
	- Measurement = $(setup.measurement)
	- Point = $(setup.point)
	- Filter = $(setup.filter)
	- Repetition = $(setup.rep)
	

	"""
end

# ╔═╡ a90c8edc-1405-4a8c-aacf-53cd130910ae
setup

# ╔═╡ b2623f9c-c0d4-4dd9-8de5-f0bae60c0560
zimg = lfi.LaserLab.select_image(xfiles, string(sfn));

# ╔═╡ 3bfd7691-10eb-4ed3-b3b0-de7e1b0d48c7
heatmap(zimg.img)

# ╔═╡ 406cc319-e7a9-4c68-b732-774b7d1a7e59
zcimg = zimg.img .- davgimg.img;

# ╔═╡ 10aa1194-17c5-40fa-af3a-c8c6db1f4488
begin
	czmx, czimx = findmax(zcimg)
	czmxx, czimxx = findmin(zcimg)
md"""
##### Image for filter $sfn  (dark current subtracted)
- maximum: value = $czmx,  coordinates = $czimx
- minimum: value = $czmxx,  coordinates = $czimxx
- average: $(mean(zcimg))    std = $(std(zcimg))
"""
end

# ╔═╡ e4a4e9d4-9351-4e28-9eb0-224001e2b295
begin
	zmx, zimx = findmax(zimg.img)
	zmxx, zimxx = findmin(zimg.img)
md"""
##### Image for filter $sfn
- maximum: value = $zmx,  coordinates = $zimx
- minimum: value = $zmxx,  coordinates = $zimxx
- average: $(mean(zimg.img))    std = $(std(zimg.img))
"""
end

# ╔═╡ 43ccd7f3-140f-4cb8-af41-e13534c454f3
nxfiles

# ╔═╡ 84170688-fbc6-4676-84f1-126ecce4f5f2
"""
Indexes are colum wise in julia, so (x,y) coordinates corresponde to (second, first) index 
"""
function get_coord_from_indx(xyindx)
	xyindx[2], xyindx[1]
end

# ╔═╡ af83c5aa-ae91-4483-a116-a6fa49993fa4
begin
	exmin, eymin = get_coord_from_indx(ecorn.minvx)
	exmax, eymax = get_coord_from_indx(ecorn.maxvx)
	vxmax, vymax = get_coord_from_indx(imximg)
	
	heatmap(img_edge)
	scatter!([vxmax], [vymax], label="img max",markersize=3)
	scatter!([exmin], [eymin], label="min vertex",markersize=3)
	scatter!([exmax], [eymax], label="max vertex",markersize=3)
end

# ╔═╡ bc4bb613-d4d9-49e3-9c12-335766e9071b
begin
	vxmaxr, vymaxr = get_coord_from_indx(imxroi)
	heatmap(iroix)
	
	scatter!([vxmax], [vymax], label="img max",markersize=3)
	scatter!([vxmaxr], [vymaxr], label="roi max",markersize=3)
	scatter!([exmin], [eymin], label="min vertex",markersize=3)
	scatter!([exmax], [eymax], label="max vertex",markersize=3)
end

# ╔═╡ 1da5474d-dabe-45f1-b856-94085052099d
begin
	heatmap(iboxx)
	scatter!([vxmaxr], [vymaxr], label="roi max",markersize=3)
end

# ╔═╡ e5c09fba-6f8c-4b56-a791-f071532e2dbd
begin
	iroizc = lfi.LaserLab.imgroi(zcimg, ecorn)
	mxxroizc, imxroizc =findmax(iroizc)
	minroizc, imimroizc =findmin(iroizc)
	xmxzc, ymxzc = get_coord_from_indx(imxroizc)
	heatmap(iroizc)
	scatter!([xmxzc], [ymxzc], label="roi max",markersize=3)
end

# ╔═╡ 0ad139e9-85bc-421b-bdd7-e61711d47454
"""
Return a tuple ((x1,x1)...(xn,yn)) with the coordinates of the points in edge
"""
function indx_from_edge(iedge::Matrix{Float64}, isize=512)
	indx = []
	for i in 1:isize
		for j in 1:isize
			if iedge[i,j] == 1
				push!(indx,(i,j))
			end
		end
	end
	indx
end

# ╔═╡ 92900aa3-c295-4def-8700-384ddbea43d9
"""
Returns the indexes of the (xmin, ymin), (xmax,ymax) corners from the edge 
"""
function edge_corners(iedge)
	indxmin = maximum([ii[1] for ii in iedge])
	lindxmin = [ii for ii in iedge if ii[1] == indxmin ]
	zmindy = minimum([ii[2] for ii in lindxmin])
	indxmax = minimum([ii[1] for ii in iedge])
	lindxmax = [ii for ii in iedge if ii[1] == indxmax ]
	zmaxy = maximum([ii[2] for ii in lindxmax])
	(minvx=(indxmin,zmindy ), maxvx = (indxmax, zmaxy))
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

# ╔═╡ d85ce87b-043b-4109-9b3b-2d562bd7aabc


# ╔═╡ 89b2979a-6e18-40df-8ac4-866e3a0639d6
"""
Given an image (img) and the vertex of the image edge (vxe), returns the image
in the roi defined by a square box around the edge
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

# ╔═╡ 52c48abb-dcc8-41cb-bae9-0bbcc1fcbb46
iroix2 = imgroi2(f1img.img, ecorn; prtlv=0);

# ╔═╡ 630b04bd-7060-4f0f-bb5b-0e671ea3fbe7
md"""
Check that the sum and maximum of both ROIS are the same

- sum(roi)  --> $(sum(iroix))
- sum(roi2)  --> $(sum(iroix2))
"""

# ╔═╡ b81df45b-e64e-4a07-982f-368ae03353c2
begin
	hroi, proi = histo_signal(iroix2)
	plot(proi)
end

# ╔═╡ 3d979530-d087-4823-9aa4-fffb78dbd0af
begin
	iroizc2 = imgroi2(zcimg, ecorn)
	mxxroizc2, imxroizc2 =findmax(iroizc2)
	minroizc2, imimroizc2 =findmin(iroizc2)
end

# ╔═╡ 2d6fba42-0e9e-4714-bc79-d47b9bab6c5c
begin
hroiz, proiz = histo_signal(iroizc2)
	plot(proiz)
end

# ╔═╡ 00dc58cd-7575-4410-a430-dde56a890e7a
md"""
##### Image for filter $sfn  (dark current subtracted) in ROI
- maximum: value = $mxxroizc2,  coordinates = $imxroizc
- minimum: value = $minroizc2,  coordinates = $imimroizc
- average: $(mean(iroizc2))    std = $(std(iroizc2))
- Total charge : $(sum(iroizc))
"""

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

# ╔═╡ 4d56de37-7398-49b1-a500-5945f693db1c
sum_ovth(iroix2, nsigma  *std(davgimg.img))

# ╔═╡ 0f736101-628b-4fc4-92fe-cafa7c07e334
sum_ovth(iroizc, nsigma  *std(davgimg.img))

# ╔═╡ 8d4429c8-c7e2-43fd-ba98-24534cca8799
if zroi
	PLTF = []
	SUM= []
	for flt in 2:11
		fltimg = lfi.LaserLab.select_image(xfiles, string(flt));
		fltcimg = fltimg.img .- davgimg.img;
		fltroi = imgroi2(fltcimg, ecorn)
		hfltroi, pfltroi = histo_signal(fltroi)
		push!(PLTF,pfltroi)
		push!(SUM, sum_ovth(fltroi, nsigma  *std(davgimg.img)))
	end
	
	plot(filtnm.center, SUM, lw=2, label=spointf1, title="Signal")
	scatter!(filtnm.center, SUM, label="")
	xlabel!("λ (nm)")
	ylabel!("counts")
end

# ╔═╡ 2761ca35-bd96-4aed-9622-1b33fff2738c
if zroi
	plot(PLTF..., layout=(5,2), titlefontsize=10)
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

# ╔═╡ 4760fdb6-5a0b-4ba2-89b7-0cc7f764d68e
begin
	image, dimage = select_image(xfiles, xfdrk, nxdrk,setup)
	cimg = get_corrected_image(image, dimage)
	mosaicview(Gray.(image.imgn), Gray.(cimg.imgn); nrow = 1)
end

# ╔═╡ 7102dae4-99f8-49a2-9dbd-c69bfeb95546
heatmap(image.img)

# ╔═╡ cdf03376-4450-4b71-a820-33564d1ed71d
begin
	imgmax, imgpos = lfi.LaserLab.signal_around_maximum(cimg.img, cimg.dark; nsigma=nsigma)
	if typeof(imgmax.imgn) == Matrix{Float64}
		Gray.(imgmax.imgn)
	end
end

# ╔═╡ 986b6b9e-3a02-4703-91d0-8088a8066810
let
	spngn = string(setup.series, "_", setup.measurement, "_",
		           setup.point,  "_Filter_", setup.filter, "_imageZoomCMOS.png")
	pxth = joinpath(pngdir, spngn)	
	if typeof(imgmax.imgn) == Matrix{Float64}
		save(pxth, colorview(Gray, map(clamp01nan, imgmax.imgn) ))
	end
end

# ╔═╡ aba1623c-1fd1-4210-b587-a76e21522397
begin
	stot = sum(imgmax.img)
md"""
- Size of the image =$(size(imgmax.img))
- Total sum around the peak = $(round(stot, sigdigits=3))
- value of max = $(imgpos.max)
- position of max: i = $(imgpos.imax)  j = $(imgpos.jmax)
"""
end

# ╔═╡ 8c1ed52f-f5c5-4276-90ac-4d28da153318
function corona(img::Matrix{Float64}, drk::Matrix{Float64}; 
	            nsigma::Float64=2.0, rpxl=15, isize=512)
	
	function crx(imgx::Matrix{Float64}, cutx=0.0)
		local maxx, indx
		maxx, indx = findmax(imgx)
		ixl = max(indx[1] - rpxl, 1)
		ixr = min(indx[1] + rpxl, isize)
		iyl = max(indx[2] - rpxl,1)
		iyr = min(indx[2] + rpxl, isize)
		
		img2 = Array{Float64}(undef,2*rpxl+1,2*rpxl+1)
		for (i, ix) in enumerate(ixl:ixr) 
			for (j, iy) in enumerate(iyl:iyr)
				if img[ix,iy] > cutx
					img2[i,j] = imgx[ix,iy]
				else
					img2[i,j] = 0.0
				end
			end
		end
		maxx, indx, img2
	end
	
	cutoff::Float64 =  nsigma * std(drk)
	cimg::Matrix{Float64} = img .- drk
	
	maxs, indxs, ximgs = crx(img)
	maxd, indxd, ximgd = crx(drk)
	maxc, indxc, ximgc = crx(cimg, 2.0)
	(max=maxs, indx=indxs, ximg=ximgs), (max=maxd, indx=indxd, ximg=ximgd), 
	(max=maxc, indx=indxc, ximg=ximgc)
end

# ╔═╡ a62b6fd5-959c-421a-a160-c420dae4ca99
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

# ╔═╡ d5540947-8b91-4cba-9738-c707e9945eab
if zroi
	sdf = spectrum_max(setup, xfn, filtnm, adctopes; nsigma=nsigma)

	psing = plot(sdf.cflt, sdf.sumpes, lw=2, label=setup.point, title="Spectrum around maximum")
	scatter!(sdf.cflt, sdf.sumpes, label="")
	xlabel!("λ (nm)")
	ylabel!("pes")
end

# ╔═╡ 03c6566e-f4d3-47b6-8a20-4744efc541a0
if zroi
	spsngn = string(setup.series, "_", setup.measurement, "_", setup.point, ".png")
	pxxth = joinpath(pngdir, spsngn)	
	png(psing, pxxth)
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

# ╔═╡ 6325910e-b376-42ab-962a-28f749bc27b2
if zfull
	sdf2 = spectrum_sum(setup,  xfn, filtnm, adctopes)
	plot(sdf2.cflt, sdf2.sumpes, lw=2, label=setup.point, title="spectrum full CMOS")
	scatter!(sdf2.cflt, sdf2.sumpes, label="")
	xlabel!("λ (nm)")
	ylabel!("pes")
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

# ╔═╡ 25219398-6903-4cb0-a336-127eedbfa902
if zrec
	spectrum_max_allpoints!(setup, sspdirs, xfn, 	
		                    filtnm, adctopes; 
                            nsigma=nsigma, odir=csvdir)
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

# ╔═╡ 1f1b334e-941b-4fbd-b964-b4c098ef3231
if zread
	dfdict = spectrum_fromfile_allpoints(setup, sspdirs, csvdir);
	PLT=[]
	for pt in sspdirs
		sdfp = dfdict[pt]
		push!(PLT, plot_spectrum_for_point(sdfp, pt, "cflt"))
	end
	pall = plot(size=(750,750), PLT[1:end]..., layout=(nx,ny), titlefontsize=8)
	
	
end

# ╔═╡ f8b65718-3e1b-454a-817f-1e78feb43225
if zread
	let
		spngn = string(setup.series, "_", setup.measurement, ".png")
		pxth = joinpath(pngdir, spngn)	
		png(pall, pxth)
	end
end

# ╔═╡ Cell order:
# ╠═f7bb5111-2fc9-49df-972a-0737182da98c
# ╠═981730a6-61fc-484b-ba3c-66920ee7cf83
# ╠═8e7ec382-c738-11ec-3aae-b50d60f15c4f
# ╠═06b8ed45-43bc-464f-89c0-dc0406312b81
# ╠═8833b198-03e4-4679-8949-0c76546cb847
# ╠═6163ba69-1237-4b49-988e-9a73cfef67f6
# ╠═5edc41bc-b912-44bf-9be5-a013f27a75ab
# ╠═c9aaf1cc-80c4-475b-8a81-e00918d91b1e
# ╠═7ce42aec-b319-4de9-b70c-84046d45a600
# ╠═58269465-ba8c-4840-bbfc-0a27897f3e2a
# ╠═0b1c5662-ec4f-486e-9ee6-7fa6ba953e45
# ╠═b98ad447-6055-46e5-bb4f-9e67f9c3176a
# ╠═88853edb-dc1f-4e7a-a0ba-1868276d1ada
# ╠═dde07319-e6a0-48ac-94d6-92e0fda4df41
# ╠═f6dcfbd5-8416-42f5-b029-c794f92ee413
# ╠═88c44435-9685-41c3-a6c8-ffaf7498d60f
# ╠═57d96432-4318-4291-8255-bfa5d6d3635c
# ╠═0b4d2c08-a677-492c-b5da-982d3d5096fc
# ╠═161b12b6-88a0-4d4d-bfc5-01310534cbdc
# ╠═87519878-fb5d-405c-9222-a71872216ce7
# ╠═e87f48e3-5e5a-44d5-83de-c520e522e33a
# ╠═50ea2ecc-970f-4630-8c7e-acf5e69cc4c9
# ╠═f26bb6e0-45ac-4419-bcb2-46e2cac1f75b
# ╠═308649b5-5c65-40dd-bc66-5b0273648341
# ╠═4218a405-5cb5-464f-9ce1-5d23daeabbef
# ╠═def14fbe-f0cf-477f-910a-e8d2ede5eeaf
# ╠═9f898e05-51da-47a9-8654-c592ff5bde01
# ╠═dfa2081b-dcc2-444c-8e6c-885affb89426
# ╠═842ba4d1-8511-4cef-aca4-f2afffeaa298
# ╠═b270c34c-177b-41ba-8024-56576770b45c
# ╠═077d7e4e-1b94-4b30-a70a-f3b5d3a6fc46
# ╠═c69c8b9d-50bf-46ce-8614-1fee1661e424
# ╠═aa5bdf46-a4db-4ff7-9e51-8d55bc6c203d
# ╠═95100b23-f017-4861-93c1-4adc571e467e
# ╠═e8b6d609-7357-4996-bda5-f1119b4b0263
# ╠═8c911d89-107b-4445-9bda-e7e9b50ff051
# ╠═d182ae16-b964-48e2-ab7f-9016e44c5d32
# ╠═81fb5293-c82f-41ef-aa55-50055f74fcf4
# ╠═ecbc5d11-3397-4495-a21f-fa7151dabcd1
# ╠═9b118763-2739-4535-99c3-da6245ba1eae
# ╠═5f8dbfec-6354-49ac-9f72-449c7280ff37
# ╠═ff46c65c-5b5c-4db5-9265-5f335eca84aa
# ╠═7c005ef7-c9f4-4eab-8f98-a002e42138b6
# ╠═57ca9236-cb82-4429-923e-51ca5feb28bd
# ╠═1c2f1eb7-cefa-4df6-8183-edbe99c1e5a7
# ╠═2864ed9a-143d-4c34-a5ba-d9e5b0f23d27
# ╠═af83c5aa-ae91-4483-a116-a6fa49993fa4
# ╠═ae87b961-58ec-4ca2-b6e9-e4f14b120b87
# ╠═e1c80222-4e9f-4dad-9bf9-d799e02b1527
# ╠═bc4bb613-d4d9-49e3-9c12-335766e9071b
# ╠═6684ea02-5feb-4328-a174-9ea8e4f2e505
# ╠═52c48abb-dcc8-41cb-bae9-0bbcc1fcbb46
# ╠═630b04bd-7060-4f0f-bb5b-0e671ea3fbe7
# ╠═cecdf185-4a1b-481e-8bc7-a8c4bb7d4990
# ╠═54664106-32d8-4ba9-b3e4-0a926a02309c
# ╠═b81df45b-e64e-4a07-982f-368ae03353c2
# ╠═4d56de37-7398-49b1-a500-5945f693db1c
# ╠═ff242019-65d1-449b-bb52-c429966c33e7
# ╠═9ce16f28-f56b-4f74-a7e8-23f909cf2fb6
# ╠═d2af1b2b-9c86-408c-b008-9e5a8ddd6c0c
# ╠═1da5474d-dabe-45f1-b856-94085052099d
# ╠═2113ea1a-b784-48fd-bcf1-7b45a85234bd
# ╠═a5dc8f3a-420b-4676-93e2-b6d947f26d4c
# ╠═d389f99a-14c2-408f-ad7b-838e00225357
# ╠═ba301d3e-3852-4449-bec0-16f53aa7a002
# ╠═95e7ab4d-68d5-483b-961b-b20eb9b37d17
# ╠═f13173e1-088c-4f5c-ae6a-f9aebcd6bc57
# ╠═479f8f86-372c-4b91-9f73-e57a85d3d194
# ╠═b2623f9c-c0d4-4dd9-8de5-f0bae60c0560
# ╠═3bfd7691-10eb-4ed3-b3b0-de7e1b0d48c7
# ╠═406cc319-e7a9-4c68-b732-774b7d1a7e59
# ╠═e4a4e9d4-9351-4e28-9eb0-224001e2b295
# ╠═10aa1194-17c5-40fa-af3a-c8c6db1f4488
# ╠═e5c09fba-6f8c-4b56-a791-f071532e2dbd
# ╠═3d979530-d087-4823-9aa4-fffb78dbd0af
# ╠═2d6fba42-0e9e-4714-bc79-d47b9bab6c5c
# ╠═dd30581f-7968-4a6e-bdbf-3f69345e372b
# ╠═0f736101-628b-4fc4-92fe-cafa7c07e334
# ╠═00dc58cd-7575-4410-a430-dde56a890e7a
# ╠═e9e8e58a-e1db-49f1-8429-420271fb1852
# ╠═54bd1f6c-2b10-47a1-838f-b428fe6b7635
# ╠═8d4429c8-c7e2-43fd-ba98-24534cca8799
# ╠═2761ca35-bd96-4aed-9622-1b33fff2738c
# ╠═8dbf64ec-5854-44b1-ac73-7cd0363a1c6d
# ╠═d762313e-1466-4724-9c35-dd89e657a11c
# ╠═7e185082-f3bf-4d95-9d5d-57101f47a684
# ╠═43ccd7f3-140f-4cb8-af41-e13534c454f3
# ╠═4760fdb6-5a0b-4ba2-89b7-0cc7f764d68e
# ╠═7102dae4-99f8-49a2-9dbd-c69bfeb95546
# ╠═cdf03376-4450-4b71-a820-33564d1ed71d
# ╠═986b6b9e-3a02-4703-91d0-8088a8066810
# ╠═aba1623c-1fd1-4210-b587-a76e21522397
# ╠═d5540947-8b91-4cba-9738-c707e9945eab
# ╠═03c6566e-f4d3-47b6-8a20-4744efc541a0
# ╠═5477b1dc-ef52-4426-aa47-f513854fcdae
# ╠═6325910e-b376-42ab-962a-28f749bc27b2
# ╠═a90c8edc-1405-4a8c-aacf-53cd130910ae
# ╠═a48af8f4-4ed2-45cd-b4e8-9b3106c885f3
# ╠═1794afb6-6ef0-46d6-b182-d54362b9a07d
# ╠═25219398-6903-4cb0-a336-127eedbfa902
# ╟─82ddd81b-8aea-4711-97d9-a645121786f8
# ╟─155d5066-935b-4643-8aad-f4c635aa7eec
# ╠═1f1b334e-941b-4fbd-b964-b4c098ef3231
# ╠═f8b65718-3e1b-454a-817f-1e78feb43225
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
# ╠═ae5cf21c-7656-4f9f-bd04-73e0e0d8fbee
# ╠═0d5a7021-6072-464e-836e-f05b1e178b80
# ╠═a759ecf7-7373-46bf-ab15-df39cbfa6814
# ╠═bc421f3e-e091-4f83-bc47-ab8572570e1b
# ╠═2498aa42-2c24-47e3-bf5b-647377af0dbc
# ╠═afeb42ca-b462-4320-b364-98a0b4730e33
# ╠═46b1d54e-4ebf-45a2-b137-9690b7a51d38
# ╠═84170688-fbc6-4676-84f1-126ecce4f5f2
# ╠═0ad139e9-85bc-421b-bdd7-e61711d47454
# ╠═92900aa3-c295-4def-8700-384ddbea43d9
# ╠═aeb87084-72fe-4e82-acfc-0a92bd534bcd
# ╠═d85ce87b-043b-4109-9b3b-2d562bd7aabc
# ╠═89b2979a-6e18-40df-8ac4-866e3a0639d6
# ╠═fda4ec0c-ca4f-48e3-ac98-85644cc2ba67
# ╠═50c2178f-eec0-4d94-b621-db383a509789
# ╠═136289e4-f0ca-48f4-8ea0-463f7c49dbbf
# ╠═0d02c70b-af32-442e-adb8-fd5a666ba574
# ╠═78de6bcd-4173-40e3-b500-499568289ba1
# ╠═f9ea012f-9aae-4a8b-89d9-0f86802bb14f
# ╠═1be45690-2ce1-45fb-bbc6-f9f4f1bf11e4
# ╠═58f3c3ac-698b-4381-bdcb-ca8a9cadc8d8
# ╠═8c1ed52f-f5c5-4276-90ac-4d28da153318
# ╠═a62b6fd5-959c-421a-a160-c420dae4ca99
# ╠═ac542727-b476-437e-9bc8-8834a0653355
# ╠═18c767aa-1461-4847-ac0d-26ad5a06dd1c
# ╠═effb1278-5896-4b19-a6ad-7f12bf5ba9b5
# ╠═f9608d49-3604-4c8d-913c-6cbf35f7a85f
# ╠═8b09554f-bf5f-4cc8-ab16-9ac34036f111
# ╠═9ad17eec-bb5c-4980-9c32-27e01c5b7fcf
