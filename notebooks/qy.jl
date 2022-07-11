### A Pluto.jl notebook ###
# v0.19.4

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

# ╔═╡ 182881dc-f9fb-11ec-0f3b-314d43eb0762
using Pkg; Pkg.activate("/Users/jj/JuliaProjects/LaserLab/")


# ╔═╡ de40f1a3-5115-4c5c-9c54-2214bf382f85
begin
	using PlutoUI
	using CSV
	using DataFrames
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
	import Glob
end

# ╔═╡ 12a361f7-6a46-428e-8c1f-a976668b92c8
import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
    A, N, mol, mmol, V, L, M

# ╔═╡ 7486f4d0-e8c5-43f9-af89-3e8232ee8b37
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

# ╔═╡ 6c4ca058-4e4d-429f-bae1-509048ad8370
lfi = ingredients("../src/LaserLab.jl")


# ╔═╡ 29095b9d-8d04-4527-acb0-a61311872e63
begin
	cmdir   = "/Users/jj/JuliaProjects/LaserLab/fluori/"
	
	odir    = "/Users/jj/JuliaProjects/LaserLab/pfluori/"
	dfiles  = "*.csv"
	dplots  = "*.png"
	md"""
	FLUORI dir = $cmdir
	"""
end

# ╔═╡ a45dd231-204d-4842-8c8d-e9c768d893b8
let
	readdir(cmdir)
	dirs = lfi.LaserLab.getbolddirs(cmdir)
	md""" Select directory : $(@bind sdir Select(dirs))"""
end

# ╔═╡ 29544e11-46e6-4868-8457-0f3347c291fa
sdir

# ╔═╡ 0316393b-0a4d-4f7c-9cc4-95c3d7826ca7
begin
	readdir(cmdir)
	xpdir = joinpath(cmdir,sdir)
	xdirs = lfi.LaserLab.getbolddirs(xpdir)
	md""" Select experiment  : $(@bind sexp Select(xdirs))"""
end

# ╔═╡ be09a685-930a-46ab-9b7b-2b83bcc87a7d
sexp

# ╔═╡ 21ec0c9d-7456-4c73-af59-a8491095ba1e
xpdir

# ╔═╡ 49ab6dde-8235-4dd0-b735-f445ece85a2f
md"## Laser data for $sexp"

# ╔═╡ b6e1ba8d-63e6-469e-ab61-afb675ac20e8
md"## Fluorimeter (discretization) and laser data compared"

# ╔═╡ 67c38779-2d01-4dca-a6d1-149c182e81ea
begin
	nmx =split(sexp,".")
	spsngn = string(odir, nmx[1], ".png")
	#pxxth = joinpath(pngdir, spsngn)	
	#png(psing, pxxth)
end

# ╔═╡ f9b8cd59-7cd5-468b-95e4-429a15c1fa51
md" # Notebook"

# ╔═╡ 1c69ab8d-0952-4282-89cd-30db5db41392
begin
	xfnm = [420.0,438.0,465.0,503.0,550.0,600.0,650.0,692.0,732.0,810.0]
	wfnm = [10.0, 24.0, 30.0, 40.0, 49.0, 52.0, 60.0, 40.0, 68.0, 10.0]
	filtnm = (center=xfnm,
	width  = wfnm,
	left = xfnm .- 0.5*wfnm,
	right = xfnm .+ 0.5*wfnm)

	println("Filter central values (nm) = ", filtnm.center)
	println("Filter width (nm) = ", filtnm.width)
end

# ╔═╡ 37a09cdb-6c63-4937-b5e5-c16b43f918bc
if sexp == "QY_450_exc_350.csv"
	csvdir   = "/Users/jj/JuliaProjects/LaserLab/data/CMOS/QY_REFERENCES_WHITE/csv"
	csvfile = "QY_REFERENCES_WHITE_QY_REFERENCES_WHITE_20220627_Point1.csv"
	ww=400.0:0.5:800.0
	wr=400.0:2.0:800.0
elseif sexp == "QY_540_exc_440.csv"
	csvdir   = "/Users/jj/JuliaProjects/LaserLab/data/CMOS/QY_REFERENCES_YELLOW/csv"
	csvfile = "QY_REFERENCES_YELLOW_QY_REFERENCES_YELLOW_20220627_Point1.csv"
	ww=400.0:0.5:900.0
	wr=460.0:2.0:810.0
elseif sexp == "QY_660_exc_510.csv"
	csvdir   = "/Users/jj/JuliaProjects/LaserLab/data/CMOS/QY_REFERENCES_ORANGE/csv"
	csvfile = "QY_REFERENCES_ORANGE_QY_REFERENCES_ORANGE_20220628_Point1.csv"
	ww=400.0:0.5:900.0
	wr=526.0:2.0:810.0
elseif sexp == "IrSlR1.csv"
	csvdir   = "/Users/jj/JuliaProjects/LaserLab/data/CMOS/BOLD_068_A1_IrSi_newsetup/csv/"
	csvfile = "BOLD_068_A1_IrSi_newsetup_IrSi_20220628_Point1.csv"
	ww=400.0:0.5:900.0
	wr=395.0:1.0:850.0
elseif sexp == "RuSlR1.csv"
	csvdir   = "/Users/jj/JuliaProjects/LaserLab/data/CMOS/BOLD_062_RuSi_newsetup/csv/"
	csvfile = "BOLD_062_RuSi_newsetup_RuSi_20220628_Point1.csv"
	ww=400.0:0.5:900.0
	wr=395.0:1.0:850.0
else
	println("WARNING, data not registered, you will get likely an error")
end

# ╔═╡ 634f66d1-5ce3-42c6-94b2-2e81b7725e49
begin
	qycsv = lfi.LaserLab.load_df_from_csv(csvdir, csvfile, lfi.LaserLab.enG)
	psqy = scatter(qycsv.cflt, qycsv.sumpes)
	plot!(psqy,qycsv.cflt, qycsv.sumpes, lw=2)

end

# ╔═╡ 1e5a04f7-cd6b-4cc8-8dce-6a935b9e7720
nqy = [q/maximum(qycsv.sumpes) for q in qycsv.sumpes]

# ╔═╡ 04cb22c9-4961-4853-b1d9-5133ee53044b
begin
	#ww=400.0:0.5:800.0
	ifile = string(sexp)
	qydf = lfi.LaserLab.load_df_from_csv(xpdir, ifile, lfi.LaserLab.spG)
	qydf[!,"W"] =Float64.(qydf[!,"W"])
	fqy = lfi.LaserLab.dftof(wr, qydf, "I")
	md"## Fluorimeter data for QY and discretization with laser filters"
end

# ╔═╡ 440e46e7-9fd8-4f72-8fb3-443680e1af18
qs = [lfi.LaserLab.qpdf(fqy, filtnm.left[l], filtnm.right[l])/filtnm.width[l] for l in 1:length(xfnm)]

# ╔═╡ fdc37ff3-3882-4250-955b-477c9f2404e0
nqs = [q/maximum(qs) for q in qs]

# ╔═╡ e67211a4-5d4c-433d-9426-25afd92c9e8c
begin
	pnqs = scatter(xfnm, nqs)
	pnqy = scatter!(pnqs,xfnm, nqy)
	pnr = plot!(pnqs,xfnm, nqs, lw=2)
	pnr2 = plot!(pnqy,xfnm, nqy, lw=2)
end

# ╔═╡ c0a4d8b5-5419-458a-b7a5-1d4e2a16f31b
begin
	pqyd = plot(qydf.W, qydf.I, label=sexp, lw=1)	
	pfy = plot!(pqyd,collect(ww), fqy.(ww))
	pqs = scatter!(pfy,xfnm, qs)
	plot!(pqs,xfnm, qs, lw=2)
end

# ╔═╡ Cell order:
# ╠═182881dc-f9fb-11ec-0f3b-314d43eb0762
# ╠═de40f1a3-5115-4c5c-9c54-2214bf382f85
# ╠═12a361f7-6a46-428e-8c1f-a976668b92c8
# ╠═7486f4d0-e8c5-43f9-af89-3e8232ee8b37
# ╠═6c4ca058-4e4d-429f-bae1-509048ad8370
# ╠═29095b9d-8d04-4527-acb0-a61311872e63
# ╠═a45dd231-204d-4842-8c8d-e9c768d893b8
# ╠═29544e11-46e6-4868-8457-0f3347c291fa
# ╠═0316393b-0a4d-4f7c-9cc4-95c3d7826ca7
# ╠═be09a685-930a-46ab-9b7b-2b83bcc87a7d
# ╠═21ec0c9d-7456-4c73-af59-a8491095ba1e
# ╠═49ab6dde-8235-4dd0-b735-f445ece85a2f
# ╠═634f66d1-5ce3-42c6-94b2-2e81b7725e49
# ╠═04cb22c9-4961-4853-b1d9-5133ee53044b
# ╠═c0a4d8b5-5419-458a-b7a5-1d4e2a16f31b
# ╠═b6e1ba8d-63e6-469e-ab61-afb675ac20e8
# ╠═e67211a4-5d4c-433d-9426-25afd92c9e8c
# ╠═67c38779-2d01-4dca-a6d1-149c182e81ea
# ╠═f9b8cd59-7cd5-468b-95e4-429a15c1fa51
# ╠═1c69ab8d-0952-4282-89cd-30db5db41392
# ╠═440e46e7-9fd8-4f72-8fb3-443680e1af18
# ╠═fdc37ff3-3882-4250-955b-477c9f2404e0
# ╠═1e5a04f7-cd6b-4cc8-8dce-6a935b9e7720
# ╠═37a09cdb-6c63-4937-b5e5-c16b43f918bc
