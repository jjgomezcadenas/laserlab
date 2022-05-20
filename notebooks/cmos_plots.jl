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

# ╔═╡ 61849c40-2a3b-4d53-a445-59492ad616e1
using Pkg; Pkg.activate("/Users/jj/JuliaProjects/LaserLab/")

# ╔═╡ 33dffd80-d690-11ec-0727-2552b82290be
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

# ╔═╡ a3e442db-eadc-470c-a0c5-a030cfdef360
import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
    A, N, mol, mmol, V, L, M

# ╔═╡ 2fadfcba-e14e-418a-939e-fc4894dc13b1
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

# ╔═╡ 0674ecbc-2166-4770-8775-52de84611c2e
lfi = ingredients("../src/LaserLab.jl")

# ╔═╡ b11df77d-5d33-4393-b8d7-96c53f431f69
begin
	cmdir   = "/Users/jj/JuliaProjects/LaserLab/data/CMOS"
	odir    = "/Users/jj/JuliaProjects/LaserLab/data/CMOS"
	dfiles  = "*.csv"
	dplots  = "*.png"
	md"""
	CMOS dir = $cmdir
	"""
end

# ╔═╡ e24137d1-213b-4968-87f8-0c83858bbf69
md"""
## Single point plot
"""

# ╔═╡ d03753c6-8ed0-4c56-aa51-7e974f944274
let
	readdir(cmdir)
	dirs = lfi.LaserLab.getbolddirs(cmdir)
	md""" Select experiment : $(@bind sexp Select(dirs))"""
end

# ╔═╡ b54e7056-29a1-4a5e-928a-00a5fb1e6ad9
begin
	md""" Select point : $(@bind pnt NumberField(1:9, default=1)
)"""
end

# ╔═╡ a385e6d5-195a-474e-8a59-777b7c647d65
md"""
## Single point comparison between experiments: G2
"""

# ╔═╡ 6da703a7-20c3-4d7c-a9b6-755d5f7252ef


# ╔═╡ 6bbf2263-858c-4460-b54b-8302215dca51
md"""
## Functions
"""

# ╔═╡ 0a623b94-459e-472f-b630-6d43116420f9
function select_files(cmdir,sexp, dfiles="*.csv")
	path = joinpath(cmdir,sexp,"csv")
	readdir(path)
	xfiles = Glob.glob(dfiles, path)
	#nxfiles = string.([split(f,"/")[end] for f in xfiles])
	xfiles
end

# ╔═╡ 7c85c42f-029e-4d8c-9f17-609063f9eade
xfiles = select_files(cmdir,sexp)

# ╔═╡ d6a0e346-0e93-4f85-9f69-28d514fb2859
xfiles

# ╔═╡ e0b66717-fe68-4e18-87dd-b607ed8312cf
xfiles[pnt]

# ╔═╡ 55096351-4df1-4776-8462-ed2c37cdc997
dfpnt = DataFrame(CSV.File(xfiles[pnt]));

# ╔═╡ cd97fced-e115-4eed-8d04-bad3c6e60649
begin
	g2dirs = ["G2_BOLD_073_B1", "G2_BOLD_078_B2", "G2_BOLD_078_C2"]
	XF = [select_files(cmdir,g2d) for g2d in g2dirs]
end

# ╔═╡ ff5589b7-08ea-4911-b395-1e9b3c455c4c
function plot_spectrum_for_point(sdfp, fscale="cflt", escale="sumpes")
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

# ╔═╡ eb1336af-6868-4a8b-a3b7-2c7738f965e3
plot_spectrum_for_point(dfpnt)

# ╔═╡ 8ce4c1f7-fd0e-4435-9e78-38940b034890
function plot_spectrum_for_points(DF, fscale="cflt", escale="sumpes")
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
# ╠═61849c40-2a3b-4d53-a445-59492ad616e1
# ╠═33dffd80-d690-11ec-0727-2552b82290be
# ╠═a3e442db-eadc-470c-a0c5-a030cfdef360
# ╠═2fadfcba-e14e-418a-939e-fc4894dc13b1
# ╠═0674ecbc-2166-4770-8775-52de84611c2e
# ╠═b11df77d-5d33-4393-b8d7-96c53f431f69
# ╠═e24137d1-213b-4968-87f8-0c83858bbf69
# ╠═d03753c6-8ed0-4c56-aa51-7e974f944274
# ╠═7c85c42f-029e-4d8c-9f17-609063f9eade
# ╠═d6a0e346-0e93-4f85-9f69-28d514fb2859
# ╠═b54e7056-29a1-4a5e-928a-00a5fb1e6ad9
# ╠═e0b66717-fe68-4e18-87dd-b607ed8312cf
# ╠═55096351-4df1-4776-8462-ed2c37cdc997
# ╠═eb1336af-6868-4a8b-a3b7-2c7738f965e3
# ╠═a385e6d5-195a-474e-8a59-777b7c647d65
# ╠═6da703a7-20c3-4d7c-a9b6-755d5f7252ef
# ╠═cd97fced-e115-4eed-8d04-bad3c6e60649
# ╠═6bbf2263-858c-4460-b54b-8302215dca51
# ╠═0a623b94-459e-472f-b630-6d43116420f9
# ╠═ff5589b7-08ea-4911-b395-1e9b3c455c4c
# ╠═8ce4c1f7-fd0e-4435-9e78-38940b034890
