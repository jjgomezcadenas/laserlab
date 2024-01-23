### A Pluto.jl notebook ###
# v0.19.36

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
using Pkg; Pkg.activate(ENV["JLaserLab"])

# ╔═╡ 33dffd80-d690-11ec-0727-2552b82290be
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
	odir    = "/Users/jj/JuliaProjects/LaserLab/pdata/CMOS"
	dfiles  = "*.csv"
	dplots  = "*.png"
	md"""
	CMOS dir = $cmdir
	"""
end

# ╔═╡ e24137d1-213b-4968-87f8-0c83858bbf69
md"""
## Comparing two samples plot
"""

# ╔═╡ d03753c6-8ed0-4c56-aa51-7e974f944274
let
	readdir(cmdir)
	dirs = lfi.LaserLab.getbolddirs(cmdir)
	md""" Select experiment 1 : $(@bind sexp Select(dirs))"""
end

# ╔═╡ 27448883-c43d-4426-ba16-4e85a5b741d8
let
	readdir(cmdir)
	dirs = lfi.LaserLab.getbolddirs(cmdir)
	md""" Select experiment 2 : $(@bind sexp2 Select(dirs))"""
end

# ╔═╡ b54e7056-29a1-4a5e-928a-00a5fb1e6ad9
begin
	md""" Select point : $(@bind pnt NumberField(1:9, default=1)
)"""
end

# ╔═╡ 6bbf2263-858c-4460-b54b-8302215dca51
md"""
## Functions
"""

# ╔═╡ c9d6232b-a457-4545-8a2a-52a60e0e74ee
function create_dir!(dir)
	if isdir(dir) == false
		mkdir(dir)
	end
end

# ╔═╡ 7a6c667b-014d-4736-99aa-b976787515ad
function output_dir2!(odir, sexp, sexp2)
	namex = string(sexp,"__", sexp2)
	od    = joinpath(odir, namex)
	create_dir!(od)
	od
end

# ╔═╡ 12730c01-16d0-4974-b387-f735c01be6be
function find_maxima(sdfp, fscale="cflt", escale="sumpes")
	max, imax = findmax(sdfp[!, escale])
	sdfp[imax, fscale], max
end

# ╔═╡ 0a623b94-459e-472f-b630-6d43116420f9
function select_points(cmdir,sexp, dfiles="*.csv")
	path = joinpath(cmdir,sexp, "csv")
	readdir(path)
	xfiles = Glob.glob(dfiles, path)
	#nxfiles = string.([split(f,"/")[end] for f in xfiles])
	xfiles
end

# ╔═╡ 7c85c42f-029e-4d8c-9f17-609063f9eade
xpoints = select_points(cmdir,sexp)

# ╔═╡ 27a06464-ff05-4f65-8016-fe6909f2d92e
xpoints2 = select_points(cmdir,sexp2)

# ╔═╡ d53d1049-686a-4ec0-bfa6-008a332bb273
begin
	dfpnt = DataFrame(CSV.File(xpoints[pnt]))
	dfpnt2 = DataFrame(CSV.File(xpoints2[pnt]))
	od = output_dir2!(odir, sexp, sexp2)
	md"""
	output plots in $od
	"""
end

# ╔═╡ 246efa59-1b37-4d57-bee1-b9a1d19c6ecd
begin
	max1 = find_maxima(dfpnt)
	max2 = find_maxima(dfpnt2)
md"""

- Maximum of $(sexp) at $(max1[1]) nm with value $(round(max1[2], sigdigits=2)) (pes)
- Maximum of $(sexp2) at $(max2[1]) nm with value $(round(max2[2], sigdigits=2)) (pes)
"""
end

# ╔═╡ ff5589b7-08ea-4911-b395-1e9b3c455c4c
function plot_spectrum_for_point(sdfp, labels; fscale="cflt", escale="sumpes", plg=:best)
	plt = plot(sdfp[!, fscale], sdfp[!, escale], lw=2, label="", legend=plg,
		       xtickfontsize=8,ytickfontsize=8)
	scatter!(sdfp[!, fscale], sdfp[!, escale], label=labels,legend=plg)
	xlabel!("λ (nm)")
	ylabel!("pes")
	#yticks!([2e+5,4e+5,6e+5])
    xticks!([0,400,600,800])
	
	#xtickfontsize=18,ytickfontsize=18,xlabel="wavelength",xguidefontsize=18,yscale=:log10,ylabel="flux",yguidefontsize=18,legendfontsize=18) here

	plt
end

# ╔═╡ 55096351-4df1-4776-8462-ed2c37cdc997
function plt12(dfpnt, sexp, dfpnt2, sexp2, od, max1, max2, plg)
	lbl1 = string(sexp, "-max =", round(max1[1],sigdigits=2))
	lbl2 = string(sexp2, "-max =", round(max2[1],sigdigits=2))
	p1 = plot_spectrum_for_point(dfpnt, lbl1; plg=plg)
	p2 = plot_spectrum_for_point(dfpnt2, lbl2; plg=plg)
	pp = plot(p1,p2, layout=(2,1), titlefontsize=8)
	pngf = joinpath(od, "plt12.png")
	png(pp, pngf)
	pp
end

# ╔═╡ f32448a1-2d5a-40de-aa65-5f5361299552
plt12(dfpnt, "RuSl-2", dfpnt2, "RuSl-1", od, max1, max2, :right)

# ╔═╡ 2c3d4812-069a-40df-bd8f-0b194998e0a5
function plot_spectrum_for_2point(sdfp, sdfp2, labels, fscale="cflt", escale="sumpes")
	plt = plot(sdfp[!, fscale], sdfp[!, escale], lw=2, label=labels[1], 
		       xtickfontsize=8,ytickfontsize=8)
	scatter!(sdfp[!, fscale], sdfp[!, escale], label="")
	plot!(sdfp2[!, fscale], sdfp2[!, escale], lw=2, label=labels[2], 
		       xtickfontsize=8,ytickfontsize=8)
	scatter!(sdfp2[!, fscale], sdfp2[!, escale], label="")
	xlabel!("λ (nm)")
	ylabel!("pes")
	#yticks!([2e+5,4e+5,6e+5])
    xticks!([0,400,600,800])
	
	#xtickfontsize=18,ytickfontsize=18,xlabel="wavelength",xguidefontsize=18,yscale=:log10,ylabel="flux",yguidefontsize=18,legendfontsize=18) here

	plt
end

# ╔═╡ 9bf23097-95d1-4fe9-9497-451782f6ac92
function plt1ovr2(dfpnt, sexp, dfpnt2, sexp2, od)
	pp = plot_spectrum_for_2point(dfpnt, dfpnt2, [sexp, sexp2])
	pngf = joinpath(od, "plt1ovr2.png")
	png(pp, pngf)
	pp
end

# ╔═╡ 6da703a7-20c3-4d7c-a9b6-755d5f7252ef
plt1ovr2(dfpnt, sexp, dfpnt2, sexp2, od)

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
# ╠═27448883-c43d-4426-ba16-4e85a5b741d8
# ╠═27a06464-ff05-4f65-8016-fe6909f2d92e
# ╠═b54e7056-29a1-4a5e-928a-00a5fb1e6ad9
# ╠═d53d1049-686a-4ec0-bfa6-008a332bb273
# ╠═f32448a1-2d5a-40de-aa65-5f5361299552
# ╠═6da703a7-20c3-4d7c-a9b6-755d5f7252ef
# ╠═246efa59-1b37-4d57-bee1-b9a1d19c6ecd
# ╠═6bbf2263-858c-4460-b54b-8302215dca51
# ╠═55096351-4df1-4776-8462-ed2c37cdc997
# ╠═9bf23097-95d1-4fe9-9497-451782f6ac92
# ╠═c9d6232b-a457-4545-8a2a-52a60e0e74ee
# ╠═7a6c667b-014d-4736-99aa-b976787515ad
# ╠═12730c01-16d0-4974-b387-f735c01be6be
# ╠═0a623b94-459e-472f-b630-6d43116420f9
# ╠═ff5589b7-08ea-4911-b395-1e9b3c455c4c
# ╠═2c3d4812-069a-40df-bd8f-0b194998e0a5
# ╠═8ce4c1f7-fd0e-4435-9e78-38940b034890
