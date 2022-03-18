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
begin
	using Peaks
	using Glob
	using FFTW
	using DSP
end
import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
    A, N, mol, mmol, V, L, M

import PhysicalConstants.CODATA2018: N_A

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
lfi = ingredients("../src/LaserLab.jl")
PlutoUI.TableOfContents(title="Table of Contents", indent=true)