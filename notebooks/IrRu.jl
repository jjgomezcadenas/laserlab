### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

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

# ╔═╡ de029a95-a76a-4ccb-bd58-28117d1f6aed
using HypothesisTests

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


# ╔═╡ 4fd4b0e0-396f-471e-a75e-0cc2170673f3
md"# Notebook"

# ╔═╡ 1c69ab8d-0952-4282-89cd-30db5db41392
begin
	xfnm = [420.0,438.0,465.0,503.0,550.0,600.0,650.0,692.0,732.0,810.0]
	wfnm = [10.0, 24.0, 30.0, 40.0, 49.0, 52.0, 60.0, 40.0, 68.0, 10.0]
	filtnm = (center=xfnm,
	width  = wfnm,
	left = xfnm .- 0.5*wfnm,
	right = xfnm .+ 0.5*wfnm,
	tleft=[415.0, 426.0, 450.0, 483.0, 523.0, 574.0, 626.0, 680.0, 712.0, 766.0],
	tright=[425.0, 450.0, 483.0,523.0,574.5,626.0,680.0,712.0,766.0, 815.0])

	#println("Filter central values (nm) = ", filtnm.center)
	#println("Filter width (nm) = ", filtnm.width)
end

# ╔═╡ 43c1ddbb-d475-4738-bc78-af355dcc6ce5
md"""
## Filter data:

- central values  = $(lfi.LaserLab.vect_to_fstr(xfnm, "%7.2f"))
- widths  = $(lfi.LaserLab.vect_to_fstr(wfnm, "%7.2f"))
- left bound  = $(lfi.LaserLab.vect_to_fstr(filtnm.left, "%7.2f"))
- right bound  = $(lfi.LaserLab.vect_to_fstr(filtnm.right, "%7.2f"))
- left bound (th)  = $(lfi.LaserLab.vect_to_fstr(filtnm.tleft, "%7.2f"))
- right bound (th)  = $(lfi.LaserLab.vect_to_fstr(filtnm.tright, "%7.2f"))

"""

# ╔═╡ 5454c96b-c857-4d88-af03-cdcd77055b12
md"""## Mixed surfaces
The plots below show the result of mixing different amounts of RuSl and IrSl, starting from a sample of pure IrSl and finishing in a ratio 1/10 between IrSl and RuSl
"""

# ╔═╡ 7ef4ec0b-7b02-4619-b0a8-58eb6261a44b
md"## IrSl analysis"

# ╔═╡ 22a81e67-02f2-4ca2-90fe-63ede59819c5
md"""
### Formulation

- We take IrSl as a model to develop the mathematical formalism of a Kolmogorov test (KT) able to asess the probability that a given laser measurement is consistent with the "theoretical" IrSL distribution (which is defined as the continuous distribution measured with the fluorimeter)
"""

# ╔═╡ a7a757fc-ef94-4808-83c5-c42c1d332a0d
md"""
## Kolmogorov test
"""

# ╔═╡ d484f83b-f010-4781-8e0a-0dccb2172467
md"""
## Comparing laser measurement with Fluorimeter for IrSl
"""

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

# ╔═╡ f0824bb2-147f-4d07-b862-a9e0b204a58d
md"# Functions"

# ╔═╡ 93f39257-3b07-47b8-aa15-95deaf272b27
function dftof(dir, file, wr, ft=lfi.LaserLab.spG)
	df = lfi.LaserLab.load_df_from_csv(dir, file, ft)
	df[!,"W"] =Float64.(df[!,"W"])
	ff = lfi.LaserLab.dftof(wr, df, "I")
	lfi.LaserLab.ftopdf(wr, ff)
end

# ╔═╡ 29095b9d-8d04-4527-acb0-a61311872e63
begin
	cmdir   = "/Users/jj/JuliaProjects/LaserLab/fluori/"
	irdir   = "/Users/jj/JuliaProjects/LaserLab/fluori/IrSl"
	rudir   = "/Users/jj/JuliaProjects/LaserLab/fluori/RuSl"
	odir    = "/Users/jj/JuliaProjects/LaserLab/pfluori/"
	irfluo  = "IrSlR1.csv"
	rufluo  = "RuSlR1.csv"
	irw=395.0:1.0:850.0
	ruw=395.0:1.0:850.0
	wp = 400.0:4.0:850.0
	dfiles  = "*.csv"
	dplots  = "*.png"
	lmin = 400.0
	lmax = 850.0
	Nir, fir = dftof(irdir, irfluo, irw)
	Nru, fru = dftof(rudir, rufluo, ruw)

	ldIrSl = "/Users/jj/JuliaProjects/LaserLab/data/CMOS/BOLD_068_A1_IrSi_newsetup/csv/"

	lfIrSl = "BOLD_068_A1_IrSi_newsetup_IrSi_20220628_Point1.csv"
	md"""
	## Read and Parameterize fluorimeter data

	- Nir = $(round(Nir, sigdigits=3))
	- Nru = $(round(Nru, sigdigits=3))
	"""
end

# ╔═╡ 81b479a3-4a79-4cbd-a77c-2b11f403d979
begin
	pir = plot(wp, fir.(wp), label="Ir", lw=2)
	pru = plot!(pir, wp, fru.(wp), label="Ru", lw=2)
	xlabel!("λ")
end

# ╔═╡ ec5c6244-d813-4762-8301-e7f10aec1c2d
begin
	pir2 = plot(wp, Nir * fir.(wp), label="Ir", lw=2)
	xlabel!("λ")
	pru2 = plot(wp, Nru* fru.(wp), label="Ru", lw=2)
	xlabel!("λ")
	plot(pir2, pru2)
end

# ╔═╡ 8cfbc972-c959-40c5-8a7f-1066b84b6bbd
begin
	qir = [lfi.LaserLab.qpdf(fir, filtnm.left[l], filtnm.right[l])/filtnm.width[l] for l in 1:length(xfnm)]
	tir = fir.(xfnm)
	
	md"""
	### Discretization and Exact Function
	- Discretization: integrates the signal en the filter inteval and divides by the bin width
	- Exact value: takes the value of the function in the point. 
	- NB: There is a difference between both cases, due to the fact that filters do not overlap perfectly and do not cover the full range.
	"""
end

# ╔═╡ 7d6d1b2e-fd9a-4e99-bbce-4739d318a550
qir

# ╔═╡ 53ee5458-c0c6-4258-b43d-33153a4301dd
tir

# ╔═╡ bc6f4ec4-4481-4b19-a58d-53370e2feb81
ks = ApproximateTwoSampleKSTest(qir, tir)

# ╔═╡ ab7efa9a-b461-4157-90a2-6f79899aef8e
begin
	pir3 = plot(wp, fir.(wp), label="Ir", lw=2)
	sir3 = scatter!(pir3,xfnm, qir, label="discretization")
	sirt = scatter!(sir3, xfnm, tir, label="exact function")
	xlabel!("λ")
end

# ╔═╡ 634f66d1-5ce3-42c6-94b2-2e81b7725e49
begin
	irslcsv = lfi.LaserLab.load_df_from_csv(ldIrSl, lfIrSl, lfi.LaserLab.enG)
	psirsl = scatter(irslcsv.cflt, irslcsv.sumpes)
	plot!(psirsl,irslcsv.cflt, irslcsv.sumpes, lw=2)

end

# ╔═╡ 17ebf4f7-c93d-4b8b-9405-2913a7b9bc00
function fcir(N, R, fir, fru, λmin=400.0, λmax=850.0, λmru=500.0)
	function firu(λ)
		if λ < λmru
			N * fir(λ)
		elseif λ > λmax
			0.0
		else
			N * (fir(λ) + R * fru(λ))
		end
	end
	firu
end

# ╔═╡ 724ac62f-5490-4b8f-a8ea-77553e471ff6
begin
	we = 0.0:2.0:10.0
	firu = [fcir(1.0, r, fir, fru) for r in 0.0:2.0:10.0]
	piru = [plot(wp, firu[i].(wp), label="", lw=2) for i in 1:length(collect(we))]
	plot(piru...)
	#xlabel!("λ")
end

# ╔═╡ 3628a47d-e100-4ccd-9e2d-1f16bdbc2bba
function irfit(df, fir, fru, fitlnm, NN)
	function chi2(p) 
		ff0 = fcir(p, 0.0, fir, fru)
		xqs = [lfi.LaserLab.qpdf(ff0, filtnm.left[l], filtnm.right[l])/filtnm.width[l] for l in 1:length(xfnm)]
		c2 = sum((xqs .-  df.sumpes).^2 ./df.sumpes)
		println("NN = ", p)
		println("xqs = ", xqs)
		println("supmes = ", df.sumpes)
		println("c2 = ", c2)
		c2
	end
	#chi2(NN)
	rr = exp10.(range(log10(NN/100.0), stop=log10(NN/10.), length=5))
	
	[chi2(p) for p in rr]
end

# ╔═╡ 638a0fe1-d1fa-4a11-a96a-84f3b48bb0f6
function CF(f, xmin::Real, xmax::Real)
	function cf(x)
		lfi.LaserLab.qpdf(f, xmin, x) / lfi.LaserLab.qpdf(f, xmin, xmax)
	end
	cf
end

# ╔═╡ b77b128a-6824-4b5f-a6e8-c2564e9891cf
function eCF(qs, filtnm)
	function fm(x)
		#println("x = ", x)
		if x < filtnm.left[1]
			cum = 0.0
		elseif x >= filtnm.right[end]
			cum = 1.0
		else
			cum=0.0
			for n in 1:length(filtnm.center)
				if x > filtnm.tleft[n] && x <= filtnm.tright[n]
					#println("n = ", n)
					for i in 1:n
						cum += qs[i]
					end
					break
				end
			end
			
			cum = cum/sum(qs)
		end
			
		#println("x =", x, " ncum = ", cum/sum(qs))
		cum
	end
	fm
end

# ╔═╡ 277bf000-2814-485a-827c-5114e199e167
begin
	cfir = CF(fir,lmin, lmax)
	ecfirq = eCF(qir, filtnm)
	ecfirt = eCF(tir, filtnm)
	md"""
	### Cumulative function (analytical) and empiric cumulative function
	- cfir is the cumulative function obtained from the "theoretical" value (fluorimeter) of the distribution by integration.
	- ecfir is the empiric cumulative function obtained from a set of points corresponding to a discretization (e.g, from laser measurements)
	"""
end

# ╔═╡ be4f8281-33da-40fc-80fc-6508508df1fd
begin
	pcfir = plot(wp, cfir.(wp), label="cfir", lw=2)
	pcfir2 = plot!(pcfir, wp, ecfirq.(wp), label="ecfir (filters)", lw=2)
	pcfir3 = plot!(pcfir2, wp, ecfirt.(wp), label="ecfir (theoretical)", lw=2)
	xlabel!("λ")
end

# ╔═╡ 6c84f26c-8cba-4804-8b32-46b80b2de675
begin
	ecfirl = eCF(irslcsv.sumpes, filtnm)
	lir = irslcsv.sumpes / sum(irslcsv.sumpes)
end

# ╔═╡ f791682c-1fad-4f8d-bb06-c5e34c2d18f7
begin
	pcfir4 = plot!(pcfir3, wp, ecfirl.(wp), label="ecfir (laser)", lw=2)
	xlabel!("λ")
end

# ╔═╡ 78e1176e-3be6-4fc8-b1f6-87c4c3d155d1
function chi2(f, X::Vector{Float64}, Y::Vector{Float64})
	@assert(length(X) == length(Y))
	sum([(f(X[i]) - Y[i])^2 for i in 1:length(X)])
end

# ╔═╡ a7d36225-9137-45d9-8401-b4af6fb44126
function chi2(X::Vector{Float64}, Y::Vector{Float64})
	@assert(length(X) == length(Y))
	sum([(X[i] - Y[i])^2 for i in 1:length(X)])
end

# ╔═╡ b4528cd1-a87b-4332-a2e3-62427d2b5e4b
md"""
### Computing χ2

χ2 test shows that both the discretization (qir) and theoretical function are obtained from theoretical curve (values close to zero)

- χ2 for qir =  $(chi2(fir, xfnm, qir))
- χ2 for tir =  $(chi2(fir, xfnm, tir))
"""

# ╔═╡ ba404f8c-401a-447b-aac7-4ea35e3da652
chi2(tir, lir)

# ╔═╡ a14210c0-7467-4df1-9dc5-41b1af6138c0
function Dn(F, Fn, X)
	maximum([abs(F(x) - Fn(x)) for x in X])
end

# ╔═╡ a0740faf-1ea3-48bb-a2e3-fa0fa6712a99
Dn(ecfirt, ecfirl, filtnm.center)

# ╔═╡ fb2cef66-2bf5-47f4-a342-e9cb5b6f9473
Dn(ecfirt, ecfirq, filtnm.center)

# ╔═╡ 225380b3-8710-442d-b731-7133c9829e97
function pvalueKS(F, Fn, X)
	d = Dn(F, Fn, X)
	n = length(X)
	α = 2 * exp(-n/d)
	1.0 - α
end

# ╔═╡ 3ac11fab-f7bd-430b-b902-2cede2aecef8
md"""

KS test pvalue (pvalue=1 confirms null hypothesis e.g, two distributions idential at 100 % CL)

- HypothesisTest function: pvalue = $(pvalue(ks))
- Simple implementation : pvalue = $(pvalueKS(ecfirt, ecfirq, filtnm.center))
"""

# ╔═╡ 8479e063-6ab9-4c38-af4b-1eca5c85ee68
function kstest(n::Int64,cl::Float64=0.95)
	sqrt(-log(cl/2.0) *(1.0/n))
end

# ╔═╡ 0a77602d-c409-4965-a4b6-70b9efe87ac9
md"""

KS: Comparing laser with fluorimeter
- Accept hypothesis (at xx % CL) if Dn < KsCL = kstest(n,CL)
- Dn = $(Dn(ecfirt, ecfirl, filtnm.center))
- Ks99 = $(kstest(length(filtnm.center), 0.99))
- pvalue = $(pvalueKS(ecfirt, ecfirl, filtnm.center))
- chi2 = $(chi2(tir, lir))
"""

# ╔═╡ Cell order:
# ╠═182881dc-f9fb-11ec-0f3b-314d43eb0762
# ╠═de40f1a3-5115-4c5c-9c54-2214bf382f85
# ╠═de029a95-a76a-4ccb-bd58-28117d1f6aed
# ╠═12a361f7-6a46-428e-8c1f-a976668b92c8
# ╠═7486f4d0-e8c5-43f9-af89-3e8232ee8b37
# ╠═6c4ca058-4e4d-429f-bae1-509048ad8370
# ╠═4fd4b0e0-396f-471e-a75e-0cc2170673f3
# ╠═43c1ddbb-d475-4738-bc78-af355dcc6ce5
# ╠═1c69ab8d-0952-4282-89cd-30db5db41392
# ╠═29095b9d-8d04-4527-acb0-a61311872e63
# ╠═81b479a3-4a79-4cbd-a77c-2b11f403d979
# ╠═ec5c6244-d813-4762-8301-e7f10aec1c2d
# ╠═5454c96b-c857-4d88-af03-cdcd77055b12
# ╠═724ac62f-5490-4b8f-a8ea-77553e471ff6
# ╠═7ef4ec0b-7b02-4619-b0a8-58eb6261a44b
# ╠═22a81e67-02f2-4ca2-90fe-63ede59819c5
# ╠═8cfbc972-c959-40c5-8a7f-1066b84b6bbd
# ╠═ab7efa9a-b461-4157-90a2-6f79899aef8e
# ╠═b4528cd1-a87b-4332-a2e3-62427d2b5e4b
# ╠═7d6d1b2e-fd9a-4e99-bbce-4739d318a550
# ╠═53ee5458-c0c6-4258-b43d-33153a4301dd
# ╠═277bf000-2814-485a-827c-5114e199e167
# ╠═be4f8281-33da-40fc-80fc-6508508df1fd
# ╠═a7a757fc-ef94-4808-83c5-c42c1d332a0d
# ╠═bc6f4ec4-4481-4b19-a58d-53370e2feb81
# ╠═3ac11fab-f7bd-430b-b902-2cede2aecef8
# ╠═d484f83b-f010-4781-8e0a-0dccb2172467
# ╠═634f66d1-5ce3-42c6-94b2-2e81b7725e49
# ╠═6c84f26c-8cba-4804-8b32-46b80b2de675
# ╠═f791682c-1fad-4f8d-bb06-c5e34c2d18f7
# ╠═a0740faf-1ea3-48bb-a2e3-fa0fa6712a99
# ╠═fb2cef66-2bf5-47f4-a342-e9cb5b6f9473
# ╠═ba404f8c-401a-447b-aac7-4ea35e3da652
# ╠═0a77602d-c409-4965-a4b6-70b9efe87ac9
# ╠═3628a47d-e100-4ccd-9e2d-1f16bdbc2bba
# ╠═37a09cdb-6c63-4937-b5e5-c16b43f918bc
# ╠═f0824bb2-147f-4d07-b862-a9e0b204a58d
# ╠═93f39257-3b07-47b8-aa15-95deaf272b27
# ╠═17ebf4f7-c93d-4b8b-9405-2913a7b9bc00
# ╠═638a0fe1-d1fa-4a11-a96a-84f3b48bb0f6
# ╠═b77b128a-6824-4b5f-a6e8-c2564e9891cf
# ╠═78e1176e-3be6-4fc8-b1f6-87c4c3d155d1
# ╠═a7d36225-9137-45d9-8401-b4af6fb44126
# ╠═225380b3-8710-442d-b731-7133c9829e97
# ╠═a14210c0-7467-4df1-9dc5-41b1af6138c0
# ╠═8479e063-6ab9-4c38-af4b-1eca5c85ee68
