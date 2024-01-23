### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ bea979a0-b2c6-11ee-201e-21e7beabe20e
using Pkg; Pkg.activate(ENV["JLaserLab"])

# ╔═╡ 7627d7f2-caa5-465f-b28c-4afd17ab50be
begin
	using PlutoUI
	using CSV
	using DataFrames
	using Plots
	using Printf
	using Markdown
	using InteractiveUtils
	using LsqFit
	using EasyFit
	using Statistics
	using StatsBase
	using LaTeXStrings
end

# ╔═╡ b8ffd5b1-42dd-48fb-97de-86e4af1fc06d
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

# ╔═╡ 752460c2-6bbb-4a25-a188-e797776508ad
lfi = ingredients("../src/LaserLab.jl")


# ╔═╡ 0a7c421d-0c76-4fc3-9f39-53217228670c
adf = DataFrame(CSV.File("AAN155_1E-5M.csv", delim=';', decimal=','))

# ╔═╡ 7321e06b-5c27-4797-9480-4c74f292664f
ntp, _ = lfi.LaserLab.p1df(adf.time, adf.AAN155Free, 50)

# ╔═╡ 53d3a69d-442a-468b-879d-5863c5bbcf59
tpf = plot(ntp[2:50,"x_mean"] ,ntp[2:50,"y_mean"], yerror=ntp[2:50,"y_std"], 
	       fmt = :png,
	       shape = :circle, color = :black, label="ANN155 free", legend=true)

# ╔═╡ 7775493a-cf75-4452-9316-de47b9ac40c0
fitf = fitexp(ntp[2:50,"x_mean"], ntp[2:50,"y_mean"],n=2)

# ╔═╡ ae560dd7-4c79-48cf-9bf7-6398596eafd7


# ╔═╡ 523bd7dc-a449-473c-b938-3e8b59fba050
ntpba, _ = lfi.LaserLab.p1df(adf.time, adf.AAN155_Ba, 50)

# ╔═╡ 0fecd06a-d77a-491d-af24-015b20339f92
tpba = plot(ntpba[2:50,"x_mean"] ,ntpba[2:50,"y_mean"], yerror=ntpba[2:50,"y_std"], 
	       fmt = :png,
	       shape = :circle, color = :black, label="ANN155 Ba2+", legend=true)

# ╔═╡ a171f871-c7c7-4143-8291-fdad8e6a5ad0
fitba = fitexp(ntpba[2:50,"x_mean"], ntpba[2:50,"y_mean"],n=2)

# ╔═╡ cac4c12e-c3d1-4630-9dbb-fbee7de52236
function plot_fit(tdata, xefit, zefit, p3dtfx, p3dtfbax)
	sexp1 = round(zefit.b[1], digits=1)
	sexp2 = round(zefit.b[2], digits=1)
	sexpb1 = round(xefit.b[1], digits=1)
	sexpb2 = round(xefit.b[2], digits=1)
	lnd1 = L"\lambda_1 = %$sexp1 \, \mu s, \lambda_2 = %$sexp2 \, \mu s"
	lnd2 = L"\lambda1 = %$sexpb1 \, \mu s, \lambda2 = %$sexpb2 \, \mu s"
	iexp1 = round(zefit.a[1]/(zefit.a[1]+zefit.a[2]), digits=2) 
	iexp2 = round(zefit.a[2]/(zefit.a[1]+zefit.a[2]), digits=2)
	iexpb1 = round(xefit.a[1]/(xefit.a[1]+xefit.a[2]), digits=2) 
	iexpb2 = round(xefit.a[2]/(xefit.a[1]+xefit.a[2]), digits=2)
	in1 = L"I(\lambda_1) = %$iexp1, I(\lambda_1) = %$iexp2"
	in2 = L"I(\lambda_2) = %$iexpb1, I(\lambda_2) = %$iexpb2"
	lbl1 = string(lnd1,"\n", in1) 
	lbl2 = string(lnd2,"\n", in2)
	
	zzpf3 = plot(p3dtfx, tdata, zefit.ypred,  ylim=(1, 3*10^2), 
		     label=lbl1)

	zzpba3 = plot(p3dtfbax, tdata, xefit.ypred, yaxis=:identity, ylim=(1, 3*10^3),label=lbl2)
	plot(zzpf3, zzpba3)
end

# ╔═╡ a4418194-a4d9-4db1-bef6-0acafde18e47
plot_fit(ntpba[2:50,"x_mean"], fitba, fitf, tpf, tpba)

# ╔═╡ 2b537d8f-fe42-44bd-bfe2-adc07364a1d7
function wmean(x1,x2, w1, w2)
	return x1*w1 + x2 * w2
end

# ╔═╡ daf25718-c6d9-477a-a245-0d18122b869a
wmean(162, 488, 0.88, 0.12)

# ╔═╡ 440fd107-619d-44ae-aeed-8281c0807628
wmean(242, 1171, 0.27, 0.73)

# ╔═╡ ccf82a7a-46d5-4049-bab5-08eb01c5fc04


# ╔═╡ Cell order:
# ╠═bea979a0-b2c6-11ee-201e-21e7beabe20e
# ╠═7627d7f2-caa5-465f-b28c-4afd17ab50be
# ╠═b8ffd5b1-42dd-48fb-97de-86e4af1fc06d
# ╠═752460c2-6bbb-4a25-a188-e797776508ad
# ╠═0a7c421d-0c76-4fc3-9f39-53217228670c
# ╠═7321e06b-5c27-4797-9480-4c74f292664f
# ╠═53d3a69d-442a-468b-879d-5863c5bbcf59
# ╠═7775493a-cf75-4452-9316-de47b9ac40c0
# ╠═ae560dd7-4c79-48cf-9bf7-6398596eafd7
# ╠═523bd7dc-a449-473c-b938-3e8b59fba050
# ╠═0fecd06a-d77a-491d-af24-015b20339f92
# ╠═a171f871-c7c7-4143-8291-fdad8e6a5ad0
# ╠═cac4c12e-c3d1-4630-9dbb-fbee7de52236
# ╠═a4418194-a4d9-4db1-bef6-0acafde18e47
# ╠═2b537d8f-fe42-44bd-bfe2-adc07364a1d7
# ╠═daf25718-c6d9-477a-a245-0d18122b869a
# ╠═440fd107-619d-44ae-aeed-8281c0807628
# ╠═ccf82a7a-46d5-4049-bab5-08eb01c5fc04
