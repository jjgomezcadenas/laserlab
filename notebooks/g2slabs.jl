### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ 7720099d-8e50-43b4-bc99-31c1092a8d08
using Pkg; Pkg.activate("/Users/jj/JuliaProjects/LaserLab/")

# ╔═╡ 4b5ffa44-bd99-45bf-b3d8-632a4840932c
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

# ╔═╡ c6a891aa-36f8-4da6-875e-f337dd486687
import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
    A, N, mol, mmol, V, L, M

# ╔═╡ 48f3499c-7603-4590-84b6-88f89a83bae3
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

# ╔═╡ b2cf92f5-2d8b-4c34-a6da-12872bb2e13f
lfi = ingredients("../src/LaserLab.jl")

# ╔═╡ de5580d2-cc12-4acf-b4f2-cce6a1bba577
PlutoUI.TableOfContents(title="G2 SL abs spectrum", indent=true)

# ╔═╡ bb08587b-2400-4431-886f-156b5282a4c3
md"# Notebook"

# ╔═╡ 1d43b062-f7bb-4e1c-97df-0323cc57ee52
md"## Absorption"

# ╔═╡ fb18d134-df34-4788-ad93-5cbb03251fff
begin
	idir   = "/Users/JJ/JuliaProjects/LaserLab/fluori/G2Sl"
	odir   = "/Users/jj/JuliaProjects/LaserLab/fluori/G2Sl"
	dfiles = "*.csv"
	dplots = "*.png"

	ifile  = "Abs_espectr_BOLD_104_G2-SIL_GB_24h_65DEG_ALL.csv"
	md"""
	### Absorption spectrum G2Sl (in quartz)
	- Input dir = $idir
	- Input file =$ifile
	"""
end

# ╔═╡ 14f4ced7-2757-402e-8c63-d3568e0ac7bd
begin
	g2df = lfi.LaserLab.load_df_from_csv(idir, ifile, lfi.LaserLab.enG; header=[1,2])
	sort!(g2df, rev = false) # sort df so that waveforms go in ascending order
	dfn = names(g2df) # get all the names in the df
	indx = findall([occursin("W", name) for name in dfn]) # indexes of columns with "W"
	for w in dfn[indx]
		g2df[!,w] =Float64.(g2df[!,w]) # convert W to float
	end
	select!(g2df, Not(:Column21_Column21)) # remove missing column 
	first(g2df, 5)
end

# ╔═╡ cf6bc5bd-82ce-4016-96ac-e42b5679bf02
begin
	sdir   = "/Users/JJ/JuliaProjects/LaserLab/fluori/G2Sl/G2_Siltrane_abs_espec_disolucion"
	
	iffrs  = "G2-SL_Libre_disolucin.csv"
	ifbas  = "G2-SL_1eqBa_disolucion.csv"
	md"""
	### Absorption spectrum G2Sl (in solution)
	- Input dir = $sdir
	- Input file (Free) =$iffrs
	- Input file (Ba) =$ifbas
	"""
end

# ╔═╡ 876f74df-e4f7-4bc1-98c4-cafff239de95
md"## Emission"

# ╔═╡ 5756d752-b52c-4f19-9452-28d356f9e5e5
begin
	edirs   = "/Users/JJ/JuliaProjects/LaserLab/fluori/G2Sl/G2_Siltrane_emision_solution"
	
	ifbase  = "G2-SIL_1eqBa_emi_340nm.csv"
	ifres  = "G2-SIL_disolucion_emi_375nm.csv"
	md"""
	### G2Sl on solution
	- Input dir = $edirs
	- Input file (Free) =$ifres
	- Input file (Ba) =$ifbase
	"""
end

# ╔═╡ 218801d5-c39f-4608-86ea-e0c4742852b5
begin
	edir   = "/Users/JJ/JuliaProjects/LaserLab/fluori/G2Sl/BOLD109"
	
	ifba  = "BOLD_109_G2SIL_Ba_on quartz_24h_65C_4h_1E-5M_Ba_exc_375_A7.csv"
	iffr  = "BOLD_109_G2SIL_on_quartz_65C_24h_exc_375A3.csv"
	md"""
	### G2Sl on quartz
	- Input dir = $idir
	- Input file (Free) =$iffr
	- Input file (Ba) =$ifba
	"""
end

# ╔═╡ d0936ab8-65ae-40f9-a1f0-b297fcfe00f8
begin
	g2fr = lfi.LaserLab.load_df_from_csv(edir, iffr, lfi.LaserLab.enG; header=1)
	select!(g2fr, Not(:Column3))
	g2ba = lfi.LaserLab.load_df_from_csv(edir, ifba, lfi.LaserLab.enG; header=1)
	select!(g2ba, Not(:Column3))
	first(g2ba, 5)
end

# ╔═╡ 85207b4d-76c0-482e-bfdf-a2efb2a0527f
plot(g2fr.W, g2fr.I, lw=2, label="G2Sl free")

# ╔═╡ b834feaa-1894-4837-9cc7-698e90be869e
plot(g2ba.W, g2ba.I, lw=2, label="G2Sl Ba2+")

# ╔═╡ 0edf5ba4-baa5-4874-aab7-ab052b58cd6b
begin
	plot(g2fr.W, g2fr.I, lw=2, label="G2Sl free")
	plot!(g2ba.W, g2ba.I, lw=2, label="G2Sl Ba2+")
end

# ╔═╡ 5824f133-7fa0-4131-98c5-713cc7ae8636
md"### Comparison with laser"

# ╔═╡ 00768985-893e-4025-a613-8023932ffa97
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

# ╔═╡ c62a82de-7e3f-436c-bc7d-8407d6bf2aeb
md"""
## Filter data:

- central values  = $(lfi.LaserLab.vect_to_fstr(xfnm, "%7.2f"))
- widths  = $(lfi.LaserLab.vect_to_fstr(wfnm, "%7.2f"))
- left bound  = $(lfi.LaserLab.vect_to_fstr(filtnm.left, "%7.2f"))
- right bound  = $(lfi.LaserLab.vect_to_fstr(filtnm.right, "%7.2f"))
- left bound (th)  = $(lfi.LaserLab.vect_to_fstr(filtnm.tleft, "%7.2f"))
- right bound (th)  = $(lfi.LaserLab.vect_to_fstr(filtnm.tright, "%7.2f"))

"""

# ╔═╡ c80ddf83-584c-4650-a663-898f11f82681
begin
	wr = 396:850
	fg2fr = lfi.LaserLab.dftof(wr, g2fr, "I")
	fg2ba = lfi.LaserLab.dftof(wr, g2ba, "I")
	qsfr = [lfi.LaserLab.qpdf(fg2fr, filtnm.left[l], filtnm.right[l])/filtnm.width[l] for l in 1:length(xfnm)]
	qsfrt = [lfi.LaserLab.qpdf(fg2fr, filtnm.tleft[l], filtnm.tright[l])/filtnm.width[l] for l in 1:length(xfnm)]
	qsba = [lfi.LaserLab.qpdf(fg2ba, filtnm.left[l], filtnm.right[l])/filtnm.width[l] for l in 1:length(xfnm)]
end

# ╔═╡ 4ef47957-b4bc-4b9e-ac24-83828fa7bb38
begin
	plot(g2fr.W, g2fr.I, lw=2, label="G2Sl free")
	scatter!(filtnm.center, qsfr, lw=2, label="G2Sl free filters")
	plot!(filtnm.center, qsfr, lw=2, label="G2Sl free filters")
end

# ╔═╡ 6954fb06-0891-43c0-8ef4-b1633781c6b1
begin
	plot(g2ba.W, g2ba.I, lw=2, label="G2Sl ba2+")
	scatter!(filtnm.center, qsba, lw=2, label="G2Sl ba2+ filters")
	#plot!(filtnm.center, qsba, lw=2, label="G2Sl ba2+ filters")
end

# ╔═╡ 8cf0430e-b6e7-4f2f-b90a-646696213ada
md"## G2SL + Cl2Ba (evaporation)"

# ╔═╡ 82fe96ed-6f82-4480-a313-8e83b03c2c07
begin
	point="Point5"
	g2sia1   = "/Users/JJ/JuliaProjects/LaserLab/data/CMOS/BOLD_104_GBG2_Si_A1/csv"
	g2sia1ba  = "/Users/JJ/JuliaProjects/LaserLab/data/CMOS/BOLD_104_A1_GBG2_Si_Bario/csv"
	a1ba  = string("BOLD_104_A1_GBG2_Si_Bario_GBG2Bario_20220707_",point,".csv")
	a1fr  = string("BOLD_104_GBG2_Si_A1_GBG2_20220630_", point, ".csv")
	md"""
	### G2Sl on quartz + Cl2Ba
	- Input dirs = $g2sia1
	- Input file (Free) =$a1ba
	- Input file (Ba) =$a1ba
	"""
end

# ╔═╡ 08e9b52a-5ea0-4ce5-8ae1-316bfea6f2e3
begin
dfa1 = lfi.LaserLab.load_df_from_csv(g2sia1, a1fr, lfi.LaserLab.enG)
dfa1ba = lfi.LaserLab.load_df_from_csv(g2sia1ba, a1ba, lfi.LaserLab.enG)
end

# ╔═╡ 4082789b-b75a-48fb-93b8-6fbb5cd73f1e
md"# Functions"

# ╔═╡ 948a6f75-e7bf-4839-a075-3e7d8bf5acb8
function loaddf(dir, file, ftype="emi")
	df = lfi.LaserLab.load_df_from_csv(dir, file, lfi.LaserLab.enG)
	select!(df, Not(:Column3))

	if ftype == "abs"
		sort!(df, rev = false) 
		df[!,"W"] =Float64.(df[!,"W"])
	end
	df
end

# ╔═╡ 35acfc65-796f-490d-8aae-0007eac2a4de
begin
	g2frs = loaddf(sdir, iffrs, "abs") 
	g2bas = loaddf(sdir, ifbas, "abs") 
	plot(g2frs.W, g2frs.Abs, lw=2, label="Abs: G2Sl free in solution")
	plot!(g2bas.W, g2bas.Abs, lw=2, label="Abs: G2Sl Ba2+ in solution")
end

# ╔═╡ 0f05f679-8910-4e1c-866c-c6845d10dc5e
begin
	g2frse = loaddf(edirs, ifres, "emi") 
	g2base = loaddf(edirs, ifbase, "emi") 
	p1 = plot(g2frse.W, g2frse.I, lw=2, label="Emission: G2Sl free in solution")
	p2 = plot(g2base.W, g2base.I, lw=2, label="Emission: G2Sl Ba2+ in solution")
	plot(p1, p2)
end

# ╔═╡ d0e6abae-cce3-4c30-a58e-453f5b240978
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

# ╔═╡ ffa9752e-73e3-4d36-869f-e55c5026336a
begin
	lblf = string("Free-", point)
	plot_spectrum_for_point(dfa1, lblf; fscale="cflt", escale="sumpes", plg=:best)
end

# ╔═╡ 98eaeac9-240d-45ab-afd7-edf3c3e1f427
begin
	lblb = string("Ba2+-", point)
	plot_spectrum_for_point(dfa1ba, lblb; fscale="cflt", escale="sumpes", plg=:best)
end

# ╔═╡ 67605946-dd92-4582-ad3a-d42727f86760
function plotg2abs(df, exp="G2"; samples=["A2", "A3"], reps =[0,1,2])

	function gplot(base)
		bs = string(base,samples[1])
		bs = string(bs,"_R",reps[1])
		ww = string(bs,"_W")
		aa = string(bs,"_Abs")
		lbl = string(samples[1], "-R", reps[1])
		p = plot(g2df[!,ww], g2df[!,aa], lw=2, label=lbl)
		
		if length(reps) > 1
			for j in 2:length(reps)
				bs = string(base,samples[1])
				bs = string(bs,"_R",reps[j])
				ww = string(bs,"_W")
				aa = string(bs,"_Abs")
				lbl = string(samples[1], "-R", reps[j])
				#println(lbl)
				p = plot!(p, g2df[!,ww], g2df[!,aa], lw=2, label=lbl)
			end
		end
		
		if length(samples) > 1
			for i in 2:length(samples)
				bs = string(base,samples[i])
				bs = string(bs,"_R",reps[1])
				ww = string(bs,"_W")
				aa = string(bs,"_Abs")
				lbl = string(samples[1], "-R", reps[1])
				p = plot!(g2df[!,ww], g2df[!,aa], lw=2, label=lbl)
				if length(reps) > 1
					for j in 2:length(reps)
						bs = string(base,samples[i])
						bs = string(bs,"_R",reps[j])
						ww = string(bs,"_W")
						aa = string(bs,"_Abs")
						lbl = string(samples[i], "-R", reps[j])
						#println(lbl)
						p = plot!(p, g2df[!,ww], g2df[!,aa], lw=2, label=lbl)
					end
				end
			end
		end
		p
	end
	
	function qplot(base)
		bs = string(base,"_R",reps[1])
		ww = string(bs,"_W")
		aa = string(bs,"_Abs")
		lbl = string(samples[1], "-R", reps[1])
		p = plot(g2df[!,ww], g2df[!,aa], lw=2, label=lbl)

		if length(reps) > 1
			for j in 2:length(reps)
				bs = string(base,"_R",reps[j])
				ww = string(bs,"_W")
				aa = string(bs,"_Abs")
				lbl = string(samples[1], "-R", reps[j])
				#println(lbl)
				p = plot!(p, g2df[!,ww], g2df[!,aa], lw=2, label=lbl)
			end
		end
		p
	end
		
	if exp== "G2"
		base = "BOLD_104_G2SIL_GB_24h_65C_"
		gplot(base)
	else
		base = "Quartz"
		qplot(base)
	end
end

# ╔═╡ 07bbf037-e4dd-4ba2-910e-49a7fc44a6ca
plotg2abs(g2df,"G2"; samples=["A2"], reps=[0,1,2])

# ╔═╡ 6b30f2a1-ffe7-4210-8a4b-74107a87ef38
plotg2abs(g2df; samples=["A3"], reps=[0,1,2])

# ╔═╡ ba44431a-c471-4ebe-8a75-46f1fcef9c7b
plotg2abs(g2df, "G2"; samples=["A2", "A3"], reps=[0,1,2])

# ╔═╡ 9539038a-4c02-44bd-950c-e34e799080a0
plotg2abs(g2df, "Q"; samples=["A2", "A3"], reps=[0,1,2])

# ╔═╡ Cell order:
# ╠═7720099d-8e50-43b4-bc99-31c1092a8d08
# ╠═4b5ffa44-bd99-45bf-b3d8-632a4840932c
# ╠═c6a891aa-36f8-4da6-875e-f337dd486687
# ╠═48f3499c-7603-4590-84b6-88f89a83bae3
# ╠═b2cf92f5-2d8b-4c34-a6da-12872bb2e13f
# ╠═de5580d2-cc12-4acf-b4f2-cce6a1bba577
# ╠═bb08587b-2400-4431-886f-156b5282a4c3
# ╠═1d43b062-f7bb-4e1c-97df-0323cc57ee52
# ╠═fb18d134-df34-4788-ad93-5cbb03251fff
# ╠═14f4ced7-2757-402e-8c63-d3568e0ac7bd
# ╠═07bbf037-e4dd-4ba2-910e-49a7fc44a6ca
# ╠═6b30f2a1-ffe7-4210-8a4b-74107a87ef38
# ╠═ba44431a-c471-4ebe-8a75-46f1fcef9c7b
# ╠═9539038a-4c02-44bd-950c-e34e799080a0
# ╠═cf6bc5bd-82ce-4016-96ac-e42b5679bf02
# ╠═35acfc65-796f-490d-8aae-0007eac2a4de
# ╠═876f74df-e4f7-4bc1-98c4-cafff239de95
# ╠═5756d752-b52c-4f19-9452-28d356f9e5e5
# ╠═0f05f679-8910-4e1c-866c-c6845d10dc5e
# ╠═218801d5-c39f-4608-86ea-e0c4742852b5
# ╠═d0936ab8-65ae-40f9-a1f0-b297fcfe00f8
# ╠═85207b4d-76c0-482e-bfdf-a2efb2a0527f
# ╠═b834feaa-1894-4837-9cc7-698e90be869e
# ╠═0edf5ba4-baa5-4874-aab7-ab052b58cd6b
# ╠═5824f133-7fa0-4131-98c5-713cc7ae8636
# ╠═00768985-893e-4025-a613-8023932ffa97
# ╠═c62a82de-7e3f-436c-bc7d-8407d6bf2aeb
# ╠═c80ddf83-584c-4650-a663-898f11f82681
# ╠═4ef47957-b4bc-4b9e-ac24-83828fa7bb38
# ╠═6954fb06-0891-43c0-8ef4-b1633781c6b1
# ╠═8cf0430e-b6e7-4f2f-b90a-646696213ada
# ╠═82fe96ed-6f82-4480-a313-8e83b03c2c07
# ╠═08e9b52a-5ea0-4ce5-8ae1-316bfea6f2e3
# ╠═ffa9752e-73e3-4d36-869f-e55c5026336a
# ╠═98eaeac9-240d-45ab-afd7-edf3c3e1f427
# ╠═4082789b-b75a-48fb-93b8-6fbb5cd73f1e
# ╠═948a6f75-e7bf-4839-a075-3e7d8bf5acb8
# ╠═d0e6abae-cce3-4c30-a58e-453f5b240978
# ╠═67605946-dd92-4582-ad3a-d42727f86760
