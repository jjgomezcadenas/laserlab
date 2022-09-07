# Model LMU-40X-NUV  : https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=3271

lmu_40x_nuv  = LaserLab.Objective("LMU-40X-NUV", 5.3mm, 5.1mm, 40.0)
lmm40xf_uvv  = LaserLab.Objective("LMM40XF-UVVV", 5.0mm, 5.1mm, 40.0)
rdir = ENV["JLaserLab"]
objdir = joinpath(rdir,"data/objectives")
df  = LaserLab.load_df_from_csv(objdir, "LMU-40X-NUV.csv", LaserLab.enG)
wl  = 200.0:20:800.0
fdf = LaserLab.dftof(wl, df, "Transmission(%)", 0.0, 0.01)

function lmu_40x_nuv_transmission(λ::Unitful.Length)
    lnm = uconvert(nm, λ)
    fdf(lnm/nm)
end

# response of the objective is very flat from 200 nm to 1000 nm, around 85 %
# the approximation of taking the function to be a constant is sufficient for now. 
function lmm40xf_uvv_transmission()
    return 0.85 
end



 