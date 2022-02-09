# Model LMU-40X-NUV  : https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=3271

lmu_40x_nuv = LaserLab.Objective("LMU-40X-NUV", 5.3mm, 5.1mm, 40.0)
df  = LaserLab.load_df_from_csv("/Users/jj/JuliaProjects/LaserLab/data/objectives", "LMU-40X-NUV.csv", LaserLab.enG)
wl  = 200.0:20:800.0
fdf = LaserLab.dftof(wl, df, "Transmission(%)", 0.0, 0.01)

function lmu_40x_nuv_transmission(λ::Unitful.Length)
    lnm = uconvert(nm, λ)
    fdf(lnm/nm)
end



 