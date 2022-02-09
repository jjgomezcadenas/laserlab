@testset "setup_test.jl" begin

    fov = LaserLab.Fov(1.0mm,1.0mm)
    luv = LaserLab.CLaser(375.0nm, 1.0μW)
    obj = LaserLab.Objective("Generic", 0.5, 40.0)
    obj2 = LaserLab.Objective("Generic", 5.0mm, 5.0mm, 40.0)
    
    epl375 = LaserLab.PulsedLaser(375.0nm, 140.0mW, 1.0MHz, 65.0ps)
    lmu_40x_nuv = LaserLab.Objective("LMU-40X-NUV", 5.3mm, 5.1mm, 40.0)
    gl  = LaserLab.GaussianLaser(epl375,lmu_40x_nuv)

    @test LaserLab.lmu_40x_nuv_transmission(0.500μm) ==LaserLab.lmu_40x_nuv_transmission(500nm)
    @test fov.a ≈ π * (fov.d/2.)^2 
    @test obj2.NA ≈ 0.5
    @test fov.v ≈ fov.a * fov.z
    @test gl.w0 ≈ luv.λ/(π * lmu_40x_nuv.NA)
    @test isapprox(LaserLab.photon_energy(500.0nm)/eV, hnu_eV(500.0), rtol=0.01)
    @test uconvert(J, LaserLab.delivered_energy(luv, 1*s)) ≈ 1.0μJ
    @test LaserLab.n_photons(luv) ≈ LaserLab.n_photons(375.0nm, 1.0μW)
    @test LaserLab.n_photons(luv) ≈ LaserLab.n_photons(375.0nm, 1.0μW)
    @test LaserLab.n_photons_int(luv, 1.0s) ≈ LaserLab.n_photons(luv)/Hz
    @test LaserLab.photon_density(375.0nm, 1.0μW, fov.a) ≈ LaserLab.photon_density(luv, fov)
end
