@testset "dffunctions_test.jl" begin
    # test reading out a csv file into a DataFrame, then interpolate it to produce a function 
    spG = LaserLab.CsvG(';',',')  # spanish: decimals represented with ',' delimited with ';'
    enG = LaserLab.CsvG(',','.')  # english
    df = LaserLab.load_df_from_csv("/Users/jj/JuliaProjects/LaserLab/data/objectives", "LMU-40X-NUV.csv", enG)
    wl=200.0:20:800.0
    fdf =LaserLab.dftof(wl, df, "Transmission(%)")
    @test all([fdf(df[1,"Wavelength(µm)"]) == df[1,"Transmission(%)"] for i in 1:31])

    green_area = LaserLab.qpdf(fdf, 400.0, 600.0)
    red_area = LaserLab.qpdf(fdf, 600.0, 800.0)
    @test isapprox(green_area/red_area, 1.88, rtol=0.01)

    total_area = LaserLab.qpdf(fdf, 200.0, 800.0)
    N, dfpdf = LaserLab.ftopdf(wl, fdf)
    @test total_area ≈ N 
end
