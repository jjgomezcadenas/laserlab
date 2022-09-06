
	@testset "pmt_test.jl" begin
        tdir = pwd()
        fword  = joinpath(tdir, "C1--Trace--00010.trc")
        fbyte  = joinpath(tdir, "C1--Trace--00108.trc")
        
        iob = open(fbyte, "r")
        iow = open(fword, "r")
        
        @test typeof(iob) == IOStream
        @test typeof(iow) == IOStream
    
        WAVEDESCb = LaserLab.wavedesc(iob)
        WAVEDESCw = LaserLab.wavedesc(iow)
        @test WAVEDESCb == WAVEDESCw
    
        @test LaserLab.readword(iob, WAVEDESCb + 32) == 0  #Byte
        @test LaserLab.readword(iow, WAVEDESCw + 32) == 1  #Word

        @test LaserLab.readlong(iow, WAVEDESCw + 60) == 2 * LaserLab.readlong(iob, WAVEDESCb + 60)
        
        @test LaserLab.readstring(iow, WAVEDESCw + 76, 14) == "LECROYWS4104HD"
        @test LaserLab.readstring(iob, WAVEDESCb + 76, 14) == "LECROYWS4104HD"

        @test LaserLab.readfloat(iow, WAVEDESCw + 160) == LaserLab.readfloat(iob, WAVEDESCb + 160)

        @test LaserLab.readdouble(iow, WAVEDESCw+ 180) ≈ LaserLab.readdouble(iob, WAVEDESCb + 180)

        tsw = LaserLab.readtimestamp(iow, WAVEDESCw + 296)
        @test tsw.year == 2022
    end

