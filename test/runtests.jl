using LaserLab
using Test
using Unitful
import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W


function hnu_eV(lambda_nm)
    return 1240.0 * (1.0/lambda_nm)
end

dspot(λ::Real, NA::Real) = 1.83*λ/(2*NA)

#include("dffunctions_test.jl")
#include("setup_test.jl")
#include("glaser_test.jl")
#include("pmt_test.jl")

function pmttest()
    rdir = pwd()
    tdir = joinpath(rdir,"test")
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

    xadd = LaserLab.scope_rawdata(iow, WAVEDESCw)
    @test xadd["NOMINAL_BITS"] == 12

    wvf   = LaserLab.xydata(iow, xadd)
    stats = LaserLab.wstats(wvf; nsigma=3.0)

    @test stats.mean + 3.0*stats.std ≈ stats.thrp
    @test stats.mean - 3.0*stats.std ≈ stats.thrn
    #end
end
