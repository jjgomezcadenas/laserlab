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

    #rdir = pwd()
    #tdir = joinpath(rdir,"test")
    #fbyte  = joinpath(tdir, "C1--Trace--00108.trc")
    #fword  = joinpath(tdir, "C1--Trace--00010.trc")


    function test_matrix(xs, ys, irng, jrng)
        tmx = zeros(xs, ys)
        indx = []
        for i in irng
            for j in jrng
                tmx[i,j] = 1.0
                push!(indx, (i,j))
            end
        end
        tmx, indx
    end


function test_edge_corners(tedge,edcorn)
    iymax = maximum([ii[2] for ii in tedge])
	ixmax = minimum([ii[1] for ii in tedge])
	iymin = minimum([ii[2] for ii in tedge])
	ixmin = maximum([ii[1] for ii in tedge])
    @test edcorn.minvx == (ixmin, iymin)
    @test edcorn.maxvx == (ixmax, iymax)
    
end


function imgtest()
    xsz=5
	irng = 2:4
	jrng = 2:4
    vx = (3,3)
    rx = 1
    
    tmrx, indx = test_matrix(xsz, xsz, irng, jrng)
    tedge = LaserLab.indx_from_edge(tmrx, xsz)
    edcorn = LaserLab.edge_corners(tedge)
	@test tedge == indx
    test_edge_corners(tedge,edcorn)
    LaserLab.imgroi(tmrx, edcorn; isize=xsz) == tmrx
    LaserLab.imgbox(tmrx, vx, rx; isize=xsz) == tmrx
end
