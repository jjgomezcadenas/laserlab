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

include("dffunctions_test.jl")
include("setup_test.jl")

