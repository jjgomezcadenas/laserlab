using Test
using Dates
using PythonStructs
using LaserLab

function test_tpu()
    
    path="/Users/jjgomezcadenas/Projects/LaserLab/test"
    fname = "test.ptu"
    zpath = joinpath(path,fname)
    nevents = 10
    df = readHH(zpath, nevents, false)
end

function test_pystructs()
    @test LaserLab.hex2dec("FFFF0008") == unpack(">i", hex2bytes("FFFF0008"))[1]
    @test LaserLab.hex2dec("00000008") == unpack(">i", hex2bytes("00000008"))[1]
    @test LaserLab.hex2dec("10000008") == unpack(">i", hex2bytes("10000008"))[1]
    @test LaserLab.hex2dec("11000008") == unpack(">i", hex2bytes("11000008"))[1]
    @test LaserLab.hex2dec("12000008") == unpack(">i", hex2bytes("12000008"))[1]
    @test LaserLab.hex2dec("20000008") == unpack(">i", hex2bytes("20000008"))[1]
    @test LaserLab.hex2dec("21000008") == unpack(">i", hex2bytes("21000008"))[1]
    @test LaserLab.hex2dec("2001FFFF") == unpack(">i", hex2bytes("2001FFFF"))[1]
    @test LaserLab.hex2dec("4001FFFF") == unpack(">i", hex2bytes("4001FFFF"))[1]
    @test LaserLab.hex2dec("4002FFFF") == unpack(">i", hex2bytes("4002FFFF"))[1]
    @test LaserLab.hex2dec("FFFFFFFF") == unpack(">i", hex2bytes("FFFFFFFF"))[1]
end
