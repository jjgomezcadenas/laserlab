### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 1459f36b-afa4-405e-a9e8-e6ae32ea752e
import Pkg

# ╔═╡ 0735b900-2548-412f-a675-c3c19945ae47
Pkg.add("PyCall")

# ╔═╡ ec3057cd-41a8-4c87-8808-5f612bfa82a5
Pkg.add("Dates")

# ╔═╡ c83e18f6-8d4b-4fd7-ad3c-950c4eb4b91d
using PyCall

# ╔═╡ e6550506-f9a3-4999-9fef-4e0cf70cc1ec
using Dates

# ╔═╡ f9470c2c-4386-41b5-9473-f88fb713a0bb

timestamp = time()  # Get current Unix timestamp


# ╔═╡ 73f68dd8-863d-49fe-9166-41fbd0ec36cb
unix2datetime(timestamp)

# ╔═╡ 6753e916-8645-4256-ab5a-191e6cbca9df
gmtime = DateTime(floor(timestamp))

# ╔═╡ b6517e62-c865-46aa-a43f-fa3075760309
pystruct =  pyimport("struct")

# ╔═╡ 8f01ac68-1cbf-4341-aabc-241bc450f869
py"bytes.fromhex"("FFFF0008")

# ╔═╡ e9f05987-8092-4164-a6a9-e0a041be123c
typeof("FFFF0008")

# ╔═╡ 4ca63a99-8679-4f3c-be4a-049442499e77
function bytes_from_hex(x)
	py"""
	import struct
	def bytesfromhex(x):
		print(type(x))
		#rx = bytes.fromhex(x)
		return struct.unpack(">i", bytes.fromhex(x))[0]
	"""
	rx = py"bytesfromhex"(x)
end

# ╔═╡ 511ee273-3401-4a4e-af28-d19b24c26642
bytes_from_hex("FFFF0008")

# ╔═╡ 9f32a075-eab9-421f-801f-eabaa9400251
x = UInt32(3345684)

# ╔═╡ e6a6aea1-93ab-4b75-98ef-c073e7310b83
y = 0x3FF

# ╔═╡ e71daa46-42dd-4b5b-ac40-27e66ffd4d6d
z = Int32(1023)

# ╔═╡ b4941210-d353-441d-8d51-d537c8fcb16e
x&y

# ╔═╡ 62e60c02-321f-47bb-94a6-c59e61824ab2
bitstring(x)

# ╔═╡ 70a97b9a-dab1-4a04-a2ac-2123cf7da83a
bitstring(y)

# ╔═╡ f7ae7992-4a66-4f55-a5db-d8dcdcaed7a3
bitstring(z)

# ╔═╡ b026d200-1763-4ea0-8336-24fc9a06dd4f
bitstring(x&y)

# ╔═╡ 28936553-406e-433f-8a4b-167ca0952af8
bitstring(x&z)

# ╔═╡ 8e7c68c7-ff2e-4895-aaf7-189f1d9eb39d
md"""

## test_tpu
"""

# ╔═╡ 0a064d7b-64da-4038-80e1-d0f8227a9052
function gotOverflow(count, recNum)
    #global outputfile, recNum
    #outputfile.write("%u OFL * %2x\n" % (recNum, count))
	println("gotOverflow: recNum  = ", recNum, " count =", count)
end

# ╔═╡ 1feb11fc-2f54-4c75-a6a0-61b9435bd411
function gotMarker(timeTag, markers, recNum)
	println("gotMarker: recNum  = ", recNum, " timeTag =", timeTag, " markers=",markers)
    #global outputfile, recNum
    #outputfile.write("%u MAR %2x %u\n" % (recNum, markers, timeTag))
end

# ╔═╡ 7f0b9b19-e92b-4f15-8524-5ddee7d18607
function gotPhoton(timeTag, channel, dtime, recNum, isT2, globRes)
    #global outputfile, isT2, recNum
	println("\n gotPhoton: isT2 = ", isT2)
    if isT2
        #outputfile.write(
        #    "%u CHN %1x %u %8.0lf\n"
        #    % (recNum, channel, timeTag, (timeTag * globRes * 1e12))
        #)
		println("recNum, = ", recNum, " channel = ", channel, " timeTag = ", timeTag)
		println("resol, = ", (timeTag * globRes * 1e12))
    else
        #outputfile.write(
        #    "%u CHN %1x %u %8.0lf %10u\n"
        #    % (recNum, channel, timeTag, (timeTag * globRes * 1e9), dtime)
        #)
		println("recNum, = ", recNum, " channel = ", channel, " timeTag = ", timeTag)
		println("resol, = ", (timeTag * globRes * 1e9), " dtime =", dtime)
	end
end

# ╔═╡ e58d12e4-b51e-4bf4-a96b-2e538dc03740
function readHT3(version, fid, numRecords, oflcorrection, isT2, globRes)
   
    T3WRAPAROUND = 1024
	OverflowCorrection = 0
	
	for recNum in 1:10
		t3record   = read(fid, UInt32)
		recordData = bitstring(t3record)
		special = parse(Int32, recordData[1], base=2)  # last bit
		channel = parse(Int32, recordData[2:7], base=2) # next 6 bits
		dtime   = parse(Int32, recordData[8:22], base=2) # next 15 bits
		nsync   = parse(Int32, recordData[23:32], base=2) #lowest 10 bits

		#println("\n recNum = ", recNum)
		#println(" recordData = ", recordData)
		#println(" special = ", special)
		#println(" channel = ", channel)
		#println(" dtime = ", dtime)
		#println(" nsync = ", nsync)

		if special == 0   # this means a regular input channel
        	true_nSync = OverflowCorrection + nsync;
           	# one nsync time unit equals to "syncperiod" which can be
           	# calculated from "SyncRate"

			#println("true_nSync =", true_nSync)
           	gotPhoton(true_nSync, channel, dtime, recNum, isT2, globRes)
			#println("Will call GotPhoton")
			
		else  # this means we have a special record
        	if channel == 63  # overflow of nsync occured
            	if (nsync == 0) || (version == 1) # if nsync is zero it is 
				                                  # an old style      
												  # single oferflow or old Version
                	OverflowCorrection = OverflowCorrection + T3WRAPAROUND
					#println("OverflowCorrection =", OverflowCorrection)
					#println("Will call GotOverflow(1)")
                	gotOverflow(1, recNum)
				else                      # otherwise nsync indicates the number of 
				                          # overflows 
				                          # THIS   IS NEW IN FORMAT V2.0
					
                	OverflowCorrection = OverflowCorrection + T3WRAPAROUND * nsync
					#println("OverflowCorrection =", OverflowCorrection)
                	gotOverflow(nsync, recNum)
					#println("Will call GotOverflow(nsync)")
              	end
            end
            if (channel >= 1) && (channel <= 15)  # these are markers
            	true_nSync = OverflowCorrection + nsync
				
              	gotMarker(true_nSync, channel, recNum)
			  	#println("Will call GotMarker(true_nSync, channel)")
			  	#println("OverflowCorrection =", OverflowCorrection)
				#println("OverflowCorrection =", OverflowCorrection)
            end
        end
    end
end
		
        

# ╔═╡ 79fd185d-0095-415a-8944-bf3eac7b64db
begin
function read_bool8!(fid, tagDataList, evalName)
    tagInt = read(fid, Int64)
    #println("case  tyBool8: TagInt =", tagInt)
    if tagInt == 0
        push!(tagDataList, (evalName, "False"))
    else
        push!(tagDataList, (evalName, "True"))
    end
end

function read_int8!(fid, tagDataList, evalName)
    tagInt = read(fid, Int64)
    #println("case  tyInt8: TagInt =", tagInt)
    push!(tagDataList, (evalName, tagInt))
end

function read_float8!(fid, tagDataList, evalName)
    tagFloat = read(fid, Float64)
    #println("case  tyFloat8: tagFloat =", tagFloat)
    push!(tagDataList, (evalName, tagFloat))
end

function read_empty8!(fid, tagDataList, evalName)
    tagInt = read(fid, Int64)
    #println("case  tyEmpty8: TagInt =", tagInt)
    push!(tagDataList, (evalName, "<empty Tag>"))
end

function read_ansistring!(fid, tagDataList, evalName)
    tagInt = read(fid, Int64)
    #println("case  tyAnsiString: TagInt =", tagInt)
    tagString = read_signature(fid,tagInt)
    #println("TagString =", tagString)
    push!(tagDataList, (evalName, tagString))
end

function read_datetime!(fid, tagDataList, evalName)
    tagFloat = read(fid, Float64)
    #println("case  tyAtyTDateTime: tagFloat =", tagFloat)
    tagTime = unix2datetime(tagFloat)
    #println("tagTime =", tagTime)
    push!(tagDataList, (evalName, tagTime))
end

function read_next(fid)
    #println("function: read_next()")
    tagIdent = read_signature(fid,32)
    #println("tagIdent =", tagIdent)

    tagIdx = read(fid, Int32)
    #println("TagIdx =", tagIdx)

    tagType = read(fid, Int32)
    #println("TagType =", tagType)

    if tagIdx > -1
        evalName = string(tagIdent, "(", string(tagIdx),")")
    else
        evalName = tagIdent
    end
    evalName, tagType, tagIdent
    #println("evalName =", evalName)

end

function read_signature(fid, lword=8)
	#fid = open(zpath, "r")
	sgnt = []
	for i in 1:lword
		append!(sgnt,read(fid, Char))
	end
	#close(fid)
	sgntx = join(string.(sgnt))
	rstrip(sgntx,'\0')
end

hex2dec(hex) = first(reinterpret(Int32, [parse(UInt32, "0x"*hex)]))
end

# ╔═╡ da993660-1647-4d2d-a12d-62e70e5669ac
function test_tpu()
    println("Function test_tpu")

    tyEmpty8      = hex2dec("FFFF0008")
    tyBool8       = hex2dec("00000008")
    tyInt8        = hex2dec("10000008")
    tyBitSet64    = hex2dec("11000008")
    tyColor8      = hex2dec("12000008")
    tyFloat8      = hex2dec("20000008")
    tyTDateTime   = hex2dec("21000008")
    tyFloat8Array = hex2dec("2001FFFF")
    tyAnsiString  = hex2dec("4001FFFF")
    tyWideString  = hex2dec("4002FFFF")
    tyBinaryBlob  = hex2dec("FFFFFFFF")

    println("tyEmpty8 =", tyEmpty8)
    println("tyBool8 =", tyBool8)
    println("tyInt8 =", tyInt8)
    println("tyBitSet64 =", tyBitSet64)
    println("tyColor8 =", tyColor8)
    println("tyFloat8 =", tyFloat8)
    println("tyTDateTime =", tyTDateTime)
    println("tyFloat8Array =", tyFloat8Array)
    println("tyAnsiString =", tyAnsiString)
    println("tyWideString =", tyWideString)
    println("tyBinaryBlob =", tyBinaryBlob)

    # RecordTypes
    rtPicoHarpT3     = hex2dec("00010303")
    rtPicoHarpT2     = hex2dec("00010203")
    rtHydraHarpT3    = hex2dec("00010304")
    rtHydraHarpT2    = hex2dec("00010204")
    rtHydraHarp2T3   = hex2dec("01010304")
    rtHydraHarp2T2   = hex2dec("01010204")
    rtTimeHarp260NT3 = hex2dec("00010305")
    rtTimeHarp260NT2 = hex2dec("00010205")
    rtTimeHarp260PT3 = hex2dec("00010306")
    rtTimeHarp260PT2 = hex2dec("00010206")
    rtMultiHarpT3    = hex2dec("00010307")
    rtMultiHarpT2    = hex2dec("00010207")

    println("\n rtPicoHarpT3 =", rtPicoHarpT3)
    println("rtPicoHarpT2 =", rtPicoHarpT2)
    println("rtHydraHarpT3 =", rtHydraHarpT3)
    println("rtHydraHarpT2 =", rtHydraHarpT2)
    println("rtHydraHarp2T3 =", rtHydraHarp2T3)
    println("rtHydraHarp2T2 =", rtHydraHarp2T2)
    println("rtTimeHarp260NT3 =", rtTimeHarp260NT3)
    println("rtTimeHarp260NT2 =", rtTimeHarp260NT2)
    println("rtTimeHarp260PT3 =", rtTimeHarp260PT3)
    println("rtTimeHarp260PT2 =", rtTimeHarp260PT2)
    println("rtMultiHarpT3 =", rtMultiHarpT3)
    println("rtMultiHarpT2 =", rtMultiHarpT2)

    path="/Users/jjgomezcadenas/LaserLab/Proyectos/data/PMT/hydraharp"
    fname = "test.ptu"
    zpath = joinpath(path,fname)
    fid = open(zpath, "r")


    magic = read_signature(fid,8)
    println("magic =", magic)

    version = read_signature(fid, 8)
    println("version =", version)

    tagDataList = []
    i = 0
    while true 
        i+=1
        
        evalName, tagType, tagIdent = read_next(fid)
        #println("\n---")
        #println("reading cycle = ", i, " evalName =", evalName, " tagType =", tagType)

        if tagType == tyAnsiString
            read_ansistring!(fid, tagDataList, evalName)
        elseif tagType == tyTDateTime
            read_datetime!(fid, tagDataList, evalName)
        elseif tagType == tyInt8
            read_int8!(fid, tagDataList, evalName)
        elseif tagType == tyColor8
            read_int8!(fid, tagDataList, evalName)
        elseif tagType == tyFloat8
            read_float8!(fid, tagDataList, evalName)
        elseif tagType == tyEmpty8
            read_empty8!(fid, tagDataList, evalName)
        elseif tagType == tyBool8
            read_bool8!(fid, tagDataList, evalName)
        
        else
            println("ERROR: Unknown tag type: tagType =", tagType)
            #exit(0)
        end
    
    
        if tagIdent == "Header_End"
            break
        end
    end

    #println("tagDataList =", tagDataList)
    
    tagNames = [tagDataList[i][1] for i in 1:length(tagDataList)]
    tagValues = [tagDataList[i][2] for i in 1:length(tagDataList)]

    tagDict = Dict(tagDataList)
    #println("tagNames =", tagNames)
    #println("tagValues =", tagValues)
    println("\n\n tagDict =", tagDict)

    numRecords = tagDict["TTResult_NumberOfRecords"]
    globRes = tagDict["MeasDesc_GlobalResolution"]
    TTTRTagRes = tagDict["MeasDesc_Resolution"]  # dtime resolution

    println("\n dtime resolution = ", TTTRTagRes)
    println("numRecords = ", numRecords) 

    oflcorrection = 0
    dlen = 0
    recordType = tagDict["TTResultFormat_TTTRRecType"]
    println("recordType = ", recordType) 

    if recordType == rtHydraHarp2T3
        isT2 = false
        println("HydraHarp V2 T3 data")
        #outputfile.write("HydraHarp V2 T3 data\n")
        #outputfile.write("\nrecord# chan   nsync truetime/ns dtime\n")
        readHT3(2, fid, numRecords, oflcorrection, isT2, globRes)

	else 
		println("Not yet implemented")
    end

    close(fid)

end



# ╔═╡ 843c6950-d629-495b-9d9e-06d920d80229
test_tpu()

# ╔═╡ Cell order:
# ╠═1459f36b-afa4-405e-a9e8-e6ae32ea752e
# ╠═0735b900-2548-412f-a675-c3c19945ae47
# ╠═c83e18f6-8d4b-4fd7-ad3c-950c4eb4b91d
# ╠═ec3057cd-41a8-4c87-8808-5f612bfa82a5
# ╠═e6550506-f9a3-4999-9fef-4e0cf70cc1ec
# ╠═f9470c2c-4386-41b5-9473-f88fb713a0bb
# ╠═73f68dd8-863d-49fe-9166-41fbd0ec36cb
# ╠═6753e916-8645-4256-ab5a-191e6cbca9df
# ╠═b6517e62-c865-46aa-a43f-fa3075760309
# ╠═8f01ac68-1cbf-4341-aabc-241bc450f869
# ╠═e9f05987-8092-4164-a6a9-e0a041be123c
# ╠═4ca63a99-8679-4f3c-be4a-049442499e77
# ╠═511ee273-3401-4a4e-af28-d19b24c26642
# ╠═9f32a075-eab9-421f-801f-eabaa9400251
# ╠═e6a6aea1-93ab-4b75-98ef-c073e7310b83
# ╠═e71daa46-42dd-4b5b-ac40-27e66ffd4d6d
# ╠═b4941210-d353-441d-8d51-d537c8fcb16e
# ╠═62e60c02-321f-47bb-94a6-c59e61824ab2
# ╠═70a97b9a-dab1-4a04-a2ac-2123cf7da83a
# ╠═f7ae7992-4a66-4f55-a5db-d8dcdcaed7a3
# ╠═b026d200-1763-4ea0-8336-24fc9a06dd4f
# ╠═28936553-406e-433f-8a4b-167ca0952af8
# ╠═8e7c68c7-ff2e-4895-aaf7-189f1d9eb39d
# ╠═843c6950-d629-495b-9d9e-06d920d80229
# ╠═da993660-1647-4d2d-a12d-62e70e5669ac
# ╠═e58d12e4-b51e-4bf4-a96b-2e538dc03740
# ╠═0a064d7b-64da-4038-80e1-d0f8227a9052
# ╠═1feb11fc-2f54-4c75-a6a0-61b9435bd411
# ╠═7f0b9b19-e92b-4f15-8524-5ddee7d18607
# ╠═79fd185d-0095-415a-8944-bf3eac7b64db
