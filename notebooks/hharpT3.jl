### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 1459f36b-afa4-405e-a9e8-e6ae32ea752e
using Pkg; Pkg.activate(ENV["JLaserLab"])

# ╔═╡ e6550506-f9a3-4999-9fef-4e0cf70cc1ec
begin
	using PlutoUI
	using CSV
	using DataFrames
	using Dates
	using Plots
	using Printf
	using Markdown
	using InteractiveUtils
	using LsqFit
	using Statistics
	using StatsBase
	using Unitful 
	using UnitfulEquivalences 
	using PhysicalConstants
end

# ╔═╡ 8e7c68c7-ff2e-4895-aaf7-189f1d9eb39d
md"""

## Analysis
"""

# ╔═╡ 0b0bb96c-0f32-41dd-95b4-9df7c5550b1e
b_range = range(0, 2000, length=200)

# ╔═╡ adef1149-cbf5-4403-88b0-e51196030cec
md"""
## Functions
"""

# ╔═╡ e58d12e4-b51e-4bf4-a96b-2e538dc03740
function readHT3(fid, version::Integer,  numRecords::Integer; verbose=false)
   
    T3WRAPAROUND = 1024
	oflcorrection = 0

	
    PRecNum = Int32[]
    PChannel = Int32[]
    PTimeTag = Int32[]
    PDtime = Int32[]
    MRecNum = Int32[]
    MChannel = Int32[]
    MTimeTag = Int32[]
    
	for recNum in 1:numRecords
		t3record   = read(fid, UInt32)                       # reads a 32 bits number
		recordData = bitstring(t3record)                     # bitstring 
		                                                     # representation
		special    = parse(Int32, recordData[1], base=2)     # last bit
		channel    = parse(Int32, recordData[2:7], base=2)   # next 6 bits
		dtime      = parse(Int32, recordData[8:22], base=2)  # next 15 bits
		nsync      = parse(Int32, recordData[23:32], base=2) #lowest 10 bits

		#println("\n recNum = ", recNum)
		#println(" recordData = ", recordData)
		#println(" special = ", special)
		#println(" channel = ", channel)
		#println(" dtime = ", dtime)
		#println(" nsync = ", nsync)

        
		if special == 0   # this means a regular input channel
        	true_nSync = oflcorrection + nsync
            
            push!(PRecNum, recNum)
            push!(PChannel, channel)
            push!(PTimeTag, true_nSync)
            push!(PDtime, dtime)

			
		else                                      # this means we have a special record
        	if channel == 63                      # overflow of nsync occured
            	if (nsync == 0) || (version == 1) # if nsync is zero it is 
				                                  # an old style      
												  # single oferflow or old Version
                	oflcorrection = oflcorrection + T3WRAPAROUND
                    count = 1

					if verbose 
                    	println(" gotOverflow: recNum  = ", recNum, 
						        " oflcorrection =", oflcorrection, 
						        " count =", count)
					end
                	
				else                              # otherwise nsync indicates the
				                                  # number of 
				                                  # overflows 
				                                  # This corresponds to format v2
					
                	oflcorrection = oflcorrection + T3WRAPAROUND * nsync

					if verbose
                   		println(" gotOverflow: recNum  = ", recNum, 
						        " oflcorrection =", oflcorrection, 
						        " count =", nsync)
					end
              	end
            end
            if (channel >= 1) && (channel <= 15)  # these are markers
            	true_nSync = oflcorrection + nsync
                
                push!(MRecNum, recNum)
                push!(MChannel, channel)
                push!(MTimeTag, true_nSync)
                  
            end
        end
    end
	PhotonDF = Dict("recNum"=>PRecNum,   "channel"=>PChannel, 
		            "timeTag"=>PTimeTag, "dtime"=>PDtime)
	LaserDF = Dict("recNum"=>MRecNum,    "channel"=>MChannel, 
		           "timeTag"=>MTimeTag)
	
	DataFrame(PhotonDF), DataFrame(LaserDF)
				  
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
function readHH(nevents, runall=false)
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

    # println("tyEmpty8 =", tyEmpty8)
    # println("tyBool8 =", tyBool8)
    # println("tyInt8 =", tyInt8)
    # println("tyBitSet64 =", tyBitSet64)
    # println("tyColor8 =", tyColor8)
    # println("tyFloat8 =", tyFloat8)
    # println("tyTDateTime =", tyTDateTime)
    # println("tyFloat8Array =", tyFloat8Array)
    # println("tyAnsiString =", tyAnsiString)
    # println("tyWideString =", tyWideString)
    # println("tyBinaryBlob =", tyBinaryBlob)

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

    # println("\n rtPicoHarpT3 =", rtPicoHarpT3)
    # println("rtPicoHarpT2 =", rtPicoHarpT2)
    # println("rtHydraHarpT3 =", rtHydraHarpT3)
    # println("rtHydraHarpT2 =", rtHydraHarpT2)
    # println("rtHydraHarp2T3 =", rtHydraHarp2T3)
    # println("rtHydraHarp2T2 =", rtHydraHarp2T2)
    # println("rtTimeHarp260NT3 =", rtTimeHarp260NT3)
    # println("rtTimeHarp260NT2 =", rtTimeHarp260NT2)
    # println("rtTimeHarp260PT3 =", rtTimeHarp260PT3)
    # println("rtTimeHarp260PT2 =", rtTimeHarp260PT2)
    # println("rtMultiHarpT3 =", rtMultiHarpT3)
    # println("rtMultiHarpT2 =", rtMultiHarpT2)

    path="/Users/jjgomezcadenas/Projects/LaserLab/test"
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
    
    #tagNames = [tagDataList[i][1] for i in 1:length(tagDataList)]
    #tagValues = [tagDataList[i][2] for i in 1:length(tagDataList)]

    tagDict = Dict(tagDataList)
    #println("tagNames =", tagNames)
    #println("tagValues =", tagValues)
    println("\n\n tagDict =", tagDict)

    numRecords = tagDict["TTResult_NumberOfRecords"]
    globRes    = tagDict["MeasDesc_GlobalResolution"]
    TTTRTagRes = tagDict["MeasDesc_Resolution"]  # dtime resolution
    recordType = tagDict["TTResultFormat_TTTRRecType"]

    println("\n global resolution (ns)= ", globRes * 1e+9)
    println("\n dtime resolution  (ns)= ", TTTRTagRes * 1e+9)
    println(" numRecords = ", numRecords)     
    println(" recordType = ", recordType) 

	if runall
		nevents = numRecords
	end
	
    if recordType == rtHydraHarp2T3
    	println("HydraHarp V2 T3 data")
        println("running number of events = ", nevents) 
        PhotonDF, LaserDF = readHT3(fid, 2, nevents)
	else 
		println("Not yet implemented")
    end

    close(fid)
	PhotonDF, LaserDF
end

# ╔═╡ 843c6950-d629-495b-9d9e-06d920d80229
PhotonDF, LaserDF = readHH(100, true)

# ╔═╡ 1896e5fb-4fcc-4e7c-9cd0-5acc4704682a
histogram(PhotonDF.dtime, label="D time", bins=b_range)

# ╔═╡ 54102c9d-2265-4d80-b896-d0fa263dc53c
histogram(PhotonDF.dtime, label="D time/cut", bins=range(300, 1500, length=40))

# ╔═╡ Cell order:
# ╠═1459f36b-afa4-405e-a9e8-e6ae32ea752e
# ╠═e6550506-f9a3-4999-9fef-4e0cf70cc1ec
# ╠═8e7c68c7-ff2e-4895-aaf7-189f1d9eb39d
# ╠═843c6950-d629-495b-9d9e-06d920d80229
# ╠═0b0bb96c-0f32-41dd-95b4-9df7c5550b1e
# ╠═1896e5fb-4fcc-4e7c-9cd0-5acc4704682a
# ╠═54102c9d-2265-4d80-b896-d0fa263dc53c
# ╠═adef1149-cbf5-4403-88b0-e51196030cec
# ╠═da993660-1647-4d2d-a12d-62e70e5669ac
# ╠═e58d12e4-b51e-4bf4-a96b-2e538dc03740
# ╠═79fd185d-0095-415a-8944-bf3eac7b64db
