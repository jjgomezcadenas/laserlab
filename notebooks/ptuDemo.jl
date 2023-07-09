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

# ╔═╡ 74a2ea89-0b20-4d71-9461-82755b6ed6d8
using PythonStructs

# ╔═╡ 04e7c8ad-d693-49b0-aa30-cbef6aa35b9b
using EasyFit

# ╔═╡ 9657e6f5-b3b9-448f-8bb5-28d5a64c4c37
function ingredients(path::String)
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol(basename(path))
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ ea788fbe-160c-40fd-b967-d49bc5403398
lfi = ingredients("../src/LaserLab.jl")

# ╔═╡ 8e7c68c7-ff2e-4895-aaf7-189f1d9eb39d
md"""

## Analysis
"""

# ╔═╡ b605ebe1-ba55-422a-9223-87948a45da4c
begin
	patht2= "/Users/jjgomezcadenas/LaserLab/Proyectos/data/APD/IrSL_VS026_B1(2.5e-6M)/Filter0"
	fnamet2 = "T3.ptu"
	patht3="/Users/jjgomezcadenas/Projects/LaserLab/test"
	fnamet3 = "test.ptu"
	zpath = joinpath(patht2,fnamet2)
end

# ╔═╡ 0b0bb96c-0f32-41dd-95b4-9df7c5550b1e
b_range = range(0, 500, length=250)

# ╔═╡ 8c68cd00-2281-4e10-8535-802cd41825b8


# ╔═╡ f769bfad-8dc0-40ce-8f49-c31aa1a94422
get_histo1d(edges::Vector{Float64}, weights::Vector{Float64}, centers::Vector{Float64}) = lfi.LaserLab.Histo1d(edges, weights, centers)

# ╔═╡ 045b2ca2-00dd-468c-bddb-c80edddc6127
function h1d(x::Vector{T}, nbins::Integer,
             xmin::T=typemin(T), xmax::T=typemax(T), norm=false) where T
	dx = (xmax - xmin) /nbins
	#edges =Vector{T}(undef, length(x)+1)
	
	edges = [xmin + i * dx for i in 0:nbins]
	h = fit(Histogram, x, edges)
	println(findall( x -> x == 0, h.weights ))
	lfi.LaserLab.Histo1d(lfi.LaserLab.edges(h), h.weights, lfi.LaserLab.centers(h))
end

# ╔═╡ 400cc8ad-f9b7-4aa7-a5c6-73daed547df2
function plot_h1d(h, xs::String; i0=1, il=-1, 
                  ylog::Bool=false, xlog::Bool=false,
                  markersize::Int64=3, fwl::Bool=false,
                  label::String="", legend::Bool=false)

    xg = h.centers
    yg = h.weights

    if il == -1 
        il = length(h.centers)
    end

    xaxis = (xs, (h.edges[i0],h.edges[il]))
	yaxis = ("a.u.", (10^3, 10^5), :log)
	marker = (:circle, markersize, 0.6, :black)
    

    p  = scatter(xg[i0:il], yg[i0:il], yerr=sqrt.(abs.(yg[i0:il])), marker = marker,
                 label = label, legend=legend, yaxis=yaxis, xaxis =xaxis)
    if fwl
        p = plot!(p, xg[i0:il], yg[i0:il], yerr=sqrt.(abs.(yg[i0:il])), fmt=:png,
                  linewidth=1, label=label, legend=legend, yaxis=:log, xaxis =(20,490))
    end
	p
end

# ╔═╡ b40c0a89-8b43-4924-a9b8-823dc8c9449f
function tfit(htime, pa0, i0, l0, ffun, mffun)
	
	tdata = htime.centers
	vdata = htime.weights
	
	#fitb = curve_fit(mffun, tdata[i0:l0], vdata[i0:l0], pa0)

	fit = curve_fit(mffun, tdata[i0:l0], vdata[i0:l0], pa0)
	coef(fit), stderror(fit), ffun.(tdata, coef(fit)...)
end

# ╔═╡ 41e1c822-de92-4744-88d6-6cfa7331b001
function tfit_2exp(htime; pa0=[1000.0, 10., 1000.0, 10.0], i0, l0)
	
	ffun(t, N1, λ1, N2, λ2) = N1 * exp(-t/λ1) + N2 * exp(-t/λ2)
	mffun(t, p) = p[1] * exp.(-t/p[2]) + p[3] * exp.(-t/p[4]) 

	tfit(htime, pa0, i0, l0, ffun, mffun)
end

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
    PTimeTag = Int64[]
    PDtime = Int32[]
    MRecNum = Int32[]
    MChannel = Int32[]
    MTimeTag = Int64[]
    
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
	

# ╔═╡ a6142bb4-e0ce-44e6-8688-e7bcaae38eb0
function readHT2(fid, version::Integer,  numRecords::Integer; verbose=false)
	"""
	Reads a 32 integer and return the bit representation
	"""
	function read_record(fid)
		try 
			t3record   = read(fid, UInt32)   # reads a 32 bits number
			bitstring(t3record)              # string representation (10001...)
		catch
			println("Failed reading: record = ", recNum)
		end 
		
	end

	"""
	Parse record
	"""
	function parse_record(recordData)
		try 
			special    = parse(Int32, recordData[1], base=2)     # last bit
			channel    = parse(Int32, recordData[2:7], base=2)   # next 6 bits
			timetag    = parse(Int32, recordData[8:32], base=2)  # next 25 bits
			special, channel, timetag
		catch
			println("Failed parsing: record = ", recNum)
			println("Failed parsing: recordData = ", recordData)
		end
	end
   
    T2WRAPAROUND_V1 = 33552000
    T2WRAPAROUND_V2 = 33554432
	oflcorrection = 0

    PRecNum = Int32[]
    PChannel = Int32[]
    PTimeTag = Int64[]
    MRecNum = Int32[]
    MChannel = Int32[]
    MTimeTag = Int64[]
	
    
	for recNum in 1:numRecords
		recordData   = read_record(fid)    
		special, channel, timetag = parse_record(recordData)
		                                        
		if special == 1
            if channel == 63  # Overflow
                                # Number of overflows in nsync. 
								# If old version, it's an
                                # old style single overflow
                if version == 1
                    oflcorrection += T2WRAPAROUND_V1
                else
                    if timetag == 0  # old style overflow, shouldn't happen
                        oflcorrection += T2WRAPAROUND_V2
                    else
                        oflcorrection += T2WRAPAROUND_V2 * timetag
					end
				end
			end
                        
            if channel >= 1 && channel <= 15  # markers
                truetime = oflcorrection + timetag
				push!(MRecNum, recNum)
                push!(MChannel, channel)
				push!(MTimeTag, truetime)
			end
				
            if channel == 0  # sync
				#if recNum > 4360 && recNum < 4380
				#	println("recNum = ", recNum)
				#	println("recordData = ", recordData)
				#	println("special = ", special, " channel=", channel, " timetag=", 
				 #  timetag)
					#println("type of tag =", typeof(timetag))
				#end
                truetime = oflcorrection + timetag
				#if recNum > 4300
				#	println("truetime =", truetime, " type =", typeof(truetime))
				#end
                push!(PRecNum, recNum)
            	push!(PChannel, channel)
            	push!(PTimeTag, truetime)
			end
        else  # regular input channel
            truetime = oflcorrection + timetag
			push!(PRecNum, recNum)
            push!(PChannel, channel + 1)
            push!(PTimeTag, truetime)
           
		end
    end
	PhotonDF = Dict("recNum"=>PRecNum,   "channel"=>PChannel, 
		            "timeTag"=>PTimeTag)
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
		tagTime = ((tagFloat - 25569) * 86400)
	    #println("case  tyAtyTDateTime: tagFloat =", tagFloat, " tagTime =", tagTime)
	    tagTime = unix2datetime(tagTime)
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
function readHH(zpath, nevents, runall=false)
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

    #path="/Users/jjgomezcadenas/Projects/LaserLab/test"
    #fname = "test.ptu"
    #zpath = joinpath(path,fname)
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

	elseif recordType == rtHydraHarp2T2
		println("HydraHarp V2 T2 data")
        println("running number of events = ", nevents) 
        PhotonDF, LaserDF = readHT2(fid, 2, nevents)
	else 
		println("Not yet implemented")
    end

    close(fid)
	PhotonDF, tagDict
end

# ╔═╡ 843c6950-d629-495b-9d9e-06d920d80229
PhotonDF, tagDict = readHH(zpath, 2, true)

# ╔═╡ 2d500c48-b3aa-492a-ad36-ae566264a1eb
tns = tagDict["MeasDesc_GlobalResolution"] * 1e+9

# ╔═╡ 1c55bf06-2813-4b9e-bee8-2b4920a897c5
ttrns = tagDict["MeasDesc_Resolution"]* 1e+9

# ╔═╡ 1cedaad8-c2ae-4fba-8326-967d9e37ea9b
h1dtns  = h1d(PhotonDF.dtime * ttrns, 244, 0.0, 500.0)

# ╔═╡ 5f143f8b-b72f-4c21-8bfc-17cf061653c9
pp = plot_h1d(h1dtns, "dt (ns)"; i0=22, il=243, markersize=2, fwl=true)

# ╔═╡ bef4551d-e6cc-4228-9713-3ac216b0741e
begin
	nmin = 100
	tdata = h1dtns.centers
	vdata = h1dtns.weights
	ft2e = fitexp(tdata[nmin:243], vdata[nmin:243],n=1)
end

# ╔═╡ 4b81c186-337e-4080-950c-52bcd099e8f8
begin
	nminx = 20
	nmax = 40
	ft1e = fitexp(tdata[nminx:nmax], vdata[nminx:nmax],n=1)
end

# ╔═╡ a3593422-04f3-45e0-9747-c73a99a27226
pp2 = plot(pp, tdata[nmin:243], ft2e.ypred)

# ╔═╡ c1242a5a-8637-4e7e-b9e3-f5ba399fc151
plot(tdata[nminx:nmax], ft1e.ypred)

# ╔═╡ 843186c5-befe-4465-9fdb-81a7879db8fd
plot(pp, tdata[nminx:nmax], ft1e.ypred)

# ╔═╡ Cell order:
# ╠═1459f36b-afa4-405e-a9e8-e6ae32ea752e
# ╠═e6550506-f9a3-4999-9fef-4e0cf70cc1ec
# ╠═74a2ea89-0b20-4d71-9461-82755b6ed6d8
# ╠═04e7c8ad-d693-49b0-aa30-cbef6aa35b9b
# ╠═9657e6f5-b3b9-448f-8bb5-28d5a64c4c37
# ╠═ea788fbe-160c-40fd-b967-d49bc5403398
# ╠═8e7c68c7-ff2e-4895-aaf7-189f1d9eb39d
# ╠═b605ebe1-ba55-422a-9223-87948a45da4c
# ╠═843c6950-d629-495b-9d9e-06d920d80229
# ╠═2d500c48-b3aa-492a-ad36-ae566264a1eb
# ╠═1c55bf06-2813-4b9e-bee8-2b4920a897c5
# ╠═0b0bb96c-0f32-41dd-95b4-9df7c5550b1e
# ╠═8c68cd00-2281-4e10-8535-802cd41825b8
# ╠═f769bfad-8dc0-40ce-8f49-c31aa1a94422
# ╠═045b2ca2-00dd-468c-bddb-c80edddc6127
# ╠═1cedaad8-c2ae-4fba-8326-967d9e37ea9b
# ╠═400cc8ad-f9b7-4aa7-a5c6-73daed547df2
# ╠═5f143f8b-b72f-4c21-8bfc-17cf061653c9
# ╠═bef4551d-e6cc-4228-9713-3ac216b0741e
# ╠═4b81c186-337e-4080-950c-52bcd099e8f8
# ╠═a3593422-04f3-45e0-9747-c73a99a27226
# ╠═c1242a5a-8637-4e7e-b9e3-f5ba399fc151
# ╠═843186c5-befe-4465-9fdb-81a7879db8fd
# ╠═41e1c822-de92-4744-88d6-6cfa7331b001
# ╠═b40c0a89-8b43-4924-a9b8-823dc8c9449f
# ╠═adef1149-cbf5-4403-88b0-e51196030cec
# ╠═da993660-1647-4d2d-a12d-62e70e5669ac
# ╠═e58d12e4-b51e-4bf4-a96b-2e538dc03740
# ╠═a6142bb4-e0ce-44e6-8688-e7bcaae38eb0
# ╠═79fd185d-0095-415a-8944-bf3eac7b64db
