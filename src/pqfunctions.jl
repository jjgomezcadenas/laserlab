import Dates
#using PythonStructs

function readHH(zpath, nevents, runall=false, verbose=false, summary=true, empty=false)
if summary
    println("Function readHH")
end

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

    if verbose
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
    end

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

    if verbose
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
    end

    fid = open(zpath, "r")

    magic = read_signature(fid,8)
    #println("magic =", magic)

    if magic != "PQTTTR"
        println("ERROR: Magic invalid, this is not a PTU file.")
        close(fid)
        exit(0)
    end

    version = read_signature(fid, 8)
    if summary 
        println("version =", version)
    end

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
            exit(0)
        end
    
    
        if tagIdent == "Header_End"
            break
        end
    end

    tagDict = Dict(tagDataList)
    
    if verbose
        println("\n\n tagDict =", tagDict)
    end

    numRecords = tagDict["TTResult_NumberOfRecords"]
    globRes    = tagDict["MeasDesc_GlobalResolution"]
    TTTRTagRes = tagDict["MeasDesc_Resolution"]  # dtime resolution
    recordType = tagDict["TTResultFormat_TTTRRecType"]

    if summary
        println("\n global resolution (ns)= ", globRes * 1e+9)
        println("\n dtime resolution  (ns)= ", TTTRTagRes * 1e+9)
        println(" numRecords = ", numRecords)     
        println(" recordType = ", recordType) 
    end

	if runall
		nevents = numRecords
	end
	
    if recordType == rtHydraHarp2T3
    	println("HydraHarp V2 T3 data")
        println("running number of events = ", nevents) 
        if empty 
            return nevents
        else

            PhotonDF, LaserDF = readHT3(fid, 2, nevents)
        end

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


"""
Reads hydraHarp (T3 format) data

"""
function readHT3(fid, version::Integer,  numRecords::Integer, verbose=false)
   
    function read_and_parse_record()
        try 
            t3record   = read(fid, UInt32)                       # reads a 32 bits number
            recordData = bitstring(t3record)                     # bitstring representation
            special    = parse(Int32, recordData[1], base=2)     # last bit
            channel    = parse(Int32, recordData[2:7], base=2)   # next 6 bits
            dtime      = parse(Int32, recordData[8:22], base=2)  # next 15 bits
            nsync      = parse(Int32, recordData[23:32], base=2) #lowest 10 bits
            special, channel, dtime, nsync
        catch
            println("Error reading data at record =>", recNum)
            exit(0) 
        end

    end
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
		special, channel, dtime, nsync = read_and_parse_record()

        if verbose
            println("\n recNum = ", recNum)
            println(" special = ", special)
            println(" channel = ", channel)
            println(" dtime = ", dtime)
            println(" nsync = ", nsync)
        end

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
                        println("gotOverflow: recNum  = ", recNum, 
                                " oflcorrection =", oflcorrection, " count =", count)
                    end
				else                              # otherwise nsync indicates the number of 
				                                  # overflows 
				                                  # This corresponds to format v2
					
                	oflcorrection = oflcorrection + T3WRAPAROUND * nsync
                    ount = 1
                    if verbose 
                        println("gotOverflow: recNum  = ", recNum, 
                                " oflcorrection =", oflcorrection, " count =", count)
                    end
              	end
            end
            if (channel >= 1) && (channel <= 15)  # these are markers
            	true_nSync = oflcorrection + nsync
                
                #printlng("Got a marker: channel =", channel)
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
	

"""
Reads hydraHarp (T2 format) data

"""
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


"""
Reads an integer and interprets it as a bool: push the result into tagDataList
"""
function read_bool8!(fid, tagDataList, evalName)
    tagInt = read(fid, Int64)
    #println("case  tyBool8: TagInt =", tagInt)
    if tagInt == 0
        push!(tagDataList, (evalName, "False"))
    else
        push!(tagDataList, (evalName, "True"))
    end
end


"""
Read an integer: push the result into tagDataList
"""
function read_int8!(fid, tagDataList, evalName)
    tagInt = read(fid, Int64)
    #println("case  tyInt8: TagInt =", tagInt)
    push!(tagDataList, (evalName, tagInt))
end


"""
Read a float: push the result into tagDataList
"""
function read_float8!(fid, tagDataList, evalName)
    tagFloat = read(fid, Float64)
    #println("case  tyFloat8: tagFloat =", tagFloat)
    push!(tagDataList, (evalName, tagFloat))
end

"""
Read an empty record: push the result into tagDataList
"""
function read_empty8!(fid, tagDataList, evalName)
    tagInt = read(fid, Int64)
    #println("case  tyEmpty8: TagInt =", tagInt)
    push!(tagDataList, (evalName, "<empty Tag>"))
end


"""
Read an ANSI string: push the result into tagDataList
"""
function read_ansistring!(fid, tagDataList, evalName)
    tagInt = read(fid, Int64)
    #println("case  tyAnsiString: TagInt =", tagInt)
    tagString = read_signature(fid,tagInt)
    #println("TagString =", tagString)
    push!(tagDataList, (evalName, tagString))
end


"""
Read a date/time: push the result into tagDataList
"""
function read_datetime!(fid, tagDataList, evalName)
    tagFloat = read(fid, Float64)
    tagTime = ((tagFloat - 25569) * 86400)
    #println("case  tyAtyTDateTime: tagFloat =", tagFloat, " tagTime =", tagTime)
    tagTime = Dates.unix2datetime(tagTime)
    #println("tagTime =", tagTime)
    push!(tagDataList, (evalName, tagTime))
end

"""
Reads the stream
"""
function read_next(fid)
    #println("function: read_next()")
    tagIdent = read_signature(fid,32)
    tagIdx = read(fid, Int32)
    tagType = read(fid, Int32)
    

    if tagIdx > -1
        evalName = string(tagIdent, "(", string(tagIdx),")")
    else
        evalName = tagIdent
    end
    evalName, tagType, tagIdent

end

"""
Reads lword characters and joins them into a string
"""
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


"""
Yields the decimal representation of a hex
"""
hex2dec(hex) = first(reinterpret(Int32, [parse(UInt32, "0x"*hex)]))