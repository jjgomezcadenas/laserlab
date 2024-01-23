
using Peaks
using Glob
using FFTW
using DSP
using DataFrames



inputfile = open(sys.argv[1], "rb")
# The following is needed for support of wide strings
outputfile = io.open(sys.argv[2], "w+", encoding="utf-16le")

# Check if inputfile is a valid PTU file
# Python strings don't have terminating NULL characters, so they're stripped
magic = inputfile.read(8).decode("utf-8").strip("\0")
if magic != "PQTTTR":
    print("ERROR: Magic invalid, this is not a PTU file.")
    inputfile.close()
    outputfile.close()
    exit(0)

version = inputfile.read(8).decode("utf-8").strip("\0")
outputfile.write("Tag version: %s\n" % version)


% start Main program
    [filename, pathname]=uigetfile('*.ptu', 'T-Mode data:');
    fid=fopen([pathname filename]);

    fprintf(1,'\n');
    Magic = fread(fid, 8, '*char');
    if not(strcmp(Magic(Magic~=0)','PQTTTR'))
        error('Magic invalid, this is not an PTU file.');
    end;
    Version = fread(fid, 8, '*char');
    fprintf(1,'Tag Version: %s\n', Version);
    
"""
Reads  waveform specifying directory and file name 
"""
function read_ptu(pmtdir::String, fname::String)
    zpath = joinpath(pmtdir,fname)
	zio = open(zpath, "r")
end

