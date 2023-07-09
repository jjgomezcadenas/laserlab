println("Hello")
path="/Users/jjgomezcadenas/LaserLab/Proyectos/data/PMT/hydraharp"
fname = "test.ptu"
zpath = joinpath(path,fname)
fid = open(zpath, "r")


magic = read_signature(8)
close(fid)

println("magic =", magic)

function read_signature(lword=8)
	#fid = open(zpath, "r")
	sgnt = []
	for i in 1:lword
		append!(sgnt,read(fid, Char))
	end
	#close(fid)
	sgntx = join(string.(sgnt))
	rstrip(sgntx,'\0')
end