import subprocess, csv
from os import walk

def process(setname, solver):
	print(setname)
	print(solver)
	csvfile = "{}.{}.dat".format(setname, solver)
	f = open(csvfile)
	reader = csv.reader(f, delimiter = " ")
	sum = 0.0
	sum100 = 0.0
	cnt = 0
	avg_sum = {}
	avg_cnt = {}

	for line in reader:
		sum += float(line[2])
		cnt += 1
		if(int(line[1]) >= 100):
			sum100 += float(line[2])
		if(avg_sum.get(int(line[1]), "") == ""):
			avg_sum[int(line[1])] = 0.0
			avg_cnt[int(line[1])] = 0
		avg_sum[int(line[1])] += float(line[2])
		avg_cnt[int(line[1])] += 1
	vs = sorted(avg_sum.keys())
	avgfilename = "{}.{}.avg".format(setname, solver)
	avgfile     = open(avgfilename,'w')
	avgfile.write("V time\n") 
	for v in vs: 
	    avgfile.write("{} {}\n".format(v, avg_sum[v] / avg_cnt[v])) 
	avgfile.close()
	sumfilename = "{}.{}.sum".format(setname, solver)
	sumfile     = open(sumfilename,'w')
	sumfile.write("{}".format(sum))
	sumfilename100 = "{}.{}.sum100".format(setname, solver)
	sumfile100     = open(sumfilename100,'w')
	sumfile100.write("{}".format(sum100))
	cntfilename = "{}.{}.cnt".format(setname, solver)
	cntfile     = open(cntfilename,'w')
	cntfile.write("{}".format(cnt))
	print(cnt)
	sumfile.close()

def process_set(setname):
	process(setname, "dejavu8")
	process(setname, "dejavu7")
	process(setname, "dejavu4")
	process(setname, "dejavu2")
	process(setname, "dejavu1")
	process(setname, "traces")
	process(setname, "nauty")

# process_set("results/lattice")
# process_set("results/grid")
# process_set("results/ag")
#process_set("results/latin")
#process_set("results/k")
#process_set("results/latin-sw")
#process_set("results/sts")
#process_set("results/sts-sw")
#process_set("results/had")
#process_set("results/had-sw")
#process_set("results/cfi")
#process_set("results/rantree")
#process_set("results/ranreg")
#process_set("results/hypercubes")
#process_set("results/dac_pipe")
process_set("results/ran2")
process_set("results/ran10")
process_set("results/ransq")
process_set("results/dac_other")
process_set("results/tran")
process_set("results/combinatorial")
process_set("results/pp16")
process_set("results/pp25")
