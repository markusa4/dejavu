import subprocess, csv
from os import walk

setname = "results/rantree"
# solver  = "dejavu7"

def process(setname, solver):
	csvfile = "{}.{}.dat".format(setname, solver)
	f = open(csvfile)
	reader = csv.reader(f, delimiter = " ")
	sum = 0.0;
	avg_sum = {}
	avg_cnt = {}

	for line in reader:
		sum += float(line[2])
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

	sumfile.close()

process(setname, "dejavu7")
process(setname, "dejavu4")
process(setname, "dejavu1")
process(setname, "traces")
process(setname, "nauty")
