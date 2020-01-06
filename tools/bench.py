import subprocess
from datetime import datetime
from random import seed
from random import randint
from os import walk

def bench(dir, setname, timeout, nauty_timeout, selector):
	f     = []
	seeds = []
	seed(datetime.now())
	for (dirpath, dirnames, filenames) in walk(dir):
	    f.extend(filenames)
	    break
	for i in range(len(f)):
		seeds += [randint(0, 100000)]
	for i in range(len(f)):
		subprocess.run("nice -15 ./bench --permute_seed {} --timeout {} --threads 7 --no_nauty --no_traces --file \"{}{}\" --stat_file {}.dejavu8.dat".format(seeds[i], timeout, dir, f[i], setname), shell=True)
	for i in range(len(f)):
		subprocess.run("nice -15 ./bench --permute_seed {} --timeout {} --threads 6 --no_nauty --no_traces --file \"{}{}\" --stat_file {}.dejavu7.dat".format(seeds[i], timeout, dir, f[i], setname), shell=True)
	for i in range(len(f)):
		subprocess.run("nice -15 ./bench --permute_seed {} --timeout {} --threads 3 --no_nauty --no_traces --file \"{}{}\" --stat_file {}.dejavu4.dat".format(seeds[i], timeout, dir, f[i], setname), shell=True)
	for i in range(len(f)):
		subprocess.run("nice -15 ./bench --permute_seed {} --timeout {} --threads 1 --force_selector {} --no_nauty --no_traces --file \"{}{}\" --stat_file {}.dejavu2.dat".format(seeds[i], timeout, selector, dir, f[i], setname), shell=True)
	for i in range(len(f)):
		subprocess.run("nice -15 ./bench --permute_seed {} --timeout {} --threads 0 --force_selector {} --no_nauty --no_traces --file \"{}{}\" --stat_file {}.dejavu1.dat".format(seeds[i], timeout, selector, dir, f[i], setname), shell=True)
	for i in range(len(f)):
		subprocess.run("nice -15 ./bench --permute_seed {} --timeout {} --no_nauty --no_dejavu --file \"{}{}\" --stat_file {}.traces.dat".format(seeds[i], timeout, dir, f[i], setname), shell=True)
	for i in range(len(f)):
		subprocess.run("nice -15 ./bench --permute_seed {} --timeout {} --no_traces --no_dejavu --file \"{}{}\" --stat_file {}.nauty.dat".format(seeds[i], nauty_timeout, dir, f[i], setname), shell=True)

bench("./graphs/undirected_dim/undirected_dim/k/k/", "results/k", -1, -1, 1)
bench("./graphs/undirected_dim/undirected_dim/grid/grid/", "results/grid", -1, -1, 1)
bench("./graphs/undirected_dim/undirected_dim/lattice/lattice/", "results/lattice", -1, -1, 1)
bench("./graphs/undirected_dim/undirected_dim/ag/ag/", "results/ag", -1, -1, 1)
bench("./graphs/undirected_dim/undirected_dim/latin/latin/", "results/latin", -1, -1, 1)
bench("./graphs/undirected_dim/undirected_dim/latin-sw/latin-sw/", "results/latin-sw", -1, -1, 2)
bench("./graphs/undirected_dim/undirected_dim/sts/sts/", "results/sts", -1, -1, 2)
bench("./graphs/undirected_dim/undirected_dim/sts-sw/sts-sw/", "results/sts-sw", -1, -1, 2)
bench("./graphs/undirected_dim/undirected_dim/had/had/", "results/had", -1, 30, 1)
bench("./graphs/undirected_dim/undirected_dim/had-sw/had-sw/", "results/had-sw", -1, 30, 1)
bench("./graphs/undirected_dim/undirected_dim/cfi/cfi/", "results/cfi", -1, -1, 1)
bench("./graphs/ran2/ran2/", "results/ran2", -1, -1, 1)
bench("./graphs/ran10/ran10/", "results/ran10", -1, -1, 1)
bench("./graphs/ransq/ransq/", "results/ransq", -1, -1, 1)
bench("./graphs/dac/dac/other/", "results/dac_other", -1,-1, 2)
bench("./graphs/tran/tran/", "results/tran", -1, -1, 1)
# double
bench("./graphs/undirected_dim/undirected_dim/k/k/", "results/k", -1, -1, 1)
bench("./graphs/undirected_dim/undirected_dim/grid/grid/", "results/grid", -1, -1, 1)
bench("./graphs/undirected_dim/undirected_dim/lattice/lattice/", "results/lattice", -1, -1, 1)
bench("./graphs/undirected_dim/undirected_dim/ag/ag/", "results/ag", -1, -1, 1)
bench("./graphs/undirected_dim/undirected_dim/latin/latin/", "results/latin", -1, -1, 1)
bench("./graphs/undirected_dim/undirected_dim/latin-sw/latin-sw/", "results/latin-sw", -1, -1, 2)
bench("./graphs/undirected_dim/undirected_dim/sts/sts/", "results/sts", -1, -1, 2)
bench("./graphs/undirected_dim/undirected_dim/sts-sw/sts-sw/", "results/sts-sw", -1, -1, 2)
bench("./graphs/undirected_dim/undirected_dim/had/had/", "results/had", -1, 30, 1)
bench("./graphs/undirected_dim/undirected_dim/had-sw/had-sw/", "results/had-sw", -1, 30, 1)
bench("./graphs/undirected_dim/undirected_dim/cfi/cfi/", "results/cfi", -1, -1, 1)
bench("./graphs/ran2/ran2/", "results/ran2", -1, -1, 1)
bench("./graphs/ran10/ran10/", "results/ran10", -1, -1, 1)
bench("./graphs/ransq/ransq/", "results/ransq", -1, -1, 1)
bench("./graphs/dac/dac/other/", "results/dac_other", -1,-1, 2)
bench("./graphs/tran/tran/", "results/tran", -1, -1, 1)
#
bench("./graphs/large-cfi/", "results/large-cfi", -1, 30, 1)
bench("./graphs/rantree/rantree/", "results/rantree", -1, 120, 1)
bench("./graphs/ranreg/ranreg/", "results/ranreg", -1, 120, 1)
bench("./graphs/hypercubes/", "results/hypercubes", -1, -1, 1)
bench("./graphs/dac/dac/pipe/", "results/dac_pipe", -1, -1, 1)
bench("./graphs/combinatorial/combinatorial/", "results/combinatorial", 60, 60, 1)
bench("./graphs/undirected_dim/undirected_dim/pp/pp/pp16/", "results/pp16", -1, 60, 1)
bench("./graphs/undirected_dim/undirected_dim/pp/pp/pp25/", "results/pp25", 60, 60, 1)
