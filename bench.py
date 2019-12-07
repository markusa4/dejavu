import subprocess
from datetime import datetime
from random import seed
from random import randint
from os import walk

def bench(dir, setname, timeout, selector):
	f     = []
	seeds = []
	seed(datetime.now())
	for (dirpath, dirnames, filenames) in walk(dir):
	    f.extend(filenames)
	    break
	for i in range(len(f)):
		seeds += [randint(0, 100000)]
	for i in range(len(f)):
		subprocess.run("./cmake-build-debug/dejavu --permute_seed {} --timeout {} --threads 7 --no_nauty --no_traces --file \"{}{}\" --stat_file {}.dejavu8.dat".format(seeds[i], timeout, dir, f[i], setname), shell=True)
	for i in range(len(f)):
		subprocess.run("./cmake-build-debug/dejavu --permute_seed {} --timeout {} --threads 6 --no_nauty --no_traces --file \"{}{}\" --stat_file {}.dejavu7.dat".format(seeds[i], timeout, dir, f[i], setname), shell=True)
	for i in range(len(f)):
		subprocess.run("./cmake-build-debug/dejavu --permute_seed {} --timeout {} --threads 3 --no_nauty --no_traces --file \"{}{}\" --stat_file {}.dejavu4.dat".format(seeds[i], timeout, dir, f[i], setname), shell=True)
	for i in range(len(f)):
		subprocess.run("./cmake-build-debug/dejavu --permute_seed {} --timeout {} --threads 1 --force_selector {} --no_nauty --no_traces --file \"{}{}\" --stat_file {}.dejavu2.dat".format(seeds[i], timeout, selector, dir, f[i], setname), shell=True)
	for i in range(len(f)):
		subprocess.run("./cmake-build-debug/dejavu --permute_seed {} --timeout {} --threads 0 --force_selector {} --no_nauty --no_traces --file \"{}{}\" --stat_file {}.dejavu1.dat".format(seeds[i], timeout, selector, dir, f[i], setname), shell=True)
	for i in range(len(f)):
		subprocess.run("./cmake-build-debug/dejavu --permute_seed {} --timeout {} --no_nauty --no_dejavu --file \"{}{}\" --stat_file {}.traces.dat".format(seeds[i], timeout, dir, f[i], setname), shell=True)
	for i in range(len(f)):
		subprocess.run("./cmake-build-debug/dejavu --permute_seed {} --timeout {} --no_traces --no_dejavu --file \"{}{}\" --stat_file {}.nauty.dat".format(seeds[i], timeout, dir, f[i], setname), shell=True)


bench("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/grid/grid/", "results/grid", 30)
bench("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/lattice/lattice/", "results/lattice", 30)
bench("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/ag/ag/", "results/ag", 30, 1)
bench("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/latin/latin/", "results/latin", 30, 1)
bench("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/latin-sw/latin-sw/", "results/latin-sw", 30, 2)
bench("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/sts/sts/", "results/sts", 30, 2)
bench("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/sts-sw/sts-sw/", "results/sts-sw", 30, 2)
bench("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/had/had/", "results/had", 30, 1)
bench("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/had-sw/had-sw/", "results/had-sw", 30, 1)
bench("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/cfi/cfi/", "results/cfi", 30, 0)
bench("/home/markus/Downloads/graphs/rantree/rantree/", "results/rantree", 30, 1)
bench("/home/markus/Downloads/graphs/rantree/ran2/", "results/ran2", 30, 1)
bench("/home/markus/Downloads/graphs/rantree/ran10/", "results/ran10", 30, 1)
bench("/home/markus/Downloads/graphs/rantree/ransq/", "results/ransq", 30, 1)
bench("/home/markus/Downloads/graphs/ranreg/ranreg/", "results/ranreg", 30, 1)
bench("/home/markus/Downloads/graphs/hypercubes/", "results/hypercubes", 30, 1)
bench("/home/markus/Downloads/graphs/dac/dac/pipe/", "results/dac_pipe", 30, 1)
bench("/home/markus/Downloads/graphs/dac/dac/other/", "results/dac_other", 30, 2)
bench("/home/markus/Downloads/graphs/tran/tran/", "results/tran", 30, 1)
bench("/home/markus/Downloads/graphs/combinatorial/combinatorial/", "results/combinatorial", 60, 2)
bench("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/pp/pp/pp16/", "results/pp16", 60, 1)
bench("/home/markus/Downloads/graphs/undirected_dim/undirected_dim/pp/pp/pp25/", "results/pp25", 60, 1)

