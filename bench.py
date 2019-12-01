import subprocess
from os import walk

timeout = 30
setname = "results/rantree"
dir     = "/home/markus//Downloads/graphs/rantree/rantree/"
# dir     = "/home/markus//Downloads/graphs/undirected_dim/undirected_dim/rnd-3-reg/"
# dir     = "/home/markus//Downloads/graphs/ranreg/ranreg/"
# dir = "/home/markus/Downloads/graphs/undirected_dim/undirected_dim/sts-sw/sts-sw/"
# dir = "/home/markus/Downloads/graphs/undirected_dim/undirected_dim/grid/grid/"
#dir = "/home/markus/Downloads/graphs/dac/dac/pipe/"
# dir = "/home/markus/Downloads/graphs/undirected_dim/undirected_dim/pp/pp/pp25/"
# dir = "/home/markus/Downloads/graphs/benchmarks/benchmarks/srg/"
# dir = "/home/markus/Downloads/graphs/undirected_dim/undirected_dim/ag/ag/"
# dir = "/home/markus/Downloads/graphs/combinatorial/combinatorial/"
# dir = "/home/markus/Downloads/graphs/cfi-rigid-t2-tar/cfi-rigid-t2/larger/"
f = []
for (dirpath, dirnames, filenames) in walk(dir):
    f.extend(filenames)
    break

for i in range(len(f)):
        subprocess.run("./cmake-build-debug/dejavu --timeout {} --threads 6 --no_nauty --no_traces --file \"{}{}\" --stat_file {}.dejavu7.dat".format(timeout, dir, f[i], setname), shell=True)

for i in range(len(f)):
        subprocess.run("./cmake-build-debug/dejavu --timeout {} --threads 3 --no_nauty --no_traces --file \"{}{}\" --stat_file {}.dejavu4.dat".format(timeout, dir, f[i], setname), shell=True)

for i in range(len(f)):
        subprocess.run("./cmake-build-debug/dejavu --timeout {} --threads 0 --force_selector 1 --no_nauty --no_traces --file \"{}{}\" --stat_file {}.dejavu1.dat".format(timeout, dir, f[i], setname), shell=True)


for i in range(len(f)):
        subprocess.run("./cmake-build-debug/dejavu --timeout {} --no_nauty --no_dejavu --file \"{}{}\" --stat_file {}.traces.dat".format(timeout, dir, f[i], setname), shell=True)


for i in range(len(f)):
        subprocess.run("./cmake-build-debug/dejavu --timeout {} --no_traces --no_dejavu --file \"{}{}\" --stat_file {}.nauty.dat".format(timeout, dir, f[i], setname), shell=True)
