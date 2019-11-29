import subprocess
from os import walk

setname = "results/pp25"
# dir     = "/home/markus//Downloads/graphs/undirected_dim/undirected_dim/cfi/cfi/"
# dir     = "/home/markus//Downloads/graphs/undirected_dim/undirected_dim/rnd-3-reg/"
# dir     = "/home/markus//Downloads/graphs/ranreg/ranreg/"
# dir = "/home/markus/Downloads/graphs/undirected_dim/undirected_dim/lattice/lattice/"
# dir = "/home/markus/Downloads/graphs/undirected_dim/undirected_dim/grid/grid/"
#dir = "/home/markus/Downloads/graphs/dac/dac/pipe/"
# dir = "/home/markus/Downloads/graphs/undirected_dim/undirected_dim/latin/latin/"
dir = "/home/markus/Downloads/graphs/undirected_dim/undirected_dim/pp/pp/pp25/"
# dir = "/home/markus/Downloads/graphs/undirected_dim/undirected_dim/ag/ag/"
# dir = "/home/markus/Downloads/graphs/cfi-rigid-t2-tar/cfi-rigid-t2/larger/"
f = []
for (dirpath, dirnames, filenames) in walk(dir):
    f.extend(filenames)
    break

for i in range(len(f)):
        subprocess.run("./cmake-build-debug/dejavu --no_nauty --no_traces --THREADS_REFINEMENT_WORKERS 3 --file \"{}{}\" --stat_file {}.dat".format(dir, f[i], setname), shell=True)

