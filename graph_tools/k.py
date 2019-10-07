n = 150
f = open("k%d.dimacs" %n,"w+")

f.write("p edge {} {}\n".format(n, int((n * (n - 1)) / 2)))

for i in range(1, n + 1):
    for j in range(i + 1, n + 1):
     f.write("e {} {}\n".format(i, j))
