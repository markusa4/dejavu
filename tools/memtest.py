# Testing memory consumption using the python library of dejavu.
import libdejavu_python

for f in range(500000):
	print("GRAPH NUMBER", f)
	graph = libdejavu_python.pygraph()
	size = 50
	graph.set_size(size)
	for i in range(50):
		graph.add_edge(i, (i + 1) % 50)

	coloring = libdejavu_python.pycoloring()
	coloring.add_colors([0] * size) 

	max_length = 2
	num        = 1

	nodes = libdejavu_python._random_paths(graph, coloring, max_length, num)

	for node in nodes:
	  base = node.get_base_points()
	  coloring = node.get_vertex_to_col()

