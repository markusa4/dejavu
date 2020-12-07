import libdejavu_python

for f in range(500000):
	print("GRAPH NUMBER", f)
	graph = libdejavu_python.pygraph()
	size = 50
	graph.set_size(size)
	for i in range(50):
		graph.add_edge(i, (i + 1) % 50)
	#graph.set_size(4)
	#graph.add_edges([0, 1, 2],
	#	        [1, 2, 3])
	# graph.set_directed_dimacs(1) # directed edge definition of DIMACS is used (CAUTION: directed graphs are not supported)

	coloring = libdejavu_python.pycoloring()
	coloring.add_colors([0] * size) 

	max_length = 2
	num        = 1

	nodes = libdejavu_python._random_paths(graph, coloring, max_length, num)
	print("Got", len(nodes), "paths.")

	for node in nodes:
	  print("Node", node)
	  print("Path invariant", node.get_invariant())
	  base = node.get_base_points()
	  print("Path length", len(base))
	  print("Path nodes", base)
	  coloring = node.get_vertex_to_col()
	  print("Coloring at the end of path", coloring)

