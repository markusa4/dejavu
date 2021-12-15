import ctypes
import ctypes.util
import importlib.machinery as impm
import os

mod_path = os.path.dirname(__file__)
for suff in impm.EXTENSION_SUFFIXES:
	dejavuc_path = str(mod_path) + "/libdejavu_api" + str(suff)
	if os.path.isfile(dejavuc_path):
		break

if os.name == 'nt':
	dejavuc = ctypes.CDLL(dejavuc_path, winmode=0)
else:
	dejavuc = ctypes.CDLL(dejavuc_path)

#dejavuc = ctypes.CDLL("./libdejavu-api.so")   

def __extract_paths(path_handle, n, return_coloring=False):
	paths = []
	path_num = dejavuc.path_get_num(ctypes.c_int(path_handle))
	for path_id in range(0, path_num):
		path_constr = []
		path_size = dejavuc.path_get_size(path_handle, path_id)
		for pos in range(0, path_size):
			path_constr += [dejavuc.path_get_point(path_handle, path_id, pos)]
		if return_coloring:
			paths += [{'base_points' : path_constr, 'invariant' : dejavuc.path_get_inv(path_handle, path_id), 
				   'coloring' : __extract_coloring(n, path_handle, path_id)}]
		else:
			paths += [{'base_points' : path_constr, 'invariant' : dejavuc.path_get_inv(path_handle, path_id)}]
	return paths
	
def __extract_generators(path_handle, n):
	gens = []
	path_num = dejavuc.path_get_num(ctypes.c_int(path_handle))
	for path_id in range(0, path_num):
		gens += [__extract_coloring(n, path_handle, path_id)]
	
	base_sz = dejavuc.path_get_base_size(ctypes.c_int(path_handle))
	base = []
	for i in range(0, base_sz):
		base += [dejavuc.path_get_base_point(ctypes.c_int(path_handle), ctypes.c_int(i))]
	
	dejavuc.path_get_grpsz1.restype = ctypes.c_double
	grp_sz1 = dejavuc.path_get_grpsz1(ctypes.c_int(path_handle))
	grp_sz2 = dejavuc.path_get_grpsz2(ctypes.c_int(path_handle))
	grp_sz = grp_sz1 * pow(10, grp_sz2) 
	
	return {'generators' : gens, 'base' : base, 'size' : grp_sz}  
	
def __extract_coloring(n, path_handle, path_id):
	vertex_to_col = []
	for v in range(0, n):
		vertex_to_col += [dejavuc.path_get_vertex_color(ctypes.c_int(path_handle), 
								  ctypes.c_int(path_id), ctypes.c_int(v))]
	return vertex_to_col

def __graph_to_handle(n, edges, vertex_labels, edge_labels, directed_dimacs=False):
	graph_handle = dejavuc.graph_create(ctypes.c_int(n))
	dejavuc.graph_set_directed_dimacs(ctypes.c_int(graph_handle), ctypes.c_bool(directed_dimacs))
	if edge_labels == []:
		for [v1, v2] in edges:
			dejavuc.graph_add_edge(ctypes.c_int(graph_handle), ctypes.c_int(v1), ctypes.c_int(v2))
	else:
		assert(len(edges) == len(edge_labels))
		for i in range(0, len(edges)):
			v1, v2 = edges[i]
			assert(edge_labels[i] < 1073741823 and edge_labels[i] >= 0)
			assert(v1 < n and v1 >= 0)
			assert(v2 < n and v1 >= 0)
			dejavuc.graph_add_edge_labelled(ctypes.c_int(graph_handle), 
							 ctypes.c_int(v1), ctypes.c_int(v2), ctypes.c_int(edge_labels[i]))
		
	if vertex_labels == []:
		vertex_labels = [0] * n
	else:
		assert(n == len(vertex_labels))
	for v in range(0, n):
		assert(vertex_labels[v] < 1073741823 and vertex_labels[v] >= 0)
		dejavuc.graph_label(ctypes.c_int(graph_handle), ctypes.c_int(v), ctypes.c_int(vertex_labels[v]))
	return graph_handle;

def random_ir_paths(n, edges, path_length, number_of_paths=1, fill_paths=False, vertex_labels=[], 
                    edge_labels=[], return_coloring=False, directed_dimacs=False):
	graph_handle = __graph_to_handle(n, edges, vertex_labels, edge_labels, directed_dimacs=directed_dimacs)
	path_handle = dejavuc.random_paths(ctypes.c_int(graph_handle), ctypes.c_int(path_length), 
					    ctypes.c_int(number_of_paths), ctypes.c_bool(fill_paths))
	paths = __extract_paths(path_handle, n, return_coloring)
	dejavuc.clean()
	return paths
	
def color_refinement(n, edges, vertex_labels=[], edge_labels=[], directed_dimacs=False):
	path = random_ir_paths(n, edges, 0, number_of_paths=1, vertex_labels=vertex_labels, 
	                       edge_labels=edge_labels, return_coloring = True, directed_dimacs=directed_dimacs)
	dejavuc.clean()
	return path[0]['coloring']
	
def are_isomorphic(n1, edges1, n2, edges2, err=8):
	graph_handle1 = __graph_to_handle(n1, edges1, [], [])
	graph_handle2 = __graph_to_handle(n2, edges2, [], [])
	is_iso = dejavuc.are_isomorphic(ctypes.c_int(graph_handle1), ctypes.c_int(graph_handle2), ctypes.c_int(err))
	dejavuc.clean()
	return is_iso

def get_automorphisms(n, edges, err=8, vertex_labels=[], directed_dimacs=False):
	graph_handle = __graph_to_handle(n, edges, vertex_labels, [], directed_dimacs=directed_dimacs)
	auto_handle = dejavuc.get_automorphisms(ctypes.c_int(graph_handle), ctypes.c_int(err))
	grp_info = __extract_generators(auto_handle, n)
	dejavuc.clean()
	return grp_info

#cycle1 = [[0,1], [1,2], [2,3], [4,3], [4,0]]
#cycle2 = [[1,0], [2,3], [3,4], [1,2], [4,0]]



#print(are_isomorphic(5, cycle1, 5, cycle2))

def _fix_digraph(edges):
	n_edges = []
	for e in edges:
		v1, v2 = e
		if not [v2, v1] in n_edges:
			n_edges += [e]
	return n_edges

#edges = [[0, 23], [0, 27], [0, 28], [1, 2], [1, 3], [1, 4], [1, 30], [2, 1], [2, 3], [2, 4], [2, 5], [3, 1], [3, 2], [3, 4], [3, 5], [3, 32], [4, 1], [4, 2], [4, 3], [4, 5], [5, 2], [5, 3], [5, 4], [5, 8], [5, 32], [6, 7], [6, 8], [6, 31], [6, 32], [7, 6], [7, 8], [7, 31], [8, 5], [8, 6], [8, 7], [8, 26], [9, 10], [9, 24], [9, 25], [9, 26], [10, 9], [10, 11], [10, 24], [11, 10], [11, 12], [11, 33], [12, 11], [12, 33], [12, 36], [13, 14], [13, 33], [13, 36], [14, 13], [14, 15], [14, 36], [15, 14], [15, 16], [15, 38], [16, 15], [16, 17], [16, 38], [17, 16], [17, 38], [17, 39], [18, 19], [18, 35], [18, 39], [19, 18], [19, 20], [19, 35], [20, 19], [20, 21], [20, 35], [21, 20], [21, 22], [21, 34], [21, 37], [22, 21], [22, 34], [22, 37], [23, 0], [23, 27], [23, 28], [24, 9], [24, 10], [24, 25], [24, 26], [25, 9], [25, 24], [25, 26], [26, 8], [26, 9], [26, 24], [26, 25], [27, 0], [27, 23], [27, 28], [27, 29], [27, 30], [28, 0], [28, 23], [28, 27], [28, 29], [28, 30], [29, 27], [29, 28], [29, 30], [30, 1], [30, 27], [30, 28], [30, 29], [31, 6], [31, 7], [31, 32], [32,3], [32, 5], [32, 6], [32, 31], [33, 11], [33, 12], [33, 13], [33, 36], [34, 21], [34, 22], [34, 37], [35, 18], [35, 19], [35, 20], [35, 39], [36, 12], [36, 13], [36, 14], [36, 33], [37, 21], [37, 22], [37, 34], [38, 15], [38, 16], [38, 17], [38, 39], [39, 17], [39, 18], [39, 35], [39, 38]]

#vlabels = [1628273134, -2132368460, 652329178, -387690958, -1435709788, 652329178, -387690958, 652329178, 1628273134, 652329178, 26827912, 1254396766, 652329178, 1628273134, -387690958, 136319548, 1097698410, -1327612506, 1628273134, -1809771368, -1179349854, -639656304, 652329178, 2146176196, 2146176196, -1639394214, -1065704762, -250097854, -250097854, -1065704762, -1639394214, 2146176196, 2146176196, 2146176196, 2146176196, -1065704762, -1065704762, -1639394214, -1065704762, -1065704762]

#vlabels_fix = list(map(lambda i: abs(i % 1073741823), vlabels))

#vlabels_fix2x = vlabels_fix + vlabels_fix

#elabels = [35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232,35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232, 35235232]

#print(get_automorphisms(5, cycle1))

#elabels_fix = [35235232]*len(fix_digraph(edges))

#for i in range(0, 2000000):
#	print(get_automorphisms(5, cycle1))
#	are_isomorphic(5, cycle1, 5, cycle2)
#	graph_to_handle(40, fix_digraph(edges), vlabels, elabels_fix)
#	graph_to_handle(40, fix_digraph(edges), vlabels, [])
#	dejavuc.clean()
#	print(random_ir_paths(40, edges, 5, vertex_labels=vlabels_fix, number_of_paths=3, edge_labels=elabels, fill_paths=True, directed_dimacs=True))
#	print(random_ir_paths(80, edges, 5, vertex_labels=vlabels_fix2x, number_of_paths=3, edge_labels=elabels, directed_dimacs=True))
	#print(random_ir_paths(5, [[0,1], [1,2], [2,3], [3,4], [4,0]], 10, edge_labels=[0,0,0,0,0], number_of_paths=3,return_coloring=True))
	#print(random_ir_paths(5, [[0,1], [1,2], [2,3], [3,4], [4,0], [3,0]], 10, number_of_paths=3,return_coloring=True))
	#print(random_ir_paths(5, [[0,1], [1,2], [2,3], [3,4], [4,0]], 10, edge_labels=[0,0,0,0,1], number_of_paths=3,return_coloring=True))
	#print(random_ir_paths(10, [[0,1], [1,2], [2,3], [3,4], [4,0]], 10, edge_labels=[0,0,0,0,0], number_of_paths=3,return_coloring=True))
	#print(color_refinement(5, [[0,1], [1,2], [2,3], [3,4], [4,0]], edge_labels=[0,0,0,0,1]))

#print(color_refinement(5, [[0,1], [1,2], [2,3], [3,4], [4,0]], edge_labels=[0,0,0,0,0]))


#random_ir_paths(5, [[0,1], [1,2], [2,3], [3,4], [4,0]], 4)
