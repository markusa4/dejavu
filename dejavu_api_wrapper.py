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

"""
    Compute random paths starting from the root of the IR tree. The graph with 'n' vertices is given as an edgelist in 'edges'. Paths are limited to length 'path_length'. Can also return multiple paths at once using 'number_of_paths'. A vertex labelling can be given using 'vertex_labels', which means that vertex 'i' is labelled with 'vertex_labels[i]'. Similarly, an edge labelling can be given using 'edge_labels'. However, note that edge labels are internally resolved using edge subdivision and vertex labels, reducing performance. 
    
    The method returns a dictionary containing the base points of the computed paths and an invariant. If the coloring of the node where a path ends is required, set 'return_coloring'.    
"""
def random_ir_paths(n, edges, path_length, number_of_paths=1, fill_paths=False, vertex_labels=[], 
                    edge_labels=[], return_coloring=False, directed_dimacs=False):
	graph_handle = __graph_to_handle(n, edges, vertex_labels, edge_labels, directed_dimacs=directed_dimacs)
	path_handle = dejavuc.random_paths(ctypes.c_int(graph_handle), ctypes.c_int(path_length), 
					    ctypes.c_int(number_of_paths), ctypes.c_bool(fill_paths))
	paths = __extract_paths(path_handle, n, return_coloring)
	dejavuc.clean()
	return paths

"""
    Computes color refinement on a given graph. The graph with 'n' vertices is given as an edgelist in 'edges'. A vertex labelling can be given using 'vertex_labels', which means that vertex 'i' is labelled with 'vertex_labels[i]'. Similarly, a edge labelling can be given using 'edge_labels'. However, note that edge labels are internally resolved using edge subdivision and vertex labels, reducing performance. 
"""
def color_refinement(n, edges, vertex_labels=[], edge_labels=[], directed_dimacs=False):
	path = random_ir_paths(n, edges, 0, number_of_paths=1, vertex_labels=vertex_labels, 
	                       edge_labels=edge_labels, return_coloring = True, directed_dimacs=directed_dimacs)
	dejavuc.clean()
	return path[0]['coloring']

"""
    Probabilistic isomorphism test for two given graphs. The first graph on 'n1' vertices is given as an edgelist in 'edges1'. The second graph is given using 'n2' and 'edges2'. The optional parameter 'err' determines the error probability. If the two graphs are non-isomorphic, the method is guaranteed to determine this. If they are isomorphic, the method returns non-isomorphic with probability at most 1 / 2^err.
"""
def are_isomorphic(n1, edges1, n2, edges2, err=8):
	graph_handle1 = __graph_to_handle(n1, edges1, [], [])
	graph_handle2 = __graph_to_handle(n2, edges2, [], [])
	is_iso = dejavuc.are_isomorphic(ctypes.c_int(graph_handle1), ctypes.c_int(graph_handle2), ctypes.c_int(err))
	dejavuc.clean()
	return is_iso

"""
    Computes a generating set for the automorphism group of the given graph probabilistically. The graph on 'n' vertices is given as an edgelist in 'edges'. A vertex labelling can be given using 'vertex_labels'. The optional parameter 'err' determines the error probability. The method will never return generators that are not symmetries of the graph (one-sided error). The method is guaranteed to return all generators with probability at least 1 / 2^err.
"""
def get_automorphisms(n, edges, err=8, vertex_labels=[], directed_dimacs=False):
	graph_handle = __graph_to_handle(n, edges, vertex_labels, [], directed_dimacs=directed_dimacs)
	auto_handle = dejavuc.get_automorphisms(ctypes.c_int(graph_handle), ctypes.c_int(err))
	grp_info = __extract_generators(auto_handle, n)
	dejavuc.clean()
	return grp_info

