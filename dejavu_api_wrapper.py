from ctypes import *

dejavuc = CDLL("./libdejavu-api.so")   

def extract_paths(path_handle, n, return_coloring=False):
	paths = []
	path_num = dejavuc.path_get_num(path_handle)
	for path_id in range(0, path_num):
		path_constr = []
		path_size = dejavuc.path_get_size(path_handle, path_id)
		for pos in range(0, path_size):
			path_constr += [dejavuc.path_get_point(path_handle, path_id, pos)]
		if return_coloring:
			paths += [{'base_points' : path_constr, 'invariant' : dejavuc.path_get_inv(path_handle, path_id), 'coloring' : 				    extract_coloring(n, path_handle, path_id)}]
		else:
			paths += [{'base_points' : path_constr, 'invariant' : dejavuc.path_get_inv(path_handle, path_id)}] 
	return paths
	
def extract_coloring(n, path_handle, path_id):
	vertex_to_col = []
	for v in range(0, n):
		vertex_to_col += [dejavuc.path_get_vertex_color(path_handle, path_id, v)]
	return vertex_to_col

def random_ir_paths(n, edges, path_length, number_of_paths=1, fill_paths=False, vertex_labels=[], edge_labels=[], return_coloring=False):
	graph_handle = dejavuc.graph_create(n)
	print(graph_handle)
	if edge_labels == []:
		for [v1, v2] in edges:
			dejavuc.graph_add_edge(graph_handle, v1, v2)
	else:
		assert(len(edges) == len(edge_labels))
		for i in range(0, len(edges)):
			v1, v2 = edges[i]
			dejavuc.graph_add_edge_labelled(graph_handle, v1, v2, edge_labels[i])
		
	if vertex_labels == []:
		vertex_labels = [0] * n
	else:
		assert(n == len(vertex_labels))
	for v in range(0, n):
		dejavuc.graph_label(graph_handle, v, vertex_labels[v])
		
	path_handle = dejavuc.random_paths(graph_handle, path_length, number_of_paths, fill_paths)
	paths = extract_paths(path_handle, n, return_coloring)
	dejavuc.clean()
	return paths
	
def color_refinement(n, edges, vertex_labels=[], edge_labels=[]):
	path = random_ir_paths(n, edges, 0, number_of_paths=1, vertex_labels=vertex_labels, edge_labels=edge_labels, return_coloring = True)
	dejavuc.clean()
	return path[0]['coloring']
	
	
	
		
		
#print(random_ir_paths(5, [[0,1], [1,2], [2,3], [3,4], [4,0]], 10, edge_labels=[0,0,0,0,0]))

#print(color_refinement(5, [[0,1], [1,2], [2,3], [3,4], [4,0]], edge_labels=[0,0,0,0,0]))


#random_ir_paths(5, [[0,1], [1,2], [2,3], [3,4], [4,0]], 4)
