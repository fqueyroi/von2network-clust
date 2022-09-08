# -*- coding: utf-8 -*-
'''
Functions used to run the python program 'infomap'
for a given network.
The script will write several temp files during its execution.
They will NOT be deleted.

Note: we did not used infomap python librairies since it procuded errors
(in version 1.3.0)
'''

import os
import subprocess
import re
from collections import defaultdict

################################################################################
def getStates(dupli_net, is_agg = False):
	'''
	Extract the representations
	as dict: representations -> locations
	(used for Infomap Clustering)
	'''
	nodes = []
	states = {}
	for u, neigh_u in dupli_net.items():
		node_u =  u[-1]
		if is_agg:
			node_u = u[0][-1]
		if node_u not in nodes:
			nodes.append(node_u)
		states[u] = node_u

		for v in neigh_u.keys():
			node_v =  v[-1]
			if is_agg:
				node_v = v[0][-1]
			if node_v not in nodes:
				nodes.append(node_v)
			states[v] = node_v
	return nodes, states

################################################################################
def uniformNodeWeights(network):
	node_weight = {}
	## loops also on target in order to
	## not miss some nodes
	for src in network.keys():
		node_weight[src] = 1.0
		for tgt in network[src].keys():
			node_weight[tgt] = 1.0
	return node_weight

################################################################################
def printGraph(graph, node_weight, filename, output_dir):
	'''
	Output the graph in output_dir in the Pajek .net format

	Parameters:
	--------
	graph: dict of (xxx -> dict of (xxx -> float))
	node_weight: dict of (xxx -> float)

	Returns:
	--------
	name_net_file: path of Pajek net file
	mapRule		 :  dict of (int -> xxx)
	'''
	mapIndex = {}
	mapRule  = {}

	name_net_file = output_dir + filename + '.net'
	# print(f'Printing net in {name_net_file}')

	file_net = open(name_net_file,'w')

	## Add vertices
	file_net.write(f'*Vertices {len(node_weight.keys())}\n')
	tmp_index = 1
	for src in node_weight.keys():
		mapIndex[src]      = tmp_index
		mapRule[tmp_index] = src
		file_net.write(f'{tmp_index} "{src}" {node_weight[src]}\n')
		tmp_index = tmp_index + 1

	## Add arcs
	file_net.write('*Arcs\n')
	for src in graph.keys():
		id_src = mapIndex[src]
		neigh_src = graph[src]
		for tgt, w in neigh_src.items():
			id_tgt = mapIndex[tgt]
			file_net.write(f'{id_src} {id_tgt} {w}\n')
	return name_net_file, mapRule

################################################################################
def runInfomap(network_file, output_dir, nb_loop = 10, is_state=False):
	'''
	Infomap options:
	----------------
	clu     				: "Write a clu file with the top cluster
							   ids for each node."
	d       				: "Assume directed links"
	to_nodes				: "Teleport to nodes instead of to links, assuming
							  uniform node weights if no such input data."
	N<n>		  		    : "Number of outer-most loops to run before
							   picking the best solution"
	'''
	options      = f'--clu -d --to-nodes -N {nb_loop}'
	# options      = f'--clu -d -N {nb_loop}'
	if is_state:
		options += ' -i states'
	res = subprocess.run(['infomap',options,network_file,output_dir], stdout=subprocess.PIPE)

	stdout = res.stdout.decode('utf-8')
	# print(stdout)
	match = re.search('Per level codelength total:.*\(sum: (\d+\.\d+)\)', stdout)
	if not match:
		print(f"Error in Infomap {stdout}")
		return False
	return True
################################################################################
def readClustering(mapRule, network_file, output_dir):
	'''
	Parameters:
	-----------
	mapRule: dict of (int -> xxx)

	Returns:
	--------
	codelength: float
	clust     : dict (node xxx -> clust id)
	'''
	name       = os.path.splitext(os.path.basename(network_file))
	clust_file = output_dir+name[0]+'.clu'
	# print(f"Reading {clust_file}")

	codelength, gain, time = 0., 0., 0.
	clust = {}
	file_seq = open(clust_file,'r')
	for line in file_seq :
		if line[0] == '#':
			match = re.search('# codelength (\d+\.\d+) bits', line)
			if match:
				codelength = float(match.group(1))
			match = re.search('# completed in (\d+\.\d+) s', line)
			if match:
				time = float(match.group(1))
			match = re.search('# relative codelength savings (\d+\.\d+)%', line)
			if match:
				gain = float(match.group(1))

		else :
			split_line = line.strip().split(' ')
			n = int(split_line[0])
			c = int(split_line[1])
			clust[mapRule[n]] = c
	return codelength, gain, clust, time

################################################################################
def infomapClustering(graph, node_weight, file_path = '', nb_loop = 10):
	'''
	Main method

	Parameters:
	--------
	graph: dict of (xxx -> dict of (xxx -> float))
	node_weight: dict of (xxx -> float)

	Returns:
	--------
	codelength: float
	clust: dict (node id -> clust id)
	'''
	output_dir = './output_infomap/'
	filename = 'out_infomap'
	if file_path != '':
		dir_path = os.path.dirname(file_path)
		filename = os.path.splitext(os.path.basename(file_path))[0]
		output_dir = dir_path + "/output_infomap/"
	os.makedirs(output_dir, exist_ok=True)
	net_filepath, mapRule = printGraph(graph, node_weight, filename, output_dir)
	run_sucess = runInfomap(net_filepath, output_dir, nb_loop)
	codelength, gain, clust, time = readClustering(mapRule, net_filepath, output_dir)
	#os.removedirs(output_dir)
	return codelength, gain, clust, time / nb_loop

################################################################################
def infomapStateClustering(nodes, states, states_network, weight, file_path = '', nb_loop = 10):
	'''
	Main method (for Network with states)

	Parameters:
	--------
	graph: dict of (xxx -> dict of (xxx -> float))
	node_weight: dict of (xxx -> float)

	Returns:
	--------
	codelength: float
	clust: dict (node id -> clust id)
	'''
	output_dir = './output_infomap/'
	filename = 'out_infomap'
	if file_path != '':
		dir_path = os.path.dirname(file_path)
		filename = os.path.splitext(os.path.basename(file_path))[0]
		output_dir = dir_path + "/output_infomap/"
	os.makedirs(output_dir, exist_ok=True)
	net_filepath, mapRule = printStateGraph(nodes, states, states_network, weight, filename, output_dir)
	run_sucess = runInfomap(net_filepath, output_dir, nb_loop, True)
	codelength, gain, clust, time = readStateClustering(mapRule, net_filepath, output_dir)
	# os.removedirs(output_dir)
	return codelength, gain, clust, time / nb_loop

################################################################################
def printStateGraph(nodes, states, states_network, weight, filename, output_dir):
	mapIndex = {}
	mapRule  = {}

	name_net_file = output_dir + filename + '.net'
	# print(f'Printing net in {name_net_file}')

	file_net = open(name_net_file,'w')

	## Add vertices
	file_net.write(f'*Vertices {len(nodes)}\n')
	tmp_index = 1
	for src in nodes:
		mapIndex[src]      = tmp_index
		mapRule[tmp_index] = src
		file_net.write(f'{tmp_index} "{src}"\n')
		tmp_index += 1

	## Add states
	mapIndexState = {}
	file_net.write(f'*States\n')
	tmp_index = 1
	for state, phys_node in states.items():
		phys_id = mapIndex[phys_node]
		mapIndexState[state] = tmp_index
		file_net.write(f'{tmp_index} {phys_id} "{state}" {weight[state]}\n')
		tmp_index += 1

	## Add arcs
	file_net.write('*Arcs\n')
	for src in states_network.keys():
		id_src = mapIndexState[src]
		neigh_src = states_network[src]
		for tgt, w in neigh_src.items():
			id_tgt = mapIndexState[tgt]
			file_net.write(f'{id_src} {id_tgt} {w}\n')
	file_net.close()

	return name_net_file, mapRule
################################################################################
def readStateClustering(mapRule, network_file, output_dir):
	name       = os.path.splitext(os.path.basename(network_file))
	clust_file = output_dir+name[0]+'.clu'
	# print(f"Reading {clust_file}")

	codelength, gain, time = 0., 0., 0.
	clust = {}
	file_seq = open(clust_file,'r')
	for line in file_seq :
		if line[0] == '#':
			match = re.search('# codelength (\d+\.\d+) bits', line)
			if match:
				codelength = float(match.group(1))
			match = re.search('# completed in (\d+\.\d+) s', line)
			if match:
				time = float(match.group(1))
			match = re.search('# relative codelength savings (\d+\.\d+)%', line)
			if match:
				gain = float(match.group(1))
		else :
			split_line = line.strip().split(' ')
			n = int(split_line[0])
			c = int(split_line[1])
			if mapRule[n] not in clust.keys():
				clust[mapRule[n]] = []
			clust[mapRule[n]].append(c)
	return codelength, gain, clust, time

