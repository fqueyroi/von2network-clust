# -*- coding: utf-8 -*-
'''
Functions to generate 1st order networks using command LFRBenchmark/benchmark
also contains functions to generate the ideal higher order network built from
the LFR generated network, clustering and using the flow dynamic described in the
paper.
'''

import subprocess
import re

def generateLFR(nb_seq=100, part_overlap = 0.01, om = 2, mu = 0.1, N=1000, k=10, kmax= 50, t1= 2, t2 = 1, minc = 10, maxc = 50):
	'''
	Generate a network using the LFR Benchmark with the given parameters.
	Outputs the network and the ground truth clustering

	WARNING!!!
	The function call the program './LFRBenchmark/benchmark' as a subprocess
	Halting the python script may not halt the subprocess
	'''
	on = int(N*part_overlap)
	flags = f'-N {N} -k {k} -maxk {kmax} -mu {mu} -t1 {t1} -t2 {t2} -minc {minc} -maxc {maxc} -on {on} -om {om}'

	cmd = f'./LFRBenchmark/benchmark {flags}'
	res = subprocess.run(cmd, shell = True, stdout=subprocess.PIPE)

	filename_net = "./network.dat"
	network = {}
	file_net = open(filename_net,'r')
	for line in file_net :
		split_line = line.strip().split('\t')
		src = split_line[0]
		tgt = split_line[1]
		if src not in network.keys():
			network[src] = []
		network[src].append(tgt)

	filename_clu = "./community.dat"
	clustering = {}
	file_clu = open(filename_clu,'r')
	communties = {}
	for line in file_clu :
		split_line = line.strip().split('\t')
		src = split_line[0]
		comms = split_line[1].strip().split(' ')
		clustering[src] = comms
		for c in comms:
			if c not in communties.keys():
				communties[c] = []
			communties[c].append(src)

	return network, clustering

def hon(network, clustering):
	'''
	Generate the ideal VON2 network
	'''
	hon={}
	index_clust = 0
	for i in network.keys():
		clusts_i = clustering[i]
		hon[tuple([i])] = {}
		if len(clusts_i)>1:
			## only add order 1 links
			for j in network[i]:
				hon[tuple([i])][tuple([j])] = 1
		else:
			c_i = clusts_i[0]
			for j in network[i]:
				clusts_j = clustering[j]
				ci_in_j = c_i in clusts_j
				if len(clusts_j)==1 or not ci_in_j:
					hon[tuple([i])][tuple([j])] = 1
				else:
					di = tuple([i,j])
					if di not in hon.keys():
						hon[di] = {}
					hon[tuple([i])][di] = 1
					for k in network[j]:
						clusts_k = clustering[k]
						if c_i in clusts_k:
							if len(clusts_k)==1:
								hon[di][tuple([k])] = 1

	return hon

def idealHON(network, clustering):
	'''
	Generate the Minimal VON2 network
	'''
	hon={}
	for i in network.keys():
		clusts_i = clustering[i]
		hon[tuple([i])] = {}
		if len(clusts_i)>1:
			## only add order 1 links
			for j in network[i]:
				hon[tuple([i])][tuple([j])] = 1
		else:
			c_i = clusts_i[0]
			for j in network[i]:
				clusts_j = clustering[j]
				ci_in_j = c_i in clusts_j
				if len(clusts_j)==1 or not ci_in_j:
					hon[tuple([i])][tuple([j])] = 1
				else:
					di = tuple([c_i,j])
					if di not in hon.keys():
						hon[di] = {}
					hon[tuple([i])][di] = 1
					for k in network[j]:
						clusts_k = clustering[k]
						if c_i in clusts_k:
							if len(clusts_k)==1:
								if tuple([k]) not in hon[di].keys():
									hon[di][tuple([k])] = 0.
								hon[di][tuple([k])] += 1.
	return(hon)

def stateNetwork(network, clustering):
	'''
	Generate the ideal FON2 network
	'''
	nodes = list(network.keys());
	states = {}
	states_network = {}
	for i in network.keys():
		clust_i = set(clustering[i])
		for j in network[i]:
			state_ij = str(i)+'_'+str(j)
			if state_ij not in states.keys():
				states[state_ij] = j
				states_network[state_ij] = {}

			clust_j = set(clustering[j])
			c_ij = clust_i & clust_j
			if len(clust_j)==1 or len(c_ij)==0:
				for k in network[j]:
					state_jk =  str(j)+'_'+str(k)
					if state_jk not in states.keys():
						states[state_jk] = k
						states_network[state_jk] = {}
					states_network[state_ij][state_jk] = 1.
			else:
				for k in network[j]:
					clust_k = set(clustering[k])
					c_ijk = c_ij & clust_k
					#### !!!!!!!!!
					if len(c_ijk)>0 and len(clust_k)==1:
					#if len(c_ijk)>0:
						state_jk =  str(j)+'_'+str(k)
						if state_jk not in states.keys():
							states[state_jk] = k
							states_network[state_jk] = {}
						states_network[state_ij][state_jk] = 1.
	return nodes, states, states_network

