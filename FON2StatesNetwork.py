# -*- coding: utf-8 -*-
'''
Functions used to generated a FON2 networks
(Fixed order model taking all subsequence of length 2)
'''

def order2Rules(sequences):
	order2 = {}
	for seq in sequences:
		if len(seq)>=3:
			for i in range(2,len(seq)):
				i, j, k = seq[i-2], seq[i-1], seq[i]
				if i not in order2.keys():
					order2[i] = {}
				o2i = order2[i]
				if j not in o2i.keys():
					o2i[j] = {}
				o2ij = o2i[j]
				if k not in o2ij.keys():
					o2ij[k] = 0.
				o2ij[k] += 1.
	for i, o2i in order2.items():
		for j, o2ij in o2i.items():
			sum_ij = sum(o2ij.values())
			for k in list(o2ij.keys()):
				o2ij[k] = o2ij[k] / sum_ij
	return order2

def buildFON2Network(sequences):
	nodes = []
	states = {}
	states_network = {}

	for seq in sequences:
		for s in seq:
			if s not in nodes:
				nodes.append(s)
		if len(seq)>=3:
			for i in range(2,len(seq)):
				i, j, k = seq[i-2], seq[i-1], seq[i]
				state_ij = str(i)+'-'+str(j)
				state_jk = str(j)+'-'+str(k)
				states[state_ij] = j
				states[state_jk] = k
				if state_ij not in states_network.keys():
					states_network[state_ij] = {}
				if state_jk not in states_network.keys():
					states_network[state_jk] = {}
				if state_jk not in states_network[state_ij].keys():
					states_network[state_ij][state_jk] = 0.
				states_network[state_ij][state_jk] += 1.
	return nodes, states, states_network

