# -*- coding: utf-8 -*-
'''
Functions used to generated a FON2 networks
(State Network taking all subsequence of length <=2)
(Fixed order model taking all subsequence of length <=2)
'''

def fon2Rules(sequences):
	fon2 = {}
	for seq in sequences:
		if len(seq)>=2:
			for pos in range(1,len(seq)):
				
				i, j, k =  None, seq[pos-1], seq[pos]
				if pos > 1:
					i, j, k = seq[pos-2], seq[pos-1], seq[pos]
				
				## add order 1 rules (j,)
				if (j,) not in fon2.keys():
					fon2[(j,)] = {}
				if k not in fon2[(j,)].keys():
					fon2[(j,)][k] = 0.
				fon2[(j,)][k] += 1.
				
				if pos > 1:
					## add order 2 rules (i,j)
					if (i, j) not in fon2.keys():
						fon2[(i,j)] = {}
					if k not in fon2[(i,j)].keys():
						fon2[(i,j)][k] = 0.
					fon2[(i,j)][k] += 1.
	
	## normalize transition probability
	for r, nc in fon2.items():
		sum_count = sum(nc.values())
		for item, count_item in nc.items():
			fon2[r][item] = count_item / sum_count
	return fon2

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
			for pos in range(2,len(seq)):
				i, j, k = seq[pos-2], seq[pos-1], seq[pos]
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

