# -*- coding: utf-8 -*-
'''
Main script to run experiments with the LFR benchmark
(Section IV, Fig. 3)
'''

import TestCasesGeneration

import InfoMapClust
import OverlappingNMI

from collections import defaultdict

########################################################################
def clusteringBySymbol(clust, isAgg = False, isDummy = False):
	flatten_clust = defaultdict(list)
	for rule, c in clust.items():
		if isAgg:
			last_symb = rule[0][-1]
			for r in rule:
				flatten_clust[last_symb].append(c)
		else:
			last_symb = rule[:rule.find('_')] if isDummy else rule[-1]
			flatten_clust[last_symb].append(c)
	return flatten_clust
########################################################################
def countClasses(clust):
	## clust possible overlap
	count_classes = defaultdict(int)
	sum_classes = 0
	for symb, classes in clust.items():
		sum_classes += len(classes)
		for c in classes:
			count_classes[c] = count_classes[c] + 1

	nb_classes     = len(count_classes.keys())
	nb_non_trivial = len([i for c,count in count_classes.items() if count > 1])
	return nb_classes, nb_non_trivial, sum_classes
########################################################################

###################################
## List of variable in the outputs
## (CSV data format output)
## The NMI similarity values used to make Fig. 3
## correspond to columns
## 10 (2-Von Diff Code)
## 14 (2-Von Same Code)
## 18 (Min 2-Von)
## 22 (2-Fon)
###################################
var_list = 'nb_nodes'
var_list+= ',om,mu,minc,maxc'
var_list+= ',nbc_base'

var_list+= ',nb_nodes_2o_ns,cl_2o_ns,nbc_2o_ns,nmi_2o_ns'
var_list+= ',nb_nodes_2o_unif,cl_2o_unif,nbc_2o_unif,nmi_2o_unif'
var_list+= ',nb_nodes_agg,cl_agg,nbc_agg,nmi_agg'
var_list+= ',nb_nodes_fon,cl_fon,nbc_fon,nmi_fon'
print(var_list)

#####################################
### List of LFR Benchmark parameters
### Entries are
### om   : ratio of overlapping nodes
### mu   : ratio of inter cluster edges
### minc : minimum size of clutsers
### maxc : maximum size of clutsers
#####################################
params = []
params.append(tuple([0.1,0.1,10,50]))
params.append(tuple([0.1,0.1,50,100]))
params.append(tuple([0.1,0.3,10,50]))
params.append(tuple([0.1,0.3,50,100]))
params.append(tuple([0.3,0.1,10,50]))
params.append(tuple([0.3,0.1,50,100]))
params.append(tuple([0.3,0.3,10,50]))
params.append(tuple([0.3,0.3,50,100]))

### Others (fixed) LFR Benchmark parameters (explain in Caption Fig. 3)
on, nb_nodes, k, maxk, t1, t2 = 2, 1000, 10, 50, 2., 1.

### Number of tests for each combinaison of parameters (should be 50)
nb_test = 50

### For each set of parameters
for p in params:
	om, mu, minc, maxc = p
	### We run 'nb test' tests
	for i in range(nb_test):
		## Generation of the 1st order network using the LFR Benchmark
		real_net, clustering = TestCasesGeneration.generateLFR(1, om, on, mu, nb_nodes,k,maxk,t1,t2,minc,maxc)

		#######################################
		## Statistics about the generated network
		comms = set()
		for n, clusts_n in clustering.items():
			for c in clusts_n:
				comms.add(c)
		nb_comm_lfr = len(comms)

		line_res = f'{nb_nodes}'
		line_res+= f',{om},{mu},{minc},{maxc}'
		line_res+= f',{nb_comm_lfr}'

		##########################################
		## Computation of Infomap
		## for the four different output
		## (see Section IV, p5)

		## VON2 using different code word for
		## representations of the same location in the same clusters
		hon2 = TestCasesGeneration.hon(real_net, clustering)
		hon2_node_weight  = InfoMapClust.uniformNodeWeights(hon2)
		hon2_codelength_ns, hon2_gain_ns, hon2_infomap_clust, time = InfoMapClust.infomapClustering(hon2, hon2_node_weight)
		hon2_clust_ns = clusteringBySymbol(hon2_infomap_clust)

		nb_nodes_hon2_ns = len(hon2_node_weight.keys())
		nb_classes_hon2_ns, nb_non_trivial_hon2_ns, sumc_hon2_ns = countClasses(hon2_clust_ns)
		nmi_2o_ns = OverlappingNMI.NMI_max(hon2_clust_ns, clustering)

		line_res+=f',{nb_nodes_hon2_ns},{hon2_codelength_ns},{nb_classes_hon2_ns},{round(nmi_2o_ns,4)}'

		## VON2 using the same code word for
		## representations of the same location in the same clusters
		nodes, states_hon2 = InfoMapClust.getStates(hon2)
		hon2_codelength, hon2_gain, hon2_clust, time = InfoMapClust.infomapStateClustering(nodes, states_hon2, hon2, hon2_node_weight)

		nb_nodes_hon2 = len(states_hon2.keys())
		nb_classes_hon2, nb_non_trivial_hon2, sumc_hon2 = countClasses(hon2_clust)
		nmi_2o = OverlappingNMI.NMI_max(hon2_clust, clustering)

		line_res+=f',{nb_nodes_hon2},{hon2_codelength},{nb_classes_hon2},{round(nmi_2o,4)}'

		## Minimal VON2 Network
		ideal_hon = TestCasesGeneration.idealHON(real_net, clustering)
		ideal_node_weight  = InfoMapClust.uniformNodeWeights(ideal_hon)
		nodes, states_ideal = InfoMapClust.getStates(ideal_hon)
		ideal_codelength, ideal_gain, ideal_clust, time = InfoMapClust.infomapStateClustering(nodes, states_ideal, ideal_hon, ideal_node_weight)

		nb_nodes_ideal = len(states_ideal.keys())
		nb_classes_ideal, nb_non_trivial_ideal, sumc_ideal = countClasses(ideal_clust)
		nmi_ideal = OverlappingNMI.NMI_max(ideal_clust, clustering)

		line_res+=f',{nb_nodes_ideal},{ideal_codelength},{nb_classes_ideal},{round(nmi_ideal,4)}'

		## FON2 Network
		nodes, states, states_network = TestCasesGeneration.stateNetwork(real_net, clustering)
		sta_node_weight  = InfoMapClust.uniformNodeWeights(states_network)
		state_codelength, state_gain, state_clust, time = InfoMapClust.infomapStateClustering(nodes, states, states_network, sta_node_weight)

		nb_states = len(states)
		nb_classes_state, nb_non_trivial_state, sumc_sta = countClasses(state_clust)
		nmi_state = OverlappingNMI.NMI_max(state_clust, clustering)

		line_res+=f',{nb_states},{state_codelength},{nb_classes_state},{round(nmi_state,4)}'

		## write output
		print(line_res)

