# -*- coding: utf-8 -*-
'''
Main script used to evaluate the clusterings of
VON2, Aggregated Von2 and FON2 Networks
'''

import HONUtils
import BuildRulesFast
import BuildNetwork

import AggOrder2Rules
import FON2StatesNetwork

import InfoMapClust
import OverlappingNMI

import time
from collections import defaultdict

##########################
## SEQUENCES FILE INPUT ##
##########################

## Computing results on the Maritime Dataset
filename = "./maritime_sequences.csv"
sep = " " ## the string separating elements in a sequence

## Computing results on the Airports Dataset
# filename = "./2011Q1_SEQ.csv"
# sep = ","

## Computing results on the Taxis Dataset
#filename = "./trajectories_PoliceStation.csv"
#sep = " "

save_cluster_file = "./clusters.csv"
save_stats_file   = "./clusters_stats.csv"

################
## PARAMETERS ##
################

## Number of execution of Infomap for each network (eg 50)
nb_loop = 50

########################
## VARIOUS FUNCTIONS  ##
########################
################################################################################
def getNbDuplication(states):
	nb_dupli = defaultdict(int)
	for d, n  in states.items():
		nb_dupli[n] = nb_dupli[n] + 1
	return nb_dupli
################################################################################
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
################################################################################
def countClusters(clust):
	clust_set = set()
	clust_count = defaultdict(int)
	symb_count_clust = {}
	for symb, id_clusts in clust.items():
		clust_set.update(id_clusts)
		symb_count_clust[symb] = len(set(id_clusts))
		for c in set(id_clusts):
			clust_count[c] = clust_count[c] + 1
	return len(clust_set), symb_count_clust, clust_count
################################################################################

################################################
## BUILDING 2nd Order NETWORK and Aggregation ##
################################################

## Read trajectories file
sequences = HONUtils.readSequenceFile(filename,True,sep)
sequences = HONUtils.removeRepetitions(sequences)

## Build rules
start_time = time.time()
rule_builder = BuildRulesFast.FastHONRulesBuilder(sequences,2,1,1.)
rules = rule_builder.ExtractRules()
time_build_rules =  time.time() - start_time

################################################
## INFOMAP CLUSTERING FOR EACH NETWORK ##
################################################

### List to keep output values
final_nb_dupli = []
final_clusts = []
final_codelengths = []
final_gain        = []
final_clusts_name = []
final_clusts_times = []
final_build_times  = []

print('#################################')
print(f'Computing clustering for VON2 ({nb_loop} infomap runs)')
## VON2 Network clustering
start_time = time.time()
network = BuildNetwork.BuildNetwork(rules)
time_build_von2 = time_build_rules + (time.time() - start_time)
von2_node_weight  = InfoMapClust.uniformNodeWeights(network)

nodes, states_von2 = InfoMapClust.getStates(network)
von2_codelength, von2_gain, von2_clust, time_2o = InfoMapClust.infomapStateClustering(nodes, states_von2, network, von2_node_weight, filename, nb_loop)

final_clusts.append(von2_clust)
final_codelengths.append(von2_codelength)
final_gain.append(von2_gain)
final_clusts_name.append("Von2")
final_clusts_times.append(time_2o)
final_nb_dupli.append(getNbDuplication(states_von2))
final_build_times.append(time_build_von2)
print('Done.')

## Agg Von-2
print('#################################')
print(f'Computing clustering for Aggregated VON2 ({nb_loop} infomap runs)')
start_time = time.time()
clusts = AggOrder2Rules.aggregateRules(rules)
flatten_agg_rules = AggOrder2Rules.flattenAgg2ndOrderRules(rules, clusts)
agg_network       = AggOrder2Rules.mergeNodes(network, flatten_agg_rules)
time_build_agg = time_build_rules + (time.time() - start_time)
agg_node_weight   = InfoMapClust.uniformNodeWeights(agg_network)

nodes, states_agg = InfoMapClust.getStates(agg_network, True)
agg_codelength, agg_gain, agg_clust, time_agg = InfoMapClust.infomapStateClustering(nodes, states_agg, agg_network, agg_node_weight, filename, nb_loop)

final_clusts.append(agg_clust)
final_codelengths.append(agg_codelength)
final_gain.append(agg_gain)
final_clusts_name.append("Agg Von2")
final_clusts_times.append(time_agg)
final_nb_dupli.append(getNbDuplication(states_agg))
final_build_times.append(time_build_agg)
print('Done.')

## FON2 Network
print('#################################')
print(f'Computing clustering for FON2 ({nb_loop} infomap runs)')
start_time = time.time()
nodes, fon2_states, fon2_state_net = FON2StatesNetwork.buildFON2Network(sequences)
time_build_fon2 = time.time() - start_time
fon2_weight = InfoMapClust.uniformNodeWeights(fon2_state_net)

state_codelength, fon2_gain, state_clust, time_state = InfoMapClust.infomapStateClustering(nodes, fon2_states, fon2_state_net, fon2_weight, filename, nb_loop)

final_clusts.append(state_clust)
final_codelengths.append(state_codelength)
final_gain.append(fon2_gain)
final_clusts_name.append("Fon2")
final_clusts_times.append(time_state)
final_nb_dupli.append(getNbDuplication(fon2_states))
final_build_times.append(time_build_fon2)
print('Done.')

####################################
## PRINT CLUSTERING RESULTS       ##
####################################
nb_classes = []
cluster_counts = []
symb_clust_counts = []
for i in range(len(final_clusts)):
	# print(f'{final_clusts_name[i]}')
	clust = final_clusts[i]
	nb_c, symb_count, counts = countClusters(clust)
	nb_non_t = 0
	for c, nb_occ in counts.items():
		if nb_occ > 1:
			nb_non_t += 1
	nb_classes.append(nb_c)
	symb_clust_counts.append(symb_count)
	cluster_counts.append(counts)

print('#################################')
print('Build / Infomap Exec Time')
### Values reported in Table III (page 9)
for i in range(len(nb_classes)):
	print(f'{final_clusts_name[i]}: {round(final_build_times[i],2)} s / {round(final_clusts_times[i],2)} s')

print('#################################')
print('Number of Clusters / Code length / Gain CL:')
for i in range(len(nb_classes)):
	print(f'{final_clusts_name[i]} : {nb_classes[i]} {final_codelengths[i]} {final_gain[i]}')

print('#################################')
print('NMI similarity between clusterings')
### Values reported in Table II (page 8)
for i in range(len(nb_classes)-1):
	for j in range(i+1,len(nb_classes)):
		nmi = OverlappingNMI.NMI_max(final_clusts[i],final_clusts[j])
		print(f'{final_clusts_name[i]} / {final_clusts_name[j]} : {round(nmi,3)}')
print()

print('#################################')
print(f'Writing Stats results in {save_stats_file}')
stats_file = open(save_stats_file,'w')
name_var = 'id'
for i in range(len(final_clusts_name)):
	name_var+= f',{final_clusts_name[i]}_nd'
for i in range(len(final_clusts_name)):
	name_var+= f',{final_clusts_name[i]}_nc'
# name_var += ',nmi'
name_var += '\n'
stats_file.write(name_var)
for symb in symb_clust_counts[0].keys():
	line_symb = f'{symb}'
	for i in range(len(final_clusts_name)):
		if symb in final_nb_dupli[i].keys():
			line_symb += f',{final_nb_dupli[i][symb]}'
	for i in range(len(final_clusts_name)):
		if symb in symb_clust_counts[i].keys():
			line_symb += f',{symb_clust_counts[i][symb]}'
	# line_symb += f',{nmi[symb]}'
	line_symb += '\n'
	stats_file.write(line_symb)
stats_file.close()
print('#################################')
print(f'Writing Clusters in {save_cluster_file}')
clusters_file = open(save_cluster_file,'w')
name_var = 'id'
for i in range(len(final_clusts_name)):
	name_var+= f';{final_clusts_name[i]}'
name_var += '\n'
clusters_file.write(name_var)
for symb in symb_clust_counts[0].keys():
	line_symb = f'{symb}'
	for i in range(len(final_clusts_name)):
		if symb in final_clusts[i].keys():
			line_symb += f';{list(set(final_clusts[i][symb]))}'
	line_symb += '\n'
	clusters_file.write(line_symb)
clusters_file.close()

