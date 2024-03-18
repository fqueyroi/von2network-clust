# -*- coding: utf-8 -*-
'''
Main script used to evaluate the aggregation of 2nd order representations
Compute the accuracy of the VON2, Agg VON2 and FON2 models 
on a testing subset of sequences
'''
from collections import defaultdict

import numpy as np

import HONUtils
import AccuracyUtils

import BuildRulesFast
import AggOrder2Rules
import FON2StatesNetwork
import time

##########################
## SEQUENCES FILE INPUT ##
##########################

## Computing results on the Maritime Dataset
# filename = "./maritime_sequences.csv"
# sep = " " ## the string separating elements in a sequence
# threshold_multi = 2.8 # Alpha, to adjust size of VON2 to Agg-VON2
# filename = "./ports_2009.ngram"
# sep = ","
# threshold_multi = 3.05

## Computing results on the Airports Dataset
#filename = "./2011Q1_SEQ.csv"
filename = "./US_flights.ngram"
sep = ","
threshold_multi = 2.8 # Alpha, to adjust size of VON2 to Agg-VON2

## Computing results on the Taxis Dataset
#filename = "./trajectories_PoliceStation.csv"
#sep = " "
# threshold_multi = 4.4 # Alpha, to adjust size of VON2 to Agg-VON2

## Computing results on the Wikispeedia dataset
# filename = "./wikispeedia_top100.ngram"
# sep = ","
# threshold_multi = 2.4 

## Computing results on the MSNBC dataset
# filename = "./msnbc990928.ngram"
# sep = " "
# threshold_multi = 2.5

## Read trajectories file
sequences = HONUtils.readSequenceFile(filename, False, sep)
sequences = HONUtils.removeRepetitions(sequences)


#########################
## Parameters          ##
#########################

nb_run = 50 ## Number of test for each network (eg 50)
testing_ratio = 0.1

## Perform 'nb_run' tests taking 'testing_ratio' of the input sequences
## as testing subset
## The last 3 cols printed are the Acc (Eq.10 in the paper)
## for the Von2, Von2(alpha*), Agg Von2 et Fon2 networks respectively
## alpha* threshold given such that Von2(alpha*) is as sparse as Agg Von2 (Eq. 2)

print('id_run,nb_symb,nb_seq_build,nb_seq_test,nn_von2,nn_sparse,nn_agg,nn_fon2,score_von2,score_sparse,score_agg,score_fon2')

for i in range(nb_run):
	## Split into Training and Testing sets
	training, testing = [], []
	training, testing = AccuracyUtils.sampleSequences(sequences, testing_ratio)

	## Build relevant order 2 extensions (Von2 network)
	rule_builder = BuildRulesFast.FastHONRulesBuilder(training, max_order=2, min_support=1, ThresholdMultiplier=1.)
	rules = rule_builder.ExtractRules()
	## Aggregate 2nd order extension (Agg Von2 network)
	clusts = AggOrder2Rules.aggregateRules(rules)
	flatten_agg_rules = AggOrder2Rules.flattenAgg2ndOrderRules(rules, clusts)
	## Build VON2 using the given threshold multiplier choosen such that
	## nb of representations is close to Agg-Von2
	sparse_rule_builder = BuildRulesFast.FastHONRulesBuilder(training, max_order=2, min_support=1, ThresholdMultiplier=threshold_multi)
	sparse_rules = sparse_rule_builder.ExtractRules()
	## Build Fix order 2 extensions (Fon2 network)
	nodes, fon2_states, fon2_state_net = FON2StatesNetwork.buildFON2Network(training)
	fon2rules = FON2StatesNetwork.order2Rules(training)


	#############
	## Results ##
	#############

	## Extract 1st rules
	## (for tests using the Agg Von2 model)
	firstOrderRules = {}
	avg = {}
	for r in rules:
		if len(r)==1:
			firstOrderRules[r] = rules[r]
			avg[r] = 0
	for r in rules:
		if len(r) == 1:
			for p in rules:
				if len(p) > 1 and p[1] == r[0]:
					avg[r] += 1
	## Compute the Accuracy capabilities of the three HON networks
	## Eq. (10) in the paper
	score_non_agg = 0.
	score_agg     = 0.
	score_fon2    = 0.
	score_sparse  = 0.
	nb_test_length3 = 0.
	for test_index in range(len(testing)):
		test_seq = testing[test_index]
		if len(test_seq)>=3:
			score_non_agg += AccuracyUtils.nonAggRulesProbSeq(rules,test_seq)
			score_agg     += AccuracyUtils.aggRulesProbSeq(firstOrderRules,clusts,test_seq)
			score_fon2    += AccuracyUtils.fonProbSeq(fon2rules,test_seq)
			score_sparse  += AccuracyUtils.nonAggRulesProbSeq(sparse_rules,test_seq)
			nb_test_length3 += 1.
	if nb_test_length3 > 0:
		score_non_agg = score_non_agg / nb_test_length3
		score_agg     = score_agg / nb_test_length3
		score_fon2    = score_fon2 / nb_test_length3
		score_sparse  = score_sparse / nb_test_length3
		
	## size Alphabet
	A = set()
	for seq in training:
		for s in seq:
			A.add(s)

	# nn_fon2 = len(fon2_state_net)
	nn_fon2 = len(A)
	for state_ij, out_neigh in fon2_state_net.items():
		if len(out_neigh)>0:
			nn_fon2 += 1


	## Outputs
	print(f'{i+1},{len(A)},{len(training)},{int(nb_test_length3)},{len(rules)},{len(sparse_rules)},{len(flatten_agg_rules)},{nn_fon2},{score_non_agg},{score_sparse},{score_agg},{score_fon2}')
