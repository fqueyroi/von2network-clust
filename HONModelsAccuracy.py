# -*- coding: utf-8 -*-
'''
Main script used to evaluate the aggregation of 2nd order representations
Compute the accuracy of the VON2, Agg VON2 and FON2 models 
on a testing subset of sequences
'''
import HONUtils
import AccuracyUtils

import BuildRulesFast
import AggOrder2Rules
import FON2StatesNetwork

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
# filename = "./trajectories_PoliceStation.csv"
# sep = " "

## Read trajectories file
sequences = HONUtils.readSequenceFile(filename,True,sep)
sequences = HONUtils.removeRepetitions(sequences)

#########################
## Parameters          ##
#########################

nb_run = 50 ## Number of test for each network (eg 50)
testing_ratio = 0.1

## Perform 'nb_run' tests taking 'testing_ratio' of the input sequences
## as testing subset
## The last 3 cols printed are the Acc (Eq.10 in the paper)
## for the Von2, Agg Von2 et Fon2 networks respectively
print('id_run,nb_seq_build,nb_seq_test,score_von2,score_agg,score_fon2')
for i in range(nb_run):
	## Split into Training and Testing sets
	training, testing = [], []
	training, testing = AccuracyUtils.sampleSequences(sequences, testing_ratio)

	## Build relevant order 2 extensions (Von2 network)
	rule_builder = BuildRulesFast.FastHONRulesBuilder(training,2,1)
	rules = rule_builder.ExtractRules()

	## Aggregate 2nd order extension (Agg Von2 network)
	clusts = AggOrder2Rules.aggregateRules(rules)
	flatten_agg_rules = AggOrder2Rules.flattenAgg2ndOrderRules(rules, clusts)

	## Build Fix order 2 extensions (Fon2 network)
	fon2rules = FON2StatesNetwork.order2Rules(training)

	#############
	## Results ##
	#############

	## Extract 1st rules
	## (for tests using the Agg Von2 model)
	firstOrderRules = {}
	for r in rules:
		if len(r)==1:
			firstOrderRules[r] = rules[r]

	## Compute the Accuracy capabilities of the three HON networks
	## Eq. (10) in the paper
	score_non_agg = 0.
	score_agg     = 0.
	score_fon2    = 0.
	nb_test_length3 = 0.
	for test_index in range(len(testing)):
		test_seq = testing[test_index]
		if len(test_seq)>=3:
			score_non_agg += AccuracyUtils.nonAggRulesProbSeq(rules,test_seq)
			score_agg     += AccuracyUtils.aggRulesProbSeq(firstOrderRules,clusts,test_seq)
			score_fon2    += AccuracyUtils.fonProbSeq(fon2rules,test_seq)
			nb_test_length3 += 1.
	score_non_agg = score_non_agg / nb_test_length3
	score_agg     = score_agg / nb_test_length3
	score_fon2    = score_fon2 / nb_test_length3

	## Outputs
	print(f'{i+1},{len(training)},{int(nb_test_length3)},{score_non_agg},{score_agg},{score_fon2}')
