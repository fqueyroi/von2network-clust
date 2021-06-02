# -*- coding: utf-8 -*-
'''
Take the results of BuildRulesFast.py (with max_order=2)
and compute the aggregation of 2nd order rules
(Alg. 1 of the paper)

main method: aggregateRules(..)

dependencies:
HeapDict https://pypi.org/project/HeapDict/
'''

import HONUtils
import math, sys
import heapdict

from collections import defaultdict

##################################################################
def lastSymbolMapping(rules):
	'''
	Associate each rule of length > 1 in 'rules' to the last symbol
	of the rule i.e. lastSymbolMapping(rules)[x] contains the rules
	that have x as last symbol

	Parameters:
	-----------
	rules: dict tuple of str -> dict (str -> float)

	Returns:
	--------
	mapping: dist str -> (dict tuple of str -> dict (str -> float))
	'''
	mapping = defaultdict(dict)
	for rule in rules.keys():
		if len(rule) > 1 :
			symb = rule[-1]
			mapping[symb][rule] = rules[rule]
	return mapping

##################################################################
def getUnionCount(count1, count2):
	'''
	Return the union of count1 and count2

	Parameters:
	-----------
	count1: map str -> float
	count2: map str -> float
	'''
	res = {}
	for target in set(count1.keys()) | set(count2.keys()):
		res[target] = 0.
		if target in count1.keys():
			res[target] += count1[target]
		if target in count2.keys():
			res[target] += count2[target]
	return res

##################################################################
def aggregationScore(c1, c2, dp): # Return bool, float
	'''
	Parameters:
	-----------
	cX: count of rules X (dict str -> float)
	dp: distribution of the parent rule of the two rules

	Return:
	-------
	can_be_merge : bool (can the two rules be merged?)
	score        : float (proximity between the two rules)
	'''

	d1 = HONUtils.getDistribution(c1)
	d2 = HONUtils.getDistribution(c2)

	c1u2 = getUnionCount(c1,c2)
	s1, s2, s1u2 = sum(c1.values()), sum(c2.values()), sum(c1u2.values())
	d1u2 = HONUtils.getDistribution(c1u2)

	kld1m, thres1m = HONUtils.KLD(d1,d1u2), HONUtils.KLDThreshold(2,s1)
	kld2m, thres2m = HONUtils.KLD(d2,d1u2), HONUtils.KLDThreshold(2,s2)
	kld1u2p, thres1u2p = HONUtils.KLD(d1u2,dp), HONUtils.KLDThreshold(2,s1u2)

	can_be_merge = kld1u2p > thres1u2p and kld1m < thres1m and kld2m < thres2m
	score = kld1m + kld2m # - 2.*kld1u2p
	return can_be_merge, score

##################################################################
def aggregate(rules, dp):
	'''
	Aggregate 2nd-order rules that have the same last symbol "symb"
	using an hierarchical clustering procedure

	Parameters:
	-----------
	rules: dict list of str -> dict (str -> float)

	Returns:
	--------
	clust: dict tuple of tuple of str -> dict (str -> float)
		   group of rules -> union count of rules in the group
	'''

	clust = defaultdict(dict) ## rules cluster -> cluster count
	for r in rules:
		## Initialisation: each rule in one cluster
		clust[tuple([r])] = rules[r]

	if len(clust.keys()) < 2:
		## No aggregation possible if only two rules
		return clust

	# Create a heap with the pair of clusters that can be merged
	# using the aggregation score
	hd = heapdict.heapdict()
	keys_clust = list(clust.keys())
	for i1 in range(len(clust)-1):
		r1 = keys_clust[i1]
		c1 = clust[r1]
		for i2 in range(i1 + 1, len(clust)):
			r2 = keys_clust[i2]
			c2 = clust[r2]
			can_merge, score = aggregationScore(c1, c2, dp)
			if can_merge:
				# print "Can merge !"
				hd[r1, r2] = score

	while len(hd) > 0:
		# while there are possible merges
		# pick the best pair and merge it
		best_pair, score = hd.popitem()
		# print 'best = '+str(best_pair)+" score : "+str(score)
		new_grp_rules = []
		r1, r2 = best_pair
		for sr1 in r1:
			new_grp_rules.append(sr1)
		for sr2 in r2:
			new_grp_rules.append(sr2)
		c1u2 = getUnionCount(clust[r1],clust[r2])

		del clust[r1]
		del clust[r2]
		## update heapdict removing entries with best_pair
		keys_hd = list(hd.keys())
		for hpair in keys_hd:
			if hpair[0] == r1 or hpair[1] == r1 or hpair[0] == r2 or hpair[1] == r2:
				del hd[hpair]

		## update heapdict with new pairs that can be merged
		for r in clust.keys():
			can_merge, score = aggregationScore(clust[r],c1u2,dp)
			if can_merge:
				hd[r, tuple(new_grp_rules)] = score
		clust[tuple(new_grp_rules)] = c1u2

	return clust

##################################################################
def aggregateRules(rules):
	'''
	Main method
	Aggregate detected 2nd order rules that are extension of the same
	1st order rule.

	Parameters:
	-----------
	rules: dict tuple of str -> dict (str -> float)

	Returns:
	--------
	clusts: dict of str -> dict tuple of tuple of str -> dict (str -> float)
		Node x -> clustering of duplications of x -> union count of rules
		in the group
	'''
	firstOrderMap = lastSymbolMapping(rules)

	clusts = defaultdict(dict)
	for firstOrderRule in firstOrderMap.keys():
		subRules = firstOrderMap[firstOrderRule]
		cFirstOrder  = rules[tuple([firstOrderRule])]
		dp = HONUtils.getDistribution(cFirstOrder)
		subClust = aggregate(subRules,dp)
		clusts[firstOrderRule] = subClust
	return clusts

##################################################################
def flattenAgg2ndOrderRules(base_rules, clusts):
	'''
	Returns:
	--------
	flatten_rules: dict tuple of tuple of str ->  dict (str -> float)
	'''
	flatten_rules = {}
	for r, distr in base_rules.items():
		if len(r) == 1:
			flatten_rules[tuple([r])] = distr

	for k, k_clusts in clusts.items():
		for subclust, distr in k_clusts.items():
			flatten_rules[subclust] = distr
	return flatten_rules

##################################################################
def mergeNodes(graph, flatten_agg_rules):
	new_graph = {}
	nodeToClust = {}
	for c in flatten_agg_rules.keys():
		new_graph[c] = {}
		for r in c:
			nodeToClust[r] = c

	for src in graph.keys():
		for tgt in graph[src].keys():
			if tuple(tgt) not in nodeToClust.keys():
				new_graph[tuple([tgt])] = {}
				nodeToClust[tgt] = tuple([tgt])
			if nodeToClust[tgt] not in new_graph[nodeToClust[src]].keys():
				new_graph[nodeToClust[src]][nodeToClust[tgt]] = 0
			new_graph[nodeToClust[src]][nodeToClust[tgt]] += graph[src][tgt]
	return new_graph
