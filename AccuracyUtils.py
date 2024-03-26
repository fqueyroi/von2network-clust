# -*- coding: utf-8 -*-
'''
Contains functions to prepare
datasets and for testing the prediction capabilities
of the models
'''
import random

################################################################################
def sampleSequences(sequences, ratio_testing):
	'''
	Cut sequences into training and testing subsets

	Parameters:
	-----------
	sequences: list of (list of str)

	Returns:
	--------
	build_seqs: list of (list of str)
	test_seqs : list of (list of str)
	'''
	build_seqs = []
	test_seqs  = []
	for i in range(len(sequences)):
		if random.random() < ratio_testing:
			test_seqs.append(sequences[i])
		else:
			build_seqs.append(sequences[i])
	return build_seqs, test_seqs

################################################################################
def numberOfRulesOfOrder(rules,order=1):
	'''
	Parameters:
	-----------
	rules: dict list of str -> dict (str -> float)
	'''
	res = 0
	for k in rules.keys():
		if len(k) == order:
			res += 1
	return res
################################################################################
def numberOfRulesAgg(agg2ndOrderRules):
	'''
	Parameters:
	-----------
	agg2ndOrderRules: dict of str -> dict list of str -> dict (str -> float)
			Node x -> clustering of duplication of x -> union count of rules
			in the group
	'''
	res = 0
	for k in agg2ndOrderRules.keys():
		res += len(agg2ndOrderRules[k].keys())
	return res


################################################################################
def probabilityNextSymb(rules, context, symb):
	'''
	Returns the probability associated with symbol 'symb' given the 'context'
	according to the 'rules'

	Parameters:
	-----------
	rules:  dict tuple of str -> dict (str -> float)
	context: list of str
	symb: str

	Returns:
	--------
	float
	'''
	t_con = tuple(context)
	if t_con not in rules.keys():
		if len(t_con)>1:
			context_suffix = context[1:]
			return probabilityNextSymb(rules, context_suffix, symb)
		else:
			return 0.
	count = rules[t_con]
	count_symb = 0
	if symb in count.keys():
		count_symb = count[symb]
	return (1.*count_symb) / sum(count.values())

################################################################################
def nonAggRulesProbSeq(rules, seq):
	'''
	Compute the average probability of each element s[i] in 'seq' using
	as context the symbols s[i-2],s[i-1]

	Returns:
	--------
	float
	'''
	if len(seq)<2:
		return 0.
	res = 0
	for i in range(1,len(seq)):
		#print "Context : "+str(seq[max(0,i-2):i])+" symb : "+str(seq[i])
		p = probabilityNextSymb(rules,seq[max(0,i-2):i],seq[i])
		#print "		p  = "+str(p)
		res += p
	return  res / (len(seq) - 1.)

################################################################################
def aggProbabilityNextSymb(firstOrderRules, agg2ndOrderRules, context, symb):
	'''
	Returns the probability associated with symbol 'symb' given the 'context'
	according to 1st order rules 'firstOrderRules'
	and aggregated 2nd order rules agg2ndOrderRules

	Parameters:
	-----------
	firstOrderRules:  dict tuple of str -> dict (str -> float)
	agg2ndOrderRules: dict of str -> dict tuple of tuple of str -> dict (str -> float)
	context: list of str
	symb: str

	Returns:
	--------
	float
	'''
	t_con = tuple(context)
	if len(t_con) == 0:
		return 0.
	count = {}
	if len(t_con) == 1:
		if t_con not in firstOrderRules.keys():
			return 0.
		count = firstOrderRules[t_con]
	if len(t_con) == 2:
		last_symb = t_con[-1]
		if last_symb not in agg2ndOrderRules.keys():
			return probabilityNextSymb(firstOrderRules, tuple([last_symb]), symb)
		cluster_found = False
		for clusts, u_counts in agg2ndOrderRules[last_symb].items():
			if t_con in clusts:
				count = u_counts
				cluster_found = True
				# print "Cluster found : "+str(clusts)
				break
		if not cluster_found:
			return probabilityNextSymb(firstOrderRules, tuple([last_symb]), symb)

	count_symb = 0
	if symb in count.keys():
		count_symb = count[symb]
	return (1.*count_symb) / sum(count.values())

################################################################################
def aggRulesProbSeq(firstOrderRules, agg2ndOrderRules, seq):
	'''
	Compute the average probability of each element s[i] in 'seq' using
	as context the symbols s[i-2],s[i-1]
	for the Aggregated 2nd order network

	Returns:
	--------
	float
	'''
	if len(seq)<2:
		return 0.
	res = 0
	for i in range(1,len(seq)):
		p = aggProbabilityNextSymb(firstOrderRules, agg2ndOrderRules, seq[max(0,i-2):i],seq[i])
		res += p
	return  res / (len(seq) - 1.)

################################################################################
def fonProbabilityNextSymbol(order2rules, context, symb):
	if len(context) < 2:
		return 0.
	i, j = context[-2], context[-1]
	if i not in order2rules.keys():
		return 0.
	o2i = order2rules[i]
	if j not in o2i.keys():
		return 0.
	o2ij = o2i[j]
	if symb not in o2ij.keys():
		return 0.
	return o2ij[symb]

################################################################################
def fonProbSeq(order2rules, seq):
	if len(seq)<3:
		return 0.
	res = 0.
	for i in range(2,len(seq)):
		context = seq[(i-2):i]
		p = fonProbabilityNextSymbol(order2rules, context, seq[i])
		res += p
	return res / (len(seq) - 2.)

