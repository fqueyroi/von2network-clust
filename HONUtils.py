# -*- coding: utf-8 -*-
'''
Various functions used in all scripts
'''
import math

##################################################################
def getDistribution(count):
	distr = {}
	support = float(sum(count.values()))
	for target in count.keys():
		distr[target] = float(count[target]) / support
	return distr

######################################################################################
def MaxDivergence(Distr):
	MaxValKey = sorted(Distr, key=Distr.__getitem__)
	d = {MaxValKey[0]: 1}
	return d

######################################################################################
def KLD(dA, dB):
	res = 0
	for target in dA.keys():
		if target in dB.keys():
			res += dA[target] * math.log(dA[target]/dB[target], 2)
	return res

######################################################################################
def KLDThreshold(order, support,ThresholdMultiplier):
	return ThresholdMultiplier * (order / math.log(1 + support, 2))

######################################################################################
def sequenceToString(seq):
	if len(seq) == 1 :
		return str(seq[0]);
	curr = seq[-1]
	res = str(curr) + '/'
	seq = seq[:-1]
	while len(seq) > 0:
		curr = seq[-1]
		res = res + str(curr) + '.'
		seq = seq[:-1]
	if res[-1] == '.':
		return res[:-1]
	else:
		return res
######################################################################################
def readSequenceFile(filename, is_line_id = True, separator=' ') :
	sequences = []
	file_seq = open(filename,'r')
	count_line = 0
	for line in file_seq :
		count_line += 1
		split_line = line.strip().split(separator)
		if len(split_line) > 1 :
			seq = []
			if is_line_id :
				id_seq = split_line[0]
				seq = split_line[1:]
				if len(seq) > 1:
					sequences.append(seq)
			else :
				if len(split_line) > 1 :
					sequences.append(split_line)
	return sequences
######################################################################################
def removeRepetitions(sequences):
	'''
	Remove repetition of symbols in each sequence

	Parameters:
	-----------
	sequences: list of (list of str)

	Returns:
	--------
	nsequences: list of (list of str)
	'''
	nsequences = []
	for index_seq in range(len(sequences)):
		seq = sequences[index_seq]
		nseq = []
		for s in seq:
			if len(nseq) == 0 or s != nseq[-1]:
				nseq.append(s)
		nsequences.append(nseq)
	return nsequences

