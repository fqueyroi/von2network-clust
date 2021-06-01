# -*- coding: utf-8 -*-
'''
Implementation of Overlapping NMI similarity measure
between clustering

Original code at https://github.com/aaronmcdaid/Overlapping-NMI

main method:  NMI_max(clust1, clust2)
returns a value in [0,1]
where 1 means clust1 and clust2 are equal
	  0 means they are completly different
'''

from collections import defaultdict
import math

def h(w,n):
	if w==0:
		return 0.
	return -(w+0.)*math.log((w+0.)/n,2)

def clusteringAsSets(clust):
	res = defaultdict(set)
	for n, clusts_n in clust.items():
		for c in clusts_n:
			res[c].add(n)
	return list(res.values())

def entropyClustering(n, X):
	## H(X)
	res = 0.
	for Xi in X:
		res += h(len(Xi),n) + h(n-len(Xi),n)
	return res

def countInters(n, Xi, Yj):
	a = n - len(Xi | Yj)
	b = len(Yj - Xi)
	c = len(Xi - Yj)
	d = len(Xi & Yj)
	return a, b, c, d

def lackInfo(n, Xi, Yj):
	## H(X_i|Y_j)
	a, b, c, d = countInters(n, Xi, Yj)
	## Apply Eq. 2
	if  h(a,n) + h(d,n) >= h(b,n) + h(c,n):
		## uncorrected Eq. 1
		lack_un = h(a,n)+h(b,n)+h(c,n)+h(d,n)-h(b+d,n)-h(a+c,n)
		return lack_un
	else:
		return h(c+d,n)+ h(a+b,n)

def lackInfoMatch(n, Xi, Y):
	## Apply Eq. 3
	min_lack = n*n
	for Yj in Y:
		lack_j = lackInfo(n, Xi, Yj)
		if lack_j < min_lack:
			minY = Yj
			min_lack = lack_j
	return min_lack

def NMI_max(clust1, clust2):
	n = len(clust1.keys())
	X = clusteringAsSets(clust1)
	Y = clusteringAsSets(clust2)

	hX = entropyClustering(n, X)
	hY = entropyClustering(n, Y)

	## Apply Eq. 4 for X and Y
	hXY = 0.
	for Xi in X:
		hXY += lackInfoMatch(n, Xi, Y)
	hYX = 0.
	for Yj in Y:
		hYX += lackInfoMatch(n, Yj, X)

	## Apply Eq. 5
	iXY = (hX - hXY + hY - hYX) / 2.

	## Apply Eq. 10
	nmi = iXY / max(hX,hY)

	return nmi
