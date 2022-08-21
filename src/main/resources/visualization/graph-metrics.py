from sys import argv
from math import sqrt
from igraph import Graph

def variance(values, mean):
	var = 0
	for a in values:
		var = var + (a - mean)**2
	return var / (len(values) - 1)

def compute_metrics(g):
	n = g.vcount()
	m = g.ecount()
	d = g.density()

	deg = g.degree(mode="in")
	mean = sum(deg) / float(n)
	norm_mean = mean / (n - 1)

	var_in = variance(deg, mean)
	std_dev_in = sqrt(var_in)
	norm_var_in = var_in * (n - 1) / (n * m * (1 - d))
	deg = g.degree(mode="out")
	var_out = variance(deg, mean)
	std_dev_out = sqrt(var_out)
	norm_var_out = var_out * (n - 1) / (n * m * (1 - d))

	print("	vertices					", n)
	print("	edges						", m)
	print("	mean degree					", mean)
	print("	norm. mean degree			", norm_mean)
	print("	degree variance	(in)		", var_in)
	print("	norm. degree variance (in)	", norm_var_in)
	print("	degree std.dev (in)			", std_dev_in)
	print("	degree variance	(out)		", var_out)
	print("	norm. degree variance (out)	", norm_var_out)
	print("	degree std.dev (out)		", std_dev_out)
	print("	density						", d)

	cc = g.transitivity_undirected()
	assort = g.assortativity_degree(directed=True)

	print("	clustering coefficient	", cc)
	print("	assortativity			", assort)

	print(n, ";", m, ";", mean, ";", norm_mean, ";", var_in, ";", norm_var_in, ";",
		std_dev_in, ";", var_out, ";", norm_var_out, ";", std_dev_out, ";", d, ";", cc,
		";", assort, sep="")


g = Graph.Read(argv[1], directed=True)
# If the graph represents one undirected edge by two directed edges, the reported metrics
# have to be adjusted as follows: edge count / 2, norm. degree variance * 2

print("metrics of whole graph:")
compute_metrics(g)

clust = g.components("strong")
print("	# strongly connected components			", len(clust))
max = 0
for c in clust:
	if len(c) > max:
		max = len(c)
print("	largest strongly connected component	", max)
clust = g.components("weak")
print("	# weakly connected components			", len(clust))
max = 0
histo = {}
for c in clust:
	l = len(c)
	if l > max:
		max = l
	if l in histo:
		histo[l] = histo[l] + 1
	else:
		histo[l] = 1
print("	largest weakly connected component		", max)
print(" histogram: ", histo)

print("metrics of largest weakly connected component:")
g = clust.giant()
compute_metrics(g)
