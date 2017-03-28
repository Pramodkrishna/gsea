#!/usr/bin/env python

# Matjaz Zganec, 2017

"""
Gene Set Enrichment Analysis

Implementation of gene set enrichment analysis as described in
Subramanian et al., Gene set enrichment analysis: A knowledge-based approach
for interpreting genome-wide expression profiles, PNAS 102 (43): 15545-15550,
2005.

The most important functions are:
(1) gsea(...), which computes normalized enrichment scores and the
    corresponding p-values based on the data in gene expression profiles.
(2) plot(...), which displays Fig. 1 from the original paper: (a) a heat map
    of gene expression profiles grouped by class, (b) gene expression
    correlation profile, (c) gene set member position in the profiles, and
    (d) the running sum (see the original paper).
(3) my_gsea(...), which computes estimates of normalized enrichment scores
    and the corresponding p-vlaues without gene expression profiles. Instead
    a user-provided set of interesting genes is used.
"""

import csv
import math
import matplotlib.pyplot as plt
import matplotlib.patches
import multiprocessing
import numpy
import operator
import scipy.stats

def get_gene_expressions(filename, delimiter='\t'):
	"""
	Parse a file with gene expression profiles. The file should contain a
	header row specifying the class label of each gene expression profile.
	After that, the first column in the file contains gene labels and the
	rest contain numeric values associated with each gene expression profile.

	Args:
		filename (str): Path to the input file.
		delimiter (str): Delimiter between fields in the input file.

	Returns:
		list: Gene expression profiles with the following structure:
			gene_expressions[gene_index][profile_index].
		list: Class labels associated with the profiles.
		list: Gene labels associated with the profiles.
	"""
	with open(filename) as f:
		reader = csv.reader(f, delimiter=delimiter)
		header = next(reader)
		class_labels = [label for label in header[1:]]
		gene_labels = []
		gene_expressions = []
		for row in reader:
			try:
				gene_expressions.append(list(map(float, row[1:])))
				gene_labels.append(row[0])
			except:
				print('Parsing error in line', row)
	return gene_expressions, class_labels, gene_labels

def get_gene_sets(filename, delimiter='\t'):
	"""
	Parse a file with gene sets e.g., genes encoding products in a metabolic
	pathway. The first column in the file should contain a gene set label,
	followed by a description in the second column. All other columns specify
	a list of gene labels in each gene set, one set per line.

	Args:
		filename (str): Path to the input file.
		delimiter (str): Delimiter between fields in the input file.

	Returns:
		list: Gene sets with the following structure:
			gene_sets[set_index].
		list: Gene set labels.
	"""
	with open(filename) as f:
		reader = csv.reader(f, delimiter=delimiter)
		gene_sets = list()
		gene_set_labels = list()
		for row in reader:
			gene_set_labels.append(row[0])
			gene_sets.append(set(row[2:]))
	return gene_sets, gene_set_labels

def get_my_genes(filename, delimiter='\t'):
	"""
	Parse a file containing a list of gene labels, one gene label per line.

	Args:
		filename (str): Path to the input file.
		delimiter (str): Delimiter between fields in the input file.

	Returns:
		set: Set of gene labels listed in the input file.
	"""
	with open(filename) as f:
		reader = csv.reader(f, delimiter=delimiter)
		my_gene_labels = set()
		for row in reader:
			my_gene_labels.add(row[0])
	return my_gene_labels

def corr(x, y):
	"""
	Compute Pearson correlation coefficient between gene expressions and
	profile class.

	Args:
		x (list): Expressions of a single gene across profiles.
		y (list): Zero-based encoded class (0.0 for profiles in the first class
			and 1.0 for profiles in the second class).

	Returns:
		Pearson correlation coefficient betweem gene expressions and profile
		class.
	"""
	return scipy.stats.pearsonr(x, y)[0]

def get_gene_ranks(gene_expressions, class_labels, gene_labels, func=corr):
	"""
	Rank genes in the profiles by correlation between their expression and
	profile class.

	Args:
		gene_expressions (list): Gene expression profiles with the following
			structure: gene_expressions[gene_index][profile_index].
		class_labels (list): Class labels associated with the profiles.
		gene_labels (list): Gene labels associated with the profiles.
		func (function): Quantifier of the relationship between gene expressions
			and class. The function accepts two arguments: (1) a list of
			expressions of a single gene accross profiles and (2) a list of
			zero-based encoded class (0.0 for profiles in the first class
			and 1.0 for profiles in the second class).

	Returns:
		list: Tuples, each containing a gene label and the corresponding
			correlation coefficient. The list is sorted by the value of
			correlation coefficient.
	"""
	classes = sorted(list(set(class_labels)))
	class_indices = list(map(classes.index, class_labels))
	gene_metrics = [func(gene_expression, class_indices) \
		for gene_expression in gene_expressions]
	return sorted(zip(gene_labels, gene_metrics), key=lambda item: item[1])

def get_enrichment_running_sum(gene_ranks, gene_set, p=1.0):
	"""
	Compute the running sum used in the calculation of enrichment score.

	Args:
		gene_ranks (list): Tuples, each containing a gene label and the
			corresponding correlation coefficient. The list is sorted by the
			value of correlation coefficient.
		gene_set (set): Labels of genes to be used in the calculation
			of enrichment running sum.
		p (float): Exponent used in the calculation of the running sum (see
			the original paper).

	Returns:
		list: Running sum used in the calculation of enrichment score.
	"""
	n_r = sum([math.pow(abs(gene_rank[1]), p) for gene_rank in gene_ranks \
		if gene_rank[0] in gene_set])
	n_intersection = len(set([gene_rank[0] \
		for gene_rank in gene_ranks]).intersection(gene_set))
	d_miss = 1.0 / (len(gene_ranks) - n_intersection)
	if n_r == 0.0:
		return [0.0 for gene_rank in gene_ranks]
	running_sum = []
	p_hit = 0.0
	p_miss = 0.0
	for gene_rank in gene_ranks:
		if gene_rank[0] in gene_set:
			p_hit += math.pow(abs(gene_rank[1]), p) / n_r
		else:
			p_miss += d_miss
		running_sum.append(p_hit - p_miss)
	return running_sum

def get_enrichment_score(running_sum):
	"""
	Compute enrichment score of a gene set from the running sum.

	Args:
		running_sum (list): The running sum obtained from
			get_enrichment_running_sum(gene_ranks, gene_set, p=1.0).

	Returns:
		float: Enrichment score value.
	"""
	max_val = max(running_sum)
	min_val = min(running_sum)
	return max_val if abs(max_val) > abs(min_val) else min_val

def get_enrichment_scores(gene_ranks, gene_sets, p=1.0):
	"""
	Compute enrichment scores of a list of gene sets.

	Args:
		gene_ranks (list): Tuples, each containing a gene label and a
			correlation coefficient. The list is sorted by the value of
			correlation coefficient.
		gene_sets (list): Sets of labels of genes to be used in the
			calculation of enrichment scores.
		p (float): Exponent used in the calculation of the running sums (see
			the original paper).

	Returns:
		list: Enrichment score values.
	"""
	return [get_enrichment_score(get_enrichment_running_sum(gene_ranks, \
		gene_set, p)) for gene_set in gene_sets]

def gsea(gene_expressions_fn, gene_sets_fn, func=corr, p=1.0, n_trials=1000, \
	n_jobs=multiprocessing.cpu_count()):
	"""
	Compute normalized enrichment scores and the corresponding p-values.
	The results are sorted by normalized enrichment score and printed out
	in three columns: (1) gene set label, (2) normalized_enrichment_score,
	and (3) statistical significance (p-value).

	Args:
		gene_expressions_fn (str): Name of the file with gene expression
			profiles data.
		gene_sets_fn (str): Name of the file with a list of gene sets.
		func (function): Quantifier of the relationship between gene expressions
			and class. The function accepts two arguments: (1) a list of
			expressions of a single gene accross profiles and (2) a list of
			zero-based encoded class (0.0 for profiles in the first class
			and 1.0 for profiles in the second class).
		p (float): Exponent used in the calculation of the running sum (see
			the original paper).
		n_trials (int): number of permutations of profile classes used in
			normalization of enrichment scores and calculation of the
			corresponding p-values.
		n_jobs (int): number of threads used in the computation of permuted
			gene ranks and permuted enrichment scores.

	Returns:
		None
	"""
	gene_expressions, class_labels, gene_labels = \
		get_gene_expressions(gene_expressions_fn)
	gene_sets, gene_set_labels = get_gene_sets(gene_sets_fn)
	gene_ranks = get_gene_ranks(gene_expressions, class_labels, gene_labels, \
		func)
	enrichment_scores = get_enrichment_scores(gene_ranks, gene_sets, p)
	permuted_class_labels = list(map(numpy.random.permutation, (class_labels \
		for _ in range(n_trials))))
	with multiprocessing.Pool(n_jobs) as pool:
		permuted_gene_ranks = pool.starmap(get_gene_ranks, \
			((gene_expressions, permuted_class_labels[i], gene_labels, func) \
			for i in range(n_trials)))
		permuted_enrichment_scores = pool.starmap(get_enrichment_scores, \
			((permuted_gene_ranks[i], gene_sets, p) for i in range(n_trials)))
	normalized_enrichment_scores = []
	significances = []
	permuted_enrichment_scores = numpy.array(permuted_enrichment_scores)
	for i in range(len(gene_set_labels)):
		signed_enrichment_scores = [score \
			for score in permuted_enrichment_scores[:,i] \
			if numpy.sign(score) == numpy.sign(enrichment_scores[i])]
		average_enrichment_score = \
			sum(signed_enrichment_scores) / len(signed_enrichment_scores) \
			if len(signed_enrichment_scores) != 0 else 0.0
		normalized_enrichment_score = \
			enrichment_scores[i] / average_enrichment_score \
			if average_enrichment_score != 0.0 else 0.0
		normalized_enrichment_scores.append(normalized_enrichment_score)
		p_value = float(len(list(filter(lambda score: \
			abs(score) >= abs(enrichment_scores[i]), signed_enrichment_scores)))) \
			/ len(signed_enrichment_scores) \
			if len(signed_enrichment_scores) != 0 else 1.0
		significances.append(p_value)
	results = sorted(zip(gene_set_labels, normalized_enrichment_scores, \
		significances), key=lambda item: item[1], reverse=True)
	print('gene_set_label', 'normalized_enrichment_score', \
		'statistical_significance', sep='\t')
	for result in results:
		print(result[0], '%.3f' % result[1], '%.3f' % result[2], sep='\t')

def plot(gene_expressions_fn, gene_set=set(), func=corr, p=1.0):
	"""
	Display all plots in Fig. 1 of the original paper:
	(1) a heat map of gene expression profiles grouped by class,
	(2) gene expression correlation profile,
	(3) gene set member position in the profiles,
	(4) the running sum (see the original paper).

	Args:
		gene_expressions_fn (str): Name of the file with gene expression
			profiles data.
		gene_set (set): Labels of genes in the gene set used in the analysis.
		func (function): Quantifier of the relationship between gene expressions
			and class; the function accepts two arguments: (1) a list of
			expressions of a single gene accross profiles and (2) a list of
			zero-based encoded class (0.0 for profiles in the first class
			and 1.0 for profiles in the second class).
		p (float): Exponent used in the calculation of the running sum (see
			the original paper).

	Returns:
		None
	"""
	gene_expressions, class_labels, gene_labels = \
		get_gene_expressions(gene_expressions_fn)
	classes = sorted(list(set(class_labels)))
	gene_ranks = get_gene_ranks(gene_expressions, class_labels, gene_labels, \
		func)
	data = []
	for gene_rank in gene_ranks:
		values = [gene_expression for gene_expression, class_label in \
			zip(gene_expressions[gene_labels.index(gene_rank[0])], class_labels) \
			if class_label == classes[0]] + [gene_expression \
			for gene_expression, class_label in \
			zip(gene_expressions[gene_labels.index(gene_rank[0])], class_labels) \
			if class_label == classes[1]]
		data.append(values)
	for i in range(len(data)):
		min_v = min(data[i])
		max_v = max(data[i])
		data[i] = [(v - min_v) / (max_v - min_v) for v in data[i]]
	f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True)
	ax1.imshow(data, cmap='hot', interpolation='nearest')
	ax1.axvline(x=class_labels.count(classes[0]) - 0.5, linewidth=2)
	ax1.set_aspect('auto')
	ax1.set_xticklabels([])
	ax1.set_yticklabels([])
	ax1.set_ylabel('Ranked Gene List')
	ax1.set_title('Gene Expression Profiles Grouped by Class')
	ax1.axis([-0.5, len(class_labels) - 0.5, -0.5, len(gene_labels) - 0.5])
	ax2.plot([gene_rank[1] for gene_rank in gene_ranks], \
		[i for i in range(len(gene_ranks))])
	ax2.set_yticklabels([])
	ax2.axis([-1.0, 1.0, -0.5, len(gene_labels) - 0.5])
	ax2.set_xlabel('Correlation Coefficient')
	ax2.set_title('Gene Expression Correlation Profile')
	indices = [i for i in range(len(gene_ranks)) \
		if gene_ranks[i][0] in gene_set]
	for i in indices:
		ax3.axhline(y=i)
	ax3.set_xticklabels([])
	ax3.set_yticklabels([])
	ax3.set_title('Gene Set Member Positions')
	running_sum = get_enrichment_running_sum(gene_ranks, gene_set, p)
	abs_running_sum = [abs(v) for v in running_sum]
	index, _ = max(enumerate(abs_running_sum), key=operator.itemgetter(1))
	ax4.plot(running_sum, range(len(running_sum)))
	e = matplotlib.patches.Ellipse(xy=(running_sum[index], index), \
		height=0.05 * len(running_sum), width=0.1 * abs_running_sum[index], \
		color='yellow')
	ax4.add_artist(e)
	ax4.set_xlabel('Running Sum')
	ax4.set_title('Random Walk')
	ax4.axis([-1.2 * abs_running_sum[index], 1.2 * abs_running_sum[index], \
		-0.5, len(gene_labels) - 0.5])
	plt.show()

def get_my_enrichment_score(my_gene_labels, gene_set, p=1.0):
	"""
	Compute an estimate of enrichment score without gene expression profiles
	data. A set of interesting genes (my_gene_labels) is provided instead.

	Args:
		my_gene_labels (set): A set of gene labels appearing near the top or
			near the bottom of the list with ranked genes.
		gene_set (set): A set of genes to calculate enrichment score for.
		p (float): Exponent used in the calculation of the running sum (see
			the original paper).

	Returns:
		float: An estimate of enrichment score value.
	"""
	if not gene_set:
		return 0.0
	common_gene_labels = my_gene_labels.intersection(gene_set)
	q = float(len(common_gene_labels)) / len(gene_set)
	if p > 0.0:
		return 1.0 if common_gene_labels else 0.5
	if p == 0.0:
		return q if q >= 0.5 else q - 0.5
	if p < 0.0:
		return 0.5 if q < 1.0 else 1.0

def get_my_enrichment_scores(my_gene_labels, gene_sets, p=1.0):
	"""
	Compute estimates of enrichment scores without gene expression profiles
	data. A set of interesting genes (my_gene_labels) is provided instead.

	Args:
		my_gene_labels (set): A set of gene labels appearing near the top or
			near the bottom of the list with ranked genes.
		gene_sets (list): Sets of labels of genes to be used in the
			calculation of enrichment scores.
		p (float): Exponent used in the calculation of the running sum (see
			the original paper).

	Returns:
		list: Estimates of enrichment score values.
	"""
	return [get_enrichment_score(my_gene_labels, gene_set, p) \
		for gene_set in gene_sets]

def my_gsea(my_genes_fn, gene_sets_fn, p=1.0):
	"""
	Compute estimates of normalized enrichment scores and the corresponding
	p-values. The results are sorted by normalized enrichment score and printed
	out in three columns: (1) gene set label, (2) normalized_enrichment_score,
	and (3) statistical significance (p-value).

	Args:
		my_genes_fn (str): Name of the file with a set of interesting genes.
		gene_sets_fn (str): Name of the file with a list of gene sets.
		p (float): Exponent used in the calculation of the running sum (see
			the original paper).

	Returns:
		None
	"""
	my_gene_labels = get_my_genes(my_genes_fn)
	gene_sets, gene_set_labels = get_gene_sets(gene_sets_fn)
	normalized_enrichment_scores = []
	significances = []
	for gene_set in gene_sets:
		common_gene_labels = my_gene_labels.intersection(gene_set)
		q = float(len(common_gene_labels)) / len(gene_set)
		if q == 0.0:
			normalized_enrichment_scores.append(1.0)
			significances.append(0.5)
		if q < 1.0:
			normalized_enrichment_scores.append(2.0)
			significances.append(0.0)
		if q == 1.0:
			normalized_enrichment_scores.append(4.0 / 3.0)
			significances.append(0.0)
	results = sorted(zip(gene_set_labels, normalized_enrichment_scores, \
		significances), key=lambda item: item[1], reverse=True)
	print('gene_set_label', 'normalized_enrichment_score', \
		'statistical_significance', sep='\t')
	for result in results:
		print(result[0], '%.3f' % result[1], '%.3f' % result[2], sep='\t')

numpy.random.seed(2017)

#gsea('leukemia.txt', 'pathways.txt')

# Fig. 1 in the original paper for MAP00480_Glutathione_metabolism, one of the
# highest scoring gene sets
#plot('leukemia.txt', set(['GSTM2', 'GPX3', 'GSTT1', 'MGST1', 'GCLC', 'GSTP1', \
#	'GSS', 'GSR', 'IDH2', 'GPX4', 'GGT1', 'GPX1', 'GSTM1', 'G6PD', 'IDH1', \
#	'GSTA2', 'GSTM5', 'GSTM3', 'GSTT2', 'GPX2', 'GSTM4',	'GCLM']))

#my_gsea('my_genes.txt', 'pathways.txt')
