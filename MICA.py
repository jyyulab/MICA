#####################################################################################################################
# Copyright (c) 2017, St. Jude Children's Research Hospital															#
# All rights reserved.																								#
#																													#
# Redistribution and use in source and binary forms, with or without												#
# modification, are permitted provided that the following conditions are met:										#
#																													#
# 1. Redistributions of source code must retain the above copyright notice, this 									#
#    list of conditions and the following disclaimer.																#
#																													#
# 2. Redistributions in binary form must reproduce the above copyright notice,										#
#    this list of conditions and the following disclaimer in the documentation										#
#    and/or other materials provided with the distribution.															#
#																													#
# 3. Neither the name of the copyright holder nor the names of its contributors 									#
#    may be used to endorse or promote products derived from this software 											#
#    without specific prior written permission.																		#
#																													#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 									#
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED										#
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 											#
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 										#
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 										#
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 										#
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 										#
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 									#
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 									#
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 												#
#####################################################################################################################
# MICA: Mutual Information-Based Clustering Analysis 																#
# This program is provided for single-cell RNA-seq analysis, in which a new approach 								#
# has been developed for cell type identification and corresponding analysis.										#
# Additional options like overlaying gene expressions, and measuring the quality of 								#
# the results are included in the package.																			#
# Author: Alireza Khatamian (akhatami@stjude.org) 																	#
# Version: 1.0.0																									#
#####################################################################################################################
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['font.family'] = 'Arial'
import sys, os, math, csv, numpy, time, pandas
import parameter as parm
import matplotlib.pyplot as plt
from sklearn import manifold 														# TSNE, Laplacian
from sklearn import cluster 														# KMeans 
from sklearn import decomposition 													# PCA
from sklearn.cluster import AgglomerativeClustering									# Hierarchical clustering
from scipy.cluster.hierarchy import dendrogram, linkage 							# Hierarchical clustering visualization
import matplotlib.gridspec as gridspec 												# Grid panel for figures
from sklearn.metrics.cluster import adjusted_rand_score 							# Clustering quality measurement
from sklearn.metrics.cluster import normalized_mutual_info_score 					# Clustering quality measurement
from sklearn.metrics import silhouette_samples, silhouette_score					# Optimal clustering score for K
from scipy import stats 
from threading import Thread, Lock															

fig_num = 1
mutex = Lock()
mutex1 = Lock()
mutex2 = Lock()

# Plotting three different graphs with respect to gene analysis
# 1. Scatter plot, overlaying the gene expression on the clustering results to show the cluster with a specific highly expressed gene
# 2. Violin plot, cross cluster comparison of a specific gene to show the cluster(s) with high expression
# 3. Bar plot, P-value and t-test related plot for cross-cluster comparison
def overlay(parameter, _object):
	print '[INFO] --> [OVLY] Overlaying gene expressions on cells ...'

	global fig_num

	path = '/'.join(parameter.infile['path'].split('/')[:-1]) + '/'
	result_path = path + 'Results_' + parameter.job_name + '/'

	pred_labels_dict = _object['mica_labels']
	biomarkers = _object['biomarkers']
	expression = _object['expression']
	expression_genes = list(expression['header'])
	expression_matrix = expression['matrix']
	expression_cells = list(expression['cells'])

	n = int(numpy.ceil(numpy.sqrt(len(biomarkers) + 1)))
	n1 = numpy.sqrt(n)
	n_cell = _object['n_cell']

	for item in pred_labels_dict:
		fig = plt.figure(fig_num, figsize=(16, 16), dpi=100)
		fig_num += 1
		fig1 = plt.figure(fig_num, figsize=(16, 16), dpi=100)
		fig_num += 1
		fig2 = plt.figure(fig_num, figsize=(16, 16), dpi=100)
		fig_num += 1

		pred_labels = []
		actual_labels = []
		for k, v in item.iteritems():
			actual_labels.append(k)
			pred_labels.append(v)

		i = numpy.max(pred_labels) + 1
		class_size = [b[1] for b in sorted(enumerate([pred_labels.count(z) for z in range(i)]), key=lambda z:z[1], reverse=True)]
		trans = _object['trans']['Y']
		trans_labels = _object['trans']['labels']
		tsne = manifold.TSNE(n_components=2, n_iter=5000, learning_rate=200, perplexity=numpy.min([numpy.max(class_size),parm.MAX_PERPLEXITY]), random_state=10)
		trans_tsne = tsne.fit_transform(numpy.array(trans[:,0:max(_object['dim_use'])]))		

		f = open(result_path + 'Overlay_Biomarkers_k' + str(i) + '.txt', 'w')
		f1 = open(result_path + 'Biomarkers_Expression_k' + str(i) + '.txt', 'w')
		f2 = open(result_path + 'Biomarkers_Cluster_PValue_k' + str(i) + '.txt', 'w')
		f.write('id\t' + '\t'.join(trans_labels) + '\n')

		grid = gridspec.GridSpec(n, n, wspace=0.2, hspace=0.4)
		grid1 = gridspec.GridSpec(n, n, wspace=0.2, hspace=0.4)
		grid2 = gridspec.GridSpec(n, n, wspace=0.2, hspace=0.4)

		subgrid = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=grid[0], wspace=0.1, hspace=0.1)
		subgrid1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=grid1[0], wspace=0.6, hspace=0.1)
		subgrid2 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=grid2[0], wspace=0.6, hspace=0.1)

		ax = plt.Subplot(fig, subgrid[0])
		ax1 = plt.Subplot(fig1, subgrid1[0])
		ax2 = plt.Subplot(fig2, subgrid2[0])
		marker_size = parm.MAXIMUM_SPECIAL_MARKER_SIZE if n_cell < parameter.MARKER_SIZE_THRESHOLD else parm.MINIMUM_SPECIAL_MARKER_SIZE
		marker_size /= n1
		colors = plt.cm.jet(numpy.linspace(0, 1, i+1))
		centers = []

		for z in range(i):
			pc = [p for p in range(n_cell) if pred_labels[actual_labels.index(trans_labels[p])] == z]
			xc = [trans_tsne[p, 0] for p in pc]
			yc = [trans_tsne[p, 1] for p in pc]
			ax.scatter(xc, yc, facecolor=colors[z+1], s=marker_size, marker='o', vmin=0, vmax=i, label=str(z+1), alpha=0.7)
			ax1.scatter(xc, yc, facecolor=colors[z+1], s=marker_size, marker='o', vmin=0, vmax=i, label=str(z+1), alpha=0.7)
			ax2.scatter(xc, yc, facecolor=colors[z+1], s=marker_size, marker='o', vmin=0, vmax=i, label=str(z+1), alpha=0.7)
			centers.append([numpy.sum(xc) / float(len(xc)), numpy.sum(yc) / float(len(yc))])

		centers = numpy.array(centers)
		for z, c in enumerate(centers):
			ax.scatter(c[0], c[1], marker='$%d$' % (z+1), alpha=0.8, s=marker_size*2, edgecolor='k')
			ax1.scatter(c[0], c[1], marker='$%d$' % (z+1), alpha=0.8, s=marker_size*2, edgecolor='k')
			ax2.scatter(c[0], c[1], marker='$%d$' % (z+1), alpha=0.8, s=marker_size*2, edgecolor='k')
		ax.set_ylabel('MICA-2', fontsize=parm.AXIS_LABEL_FONTSIZE/n1)	
		ax.set_xlabel('MICA-1', fontsize=parm.AXIS_LABEL_FONTSIZE/n1)
		ax1.set_ylabel('MICA-2', fontsize=parm.AXIS_LABEL_FONTSIZE/n1)
		ax1.set_xlabel('MICA-1', fontsize=parm.AXIS_LABEL_FONTSIZE/n1)
		ax2.set_ylabel('MICA-2', fontsize=parm.AXIS_LABEL_FONTSIZE/n1)
		ax2.set_xlabel('MICA-1', fontsize=parm.AXIS_LABEL_FONTSIZE/n1)
		ax.set_xticks([])
		ax.set_yticks([])
		ax1.set_xticks([])
		ax1.set_yticks([])
		ax2.set_xticks([])
		ax2.set_yticks([])
		fig.add_subplot(ax)
		fig1.add_subplot(ax1)
		fig2.add_subplot(ax2)

		t = 1
		for bm in biomarkers:
			subgrid = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=grid[t], wspace=0.1, hspace=0.1)
			subgrid1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=grid1[t], wspace=0.6, hspace=0.1)
			subgrid2 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=grid2[t], wspace=0.6, hspace=0.1)

			ax = plt.Subplot(fig, subgrid[0])
			if bm not in expression_genes:
				print '[WARN] --> [OVLY] Gene ' + str(bm) + ' not found, ignored.'
				continue
			gene_index = expression_genes.index(bm)
			gene_expression = [expression_matrix[expression_cells.index(trans_labels[j])][gene_index] for j in range(n_cell)]
			ax.scatter(trans_tsne[:,0], trans_tsne[:,1], s=marker_size, cmap=plt.cm.jet, c=gene_expression, marker='o')
			ax.set_title(bm, fontsize=parm.SPECIAL_TITLE_FONTSIZE/n1)
			ax.set_ylabel('MICA-2', fontsize=parm.AXIS_LABEL_FONTSIZE/n1)	
			ax.set_xlabel('MICA-1',fontsize=parm.AXIS_LABEL_FONTSIZE/n1)
			ax.set_xticks([])
			ax.set_yticks([])
			fig.add_subplot(ax)
			f.write(bm + '\t' + '\t'.join([str(p) for p in gene_expression]) + '\n')

			ax1 = plt.Subplot(fig1, subgrid1[0])
			gene_expression_class_wise = [[]] * i
			gene_expression_class_wise_complement = [[]] * i

			for k in range(i):
				exp = [gene_expression[p] for p in range(n_cell) if pred_labels[actual_labels.index(trans_labels[p])] == k]
				exp_complement = [gene_expression[p] for p in range(n_cell) if pred_labels[actual_labels.index(trans_labels[p])] != k]
				gene_expression_class_wise[k] = exp
				gene_expression_class_wise_complement[k] = exp_complement
				f1.write(bm + '\t' + '\t'.join([trans_labels[p] for p in range(n_cell) if pred_labels[actual_labels.index(trans_labels[p])] == k]) + '\n')
				f1.write(str(k) + '\t' + '\t'.join(str(p) for o in exp) + '\n')
				
				jitter = numpy.random.rand(len(exp)) - 0.5
				density = 1.0
				if not all(v == 0 for v in exp):
					density = stats.gaussian_kde(exp)(exp)
				xc = (k+1) + (density * jitter * 0.4)
				ax1.scatter(xc, exp, s=0.5, marker='.', c='k')
				

			#ax1.boxplot(gene_expression_class_wise, flierprops=dict(marker='o', markersize=4))
			parts = ax1.violinplot(gene_expression_class_wise, showmedians=True, showmeans=False, showextrema=False)
			for pc in parts['bodies']:
				pc.set_facecolor('none')
				pc.set_edgecolor(colors[parts['bodies'].index(pc) + 1])
			ax1.set_xticks([p+1 for p in range(i)])
			ax1.set_xticklabels([str(p+1) for p in range(i)], fontsize=parm.TICK_LABEL_FONTSIZE/n1)
			ax1.tick_params(axis='y', which='major', labelsize=parm.TICK_LABEL_FONTSIZE/n1)
			ax1.set_title(bm, fontsize=parm.SPECIAL_TITLE_FONTSIZE/n1)
			ax1.set_ylabel('Gene Expression', fontsize=parm.AXIS_LABEL_FONTSIZE/n1)
			ax1.set_xlabel('Clusters', fontsize=parm.AXIS_LABEL_FONTSIZE/n1)
			fig1.add_subplot(ax1)

			xd = [p+1 for p in range(i)]
			xd_label = [str(p+1) for p in range(i)]
			yd = [stats.ttest_ind(numpy.array(gene_expression_class_wise[p]), numpy.array(gene_expression_class_wise_complement[p]), equal_var=False) for p in range(i)]
			yd = [p[1] if p[0] > 0 else 1 for p in yd]
			yd_threshold = numpy.min([p for p in yd if p > 0])
			yd = [p if p >= yd_threshold else yd_threshold for p in yd]
			yd = [-numpy.log10(p) for p in yd]
			ax21 = plt.Subplot(fig2, subgrid2[0])
			ax21.bar(xd, yd, align='center', color=colors[1:])
			ax21.set_xticks(xd)
			ax21.set_xticklabels(xd_label)
			ax21.tick_params(axis='y', which='major', labelsize=parm.TICK_LABEL_FONTSIZE/n1)
			ax21.set_title(bm, fontsize=parm.SPECIAL_TITLE_FONTSIZE/n1)
			ax21.set_ylabel('-log10(p-value)', fontsize=parm.AXIS_LABEL_FONTSIZE/n1)
			ax21.set_xlabel('Clusters', fontsize=parm.AXIS_LABEL_FONTSIZE/n1)
			fig2.add_subplot(ax21)

			f2.write('Cluster\t' + '\t'.join(xd_label) + '\n')
			f2.write(bm + '\t' + '\t'.join([str(p) for p in yd]) + '\n')

			t += 1
		fig.savefig(result_path + 'Overlay_Biomarkers_k' + str(i) + '.png', bbox_inches='tight')
		fig.savefig(result_path + 'Overlay_Biomarkers_k' + str(i) + '.eps', bbox_inches='tight')
		fig.savefig(result_path + 'Overlay_Biomarkers_k' + str(i) + '.pdf', bbox_inches='tight')

		fig1.savefig(result_path + 'Biomarkers_Expression_k' + str(i) + '.png', bbox_inches='tight')
		fig1.savefig(result_path + 'Biomarkers_Expression_k' + str(i) + '.eps', bbox_inches='tight')
		fig1.savefig(result_path + 'Biomarkers_Expression_k' + str(i) + '.pdf', bbox_inches='tight')

		fig2.savefig(result_path + 'Biomarkers_Cluster_PValue_k' + str(i) + '.png', bbox_inches='tight')
		fig2.savefig(result_path + 'Biomarkers_Cluster_PValue_k' + str(i) + '.eps', bbox_inches='tight')
		fig2.savefig(result_path + 'Biomarkers_Cluster_PValue_k' + str(i) + '.pdf', bbox_inches='tight')

		f.close()
		f1.close()
		f2.close()

	plt.close()
			
	print '[INFO] --> [OVLY] Done.'

# Estimating an optimal K for clustering based on silhouette score of consensus clustering results over multiple given Ks
def optimal_k_2(parameter, _object):
	print '[INFO] --> [OPTK] Estimating optimal K ...'

	global fig_num

	path = '/'.join(parameter.infile['path'].split('/')[:-1]) + '/'
	result_path = path + 'Results_' + parameter.job_name + '/'

	kn = _object['kn']
	n_cell = _object['n_cell']
	n = int(len(kn))
	silhouette_below_avg_count = [0] * n
	cluster_size_below_avg_count = [0] * n
	silhouette_avg_list = [0] * n

	fig = plt.figure(fig_num, figsize=(16, 8 * len(kn)), dpi=300)
	fig_num += 1
	hratios = [20] * (n+1)
	hratios[0] = 1
	grid = gridspec.GridSpec(n+1, 1, wspace=0.1, hspace=0.2, height_ratios=hratios)
	t = 1
	for k in kn:
		trans = _object['trans']['Y'][:,0:numpy.max(_object['dim_use'])]
		con_result = _object['con_result'][str(k)]
		km_result = numpy.array([con_result['sorted_labels'][con_result['consensus_labels'].index(p)] for p in _object['trans']['labels']])
		silhouette_avg = silhouette_score(trans, km_result)
		silhouette_avg_list[kn.index(k)] = silhouette_avg
		sample_silhouette_values = silhouette_samples(trans, km_result)
		y_lower = 10
		
		subgrid = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=grid[t], wspace=0.1, hspace=0.4)
		ax = plt.Subplot(fig, subgrid[0])
		ax1 = plt.Subplot(fig, subgrid[1])
		ax.set_xlim([-0.1, 1])
		ax.set_ylim([0, len(trans) + (k + 1) * 10])

		for i in range(k):	
			ith_cluster_silhouette_values = sample_silhouette_values[km_result == i]
			ith_cluster_silhouette_values.sort()
			size_cluster_i = ith_cluster_silhouette_values.shape[0]
			silhouette_below_avg_count[kn.index(k)] += 1 if numpy.max(ith_cluster_silhouette_values) < silhouette_avg else 0
			cluster_size_below_avg_count[kn.index(k)] += 1 if size_cluster_i < float(n_cell) / k else 0
			y_upper = y_lower + size_cluster_i
			color = plt.cm.jet(float(i) / k)
			ax.fill_betweenx(numpy.arange(y_lower, y_upper), 0, ith_cluster_silhouette_values, facecolor=color, edgecolor=color, alpha=0.7)
			ax.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i), fontsize=parm.TEXT_FONTSIZE)
			y_lower = y_upper + 10
		ax.set_title('Silhouette Score', fontsize=parm.TITLE_FONTSIZE)
		ax.set_xlabel('Silhouette coefficients', fontsize=parm.AXIS_LABEL_FONTSIZE)
		ax.set_ylabel('Clusters', fontsize=parm.AXIS_LABEL_FONTSIZE)
		ax.axvline(x=silhouette_avg, color='r', linestyle='--')
		ax.set_yticks([])
		ax.set_xticks(numpy.arange(-0.2, 1.1, 0.2))

		colors = plt.cm.jet(km_result.astype(float) / k)
		marker_size = parm.MAXIMUM_OPT_MARKER_SIZE if len(km_result) < parameter.MARKER_SIZE_THRESHOLD else parm.MINIMUM_OPT_MARKER_SIZE
		tsne = None
		trans_tsne = None
		if k not in _object['tsne_trans']:
			class_size = [b[1] for b in sorted(enumerate([km_result.tolist().count(z) for z in range(k)]), key=lambda z:z[1], reverse=True)]
			tsne = manifold.TSNE(n_components=2, n_iter=5000, learning_rate=200, perplexity=numpy.min([numpy.max(class_size),parm.MAX_PERPLEXITY]), random_state=10)
			trans_tsne = tsne.fit_transform(trans)
			_object['tsne_trans'][k] = trans_tsne
		else:
			trans_tsne = _object['tsne_trans'][k]
		ax1.scatter(trans_tsne[:,0], trans_tsne[:,1], marker='o', s=marker_size, c=colors, edgecolor='k', alpha=0.7)
		centers = []
		for l in numpy.unique(km_result):
			xc = [trans_tsne[p, 0] for p in range(n_cell) if km_result[p] == l]
			yc = [trans_tsne[p, 1] for p in range(n_cell) if km_result[p] == l]
			centers.append([numpy.sum(xc) / float(len(xc)), numpy.sum(yc) / float(len(yc))])
		centers = numpy.array(centers)
		ax1.scatter(centers[:,0], centers[:,1], marker='o', c='white', alpha=0.5, s=marker_size*3, edgecolor='k')
		for i, c in enumerate(centers):
			ax1.scatter(c[0], c[1], marker='$%d$' % i, alpha=0.7, s=marker_size*2, edgecolor='k')
		ax1.set_title('Clustered Data.', fontsize=parm.TITLE_FONTSIZE)
		ax1.set_xlabel('1st feature', fontsize=parm.AXIS_LABEL_FONTSIZE)
		ax1.set_ylabel('2nd feature', fontsize=parm.AXIS_LABEL_FONTSIZE)
		ax1.set_xticks([])
		ax1.set_yticks([])
		fig.add_subplot(ax)
		fig.add_subplot(ax1)
		t += 1
	score_filter = [i for i, s in enumerate(silhouette_below_avg_count) if s == numpy.min(silhouette_below_avg_count)]
	#size_filter = [i for i, s in enumerate(cluster_size_below_avg_count) if i in score_filter]
	#silhouette_filter = [i for i in size_filter if cluster_size_below_avg_count[i] == numpy.min(cluster_size_below_avg_count)]
	#silhouette_filter = silhouette_filter if len(silhouette_filter) > 0 else size_filter
	#silhouette_filter = silhouette_filter if len(silhouette_filter) > 0 else score_filter
	silhouette_filter = score_filter
	optimal_k = [kn[i] for i in silhouette_filter]
	subgrid = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=grid[0], wspace=0.1, hspace =0.0)
	print '[INFO] --> [OPTK] Optimal K for clustering is: ' + str(optimal_k)
	optimal_k_2 = None
	if len(optimal_k) > 1:
		silhouette_filter = [silhouette_avg_list.index(numpy.max(numpy.array(silhouette_avg_list)[numpy.array(silhouette_filter)]))]
		optimal_k_2 = [kn[i] for i in silhouette_filter]
		print '[INFO] --> [OPTK] Optimal K with highest average score for clustering is: ' + str(optimal_k_2)
	ax = plt.Subplot(fig, subgrid[0])
	if optimal_k_2 == None:
		ax.set_title(('Silhouette Analysis with optimal K = ' + str(optimal_k)), fontsize=parm.HEADER_FONTSIZE, fontweight='bold')
	else:
		ax.set_title(('Silhouette Analysis with optimal K = ' + str(optimal_k) + '\nRecommended optimal K = ' + str(optimal_k_2)), fontsize=parm.HEADER_FONTSIZE, fontweight='bold')
	ax.axis('off')
	fig.add_subplot(ax)
	
	fig.savefig(result_path + 'Optimal_K_Analysis_2.png', bbox_inches='tight')
	fig.savefig(result_path + 'Optimal_K_Analysis_2.eps', bbox_inches='tight')
	fig.savefig(result_path + 'Optimal_K_Analysis_2.pdf', bbox_inches='tight')
	plt.close()
	print '[INFO] --> [OPTK] Done.'
	return optimal_k_2 if optimal_k_2 != None else optimal_k

# Estimating an optimal K for clustering based on silhouette score of simple K-Means results over multiple given Ks
def optimal_k(parameter, _object):
	print '[INFO] --> [OPTK] Estimating optimal K ...'

	global fig_num

	path = '/'.join(parameter.infile['path'].split('/')[:-1]) + '/'
	result_path = path + 'Results_' + parameter.job_name + '/'

	kn = _object['kn']
	n_cell = _object['n_cell']
	n = int(len(kn))
	silhouette_below_avg_count = [0] * n
	cluster_size_below_avg_count = [0] * n
	silhouette_avg_list = [0] * n

	fig = plt.figure(fig_num, figsize=(16, 8 * len(kn)), dpi=300)
	fig_num += 1
	hratios = [20] * (n+1)
	hratios[0] = 1
	grid = gridspec.GridSpec(n+1, 1, wspace=0.1, hspace=0.2, height_ratios=hratios)
	t = 1
	for k in kn:
		trans = _object['trans']['Y'][:,0:numpy.max(_object['dim_use'])]
		km = cluster.KMeans(n_clusters=k, random_state=10)
		km_result = km.fit_predict(trans)
		silhouette_avg = silhouette_score(trans, km_result)
		silhouette_avg_list[kn.index(k)] = silhouette_avg
		sample_silhouette_values = silhouette_samples(trans, km_result)
		y_lower = 10
		
		subgrid = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=grid[t], wspace=0.1, hspace=0.4)
		ax = plt.Subplot(fig, subgrid[0])
		ax1 = plt.Subplot(fig, subgrid[1])
		ax.set_xlim([-0.1, 1])
		ax.set_ylim([0, len(trans) + (k + 1) * 10])

		for i in range(k):	
			ith_cluster_silhouette_values = sample_silhouette_values[km_result == i]
			ith_cluster_silhouette_values.sort()
			size_cluster_i = ith_cluster_silhouette_values.shape[0]
			silhouette_below_avg_count[kn.index(k)] += 1 if numpy.max(ith_cluster_silhouette_values) < silhouette_avg else 0
			cluster_size_below_avg_count[kn.index(k)] += 1 if size_cluster_i < float(n_cell) / k else 0
			y_upper = y_lower + size_cluster_i
			color = plt.cm.jet(float(i) / k)
			ax.fill_betweenx(numpy.arange(y_lower, y_upper), 0, ith_cluster_silhouette_values, facecolor=color, edgecolor=color, alpha=0.7)
			ax.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i), fontsize=parm.TEXT_FONTSIZE)
			y_lower = y_upper + 10
		ax.set_title('Silhouette Score', fontsize=parm.TITLE_FONTSIZE)
		ax.set_xlabel('Silhouette coefficients', fontsize=parm.AXIS_LABEL_FONTSIZE)
		ax.set_ylabel('Clusters', fontsize=parm.AXIS_LABEL_FONTSIZE)
		ax.axvline(x=silhouette_avg, color='r', linestyle='--')
		ax.set_yticks([])
		ax.set_xticks(numpy.arange(-0.2, 1.1, 0.2))

		colors = plt.cm.jet(km_result.astype(float) / k)
		marker_size = parm.MAXIMUM_OPT_MARKER_SIZE if len(km_result) < parameter.MARKER_SIZE_THRESHOLD else parm.MINIMUM_OPT_MARKER_SIZE
		class_size = [b[1] for b in sorted(enumerate([km_result.tolist().count(z) for z in range(k)]), key=lambda z:z[1], reverse=True)]
		tsne = manifold.TSNE(n_components=2, n_iter=5000, learning_rate=200, perplexity=numpy.min([numpy.max(class_size),parm.MAX_PERPLEXITY]), random_state=10)
		trans_tsne = tsne.fit_transform(numpy.array(trans[:,0:max(_object['dim_use'])]))
		ax1.scatter(trans_tsne[:,0], trans_tsne[:,1], marker='o', s=marker_size, c=colors, edgecolor='k', alpha=0.7)
		centers = []
		for l in numpy.unique(km.labels_):
			xc = [trans_tsne[p, 0] for p in range(n_cell) if km.labels_[p] == l]
			yc = [trans_tsne[p, 1] for p in range(n_cell) if km.labels_[p] == l]
			centers.append([numpy.sum(xc) / float(len(xc)), numpy.sum(yc) / float(len(yc))])
		centers = numpy.array(centers)
		ax1.scatter(centers[:,0], centers[:,1], marker='o', c='white', alpha=0.5, s=marker_size * 3, edgecolor='k')
		for i, c in enumerate(centers):
			ax1.scatter(c[0], c[1], marker='$%d$' % i, alpha=0.7, s=marker_size*2, edgecolor='k')
		ax1.set_title('Clustered Data', fontsize=parm.TITLE_FONTSIZE)
		ax1.set_xlabel('1st feature', fontsize=parm.AXIS_LABEL_FONTSIZE)
		ax1.set_ylabel('2nd feature', fontsize=parm.AXIS_LABEL_FONTSIZE)
		ax1.set_xticks([])
		ax1.set_yticks([])
		fig.add_subplot(ax)
		fig.add_subplot(ax1)
		t += 1
	score_filter = [i for i, s in enumerate(silhouette_below_avg_count) if s == numpy.min(silhouette_below_avg_count)]
	#size_filter = [i for i, s in enumerate(cluster_size_below_avg_count) if i in score_filter]
	#silhouette_filter = [i for i in size_filter if cluster_size_below_avg_count[i] == numpy.min(cluster_size_below_avg_count)]
	#silhouette_filter = silhouette_filter if len(silhouette_filter) > 0 else size_filter
	#silhouette_filter = silhouette_filter if len(silhouette_filter) > 0 else score_filter
	silhouette_filter = score_filter
	optimal_k = [kn[i] for i in silhouette_filter]
	subgrid = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=grid[0], wspace=0.1, hspace =0.0)
	print '[INFO] --> [OPTK] Optimal K for clustering is: ' + str(optimal_k)
	optimal_k_2 = None
	if len(optimal_k) > 1:
		silhouette_filter = [silhouette_avg_list.index(numpy.max(numpy.array(silhouette_avg_list)[numpy.array(silhouette_filter)]))]
		optimal_k_2 = [kn[i] for i in silhouette_filter]
		print '[INFO] --> [OPTK] Optimal K with highest average score for clustering is: ' + str(optimal_k_2)
	ax = plt.Subplot(fig, subgrid[0])
	if optimal_k_2 == None:
		ax.set_title(('Silhouette Analysis with optimal K = ' + str(optimal_k)), fontsize=parm.HEADER_FONTSIZE, fontweight='bold')
	else:
		ax.set_title(('Silhouette Analysis with optimal K = ' + str(optimal_k) + '\nRecommended optimal K = ' + str(optimal_k_2)), fontsize=parm.HEADER_FONTSIZE, fontweight='bold')
	ax.axis('off')
	fig.add_subplot(ax)
	fig.savefig(result_path + 'Optimal_K_Analysis.png', bbox_inches='tight')
	fig.savefig(result_path + 'Optimal_K_Analysis.eps', bbox_inches='tight')
	fig.savefig(result_path + 'Optimal_K_Analysis.pdf', bbox_inches='tight')
	plt.close()
	print '[INFO] --> [OPTK] Done.'
	return optimal_k_2 if optimal_k_2 != None else optimal_k

# Plotting box chart with respect to c type quality measurement, this plot shows the stability of the method over the various dimensions
def cluster_cquality_plot(cquality, n_components):
	print '[INFO] --> [PLOT] Plotting c-quality ...'

	global fig_num

	path = '/'.join(parameter.infile['path'].split('/')[:-1]) + '/'
	result_path = path + 'Results_' + parameter.job_name + '/'

	purity = []
	ari = []
	nmi = []

	cross_component_mean_purity = []
	cross_component_mean_ari = []
	cross_component_mean_nmi = []
	for cc in cquality:
		component_purity = [x[2] for x in cc]
		component_ari = [x[0] for x in cc]
		component_nmi = [x[1] for x in cc]
		cross_component_mean_purity.append(numpy.mean(component_purity))
		cross_component_mean_ari.append(numpy.mean(component_ari))
		cross_component_mean_nmi.append(numpy.mean(component_nmi))
		purity.append(component_purity)
		ari.append(component_ari)
		nmi.append(component_nmi)

	fig = plt.figure(fig_num, figsize=(16, 16), dpi=300)
	fig_num += 1
	grid = gridspec.GridSpec(1, 1, wspace=0.2, hspace=0.2)
	subgrid = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=grid[0], wspace=0.1, hspace=0.1)
	metrics = ('ARI', 'NMI', 'Purity')
	x = numpy.arange(len(n_components)) + 1
	colors = plt.cm.Accent(numpy.linspace(0, 1, len(metrics)))
	ax = plt.Subplot(fig, subgrid[0])
	ax1 = plt.Subplot(fig, subgrid[1])
	ax2 = plt.Subplot(fig, subgrid[2])
	ax.plot(x, cross_component_mean_purity, marker='.', color=colors[2], label=metrics[2], ms=0, alpha=0.7)
	ax1.plot(x, cross_component_mean_ari, marker='.', color=colors[0], label=metrics[0], ms=0, alpha=0.7)
	ax2.plot(x, cross_component_mean_nmi, marker='.', color=colors[1], label=metrics[1], ms=0, alpha=0.7)

	purity_data = [[]] * len(x)
	ari_data = [[]] * len(x)
	nmi_data = [[]] * len(x)
	for c, xx in enumerate(x):
		purity_data[c] = numpy.concatenate((purity[c], [cross_component_mean_purity[c]], [p for p in purity[c] if p >= cross_component_mean_purity[c]], [p for p in purity[c] if p <= cross_component_mean_purity[c]]))
		ari_data[c] = numpy.concatenate((ari[c], [cross_component_mean_ari[c]], [p for p in ari[c] if p >= cross_component_mean_ari[c]], [p for p in ari[c] if p <= cross_component_mean_ari[c]]))
		nmi_data[c] = numpy.concatenate((nmi[c], [cross_component_mean_nmi[c]], [p for p in nmi[c] if p >= cross_component_mean_nmi[c]], [p for p in nmi[c] if p <= cross_component_mean_nmi[c]]))
	
	ax.boxplot(purity_data, boxprops=dict(color=colors[2]), flierprops=dict(marker='o', markerfacecolor=colors[2], markersize=4), medianprops=dict(color=colors[2]))
	ax1.boxplot(ari_data, boxprops=dict(color=colors[0]), flierprops=dict(marker='o', markerfacecolor=colors[0], markersize=4), medianprops=dict(color=colors[0]))
	ax2.boxplot(nmi_data, boxprops=dict(color=colors[1]), flierprops=dict(marker='o', markerfacecolor=colors[1], markersize=4), medianprops=dict(color=colors[1]))
	plt.setp(ax.get_xticklabels(), visible=False)
	plt.setp(ax1.get_xticklabels(), visible=False)
	ax.set_xticks([p for p in x if p % 5 == 0])
	ax1.set_xticks([p for p in x if p % 5 == 0])
	ax2.set_xticks([p for p in x if p % 5 == 0])
	ax2.set_xticklabels([p for p in n_components if p % 5 == 0], fontsize=parm.TICK_LABEL_FONTSIZE, rotation=40, ha='right')
	ax.set_yticks(numpy.arange(0, 1.02, 0.2))
	ax1.set_yticks(numpy.arange(0, 1.02, 0.2))
	ax2.set_yticks(numpy.arange(0, 1.02, 0.2))
	ax1.set_ylabel('Quality', fontsize=parm.AXIS_LABEL_FONTSIZE)
	ax2.set_xlabel('#Components', fontsize=parm.AXIS_LABEL_FONTSIZE)
	ax.set_title('Clustering Quality Components (' + str(len(cquality)) + ')', fontsize=parm.TITLE_FONTSIZE)

	h, l = ax.get_legend_handles_labels()
	h1, l1 = ax1.get_legend_handles_labels()
	h2, l2 = ax2.get_legend_handles_labels()
	ax1.legend(h+h1+h2, l+l1+l2, loc='center left', bbox_to_anchor=(1, 0.5))
	fig.add_subplot(ax, sahrex=ax2, sharey=ax1)
	fig.add_subplot(ax1, sharex=ax2)
	fig.add_subplot(ax2, sharey=ax1)

	fig.savefig(result_path + 'Clustering_C_Quality.png', bbox_inches='tight')
	fig.savefig(result_path + 'Clustering_C_Quality.eps', bbox_inches='tight')
	fig.savefig(result_path + 'Clustering_C_Quality.pdf', bbox_inches='tight')

	with open(result_path + 'Clustering_C_Quality.txt', 'w') as f:
		f.write('#Components\t' + '\t'.join([str(b) for b in n_components]) + '\n')
		f.write('Mean Purity\t' + '\t'.join([str(b) for b in cross_component_mean_purity]) + '\n')
		f.write('Mean ARI\t' + '\t'.join([str(b) for b in cross_component_mean_ari]) + '\n')
		f.write('Mean NMI\t' + '\t'.join([str(b) for b in cross_component_mean_nmi]) + '\n')
		purity_trans = numpy.transpose(numpy.array(purity))
		ari_trans = numpy.transpose(numpy.array(ari))
		nmi_trans = numpy.transpose(numpy.array(nmi))
		for i in range(purity_trans.shape[0]):
			f.write('Purity #' + str(i) + '\t' + '\t'.join([str(b) for b in purity_trans[i].tolist()]) + '\n')
			f.write('ARI #' + str(i) + '\t' + '\t'.join([str(b) for b in ari_trans[i].tolist()]) + '\n')
			f.write('NMI #' + str(i) + '\t' + '\t'.join([str(b) for b in nmi_trans[i].tolist()]) + '\n')

	plt.close()
	print '[INFO] --> [PLOT] Done.'

# Plotting bar chart with respect to b type quality measurement, this plot shows the confidence factor of the method
def cluster_bquality_plot(bquality):
	print '[INFO] --> [PLOT] Plotting b-quality ...'

	global fig_num

	path = '/'.join(parameter.infile['path'].split('/')[:-1]) + '/'
	result_path = path + 'Results_' + parameter.job_name + '/'

	purity = [x[2] for x in bquality]
	ari = [x[0] for x in bquality]
	nmi = [x[1] for x in bquality]

	fig = plt.figure(fig_num, figsize=(16, 16), dpi=300)
	fig_num += 1
	grid = gridspec.GridSpec(1, 1, wspace=0.2, hspace=0.2)
	subgrid = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=grid[0], wspace=0.1, hspace=0.1)
	metrics = ('ARI', 'NMI', 'Purity')
	x = numpy.arange(len(metrics))
	y = [0] * len(metrics)
	yerr_high = [0] * len(metrics)
	yerr_low = [0] * len(metrics)
	y[2] = numpy.mean(purity)
	yerr_high[2] = numpy.max(purity) - y[2]
	yerr_low[2] = -numpy.min(purity) + y[2]
	y[0] = numpy.mean(ari)
	yerr_high[0] = numpy.max(ari) - y[0]
	yerr_low[0] = -numpy.min(ari) + y[0]
	y[1] = numpy.mean(nmi)
	yerr_high[1] = numpy.max(nmi) - y[1]
	yerr_low[1] = -numpy.min(nmi) + y[1]
	colors = plt.cm.Accent(numpy.linspace(0, 1, len(metrics)))
	ax = plt.Subplot(fig, subgrid[0])
	ax.bar(x, y, align='center', alpha=0.5, color=colors, width=0.3)
	ax.errorbar(x, y, yerr=[yerr_low, yerr_high], xerr=0.15, fmt='o', ecolor='r', mfc='r', ms=5, capsize=10)
	ax.set_xticks(x)
	ax.set_xticklabels(metrics, fontsize=parm.TICK_LABEL_FONTSIZE)
	ax.set_yticks(numpy.arange(0, 1.01, 0.05))
	ax.set_ylabel('Quality', fontsize=parm.AXIS_LABEL_FONTSIZE)
	ax.set_xlabel('Metrics', fontsize=parm.AXIS_LABEL_FONTSIZE)
	ax.set_title('Clustering Quality Bootstrap (' + str(len(bquality)) + ')', fontsize=parm.TITLE_FONTSIZE)
	fig.add_subplot(ax)

	fig.savefig(result_path + 'Clustering_B_Quality_c' + str(parameter.infile['dmin']) + '-' + str(parameter.infile['dmax']) + '.png', bbox_inches='tight')
	fig.savefig(result_path + 'Clustering_B_Quality_c' + str(parameter.infile['dmin']) + '-' + str(parameter.infile['dmax']) + '.eps', bbox_inches='tight')
	fig.savefig(result_path + 'Clustering_B_Quality_c' + str(parameter.infile['dmin']) + '-' + str(parameter.infile['dmax']) + '.pdf', bbox_inches='tight')

	with open(result_path + 'Clustering_B_Quality_c' + str(parameter.infile['dmin']) + '-' + str(parameter.infile['dmax']) + '.txt', 'w') as f:
		bootstrap = numpy.arange(len(bquality))
		f.write('#Bootstrap\t' + '\t'.join([str(b) for b in bootstrap]) + '\tMean\n')
		f.write('Purity\t' + '\t'.join([str(b) for b in purity]) + '\t' + str(y[2]) + '\n')
		f.write('ARI\t' + '\t'.join([str(b) for b in ari]) + '\t' + str(y[0]) + '\n')
		f.write('NMI\t' + '\t'.join([str(b) for b in nmi]) + '\t' + str(y[1]) + '\n')

	plt.close()
	print '[INFO] --> [PLOT] Done.'

# Measuring Purity, ARI and NMI in the quality mode, b and c quality measurements are two rigorous methods of quality measurements
# In type b, the number of dimensions are fixed, but the consensus clustering is done multiple time to measure the confidence factor of the method
# In type c, the number of selected dimensions is varios from 1 to at most the number of data points to measure the stability of the method
def cluster_bc_quality(parameter, _object):
	print '[INFO] --> [QUAL] Cluster quality measurement ...'

	kn = _object['kn']
	true_labels_dict = _object['true_labels']
	n_cell = _object['n_cell']

	consensus_cluster = _object['con_result'][str(kn[0])]
	pred_labels = consensus_cluster['sorted_labels']
	actual_labels = consensus_cluster['consensus_labels']
	true_labels = [true_labels_dict[actual_labels[k]] for k in range(n_cell)]
	metrics = ('ARI', 'NMI', 'Purity')
	y = [0] * len(metrics)
	y[2] = cluster_purity(true_labels, pred_labels, actual_labels, kn[0])['score']
	y[0] = adjusted_rand_score(true_labels, pred_labels)
	y[1] = normalized_mutual_info_score(true_labels, pred_labels)
	print '[INFO] --> [QUAL] Done.'
	return y

# Implementation of Purity formula to measure this quality factor
def cluster_purity(true_labels, pred_labels, actual_labels, k):
	n1 = len(true_labels)
	n2 = len(pred_labels)
	if n1 != n2:
		return 0.0
	score = 0.0
	cluster_mapping = {}
	true_labels_unique = numpy.unique(true_labels)
	false_list = [False] * len(true_labels_unique)
	assigned = dict(zip(true_labels_unique, false_list))
	for i in range(k):
		ith_pred_cluster = [actual_labels[x] for x in range(n2) if pred_labels[x] == i]
		jth_best_match = []
		for j in true_labels_unique:
			jth_true_cluster = [actual_labels[x] for x in range(n1) if true_labels[x] == j]
			jth_best_match.append(len(numpy.intersect1d(ith_pred_cluster, jth_true_cluster)))
		score += numpy.max(jth_best_match)
		cluster_mapping[i] = true_labels_unique[sorted([x for x in enumerate(jth_best_match) if assigned[true_labels_unique[x[0]]] == False], key=lambda z:z[1], reverse=True)[0][0]]
		assigned[cluster_mapping[i]] = True
	return dict(score=score / float(n1), map=cluster_mapping)

# Measuring the Purity, ARI and NMI quality factor of the consensus clustering results based on the given true labels and plotting the bar chart
# Plotting the miscalssification scatter plot based on purity class matching
# Running in the validate mode
def cluster_quality(parameter, _object):
	print '[INFO] --> [QUAL] Cluster quality measurement ...'
	global fig_num

	path = '/'.join(parameter.mica_labels.split('/')[:-1]) + '/'
	name = parameter.mica_labels.split('/')[-1].split('.')[0]
	result_path = path
	
	pred_labels_dict = _object['mica_labels']
	true_labels_dict = _object['true_labels']

	n = int(numpy.ceil(numpy.sqrt(len(pred_labels_dict))))
	fig = plt.figure(fig_num, figsize=(16, 16), dpi=300)
	fig_num += 1
	fig1 = plt.figure(fig_num, figsize=(16, 16), dpi=300)
	fig_num += 1
	grid = gridspec.GridSpec(n, n, wspace=0.2, hspace=0.2)
	grid1 = gridspec.GridSpec(n, n, wspace=0.6, hspace=0.2)
	t = 0

	for item in pred_labels_dict:
		pred_labels = []
		true_labels = []
		actual_labels = []
		n_cell = _object['n_cell']

		subgrid = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=grid[t], wspace=0.1, hspace=0.1)
		subgrid1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=grid1[t], wspace=0.1, hspace=0.1)

		for k, v in item.iteritems():
			actual_labels.append(k)
			pred_labels.append(v)
			true_labels.append(true_labels_dict[k])
		i = numpy.max(pred_labels) + 1
		metrics = ('ARI', 'NMI', 'Purity')
		x = numpy.arange(len(metrics))
		y = [0] * len(metrics)
		purity = cluster_purity(true_labels, pred_labels, actual_labels, i)
		y[2] = purity['score']
		y[0] = adjusted_rand_score(true_labels, pred_labels)
		y[1] = normalized_mutual_info_score(true_labels, pred_labels)
		colors = plt.cm.Accent(numpy.linspace(0, 1, len(metrics)))
		ax = plt.Subplot(fig, subgrid[0])
		ax.bar(x, y, align='center', alpha=0.5, color=colors, width=0.3)
		print '[INFO] --> [QUAL] ' + str(metrics) + ' : ' + str(y)
		ax.set_xticks(x)
		ax.set_xticklabels(metrics)
		ax.set_yticks(numpy.arange(0, 1.01, 0.05))
		ax.set_ylabel('Quality', fontsize=parm.AXIS_LABEL_FONTSIZE)
		ax.set_xlabel('Metrics', fontsize=parm.AXIS_LABEL_FONTSIZE)
		ax.set_title('Clustering Quality Measurement (' + str(i) + ')', fontsize=parm.TITLE_FONTSIZE)
		fig.add_subplot(ax)
		
		trans_labels = _object['trans']['labels']
		tsne = None
		trans_tsne = None
		if i not in _object['tsne_trans']:
			trans = _object['trans']['Y'][:,0:numpy.max(_object['dim_use'])]
			class_size = [b[1] for b in sorted(enumerate([pred_labels.count(z) for z in range(i)]), key=lambda z:z[1], reverse=True)]
			tsne = manifold.TSNE(n_components=2, n_iter=5000, learning_rate=200, perplexity=numpy.min([numpy.max(class_size),parm.MAX_PERPLEXITY]), random_state=10)
			trans_tsne = tsne.fit_transform(trans)
			_object['tsne_trans'][i] = trans_tsne
		else:
			trans_tsne = _object['tsne_trans'][i]
			
		ax1 = plt.Subplot(fig1, subgrid1[0])
		colors = plt.cm.jet(numpy.linspace(0, 1, i+1))	
		marker_size = parm.MAXIMUM_NORMAL_MARKER_SIZE if n_cell < parameter.MARKER_SIZE_THRESHOLD else parm.MINIMUM_NORMAL_MARKER_SIZE
		for z in range(i):
			xc = [trans_tsne[p, 0] for p in range(n_cell) if pred_labels[actual_labels.index(trans_labels[p])] == z and true_labels[actual_labels.index(trans_labels[p])] == purity['map'][z]]
			yc = [trans_tsne[p, 1] for p in range(n_cell) if pred_labels[actual_labels.index(trans_labels[p])] == z and true_labels[actual_labels.index(trans_labels[p])] == purity['map'][z]]
			ax1.scatter(xc, yc, facecolor=colors[z+1], s=marker_size, marker='o', vmin=0, vmax=i, label='T ' + str(z+1) + ' (' + str(len(xc)) + ')')
			xc = [trans_tsne[p, 0] for p in range(n_cell) if pred_labels[actual_labels.index(trans_labels[p])] == z and true_labels[actual_labels.index(trans_labels[p])] != purity['map'][z]]
			yc = [trans_tsne[p, 1] for p in range(n_cell) if pred_labels[actual_labels.index(trans_labels[p])] == z and true_labels[actual_labels.index(trans_labels[p])] != purity['map'][z]]
			ax1.scatter(xc, yc, facecolor=colors[z+1], s=marker_size, marker='^', vmin=0, vmax=i, label='F ' + str(z+1) + ' (' + str(len(xc)) + ')')
		ax1.set_title('Misclassified MICA plot (' + str(i) + ')', fontsize=parm.TITLE_FONTSIZE)
		ax1.set_ylabel('MICA-2', fontsize=parm.AXIS_LABEL_FONTSIZE)	
		ax1.set_xlabel('MICA-1', fontsize=parm.AXIS_LABEL_FONTSIZE)
		ax1.set_xticks([])
		ax1.set_yticks([])
		ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
		fig1.add_subplot(ax1)

		t += 1
	fig.savefig(result_path + 'Quality_' + name + '.png', bbox_inches='tight')
	fig.savefig(result_path + 'Quality_' + name + '.eps', bbox_inches='tight')
	fig.savefig(result_path + 'Quality_' + name + '.pdf', bbox_inches='tight')

	fig1.savefig(result_path + 'Misclassified_' + name + '.png', bbox_inches='tight')
	fig1.savefig(result_path + 'Misclassified_' + name + '.eps', bbox_inches='tight')
	fig1.savefig(result_path + 'Misclassified_' + name + '.pdf', bbox_inches='tight')	

	plt.close()

	print '[INFO] --> [QUAL] Done.'

# Plotting the flat cluster heatmap based on the consensus clustering results
# Running in the clustering mode
def heatmap_clusters(parameter, _object):
	print '[INFO] --> [PLOT] Plotting heatmap ...'
	
	global fig_num

	path = '/'.join(parameter.infile['path'].split('/')[:-1]) + '/'
	result_path = path + 'Results_' + parameter.job_name + '/'
	
	kn = _object['kn']
	n = int(numpy.ceil(numpy.sqrt(len(kn))))
	n1 = numpy.sqrt(n)
	t = 0

	fig = plt.figure(fig_num, figsize=(16, 16), dpi=300)
	fig_num += 1
	fig1 = plt.figure(fig_num, figsize=(16, 16), dpi=300)
	fig_num += 1
	grid = gridspec.GridSpec(n, n, wspace=0.2, hspace=0.2)
	grid1 = gridspec.GridSpec(n, n, wspace=0.2, hspace=0.2)
	
	for i in kn:
		#f = open(result_path + 'MICA_Clusters_k' + str(i) + '.txt', 'w')
		flat_clusters = _object['con_result'][str(i)]['flat_clusters']
		sorted_labels = _object['con_result'][str(i)]['sorted_labels']
		distance_matrix = _object['matrix']
		consensus_labels = list(_object['con_result'][str(i)]['consensus_labels'])
		mds_labels = list(_object['trans']['labels'])
		hc1_labels = list(_object['hc_trans']['hc_labels'])
		original_labels = list(_object['header'][1:])
		n_cells = _object['n_cell']
		hc_trans = numpy.zeros((n_cells, n_cells))
		hc_mi_trans = numpy.zeros((n_cells, n_cells))
		hc_labels = []
		
		#linkage_ = linkage(flat_clusters, method='ward', metric='euclidean')				# Create dendrogram linkage on ordered consensus clustering matrix

		subgrid = gridspec.GridSpecFromSubplotSpec(3, 2, subplot_spec=grid[t], wspace=0.001, hspace=0.1, height_ratios=[1,1,20], width_ratios=[21,1])
		subgrid1 = gridspec.GridSpecFromSubplotSpec(3, 2, subplot_spec=grid1[t], wspace=0.001, hspace=0.1, height_ratios=[1,1,20], width_ratios=[21,1])
		
		#ax1 = plt.Subplot(fig, subgrid[0])
		#ax1.set_title('MICA Result with ' + str(i) + ' clusters', fontsize=parm.TITLE_FONTSIZE/n1)
		#sys.setrecursionlimit(5000)
		#leaves = dendrogram(linkage_, no_plot=True, orientation='top', show_leaf_counts=True, labels=consensus_labels, leaf_font_size=2, color_threshold=0, above_threshold_color='b')['leaves']
		#ax1.set_xticks([])
		#ax1.set_yticks([])
		#ax1.axis('off')
		#fig.add_subplot(ax1)
		
		colorbar = []
		for k in range(n_cells):
			for j in range(n_cells):
				hc_mi_trans[k][j] = distance_matrix[original_labels.index(hc1_labels[hc1_labels.index(mds_labels[mds_labels.index(consensus_labels[k])])])][original_labels.index(hc1_labels[hc1_labels.index(mds_labels[mds_labels.index(consensus_labels[j])])])]
		hc_labels = consensus_labels
		colorbar = numpy.array(sorted_labels) + 1
		bar = numpy.array(colorbar)
		for k in range(0, 5):
			bar = numpy.vstack((bar, colorbar))

		ax2 = plt.Subplot(fig, subgrid[2])
		ax2.set_title('MICA Result (' + str(i) + ')', fontsize=parm.TITLE_FONTSIZE/n1)
		ax2.matshow(bar, cmap='jet', interpolation=None, aspect='auto', vmin=1, vmax=i+1)
		ax2.set_xticks([])
		ax2.set_yticks([])
		fig.add_subplot(ax2)
		
		ax3 = plt.Subplot(fig, subgrid[4])
		cax = ax3.matshow(flat_clusters, cmap='YlOrRd', interpolation=None, aspect='auto')
		ax3.set_xticks([])
		ax3.set_yticks([])
		#ax3.set_yticks(numpy.arange(len(hc_labels)))
		#ax3.set_yticklabels(hc_labels, fontsize=2)
		fig.add_subplot(ax3)
		
		ax4 = plt.Subplot(fig, subgrid[5])
		ax4.axis('off')
		cbar = fig.colorbar(cax, ax=ax4, shrink=1.5, aspect=20, fraction=.5, pad=0)
		cbar.ax.tick_params(labelsize=parm.TICK_LABEL_FONTSIZE/n1) 

		ax11 = plt.Subplot(fig1, subgrid1[2])
		ax11.set_title('Clustered Mutual Information (' + str(i) + ')', fontsize=parm.TITLE_FONTSIZE/n1)
		ax11.matshow(bar, cmap='jet', interpolation=None, aspect='auto', vmin=1, vmax=i+1)
		ax11.set_xticks([])
		ax11.set_yticks([])
		fig1.add_subplot(ax11)

		ax21 = plt.Subplot(fig1, subgrid1[4])
		cax1 = ax21.matshow(hc_mi_trans, cmap='YlOrRd', interpolation=None, aspect='auto')
		ax21.axis('off')
		fig1.add_subplot(ax21)

		ax31 = plt.Subplot(fig1, subgrid1[5])
		ax31.axis('off')
		ccbar = fig1.colorbar(cax1, ax=ax31, shrink=1.5, aspect=20, fraction=.5, pad=0)
		ccbar.ax.tick_params(labelsize=parm.TICK_LABEL_FONTSIZE/n1) 


		'''
		f.write('Cell\t' + '\t'.join(hc_labels) + '\n')
		for k in range(n_cells):
			f.write(hc_labels[k] + '\t' + '\t'.join(str(p) for p in flat_clusters[k]) + '\n')
		f.close()
		'''

		t += 1
	
	fig.savefig(result_path + 'MICA_Clusters.png', bbox_inches='tight')
	fig.savefig(result_path + 'MICA_Clusters.eps', bbox_inches='tight')
	fig.savefig(result_path + 'MICA_Clusters.pdf', bbox_inches='tight')

	fig1.savefig(result_path + 'MICA-MI.png', bbox_inches='tight')
	fig1.savefig(result_path + 'MICA-MI.eps', bbox_inches='tight')
	fig1.savefig(result_path + 'MICA-MI.pdf', bbox_inches='tight')

	plt.close()
	
	print '[INFO] --> [PLOT] Done.'

# Construct membership matrix for each KMeans result, compute the average membership matrix and apply hierarchical clustering for the final clustering (consensus)
# Running in the clustering, opt-k-2 and any quality mode, in the quality and opt-k-2, the plotting will remain off
def consensus_clustering(parameter, _object):
	print '[INFO] --> [CCLS] Consensus clustering ...'
	
	global fig_num
	
	fig = None
	if parameter.intermediate_plotting:
		fig = plt.figure(fig_num, figsize=(16, 16), dpi=300)
		fig_num += 1
	
	kn = _object['kn']
	k_len = len(kn)
	n = int(numpy.ceil(numpy.sqrt(k_len)))
	n1 = numpy.sqrt(n)
	grid = gridspec.GridSpec(n, n, wspace=0.5, hspace=0.2)

	path = '/'.join(parameter.infile['path'].split('/')[:-1]) + '/'
	result_path = path + 'Results_' + parameter.job_name + '/'
	
	min_dim = numpy.min(_object['dim_use'])
	max_dim = numpy.max(_object['dim_use'])
	max_bs = numpy.max(_object['bs'])
	f = open(result_path + 'Consensus_Clustering_c' + str(min_dim) + '-' + str(max_dim) + '_b' + str(max_bs + 1) + '.txt', 'w') if parameter.intermediate_plotting else None
	d = _object['d']
	n_cell = _object['n_cell']
	consensus_clustering_results = {}
	t = 0
	
	for i in kn:
		f1 = open(result_path + 'Consensus_Clustering_c' + str(min_dim) + '-' + str(max_dim) + '_b' + str(max_bs + 1) + '_k' + str(i) + '.txt', 'w') if parameter.intermediate_plotting else None
		km_result = _object['km_result'][str(i)]
		# Building consensus clustering matrix by averaging across all k-means results with same k and different number of eigen vectors and seed numbers
		clustering_matrix = numpy.zeros((n_cell, n_cell))
		for km in km_result.values():
			km = numpy.transpose(km)
			for j in range(0, d):
				for k in range(0, n_cell):
					clustering_matrix[k] = clustering_matrix[k] + [(1 if km[j][k] == km[j][l] else 0) for l in range(0, n_cell)]
		clustering_matrix /= (float(d) * len(km_result))																			# averaging the k-means matrices

		clustering = AgglomerativeClustering(linkage='ward', n_clusters=i, affinity='euclidean')									# hierarchical clustring
		clustering.fit(clustering_matrix)																							# fit the model into the consensus clustering matrix
		
		labels = [int(z) for z in clustering.labels_]																				# labels from hierarchical clustreing on consensus clusterint
		km_labels_count = sorted(enumerate([labels.count(z) for z in range(i)]), key=lambda p:p[1], reverse=True)
		class_orders = [b[0] for b in km_labels_count]																				# order the clusters based on the size of the cluster
		relabel = [class_orders.index(z) for z in labels]																			# labels after re-ordering based on size of the clusters
		sorted_labels = [b for b in sorted(enumerate(relabel), key=lambda p:p[1])]													# sorting labels based on the value of the cluster name (descending size)	
		consensus_labels = [_object['trans']['labels'][sorted_labels[z][0]] for z in range(n_cell)]							# preserve actual labels of cells after re-ordering and sorting

		# Create a flat clustering matrix with reordering the clustering matrix based on the re-ordered labels
		flat_clusters = numpy.zeros((n_cell, n_cell))
		for j in range(n_cell):
			for k in range(n_cell):
				flat_clusters[j][k] = clustering_matrix[sorted_labels[j][0]][sorted_labels[k][0]]
				#flat_clusters[j][k] = int(sorted_labels[j][1]) + 1 if sorted_labels[j][1] == sorted_labels[k][1] else 0
		sorted_labels = [z[1] for z in sorted_labels]
		
		if parameter.intermediate_plotting:
			# Scatter plot of first two components of MDS-TSNE with class labeled based on size
			tsne = None
			trans_tsne = None
			trans = _object['trans']['Y']
			if i not in _object['tsne_trans']:
				class_size = [b[1] for b in km_labels_count]
				tsne = manifold.TSNE(n_components=2, n_iter=5000, learning_rate=200, perplexity=numpy.min([numpy.max(class_size),parm.MAX_PERPLEXITY]), random_state=10)
				trans_tsne = tsne.fit_transform(trans[:,0:max_dim])
				_object['tsne_trans'] = trans_tsne
			else:
				trans_tsne = _object['tsne_trans'][i]
			
			subgrid = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=grid[t], wspace=0.1, hspace=0.1)
			ax = plt.Subplot(fig, subgrid[0])
			colors = plt.cm.jet(numpy.linspace(0, 1, i+1))
			marker_size = parm.MAXIMUM_NORMAL_MARKER_SIZE if len(relabel) < parameter.MARKER_SIZE_THRESHOLD else parm.MINIMUM_NORMAL_MARKER_SIZE
			marker_size /= n1
			xx = []
			yy = []
			xmds = [[]] * max_dim
			ids = []
			classes = []
			centers = []
			for z in range(i):
				xc = [trans_tsne[p, 0] for p in range(n_cell) if relabel[p] == z]
				yc = [trans_tsne[p, 1] for p in range(n_cell) if relabel[p] == z]
				ax.scatter(xc, yc, facecolor=colors[z+1], s=marker_size, marker='o', vmin=0, vmax=i, label=str(z+1) + ' (' + str(len(xc)) + ')')
				xx = xx + xc
				yy = yy + yc
				classes = classes + [z+1 for p in xc]
				ids = ids + [_object['trans']['labels'][p] for p in range(n_cell) if relabel[p] == z]
				centers.append([numpy.sum(xc) / float(len(xc)), numpy.sum(yc) / float(len(yc))])
				for kk in range(max_dim):
					xmds[kk] = xmds[kk] + [trans[p, kk] for p in range(n_cell) if relabel[p] == z]

			centers = numpy.array(centers)
			ax.scatter(centers[:,0], centers[:,1], marker='o', c='white', alpha=0.5, s=marker_size*3, edgecolor='k')
			for ii, cc in enumerate(centers):
				ax.scatter(cc[0], cc[1], marker='$%d$' % (ii+1), alpha=0.7, s=marker_size*2, edgecolor='k')
			f1.write('ids\t' + '\t'.join(ids) + '\n')
			f1.write('cluster\t' + '\t'.join([str(p) for p in classes]) + '\n')
			f1.write('TSNE-x-component\t' + '\t'.join([str(p) for p in xx]) + '\n')
			f1.write('TSNE-y-component\t' + '\t'.join([str(p) for p in yy]) + '\n')
			for kk in range(max_dim):
				f1.write('MDS-' + str(kk) + '-component\t' + '\t'.join([str(p) for p in xmds[kk]]) + '\n')
			ax.set_title('Consensus Clustering (' + str(i) + ')', fontsize=parm.TITLE_FONTSIZE/n1)
			ax.set_ylabel('MICA-2', fontsize=parm.AXIS_LABEL_FONTSIZE/n1)	
			ax.set_xlabel('MICA-1', fontsize=parm.AXIS_LABEL_FONTSIZE/n1)
			ax.set_xticks([])
			ax.set_yticks([])
			ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=parm.TICK_LABEL_FONTSIZE/n1)
			fig.add_subplot(ax)
		
			f.write('\t'.join(str(x) for x in consensus_labels) + '\n')
			f.write('\t'.join([str(x) for x in sorted_labels]) + '\n')
			f1.close()
		consensus_clustering_results[str(i)] = dict(consensus=clustering_matrix, consensus_labels=consensus_labels, sorted_labels=sorted_labels, flat_clusters=flat_clusters, leaves=clustering.children_)
		t += 1

	if parameter.intermediate_plotting:
		plt.savefig(result_path + 'Consensus_Clustering_c' + str(numpy.min(_object['dim_use'])) + '-' + str(numpy.max(_object['dim_use'])) + '_b' + str(numpy.max(_object['bs']) + 1) + '.png', bbox_inches='tight')
		plt.savefig(result_path + 'Consensus_Clustering_c' + str(numpy.min(_object['dim_use'])) + '-' + str(numpy.max(_object['dim_use'])) + '_b' + str(numpy.max(_object['bs']) + 1) + '.eps', bbox_inches='tight')
		plt.savefig(result_path + 'Consensus_Clustering_c' + str(numpy.min(_object['dim_use'])) + '-' + str(numpy.max(_object['dim_use'])) + '_b' + str(numpy.max(_object['bs']) + 1) + '.pdf', bbox_inches='tight')
		plt.close()
		f.close()
	
	_object['con_result'] = consensus_clustering_results
	
	print '[INFO] --> [CCLS] Done.'

# Single kmeans running on thread
def single_kmeans(parameter, _object, fig, grid, n_cell, i, k, km_result, relabel):
	trans = _object['trans']['Y'][:,0:_object['dim_use'][i]]				# Retrieving MDS-transformed distance matrix
	n = int(numpy.ceil(numpy.sqrt(_object['d'])))
	n /= numpy.sqrt(n)
	km = cluster.KMeans(n_clusters=k, max_iter=1000, n_init=1000)		# Running KMeans
	km_res = km.fit_predict(trans)
	mutex.acquire()
	try:
		km_result[i] = km_res
		km_labels = [int(z) for z in km_result[i]]
		# Reordering the clusters based on size and relabel independent variables based on the new class label
		km_labels_count = sorted(enumerate([km_labels.count(z) for z in range(k)]), key=lambda z:z[1], reverse=True)
		class_orders = [b[0] for b in km_labels_count]
		relabel[i] = [class_orders.index(z) for z in km_result[i]]
	
		if parameter.intermediate_plotting_:
			tsne = None
			trans_tsne = None
			if k not in _object['tsne_trans']:
				class_size = [b[1] for b in km_labels_count]
				trans = _object['trans']['Y'][:,0:numpy.max(_object['dim_use'])]
				tsne = manifold.TSNE(n_components=2, n_iter=5000, learning_rate=200, perplexity=numpy.min([numpy.max(class_size),parm.MAX_PERPLEXITY]), random_state=10)
				trans_tsne = tsne.fit_transform(trans)
				_object['tsne_trans'][k] = trans_tsne
			else:
				trans_tsne = _object['tsne_trans'][k]

			subgrid = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=grid[i], wspace=0.1, hspace=0.1)
			ax = plt.Subplot(fig, subgrid[0])
			colors = plt.cm.jet(numpy.linspace(0, 1, k+1))
			marker_size = parm.MAXIMUM_NORMAL_MARKER_SIZE if len(relabel[i]) < parameter.MARKER_SIZE_THRESHOLD else parm.MINIMUM_NORMAL_MARKER_SIZE
			marker_size /= n
			centers = []
			for z in range(k):
				xc = [trans_tsne[p, 0] for p in range(n_cell) if relabel[i][p] == z]
				yc = [trans_tsne[p, 1] for p in range(n_cell) if relabel[i][p] == z]
				ax.scatter(xc, yc, facecolor=colors[z+1], s=marker_size, marker='o', label=str(z+1) + ' (' + str(len(xc)) + ')', vmin=0, vmax=k)
				centers.append([numpy.sum(xc) / float(len(xc)), numpy.sum(yc) / float(len(yc))])
			centers = numpy.array(centers)
			ax.scatter(centers[:,0], centers[:,1], marker='o', c='white', alpha=0.5, s=marker_size*3, edgecolor='k')
			for ii, cc in enumerate(centers):
				ax.scatter(cc[0], cc[1], marker='$%d$' % (ii+1), alpha=0.7, s=marker_size*2, edgecolor='k')
			ax.set_title('K-Means (' + str(k) + ',' + str(_object['dim_use'][i]) + ')', fontsize=parm.TITLE_FONTSIZE/n)
			ax.set_ylabel('MICA-2', fontsize=parm.AXIS_LABEL_FONTSIZE/n)	
			ax.set_xlabel('MICA-1', fontsize=parm.AXIS_LABEL_FONTSIZE/n)
			ax.set_xticks([])
			ax.set_yticks([])
			ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
			fig.add_subplot(ax)
	finally:
		mutex.release()

# Single bootstrap, single k, combined components running on threads
def single_kn_clustering(parameter, _object, result_path, k, bsi, n_cell, km_results, fig):
	n = int(numpy.ceil(numpy.sqrt(_object['d'])))
	km_result = [[]] * _object['d']
	relabel = [[]] * _object['d']
	d_threads = []

	grid = gridspec.GridSpec(n, n, wspace=0.4, hspace=0.2)
	for i in range(0, _object['d']):
		print '[INFO] --> [KMNS] Running #' + str(len(d_threads) + 1) + ' kmeans ...'
		d_threads.append(Thread(target=single_kmeans, args=(parameter, _object, fig, grid, n_cell, i, k, km_result, relabel)))
		d_threads[-1].start()
		time.sleep(1)
	for d_thread in d_threads:
		d_thread.join()

	if parameter.intermediate_plotting_:
		fig.savefig(result_path + 'KMeans_k' + str(k) + '_b' + str(bsi+1) + '.png', bbox_inches='tight')
		fig.savefig(result_path + 'KMeans_k' + str(k) + '_b' + str(bsi+1) + '.eps', bbox_inches='tight')
		fig.savefig(result_path + 'KMeans_k' + str(k) + '_b' + str(bsi+1) + '.pdf', bbox_inches='tight')
	
		with open(result_path + 'KMeans_k' + str(k) + '_b' + str(bsi+1) + '.txt', 'w') as f:
			f.write('\t'.join([str(x) for x in _object['trans']['labels']]) + '\n')
			for v in relabel:
				f.write('\t'.join([str(x+1) for x in v]) + '\n')
	
	mutex1.acquire()
	try:
		km_result = numpy.transpose(numpy.array(km_result))
		if str(k) not in km_results:
			km_results[str(k)] = {}
		km_results[str(k)][str(bsi)] = km_result
	finally:
		mutex1.release()

# Single bootstrap running on thread
def single_bootstrap_clustering(parameter, _object, result_path, kn, bsi, n_cell, km_results):
	k_threads = []																				# Running KMeans over number of bootstraps
	for k in kn:																			# Running KMeans for provided Ks 
		print '[INFO] --> [KMNS] Running kmeans for k = ' + str(k) + ' ...'
		mutex2.acquire()
		try:
			global fig_num
			fig = None
			if parameter.intermediate_plotting_:
				fig = plt.figure(fig_num, figsize=(16, 16), dpi=300)
				fig_num += 1
			k_threads.append(Thread(target=single_kn_clustering, args=(parameter, _object, result_path, k, bsi, n_cell, km_results, fig)))
			k_threads[-1].start()
			time.sleep(_object['d'])
		finally:
			mutex2.release()
	for k_thread in k_threads:
		k_thread.join()

# Clustering the MDS-transformed distance matrix over selected dimensions and with random seed for multiple times 
# Running in the clustering, opt-k-2 and any quality mode, in the quality and opt-k-2, the plotting will remain off
def kmeans(parameter, _object):
	print '[INFO] --> [KMNS] K-Means clustering ...'

	path = '/'.join(parameter.infile['path'].split('/')[:-1]) + '/'
	result_path = path + 'Results_' + parameter.job_name + '/'

	km_results = {}
	kn = _object['kn']
	bs = _object['bs']
	n_cell = _object['n_cell']
	b_threads = []

	for bsi in bs:
		print '[INFO] --> [KMNS] Running kmeans for b = ' + str(bsi) + ' ...'
		b_threads.append(Thread(target=single_bootstrap_clustering, args=(parameter, _object, result_path, kn, bsi, n_cell, km_results)))
		b_threads[-1].start()
		time.sleep(_object['d'] * len(kn))
	for b_thread in b_threads:
		b_thread.join()

	_object['km_result'] = km_results

	print '[INFO] --> [KMNS] Done.'

# MDS transformation of the distance matrix
# Running if MDS decomposition method is selected in all modes
def cmdscale(parameter, _object):
	print '[INFO] --> [CMDS] Classic Multidimensional Scaling ...'
	global fig_num

	path = '/'.join(parameter.infile['path'].split('/')[:-1]) + '/'
	result_path = path + 'Results_' + parameter.job_name + '/'

	distance = _object['hc_trans']['hc_matrix']
	n = distance.shape[0]															# Number of points
	
	H = numpy.eye(n) - numpy.ones((n, n))/n 										# Centering matrix
	B = -H.dot(distance ** 2).dot(H)/2									# YY^T
	evals, evecs = numpy.linalg.eigh(B)												# Diagonalize
	idx   = numpy.argsort(evals)[::-1]												# Sort by eigenvalue in descending order
	evals = evals[idx]
	evecs = evecs[:,idx]	
	# Compute the coordinates using positive-eigenvalued components only
	w, = numpy.where(evals > 0)
	L  = numpy.diag(numpy.sqrt(evals[w]))
	V  = evecs[:,w]
	Y  = V.dot(L)

	#mds = manifold.MDS(n_components=n, dissimilarity='precomputed', random_state=10)#, max_iter=3000, eps=1e-9)
	#Y = mds.fit_transform(numpy.array(distance))#.embedding_
	mds_labels = _object['hc_trans']['hc_labels'] 									# Reordering labels as mds_labels

	if parameter.intermediate_plotting:
		fig = plt.figure(fig_num, figsize=(16, 16), dpi=300)
		fig_num += 1
		
		tsne = manifold.TSNE(n_components=2, n_iter=5000, learning_rate=200, perplexity=50, random_state=10)
		Y_tsne = tsne.fit_transform(Y[:,0:numpy.max(_object['dim_use'])])
		marker_size = parm.MAXIMUM_NORMAL_MARKER_SIZE if n < parameter.MARKER_SIZE_THRESHOLD else parm.MINIMUM_NORMAL_MARKER_SIZE
		plt.scatter(Y_tsne[:,0], Y_tsne[:,1], s=marker_size, marker='o', facecolors='none', edgecolors='r')
		plt.title('MI-MDS', fontsize=parm.TITLE_FONTSIZE)
		plt.ylabel('MICA-2', fontsize=parm.AXIS_LABEL_FONTSIZE)	
		plt.xlabel('MICA-1', fontsize=parm.AXIS_LABEL_FONTSIZE)
		plt.xticks([])
		plt.yticks([])
		
		plt.savefig(result_path + 'MI-CMDS.png', bbox_inches='tight')
		plt.savefig(result_path + 'MI-CMDS.eps', bbox_inches='tight')
		plt.savefig(result_path + 'MI-CMDS.pdf', bbox_inches='tight')
		plt.close()
	
		''' Writing out the transformation/decomposition matrix
		with open(result_path + 'MI-CMDS.txt', 'w') as f:
			for v in Y.tolist():
				f.write(mds_labels[Y.tolist().index(v)] + '\t' + '\t'.join([str(x) for x in v]) + '\n')
		'''
	
	_object['mds_trans'] = dict(Y=Y, labels=mds_labels)							# Updating MICA object 
	
	print '[INFO] --> [CMDS] Done.'

# Laplacian decompositin method which uses EigenMap implementation
# Running if LPL decomposition method is selected in all modes or LPCA2 in clustering, opt-k-2 or any quality mode
def laplacian(parameter, _object):
	print '[INFO] --> [CLPL] Laplacian Decomposition ...'
	global fig_num

	path = '/'.join(parameter.infile['path'].split('/')[:-1]) + '/'
	result_path = path + 'Results_' + parameter.job_name + '/'

	distance = _object['hc_trans']['hc_matrix']
	n = distance.shape[0]															# Number of points
	
	lpl = manifold.SpectralEmbedding(n_components=n, eigen_solver='lobpcg', n_jobs=-1)
	Y = lpl.fit_transform(distance)
	lpl_labels = _object['hc_trans']['hc_labels'] 									# Reordering labels as lpl_labels

	if parameter.intermediate_plotting:
		fig = plt.figure(fig_num, figsize=(16, 16), dpi=300)
		fig_num += 1
		
		tsne = manifold.TSNE(n_components=2, n_iter=5000, learning_rate=200, perplexity=50, random_state=10)
		Y_tsne = tsne.fit_transform(Y[:,0:numpy.max(_object['dim_use'])])
		marker_size = parm.MAXIMUM_NORMAL_MARKER_SIZE if n < parameter.MARKER_SIZE_THRESHOLD else parm.MINIMUM_NORMAL_MARKER_SIZE
		plt.scatter(Y_tsne[:,0], Y_tsne[:,1], s=marker_size, marker='o', facecolors='none', edgecolors='r')
		plt.title('MI-LPL', fontsize=parm.TITLE_FONTSIZE)
		plt.ylabel('MICA-2', fontsize=parm.AXIS_LABEL_FONTSIZE)	
		plt.xlabel('MICA-1', fontsize=parm.AXIS_LABEL_FONTSIZE)
		plt.xticks([])
		plt.yticks([])
		
		plt.savefig(result_path + 'MI-LPL.png', bbox_inches='tight')
		plt.savefig(result_path + 'MI-LPL.eps', bbox_inches='tight')
		plt.savefig(result_path + 'MI-LPL.pdf', bbox_inches='tight')
		plt.close()
		
		''' Writing out the transformation/decomposition matrix
		with open(result_path + 'MI-LPL.txt', 'w') as f:
			for v in Y.tolist():
				f.write(lpl_labels[Y.tolist().index(v)] + '\t' + '\t'.join([str(x) for x in v]) + '\n')
		'''
	
	_object['lpl_trans'] = dict(Y=Y, labels=lpl_labels)							# Updating MICA object 
	
	print '[INFO] --> [CLPL] Done.'

# PCA transformation of the distance matrix --> Using the rotation matrix not the transformed (scaled) data points
# Running if PCA decomposition method is selected in all modes or LPCA2 in clustering, opt-k-2 or any quality mode
def pca(parameter, _object):
	print '[INFO] --> [CPCA] Principal Component Analysis ...'
	global fig_num

	path = '/'.join(parameter.infile['path'].split('/')[:-1]) + '/'
	result_path = path + 'Results_' + parameter.job_name + '/'

	distance = _object['hc_trans']['hc_matrix']
	n = distance.shape[0]															# Number of points
	
	pca = decomposition.PCA(n_components=n)
	Y = numpy.transpose(pca.fit(distance).components_)								# PCA model or the rotation matrix based on the distance matrix --> This is a wrong way of using PCA from the author's point of view
	pca_labels = _object['hc_trans']['hc_labels'] 									# Reordering labels as pca_labels

	if parameter.intermediate_plotting:
		fig = plt.figure(fig_num, figsize=(16, 16), dpi=300)
		fig_num += 1
		
		tsne = manifold.TSNE(n_components=2, n_iter=5000, learning_rate=200, perplexity=50, random_state=10)
		Y_tsne = tsne.fit_transform(Y[:,0:numpy.max(_object['dim_use'])])
		marker_size = parm.MAXIMUM_NORMAL_MARKER_SIZE if n < parameter.MARKER_SIZE_THRESHOLD else parm.MINIMUM_NORMAL_MARKER_SIZE
		plt.scatter(Y_tsne[:,0], Y_tsne[:,1], s=marker_size, marker='o', facecolors='none', edgecolors='r')
		plt.title('MI-PCA', fontsize=parm.TITLE_FONTSIZE)
		plt.ylabel('MICA-2', fontsize=parm.AXIS_LABEL_FONTSIZE)	
		plt.xlabel('MICA-1', fontsize=parm.AXIS_LABEL_FONTSIZE)
		plt.xticks([])
		plt.yticks([])
		
		plt.savefig(result_path + 'MI-PCA.png', bbox_inches='tight')
		plt.savefig(result_path + 'MI-PCA.eps', bbox_inches='tight')
		plt.savefig(result_path + 'MI-PCA.pdf', bbox_inches='tight')
		plt.close()
	
		''' Writing out the transformation/decomposition matrix
		with open(result_path + 'MI-PCA.txt', 'w') as f:
			for v in Y.tolist():
				f.write(pca_labels[Y.tolist().index(v)] + '\t' + '\t'.join([str(x) for x in v]) + '\n')
		'''
	
	_object['pca_trans'] = dict(Y=Y, labels=pca_labels)							# Updating MICA object 
	
	print '[INFO] --> [CPCA] Done.'

# Sequential transformation of the distance matrix using PCA and Laplacian (eigenmap) in order
# Running if LPCA decomposition method is selected in all modes
def laplacian_pca(parameter, _object):
	print '[INFO] --> [LPCA] Laplacian+PCA Decomposition ...'
	global fig_num

	path = '/'.join(parameter.infile['path'].split('/')[:-1]) + '/'
	result_path = path + 'Results_' + parameter.job_name + '/'

	distance = _object['hc_trans']['hc_matrix']
	n = distance.shape[0]															# Number of points
	
	pca = decomposition.PCA(n_components=n)
	Y = numpy.transpose(pca.fit(distance).components_)				# PCA model or the rotation matrix based on the distance matrix --> This is a wrong way of using PCA from the author's point of view
	lpl = manifold.SpectralEmbedding(n_components=n, eigen_solver='lobpcg', n_jobs=-1)
	Y = lpl.fit_transform(Y)
	lpca_labels = _object['hc_trans']['hc_labels'] 									# Reordering labels as lpca_labels

	if parameter.intermediate_plotting:
		fig = plt.figure(fig_num, figsize=(16, 16), dpi=300)
		fig_num += 1
		
		tsne = manifold.TSNE(n_components=2, n_iter=5000, learning_rate=200, perplexity=50, random_state=10)
		Y_tsne = tsne.fit_transform(Y[:,0:numpy.max(_object['dim_use'])])
		marker_size = parm.MAXIMUM_NORMAL_MARKER_SIZE if n < parameter.MARKER_SIZE_THRESHOLD else parm.MINIMUM_NORMAL_MARKER_SIZE
		plt.scatter(Y_tsne[:,0], Y_tsne[:,1], s=marker_size, marker='o', facecolors='none', edgecolors='r')
		plt.title('MI-LPCA', fontsize=parm.TITLE_FONTSIZE)
		plt.ylabel('MICA-2', fontsize=parm.AXIS_LABEL_FONTSIZE)	
		plt.xlabel('MICA-1', fontsize=parm.AXIS_LABEL_FONTSIZE)
		plt.xticks([])
		plt.yticks([])
		
		plt.savefig(result_path + 'MI-LPCA.png', bbox_inches='tight')
		plt.savefig(result_path + 'MI-LPCA.eps', bbox_inches='tight')
		plt.savefig(result_path + 'MI-LPCA.pdf', bbox_inches='tight')
		plt.close()
	
		''' Writing out the transformation/decomposition matrix
		with open(result_path + 'MI-LPCA.txt', 'w') as f:
			for v in Y.tolist():
				f.write(lpca_labels[Y.tolist().index(v)] + '\t' + '\t'.join([str(x) for x in v]) + '\n')
		'''
	
	_object['lpca_trans'] = dict(Y=Y, labels=lpca_labels)							# Updating MICA object 
	
	print '[INFO] --> [LPCA] Done.'

# Plotting original distance matrix and hierarchical-transformed distance matrix
def heatmap(parameter, _object):
	print '[INFO] --> [HCLS] Hierarchical Clustering ...'
	
	global fig_num

	path = '/'.join(parameter.infile['path'].split('/')[:-1]) + '/'
	result_path = path + 'Results_' + parameter.job_name + '/'

	# Initialization
	n = _object['n_cell']
	distance_matrix = _object['matrix']
	original_labels = _object['header'][1:]
	hc_trans = numpy.zeros((n, n))
	leaves_labels = []

	linkage_ = linkage(distance_matrix, method='ward')											# Applying hierarchical clustering on distance matrix

	if parameter.intermediate_plotting:															# In clustering mode all figures will be plotted
		# Original distance matrix plot
		fig = plt.figure(fig_num, figsize=(16, 16), dpi=300)
		fig_num += 1
		grid = gridspec.GridSpec(1, 1, wspace=0.2, hspace=0.2)
		subgrid = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=grid[0], wspace=0.001, hspace=0.01, height_ratios=[1,10], width_ratios=[10,1])

		ax1 = plt.Subplot(fig, subgrid[2])
		ax1.set_title('Mutual Information-Based Dissimilarity Matrix', fontsize=parm.TITLE_FONTSIZE)
		dax = ax1.matshow(distance_matrix, cmap='YlOrRd', interpolation=None, aspect='auto')
		ax1.set_xticks([])
		ax1.set_yticks([])
		#ax1.set_yticks(numpy.arange(original_labels.shape[0]))
		#ax1.set_yticklabels(original_labels, fontsize=2)
		fig.add_subplot(ax1)

		ax2 = plt.Subplot(fig, subgrid[3])
		ax2.axis('off')
		cbar = fig.colorbar(dax, ax=ax2, shrink=1.5, aspect=20, fraction=.5, pad=0)
		cbar.ax.tick_params(labelsize=parm.TICK_LABEL_FONTSIZE) 

		plt.savefig(result_path + 'dissimilarity_matrix.png', bbox_inches='tight')
		plt.savefig(result_path + 'dissimilarity_matrix.eps', bbox_inches='tight')
		plt.savefig(result_path + 'dissimilarity_matrix.pdf', bbox_inches='tight')
		plt.close()
		
		# Hierarchical-transformed distance matrix plot with dendrogram
		fig = plt.figure(fig_num, figsize=(16,16), dpi=100)
		fig_num += 1
		grid = gridspec.GridSpec(1, 1, wspace=0.2, hspace=0.2)
		subgrid = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=grid[0], wspace=0.01, hspace=0.01, height_ratios=[1,10], width_ratios=[10,1])
		
		ax1 = plt.Subplot(fig, subgrid[0])
		ax1.set_title('Mutual Information-Based Dissimilarity Matrix', fontsize=parm.TITLE_FONTSIZE)
		leaves = dendrogram(linkage_, ax=ax1, orientation='top', labels=original_labels, leaf_font_size=2, color_threshold=0)['leaves']
		ax1.axis('off')
		fig.add_subplot(ax1)

		ax3 = plt.Subplot(fig, subgrid[2])
		leaves = numpy.asarray(leaves)
		leaves = leaves[numpy.where(leaves < n)[0]]
		for i in range(n):																		# Put all independent variables together to construct a flat matrix based on hierarchical clustering result
			for j in range(n):
				hc_trans[i][j] = distance_matrix[leaves[i]][leaves[j]]
			leaves_labels.append(original_labels[leaves[i]])
		cax = ax3.matshow(hc_trans, cmap='YlOrRd', interpolation=None, aspect='auto')
		ax3.set_xticks([])
		ax3.set_yticks([])
		#ax3.set_yticks(numpy.arange(len(leaves_labels)))
		#ax3.set_yticklabels(leaves_labels, fontsize=2)
		fig.add_subplot(ax3)
		
		ax4 = plt.Subplot(fig, subgrid[3])
		ax4.axis('off')
		ccbar = fig.colorbar(cax, ax=ax4, shrink=1.5, aspect=20, fraction=.5, pad=0)
		ccbar.ax.tick_params(labelsize=parm.TICK_LABEL_FONTSIZE)
		
		plt.savefig(result_path + 'dissimilarity_matrix_dendrogram.png', bbox_inches='tight')
		plt.savefig(result_path + 'dissimilarity_matrix_dendrogram.eps', bbox_inches='tight')
		plt.savefig(result_path + 'dissimilarity_matrix_dendrogram.pdf', bbox_inches='tight')
		plt.close()
		'''
		with open(result_path + 'dissimilarity_matrix.txt', 'w') as f:
			f.write('\t'.join(leaves_labels) + '\n')
			for row in distance_matrix:
				f.write('\t'.join([str(x) for x in row]) + '\n')
		'''
	else:																						# In modes other than clustering the intermediate plotting will be discarded
		leaves = dendrogram(linkage_, no_plot=True, orientation='top', labels=original_labels, leaf_font_size=2, color_threshold=0)['leaves']
		leaves = numpy.asarray(leaves)
		leaves = leaves[numpy.where(leaves < n)[0]]
		for i in range(n):																		# Put all independent variables together to construct a flat matrix based on hierarchical clustering result
			for j in range(n):
				hc_trans[i][j] = distance_matrix[leaves[i]][leaves[j]]
			leaves_labels.append(original_labels[leaves[i]])

	_object['hc_trans'] = dict(hc_matrix=hc_trans, hc_labels=numpy.asarray(leaves_labels))						# Updating hierarchical-transformed matrix and the re-ordred labels on MICA object
	
	print '[INFO] --> [HCLS] Done.'

# Constructing the MICA object based on the given parameter which will be completed over the process
# First step to read the input parameters and files
def setup(parameter):
	print '[INFO] --> [PREP] Preparing the input file with the given parameters ...'
	# Initial parameters of MICA, given by user or set as default
	dmin = parameter.infile['dmin']
	dmax = parameter.infile['dmax']
	kn = parameter.infile['kn']
	bs = parameter.infile['bs']
	true_labels = None
	mica_labels = None
	expression = None
	biomarkers = None
	
	
	infile = parameter.infile['path']
	delimiter = parameter.infile['delimiter']
	reader = pandas.read_csv(infile, sep=delimiter, header=None, dtype=str)									# Reading the Mutual Information matrix calculated by MIE or any other packages, NxN matrix, each row and column represent an independent variable
	header = reader.iloc[0,:].values.astype(str)															# Header: Independent variables names
	n_cell = reader.shape[0] - 1																			# Number of independent variables
	n_dim = numpy.arange(int(dmin), int(dmax) + 1)															# Selected dimension set for MICA
	distance_matrix = 1 - reader.iloc[1:,1:].values.astype(float)											# The distance matrix built upon Mutual Information
	'''
	with open(parameter.infile['path'], 'r') as infile:														# Reading the Mutual Information matrix calculated by MIE or any other packages, NxN matrix, each row and column represent an independent variable
		reader = csv.reader(infile, delimiter=parameter.infile['delimiter'])								# CSV reader
		header = reader.next()																				# Header: Independent variables names
		n_cell = len(header) - 1																			# Number of independent variables
		n_dim = [x for x in range(int(dmin), int(dmax) + 1)] 												# Selected dimension set for MICA
		distance_matrix = []																				# The distance matrix built upon Mutual Information
		for row in reader:
			distance_matrix.append(1.0 - numpy.array([float(x) for x in row[1:]]))							# The distance or dissimilarity is 1 - normalized_mutual_information
	'''
	# Making a directory for the MICA results if not exists
	path = '/'.join(parameter.infile['path'].split('/')[:-1]) + '/'											
	result_path = path + 'Results_' + parameter.job_name
	if not os.path.exists(result_path):
		os.mkdir(result_path)
	if parameter.mode == 'validate' or parameter.mode == 'b-quality' or parameter.mode == 'c-quality':		# Reading ground truth if running mode is validate, b-quality or c-quality
		if parameter.true_labels != None:
			
			t_reader = pandas.read_csv(parameter.true_labels, sep="\t", header=None, dtype=str)
			actual = t_reader.iloc[0,:].values.astype(str)
			labels = t_reader.iloc[1,:].values.astype(int)
			true_labels = dict(zip(actual, labels))
			'''
			with open(parameter.true_labels, 'r') as f:														# Reading ground truth 
				reader = csv.reader(f, delimiter='\t')
				actual = reader.next()
				labels = [int(x) for x in reader.next()]
				true_labels = dict(zip(actual, labels))
			'''
	if parameter.mode == 'validate' or parameter.mode == 'overlay':											# Reading clustering result if running mode is validate or overlay
		mica_labels = []
		if parameter.mica_labels != None:
			
			m_reader = pandas.read_csv(parameter.mica_labels, sep="\t", header=None, dtype=str)
			for i in range(0, m_reader.shape[0], 2):
				actual = m_reader.iloc[i,:].values.astype(str)
				labels = m_reader.iloc[i+1,:].values.astype(int)
				mica_labels.append(dict(zip(actual, labels)))
			'''
			with open(parameter.mica_labels, 'r') as f:														# Reading the clustering result
				reader = list(csv.reader(f, delimiter='\t'))
				for i in range(0, len(reader), 2):
					actual = reader[i]
					labels = [int(x) for x in reader[i + 1]]
					mica_labels.append(dict(zip(actual, labels)))		
			'''
	if parameter.mode == 'overlay':																			# Reading expression matrix if running mode is overlay
		biomarkers = parameter.biomarkers
		
		exp = parameter.expression['path']
		exp_delim = parameter.expression['delimiter']
		exp_ncol = parameter.expression['ncol']
		e_reader = pandas.read_csv(exp, sep=exp_delim, header=None, dtype=str)
		exp_header = e_reader.iloc[0,exp_ncol:].values.astype(str)
		expression_matrix = e_reader.iloc[1:,exp_ncol:].values.astype(float)
		cell_ids = e_reader.iloc[1:,0].values.astype(str)
		expression = dict(header=exp_header, cells=cell_ids, matrix=expression_matrix)
		'''
		with open(parameter.expression['path'], 'r') as exp:
			reader = csv.reader(exp, delimiter=parameter.expression['delimiter'])
			exp_header = reader.next()[parameter.expression['ncol']:]
			expression_matrix = []
			cell_ids = []
			for row in reader:
				entry = []
				cell_ids.append(row[0])
				for i in range(parameter.expression['ncol'], len(row)):
					entry.append(float(row[i]))
				expression_matrix.append(entry)
			expression = dict(header=exp_header, cells=cell_ids, matrix=expression_matrix)
		'''
	_object = dict(header=header, matrix=distance_matrix, n_cell=n_cell, kn=kn, bs=bs, dim_use=n_dim, d=len(n_dim), hc_trans=None, mds_trans=None, km_result=None, con_result=None, mica_labels=mica_labels, true_labels=true_labels, expression=expression, biomarkers=biomarkers, tsne_trans={})

	print '[INFO] --> [PREP] Done.'
	return _object

# Managing the different modes of MICA and building the workflow based on each mode
def run(parameter):
	# Program will be running on local machine 
	if parameter.host['host'] == 'local':	
		print '[INFO] --> [MAIN] Running on local machine ...'
		_object = setup(parameter)							 							# Construct the object for the entire process
		kn_original = _object['kn']
		if _object['kn'] == None and (parameter.mode == 'clustering' or parameter.mode == 'opt-k-2'):
			plotting = parameter.intermediate_plotting
			bs = _object['bs']
			parameter.intermediate_plotting = False
			_object['bs'] = [1]
			_object['kn'] = parm.DEFAULT_KN
			heatmap(parameter, _object)													# Plotting heatmap of the Mutual Information-based distance matrix
			if parameter.decomposition == 'MDS':
				cmdscale(parameter, _object)											# Computing MDS transformation of distance matrix and plotting TSNE transformation of MDS result
				_object['trans'] = _object['mds_trans']
				kmeans(parameter, _object)												# Applying bootstraps of KMeans on various dimensions of MDS matrix 
			elif parameter.decomposition == 'LPL':
				laplacian(parameter, _object)
				_object['trans'] = _object['lpl_trans']
				kmeans(parameter, _object)												# Applying bootstraps of KMeans on various dimensions of MDS matrix 
			elif parameter.decomposition == 'PCA':
				pca(parameter, _object)
				_object['trans'] = _object['pca_trans']
				kmeans(parameter, _object)												# Applying bootstraps of KMeans on various dimensions of MDS matrix 
			elif parameter.decomposition == 'LPCA':
				laplacian_pca(parameter, _object)
				_object['trans'] = _object['lpca_trans']
				kmeans(parameter, _object)												# Applying bootstraps of KMeans on various dimensions of MDS matrix 
			elif parameter.decomposition == 'LPCA2':
				_object['bs'] = numpy.array(_object['bs']) * 2
				double_bs = _object['bs']
				laplacian(parameter, _object)
				_object['trans'] = _object['lpl_trans']
				kmeans(parameter, _object)
				_object['bs'] = numpy.array(_object['bs']) + 1
				double_bs = double_bs + _object['bs']
				pca(parameter, _object)
				_object['trans'] = _object['pca_trans']
				kmeans(parameter, _object)
				_object['bs'] = double_bs
				laplacian_pca(parameter, _object)
				_object['trans'] = _object['lpca_trans']
			consensus_clustering(parameter, _object)									# Averaging all KMeans results and applying hierarchical clustring			
			_object['kn'] = optimal_k_2(parameter, _object)
			#_object['kn'] = [_object['kn'][len(_object['kn']) / 2]]
			_object['bs'] = bs
			parameter.intermediate_plotting = plotting
		if parameter.mode == 'clustering' or (parameter.mode == 'opt-k-2' and kn_original != None):				# Clustering mode 
			heatmap(parameter, _object)													# Plotting heatmap of the Mutual Information-based distance matrix
			if parameter.decomposition == 'MDS':
				cmdscale(parameter, _object)											# Computing MDS transformation of distance matrix and plotting TSNE transformation of MDS result
				_object['trans'] = _object['mds_trans']
				kmeans(parameter, _object)												# Applying bootstraps of KMeans on various dimensions of MDS matrix 
			elif parameter.decomposition == 'LPL':
				laplacian(parameter, _object)
				_object['trans'] = _object['lpl_trans']
				kmeans(parameter, _object)												# Applying bootstraps of KMeans on various dimensions of MDS matrix 
			elif parameter.decomposition == 'PCA':
				pca(parameter, _object)
				_object['trans'] = _object['pca_trans']
				kmeans(parameter, _object)												# Applying bootstraps of KMeans on various dimensions of MDS matrix 
			elif parameter.decomposition == 'LPCA':
				laplacian_pca(parameter, _object)
				_object['trans'] = _object['lpca_trans']
				kmeans(parameter, _object)												# Applying bootstraps of KMeans on various dimensions of MDS matrix 
			elif parameter.decomposition == 'LPCA2':
				_object['bs'] = numpy.array(_object['bs']) * 2
				double_bs = _object['bs']
				laplacian(parameter, _object)
				_object['trans'] = _object['lpl_trans']
				kmeans(parameter, _object)
				_object['bs'] = numpy.array(_object['bs']) + 1
				double_bs = double_bs + _object['bs']
				pca(parameter, _object)
				_object['trans'] = _object['pca_trans']
				kmeans(parameter, _object)
				_object['bs'] = double_bs
				laplacian_pca(parameter, _object)
				_object['trans'] = _object['lpca_trans']
			consensus_clustering(parameter, _object)									# Averaging all KMeans results and applying hierarchical clustring			
			if parameter.intermediate_plotting:
				heatmap_clusters(parameter, _object)									# Plotting the final membership matrix 
			if parameter.mode == 'opt-k-2':
				optimal_k_2(parameter, _object)
		elif parameter.mode == 'validate':												# Validation mode: If the clustering result and ground truth is available, checks the result with the ground truth
			parameter.intermediate_plotting = False 			
			heatmap(parameter, _object)													# Hierarchical transformation of the distance matrix
			if parameter.decomposition == 'MDS':
				cmdscale(parameter, _object)											# Computing MDS transformation of hierarchical-transformed distance matrix
				_object['trans'] = _object['mds_trans']
			elif parameter.decomposition == 'LPL':
				laplacian(parameter, _object)
				_object['trans'] = _object['lpl_trans']
			elif parameter.decomposition == 'PCA':
				pca(parameter, _object)
				_object['trans'] = _object['pca_trans']
			elif parameter.decomposition == 'LPCA' or parameter.decomposition == 'LPCA2':
				laplacian_pca(parameter, _object)
				_object['trans'] = _object['lpca_trans']
			cluster_quality(parameter, _object)											# Validate and measure the clustering quality with ground truth
		elif parameter.mode == 'b-quality':												# Quality type B mode: Runs the clustering mode multiple time and returns the average clustering quality measurements
			parameter.intermediate_plotting = False
			it2 = parameter.run
			it = _object['dim_use']
			bs = _object['bs']
			kn = _object['kn']
			bquality = []
			heatmap(parameter, _object)													# Hierarchical transformation of the distance matrix
			if parameter.decomposition == 'MDS':
				cmdscale(parameter, _object)											# Computing MDS transformation of hierarchical-transformed distance matrix
				_object['trans'] = _object['mds_trans']
			elif parameter.decomposition == 'LPL':
				laplacian(parameter, _object)
				_object['trans'] = _object['lpl_trans']
			elif parameter.decomposition == 'PCA':
				pca(parameter, _object)
				_object['trans'] = _object['pca_trans']
			elif parameter.decomposition == 'LPCA':
				laplacian_pca(parameter, _object)
				_object['trans'] = _object['lpca_trans']
			elif parameter.decomposition == 'LPCA2':
				for i in range(it2):
					_object['bs'] = numpy.array(_object['bs']) * 2
					double_bs = _object['bs']
					laplacian(parameter, _object)
					_object['trans'] = _object['lpl_trans']
					kmeans(parameter, _object)
					_object['bs'] = numpy.array(_object['bs']) + 1
					double_bs = double_bs + _object['bs']
					pca(parameter, _object)
					_object['trans'] = _object['pca_trans']
					kmeans(parameter, _object)
					_object['bs'] = double_bs
					laplacian_pca(parameter, _object)
					_object['trans'] = _object['lpca_trans']
					consensus_clustering(parameter, _object)
					bquality.append(cluster_bc_quality(parameter, _object))
			if parameter.decomposition != 'LPCA2':
				for i in range(it2):														# Clustering with random seed KMeans multiple times and report the clustering quality with error bar, here the randomness of consensus clustering comes from the dimension range that is used by default
					_object['bs'] = bs													# Setting the bootstrap to a single number; it was used as seed for KMeans in the earlier stage
					kmeans(parameter, _object)												# Running KMeans on multiple dimensions
					consensus_clustering(parameter, _object)								# Building the consensus cluster from the KMeans results
					bquality.append(cluster_bc_quality(parameter, _object))					# Measuring the clustering quality 
			_object['bs'] = bs 															# Setting the bootstrap to the original value
			cluster_bquality_plot(bquality) 											# Plotting the bar chart with error bar of the clustering quality measurements
		elif parameter.mode == 'c-quality':												# Quality type C mode: Runs the clustering mode across different number of dimensions and returns the average clustering quality measurements for each dimension
			parameter.intermediate_plotting = False
			it = _object['dim_use']
			bs = _object['bs']
			kn = _object['kn']
			it2 = parameter.run
			cquality = []
			heatmap(parameter, _object)													# Hierarchical transformation of the distance matrix
			if parameter.decomposition == 'MDS':
				cmdscale(parameter, _object)											# Computing MDS transformation of hierarchical-transformed distance matrix
				_object['trans'] = _object['mds_trans']
			elif parameter.decomposition == 'LPL':
				laplacian(parameter, _object)
				_object['trans'] = _object['lpl_trans']
			elif parameter.decomposition == 'PCA':
				pca(parameter, _object)
				_object['trans'] = _object['pca_trans']
			elif parameter.decomposition == 'LPCA':
				laplacian_pca(parameter, _object)
				_object['trans'] = _object['lpca_trans']
			elif parameter.decomposition == 'LPCA2':
				for i in it:																# Clustering with one dimension over a range of selected dimensions
					ccquality = []
					_object['dim_use'] = [i]
					_object['d'] = 1
					_object['bs'] = bs
					for k in range(it2):													# Clustering based on one dimension multiple times to measure the confidence intervals of the clustering result for each dimension
						_object['bs'] = numpy.array(_object['bs']) * 2
						double_bs = _object['bs']
						laplacian(parameter, _object)
						_object['trans'] = _object['lpl_trans']
						kmeans(parameter, _object)											# Running KMeans on single dimension, here the randomness of consensus clustering comes from running KMeans for each dimension with mutilple random seed
						_object['bs'] = numpy.array(_object['bs']) + 1
						double_bs = double_bs + _object['bs']
						pca(parameter, _object)
						_object['trans'] = _object['pca_trans']
						kmeans(parameter, _object)
						_object['bs'] = double_bs
						laplacian_pca(parameter, _object)
						_object['trans'] = _object['lpca_trans']
						consensus_clustering(parameter, _object)							# Building the consensus clustering based on one selected dimension, and running KMeans with random seeds
						ccquality.append(cluster_bc_quality(parameter, _object))			# Measure the clustering quality 
					cquality.append(ccquality)
			if parameter.decomposition != 'LPCA2':
				for i in it:																# Clustering with one dimension over a range of selected dimensions
					ccquality = []
					_object['dim_use'] = [i]
					_object['d'] = 1
					_object['bs'] = bs
					for k in range(it2):													# Clustering based on one dimension multiple times to measure the confidence intervals of the clustering result for each dimension
						kmeans(parameter, _object)											# Running KMeans on single dimension, here the randomness of consensus clustering comes from running KMeans for each dimension with mutilple random seed
						consensus_clustering(parameter, _object)							# Building the consensus clustering based on one selected dimension, and running KMeans with random seeds
						ccquality.append(cluster_bc_quality(parameter, _object))			# Measure the clustering quality 
					cquality.append(ccquality)												
			cluster_cquality_plot(cquality, it)											# Plotting the box-plot for all dimensions
		elif parameter.mode == 'opt-k':													# Optimal K Analysis: Based on Silhouette score of various KMeans with different Ks, returns all optimal Ks for clustering
			parameter.intermediate_plotting = False
			heatmap(parameter, _object)													# Hierarchical transformation of the distance matrix
			if parameter.decomposition == 'MDS':
				cmdscale(parameter, _object)											# Computing MDS transformation of hierarchical-transformed distance matrix
				_object['trans'] = _object['mds_trans']
			elif parameter.decomposition == 'LPL':
				laplacian(parameter, _object)
				_object['trans'] = _object['lpl_trans']
			elif parameter.decomposition == 'PCA':
				pca(parameter, _object)
				_object['trans'] = _object['pca_trans']
			elif parameter.decomposition == 'LPCA' or parameter.decomposition == 'LPCA2':
				laplacian_pca(parameter, _object)
				_object['trans'] = _object['lpca_trans']
			optimal_k(parameter, _object)												# Plotting optimal k graph based on silhouette score and KMeans
		elif parameter.mode == 'overlay':												# Overlay mode: Gene expression overlaying and marker gene analysis
			parameter.intermediate_plotting = False
			heatmap(parameter, _object)													# Hierarchical transformation of the distance matrix
			if parameter.decomposition == 'MDS':
				cmdscale(parameter, _object)											# Computing MDS transformation of hierarchical-transformed distance matrix
				_object['trans'] = _object['mds_trans']
			elif parameter.decomposition == 'LPL':
				laplacian(parameter, _object)
				_object['trans'] = _object['lpl_trans']
			elif parameter.decomposition == 'PCA':
				pca(parameter, _object)
				_object['trans'] = _object['pca_trans']
			elif parameter.decomposition == 'LPCA' or parameter.decomposition == 'LPCA2':
				laplacian_pca(parameter, _object)
				_object['trans'] = _object['lpca_trans']
			overlay(parameter, _object)													# Overlaying marker genes expressions on the clustering results, plotting violin plot and p-value analysis to show biological meaningfulness of the clustering results
	# Program will be running on LSF -> Not Supported
	elif parameter.host['host'] == 'LSF':
		print '[INFO] --> [MAIN] Running on LSF cluster not supported.'

if __name__ == '__main__':
	run_time = time.time()
	parameter = parm.Parameter(sys.argv[1:])										# Load parameters provided in the command line
	check = parameter.validate()													# Check if the required parameters are provided they are valid
	if check:
		parameter.print_parameters()
		run(parameter)																# Run MICA
		print 'DONE'
	else:
		print 'EXIT'
	print '[INFO] --> [MAIN] Run time: ' + str(time.time() - run_time)