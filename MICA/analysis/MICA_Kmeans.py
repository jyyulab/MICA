#!/hpcf/authorized_apps/rhel7_apps/mica/openmpi/mica_env/bin/python
# coding: utf-8

# # Use MPI and Scalapy to distribute all of MICA workflow to work on multiple nodes

#load all necessary libraries onto each rank
from mpi4py import MPI
from scipy.sparse import csr_matrix
import scipy.linalg as la
import sys
import numba
import pandas as pd
import scanpy as sc
import scipy as sci
import numpy as np
import anndata
import time
#from sklearn.decomposition import PCA
#import fast_histogram
import logging
logging.basicConfig(level=logging.INFO)
import os
import collections.abc as collections
import socket
import gc
import math
from memory_profiler import profile
from MICA.lib import utils

#========================================================================
#Global variables
comm = MPI.COMM_WORLD

#============== function definitions ========================
# Post dimensionality reduction Clustering and Visualization

import itertools
import argparse
from multiprocessing import Pool
from functools import partial
#from sklearn import cluster        # for kmeanmerge_dist_matss
#from sklearn import manifold       # for tsne
import pickle

#@profile
def clustering_dist(in_file, dr, k, n_bootstrap, out_name, plot_method, umap_min_dist, tsne_perplexity, plot_dim, n_processes, dim_km):
    rank=comm.Get_rank()
    start_total_time = time.time()
    dim_km = map(int, dim_km) #converts string to int 
    start_km_time = time.time()

    #run in parallel
    result=[]
    #Read in completed results if they exist
    kmeans_file=out_name+".kmeanslist"+".pkl"
    if os.path.exists(kmeans_file) : #check for kmeans list object file
        if rank==0:
            print("reading kmeans results from file: "+kmeans_file,flush=True)        
            #read computed list of results
            inp = open(kmeans_file,'rb')  
            result=pickle.load(inp)
    else: #compute new results
        if rank==0:
            print("computing kmeans with %s bootstraps" % (n_bootstrap),flush=True)        
        result = utils.cluster_method_distributed(in_file, n_cluster=k, n_bootstraps=n_bootstrap,
                         common_name=out_name, dims=dim_km, num_processes=n_processes)
        if rank==0:#write to file
            print("cluster_method_multiprocess %s seconds" % (time.time() - start_km_time),flush=True)
            #save list of dataframes (kmeans) as object to file
            outp = open(kmeans_file,'wb')
            pickle.dump(result, outp, pickle.HIGHEST_PROTOCOL)

    start_agg_time = time.time()

    #This function should currently only be run on a single rank
    #Graph cluster consensus of Kmeans results
    #Currently loops over n_bootstraps kmeans results. 
    # Could be parallelized to process n_bootstraps results simultaneously
    if rank==0:
        agg, out_f = utils.consensus_sc3_graph(result, k, out_name)
        print("consensus clustering time %s seconds" % (time.time() - start_agg_time),flush=True)

    #Compute tSNE to visualize clusters
    if rank==0:
        start_tsne_time = time.time()
        #This function should currently only be run on a single rank
        # calls mctsne which is multithreaded
        utils.visualization_dist(agg,  # consensus clustering result
                        in_file,  # reduced_mi_file
                        dr,  # transformation
                        out_f,  # output file name
                        max_dim=plot_dim,
                        visualize=plot_method.lower(),
                        min_dist=umap_min_dist,
                        perplexity=tsne_perplexity,
                        )
        print("tsne %s seconds" % (time.time() - start_tsne_time),flush=True)
        print("--- %s seconds ---" % (time.time() - start_total_time),flush=True)

#=========================================== end funtion definintions ===============================

def main():
    ## Check to make sure MPI (mpi4py) is working
    #comm = MPI.COMM_WORLD #defined globally
    size = comm.Get_size()
    nranks=size
    rank = comm.Get_rank()
    name = MPI.Get_processor_name()

    ## Begin execution of code
    cwd=os.getcwd()
    if rank==0:
        print(cwd,flush=True)

    ## Parse command line arguments
    input_file_name=""
    project_name=""
    data_file_path=""
    nbootstraps=0
    nclusters=0
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("--input", action="store", type="string")
    parser.add_option("--project", action="store", type="string", default="default_proj")
    parser.add_option("--outdir", action="store", type="string", default="outdir")
    parser.add_option("--bootstraps", action="store", type="int", default=1)
    parser.add_option("--clusters", action="store", type="int", default=4)
    (options, args) = parser.parse_args()

    input_file_name = options.input
    project_name = options.project
    data_file_path = options.outdir
    nbootstraps=options.bootstraps
    nclusters=options.clusters

    ## define file paths and names
    data_file_path= data_file_path +'/'+ project_name +'/'
    generated_file_path = data_file_path + 'generated_files/'
    plot_file_path = data_file_path + 'plots/'
    output_file_prefix = generated_file_path+project_name

    ## create directories if necessary
    if rank==0:
        if not os.path.exists(generated_file_path):
            os.makedirs(generated_file_path)
        if not os.path.exists(plot_file_path):
            os.makedirs(plot_file_path)

    if rank==0:
        print("Begin Clustering",flush=True)

    #ceb read reduced matrix
    reduced_mi_filename = output_file_prefix+'_mi_reduced.h5'
    plot_file_name = plot_file_path+project_name
    out_file_name = generated_file_path+project_name
    if os.path.isfile(reduced_mi_filename):
        block_start = time.time()
        clustering_dist(reduced_mi_filename, 
           "mds",#dr 
           nclusters,    #k
           nbootstraps,  #test  #n_bootstrap
           #plot_file_name, #outfile name
           out_file_name, #outfile name
           "tsne", #plot method
           #"umap", #plot method
           0.1,    #umap_min_dist
           30,     #tsne_perplexity
           19,     #plot_dim
           comm.size,     #n_processes
           [19])   #dim_km
        block_end = time.time()
        if rank==0:
            print("Clustering Elapsed = %s" % (block_end - block_start),flush=True)
    else:
        if rank==0:
            print("Reduced similarity matrix file not found. Exiting",flush=True)

    #==============================================================
    comm.Barrier()
    #gracefully end program
    MPI.Finalize()
    exit()


#==============================================================
if __name__ == "__main__":
    main()


