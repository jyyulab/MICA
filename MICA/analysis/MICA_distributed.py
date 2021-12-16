#!/usr/bin/env python
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
from sklearn.decomposition import PCA
#import fast_histogram
import logging
logging.basicConfig(level=logging.INFO)
from MICA.lib import utils
from scalapy import *
from scalapy import blacs
from scalapy import core
import scalapy.routines as rt
import os
import collections.abc as collections
import socket
import gc
import math

#ceb 
from memory_profiler import profile


#========================================================================

#============== function definitions ========================

#@profile
def obtain_distributed_MI_matrix(nslices, generated_file_path, project_name, g_nrows, ncols):
    ## Read in anndata preprocessed files (in distributed mode, by node number) and calculate distance metrics between all row pairs

    #create a 2d list to hold blocks of similarity matrix
    #this should be stored in a distributed scalapack matrix
    b=nslices #row blocks

    if rank==0:
        print("Computing distance metric matrix (MI)",flush=True)

    n_jobs_per_rank= math.ceil(int((b * (b + 1)) / 2)/comm.Get_size())
    #check that n_jobs_per_rank >= 1
    assert(n_jobs_per_rank >= 1), print("Error, Discretization too small. must have at least 1 job per rank.")
    if rank==0:
        print("n_jobs_per_rank=",n_jobs_per_rank,flush=True)


    SM = [[None for j in range(b)] for i in range(b)]

    start = time.time()
    utils.calc_distance_metric_distributed(generated_file_path, project_name, g_nrows, ncols, nslices, SM)
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("calc distance metric Elapsed = %s" % (end - start),flush=True)

    # ## Create distributed matrix for scalapack and copy distributed blocks into object
    # ### This matrix needs to be dense for use in scalapack functions, so we will copy the symmetric data into both upper and lower triangular sections of the MI matrix
    block_start = time.time()
    # ## copy lower triangular transpose to upper triangular for diagonal blocks

    ##from scipy.sparse import csr_matrix #may not use csr as it complicates copy to distributed scalapack and is not used in scalapack apparently
    for i in range(b):
        for j in range(i,b):
            if isinstance(SM[i][j], collections.Iterable):
                if i==j: #copy lower triangular transpose to upper triangular
                    for ii in range(SM[i][j].shape[0]):
                        for jj in range(ii+1,SM[i][j].shape[1]):
                            (SM[i][j])[ii,jj]=(SM[i][j])[jj,ii]

    #copy SM data into global distributed matrix and then write to file?
    #then we can read that file into the Scalapack block cyclic matrix form
    start=time.time()
    #create empty distributed MI matrix
    dMI=core.DistributedMatrix(global_shape=[g_nrows,g_nrows],dtype=np.float64)
    end=time.time()
    comm.Barrier()
    if rank==0:
        print("Create empty distributed MI matrix Elapsed = %s" % (end - start),flush=True)

    #=====================================================================
    ## Copy each SM block submatrix to distributed block cyclic matrix
    start=time.time()
    blocksize=slice_size

    for i in range(b):
        for j in range(i,b): 
            idx = int(i * b + j - (i * (i + 1)) / 2)
            srank = int(idx//n_jobs_per_rank)
            #create dummy block
            lA=np.zeros(shape=(2,2))
            s_block_shape=np.shape(lA)
            if isinstance(SM[i][j], collections.Iterable):
                lA=SM[i][j]
                s_block_shape=np.shape(lA)
            #broadcast sending ranks block shape to all
            s_block_shape = comm.bcast(s_block_shape, root=srank)
            dMI.np2self(lA, srow=i*blocksize, scol=j*blocksize, block_shape=s_block_shape, rank=srank )
    end=time.time()
    #comm.Barrier()
    if rank==0:
        print("Copy MI block matrices to distributed matrix form Elapsed = %s" % (end - start),flush=True)

    ## copy transpose of blocks to fill upper triangular distributed matrix (needed for scalapack computation)
    start=time.time()
    for i in range(b):
        for j in range(i+1,b): #only upper triangular components 
            idx = int(i * b + j - (i * (i + 1)) / 2)
            srank = int(idx//n_jobs_per_rank)
            #create dummy block
            lA=np.zeros(shape=(2,2))
            s_block_shape=np.shape(lA)
            if isinstance(SM[i][j], collections.Iterable):
                lA=np.transpose(SM[i][j])
                s_block_shape=np.shape(lA)
            #broadcast sending ranks block shape to all
            s_block_shape = comm.bcast(s_block_shape, root=srank)
            #print("rank %s checkpt 6.1" % (rank),flush=True)
            dMI.np2self(lA, srow=j*blocksize, scol=i*blocksize, block_shape=s_block_shape, rank=srank )
            #print("rank %s checkpt 6.2" % (rank),flush=True)

    end=time.time()
    #comm.Barrier()
    if rank==0:
        print("Copy transpose of MI block matrices to distributed matrix form Elapsed = %s" % (end - start),flush=True)

    del SM
    #force garbage collection
    gc.collect()

    block_end = time.time()
    if rank==0:
        print("Copy distributed MI matrix Elapsed = %s" % (block_end - block_start),flush=True)

    ## need to also fill in empty symmetric upper triangular portion
    # Even though this is a symmetric matrix, for further processing, we need to copy block data to rest of matrix

    ## Write distributed MI matrix to file
    ### So we can read this in to Scalapack later on
    start=time.time()
    #Write MI matrix to file
    dMI.to_file(mi_filename)
    end=time.time()
    comm.Barrier()
    if rank==0:
        print("Write distributed MI matrix to file Elapsed = %s" % (end - start),flush=True)
    # End compute MI matrix and write
    return dMI
#=================================================================================



#@profile
def normalize_distributed_matrix(dMI):

    ## Now we need to create a normalization matrix
    ### We start with an empty matrix but add the the diagonal as the first column
    ### Then we multiply by its transpose to get a dense matrix
    block_start = time.time()
    #get global indices for diagonal
    gi, lri, lci = dMI.local_diagonal_indices()

    #create matrix to store diagonal row
    dMI_diag=core.DistributedMatrix.empty_like(dMI)

    dMI_row1=core.DistributedMatrix.empty_like(dMI)

    dgi, dlri, dlci = dMI_diag.local_diagonal_indices()

    dMI_diag.local_array[dlri,dlci]=dMI.local_array[lri,lci]

    ## Create a matrix with ones in the first row and zeros elsewhere

    ri, ci = dMI_row1.indices()
    dMI_row1.local_array[:]= ((ri==0).astype(int)).astype(float)

    start = time.time()
    ## Multiply the matrices to get diagonal values on first row of distributed matrix
    dMI_norm = rt.dot(dMI_row1,dMI_diag)
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("checkpt post dot row1 with diag Elapsed = %s" % (end - start),flush=True)

    #Have other ranks wait until prep_dist has completed
    comm.Barrier()
    del dMI_row1
    del dMI_diag
    gc.collect()

    start = time.time()
    ## Multiply the matrix with its transpose to get a dense matrix for normalization
    dMI_normT = dMI_norm.transpose()
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("checkpt test transpose Elapsed = %s" % (end - start),flush=True)
    start = time.time()
    dMI_norm2 = rt.dot(dMI_normT,dMI_norm)
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("checkpt post dot diag with diag Elapsed = %s" % (end - start),flush=True)

    #Have other ranks wait until prep_dist has completed
    comm.Barrier()
    del dMI_norm
    del dMI_normT
    gc.collect()

    ## Use scalapack to compute distributed GEMM

    ## Compute the square root of the normalization matrix

    #compute sqrt of each element
    dMI_norm_root=core.DistributedMatrix.empty_like(dMI)

    dMI_norm_root.local_array[:] = np.sqrt(dMI_norm2.local_array[:])

    ## Now we can finally compute the norm of the MI matrix
    dMI_normed=core.DistributedMatrix.empty_like(dMI)

    dMI_normed.local_array[:] = dMI.local_array[:] / dMI_norm_root.local_array[:]

    #Have other ranks wait until prep_dist has completed
    comm.Barrier()
    del dMI_norm2
    del dMI_norm_root
    gc.collect()

    #ceb write normed matrix
    dMI_normed.to_file(mi_normed_filename)

    block_end = time.time()
    comm.Barrier()
    if rank==0:
        print("Compute matrix norm Elapsed = %s" % (block_end - block_start),flush=True)

    return dMI_normed
#===============================




#@profile
def process_dissimilarity_matrix(rank, processed_dissimilarity_matrix_filename, dMI, dMI_normed, g_nrows, block_size):
    if rank==0:
        print("Preprocessing Dissimilarity Matrix")

    ## Now compute eigenvalues and eigenvectors of dissimmilarity matrix
    block_start = time.time()

    #convert similarity matrix to dissimilarity matrix
    #df= 1-df
    subblock_start = time.time()

    MDS= core.DistributedMatrix.empty_like(dMI)
    MDS.local_array[:]=1.0-dMI_normed.local_array[:]

    # H = I-Ones/n
    I = core.DistributedMatrix.identity(n=g_nrows)

    Ones= core.DistributedMatrix.empty_like(dMI)

    recipn=1.0/g_nrows
    Ones.local_array[:]=recipn

    H = core.DistributedMatrix.empty_like(dMI)

    H.local_array[:] = I.local_array[:] - Ones.local_array[:]

    #Have other ranks wait until prep_dist has completed
    comm.Barrier()

    #remove I, Ones
    del I
    del Ones
    gc.collect()
    #Have other ranks wait until prep_dist has completed

    # B = -H.dot(MDS**2).dot(H)/2
    negH= core.DistributedMatrix.empty_like(dMI)

    negH.local_array[:]= -H.local_array[:]

    MDS2= core.DistributedMatrix.empty_like(dMI)

    MDS2.local_array[:] = MDS.local_array[:]**2

    #Have other ranks wait until prep_dist has completed
    comm.Barrier()
    del MDS
    gc.collect

    #Have other ranks wait until prep_dist has completed
    comm.Barrier()

    start = time.time()
    if rank==0:
        print("checkpt 10.12.a ",flush=True)
    C = rt.dot(negH,MDS2)
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("checkpt 10.12b [dot(negH,MDS2)] Elapsed = %s" % (end - start),flush=True)    


    start = time.time()
    B = rt.dot(C,H)
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("checkpt 10.13 [dot(C,H)] Elapsed = %s" % (end - start),flush=True)

    B.local_array[:]=B.local_array[:]/2.0

    #Have other ranks wait until prep_dist has completed
    comm.Barrier()
    del H
    del C
    del negH
    del MDS2
    del dMI
    gc.collect()

    #Have other ranks wait until prep_dist has completed
    comm.Barrier()

    #dMI_norm=dMI_diag.T*dMI_diag
    #dMI_norm = rt.dot(dMI_diag,dMI_diag,transA='T')

    subblock_end = time.time()
    comm.Barrier()
    if rank==0:
        print("Prep dissimilarity matrix Elapsed = %s" % (subblock_end - subblock_start),flush=True)

    #ceb output dissimilarity matrix here prior to eigenvalue calculation as a checkpoint
    B.to_file(processed_dissimilarity_matrix_filename)

    return B
#================================





#@profile
def compute_eigenvales(reduced_mi_filename, B, g_nrows, rank):
    if rank==0:
        print("Computing Eigenvalues of Preprocessed Dissimilarity Matrix")

    block_start = time.time()

    #compute eigh(B,)
    #we want to pick out the top 200 eigenvalues/vectors from the matrix
    start = time.time()
    n= g_nrows
    evals, dZd = rt.eigh(B,eigvals=(n - np.min([n, 200]), n - 1))
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("checkpt 10.15 [eigh(B,eigvals)] Elapsed = %s" % (end - start),flush=True)

    #copy evecs to root
    evecs = dZd.to_global_array(rank=0)

    block_end = time.time()
    comm.Barrier()
    if rank==0:
        print("Compute Eigenvalues: Elapsed = %s" % (block_end - block_start),flush=True)

    ## gather the top 200 eigenvalues on a single rank
    if rank==0:
        #if os.path.isfile(reduced_mi_filename):
        #    Y = pd.read_hdf(reduced_mi_filename)
        #else:

        ## Read in original dataframe to get labels to attach to results
        #get index names from original dataframe
        adf=utils.read_anndata_file(input_file_name)
        index=adf.obs.index

        ## Postprocess the eigenvalues by sorting and removing negative vals
        idx = np.argsort(evals)[::-1]
        print(len(idx))
        evals = evals[idx]
        evecs = evecs[:, idx]
        evals_pos = evals > 0
        L = np.diag(np.sqrt(evals[evals_pos]))
        V = evecs[:, evals_pos]
        Y = pd.DataFrame(
            data=V.dot(L),
            index=index, #need to reattach index names to eigenvectors
            columns=["mds_" + str(x) for x in np.arange(1, L.shape[0] + 1)],
        )

        ## Write reduced data to file
        Y.to_hdf(reduced_mi_filename, "mds")  # save reduced mi in mds
        print(reduced_mi_filename+" written to disk")
        return Y
    else:
        return None

#====================================================================================







#        Post dimensionality reduction Clustering and Visualization

import itertools
import argparse
from multiprocessing import Pool
from functools import partial

from sklearn import cluster        # for kmeanmerge_dist_matss
from sklearn import manifold       # for tsne


import pickle
#@profile
def clustering_dist(in_file, dr, k, n_bootstrap, out_name, plot_method, umap_min_dist, tsne_perplexity, plot_dim, n_processes, dim_km):

    start_total_time = time.time()
    dim_km = map(int, dim_km) #converts string to int 
    start_km_time = time.time()

    #run in parallel
    result=[]
    #need to be able to read in completed results from outfile if exists
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





## Check to make sure MPI (mpi4py) is working
comm = MPI.COMM_WORLD
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
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--input", action="store", type="string")
parser.add_option("--project", action="store", type="string")
parser.add_option("--outdir", action="store", type="string")
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

#set block size for distributed cyclic matrix
block_size=64

#begin whole code timer
total_time_start = time.time()

#comm.Barrier()
#exit()


#==================================================================================

## Run Prep_dist() to split file into slices ================================
start = time.time()
g_nrows=0 #global number of rows (cells)
ncols=0
nslices=0
slice_size=0
#Run prep.py only on one processor to create the slice files
if rank==0:
    #g_nrows, ncols, nslices = utils.prep_dist(input_file_name, output_file_name, slice_size)
    g_nrows, ncols, nslices, slice_size = utils.prep_dist(input_file_name, output_file_prefix, nranks)
end = time.time()

sys.stdout.flush()#ceb
comm.Barrier()

if rank==0:
    print("prep_dist completed Elapsed = %s" % (end - start),flush=True)

#Have other ranks wait until prep_dist has completed
#comm.Barrier()

#broadcast resultant variables from root to the other ranks
g_nrows = comm.bcast(g_nrows, root=0)
ncols = comm.bcast(ncols, root=0)
nslices = comm.bcast(nslices, root=0)
slice_size = comm.bcast(slice_size, root=0)
##=============================================================================

#comm.Barrier()
#exit()

b=nslices
if slice_size<block_size:
    block_size=slice_size

if rank==0:
    print("global nrows, ncols, slices: ",g_nrows, ncols, nslices,flush=True)

#comm.Barrier()
#exit()#ceb

#Define process grid with process rows and process cols
PR,PC = utils.compute_process_grid_dimensions(comm.size)
if rank==0:
   print("PR=",PR," PC=",PC," block_size=",block_size,flush=True)
comm.Barrier()
#sets default context and block_shape
#@profile
core.initmpi([PR, PC],block_shape=[block_size,block_size])

#comm.Barrier()
#exit()#ceb

#=========================================================================
#if MI file already exists, then read it
mi_filename = output_file_prefix+'_mi_distributed.scalapack'

#==================================================================================

if os.path.isfile(mi_filename):
    ## The following code snippet reads MI matrix from a file and loads it into a distributed Scalapack matrix
    #Read MI matrix from file
    if rank==0:
        print("reading computed MI matrix from file: ",mi_filename)
    dMI=core.DistributedMatrix.from_file(mi_filename, global_shape=[g_nrows,g_nrows], dtype=np.float64, block_shape=[block_size,block_size])
else: #file does not exist, so compute MI matrix and write
    dMI = obtain_distributed_MI_matrix(nslices, generated_file_path, project_name, g_nrows, ncols)

## End compute MI matrix and write

#comm.Barrier()
#exit()#ceb

#==================================================================================

#if MI file already exists, then read it
mi_normed_filename = output_file_prefix+'_mi_normed_distributed.scalapack'

if os.path.isfile(mi_normed_filename):
    start=time.time()
    dMI_normed=core.DistributedMatrix.from_file(mi_normed_filename, global_shape=[g_nrows,g_nrows], dtype=np.float64, block_shape=[block_size,block_size])
    end=time.time()
    if rank==0:
        print("Read distributed normed MI matrix from file Elapsed = %s" % (end - start),flush=True)

else:

    dMI_normed=normalize_distributed_matrix(dMI)


#    ## Now we need to create a normalization matrix
#    ### We start with an empty matrix but add the the diagonal as the first column
#    ### Then we multiply by its transpose to get a dense matrix

#==================================================================================

#Read eigenvalue matrix here if exists, else compute it
#if MI file already exists, then read it
processed_dissimilarity_matrix_filename = output_file_prefix+'_dissimilarity_matrix_distributed.scalapack'
if os.path.isfile(processed_dissimilarity_matrix_filename):
    if rank==0:
        print("Reading preprocessed dissimilarity matrix from file: "+processed_dissimilarity_matrix_filename)
    B=core.DistributedMatrix.from_file(processed_dissimilarity_matrix_filename, global_shape=[g_nrows,g_nrows], dtype=np.float64, block_shape=[block_size,block_size])
else:

    B=process_dissimilarity_matrix(rank, processed_dissimilarity_matrix_filename, dMI, dMI_normed, g_nrows, block_size)

#==================================================================================

#ceb read or compute reduced matrix
reduced_mi_filename = output_file_prefix+'_mi_reduced.h5'

if os.path.isfile(reduced_mi_filename):
    if rank==0:
        print("Reading existing eigenvalue file")
        Y = pd.read_hdf(reduced_mi_filename)
else:
    Y = compute_eigenvales(reduced_mi_filename, B, g_nrows, rank)

#print reduced data to screen
if rank==0:
    print("Reduced Data Matrix")
    print(Y)
#==================================================================================

plot_file_name = plot_file_path+project_name
out_file_name = generated_file_path+project_name


#run entire clustering pipeline
block_start = time.time()

if rank==0:
    print("Begin Clustering",flush=True)

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

#==============================================================
block_end = time.time()
comm.Barrier()
if rank==0:
    print("Clustering Elapsed = %s" % (block_end - block_start),flush=True)

total_time_end = time.time()
if rank==0:
    print("Total Time Elapsed = %s" % (total_time_end - total_time_start),flush=True)

#gracefully end program
exit()




