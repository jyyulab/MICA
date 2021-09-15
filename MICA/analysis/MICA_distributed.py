#!/usr/bin/env python
# coding: utf-8

# # Use MPI and Scalapy to distribute all of MICA workflow to work on multiple nodes

# prep.py component \
# Read input file and slice into manageable sizes

# # Launch an ipython parallel cluster
# Run this on hpc node to launch a cluster with mpi engines

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
import collections
import socket
import gc
import math

## Check to make sure MPI (mpi4py) is working

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()
#print("MPI Check: {host}[{pid}]: {rank}/{size}".format(
#    host=socket.gethostname(),
#    pid=os.getpid(),
#    rank=comm.rank,
#    size=comm.size,),flush=True)


## Begin execution of code
cwd=os.getcwd()
if rank==0:
    print(cwd,flush=True)

#SMALL_TEST=False
SMALL_TEST=True

if SMALL_TEST:    
    project_name = 'pbmc3k'
    data_file_path = "/home/cburdysh/MICA_Project/MICA_distributed/MICA/test_data/pbmc3k/"
    input_file_name = data_file_path + 'input/' + 'pbmc3k_preprocessed.h5ad'
    generated_file_path = data_file_path + 'generated_files/'
    plot_file_path = data_file_path + 'plots/'
else:
    project_name = 'ecoCART'
    data_file_path = "/home/cburdysh/MICA_Project/MICA_distributed/MICA/test_data/ecoCART/"
    input_file_name = data_file_path + 'input/' + 'ecoCART_double_construct_bothdonor_MICAfilt_except_week19.h5ad'
    generated_file_path = data_file_path + 'generated_files/'
    plot_file_path = data_file_path + 'plots/'


output_file_name = generated_file_path+project_name

#We should try to automate or suggest the appropriate slice size based on the number of cells and number of ranks
    #g_nrows, ncols, nslices = prep_dist(input_file_name, output_file_name, slice_size)
#set slice size (max size of row blocks)mpirun_rsh -n 4 -export-all -hostfile $LSB_DJOB_HOSTFILE MV2_IBA_HCA=mlx5_0 MV2_ON_DEMAND_THRESHOLD=1 MV2_NUM_HCAS=1 MV2_DEFAULT_TIME_OUT=230 MV2_USE_SRQ=0 MV2_ENABLE_AFFINITY=0 MV2_NDREG_ENTRIES_MAX=32000 MV2_NDREG_ENTRIES=16000 python ./MICA/analysis/MICA_distributed.py


#ceb
#maybe determine slice size by the number or ranks available
#need to understand whay some values are causing segfault failures


#ceb
#b = slices = nrows/slice_size
# N Comparisons = b(b+1)/2
# Njobs/rank = Ncomparisons/nranks (ideally 1/rank)
# So we want ncomparisons==nranks 
# nranks = b(b+1)/2


# Maybe, this might not be best for eigenval calculation but should be irrevelant as we can use (mostly) arbitrary ranks to create a process grid
 

#there is a minimum slice size for selected ranks.
#not sure what this value is
if SMALL_TEST :
    #slice_size = 32
    #slice_size = 50
    #slice_size = 64
    #slice_size = 100
    #slice_size = 120
    #slice_size = 128 #min for 3480? Why?
    #slice_size = 512 #min for 3480? ok rank=16 1 job per rank?
    slice_size = 256 #min for 3480? ok rank=32 1 job per rank?
    #slice_size = 1024 #min for 3480? ok rank=8 1 job per rank?
    #slice_size = 1000 #min for 3480? ok rank=8 1. FAIL. Slice size must be power of 2?
    #block_size=32
    block_size=8
else:
    slice_size = 1000 #min for 50497 ?
    #slice_size = 8192 #ok for 50497?
    #slice_size = 12500 #ok for 50497 with ranks=16 (failure at 64 ranks)
    #slice_size = 12500 #ok for 50497 with ranks=32, 1 job per rank?
    if size==32:
        slice_size = 12500 #ok for 50497 with ranks=64, 1 job per rank?
    if size==64:
        slice_size = 6250 #ok for 50497 with ranks=64, 1 job per rank?
    if size==128:
        slice_size = 3125 # for 50497 with ranks=128, 1 job per rank?
    if size==256:
        #slice_size = 3125 # for 50497 with ranks=128, 1 job per rank?
        #slice_size = 1500 # for 50497 with ranks=256, 1 job per rank?
        #slice_size = 1000 # for 50497 with ranks=256, 1 job per rank?
        slice_size = 512 # for 50497 with ranks=256, 1 job per rank?
        #slice_size = 197 # for 50497 with ranks=256, 1 job per rank?

    #Maybe divide rows by ranks to make sure we have a diagonal block for each rank
    #block_size=32
    #block_size=100
    block_size=64




if rank==0:
    print (input_file_name,flush=True)

#ceb should test to see if 
#read initial file here to get dimensions
#then check to see if mi file exists.
#if so we can skip the prep_dist 


total_time_start = time.time()


# ## Run Prep_dist() to split file into slices
start = time.time()
#Run prep.py only on one processor to create the slice files
g_nrows=0 #global number of rows (cells)
ncols=0
nslices=0
if rank==0: 
    g_nrows, ncols, nslices = utils.prep_dist(input_file_name, output_file_name, slice_size)
end = time.time()
comm.Barrier()

if rank==0:
    print("prep_dist completed Elapsed = %s" % (end - start),flush=True)
    
#Have other ranks wait until prep_dist has completed
comm.Barrier()

#broadcast resultant variables from root to the other ranks
g_nrows = comm.bcast(g_nrows, root=0)
ncols = comm.bcast(ncols, root=0)
nslices = comm.bcast(nslices, root=0)
b=nslices    

if rank==0:
    print("global nrows, ncols, slices: ",g_nrows, ncols, nslices)



#function to compute process grid rows and colums
def compute_process_grid_dimensions(nranks):
    #Define process grid with process rows and process cols
    #ideally we would like BR and BC to the square root of the num_ranks to get a square process matrix
    PR=int(np.sqrt(nranks))
    PC=PR
    #if we can't create a square matrix, get next best dimensions
    if PR*PR!=size:
        PC=size//PR
    #if that doesn't match then create a 1D process array
    if PR*PC != nranks:
        PR=size
        PC=1

    #Special cases that are not optimal or valid under the general rules
    #ultimately we want to compute the squares matrix possible

    #ceb best for 32
    # Need to write code to select the greatest 2 factors of comm.size for optimal process grid.
    if comm.size==32:
        #ceb eigenvalue solve fails (MPI errors) at pr=16 pc=8. Not sure why
        PR=8
        PC=4
    #ceb best for 32
    # Need to write code to select the greatest 2 factors of comm.size for optimal process grid.
    if comm.size==64:
        #ceb eigenvalue solve fails (MPI errors) at pr=16 pc=8. Not sure why
        PR=8
        PC=8
    #ceb best for 128
    # Need to write code to select the greatest 2 factors of comm.size for optimal process grid.
    if nranks==128:
        #PR=128
        #PC=1
        #ceb eigenvalue solve fails (MPI errors) at pr=16 pc=8. Not sure why
        PR=64
        PC=2
    #getting error at PR=16 PC=16
    if nranks==256:
        PR=16
        PC=16
    #ceb Check to make sure that PR*PC == num_ranks
    assert(PR*PC == nranks), print("Error, Size of process grid must match total number of ranks.")
    #we want the greatest common factors here for efficiency
    if rank==0:
        print("PR=",PR, "PC=",PC,flush=True)
    return PR,PC
    


#=========================================================================
#if MI file already exists, then read it
#mi_filename = data_file_path+project_name+'_mi_distributed.scalapack'
mi_filename = output_file_name+'_mi_distributed.scalapack'

if os.path.isfile(mi_filename): 

    #Define process grid with process rows and process cols
    PR,PC = compute_process_grid_dimensions(comm.size)

    #sets default context and block_shape
    core.initmpi([PR, PC],block_shape=[block_size,block_size])

    ## The following code snippet reads MI matrix from a file and loads it into a distributed Scalapack matrix
    #Read MI matrix from file
    #mi_filename = data_file_path+project_name+'_mi_distributed.scalapack'
    if rank==0:
        print("reading computed MI matrix from file: ",mi_filename)
    dMI=core.DistributedMatrix.from_file(mi_filename, global_shape=[g_nrows,g_nrows], dtype=np.float64, block_shape=[block_size,block_size])

#============================================================================
else: #file does not exist, so compute MI matrix and write


    ## Read in anndata preprocessed files (in distributed mode, by node number) and calculate distance metrics between all row pairs

    #create a 2d list to hold blocks of similarity matrix
    #this should be stored in a distributed scalapack matrix
    b=nslices #row blocks

    if rank==0:
        print("Computing distance metric matrix (MI)",flush=True)

    SM = [[None for j in range(b)] for i in range(b)] 

    start = time.time()
    #utils.calc_distance_metric_distributed(data_file_path, project_name, g_nrows, ncols, nslices, SM)
    utils.calc_distance_metric_distributed(generated_file_path, project_name, g_nrows, ncols, nslices, SM)
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("calc distance metric Elapsed = %s" % (end - start),flush=True)



    # ## Create distributed matrix for scalapack and copy distributed blocks into object
    # ### This matrix needs to be dense for use in scalapack functions, so we will copy the symmetric data into both upper and lower triangular sections of the MI matrix
    block_start = time.time()
    # ## copy lower triangular transpose to upper triangular for diagonal blocks

    #start = time.time()
    ##from scipy.sparse import csr_matrix #may not use csr as it complicates copy to distributed scalapack and is not used in scalapack apparently
    for i in range(b):
        for j in range(i,b):
            if isinstance(SM[i][j], collections.Iterable):
                if i==j: #copy lower triangular transpose to upper triangular 
                    for ii in range(SM[i][j].shape[0]):
                        for jj in range(ii+1,SM[i][j].shape[1]):
                            (SM[i][j])[ii,jj]=(SM[i][j])[jj,ii]
                            #print("Rank:",rank, " SM[",i,"][",j,"]=",SM[i][j])

    #copy SM data into global distributed matrix and then write to file?
    #then we can read that file into the Scalapack block cyclic matrix form

    ##Define process grid with process rows and process cols
    PR,PC = compute_process_grid_dimensions(comm.size)

    #sets default context and block_shape
    core.initmpi([PR, PC],block_shape=[block_size,block_size])

    start=time.time()
    #create empty distributed MI matrix
    dMI=core.DistributedMatrix(global_shape=[g_nrows,g_nrows],dtype=np.float64)
    end=time.time()
    comm.Barrier()
    if rank==0:
        print("checkpt create empty distributed MI matrix Elapsed = %s" % (end - start),flush=True)



    #=====================================================================


    ## Copy each SM block submatrix to distributed block cyclic matrix
    start=time.time()
    blocksize=slice_size
    n_jobs_per_rank= math.ceil(int((b * (b + 1)) / 2)/comm.Get_size())

    #ceb check that n_jobs_per_rank >= 1
    assert(n_jobs_per_rank >= 1), print("Error, Discretization too small. must have at least 1 job per rank.") 

    if rank==0:
        print("n_jobs_per_rank=",n_jobs_per_rank,flush=True)

    comm.Barrier()

    for i in range(b):
        for j in range(i,b): # j in range [i,b]
            idx = int(i * b + j - (i * (i + 1)) / 2)
            srank = int(idx//n_jobs_per_rank)
            #print("rank: ",rank," idx: ",idx," srank: ",srank,flush=True)

            #create dummy block
            lA=np.zeros(shape=(2,2))
            s_block_shape=np.shape(lA)
            #print("rank: ",rank," s_block_shape ",s_block_shape, flush=True)

            if isinstance(SM[i][j], collections.Iterable):
                lA=SM[i][j]
                s_block_shape=np.shape(lA)
                #print("rank: ",rank," copy SM[",i,j,"] to srank: ",srank," shape: ",s_block_shape," to global i,j:",i*blocksize,j*blocksize,flush=True)

            #broadcast sending ranks block shape to all
            s_block_shape = comm.bcast(s_block_shape, root=srank) 

            #print("rank: ",rank," post comm.bcast i,j:",i,j," srank=",srank,flush=True)
            #print("rank: ",rank," copy SM[",i,j,"] to srank:",srank,"shape: ",s_block_shape," to global i,j:",i*blocksize,j*blocksize,flush=True)
            #ceb
            dMI.np2self(lA, srow=i*blocksize, scol=j*blocksize, block_shape=s_block_shape, rank=srank )      
            #print("rank: ",rank," post np2self ",flush=True)
 
    #print("rank: ",rank," post copy SM to scalapack ",flush=True)
    end=time.time()
    comm.Barrier()

    if rank==0:
        print("checkpt copy MI block matrices to distributed matrix form Elapsed = %s" % (end - start),flush=True)




    ## copy transpose of blocks to fill upper triangular distributed matrix (needed for scalapack computation)
    start=time.time()

    for i in range(b):
        for j in range(i+1,b): # j in range [i,b]
            idx = int(i * b + j - (i * (i + 1)) / 2)
            srank = int(idx//n_jobs_per_rank)
            #print("rank: ",rank," idx: ",idx," srank: ",srank,flush=True)
            #create dummy block
            lA=np.zeros(shape=(2,2))
            s_block_shape=np.shape(lA) 

            if isinstance(SM[i][j], collections.Iterable):
                #print(rank," chekcpt 1",flush=True)
                lA=np.transpose(SM[i][j])
                #print(rank," chekcpt 2",flush=True)
                s_block_shape=np.shape(lA)
                #print("copy SM[",j,i,"] shape: ",s_block_shape,flush=True)
                #print("rank: ",rank," copy SM[",i,j,"] shape: ",s_block_shape," to global i,j:",i*blocksize,j*blocksize,flush=True)

            #broadcast sending ranks block shape to all
            #print(rank," chekcpt 3",flush=True)
            s_block_shape = comm.bcast(s_block_shape, root=srank)   
            #print("rank: ",rank," post comm.bcast i,j:",i,j,flush=True)
            #print("rank: ",rank," copy SM[",i,j,"] shape: ",s_block_shape," to global i,j:",i*blocksize,j*blocksize,flush=True)
            #This line is causing mpi errors
            #ceb
            dMI.np2self(lA, srow=j*blocksize, scol=i*blocksize, block_shape=s_block_shape, rank=srank )      
            #print(rank," chekcpt 4",flush=True)
            #print("rank: ",rank," post np2self ",flush=True)

    end=time.time()
    comm.Barrier()
    if rank==0:
        print("checkpt copy transpose of MI block matrices to distributed matrix form Elapsed = %s" % (end - start),flush=True)

    del SM
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



#============================================================================
# End compute MI matrix and write



#if MI file already exists, then read it
mi_normed_filename = output_file_name+'_mi_normed_distributed.scalapack'

if os.path.isfile(mi_normed_filename): 
    start=time.time()
    dMI_normed=core.DistributedMatrix.from_file(mi_normed_filename, global_shape=[g_nrows,g_nrows], dtype=np.float64, block_shape=[block_size,block_size])
    end=time.time()
    if rank==0:
        print("Read distributed normed MI matrix from file Elapsed = %s" % (end - start),flush=True)

else:
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



#================================================================================




#ceb
#Read eigenvalue matrix here if exists, else compute it





#if MI file already exists, then read it
processed_dissimilarity_matrix_filename = output_file_name+'_dissimilarity_matrix_distributed.scalapack'
if os.path.isfile(processed_dissimilarity_matrix_filename): 
    if rank==0:
        print("Reading preprocessed dissimilarity matrix from file: "+processed_dissimilarity_matrix_filename)
    B=core.DistributedMatrix.from_file(processed_dissimilarity_matrix_filename, global_shape=[g_nrows,g_nrows], dtype=np.float64, block_shape=[block_size,block_size])
else:
    if rank==0:
        print("Preprocessing Dissimilarity Matrix")

    ## Now compute eigenvalues and eigenvectors of dissimmilarity matrix
    block_start = time.time()

    #convert similarity matrix to dissimilarity matrix
    #df= 1-df
    subblock_start = time.time()

    MDS= core.DistributedMatrix.empty_like(dMI)
    #print("rank: %s, checkpt 10.1 Elapsed = %s" % (rank,end - start),flush=True)
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
##------------------------------------------------------------------------------------------------------




#ceb read or compute reduced matrix
reduced_mi_filename = output_file_name+'_mi_reduced.h5'

if os.path.isfile(reduced_mi_filename): 
    if rank==0:
        print("Reading existing eigenvalue file")
        Y = pd.read_hdf(reduced_mi_filename)
else:
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

        #data_file_path = cwd+'/test_data/inputs/10x/PBMC/3k/pre-processed/'
        #input_file_name = data_file_path + 'pbmc3k_preprocessed.h5ad' #from initial file
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


if rank==0:
    print(Y)




#====================================================================================
#        Post dimensionality reduction Clustering and Visualization

import itertools
#import time
import argparse
#import pandas as pd
from multiprocessing import Pool
from functools import partial
#from MICA.lib import utils

#if rank==0:
#    reduced_mi_matrix = Y
#    #reduced_mi_matrix = pd.read_hdf(reduced_mi_filename)

from sklearn import cluster        # for kmeanmerge_dist_matss
from sklearn import manifold       # for tsne


def kmeans(in_mat, n_clusters, project_name, dim, bootstrap_id):
    out_file_name = project_name + "_kmeans_k" + str(n_clusters) + "_d" + str(dim) + ".h5.tmp." + str(bootstrap_id)
    km = cluster.KMeans(n_clusters=n_clusters, max_iter=1000, n_init=1000)
    km_res = pd.DataFrame(
        data=np.transpose(km.fit_predict(in_mat.iloc[:, 0:dim])),
        index=in_mat.index,
        columns=["label"],
    )
    # km_res.to_hdf(out_file_name, "kmeans")
    print("Executing kmeans:" + out_file_name,flush=True)
    return km_res


#import faiss
#https://github.com/facebookresearch/faiss
def faiss_kmeans(in_mat, n_clusters, project_name, dim, bootstrap_id):
    out_file_name = project_name + "_faiss_k" + str(n_clusters) + "_d" + str(dim) + ".h5.tmp." + str(bootstrap_id)
    max_iter=500
    n_init=1
    verbosity=False
    x=in_mat.values
    x=np.ascontiguousarray(x, dtype=np.float32)
    km = faiss.Kmeans( x.shape[1] , n_clusters, niter=max_iter, nredo=n_init, verbose=verbosity)
    km.train(x)
    D, membership = km.index.search(x,1)
    #km_res = pd.DataFrame(
    #    data=I,
    #    index=in_mat.index,
    #    columns=["label"],
    #)
    # km_res.to_hdf(out_file_name, "kmeans")
    print("Executing faiss kmeans:" + out_file_name,flush=True)
    #return km_res
    return np.array(membership)


#shared memory kmeans 
def cluster_method_multiprocess(mi_file, n_cluster, n_iter, common_name, dims=[19], num_processes=1):
   #Create a thread pool of num_processes size
    pool = Pool(processes=num_processes)
    hdf = pd.HDFStore(mi_file)
    r = []
    #ceb not sure what these "keys" are
    # but kmeans is being run independently on each key
    for trans in hdf.keys():
        print("mc clustering trans=[%s]" % (trans),flush=True)
        df = hdf[trans]
        method_iterable = partial(utils.kmeans, df, n_cluster, common_name)
        #method_iterable = partial(faiss_kmeans, df, n_cluster, common_name)
        #n_iter is the number of bootstrap cases we want to execute for kmeans
        iterations = itertools.product(dims, range(n_iter))
        print("cluster_method_multiprocess checkpt 3",flush=True)
        #ceb run n_iter kmeans functions on thread pool 
        res = pool.starmap(method_iterable, iterations)
        r = r + res
    pool.close()
    hdf.close()
    print("cluster_method_multiprocess checkpt 6",flush=True)
    return r


#https://github.com/DmitryUlyanov/Multicore-TSNE
#shared memory kmeans
def cluster_method_distributed(mi_file, n_cluster, n_bootstraps, common_name, dims=[19], num_processes=1):
    #Run on all mpi ranks
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    nranks = comm.Get_size()
    if myrank < n_bootstraps:
        color = 1
        key = myrank
    else:
        color = 2 #MPI.UNDEFINED
        key = myrank

    subcomm = MPI.COMM_WORLD.Split(color,key)
    #print("myrank=%s subcomm=%s" % (myrank,subcomm))
    #need to use n_iter to set number of bootstraps
    root = 0
    r = []
    #only use ranks in subcomm
    if color==1:
        mysubrank = subcomm.Get_rank()
        print("myrank=%s my_subrank= %s subcomm.Get_size()=%s" % (myrank,mysubrank,subcomm.Get_size()),flush=True)
        #All ranks read input file
        hdf = pd.HDFStore(mi_file)
        #ceb not sure what these "keys" are
        # but kmeans is being run independently on each key
        for trans in hdf.keys():
            #print("mc clustering trans=[%s]" % (trans),flush=True)
            df = hdf[trans]
            ncells = df.shape[0]
       
            print("cluster_method_distributed checkpt 0",flush=True)
            #n_iter is the number of bootstrap cases we want to execute for kmeans
            #we will distribute a kmeans to each rank 
            # we will send numpy array clusters to root rank
            #clusters = faiss_kmeans(df, n_cluster, common_name, dim=ncells ,bootstrap_id=myrank)
            clusters = np.asarray((utils.kmeans(df, n_cluster, common_name, dim=ncells ,bootstrap_id=mysubrank)).label,dtype=int)
            print("cluster_method_distributed checkpt 1",flush=True)
      
            clusters_buffer=None 
            if mysubrank == root:
                #clusters_buffer = np.empty( nranks*ncells, dtype=int )
                clusters_buffer = np.empty( n_bootstraps*ncells, dtype=int )

            #print("cluster_method_distributed checkpt 2",flush=True)
            #comm.Gatherv(clusters, clusters_buffer, root=root)
            subcomm.Gatherv(clusters, clusters_buffer, root=root)

            #print("cluster_method_distributed checkpt 3",flush=True)

            #collect clusters data on root rank
            if mysubrank == root:
                #for i in range(nranks):
                for i in range(subcomm.Get_size()):
                    #add dataframe to global list
                    res = pd.DataFrame(
                        data=clusters_buffer[i*ncells:(i+1)*ncells],
                        index=df.index,
                        columns=["label"],
                        ) 
                    #print(res)
                    r.append(res)
        hdf.close()

    #print("rank %s is waiting at barrier" % (myrank),flush=True)
    comm.Barrier()
    if color==1:
        subcomm.Free()
    #print("cluster_method_distributed checkpt 6",flush=True)
    return r



def visualization(
        agg_mtx,  
        reduced_mi_file, 
        transformation,  
        out_file_name,
        max_dim=0,
        visualize="umap", 
        min_dist=0.25, 
        perplexity=30, 
        nprocesses=1,
):
    cclust = agg_mtx
    hdf = pd.HDFStore(reduced_mi_file)
    transformation = "pca" if transformation == "lpca" else transformation
    hdf_trans = hdf[transformation.lower()]
    hdf.close()

    if visualize == "umap":
        embed_2d = utils.umap(hdf_trans, max_dim, min_dist)

    elif visualize == "tsne":
        perplexity = np.min([perplexity, np.max(cclust.groupby(["label"]).size())])
        #embed_2d = utils.tsne(hdf_trans, max_dim, "", None, perplexity, "False")
        embed_2d = utils.mctsne(hdf_trans, max_dim, "", None, perplexity, "False",nprocesses)
    cclust = pd.concat([cclust, embed_2d.loc[cclust.index, :]], axis=1)
    res = cclust.loc[:, ["X", "Y", "label"]]
    # save 2D embedding to txt file
    out_file_name = out_file_name + "_" + visualize
    utils.scatter2(res, out_file_name + '.png')
    res.to_csv(out_file_name + "_ClusterMem.txt", sep="\t")


def consensus_sc3(km_results, n_clusters, common_name=None):
    """ Implement SC3's consensus clustering. (https://www.nature.com/articles/nmeth.4236)
    Args:
        km_results (list of dataframes): each dataframe is a clustering results with two columns (cell_index,
                                         clustering_label)
        n_clusters (int): number of clusters
        common_name (name): common string to name the output files
    Returns:
        Clustering results
    """
    print("consensus checkpt 1",flush=True)
    n_iter = len(km_results)
    print('Number of k-mean results: {}'.format(n_iter))
    if n_iter == 0:
        return None
    if len(km_results) == 0:
        return None

    print("consensus checkpt 1",flush=True)
    conss_binary_mat = np.zeros((km_results[0].shape[0], km_results[0].shape[0]))
    print("consensus checkpt 2",flush=True)

    #loop over list of clusters generated by parallel kmeans.
    #This loop can be done in parallel
    #The concensus matrix will be sized [n_cells x n_cells], so will grow as square of n_cells
    #We are going to perform a hierarchical agglomerative clustering algorithm on this matrix.
    #This will impose a computational (and eventually memory) bottleneck

    for i in range(n_iter):
        arr = km_results[i].to_numpy()
        #compare arr with it's transpose, create matrix of bools where true if .
        # arr[:,None] is transpose of arr
        #creates an NxN matrix of bools
        mask = arr[:, None] == arr
        #convert binary matrix to integers
        #We could convert this new matrix to Anndata CSR
        binary_mat = mask[:,:,0].astype(int)

        #sum values where kmeans clusters agree
        conss_binary_mat += binary_mat

    print("consensus checkpt 3",flush=True)
    #divide summations by total number of clusterings compared
    conss__mat = conss_binary_mat / n_iter
    print("consensus checkpt 4",flush=True)
    #ceb differs from previous aggregate function which uses ward linkage
    clust = cluster.AgglomerativeClustering(
        linkage="complete", n_clusters=n_clusters, affinity="euclidean"
    )

    #this agglomerative clustering function takes most time 
    print("aggregate: fit agglomerative cluster",flush=True)
    clust.fit(conss_mat)
    print("aggregate: finish agglomerative cluster",flush=True)

    cclust = pd.DataFrame(data=clust.labels_, index=km_results[0].index, columns=["label"])
    index = cclust.groupby(["label"]).size().sort_values(ascending=False).index
    label_map = {index[i]: i + 1000 for i in range(len(index))}
    cclust = cclust.replace(to_replace={"label": label_map})
    index = cclust.groupby(["label"]).size().sort_values(ascending=False).index
    map_back = {index[i]: i + 1 for i in range(len(index))}
    cclust = cclust.replace(to_replace={"label": map_back})
    out_file = common_name + "_k" + str(n_clusters)
    out_file_hdf = out_file + "_cclust.h5"
    cclust.to_hdf(out_file_hdf, "cclust")
    adjacency = edgelist2adjacency(edge_list)# conss_binary_mat.to_hdf(out_file_hdf, "membership")
    print('consensus_sc3 is done')
    return cclust, out_file


from scipy import sparse
from sknetwork.utils import edgelist2adjacency, edgelist2biadjacency
from sknetwork.data import convert_edge_list, load_edge_list, load_graphml
from sknetwork.visualization import svg_graph, svg_digraph, svg_bigraph
from sknetwork.clustering import Louvain

def consensus_sc3_graph(km_results, n_clusters, common_name=None):
    """ Implement SC3's consensus clustering. (https://www.nature.com/articles/nmeth.4236)
    Args:
        km_results (list of dataframes): each dataframe is a clustering results with two columns (cell_index,
                                         clustering_label)
        n_clusters (int): number of clusters
        common_name (name): common string to name the output files
    Returns:
        Clustering results
    """
    print("consensus checkpt 1",flush=True)
    n_iter = len(km_results)
    print('Number of k-mean results: {}'.format(n_iter))
    if n_iter == 0:
        return None
    if len(km_results) == 0:
        return None

    print("consensus checkpt 1",flush=True)
    #create an empty csr matrix that we will be adding to
    conss_adjacency_mat = csr_matrix( np.zeros((km_results[0].shape[0], km_results[0].shape[0])) )
    print("consensus checkpt 2",flush=True)

    #loop over list of clusters generated by parallel kmeans.
    #This loop can be done in parallel but this should not be computational costly so should probably leave as is.

    #The concensus matrix will be sized [n_cells x n_cells], so will grow as square of n_cells
    #We are going to perform a hierarchical agglomerative clustering algorithm on this matrix.
    
    #ceb This might be able to be run in parallel, particularly the summation of the csr adjacency matrix
    for i in range(n_iter):
        print("processing kmeans result %s" % i,flush=True)
        #Get array of cluster labels from kmeans solution
        arr = km_results[i].to_numpy()

        #compare arr with it's transpose, create matrix of bools where true if .
        # arr[:,None] is transpose of arr
        #creates an NxN matrix of bools
        #not sure how efficient this is. Really only need the upper triangular portion
        mask = arr[:, None] == arr
        print("consensus checkpt 2.1",flush=True)

        #convert binary matrix to integers (adjacency matrix)
        adjacency_mat = mask[:,:,0].astype(int) 
        print("consensus checkpt 2.2",flush=True)
        #ceb for graph clustering
        #convert this new matrix to Anndata CSR
        adjacency_mat = csr_matrix(adjacency_mat)
        print("consensus checkpt 2.3",flush=True)

        #sum values where kmeans clusters agree
        conss_adjacency_mat += adjacency_mat
        print("consensus checkpt 2.4",flush=True)


    print("consensus checkpt 3",flush=True)
    #Normalize by dividing summations by total number of clusterings compared
    conss_adjacency_mat = conss_adjacency_mat / n_iter
    #print("conss_adjacency_mat",flush=True)
    #print(conss_adjacency_mat)

    if 0:
        print("consensus checkpt 4",flush=True)
        #ceb differs from previous aggregate function which uses ward linkage
        clust = cluster.AgglomerativeClustering( linkage="complete", n_clusters=n_clusters, affinity="euclidean" )
        #this agglomerative clustering function takes most time 
        print("aggregate: fit agglomerative cluster",flush=True)
        clust.fit(conss_adjacency_mat) 
        labels=clust.labels_
        print("consensus labels",flush=True)
        print(labels,flush=True)
        print("aggregate: finish agglomerative cluster",flush=True)
    else:
        print("consensus checkpt 4",flush=True)
        louvain=Louvain()
        labels = louvain.fit_transform(conss_adjacency_mat)
        print("consensus labels",flush=True)
        print(labels,flush=True)

    cclust = pd.DataFrame(data=labels, index=km_results[0].index, columns=["label"])
    index = cclust.groupby(["label"]).size().sort_values(ascending=False).index
    label_map = {index[i]: i + 1000 for i in range(len(index))}
    cclust = cclust.replace(to_replace={"label": label_map})
    index = cclust.groupby(["label"]).size().sort_values(ascending=False).index
    map_back = {index[i]: i + 1 for i in range(len(index))}
    cclust = cclust.replace(to_replace={"label": map_back})
    out_file = common_name + "_k" + str(n_clusters)
    out_file_hdf = out_file + "_cclust.h5"
    cclust.to_hdf(out_file_hdf, "cclust")
    # conss_binary_mat.to_hdf(out_file_hdf, "membership")
    print('consensus_sc3 is done')
    return cclust, out_file


def clustering_dist(in_file, dr, k, n_bootstrap, out_name, plot_method, umap_min_dist, tsne_perplexity, plot_dim, n_processes, dim_km):

    start_total_time = time.time()
    dim_km = map(int, dim_km) #converts string to int 
    start_km_time = time.time()
    #This function should currently only be run on a single rank
    #if rank == 0:
    #    result = cluster_method_multiprocess(in_file, n_cluster=k, n_iter=n_bootstrap,
    #                     common_name=out_name, dims=dim_km, num_processes=n_processes)
    #run in parallel
    #need to be able to read in completed results from outfile if exists
    result = cluster_method_distributed(in_file, n_cluster=k, n_bootstraps=n_bootstrap,
                         common_name=out_name, dims=dim_km, num_processes=n_processes)
    if rank==0:
        print("cluster_method_multiprocess %s seconds" % (time.time() - start_km_time),flush=True)
    start_agg_time = time.time()
    #print("pre create consensus clustering checkpt ",flush=True)
    #agg, out_f = utils.consensus_sc3(result, k, out_name)
    #agg, out_f = consensus_sc3(result, k, out_name)
    #This function should currently only be run on a single rank
    #Graph cluster consensus of Kmeans results
    if rank==0:
        agg, out_f = consensus_sc3_graph(result, k, out_name)
        print("consensus clustering time %s seconds" % (time.time() - start_agg_time),flush=True)
    #Compute tSNE to visualize clusters
    if rank==0:
        start_tsne_time = time.time()
        #This function should currently only be run on a single rank
        visualization(agg,  # consensus clustering result
                        in_file,  # reduced_mi_file
                        dr,  # transformation
                        out_f,  # output file name
                        max_dim=plot_dim,
                        visualize=plot_method.lower(),
                        min_dist=umap_min_dist,
                        perplexity=tsne_perplexity,
                        #nprocesses=n_processes,
                        nprocesses=1,
                        )
        print("tsne %s seconds" % (time.time() - start_tsne_time),flush=True)

        print("--- %s seconds ---" % (time.time() - start_total_time),flush=True)




#==============================================================
#compute tsne
plot_file_name = plot_file_path+project_name+"scatter"


#if 0: #run multiproces clustering
#    #if rank==0:
#    start_km_time = time.time()
#    dim_km=[19]
#    dim_km = map(int, dim_km) #converts string to int 
#    #result = cluster_method_multiprocess(reduced_mi_filename, n_cluster=5, n_iter=2,
#    #                     common_name=plot_file_name, dims=dim_km, num_processes=comm.size)
#    result = cluster_method_distributed(reduced_mi_filename, n_cluster=5, n_iter=2,
#                         common_name=plot_file_name, dims=dim_km, num_processes=comm.size)
#    print("cluster_method_multiprocess %s seconds" % (time.time() - start_km_time),flush=True)




#if 0: #compute tSNE and plot
#    max_dim=200
#    block_start = time.time()
#    if rank==0:
#        #vis is dataframe    
#        vis = utils.tsne( Y, max_dim, plot_file_name, "mds",)
#        #vis = TSNE(Y, max_dim, plot_file_name, "mds", plot="False", n_jobs=comm.size)
#        vis.to_hdf(plot_file_name + "_reduced.h5", "mds_tsne")  # save preview in key "mds_tsne"
#        ##utils.visualization(Y, reduced_mi_filename, "mds", plot_file_name, max_dim, "tsne")
#    block_end = time.time()
#    comm.Barrier()
#    if rank==0:
#        print("Plot TSNE Elapsed = %s" % (block_end - block_start),flush=True)

if 1: #run entire clustering pipeline
    block_start = time.time()

    if rank==0:
        print("Begin Clustering",flush=True)

    clustering_dist(reduced_mi_filename, 
               "mds",#dr 
               5,    #k
               comm.size,  #test  #n_bootstrap
               #16,  #test  #n_bootstrap
               plot_file_name, #outfile name
               "tsne", #plot method
               #"umap", #plot method
               0.1,    #umap_min_dist
               30,     #tsne_perplexity
               19,     #plot_dim
               comm.size,     #n_processes
               [19])   #dim_km

    block_end = time.time()
    comm.Barrier()
    if rank==0:
        print("Clustering Elapsed = %s" % (block_end - block_start),flush=True)

total_time_end = time.time()
if rank==0:
    print("Total Time Elapsed = %s" % (total_time_end - total_time_start),flush=True)

#gracefully end program
exit()




