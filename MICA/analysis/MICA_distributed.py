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
import sys
import numba
import pandas as pd
import scanpy as sc
import scipy as sci
import numpy as np
import anndata
import time
from sklearn.decomposition import PCA
import fast_histogram
import logging
logging.basicConfig(level=logging.INFO)
from MICA.lib import utils
from scalapy import *
from scalapy import blacs
from scalapy import core
import scalapy.routines as rt
import os
import scipy.linalg as la
import collections
import socket
import gc
import math

## Check to make sure MPI (mpi4py) is working

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()
print("MPI Check: {host}[{pid}]: {rank}/{size}".format(
    host=socket.gethostname(),
    pid=os.getpid(),
    rank=comm.rank,
    size=comm.size,
),flush=True)


## Begin execution of code
cwd=os.getcwd()
if rank==0:
    print(cwd)

SMALL_TEST=False

if SMALL_TEST:    
    data_file_path = cwd+'/test_data/inputs/10x/PBMC/3k/pre-processed/'
    input_file_name = data_file_path + 'pbmc3k_preprocessed.h5ad'
    project_name = 'pbmc3k'
else:
    data_file_path = cwd+'/test_data/inputs/80k/'
    input_file_name = data_file_path + '/ecoCART_double_construct_bothdonor_MICAfilt_except_week19.h5ad'
    project_name = 'ecoCART'


output_file_name = data_file_path+project_name

#We should try to automate or suggest the appropriate slice size based on the number of cells and number of ranks
    #g_nrows, ncols, nslices = prep_dist(input_file_name, output_file_name, slice_size)
#set slice size (max size of row blocks)mpirun_rsh -n 4 -export-all -hostfile $LSB_DJOB_HOSTFILE MV2_IBA_HCA=mlx5_0 MV2_ON_DEMAND_THRESHOLD=1 MV2_NUM_HCAS=1 MV2_DEFAULT_TIME_OUT=230 MV2_USE_SRQ=0 MV2_ENABLE_AFFINITY=0 MV2_NDREG_ENTRIES_MAX=32000 MV2_NDREG_ENTRIES=16000 python ./MICA/analysis/MICA_distributed.py


#ceb
#maybe determine slice size by the number or ranks available
#need to understand whay some values are causing segfault failures




#there is a minimum slice size for selected ranks.
#not sure what this value is
if SMALL_TEST :
    #slice_size = 32
    #slice_size = 50
    #slice_size = 64
    #slice_size = 100
    #slice_size = 120
    #slice_size = 128 #min for 3480? Why?
    slice_size = 512 #min for 3480? ok rank=16 1 job per rank?
    #slice_size = 1024 #min for 3480? ok rank=8 1 job per rank?
    #slice_size = 1000 #min for 3480? ok rank=8 1. FAIL. Slice size must be power of 2?
    block_size=32
else:
    #slice_size = 16385 #min for 50497 ?
    #slice_size = 8192 #ok for 50497?
    #slice_size = 12500 #ok for 50497 with ranks=16 (failure at 64 ranks)
    #slice_size = 12500 #ok for 50497 with ranks=32, 1 job per rank?
    slice_size = 6250 #ok for 50497 with ranks=64, 1 job per rank?
    #block_size=32
    #block_size=100
    block_size=64

if rank==0:
    print (input_file_name)



#ceb should test to see if 
#read initial file here to get dimensions
#then check to see if mi file exists.
#if so we can skip the prep_dist 



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
    print("Rank:%s prep_dist completed Elapsed = %s" % (rank,end - start),flush=True)
    
#Have other ranks wait until prep_dist has completed
comm.Barrier()

#broadcast resultant variables from root to the other ranks
g_nrows = comm.bcast(g_nrows, root=0)
ncols = comm.bcast(ncols, root=0)
nslices = comm.bcast(nslices, root=0)
b=nslices    

if rank==0:
    print("global nrows, ncols, slices: ",g_nrows, ncols, nslices)


#=========================================================================
#if MI file already exists, then read it
mi_filename = data_file_path+project_name+'_mi_distributed.scalapack'

if os.path.isfile(mi_filename): 

    #Define process grid with process rows and process cols

    #ideally we would like BR and BC to the square root of the num_ranks to get a square process matrix
    PR=int(np.sqrt(comm.size))
    PC=PR
    #if we can't create a square matrix, get next best dimensions
    if PR*PR!=size:
        PC=size//PR
    #if that doesn't match then create a 1D process array
    if PR*PC != comm.size:
        PR=size
        PC=1

    #ceb Check to make sure that PR*PC == num_ranks
    assert(PR*PC == comm.size), print("Error, Size of process grid must match total number of ranks.")

    #we want the greatest common factors here for efficiency
    if rank==0:
        print("PR=",PR, "PC=",PC,flush=True)

    start=time.time()
    #sets default context and block_shape
    core.initmpi([PR, PC],block_shape=[block_size,block_size])
    end=time.time()
    #print("rank:%s checkpt initmpi scalapack context Elapsed = %s" % (rank,end - start),flush=True)

    start=time.time()
    #create empty distributed MI matrix
    #dMI=core.DistributedMatrix(global_shape=[g_nrows,g_nrows],dtype=np.float64)
    end=time.time()
    comm.Barrier()
    if rank==0:
        print("rank:%s checkpt create empty distributed MI matrix Elapsed = %s" % (rank,end - start),flush=True)


    ## The following code snippet reads MI matrix from a file and loads it into a distributed Scalapack matrix
    start=time.time()
    #Read MI matrix from file
    #mi_filename = data_file_path+project_name+'_mi_distributed.scalapack'
    if rank==0:
        print("reading computed MI matrix from file: ",mi_filename)
    dMI=core.DistributedMatrix.from_file(mi_filename, global_shape=[g_nrows,g_nrows], dtype=np.float64, block_shape=[block_size,block_size])
    end=time.time()
    comm.Barrier()
    #print("Rank:%s Read distributed MI matrix from file Elapsed = %s" % (rank,end - start),flush=True)

    #print(dir(dMI))
    #print("rank=",rank," dMI_local: ",dMI.local_array,flush=True)

#============================================================================
else: #file does not exist, so compute MI matrix and write

    start = time.time()

    ## Read in anndata preprocessed files (in distributed mode, by node number) and calculate distance metrics between all row pairs

    #create a 2d list to hold blocks of similarity matrix
    #this should be stored in a distributed scalapack matrix
    b=nslices #row blocks


    SM = [[None for j in range(b)] for i in range(b)] 

    start = time.time()
    utils.calc_distance_metric_distributed(data_file_path, project_name, g_nrows, ncols, nslices, SM)
    end = time.time()
    comm.Barrier()

    if rank==0:
        print("Rank:%s calc distance metric Elapsed = %s" % (rank,end - start),flush=True)



    # ## Create distributed matrix for scalapack and copy distributed blocks into object
    # ### This matrix needs to be dense for use in scalapack functions, so we will copy the symmetric data into both upper and lower triangular sections of the MI matrix
    block_start = time.time()
    # ## copy lower triangular transpose to upper triangular for diagonal blocks

    start = time.time()
    ##from scipy.sparse import csr_matrix #may not use csr as it complicates copy to distributed scalapack and is not used in scalapack apparently
    for i in range(b):
        for j in range(i,b):
            if isinstance(SM[i][j], collections.Iterable):
                if i==j: #copy lower triangular transpose to upper triangular 
                    for ii in range(SM[i][j].shape[0]):
                        for jj in range(ii+1,SM[i][j].shape[1]):
                            (SM[i][j])[ii,jj]=(SM[i][j])[jj,ii]
                            #print("Rank:",rank, " SM[",i,"][",j,"]=",SM[i][j])
    end=time.time()
    comm.Barrier()
    if rank==0:
        print("rank:%s checkpt copy to upper triangluar for symmetric matrix Elapsed = %s" % (rank,end - start),flush=True)

    # ## Populate a global array with all of the MI data from each rank
    # 
    # Preferably, comm.Barrier()we would like each rank to contribute of their block MI matrices to the global matrix,
    # but currently the distributed global matrix has to be constructed from a global (not distributed) array

    #copy SM data into global distributed matrix and then write to file?

    #then we can read that file into the Scalapack block cyclic matrix form




    #Define process grid with process rows and process cols

    #ideally we would like BR and BC to the square root of the num_ranks to get a square process matrix
    PR=int(np.sqrt(comm.size))
    PC=PR
    comm.Barrier()
    #if we can't create a square matrix, get next best dimensions
    if PR*PR!=size:
        PC=size//PR
    #if that doesn't match then create a 1D process array
    if PR*PC != comm.size:
        PR=size
        PC=1

    #ceb Check to make sure that PR*PC == num_ranks
    assert(PR*PC == comm.size), print("Error, Size of process grid must match total number of ranks.") 

    #we want the greatest common factors here for efficiency
    if rank==0:
        print("PR=",PR, "PC=",PC,flush=True)

    start=time.time()
    #sets default context and block_shape
    core.initmpi([PR, PC],block_shape=[block_size,block_size])
    end=time.time()
    #print("rank:%s checkpt initmpi scalapack context Elapsed = %s" % (rank,end - start),flush=True)

    start=time.time()
    #create empty distributed MI matrix
    dMI=core.DistributedMatrix(global_shape=[g_nrows,g_nrows],dtype=np.float64)
    end=time.time()
    comm.Barrier()
    if rank==0:
        print("rank:%s checkpt create empty distributed MI matrix Elapsed = %s" % (rank,end - start),flush=True)








    ## Copy each SM block submatrix to distributed block cyclic matrix
    start=time.time()
    blocksize=slice_size
    n_jobs_per_rank= math.ceil(int((b * (b + 1)) / 2)/comm.Get_size())

    #ceb check that n_jobs_per_rank >= 1
    assert(n_jobs_per_rank >= 1), print("Error, Discretization too small. must have at least 1 job per rank.") 

    if rank==0:
        print("rank=",rank," n_jobs_per_rank=",n_jobs_per_rank,flush=True)

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
                #print("rank: ",rank," copy SM[",i,j,"] shape: ",s_block_shape," to global i,j:",i*blocksize,j*blocksize,flush=True)

            #broadcast sending ranks block shape to all
            s_block_shape = comm.bcast(s_block_shape, root=srank) 
            #print("rank: ",rank," post comm.bcast i,j:",i,j," srank=",srank,flush=True)
            #print("rank: ",rank," copy SM[",i,j,"] shape: ",s_block_shape," to global i,j:",i*blocksize,j*blocksize,flush=True)
            #ceb
            dMI.np2self(lA, srow=i*blocksize, scol=j*blocksize, block_shape=s_block_shape, rank=srank )      
            #print("rank: ",rank," post np2self ",flush=True)

    end=time.time()
    comm.Barrier()
    if rank==0:
        print("rank:%s checkpt copy MI block matrices to distributed matrix form Elapsed = %s" % (rank,end - start),flush=True)

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
        print("rank:%s checkpt copy transpose of MI block matrices to distributed matrix form Elapsed = %s" % (rank,end - start),flush=True)

    del SM
    gc.collect()

    block_end = time.time()
    if rank==0:
        print("Rank:%s Copy distributed MI matrix Elapsed = %s" % (rank,block_end - block_start),flush=True)

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
        print("Rank:%s Write distributed MI matrix to file Elapsed = %s" % (rank,end - start),flush=True)

    ### The following code snippet reads MI matrix from a file and loads it into a distributed Scalapack matrix
    #istart=time.time()
    ##Read MI matrix from file
    #dMI=core.DistributedMatrix.from_file(mi_filename, global_shape=[g_nrows,g_nrows], dtype=np.float64, block_shape=[block_size,block_size])
    #end=time.time()
    #comm.Barrier()

    print("rank=",rank," dMI_local: ",dMI.local_array,flush=True)

    if rank==0:
        print("Rank:%s Read distributed MI matrix from file %s Elapsed = %s" % (rank, mi_filename, end - start),flush=True)

##------------------------------------------------------------------------------------------------------



### The following code snippet reads MI matrix from a file and loads it into a distributed Scalapack matrix
#istart=time.time()
##Read MI matrix from file
#mi_filename = data_file_path+project_name+'_mi_distributed.scalapack'
#dMI.from_file(mi_filename, global_shape=[g_nrows,g_nrows], dtype=np.float64, block_shape=[block_size,block_size])
#end=time.time()
#comm.Barrier()
#if rank==0:
#    print("Rank:%s Read distributed MI matrix from file Elapsed = %s" % (rank,end - start),flush=True)




#if MI file already exists, then read it
mi_normed_filename = data_file_path+project_name+'_mi_normed_distributed.scalapack'

if os.path.isfile(mi_normed_filename): 

    dMI_normed=core.DistributedMatrix.from_file(mi_normed_filename, global_shape=[g_nrows,g_nrows], dtype=np.float64, block_shape=[block_size,block_size])

else:
    ## Now we need to create a normalization matrix
    ### We start with an empty matrix but add the the diagonal as the first column
    ### Then we multiply by its transpose to get a dense matrix
    block_start = time.time()
    start = time.time()
    #get global indices for diagonal
    gi, lri, lci = dMI.local_diagonal_indices()
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("rank:%s checkpt post get diagonal indicies Elapsed = %s" % (rank,end - start),flush=True)

    start = time.time()
    #create matrix to store diagonal row
    dMI_diag=core.DistributedMatrix.empty_like(dMI)
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("rank:%s checkpt post create empty matrix for diag Elapsed = %s" % (rank,end - start),flush=True)

    start = time.time()
    dMI_row1=core.DistributedMatrix.empty_like(dMI)
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("rank:%s checkpt post create empty matrix for row1 Elapsed = %s" % (rank,end - start),flush=True)

    start = time.time()
    dgi, dlri, dlci = dMI_diag.local_diagonal_indices()
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("rank:%s checkpt post get local indices Elapsed = %s" % (rank,end - start),flush=True)

    start = time.time()
    dMI_diag.local_array[dlri,dlci]=dMI.local_array[lri,lci]
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("rank:%s checkpt post get local array Elapsed = %s" % (rank,end - start),flush=True)
    #dMI_diag.local_array[0,dlci]=dMI.local_array[lri,lci]
    #my_diag[comm.rank]=dMI.local_array[lri,lci]


    ## Create a matrix with ones in the first row and zeros elsewhere

    start = time.time()
    ri, ci = dMI_row1.indices()
    dMI_row1.local_array[:]= ((ri==0).astype(int)).astype(float)
    end = time.time()
    #print(dMI_row1.local_array)
    comm.Barrier()
    if rank==0:
        print("rank:%s checkpt post convert array to float Elapsed = %s" % (rank,end - start),flush=True)

    start = time.time()
    ## Multiply the matrices to get diagonal values on first row of distributed matrix
    dMI_norm = rt.dot(dMI_row1,dMI_diag)
    #dMI_norm = dMI_diag.dot(dMI_row1)
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("rank:%s checkpt post dot row1 with diag Elapsed = %s" % (rank,end - start),flush=True)

    #Have other ranks wait until prep_dist has completed
    comm.Barrier()
    del dMI_row1
    del dMI_diag
    gc.collect()

    start = time.time()
    ## Multiply the matrix with its transpose to get a dense matrix for normalization
    #CEB taking long time, probably due to transpose
    #dMI_norm=dMI_diag.T*dMI_diag
    #dMI_normT = dMI_norm.copy()
    dMI_normT = dMI_norm.transpose()
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("rank:%s checkpt test transpose Elapsed = %s" % (rank,end - start),flush=True)
    start = time.time()
    #dMI_norm2 = rt.dot(dMI_norm,dMI_norm,transA='T')
    dMI_norm2 = rt.dot(dMI_normT,dMI_norm)
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("rank:%s checkpt post dot diag with diag Elapsed = %s" % (rank,end - start),flush=True)

    #Have other ranks wait until prep_dist has completed
    comm.Barrier()
    del dMI_norm
    del dMI_normT
    gc.collect()

    ## Use scalapack to compute distributed GEMM

    ## Compute the square root of the normalization matrix

    #compute sqrt of each element
    start = time.time()
    dMI_norm_root=core.DistributedMatrix.empty_like(dMI)
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("rank:%s checkpt post create empty squared matrix Elapsed = %s" % (rank,end - start),flush=True)

    start = time.time()
    dMI_norm_root.local_array[:] = np.sqrt(dMI_norm2.local_array[:])
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("rank:%s checkpt post compute square of matrix Elapsed = %s" % (rank,end - start),flush=True)

    ## Now we can finally compute the norm of the MI matrix
    start = time.time()
    dMI_normed=core.DistributedMatrix.empty_like(dMI)
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("rank:%s checkpt post create empty norm matrix Elapsed = %s" % (rank,end - start),flush=True)

    start = time.time()
    dMI_normed.local_array[:] = dMI.local_array[:] / dMI_norm_root.local_array[:]
    #ceb
    #dMI_normed.local_array[:] = dMI.local_array[:] / 1.0 #ceb test
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("rank:%s checkpt post divide local array by norm array Elapsed = %s" % (rank,end - start),flush=True)

    #Have other ranks wait until prep_dist has completed
    comm.Barrier()
    del dMI_norm2
    del dMI_norm_root
    gc.collect()

    start = time.time()
    #mi_normed_filename = data_file_path+project_name+'_mi_normed_distributed.scalapack'
    dMI_normed.to_file(mi_normed_filename)
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("rank:%s checkpt post write normed MI matrix to file Elapsed = %s" % (rank,end - start),flush=True)

    block_end = time.time()
    comm.Barrier()
    if rank==0:
        print("rank: %s, Compute matrix norm Elapsed = %s" % (rank,block_end - block_start),flush=True)
#================================================================================







## Now compute eigenvalues and eigenvectors of dissimmilarity matrix

block_start = time.time()
n= g_nrows

#convert similarity matrix to dissimilarity matrix
#df= 1-df
subblock_start = time.time()

start = time.time()
MDS= core.DistributedMatrix.empty_like(dMI)
#end = time.time()
#print("rank: %s, checkpt 10.1 Elapsed = %s" % (rank,end - start),flush=True)
#start = time.time()
MDS.local_array[:]=1.0-dMI_normed.local_array[:]
end = time.time()
comm.Barrier()
#if rank==0:
#    print("rank: %s, checkpt 10.2 Elapsed = %s" % (rank,end - start),flush=True)

# H = I-Ones/n
start = time.time()
I= core.DistributedMatrix.identity(n=g_nrows)
end = time.time()
comm.Barrier()
#if rank==0:
#    print("rank: %s, checkpt 10.3 Elapsed = %s" % (rank,end - start),flush=True)

start = time.time()
Ones= core.DistributedMatrix.empty_like(dMI)
end = time.time()
comm.Barrier()
#if rank==0:
#    print("rank: %s, checkpt 10.4 Elapsed = %s" % (rank,end - start),flush=True)

start = time.time()
recipn=1.0/n
Ones.local_array[:]=recipn
end = time.time()
comm.Barrier()
#if rank==0:
#    print("rank: %s, checkpt 10.5 Elapsed = %s" % (rank,end - start),flush=True)

start = time.time()
H = core.DistributedMatrix.empty_like(dMI)
end = time.time()
comm.Barrier()
#if rank==0:
#    print("rank: %s, checkpt 10.6 Elapsed = %s" % (rank,end - start),flush=True)

start = time.time()
H.local_array[:] = I.local_array[:] - Ones.local_array[:]
end = time.time()
comm.Barrier()
#if rank==0:
#    print("rank: %s, checkpt 10.7 Elapsed = %s" % (rank,end - start),flush=True)

#Have other ranks wait until prep_dist has completed
comm.Barrier()

#remove I, Ones
del I
del Ones
gc.collect()
#Have other ranks wait until prep_dist has completed
comm.Barrier()

# B = -H.dot(MDS**2).dot(H)/2
start = time.time()
negH= core.DistributedMatrix.empty_like(dMI)
end = time.time()
comm.Barrier()
#if rank==0:
#    print("rank: %s, checkpt 10.8 Elapsed = %s" % (rank,end - start),flush=True)

start = time.time()
negH.local_array[:]= -H.local_array[:]
end = time.time()
comm.Barrier()
#if rank==0:
#    print("rank: %s, checkpt 10.9 Elapsed = %s" % (rank,end - start),flush=True)

start = time.time()
MDS2= core.DistributedMatrix.empty_like(dMI)
end = time.time()
comm.Barrier()
if rank==0:
    print("rank: %s, checkpt 10.10 Elapsed = %s" % (rank,end - start),flush=True)

start = time.time()
MDS2.local_array[:] = MDS.local_array[:]**2
end = time.time()
comm.Barrier()
if rank==0:
    print("rank: %s, checkpt 10.11 Elapsed = %s" % (rank,end - start),flush=True)

#Have other ranks wait until prep_dist has completed
comm.Barrier()
del MDS
gc.collect

#Have other ranks wait until prep_dist has completed
comm.Barrier()

start = time.time()
#ceb test
#C = core.DistributedMatrix.empty_like(dMI)#not necessary to create space here
if rank==0:
    print("rank: %s, checkpt 10.12.a " % (rank),flush=True)
C = rt.dot(negH,MDS2)
end = time.time()
comm.Barrier()
if rank==0:
    print("rank: %s, checkpt 10.12b Elapsed = %s" % (rank,end - start),flush=True)

start = time.time()
#ceb test
#B = core.DistributedMatrix.empty_like(dMI)#not necessary to create space here
B = rt.dot(C,H)
end = time.time()
comm.Barrier()
if rank==0:
    print("rank: %s, checkpt 10.13 Elapsed = %s" % (rank,end - start),flush=True)

start = time.time()
B.local_array[:]=B.local_array[:]/2.0
end = time.time()
comm.Barrier()
if rank==0:
    print("rank: %s, checkpt 10.14 Elapsed = %s" % (rank,end - start),flush=True)

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
    print("rank: %s, Prep dissimilarity matrix Elapsed = %s" % (rank,subblock_end - subblock_start),flush=True)


#ceb breaks after this point ----------

#start = time.time()

#compute eigh(B,)
#we want to pick out the top 200 eigenvalues/vectors from the matrix

start = time.time()
evals, dZd = rt.eigh(B,eigvals=(n - np.min([n, 200]), n - 1)) 
end = time.time()
comm.Barrier()
if rank==0:
    print("rank: %s, checkpt 10.15 Elapsed = %s" % (rank,end - start),flush=True)

#copy evecs to root
start = time.time()
evecs = dZd.to_global_array(rank=0)
end = time.time()
comm.Barrier()
if rank==0:
    print("rank: %s, checkpt 10.16 Elapsed = %s" % (rank,end - start),flush=True)
#gZd = dZd.to_global_array(rank=0)

block_end = time.time()
comm.Barrier()
if rank==0:
    print("rank: %s, Compute Eigenvalues: Elapsed = %s" % (rank,block_end - block_start),flush=True)  







## gather the top 200 eigenvalues on a single rank

## Read in original dataframe to get labels to attach to results

#get index names from original dataframe
if rank==0:
    #data_file_path = cwd+'/test_data/inputs/10x/PBMC/3k/pre-processed/'
    #input_file_name = data_file_path + 'pbmc3k_preprocessed.h5ad'
    adf=utils.read_anndata_file(input_file_name)
    index=adf.obs.index

    ## Postprocess the eigenvalues by sorting and removing negative vals

    #if rank==0:
    idx = np.argsort(evals)[::-1]
    print(len(idx))
    
    evals = evals[idx]
    evecs = evecs[:, idx]
    evals_pos = evals > 0
    L = np.diag(np.sqrt(evals[evals_pos]))
    V = evecs[:, evals_pos]
    #print(V)
    
    Y = pd.DataFrame(
        data=V.dot(L),
        index=index, #need to reattach index names to eigenvectors
        columns=["mds_" + str(x) for x in np.arange(1, L.shape[0] + 1)],
    )
    
    #Y.to_hdf(out_file_name + "_reduced.h5", "mds")  # save reduced mi in mds

    ## Write reduced data to file

    #if rank==0:
    data_file_path = cwd+'/test_data/inputs/10x/PBMC/3k/pre-processed/'
    out_file_name = data_file_path + 'pbmc3k_preprocessed'    
    Y.to_hdf(out_file_name + "_reduced.h5", "mds")  # save reduced mi in mds

    #if rank==0:
    print(Y)


#%%px
#perplexity=30
#max_dim=200
#    Y.to_hdf(out_file_name + "_reduced.h5", "mds")  # save reduced mi in mds
#
#    if print_plot == "True":
#        vis = tsne(
#            Y,
#            max_dim,
#            out_file_name,
#            "mds",
#            perplexity,
#            print_plot,
#        )
#        vis.to_hdf(out_file_name + "_reduced", "mds_tsne")  # save preview in key "mds_tsne"






