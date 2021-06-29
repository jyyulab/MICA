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


# ## Check to make sure MPI (mpi4py) is working




import os
import socket
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()
print("{host}[{pid}]: {rank}/{size}".format(
    host=socket.gethostname(),
    pid=os.getpid(),
    rank=comm.rank,
    size=comm.size,
))


# ## Begin execution of code




cwd=os.getcwd()
if rank==0:
    print(cwd)
    
data_file_path = cwd+'/test_data/inputs/10x/PBMC/3k/pre-processed/'
input_file_name = data_file_path + 'pbmc3k_preprocessed.h5ad'
project_name = 'pbmc3k'
output_file_name = data_file_path+project_name






#set slice size (max size of row blocks)
slice_size = 500






if rank==0:
    print (input_file_name)


# ## Run Prep_dist() to split file into slices





#Run prep.py only on one processor to create the slice files
g_nrows=0 #global number of rows (cells)
ncols=0
nslices=0
if rank==0: 
    #g_nrows, ncols, nslices = prep_dist(input_file_name, output_file_name, slice_size)
    g_nrows, ncols, nslices = utils.prep_dist(input_file_name, output_file_name, slice_size)
    
#broadcast resultant variables from root to the other ranks
g_nrows = comm.bcast(g_nrows, root=0)
ncols = comm.bcast(ncols, root=0)
nslices = comm.bcast(nslices, root=0)






if rank==0:
    print("global nrows, ncols, slices: ",g_nrows, ncols, nslices)


# ## Read in anndata preprocessed files (in distributed mode, by node number) and calculate distance metrics between all row pairs
# 





#create a 2d list to hold blocks of similarity matrix
#this should be stored in a distributed scalapack matrix
b=nslices #row blocks
SM = [[None for j in range(b)] for i in range(b)] 

start = time.time()
utils.calc_distance_metric_distributed(data_file_path, project_name, g_nrows, ncols, nslices, SM)
end = time.time()
print("Elapsed = %s" % (end - start))


# ## Create distributed matrix for scalapack and copy distributed blocks into object
# ### This matrix needs to be dense for use in scalapack functions, so we will copy the symmetric data into both upper and lower triangular sections of the MI matrix

# ## copy lower triangular transpose to upper triangular for diagonal blocks





##from scipy.sparse import csr_matrix #may not use csr as it complicates copy to distributed scalapack and is not used in scalapack apparently
import collections
for i in range(b):
    for j in range(i,b):
        if isinstance(SM[i][j], collections.Iterable):
            if i==j: #copy lower triangular transpose to upper triangular 
                for ii in range(SM[i][j].shape[0]):
                    for jj in range(ii+1,SM[i][j].shape[1]):
                        (SM[i][j])[ii,jj]=(SM[i][j])[jj,ii]
                #print("Rank:",rank, " SM[",i,"][",j,"]=",SM[i][j])


# ## Populate a global array with all of the MI data from each rank
# 
# Preferably, we would like each rank to contribute of their block MI matrices to the global matrix,
# but currently the distributed global matrix has to be constructed from a global (not distributed) array




#copy SM data into global distributed matrix and then write to file?

#then we can read that file into the Scalapack block cyclic matrix form






#test to distribute matrix from local blocks rather than global array
from scalapy import blacs
import os
import numpy as np
import scipy.linalg as la
from mpi4py import MPI
from scalapy import core
import scalapy.routines as rt

#distribute MI components to ranks as scalapack distributed matrix
comm = MPI.COMM_WORLD
rank = comm.rank
size = comm.size #total number of ranks

global_num_rows =g_nrows
global_num_cols =g_nrows
local_num_rows =g_nrows/b

block_size=64 #default is 32

#Define process grid with process rows and process cols
#We'll use a 2d process grid to distribute blocks so we want to have num_ranks divisible by 2
assert((size % 2)==0)
#ideally we would like BR and BC to the square root of the num_ranks to get a square process matrix
PR=int(np.sqrt(size))
PC=PR

#if we can't create a square matrix, get next best dimensions
if PR*PR!=size:
    PC=size//PR
if rank==0:
    print("PR=",PR, "PC=",PC)

#sets default context and block_shape
core.initmpi([PR, PC],block_shape=[block_size,block_size])
#convert to fortran array indexing to match scalapack functions
#create global matrix from array on rank0
dMI=core.DistributedMatrix(global_shape=[g_nrows,g_nrows],dtype=np.float64)


# ## Copy each SM block submatrix to distributed block cyclic matrix





blocksize=slice_size
n_jobs_per_rank= (int((b * (b + 1)) / 2))/comm.Get_size()
import collections
for i in range(b):
    for j in range(i,b): # j in range [i,b]
        idx = int(i * b + j - (i * (i + 1)) / 2)
        srank = idx//n_jobs_per_rank
        lA=np.zeros(shape=(2,2))
        s_block_shape=np.shape(lA)
        if isinstance(SM[i][j], collections.Iterable):
            lA=SM[i][j]
            s_block_shape=np.shape(lA)
            print("copy SM[",i,j,"] shape: ",s_block_shape," to global i,j:",i*blocksize,j*blocksize)
        #broadcast sending ranks block shape to all
        s_block_shape = comm.bcast(s_block_shape, root=srank)   
        dMI.np2self(lA, srow=i*blocksize, scol=j*blocksize, block_shape=s_block_shape, rank=srank )      


# ## copy transpose of blocks to fill upper triangular distributed matrix (needed for scalapack computation)




blocksize=slice_size
n_jobs_per_rank= (int((b * (b + 1)) / 2))/comm.Get_size()
import collections

for i in range(b):
    for j in range(i+1,b): # j in range [i,b]
        idx = int(i * b + j - (i * (i + 1)) / 2)
        srank = idx//n_jobs_per_rank
        lA=np.zeros(shape=(2,2))
        s_block_shape=np.shape(lA)
        if isinstance(SM[i][j], collections.Iterable):
            lA=np.transpose(SM[i][j])
            s_block_shape=np.shape(lA)
            #print("copy SM[",j,i,"] shape: ",s_block_shape)
        #broadcast sending ranks block shape to all
        s_block_shape = comm.bcast(s_block_shape, root=srank)   
        dMI.np2self(lA, srow=j*blocksize, scol=i*blocksize, block_shape=s_block_shape, rank=srank )      





## need to also fill in empty symmetric upper triangular portion





# Even though this is a symmetric matrix, for further processing, we need to copy block data to rest of matrix


# ## Write distributed MI matrix to file
# ### So we can read this in to Scalapack later on





#Write MI matrix to file
mi_filename = data_file_path+project_name+'_mi_distributed.scalapack'
dMI.to_file(mi_filename)


# ## The following code snippet reads MI matrix from a file and loads it into a distributed Scalapack matrix





#Read MI matrix from file
mi_filename = data_file_path+project_name+'_mi_distributed.scalapack'
dMI.from_file(mi_filename, global_shape=[g_nrows,g_nrows], dtype=np.float64, block_shape=[block_size,block_size])


# ## Now we need to create a normalization matrix
# ### We start with an empty matrix but add the the diagonal as the first column
# ### Then we multiply by its transpose to get a dense matrix





#get global indices for diagonal
gi, lri, lci = dMI.local_diagonal_indices()





#create matrix to store diagonal row
dMI_diag=core.DistributedMatrix.empty_like(dMI)
dMI_row1=core.DistributedMatrix.empty_like(dMI)





dgi, dlri, dlci = dMI_diag.local_diagonal_indices()





dMI_diag.local_array[dlri,dlci]=dMI.local_array[lri,lci]
#dMI_diag.local_array[0,dlci]=dMI.local_array[lri,lci]
#my_diag[comm.rank]=dMI.local_array[lri,lci]


# ## Create a matrix with ones in the first row and zeros elsewhere



ri, ci = dMI_row1.indices()
dMI_row1.local_array[:]= ((ri==0).astype(int)).astype(float)
#print(dMI_row1.local_array)


# ## Multiply the matrices to get diagonal values on first row of distributed matrix



#dMI_norm = rt.dot(dMI_diag,dMI_row1,transA='T')
dMI_norm = rt.dot(dMI_row1,dMI_diag)
#dMI_norm = dMI_diag.dot(dMI_row1)


# ## Multiply the matrix with its transpose to get a dense matrix for normalization



import scalapy.routines as rt
#dMI_norm=dMI_diag.T*dMI_diag
dMI_norm2 = rt.dot(dMI_norm,dMI_norm,transA='T')


# ## Use scalapack to compute distributed GEMM

# ## Compute the square root of the normalization matrix



#compute sqrt of each element
dMI_norm_square=core.DistributedMatrix.empty_like(dMI)
dMI_norm_square.local_array[:] = np.sqrt(dMI_norm2.local_array[:])


# ## Now we can finally compute the norm of the MI matrix



dMI_normed=core.DistributedMatrix.empty_like(dMI)
dMI_normed.local_array[:] = dMI.local_array[:] / dMI_norm_square.local_array[:]




mi_normed_filename = data_file_path+project_name+'_mi_normed_distributed.scalapack'
dMI_normed.to_file(mi_normed_filename)


# ## Now compute eigenvalues and eigenvectors of dissimmilarity matrix




import time
start = time.time()
import scalapy.routines as rt

n= g_nrows

#convert similarity matrix to dissimilarity matrix
#df= 1-df
MDS= core.DistributedMatrix.empty_like(dMI)
MDS.local_array[:]=1.0-dMI_normed.local_array[:]

# H = I-Ones/n
I= core.DistributedMatrix.identity(n=g_nrows)
Ones= core.DistributedMatrix.empty_like(dMI)
Ones.local_array[:]=1.0/n
H = core.DistributedMatrix.empty_like(dMI)
H.local_array[:] = I.local_array[:] - Ones.local_array[:]

# B = -H.dot(MDS**2).dot(H)/2
negH= core.DistributedMatrix.empty_like(dMI)
negH.local_array[:]= -H.local_array[:]
MDS2= core.DistributedMatrix.empty_like(dMI)
MDS2.local_array[:] = MDS.local_array[:]**2
C=rt.dot(negH,MDS2)
B = rt.dot(C,H)
B.local_array[:]=B.local_array[:]/2.0
#dMI_norm=dMI_diag.T*dMI_diag
#dMI_norm = rt.dot(dMI_diag,dMI_diag,transA='T')

end = time.time()
print("Elapsed = %s" % (end - start))




import time
start = time.time()

import scalapy.routines as rt
#compute eigh(B,)
#we want to pick out the top 200 eigenvalues/vectors from the matrix

evals, dZd = rt.eigh(B,eigvals=(n - np.min([n, 200]), n - 1)) 

#copy evecs to root
evecs = dZd.to_global_array(rank=0)
#gZd = dZd.to_global_array(rank=0)

end = time.time()
print("Elapsed = %s" % (end - start))  


# ## gather the top 200 eigenvalues on a single rank

# ## Read in original dataframe to get labels to attach to results



#get index names from original dataframe
import pandas as pd
if rank==0:
    data_file_path = cwd+'/test_data/inputs/10x/PBMC/3k/pre-processed/'
    input_file_name = data_file_path + 'pbmc3k_preprocessed.h5ad'
    adf=utils.read_anndata_file(input_file_name)
    index=adf.obs.index


# ## Postprocess the eigenvalues by sorting and removing negative vals



if rank==0:
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


# ## Write reduced data to file



if rank==0:
    data_file_path = cwd+'/test_data/inputs/10x/PBMC/3k/pre-processed/'
    out_file_name = data_file_path + 'pbmc3k_preprocessed'    
    Y.to_hdf(out_file_name + "_reduced.h5", "mds")  # save reduced mi in mds




if rank==0:
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






