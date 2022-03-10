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
from sklearn.decomposition import PCA
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
from scalapy import *
from scalapy import core
import scalapy.routines as rt
from scalapy import blacs

#========================================================================
#global vars
comm = MPI.COMM_WORLD
#slice_size=0
WRITE_CHECKPOINT=False

#============== function definitions ========================

#@profile
def obtain_distributed_MI_matrix(nslices, slice_size, generated_file_path, project_name, g_nrows, ncols):
    rank = comm.Get_rank()
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
                #print(SM[i][j])
                if i==j: #copy lower triangular transpose to upper triangular
                    for ii in range(SM[i][j].shape[0]):
                        for jj in range(ii+1,SM[i][j].shape[1]):
                            (SM[i][j])[ii,jj]=(SM[i][j])[jj,ii]

    comm.Barrier()
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
    comm.Barrier()
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

    #del SM
    #force garbage collection
    #gc.collect()

    block_end = time.time()
    if rank==0:
        print("Copy distributed MI matrix Elapsed = %s" % (block_end - block_start),flush=True)

    ## need to also fill in empty symmetric upper triangular portion
    # Even though this is a symmetric matrix, for further processing, we need to copy block data to rest of matrix
    comm.Barrier()

    ## Write distributed MI matrix to file
    ### So we can read this in to Scalapack later on
    if rank==0:
        print("Begin write distributed MI matrix to file Elapsed = %s" % (end - start),flush=True)
    start=time.time()
    #Write MI matrix to file
    if WRITE_CHECKPOINT:
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
    rank = comm.Get_rank()

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
    #comm.Barrier()
    #del dMI_row1
    #del dMI_diag
    #gc.collect()

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
    #comm.Barrier()
    #del dMI_norm
    #del dMI_normT
    #gc.collect()

    ## Use scalapack to compute distributed GEMM

    ## Compute the square root of the normalization matrix

    #compute sqrt of each element
    dMI_norm_root=core.DistributedMatrix.empty_like(dMI)

    dMI_norm_root.local_array[:] = np.sqrt(dMI_norm2.local_array[:])

    ## Now we can finally compute the norm of the MI matrix
    dMI_normed=core.DistributedMatrix.empty_like(dMI)

    dMI_normed.local_array[:] = dMI.local_array[:] / dMI_norm_root.local_array[:]

    #Have other ranks wait until prep_dist has completed
    #comm.Barrier()
    #del dMI_norm2
    #del dMI_norm_root
    #gc.collect()

    #ceb write normed matrix
    if WRITE_CHECKPOINT:
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
    #del I
    #del Ones
    #gc.collect()
    #Have other ranks wait until prep_dist has completed

    # B = -H.dot(MDS**2).dot(H)/2
    negH= core.DistributedMatrix.empty_like(dMI)

    negH.local_array[:]= -H.local_array[:]

    MDS2= core.DistributedMatrix.empty_like(dMI)

    MDS2.local_array[:] = MDS.local_array[:]**2

    #Have other ranks wait until prep_dist has completed
    comm.Barrier()

    start = time.time()
    if rank==0:
        print("computing dot(negH,MDS2) ",flush=True)
    C = rt.dot(negH,MDS2)
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("finished computing dot(negH,MDS2) Elapsed = %s" % (end - start),flush=True)    


    start = time.time()
    if rank==0:
        print("computing dot(C,H) ",flush=True)
    B = rt.dot(C,H)
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("finished computing dot(C,H) Elapsed = %s" % (end - start),flush=True)

    B.local_array[:]=B.local_array[:]/2.0

    #Have other ranks wait until prep_dist has completed
    comm.Barrier()
    subblock_end = time.time()
    if rank==0:
        print("Prep dissimilarity matrix Elapsed = %s" % (subblock_end - subblock_start),flush=True)

    #ceb output dissimilarity matrix here prior to eigenvalue calculation as a checkpoint
    if WRITE_CHECKPOINT:
        B.to_file(processed_dissimilarity_matrix_filename)

    return B
#================================





#@profile
def compute_eigenvalues(input_file_name,reduced_mi_filename, B, g_nrows, rank):
    if rank==0:
        print("Computing Eigenvalues of Preprocessed Dissimilarity Matrix")

    block_start = time.time()

    #compute eigh(B,) in parallel
    #we want to pick out the top 200 eigenvalues/vectors from the matrix
    start = time.time()
    n= g_nrows
    evals, dZd = rt.eigh(B,eigvals=(n - np.min([n, 200]), n - 1))
    end = time.time()
    comm.Barrier()
    if rank==0:
        print("computed eigh(B,eigvals) Elapsed = %s" % (end - start),flush=True)

    #copy evecs to root
    evecs = dZd.to_global_array(rank=0)

    block_end = time.time()
    comm.Barrier()
    if rank==0:
        print("Compute Eigenvalues: Elapsed = %s" % (block_end - block_start),flush=True)

    ## gather the top 200 eigenvalues on a single rank
    if rank==0:
        # Read in original dataframe to get labels to attach to results
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
        # Return Y to rank 0
        return Y
    else:
        # Return None to all other ranks
        return None


#=========================================== end function definintions ===============================



def main():
    ## Check to make sure MPI (mpi4py) is working
    #comm = MPI.COMM_WORLD
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
    #if rank==0:
    parser = OptionParser()
    parser.add_option("--input", action="store", type="string")
    parser.add_option("--project", action="store", type="string", default="default_proj")
    parser.add_option("--outdir", action="store", type="string", default="outdir")
    parser.add_option("--write_checkpoints", action="store_true")
    (options, args) = parser.parse_args()

    input_file_name = options.input
    project_name = options.project
    data_file_path = options.outdir

    #set global variable for writing checkpoint files
    WRITE_CHECKPOINT=options.write_checkpoints

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

    #begin whole code timer
    total_time_start = time.time()

    #==================================================================================

    ## Run Prep_dist() to split file into slices ================================
    start = time.time()
    g_nrows=0 #global number of rows (cells)
    ncols=0
    nslices=0
    slice_size=0
    #Run prep.py only on one processor to create the slice files

    if rank==0:
        g_nrows, ncols, nslices, slice_size = utils.prep_dist(input_file_name, output_file_prefix, nranks)
    end = time.time()

    if rank==0:
        print("prep_dist completed Elapsed = %s" % (end - start),flush=True)

    #Have other ranks wait until prep_dist has completed
    comm.Barrier()

    #broadcast resultant variables from root to the other ranks
    ncols = comm.bcast(ncols, root=0)
    g_nrows = comm.bcast(g_nrows, root=0)
    nslices = comm.bcast(nslices, root=0)
    slice_size = comm.bcast(slice_size, root=0)
    ##=============================================================================

    comm.Barrier()

    #set default block size for distributed cyclic matrix
    block_size=32
    #if block_size is larger than slice size, adjust
    b=nslices
    if slice_size<block_size:
        block_size=slice_size

    if rank==0:
        print("rank, global nrows, ncols, slices: ",rank, g_nrows, ncols, nslices,flush=True)

    comm.Barrier()

    #Define process grid with process rows and process cols
    PR=PC=0
    #compute the process grid dimensions from nranks, optimizing to acheive as square a matrix as possible
    PR,PC = utils.compute_process_grid_dimensions(comm.size)

    if rank==0:
        print("PR=",PR," PC=",PC," block_size=",block_size,flush=True)

    comm.Barrier()

    start = time.time()
    if rank==0: 
        print("Initializing distributed MI matrix",flush=True)

    core.initmpi([PR, PC])
    end = time.time()

    if rank==0:
        print("Completed Initializing distributed MI matrix elapsed %s" % (end - start),flush=True)

    #=========================================================================
    #if MI file already exists, then read it
    mi_filename = output_file_prefix+'_mi_distributed.scalapack'

    #==================================================================================
    dMI=None
    if os.path.isfile(mi_filename):
        ## The following code snippet reads MI matrix from a file and loads it into a distributed Scalapack matrix
        #Read MI matrix from file
        if rank==0:
            print("reading computed MI matrix from file: ",mi_filename,flush=True)
        #dMI=core.DistributedMatrix.from_file(mi_filename, global_shape=[g_nrows,g_nrows], dtype=np.float64, block_shape=[block_size,block_size])
        dMI=core.DistributedMatrix.from_file(mi_filename, global_shape=[g_nrows,g_nrows], dtype=np.float64)
    else: #file does not exist, so compute MI matrix and write
        if rank==0:
            print("Computed distributed MI matrix: ",mi_filename,flush=True)
        dMI = obtain_distributed_MI_matrix(nslices, slice_size, generated_file_path, project_name, g_nrows, ncols)

    ## End compute MI matrix and write
    comm.Barrier()

    #==================================================================================

    #if MI file already exists, then read it
    mi_normed_filename = output_file_prefix+'_mi_normed_distributed.scalapack'
    dMI_normed=None
    if os.path.isfile(mi_normed_filename):
        start=time.time()
        if rank==0:
            print("Reading distributed normed MI matrix from file",flush=True)
        dMI_normed=core.DistributedMatrix.from_file(mi_normed_filename, global_shape=[g_nrows,g_nrows], dtype=np.float64, block_shape=[block_size,block_size])
        end=time.time()
        if rank==0:
            print("Completed Reading distributed normed MI matrix from file Elapsed = %s" % (end - start),flush=True)
    else:
        if rank==0:
            print("Computing distributed normed MI matrix",flush=True)
        dMI_normed=normalize_distributed_matrix(dMI)

    # Now we need to create a normalization matrix
    # We start with an empty matrix but add the the diagonal as the first column
    # Then we multiply by its transpose to get a dense matrix

    #==================================================================================

    #Read eigenvalue matrix here if exists, else compute it
    #if MI file already exists, then read it
    processed_dissimilarity_matrix_filename = output_file_prefix+'_dissimilarity_matrix_distributed.scalapack'
    B=None
    if os.path.isfile(processed_dissimilarity_matrix_filename):
        if rank==0:
            print("Reading preprocessed dissimilarity matrix from file: "+processed_dissimilarity_matrix_filename,flush=True)
        B=core.DistributedMatrix.from_file(processed_dissimilarity_matrix_filename, global_shape=[g_nrows,g_nrows], dtype=np.float64, block_shape=[block_size,block_size])
    else: 
        if rank==0:
            print("Preprocessing dissimilarity matrix ",flush=True)
        B=process_dissimilarity_matrix(rank, processed_dissimilarity_matrix_filename, dMI, dMI_normed, g_nrows, block_size)
    #==================================================================================

    #ceb read if already existing, or compute reduced matrix
    reduced_mi_filename = output_file_prefix+'_mi_reduced.h5'
    Y=None
    if os.path.isfile(reduced_mi_filename):
        if rank==0:
            print("Reading existing eigenvalue file",flush=True)
            Y = pd.read_hdf(reduced_mi_filename)
    else:
        if rank==0:
            print("Computing eigenvalues",flush=True)
        Y = compute_eigenvalues(input_file_name,reduced_mi_filename, B, g_nrows, rank)

    #print reduced data to screen
    if rank==0:
        print("Reduced Data Matrix",flush=True)
        print(Y)
    #==================================================================================

    total_time_end = time.time()
    if rank==0:
        print("Total Time Elapsed = %s" % (total_time_end - total_time_start),flush=True)

    comm.Barrier()
    #gracefully end program
    MPI.Finalize()
    exit()



#==================================================================================
if __name__ == "__main__":
    main()



