#!/usr/bin/env python3
"""
This module contains helper functions, essential to the execution of MICA (Mutual Information-based
clustering algorithm)
"""

import sys
import time

# import umap as ump  # will pop up "joblib" deprecation warning message
import numpy as np
import pandas as pd
import h5py
import matplotlib  # for plotting
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import anndata

import numba
from numba import jit

from scipy.sparse import csr_matrix
from mpi4py import MPI

# from sklearn.metrics import mutual_info_score
from sklearn import cluster        # for kmeanmerge_dist_matss
from sklearn import manifold       # for tsne
from sklearn import decomposition  # for PCA, etc
from scipy.cluster.hierarchy import dendrogram  # for heatmap
from scipy.linalg import eigh
from scipy.spatial import distance  # for euclidean distance


def read_file(in_file, out_file_name):
    """ Reads text file and stores data in a temporary HDF5-format file.

    Args:
        in_file_name  (str): path to input text file
        out_file_name (str): user-defined name for output file
    """
    if in_file.endswith('.txt'):
        #read in csv text file and convert to pandas dataframe
        frame = pd.read_csv(in_file, sep="\t", index_col=0).iloc[:, 0:]
    if in_file.endswith('.h5ad'):
        adata = anndata.read_h5ad(in_file)
        #if read in as anndata form, convert to pandas dataframe
        frame = adata.to_df()
    #write dataframe to hdf
    frame.to_hdf(out_file_name + ".h5.tmp", "slice_0")
    
    
def read_anndata_file(in_file):
    """ Reads text file and stores data in a temporary HDF5-format file.

    Args:
        in_file_name  (str): path to input text file
        out_file_name (str): user-defined name for output file
    """
    if in_file.endswith('.h5ad'):
        adata = anndata.read_h5ad(in_file)
    #write dataframe to hdf
    #adata.to_hdf(out_file_name + ".h5.tmp", "slice_0")
    return adata


def slice_file(df_file,  out_file_name, slice_size="1000"):
    """ Slices the HDF5 file.

    Determines the number of slices for the data, based on slice_size.
    Calculates start and end indices for slicing, creating a new dataframe
    based on those indices and appends them to the sliced output file using
    a unique-identifier in the format 'slice_00x' as key.

    Args:
        df_file       (str): path to HDF5-format file
        out_file_name (str): path to sliced output HDF5 file
        slice_size    (str): number of items in each slice
    """

    slice_size = int(slice_size)
    df = pd.HDFStore(df_file)["slice_0"]
    b = int(np.ceil(float(df.shape[0]) / float(slice_size)))
    digit = int(np.floor(np.log10(b)) + 1)
    for i in range(b):
        slice_name = str(i).zfill(digit)
        start = i * slice_size
        end = np.min([(i + 1) * slice_size, df.shape[0]])
        #create dataframe from selected set of rows, and all columns
        slice_ = pd.DataFrame(
            data=df.iloc[start:end, :].values,
            index=df.index[start:end],
            columns=df.columns,
        )
        # write slice of original dataframe to file
        slice_.to_hdf(out_file_name + ".sliced.h5", "slice_" + slice_name)
        
        
    pd.DataFrame(data=np.array(df.shape + (b,)), index=["row", "col", "slice"]).to_hdf(
        out_file_name + ".sliced.h5", "attr"
    )
    pd.DataFrame(data=df.columns, index=df.columns).to_hdf(
        out_file_name + ".sliced.h5", "cols"
    )
    pd.DataFrame(data=df.index, index=df.index).to_hdf(
        out_file_name + ".sliced.h5", "rows"
    )


def patch_file(df_file,  out_file_name):
    """ Prepares the HDF5 file for slicing. Completes the "temporary" HDF5-format file.

    Reads input file into several data frames. Indexes attributes as row, col and slice.
    Indexes columns and rows of data. Converts all data frames into an output HDF5 file
    with separate keys for each piece of data.
    
    Args:
        df_file       (str): path to HDF5-format file
        out_file_name (str): path to complete HDF5-format file
    """

    df = pd.HDFStore(df_file)["slice_0"]
    df.to_hdf(out_file_name + ".whole.h5", "slice_0")
    pd.DataFrame(data=np.array(df.shape + (1,)), index=["row", "col", "slice"]).to_hdf(
        out_file_name + ".whole.h5", "attr"
    )
    pd.DataFrame(data=df.columns, index=df.columns).to_hdf(
        out_file_name + ".whole.h5", "cols"
    )
    pd.DataFrame(data=df.index, index=df.index).to_hdf(
        out_file_name + ".whole.h5", "rows"
    )


def calc_prep(in_file, project_name):
    """ Prepares the already sliced input file for further calculation in MICA.
    
    Enters pairs of slices (matrices) into temporary HDF5-format files. It enters them
    individually, using their unique key. It also enters the parameter data for every single 
    pair into the key "params", which consists of: [key1, key2, num_bins, num_genes,
    pair_index, project_name, num_slices]
    
    Args:
        in_file      (str): path to sliced HDF5-format file
        project_name (str): project name used to generate path for final outputs
    """

    #input file is full input that has been segmented into b blocks of rows
    
    in_ = pd.HDFStore(in_file, "r")  # slice.h5
    bins = int(np.floor(in_["attr"].loc["row"] ** (1 / 3.0)))  # number of bins
    b = in_["attr"].loc["slice", 0]  # number of sliced matrix
    m = in_["attr"].loc["col", 0]  # number of genes
    digit = int(np.floor(np.log10(b)) + 1)  # some unique identifier
    total = int((b * (b + 1)) / 2)  # total number of calc jobs execute
    digit1 = int(np.floor(np.log10(total)) + 1)
    
    for i in range(b):
        key1 = "slice_" + str(i).zfill(digit)  # location of sliced matrix 1
        mat1 = in_[key1]  # sliced matrix 1def
        for j in range(i, b):
            key2 = "slice_" + str(j).zfill(digit)
            mat2 = in_[key2]
            idx = int(i * b + j - (i * (i + 1)) / 2)
            name = "mi_" + str(idx).zfill(digit1)  # name of MI pair
            mat_tmp = (project_name + "_" + name + ".h5.tmp")  # tmp directory as MIE_out/.tmp
            mat1.to_hdf(mat_tmp, key1)  # key: slice_00x
            mat2.to_hdf(mat_tmp, key2)
            
            pd.DataFrame(data=np.array([key1, key2, bins, m, name,  project_name, b]),
                         index=["key1", "key2",
                                "num_bins", "n_genes",
                                "MI_indx", "project_name", "num_slices"]).to_hdf(mat_tmp, "params")
    in_.close()

    
def prep_dist(input_file, out_name, slice_unit):
    import logging
    """ Preprocess input file to create sliced matrices.

    Reads file into HDF5-format, adds parameters, slices data in file, generates several files with
    different combinations of slices.

    Args:
        input_file (str): path to input text-file
        out_name   (str): common rootname of generated output files
        slice_unit (int): size of each slice of cell data in input text-file
    """
    logging.basicConfig(level=logging.INFO)
    #Read in whole file stored in anndata csr format
    #adf=utils.read_anndata_file(input_file)
    adf=read_anndata_file(input_file)
    #make sure the matrix is in csr format
    from scipy.sparse import csr_matrix
    adf.X = csr_matrix(adf.X)
    
    #if adf==None :
    #    raise Exception("Input file ",input_file," not found.")
    #print("prep_dist: initial data size=",adf.shape,flush=True)
    slice_size = int(slice_unit)
    #compute number of slices needed to break dataset in to slice_size row blocks
    b = int(np.ceil(float(adf.shape[0]) / float(slice_size)))
    #determine how many digits are in b so we can pad spaces for the string output
    digit = int(np.floor(np.log10(b)) + 1)
    #print("prep_dist: slices=",b,flush=True)

    #loop over slice numbers
    for i in range(b):
        #slice name is equal to batch index
        slice_name = str(i).zfill(digit)
        #compute batch row indices
        start = i * slice_size
        end = np.min([(i + 1) * slice_size, adf.shape[0]])
        #copy slice to array of slices
        adf_sub=adf[start:end,:]
        #write to file so we don't have to keep each slice in memory
        output_file_name = out_name + ".slice_" + slice_name +".h5ad"
        #print("output_file_name: ",output_file_name)
        adf_sub.write(output_file_name)
    #return nrows, ncols, and nslices
    return adf.shape[0], adf.shape[1], b




def vpearson(X, y):
    X_mean = np.reshape(np.mean(X, axis=1), (X.shape[0], 1))
    y_mean = np.mean(y)
    r_num = np.sum((X-X_mean)*(y-y_mean), axis=1)
    r_den = np.sqrt(np.sum((X-X_mean)**2, axis=1)*np.sum((y-y_mean)**2))
    r = r_num/r_den
    return r


def calc_mi(arr1, arr2, bins, m):
    """ Calculates mutual information in between two cells, considering their gene expression levels
    
    This function is called by calc_distance_mat. It takes gene expression data from single cells,
    and compares them using standard calculation for mutual information. It builds a 2d histogram,
    which is used to calculate P(arr1, arr2)

    Args:
        arr1 (pandas series): gene expression data for a given cell in matrix_1
        arr2 (pandas series):
        bins           (int):
        m              (int):
    
    """
    fq = np.histogram2d(arr1.values, arr2.values, bins=(bins, bins))[0] / float(m)
    sm = np.sum(fq * float(m), axis=1)
    tm = np.sum(fq * float(m), axis=0)
    sm = np.asmatrix(sm / float(sm.sum()))
    tm = np.asmatrix(tm / float(tm.sum()))
    sm_tm = np.matmul(np.transpose(sm), tm)
    div = np.divide(fq, sm_tm, where=sm_tm != 0, out=np.zeros_like(fq))
    ent = np.log(div, where=div != 0, out=np.zeros_like(div))
    agg = np.multiply(fq, ent, out=np.zeros_like(fq), where=fq != 0)
    return agg.sum()







#ceb numba processed functions for distributed version of code
@numba.jit(nopython=True)
def numba_nan_fill(x):
    shape = x.shape
    x = x.ravel()
    x[np.isnan(x)] = 0.0
    x = x.reshape(shape)
    return x

#replace inf with 0
@numba.jit(nopython=True)
def numba_inf_fill(x):
    shape = x.shape
    x = x.ravel()
    x[np.isinf(x)] = 0.0
    x = x.reshape(shape)
    return x

@numba.jit(nopython=True, fastmath=True)
def compute_bin(x, min, max, num_bins):
    """ Compute bin index for a give number.
    """
    # special case to mirror NumPy behavior for last bin
    if x == max:
        return num_bins - 1 # a_max always in last bin
    bin = int(num_bins * (x - min) / (max - min))
    if bin < 0 or bin >= num_bins:
        return None
    else:
        return bin
    
    
@numba.jit(nopython=True, fastmath=True)
def compute_bin_upperbound(x, max, num_bins):
    """ Compute bin index for a give number.
        Assume that min is always zero
    """
    # special case to mirror NumPy behavior for last bin
    if x == max:
        return num_bins - 1 # a_max always in last bin
    bin = int(num_bins * x / max)
    if bin >= num_bins:
        return None
    else:
        return bin

    
@numba.jit(nopython=True, fastmath=True)
def numba_histogram2d(arr1, arr2, num_bins):
    """ Compute the bi-dimensional histogram of two data samples.
    Args:
        arr1 (array_like, shape (N,)): An array containing the x coordinates of the points to be histogrammed.
        arr2 (array_like, shape (N,)): An array containing the y coordinates of the points to be histogrammed.
        num_bins (int): int
    Return:
        hist (2D ndarray)
    """
    bin_indices1 = np.zeros((arr1.shape[0],), dtype=np.int16)
    min1 = arr1.min()
    max1 = arr1.max()
    #note that bin_indices has same size/indices as full array x and y
    for i, x in enumerate(arr1.flat):
        bin_indices1[i] = compute_bin(x, min1, max1, num_bins)

    bin_indices2 = np.zeros((arr2.shape[0],), dtype=np.int16)
    min2 = arr2.min()
    max2 = arr2.max()
    for i, y in enumerate(arr2.flat):
        bin_indices2[i] = compute_bin(y, min2, max2, num_bins)

    hist = np.zeros((num_bins, num_bins), dtype=np.int16)
    for i, b in enumerate(bin_indices1):
        hist[b, bin_indices2[i]] += 1
    return hist
    
    
#ceb create csr version of numba_histogram2d, also compute_bin with knowledge that minx will always be zero
@numba.jit(nopython=True, fastmath=True)
def numba_histogram2d_csr(arr1, cols1, arr2, cols2, ncols, num_bins):
    """ Compute the bi-dimensional histogram of two data samples.
    Args:
        arr1 (array_like, shape (N,)): An array containing the x coordinates of the points to be histogrammed.
        arr2 (array_like, shape (N,)): An array containing the y coordinates of the points to be histogrammed.
        num_bins (int): int
    Return:
        hist (2D ndarray)
    """
    #ceb short circuit test
    #return np.zeros((num_bins, num_bins), dtype=np.int16)

    #for csr arrays we have to compute zero bins ahead of time 
    bin_indices1 = np.zeros((ncols,), dtype=np.int16)
    max1 = arr1.max()
    #note that bin_indices has same size/indices as full array x and y
    for i, x in enumerate(arr1.flat):
        #assume zero min
        bin_indices1[cols1[i]] = compute_bin_upperbound(x, max1, num_bins)
        #bin_indices1[cols1[i]] = compute_bin(x, 0, max1, num_bins)
    bin_indices2 = np.zeros((ncols,), dtype=np.int16)
    max2 = arr2.max()
    for i, y in enumerate(arr2.flat):
        #assume zero min
        bin_indices2[cols2[i]] = compute_bin_upperbound(y, max2, num_bins)
        #bin_indices2[cols2[i]] = compute_bin(y, 0, max2, num_bins)
    hist = np.zeros((num_bins, num_bins), dtype=np.int16)
    for i, b in enumerate(bin_indices1):
        hist[b, bin_indices2[i]] += 1
    return hist



@numba.jit(nopython=True, fastmath=True)
def numba_calc_mi_dis(arr1, arr2, bins, m):
    """ Calculates a mutual information distance D(X, Y) = H(X, Y) - I(X, Y) using bin-based method

    It takes gene expression data from single cells, and compares them using standard calculation for
    mutual information and joint entropy. It builds a 2d histogram, which is used to calculate P(arr1, arr2).

    Args:
        arr1 (pandas series): gene expression data for cell 1
        arr2 (pandas series): gene expression data for cell 2
        marginals  (ndarray): marginal probability matrix
        index1         (int): index of cell 1
        index2         (int): index of cell 2
        bins           (int): number of bins
        m              (int): number of genes
    Returns:
        a float between 0 and 1
    """
    hist = numba_histogram2d(arr1, arr2, bins)
    sm = np.sum(hist, axis=1)
    tm = np.sum(hist, axis=0)
    sm = sm / float(sm.sum())
    tm = tm / float(tm.sum())

    sm_tm = np.zeros((bins, bins), dtype=np.float32)
    for i, s in enumerate(sm):
        for j, t in enumerate(tm):
            sm_tm[i, j] = s * t

    fq = hist / float(m)
    div = np.true_divide(fq, sm_tm)
    numba_nan_fill(div)
    ent = np.log(div)
    numba_inf_fill(ent)
    agg = np.multiply(fq, ent)
    joint_ent = -np.multiply(fq, numba_inf_fill(np.log(fq))).sum()
    return joint_ent - agg.sum()



@numba.jit(nopython=True, fastmath=True)
def numba_calc_mi_dis_csr(arr1, cols1, arr2, cols2, bins, m):
    """ Calculates a mutual information distance D(X, Y) = H(X, Y) - I(X, Y) using bin-based method

    It takes gene expression data from single cells, and compares them using standard calculation for
    mutual information and joint entropy. It builds a 2d histogram, which is used to calculate P(arr1, arr2).

    Args:
        arr1 (float) nparray of csr: gene expression data for cell 1
        cols1 (int): column indices for csr arr1
        arr2 (float) nparray of csr: gene expression data for cell 2
        cols2 (int): column indices for csr arr
        bins           (int): number of bins
        m              (int): number of genes
    Returns:
        a float between 0 and 1
    """
    hist = numba_histogram2d_csr(arr1, cols1, arr2, cols2, m, bins)
    
    sm = np.sum(hist, axis=1)
    tm = np.sum(hist, axis=0)
    sm = sm / float(sm.sum())
    tm = tm / float(tm.sum())

    sm_tm = np.zeros((bins, bins), dtype=np.float32)
    for i, s in enumerate(sm):
        for j, t in enumerate(tm):
            sm_tm[i, j] = s * t

    fq = hist / float(m)
    div = np.true_divide(fq, sm_tm)
    numba_nan_fill(div)
    ent = np.log(div)
    numba_inf_fill(ent)
    agg = np.multiply(fq, ent)
    #joint_ent = -np.multiply(fq, numba_inf_fill(np.log(fq))).sum()
    #return joint_ent - agg.sum()
    return agg.sum()



@numba.jit(nopython=True, fastmath=True)
def numba_calc_mi_csr(arr1, cols1, arr2, cols2, bins, m):
    """ Calculates a mutual information distance D(X, Y) = H(X, Y) - I(X, Y) using bin-based method

    It takes gene expression data from single cells, and compares them using standard calculation for
    mutual information and joint entropy. It builds a 2d histogram, which is used to calculate P(arr1, arr2).

    Args:
        arr1 (float) nparray of csr: gene expression data for cell 1
        cols1 (int): column indices for csr arr1
        arr2 (float) nparray of csr: gene expression data for cell 2
        cols2 (int): column indices for csr arr
        bins           (int): number of bins
        m              (int): number of genes
    Returns:
        a float between 0 and 1
    """
    hist = numba_histogram2d_csr(arr1, cols1, arr2, cols2, m, bins)
    sm = np.sum(hist, axis=1)
    tm = np.sum(hist, axis=0)
    sm = sm / float(sm.sum())
    tm = tm / float(tm.sum())
    sm_tm = np.zeros((bins, bins), dtype=np.float32)
    for i, s in enumerate(sm):
        for j, t in enumerate(tm):
            sm_tm[i, j] = s * t

    fq = hist / float(m)
    div = np.true_divide(fq, sm_tm)
    numba_nan_fill(div)
    ent = np.log(div)
    numba_inf_fill(ent)
    agg = np.multiply(fq, ent)
    #joint_ent = -np.multiply(fq, numba_inf_fill(np.log(fq))).sum()
    #return joint_ent - agg.sum()
    return agg.sum()

    
#numba compilation cannot interpret the creation of a 2d array inside of this function so we pass in and out SM_block instead of returning it
@numba.jit(nopython=True, fastmath=True)
def process_matrices(Arows,Amat_data,Amat_indptr,Amat_indices,
                     Brows,Bmat_data,Bmat_indptr,Bmat_indices,
                     num_bins,num_genes,
                     SM_block, symmetric=False):
    for i in range(Arows):
        Arowstart = Amat_indptr[i]
        Arowend   = Amat_indptr[i+1]
        Arow_cols = Amat_indices[Arowstart:Arowend]
        Arow_data = Amat_data[Arowstart:Arowend]
        Bstart=0
        Bend=Brows
        #if(symmetric):Bstart=i #upper triangluar
        if(symmetric):Bend=i+1 #lower triangular
        for j in range(Bstart,Bend): 
            Browstart = Bmat_indptr[j]
            Browend   = Bmat_indptr[j+1]
            Brow_cols = Bmat_indices[Browstart:Browend]
            Brow_data = Bmat_data[Browstart:Browend]               
            SM_block[i,j] = numba_calc_mi_dis_csr(Arow_data, Arow_cols, Brow_data, Brow_cols, num_bins, num_genes)
    return #SM_block

    

def calc_distance_metric_distributed(in_file_path, project_name, nrows, ncols, nslices, SM):
    
    """ Prepares the already sliced input file for further calculation in MICA.
    Enters pairs of slices (matrices) into temporary HDF5-format files. It enters them
    individually, using their unique key. It also enters the parameter data for every single 
    pair into the key "params", which consists of: [key1, key2, num_bins, num_genes,
    pair_index, project_name, num_slices]
    Args:
        in_file      (str): path to sliced HDF5-format file
        project_name (str): project name used to generate path for final outputs
        nrows : number of rows in global matrix
        ncols : number of vars in global matrix
    """

    #create a 2d list to hold blocks of similarity matrix
    #this should be a global var
    #SM = [[None for j in range(b)] for i in range(b)] 
    
    #input file is full input that has been segmented into b blocks of rows
    
    #nranks would ideally be equal to  b(b+1)/2
    import math
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    myrank = comm.Get_rank()
    name = MPI.Get_processor_name()
    
    digit = int(np.floor(np.log10(nslices)) + 1)

    num_bins = int(np.floor(nrows ** (1 / 3.0)))  # number of bins
    b = nslices #number of row blocks (cells)
    m = ncols  # number of cols per row (genes)

    n_block_comparisons = int((b * (b + 1)) / 2)  # total number of row block comparisons needed to compute entire global similarity matrix
    n_jobs_per_rank = math.ceil(n_block_comparisons/size)
    if (myrank == 0): print("block comparsons = %d. jobs per rank = %d\n" % (n_block_comparisons, n_jobs_per_rank))

    #build list of row block comparisons that current mpi rank will process
    myslices=[]
    for i in range(b):
        for j in range(i,b): # j in range [i,b]
            idx = int(i * b + j - (i * (i + 1)) / 2)
            targetrank = idx//n_jobs_per_rank
            if (targetrank == myrank): myslices.append((i,j))            

    #from list just generated, only do work assigned to your rank
    for index, tuple in enumerate(myslices):
        #print("tuple=",tuple)
        i=tuple[0] #block row
        j=tuple[1] #block col      

        #get 1st slice (row block) file
        slice_name = str(i).zfill(digit)
        ##ceb we only want to read this once per i,j combination
        input_file = in_file_path + project_name + ".slice_" + slice_name +".h5ad"
        #print("infile seg1: ",input_file)
        mat1 = read_anndata_file(input_file)
    
        #get 2nd slice (row block) file
        slice_name = str(j).zfill(digit)
        input_file = in_file_path + project_name + ".slice_" + slice_name +".h5ad"
        #print("infile seg2: ",input_file)
        mat2 = read_anndata_file(input_file)    

        #check to see if block comparison will result in a symmetric SM matrix
        # so that we can reduce the number of computations in half
        symmetric=False
        if i==j: symmetric=True
        
        #print("rank: ",myrank," comparison between segs:",i," x ",j," symmetric=",symmetric,flush=True)
            
        #compute distance metrics between row blocks
        Arows = mat1.n_obs
        Brows = mat2.n_obs
        num_genes = mat1.n_vars #we will assume Acols==Bcols==num_genes
        
        #need SM for each block pair    
        #creates local array of zeros and assigns to global 2d list
        #create matrix of zeros with row order indexing
        #SM[i][j] = np.ones(shape=(Arows,Brows), dtype = float, order = 'C') #for testing only
        SM[i][j] = np.zeros(shape=(Arows,Brows), dtype = float, order = 'C')
 
        #This numba function cannot create a numpy array internally so we return SM[i,j] as a variable
        process_matrices(mat1.n_obs,mat1.X.data,mat1.X.indptr,mat1.X.indices, 
                         mat2.n_obs,mat2.X.data,mat2.X.indptr,mat2.X.indices,
                         num_bins, num_genes, SM[i][j],
                         symmetric #if i==j we can eliminate half of computations
                        )
    return    
    
    
    
    
    

def calc_distance_mat(mat1, mat2, paras, method):
    """ Calculates a distance metric in between two matrices (slices)

    Calculates a distance metric using the preferred method of comparison. Iterates over each cell's
    gene expression data and populates a new matrix with rows and columns as cells from the input
    matrices. The resulting matrix is then converted to an HDF5-format file.

    Args:
        mat1  (pandas dataframe): a sliced part of the original matrix, with some fraction of the
                                  total cells as rows from original file and all gene expression
                                  attributes as columns
        mat2  (pandas dataframe): similar to mat1
        paras (pandas dataframe): a dataframe that holds an array of parameters from the whole dataset
        method             (str): the method to be used for the distance calculation (
                                        mutual information: "mi"
                                        euclidean distance: "euclidean"
                                        pearson correlation: "pearson"
                                        spearman correlation: "spearman")
    """

    bins = int(paras.loc["num_bins", 0])
    m = int(paras.loc["n_genes", 0])
    key = paras.loc["MI_indx", 0]

    project_name = paras.loc["project_name", 0]
    out_file_name = project_name + "_" + key + ".h5"
    print(out_file_name)

    if method == "mi":
        df = pd.DataFrame(data=0, index=mat1.index, columns=mat2.index) 
        start = time.time()
        for c in mat2.index:
            df.loc[mat1.index, c] = mat1.apply(
                calc_mi, axis=1, args=(mat2.loc[c, :], bins, m)
            )
        end = time.time()
    elif method == "euclidean":
        dist = distance.cdist(mat1, mat2, method)
        df = pd.DataFrame(data=dist, index=mat1.index, columns=mat2.index, dtype="float")
    elif method == "pearson":
        df = pd.DataFrame(data=0, index=mat1.index, columns=mat2.index)
        for c in mat2.index:
            df.loc[:, c] = vpearson(mat1.values, mat2.loc[c, :].values)
    elif method == "spearman":
        df = pd.DataFrame(data=0, index=mat1.index, columns=mat2.index)
        for c in mat2.index:
            df.loc[:, c] = vpearson(mat1.rank(axis=1).values, mat2.loc[c, :].values)
    else:
        sys.exit("Distance Metrics not supported!\n")

    df.to_hdf(out_file_name, str(key))  # key: MI_indx
    paras.to_hdf(out_file_name, "params")
    
    
    
    

    


def merge_dist_mats(mi_slices, in_common_name, metrics):
    """ Iterates over and merges all distance matrices in one HDF5-format file

    Args:
        mi_slices    (str[]): a list of paths of sliced distance matrix files
        in_common_name (str): project name
        metrics        (str): metric used to calculate distance
    """
    mi_0 = pd.HDFStore(mi_slices[0])
    n_slice = int(mi_0["params"].loc["num_slices", 0])  # number of chopped dfs
    n = ((n_slice + 1) * n_slice) / 2

    if len(mi_slices) < n:
        sys.exit("Error - missing sliced MI file(s) for merging.")

    k = mi_0.keys()[0]
    mi_0[k].to_hdf(in_common_name + "_mi_whole.h5", key=k)
    mi_0.close()

    digit = int(np.floor(np.log10(n)) + 1)
    rows = []

    for i in range(len(mi_slices))[1:]:
        mi_k = pd.HDFStore(mi_slices[i])
        k = mi_k.keys()[0]
        mi_k[k].to_hdf(in_common_name + "_mi_whole.h5", key=k)
        mi_k.close()

    mi_whole = pd.HDFStore(in_common_name + "_mi_whole.h5")

    for i in range(n_slice):
        row_cols = []
        for j in range(i):
            idx = int(j * n_slice + i - (j * (j + 1)) / 2)
            key = "mi_" + str(idx).zfill(digit)  # mi_idx (name)
            mat_t = mi_whole[key]
            row_cols.append(mat_t.T)  # transposed MI add to merged file
        for j in range(i, n_slice):
            idx = int(i * n_slice + j - (i * (i + 1)) / 2)
            key = "mi_" + str(idx).zfill(digit)  # mi_idx (name)
            mat = mi_whole[key]
            row_cols.append(mat)
        rows.append(pd.concat(row_cols, axis=1))

    df = pd.concat(rows, axis=0)
    df.to_hdf(in_common_name + "_dist.h5", metrics)


def norm_mi_mat(in_mat_file, out_file_name):
    """Normalizes mutual information metric in the merged matrix
    
    Args:
        in_mat_file   (str): path to merged matrix
        out_file_name (str): name of output file
    """
    hdf = pd.HDFStore(in_mat_file)
    df = hdf["mi"]
    diag = np.asmatrix(np.diag(df))
    if in_mat_file == out_file_name + "_dist.h5":
        hdf.put("norm_mi", df / np.sqrt(np.matmul(diag.T, diag)))
    else:
        df.to_hdf(out_file_name + "_dist.h5", "mi")
        (df / np.sqrt(np.matmul(diag.T, diag))).to_hdf(
            out_file_name + "_dist.h5", "norm_mi"
        )
    hdf.close()


def tsne(
         data, max_dim, out_file_name, tag, perplexity=30, plot="True"
):
    embed = manifold.TSNE(
        n_components=2,
        n_iter=5000,
        learning_rate=200,
        perplexity=perplexity,
        random_state=10,
        early_exaggeration=12.0).fit_transform(data.iloc[:, 0:max_dim])
    res = pd.DataFrame(data=embed, index=data.index, columns=["X", "Y"])
    if plot == "True":
        scatter(embed, out_file_name, tag)
    return res


def umap(
        data, max_dim, min_dist=0.25
):
    embed = ump.UMAP(
        random_state=30,
        metric="euclidean",
        n_neighbors=10,
        min_dist=min_dist).fit_transform(data.iloc[:, 0:max_dim])
    res = pd.DataFrame(data=embed, index=data.index, columns=["X", "Y"])
    return res


def mds(
        in_mat_file, max_dim, out_file_name, perplexity=30, print_plot="True", dist_method="mi",
):
    hdf = pd.HDFStore(in_mat_file)
    if dist_method == "mi":
        df = 1 - hdf["norm_mi"]
    elif dist_method == "euclidean":
        df = hdf[dist_method]
    else:
        df = 1 - hdf[dist_method]

    hdf.close()
    n = df.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    B = -H.dot(df ** 2).dot(H) / 2
    evals, evecs = eigh(B, eigvals=(n - np.min([n, 200]), n - 1))
    idx = np.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:, idx]
    evals_pos = evals > 0
    L = np.diag(np.sqrt(evals[evals_pos]))
    V = evecs[:, evals_pos]
    Y = pd.DataFrame(
        data=V.dot(L),
        index=df.index,
        columns=["mds_" + str(x) for x in np.arange(1, L.shape[0] + 1)],
    )

    Y.to_hdf(out_file_name + "_reduced.h5", "mds")  # save reduced mi in mds

    if print_plot == "True":
        vis = tsne(
            Y,
            max_dim,
            out_file_name,
            "mds",
            perplexity,
            print_plot,
        )
        vis.to_hdf(out_file_name + "_reduced", "mds_tsne")  # save preview in key "mds_tsne"


def lpl(
         in_mat_file, max_dim, out_file_name, perplexity=30, plot="True", dist_method="mi",
):
    hdf = pd.HDFStore(in_mat_file)

    if dist_method == "mi":
        df = 1 - hdf["norm_mi"]
    elif dist_method == "euclidean":
        df = hdf[dist_method]
    else:
        df = 1 - hdf[dist_method]

    hdf.close()
    n = np.min(df.shape[0], 200)
    laplacian = manifold.SpectralEmbedding(
        n_components=n, eigen_solver="lobpcg", random_state=10
    )
    Y = pd.DataFrame(
        data=laplacian.fit_transform(df),
        index=df.index,
        columns=["lpl_" + str(x) for x in np.arange(1, n)],
    )
    Y.to_hdf(out_file_name + "_reduced.h5", "lpl")
    if plot == "True":
        tsne(
            Y,
            max_dim,
            out_file_name,
            "lpl",
            perplexity,
            plot,
        )


def pca(
         in_mat_file, max_dim, out_file_name, perplexity=30, plot="True", dist_method="mi",
):
    hdf = pd.HDFStore(in_mat_file)

    if dist_method == "mi":
        df = 1 - hdf["norm_mi"]
    elif dist_method == "euclidean":
        df = hdf[dist_method]
    else:
        df = 1 - hdf[dist_method]

    hdf.close()
    n = np.min(df.shape[0], 200)
    pca_ = decomposition.PCA(n_components=n, random_state=10)
    Y = pd.DataFrame(
        data=np.transpose(pca_.fit(df).components_),
        index=df.index,
        columns=["pca_" + str(x) for x in np.arange(1, n + 1)],
    )
    Y.to_hdf(out_file_name + "_reduced.h5", "pca")
    if plot == "True":
        tsne(
            Y,
            max_dim,
            out_file_name,
            "pca",
            perplexity,
            plot,
        )


def kmeans(in_mat, n_clusters, project_name, dim, bootstrap_id):
    out_file_name = project_name + "_kmeans_k" + str(n_clusters) + "_d" + str(dim) + ".h5.tmp." + str(bootstrap_id)
    km = cluster.KMeans(n_clusters=n_clusters, max_iter=1000, n_init=1000)

    km_res = pd.DataFrame(
        data=np.transpose(km.fit_predict(in_mat.iloc[:, 0:dim])),
        index=in_mat.index,
        columns=["label"],
    )
    # km_res.to_hdf(out_file_name, "kmeans")
    print("Executing kmeans:" + out_file_name)
    return km_res


def aggregate(km_results, n_clusters, common_name):
    n_iter = len(km_results)
    if n_iter == 0:
        return None
    mem = None
    for i in range(n_iter):
        df = km_results[i]
        dff = pd.DataFrame(data=df.values@df.T.values, index=df.index, columns=df.index)
        dff_div = pd.DataFrame(
            data=np.array((np.diag(dff),) * dff.shape[0]).T,
            index=dff.index,
            columns=dff.columns,
        )
        mem_mat = pd.DataFrame(
            data=dff / dff_div == 1,
            index=dff.index,
            columns=dff.columns,
            dtype=np.float32,
        )
        mem = mem_mat if i == 0 else mem + mem_mat.loc[mem.index, mem.columns]
        mem = mem / n_iter if i == n_iter - 1 else mem

    clust = cluster.AgglomerativeClustering(
        linkage="ward", n_clusters=n_clusters, affinity="euclidean"
    )
    clust.fit(mem)

    cclust = pd.DataFrame(data=clust.labels_, index=mem.index, columns=["label"])
    index = cclust.groupby(["label"]).size().sort_values(ascending=False).index
    label_map = {index[i]: i + 1000 for i in range(len(index))}
    cclust = cclust.replace(to_replace={"label": label_map})
    index = cclust.groupby(["label"]).size().sort_values(ascending=False).index
    map_back = {index[i]: i + 1 for i in range(len(index))}
    cclust = cclust.replace(to_replace={"label": map_back})

    out_file = common_name + "_k" + str(n_clusters)
    out_file_hdf = out_file + "_cclust.h5"
    cclust.to_hdf(out_file_hdf, "cclust")
    mem.to_hdf(out_file_hdf, "membership")
    return cclust, out_file


def visualization(
        agg_mtx,   # 1
        reduced_mi_file,  # 2
        transformation,   # 3
        out_file_name,
        max_dim=0,
        visualize="umap",  # 5
        min_dist=0.25,  # 6
        perplexity=30,  # 7
):

    cclust = agg_mtx
    hdf = pd.HDFStore(reduced_mi_file)
    transformation = "pca" if transformation == "lpca" else transformation
    hdf_trans = hdf[transformation.lower()]
    hdf.close()

    if visualize == "umap":
        embed_2d = umap(hdf_trans, max_dim, min_dist)

    elif visualize == "tsne":
        perplexity = np.min([perplexity, np.max(cclust.groupby(["label"]).size())])
        embed_2d = tsne(hdf_trans, max_dim, "", None, perplexity, "False")

    cclust = pd.concat([cclust, embed_2d.loc[cclust.index, :]], axis=1)
    res = cclust.loc[:, ["X", "Y", "label"]]
    # save 2D embedding to txt file
    out_file_name = out_file_name + "_" + visualize

    scatter2(res, out_file_name + '.png')
    res.to_csv(out_file_name + "_ClusterMem.txt", sep="\t")


def heatmap(df, linkage_, matrix, out_dir, out_file_name, tag, dendro="True"):
    fig = plt.figure(figsize=(16, 16), dpi=300)
    grid = gs.GridSpec(1, 1, wspace=0.2, hspace=0.2)
    subgrid = gs.GridSpecFromSubplotSpec(
        2,
        2,
        subplot_spec=grid[0],
        wspace=0.01,
        hspace=0.01,
        height_ratios=[1, 10],
        width_ratios=[10, 1],
    )
    if dendro == "True":
        ax = plt.Subplot(fig, subgrid[0])
        dendrogram(
            linkage_,
            ax=ax,
            orientation="top",
            labels=df.index,
            leaf_font_size=2,
            color_threshold=0,
        )
        ax.axis("off")
        fig.add_subplot(ax)
    ax1 = plt.Subplot(fig, subgrid[2])
    cax = ax1.matshow(matrix, cmap="YlOrRd", interpolation=None, aspect="auto")
    ax1.set_xticks([])
    ax1.set_yticks([])
    fig.add_subplot(ax1)
    ax3 = plt.Subplot(fig, subgrid[3])
    ax3.axis("off")
    cbar = fig.colorbar(cax, ax=ax3, shrink=1.5, aspect=20, fraction=0.5, pad=0)
    plt.savefig(out_dir + out_file_name + "_" + tag + ".pdf", bbox_inches="tight")


def heatmap2(data, labels, out_dir, out_file_name, tag):
    if data is None:
        return
    fig = plt.figure( figsize=(16, 16), dpi=300)
    grid = gs.GridSpec(1, 1, wspace=0.2, hspace=0.2)
    subgrid = gs.GridSpecFromSubplotSpec(
        2,
        2,
        subplot_spec=grid[0],
        wspace=0.001,
        hspace=0.1,
        height_ratios=[1, 20],
        width_ratios=[20, 1],
    )
    sorted_labels = labels.sort_values()
    bar = np.vstack((np.array(sorted_labels) + 1,) * 2) / np.max(np.unique(labels) + 1)
    ax1 = plt.Subplot(fig, subgrid[0])
    ax1.matshow(bar, cmap="jet", interpolation=None, aspect="auto", vmin=0, vmax=1)
    ax1.set_xticks([])
    ax1.set_yticks([])
    fig.add_subplot(ax1)
    ax2 = plt.Subplot(fig, subgrid[2])
    cax = ax2.matshow(
        data.loc[sorted_labels.index, sorted_labels.index],
        cmap="YlOrRd",
        interpolation=None,
        aspect="auto",
    )
    ax2.set_xticks([])
    ax2.set_yticks([])
    fig.add_subplot(ax2)
    ax3 = plt.Subplot(fig, subgrid[3])
    ax3.axis("off")
    cbar = fig.colorbar(cax, ax=ax3, shrink=1.5, aspect=20, fraction=0.5, pad=0)
    plt.savefig(out_dir + out_file_name + "_" + tag + ".pdf", bbox_inches="tight")


def scatter(
        df,
        out_file_name,
        tag,
        facecolor="none",
        edgecolor="r",
        marker="o",
        marker_size=20,
):
    fig = plt.figure(figsize=(16, 16), dpi=300)
    plt.scatter(
        df[:, 0],
        df[:, 1],
        facecolor=facecolor,
        edgecolor=edgecolor,
        marker=marker,
        s=marker_size,
    )
    plt.ylabel("MICA-2")
    plt.xlabel("MICA-1")
    plt.xticks([])
    plt.yticks([])
    plt.savefig(out_file_name + "_" + tag + ".pdf", bbox_inches="tight")


def scatter2(data, out_file_name, marker_size=20, marker="o"):
    if data is None:
        return
    fig = plt.figure(figsize=(10, 10), dpi=300)
    lab = np.unique(data.loc[:, "label"])
    colors = plt.cm.jet(np.linspace(0, 1, len(lab)))
    for z in lab:
        df = data.loc[data.loc[:, "label"] == z, ["X", "Y"]]
        plt.scatter(
            df.loc[:, "X"],
            df.loc[:, "Y"],
            facecolor=colors[z-1],
            s=marker_size,
            marker=marker,
            vmin=0,
            vmax=len(lab),
            label=str(z) + "(" + str(df.shape[0]) + ")",
            alpha=0.7,
        )
        center = np.mean(df, axis=0)
        plt.scatter(
            center.loc["X"],
            center.loc["Y"],
            marker="o",
            c="white",
            alpha=0.7,
            s=100,
            edgecolor="k",
        )
        plt.scatter(
            center.loc["X"],
            center.loc["Y"],
            marker="$%d$" % z,
            c="black",
            alpha=0.7,
            s=80,
            edgecolor="k",
        )
    plt.ylabel("MICA-2")
    plt.xlabel("MICA-1")
    plt.xticks([])
    plt.yticks([])
    plt.legend(
        loc="center left",
        bbox_to_anchor=(1, 0.5),
        title="Clusters(" + str(data.shape[0]) + ")",
    )
    plt.savefig(
        out_file_name,
        bbox_inches="tight",
    )
