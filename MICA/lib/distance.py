#!/usr/bin/env python3

import numpy as np
import numba
import fast_histogram
import pandas as pd


@numba.jit(nopython=True, fastmath=True)
def compute_bin(x, min, max, num_bins):
    """ Compute bin index for a give number.
    """
    # special case to mirror NumPy behavior for last bin
    if x == max:
        return num_bins - 1     # a_max always in last bin
    bin = int(num_bins * (x - min) / (max - min))
    if bin < 0 or bin >= num_bins:
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


@numba.jit(nopython=True)
def numba_nan_fill(x):
    """ Fill nan with 0.0 """
    shape = x.shape
    x = x.ravel()
    x[np.isnan(x)] = 0.0
    x = x.reshape(shape)
    return x


@numba.jit(nopython=True)
def numba_inf_fill(x):
    """ Fill inf or -inf with 0.0 """
    shape = x.shape
    x = x.ravel()
    x[np.isinf(x)] = 0.0
    x = x.reshape(shape)
    return x


@numba.jit(nopython=True, fastmath=True)
def numba_calc_mi_dis(arr1, arr2, bins, m):
    """ Calculates a normalized mutual information distance D(X, Y) = 1 - I(X, Y)/H(X, Y) using bin-based method

    It takes gene expression data from single cells, and compares them using standard calculation for
    mutual information and joint entropy. It builds a 2d histogram, which is used to calculate P(arr1, arr2).

    Args:
        arr1 (pandas series): gene expression data for cell 1
        arr2 (pandas series): gene expression data for cell 2
        bins           (int): number of bins
        m              (int): number of genes
    Returns:
        D(arr1, arr2) see above, a float between 0 and 1
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
    return (joint_ent - agg.sum()) / joint_ent    # normalized: D(X, Y) = 1 - I(X, Y)/H(X, Y)
    # return joint_ent - agg.sum()                # unnormalized: D(X, Y) = H(X, Y) - I(X, Y)


def calc_norm_mi(arr1, arr2, bins, m):
    """ Calculates a normalized mutual information distance D(X, Y) = 1 - I(X, Y)/H(X, Y) using bin-based method

    It takes gene expression data from single cells, and compares them using standard calculation for
    mutual information and joint entropy. It builds a 2d histogram, which is used to calculate P(arr1, arr2).

    Args:
        arr1 (pandas series): gene expression data for cell 1
        arr2 (pandas series): gene expression data for cell 2
        bins           (int): number of bins
        m              (int): number of genes
    Returns:
        D(arr1, arr2) see above, a float between 0 and 1
    """
    fq = fast_histogram.histogram2d(arr1, arr2, range=[[arr1.min(), arr1.max()+1e-9], [arr2.min(), arr2.max()+1e-9]],
                                    bins=(bins, bins)) / float(m)
    sm = np.sum(fq * float(m), axis=1)
    tm = np.sum(fq * float(m), axis=0)
    sm = np.asmatrix(sm / float(sm.sum()))
    tm = np.asmatrix(tm / float(tm.sum()))
    sm_tm = np.matmul(np.transpose(sm), tm)
    div = np.divide(fq, sm_tm, where=sm_tm != 0, out=np.zeros_like(fq))
    ent = np.log(div, where=div != 0, out=np.zeros_like(div))
    agg = np.multiply(fq, ent, out=np.zeros_like(fq), where=fq != 0)
    joint_ent = -np.multiply(fq, np.log(fq, where=fq != 0, out=np.zeros_like(fq)),
                             out=np.zeros_like(fq), where=fq != 0).sum()
    return (joint_ent - agg.sum()) / joint_ent


def calc_dis_mat_numba(mat1, mat2, bins, m):
    """ Wrapper of numba_calc_mi_dis for calculating mutual information for two matrices
    Args:
        mat1 (pandas dataframe): exp matrix of a slice of cells, with cells as rows from original file
                                  and all gene expression attributes as columns
        mat2 (pandas dataframe): exp matrix of another slice of cells
        bins              (int): number of bins
        m                 (int): number of genes
    Returns:
        df (pandas dataframe with dimension mat1.index * mat2.index)
    """
    df = pd.DataFrame(data=0, index=mat1.index, columns=mat2.index)
    for c1 in mat1.index:
        for c2 in mat2.index:
            # Can speed up calculation by only calculating half of the matrix
            df.loc[c1, c2] = numba_calc_mi_dis(mat1.loc[c1,:].to_numpy(), mat2.loc[c2, :].to_numpy(), bins, m)
    return df


def calc_dis_mat_paras(mat1, mat2, paras):
    """ Wrapper of calc_dis_mat_numba """
    bins = int(paras.loc["num_bins", 0])
    m = int(paras.loc["n_genes", 0])
    key = paras.loc["MI_indx", 0]

    project_name = paras.loc["project_name", 0]
    out_file_name = project_name + "_" + key + ".h5"
    print(out_file_name)

    df = calc_dis_mat_numba(mat1, mat2, bins, m)
    df.to_hdf(out_file_name, str(key))          # key: MI_indx
    paras.to_hdf(out_file_name, "params")


def calc_dis_mat_np(mat1, mat2, bins, m, index):
    """ Wrapper of calc_mi for calculating mutual information for two matrices(numpy.ndarray)
    Args:
        mat1 (numpy.ndarray): exp matrix of a slice of cells, with cells as rows from original file
                              and all gene expression attributes as columns
        mat2 (numpy.ndarray): exp matrix of another slice of cells
        bins           (int): number of bins
        m              (int): number of genes
        index       (series): cell index
    Returns:
        df (pandas dataframe with dimension mat1.index * mat2.index)
    """
    df = pd.DataFrame(data=0, index=index, columns=index)
    for i, c in enumerate(index):
        df.loc[index, c] = np.apply_along_axis(calc_norm_mi, 1, mat1, mat2[i, :], bins, m)
    return df
