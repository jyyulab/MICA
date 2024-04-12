#!/usr/bin/env python3
import pandas as pd
from MICA.lib import utils
from MICA.lib import preprocessing
from MICA.lib import distance
import numpy


#%% Read Buettner input matrix
h5_file = '/Users/lding/Documents/MICA/Datasets/with_true_labels/GoldernStd/Buettner/Buettner_MICA_input.h5ad'
adata = preprocessing.read_preprocessed_mat(h5_file)
frame = adata.to_df()
ndarray = frame.to_numpy()


#%%
def calc_distance_mat(mat1, mat2):
    """ Calculates a distance metric in between two matrices (slices)

    Calculates a distance metric using the preferred method of comparison. Iterates over each cell's
    gene expression data and populates a new matrix with rows and columns as cells from the input
    matrices. The resulting matrix is then converted to an HDF5-format file.

    Args:
        mat1  (pandas dataframe): a sliced part of the original matrix, with some fraction of the
                                  total cells as rows from original file and all gene expression
                                  attributes as columns
        mat2  (pandas dataframe): similar to mat1
    """
    # num_bins = int((mat1.shape[0]) ** (1 / 3.0))
    num_bins = 20
    num_genes = mat1.shape[1]

    df = pd.DataFrame(data=0, index=mat1.index, columns=mat2.index)
    for c in mat2.index:
        df.loc[mat1.index, c] = mat1.apply(
            utils.calc_mi, axis=1, args=(mat2.loc[c, :], num_bins, num_genes)
        )

    return df


#%% Check if MIE output matrix the same as the MI matrix created by Tracy or Ali
dist_mat = calc_distance_mat(frame, frame)


#%%
def calculate_MI_for_pairs(bins, matrix1, matrix2, i, j, n):
    source = matrix1[i, :]
    target = matrix2[j, :]
    frequency_table = numpy.histogram2d(source, target, bins=(bins, bins))[0]/float(n)
    source_marginal = numpy.sum(frequency_table*float(n), axis=1)
    target_marginal = numpy.sum(frequency_table*float(n), axis=0)
    source_marginal = source_marginal/float(source_marginal.sum())
    target_marginal = target_marginal/float(target_marginal.sum())
    mi = 0.0
    for l in range(0, bins):
        for m in range(0, bins):
            if source_marginal[l] != 0 and target_marginal[m] != 0 and frequency_table[l][m] != 0:
                mi += frequency_table[l][m] * numpy.log(frequency_table[l][m] / (source_marginal[l] * target_marginal[m]))
    return mi

#%%
num_bins = int((frame.shape[0]) ** (1 / 3.0))
num_genes = frame.shape[1]

#%%
calculate_MI_for_pairs(20, ndarray, ndarray, 0, 1, num_genes)
