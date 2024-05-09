#!/usr/bin/env python3

import numpy as np
import sklearn 


def mds_trans(df):
    n = df.shape[0]
    H = np.eye(n) - np.ones((n, n)) / n
    B = -H.dot(df ** 2).dot(H) / 2
    evals, evecs = np.linalg.eigh(B)
    idx = np.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:, idx]
    evals_pos = evals > 0
    L = np.diag(np.sqrt(evals[evals_pos]))
    V = evecs[:, evals_pos]
    mds_res = np.dot(V, L)
    
    return mds_res