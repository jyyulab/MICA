#!/usr/bin/env python3

import numpy as np
import sklearn 


def entropy_calc(x):
    res = np.sum(np.where(x>0, -x*np.log(x), 0))
    return res


def mi_calc(data, bins):
    
    samples, features = data.shape
    res = np.zeros((samples, samples))
    entropy_single = []
    
    for i in data:
        hist_temp, _ = np.histogram(i, bins=bins)
        entropy_temp = entropy_calc(hist_temp/features)
        entropy_single.append(entropy_temp)
        
    for i in range(samples):
        for j in range(samples):
            joint_hist, _, _ = np.histogram2d(data[i], data[j], bins=bins)
            joint_prob = joint_hist/features
            joint_entropy = entropy_calc(joint_prob.flatten())
            midis = (entropy_single[i] + entropy_single[j]) - joint_entropy
            
            res[i,j] = midis
    
    return res


def mi_norm(mi_ori):
    diag = np.diag(mi_ori)
    return 1 - mi_ori/np.sqrt(np.dot(diag.reshape(-1, 1), diag.reshape(1, -1)))