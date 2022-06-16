#!/usr/bin/env Rscript

act_rdata_path <- '/Users/lding/Documents/scMINER/PBMC20k_CD4/ac_mat_CD4T_wightedmean.RData'
act_rdata <- load(act_rdata_path)
write.table(t(ac_mat), file = '/Users/lding/Documents/scMINER/PBMC20k_CD4/ac_mat_CD4T_wightedmean.txt', sep = '\t')

exp_rdata_path <- '/Users/lding/Documents/scMINER/PBMC20k_CD4/eset.PBMC.CD4T.RData'
exp_rdata <- load(exp_rdata_path)
write.table(t(exprs(eset.PBMC1)), file = '/Users/lding/Documents/scMINER/PBMC20k_CD4/exp_mat_CD4T.txt', sep = '\t')







act_rdata_path_tf <- '/Users/lding/Documents/scMINER/PBMC20k_CD4/ac_mat_CD4T_weightedmean_cell_type_specific_net/master_tf.RData'
act_rdata_tf <- load(act_rdata_path_tf)
acs_master_tf <- acs_master

act_rdata_path_sig <- '/Users/lding/Documents/scMINER/PBMC20k_CD4/ac_mat_CD4T_weightedmean_cell_type_specific_net/master_sig.RData'
act_rdata_sig <- load(act_rdata_path_sig)
acs_master_sig <- acs_master

acs_master <- rbind(acs_master_tf, acs_master_sig)
acs_master_tmp <- acs_master[,-1]
acs_master_mat <- as.matrix(acs_master_tmp)
rownames(acs_master_mat) <- acs_master$ID


# Replace NAs with row means
k <- which(is.na(acs_master_mat), arr.ind=TRUE)
# acs_master_mat[k] <- rowMeans(acs_master_mat, na.rm=TRUE)[k[,1]]
library(Rfast)
acs_master_mat[k] <- colMins(acs_master_mat, value = TRUE)[k[,1]]


write.table(t(acs_master_mat), file = '/Users/lding/Documents/scMINER/PBMC20k_CD4/ac_mat_CD4T_wightedmean_cts.txt', sep = '\t', col.names=NA)
