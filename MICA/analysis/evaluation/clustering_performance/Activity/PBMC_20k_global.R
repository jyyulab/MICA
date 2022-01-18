library(scMINER)

# Read PBMC preprocessed matrix
exp.log2 <- read.table(file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_MICA_input.txt', 
                       sep = '\t', row.names = 1, header = TRUE)
exp.log2.mat <- data.matrix(exp.log2)
exp.log2.mat.T <- t(exp.log2.mat)


# Create sparse eset object
eset <- CreateSparseEset(data=exp.log2.mat.T, meta.data = NULL, feature.data = NULL, add.meta = F)
# Read MICA clustering outputs
eset_mica <- readMICAoutput(eset = eset, load_ClusterRes = TRUE, 
                            output_file = "/Users/lding/Documents/MICA/Activity/Clustering/PBMC_20k/clustering_UMAP_euclidean_24_1.8.txt")


indx<-factor(x=c("1","2","3","4","5","6","7","8","9","10"),
             levels=c("1","2","3","4","5","6","7","8","9","10"))
eset_mica$celltype <- eset_mica$ClusterRes



act_mat <- GetActivityFromSJARACNe(SJARACNe_output_path = '/Users/lding/Documents/MICA/Activity/SJARACNe/PBMC_20k/SJARACNE_out.final/TF',
                                   SJARACNe_input_eset = eset_mica,
                                   activity.method="unweighted", # we highly recommend using 'unweighted' as activity calculation method
                                   activity.norm=TRUE, 
                                   group_name = "celltype", # which group was used to partition expression profiles
                                   save_network_file=TRUE, # whether or not save network for each group
                                   functype = 'tf',
                                   save_path="/Users/lding/Documents/MICA/Activity/SJARACNe/PBMC_20k/SJARACNE_out.final/TF") #default as false, but recommended to be TRUE

# Error out...
# Retrieve Network from  1 /consensus 
# Network saved for  /consensus 
# Calculate Activity for  /consensus ! 
#   Using geneSymbol to match targets! 
#   normalized!
#   Using geneSymbol to match targets! 
#   Activity Done!! 
#   ============================================== 
#   Error in acs.mtx[, rownames(pd)] : subscript out of bounds









# Go back to NetBID2 functions for calculating activity
library(NetBID2)
tf.network <- get.SJAracne.network(network_file = '/Users/lding/Documents/MICA/Activity/SJARACNe/PBMC_20k/SJARACNE_out.final/TF/consensus_network_ncol_.txt')
sig.network <- get.SJAracne.network(network_file = '/Users/lding/Documents/MICA/Activity/SJARACNe/PBMC_20k/SJARACNE_out.final/SIG/consensus_network_ncol_.txt')



# Merge network first
merge.network <- merge_TF_SIG.network(TF_network=tf.network, SIG_network=sig.network)



# Read PBMC preprocessed matrix
exp.log2 <- read.table(file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_MICA_input.txt', 
                       sep = '\t', row.names = 1, header = TRUE)
exp.log2.mat <- data.matrix(exp.log2)
exp.log2.mat.T <- t(exp.log2.mat)


exp.log2.mat.T.downsample <- exp.log2.mat.T[,0:1000]

ac_mat <- cal.Activity(target_list = merge.network$target_list,
                       cal_mat = exp.log2.mat.T.downsample, es.method='weightedmean',
                       std = TRUE)

ac_mat_rownames <- rownames(ac_mat)


# has problems in calculating activity, scRNA global network quanlity is not that good







# try Jingjing's T-ALL bulk network 
load('/Users/lding/Documents/MICA/Activity/SJARACNe/PBMC_20k/TALL_Jingjing/addPTCRA_TARGET_net.RData')

# Network QC
out.dir.QC.tf <- '/Users/lding/Documents/MICA/Activity/SJARACNe/PBMC_20k/TALL_Jingjing'
out.dir.QC.sig <- '/Users/lding/Documents/MICA/Activity/SJARACNe/PBMC_20k/TALL_Jingjing'

draw.network.QC(addPTCRA_TALLnet$tf.network$igraph_obj, outdir=out.dir.QC.tf, prefix='TF_net_', html_info_limit=FALSE)
draw.network.QC(addPTCRA_TALLnet$sig.network$igraph_obj, outdir=out.dir.QC.sig, prefix='SIG_net_', html_info_limit=FALSE)



# Merge network first
merge.network <- merge_TF_SIG.network(TF_network=addPTCRA_TALLnet$tf.network, SIG_network=addPTCRA_TALLnet$sig.network)


# Read PBMC preprocessed matrix
exp.log2 <- read.table(file = '/Users/lding/Documents/MICA/Datasets/HPC/SilverStd/PBMC_20k/PBMC_20k_MICA_input.txt', 
                       sep = '\t', row.names = 1, header = TRUE)
exp.log2.mat <- data.matrix(exp.log2)
exp.log2.mat.T <- t(exp.log2.mat)

# exp.log2.mat.T.downsample <- exp.log2.mat.T[,0:1000]

ac_mat <- cal.Activity(target_list = merge.network$target_list,
                       cal_mat = exp.log2.mat.T, es.method='weightedmean',
                       std = TRUE)

ac_mat.T <- t(ac_mat)

write.table(ac_mat.T, file = '/Users/lding/Documents/MICA/Activity/SJARACNe/PBMC_20k/TALL_Jingjing/PBMC_20k_MICA_input_ac_mat_TALL.txt',
            sep = '\t')
