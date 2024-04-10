library(NetBID2)

tf.network.file <- '/Users/lding/Documents/MICA/Activity/SJARACNe/Goolam/SJARACNE_out.final/TF/consensus_network_ncol_.txt'
sig.network.file <- '/Users/lding/Documents/MICA/Activity/SJARACNe/Goolam/SJARACNE_out.final/SIG/consensus_network_ncol_.txt'

out.dir.QC.tf <- '/Users/lding/Documents/MICA/Activity/SJARACNe/Goolam/SJARACNE_out.final/TF'
out.dir.QC.sig <- '/Users/lding/Documents/MICA/Activity/SJARACNe/Goolam/SJARACNE_out.final/SIG'

tf.network <- get.SJAracne.network(network_file = tf.network.file)
sig.network <- get.SJAracne.network(network_file = sig.network.file)

draw.network.QC(tf.network$igraph_obj, outdir=out.dir.QC.tf, prefix='TF_net_', html_info_limit=FALSE)
draw.network.QC(sig.network$igraph_obj, outdir=out.dir.QC.sig, prefix='SIG_net_', html_info_limit=FALSE)






tf.network.file <- '/Users/lding/Documents/MICA/Activity/SJARACNe/PBMC_20k/SJARACNE_out.final/TF/consensus_network_ncol_.txt'
sig.network.file <- '/Users/lding/Documents/MICA/Activity/SJARACNe/PBMC_20k/SJARACNE_out.final/SIG/consensus_network_ncol_.txt'

out.dir.QC.tf <- '/Users/lding/Documents/MICA/Activity/SJARACNe/PBMC_20k/SJARACNE_out.final/TF'
out.dir.QC.sig <- '/Users/lding/Documents/MICA/Activity/SJARACNe/PBMC_20k/SJARACNE_out.final/SIG'

tf.network <- get.SJAracne.network(network_file = tf.network.file)
sig.network <- get.SJAracne.network(network_file = sig.network.file)

draw.network.QC(tf.network$igraph_obj, outdir=out.dir.QC.tf, prefix='TF_net_', html_info_limit=FALSE)
draw.network.QC(sig.network$igraph_obj, outdir=out.dir.QC.sig, prefix='SIG_net_', html_info_limit=FALSE)

