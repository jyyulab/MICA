library(NetBID2)

transfer_tab <- get_IDtransfer2symbol2type(from_type='ensembl_gene_id',
                                           dataset = 'mmusculus_gene_ensembl',
                                           use_level='gene', ignore_version = T)


goolam <- read.csv('/Users/lding/Documents/MICA/Activity/SJARACNe/Goolam/Goolam_SJARACNe_input.txt', sep = '\t')


# Get gene IDs for goolam emsenble IDs
transfer_tab[Filter(function(x) !is.na(x), match(transfer_tab$ensembl_gene_id, goolam$isoformId)),]$external_gene_name
