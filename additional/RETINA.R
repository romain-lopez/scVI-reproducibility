
# ORIGINAL PIPELINE FROM NATURE PAPER
source("~/R/class.R")
load("~/bipolar/bipolar_data_Cell2016.Rdata")
print(dim(bipolar_dge))

mt.genes = grep("mt-", rownames(bipolar_dge), value = TRUE)
cells.use = colnames(bipolar_dge)[colSums(bipolar_dge[mt.genes, ])/colSums(bipolar_dge) <
                                    0.1]
bipolar_dge = bipolar_dge[, cells.use]
dim(bipolar_dge)

#Execution time ~ 10-15 minutes
dsq.bip=scDrop(count.data=bipolar_dge)
dsq.bip=initialize(dsq.bip, min.genes = 500, min.cells = 30, min.counts=60)

dim(dsq.bip@count.data)



batchname = as.character(dsq.bip@meta$sample)
batchid = rep(1,length(batchname))
batchid[batchname=="Bipolar5"] = 2
batchid[batchname=="Bipolar6"] = 2

dsq.bip=doBatchCorrection(dsq.bip, batch.cov=batchid)

# NOW look at the clustering
dsq.bip@pca.scores = pca.scores
dsq.bip@pca.load = pca.load
data.plot = dsq.bip@pca.scores
data.plot$group = dsq.bip@group
ggplot(data.plot, aes(x = PC1, y = PC2)) + geom_point(aes(colour = factor(group),
                                                          size = 1)) + scale_color_hue(l = 55) + scale_size(range = c(1, 1)) + theme_bw()

dsq.bip = doGraph_clustering(dsq.bip, pcs.use = 1:37, num.nn = 30, do.jaccard = TRUE,
                             method = "Louvain")

TPM.mat = exp(dsq.bip@data) - 1 #Recompute normalized TPM counts from log-transformed values
Count.mat = dsq.bip@count.data
markers.6 = markers.binom(dsq.bip, clust.1 = 6, effect.size = log(2), TPM.mat = TPM.mat,
                          Count.mat = Count.mat)
markers.1vs2 = markers.binom(dsq.bip, clust.1 = 1, clust.2 = 2, effect.size = log(2),
                             TPM.mat = TPM.mat, Count.mat = Count.mat)
dsq.bip = merge.clusters.DE(dsq.bip, min.de.genes = 50, pcs.use = 1:37, TPM.mat = TPM.mat,
                            Count.mat = Count.mat)
ident.louvain = dsq.bip@group


# RETAIN BIO INFORMATION FROM MERGED CLUSTERS
write.csv(ident.louvain, "~/bipolar/clustering.csv")

#SAVE FULL SPARSE EXPRESSION MATRIX
write.csv(batchid, "~/bipolar/batch_id.csv")
write.csv(rownames(dsq.bip@count.data), "~/bipolar/gene_names.csv")
write.csv(colnames(dsq.bip@count.data), "~/bipolar/barcodes.csv")
library(Matrix)
counts = Matrix(as.matrix(dsq.bip@count.data), sparse = TRUE)
writeMM(counts,file='~/bipolar/count.mtx')

dsq.bip=doBatchCorrection(dsq.bip, batch.cov=batchid)
dsq.bip=doPCA(dsq.bip,pcs.store=100)

#SAVE PCs of COMBAT OUTPUT
write.csv(dsq.bip@pca.scores, "~/bipolar/pca_combat.csv")

