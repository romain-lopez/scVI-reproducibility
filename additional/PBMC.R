
library(Matrix)
library(scone)
library(plyr) 


# get cleaned data by merging genes from original scone pipeline and microarray
load("/data/yosef2/pbmc8k/scVI_scone.rda")

barcodes = scone_obj@colData@rownames
list_qc = scone_obj@colData@listData[names(scone_obj@colData@listData)[1:9]]
qc.df = do.call("cbind", lapply(list_qc, as.data.frame))
colnames(qc.df) = names(scone_obj@colData@listData)[1:9]
batch = get_batch(scone_obj)
gene_names = scone_obj@NAMES
design = get_design(scone_obj, method = "none,fq,qc_k=8,no_bio,no_batch" )

write.csv(barcodes, "~/data_PBMCs/barcodes.csv")
write.csv(batch, "~/data_PBMCs/batch.csv")
write.csv(qc.df, "~/data_PBMCs/full_qc.csv")
write.csv(gene_names, "~/data_PBMCs/gene_names.csv")
write.csv(design, "~/data_PBMCs/design.csv")

# load cells information from SEURAT, included in the original scone object
load("/data/yosef2/pbmc8k/scone_all_wposcon_extendim_biostrat2_q.rda")
bio = get_bio(scone_obj)
write.csv(bio, "~/data_PBMCs/bio.csv")


