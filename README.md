# scVI-reproducibility
Reproducing the experiments of the scVI paper
+ Current development version in PyTorch is available at https://github.com/YosefLab/scVI 
+ TensorFlow Demonstration code has been moved from YosefLab/scVI to the demo_code folder in this repo
+ NIPS MLCB submission: https://arxiv.org/abs/1709.02082
+ biorXiv preprint: https://www.biorxiv.org/content/early/2018/03/30/292037

# Contact
romain [underscore] lopez [at] berkeley [dot] edu

## Datasets
+ CORTEX

Original file expression_mRNA_17-Aug-2014.txt can be dowloaded here https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt

+ PBMC

An additional R script to extract the quantity of interest from SCONE objects is needed (additional/PBMC.R)
scVI_scone.rda, scone_all_wposcon_extendim_biostrat2_q.rda, molecule_qc_8k.txt, molecule_qc_4k.txt and gene_info.csv (containing the merged p-values from the microarray experiments) are in additional/data.zip

+ BRAIN LARGE

An additional Jupyter Notebook to format the data for scVI is needed (additional/BRAIN-LARGE.ipynb)
1M_neurons_filtered_gene_bc_matrices_h5.h5 can be downloaded here http://cf.10xgenomics.com/samples/cell-exp/1.3.0/1M_neurons/1M_neurons_filtered_gene_bc_matrices_h5.h5

+ RETINA

An additional R script, based on the original analysis is needed (additional/RETINA.R)
We modified the file class.R so that the gene filtering also apply to the raw matrix, use ours in this repo (additional/class.R)
The file bipolar_data_Cell2016.Rdata can be downloaded here https://github.com/broadinstitute/BipolarCell2016
The file ClustAssignFile.txt was obtained from the author, uploaded in additional/data.zip

+ HEMATO

Original files bBM.raw_umifm_counts.csv can be downloaded here https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2388072 bBM.spring_and_pba.csv and bBM.filtered_gene_list.paper.txt were sent by authors, uploaded in additional/data.zip

+ CBMC

Original files GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv, GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv, GSE100866_CBMC_8K_13AB_10X-ADT_clr-transformed.csv can be downloaded here https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100866

+ BRAIN SMALL

Original files genes.tsv, barcodes.tsv, matrix.mtx and clusters.csv can be downloaded here https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/neuron_9k
molecule_qc.txt was extracted from the original file http://cf.10xgenomics.com/samples/cell-exp/2.1.0/neuron_9k/neuron_9k_molecule_info.h5 using the additional python script additional/molecule_info.py (credit David Detomaso)

