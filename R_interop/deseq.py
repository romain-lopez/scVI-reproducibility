import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import warnings
import numpy as np
from rpy2.rinterface import RRuntimeWarning
import scipy

class DESeq2(object):
    
    def __init__(self, A, B, data, labels, cluster):
        """
        A: number of cells in the first cluster
        B: number of cells in the second cluster
        data: dataset to look at
        labels: clusters
        cluster: list that tells which cluster to test ex. (0, 4)
        """
        self.A = A
        self.B = B
        self.data = data
        self.labels = labels
        self.cluster = cluster
        #loading libraries
        warnings.filterwarnings("ignore", category=RRuntimeWarning)
        rpy2.robjects.numpy2ri.activate()
        ro.r["library"]("RcppCNPy")
        ro.r["library"]("DESeq2")
        ro.r["library"]("BiocParallel")
        ro.r("BiocParallel::register(BiocParallel::MulticoreParam())")
        self.X_train = np.load(self.data)
        self.c_train = np.loadtxt(self.labels)
                
        # loading data
        ro.r(str("""fmat <- npyLoad("*""")[:-1] + self.data + str("""*", "integer")""")[1:])
        #ro.r(str("""fmat <- npyLoad("*""")[:-1] + self.data + str("""*", "numeric")""")[1:])
        ro.r(str("""cmat <- read.table("*""")[:-1] + self.labels + str("""*")""")[1:])
        ro.r("cmat$V2 <- factor(cmat$V1)")
        

    def fit(self, return_fc=False):
        
        # computing data mask
        set_a = np.where(self.c_train == self.cluster[0])[0]
        subset_a = np.random.choice(set_a, self.A)
        set_b = np.where(self.c_train == self.cluster[1])[0]
        subset_b = np.random.choice(set_b, self.B)

        stochastic_set = np.hstack((subset_a, subset_b))
        f = np.array([a in stochastic_set for a in np.arange(self.X_train.shape[0])])

        nr, nc = f[:, np.newaxis].shape
        f_r = ro.r.matrix(f[:, np.newaxis], nrow=nr, ncol=nc)
        ro.r.assign("f_", f_r)
        ro.r("f <- as.integer(rownames(cmat[f_,]))")
        
        ro.r("local_fmat <- fmat[f, ]")
        ro.r("local_cmat <- cmat[f, ]")
        ro.r("local_cmat$V3 <- factor(local_cmat$V1)")
        ro.r("dds <- DESeqDataSetFromMatrix(countData = t(local_fmat), colData = local_cmat, design = ~ V3)")
        ro.r("dds <- DESeq(dds)")
        ro.r("res <- results(dds)")
        if return_fc:
            return ro.r("res@listData$pvalue"), ro.r("res@listData$log2FoldChange")
        else:
            return ro.r("res@listData$pvalue")  
    
    
class Weighted_edgeR(object):
    
    def __init__(self, A, B, data, labels, cluster, weights):
        """
        A: number of cells in the first cluster
        B: number of cells in the second cluster
        data: dataset to look at
        labels: clusters
        cluster: list that tells which cluster to test ex. (0, 4)
        """
        self.A = A
        self.B = B
        self.data = data
        self.labels = labels
        self.cluster = cluster
        self.weights = weights
        #loading libraries
        warnings.filterwarnings("ignore", category=RRuntimeWarning)
        rpy2.robjects.numpy2ri.activate()
        ro.r["library"]("RcppCNPy")
        ro.r["library"]("DESeq2")
        ro.r["library"]("BiocParallel")
        ro.r("BiocParallel::register(BiocParallel::MulticoreParam())")
        self.X_train = np.load(self.data)
        self.c_train = np.loadtxt(self.labels)
                
        # loading data
        ro.r(str("""fmat <- npyLoad("*""")[:-1] + self.data + str("""*", "integer")""")[1:])
        ro.r(str("""weight <- npyLoad("*""")[:-1] + self.weights + str("""*")""")[1:])
        ro.r(str("""cmat <- read.table("*""")[:-1] + self.labels + str("""*")""")[1:])
        ro.r("cmat$V2 <- factor(cmat$V1)")

    def fit(self):
        
        # computing data mask
        set_a = np.where(self.c_train == self.cluster[0])[0]
        subset_a = np.random.choice(set_a, self.A)
        set_b = np.where(self.c_train == self.cluster[1])[0]
        subset_b = np.random.choice(set_b, self.B)

        stochastic_set = np.hstack((subset_a, subset_b))
        f = np.array([a in stochastic_set for a in np.arange(self.X_train.shape[0])])

        nr, nc = f[:, np.newaxis].shape
        f_r = ro.r.matrix(f[:, np.newaxis], nrow=nr, ncol=nc)
        ro.r.assign("f_", f_r)
        ro.r("f <- as.integer(rownames(cmat[f_,]))")
        
        ro.r("local_fmat <- fmat[f, ]")
        ro.r("local_cmat <- cmat[f, ]")
        ro.r("local_weight <- weight[f, ]")
        ro.r("local_cmat$V3 <- factor(local_cmat$V1)")

        ro.r["library"]("zinbwave")
        ro.r["library"]("edgeR")
        ro.r("""dge <- DGEList(counts = t(local_fmat))""")
        ro.r("""dge <- suppressWarnings(edgeR::calcNormFactors(dge))""")
        ro.r("""design <- model.matrix(~V3, data = local_cmat)""")
        ro.r("""dge$weights <- t(local_weight)""")
        ro.r("""dge <- estimateDisp(dge, design)""")
        ro.r("""fit <- glmFit(dge, design)""")
        ro.r("""lrt <- glmWeightedF(fit, coef = 2)""")
        #print ro.r("lrt")
        return ro.r("lrt$table$PValue")
    
