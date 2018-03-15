import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import warnings
import numpy as np
from rpy2.rinterface import RRuntimeWarning
import scipy.sparse
import pandas as pd

class MAST(object):
    
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
        warnings.filterwarnings("ignore", category=RRuntimeWarning)
        rpy2.robjects.numpy2ri.activate()
        ro.r["library"]("RcppCNPy")
        ro.r["library"]("MAST")
        ro.r["library"]("BiocParallel")
        ro.r("BiocParallel::register(BiocParallel::MulticoreParam())")
        
        self.X_train = np.load(self.data)
        self.c_train = np.loadtxt(self.labels)
        
        # loading data
        ro.r(str("""fmat <- npyLoad("*""")[:-1] + self.data + str("""*", "integer")""")[1:])
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

        ro.r("local_fmat <- log2(fmat[f, ] + 1)")
        ro.r("local_cmat <- cmat[f, ]")
        ro.r("local_cmat$V3 <- factor(local_cmat$V1)")
        
        
        ro.r("sca <- FromMatrix(t(data.frame(local_fmat)), data.frame(local_cmat$V3))")
        ro.r("zlmCond <- zlm(~local_cmat.V3, sca)")
        ro.r("""summaryCond <- summary(zlmCond, doLRT='local_cmat.V34')""") 
        ro.r("summaryDt <- summaryCond$datatable")
        ro.r("""fcHurdle <- merge(
                summaryDt[contrast=='local_cmat.V34' & component=='H',.(primerid, `Pr(>Chisq)`)],
                    #hurdle P values
                summaryDt[contrast=='local_cmat.V34' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)],
                by='primerid') #logFC coefficients""")
        data = pd.DataFrame([ro.r("fcHurdle$primerid"), ro.r("""fcHurdle$'Pr(>Chisq)'"""), ro.r("fcHurdle$coef")]).T
        data.columns = ["gene_index", "p_value", "coeff"]
        data["gene_index"] = data["gene_index"].apply(lambda x: int(str(x)[1:]))
        data.sort_values("gene_index", inplace=True)
        if return_fc:
            return data["p_value"].values.astype(np.float), data["coeff"].values.astype(np.float)
        return data["p_value"].values.astype(np.float)