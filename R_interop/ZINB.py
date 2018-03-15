import numpy as np
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils.validation import check_array, check_is_fitted
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import warnings
from rpy2.rinterface import RRuntimeWarning

class ZINB(BaseEstimator, TransformerMixin):
    
    def __init__(self, n_components=10, learn_V=True):
        self.n_components = n_components
        self.learn_V = learn_V
        warnings.filterwarnings("ignore", category=RRuntimeWarning)
        rpy2.robjects.numpy2ri.activate()
        ro.r["library"]("zinbwave")
        ro.r["library"]("BiocParallel")
        ro.r["library"]("scRNAseq")
        ro.r["library"]("matrixStats")
        ro.r["library"]("magrittr")
        ro.r["library"]("ggplot2")
        ro.r["library"]("biomaRt")
        ro.r["library"]("tibble")
        ro.r["library"]("SIMLR")
        ro.r("BiocParallel::register(BiocParallel::MulticoreParam())")
        ro.r.assign("K", n_components)
        ro.r['source']('~/R/zinboptimizeval.R')
        
    def fit(self, X):

        self.X_ = X
        nr,nc = X.shape
        X_trainr = ro.r.matrix(X, nrow=nr, ncol=nc)
        ro.r.assign("matrix_train", X_trainr)

        ro.r("m_train <- zinbModel(n=NROW(matrix_train), J=NCOL(matrix_train), K=K)")
        ro.r("m_train <- zinbInitialize(m_train, matrix_train)")
        if not self.learn_V:
            ro.r("m_train <- zinbOptimizeTrick(m_train, matrix_train, )")
        else:
            ro.r("m_train <- zinbOptimize(m_train, matrix_train)")
        self.params = ro.r("m_train")
        # Return the classifier
        return self
    
    def transform(self, X):

        # Check is fit had been called
        check_is_fitted(self, ['X_'])

        # Input validation
        X = check_array(X)
        nr,nc = X.shape
        X_testr = ro.r.matrix(X, nrow=nr, ncol=nc)
        ro.r.assign("matrix_val", X_testr)
        ro.r("m_val <- zinbModel(n=NROW(matrix_val), J=NCOL(matrix_val), K=K)")
        ro.r("m_val <- zinbInitialize(m_val, matrix_val)") 
        if not self.learn_V:
            ro.r("m_val <- zinbOptimizeValTrick(m_train, m_val, matrix_val)")
        else:
            ro.r("m_val <- zinbOptimizeVal(m_train, m_val, matrix_val)")
        return ro.r("getW(m_val)")
    
    def score(self, X):
        """
        Mean per sample likelihood of the data
        """
        X = check_array(X)
        nr,nc = X.shape
        X_testr = ro.r.matrix(X, nrow=nr, ncol=nc)
        ro.r.assign("matrix_val", X_testr)
        ro.r("m_val <- zinbModel(n=NROW(matrix_val), J=NCOL(matrix_val), K=K)")
        ro.r("m_val <- zinbInitialize(m_val, matrix_val)") 
        if not self.learn_V:
            ro.r("m_val <- zinbOptimizeValTrick(m_train, m_val, matrix_val)")
        else:
            ro.r("m_val <- zinbOptimizeVal(m_train, m_val, matrix_val)")
        return ro.r("zinbLikelihood(m_val, matrix_val)")[0] / X.shape[0]
    
    def output_estimation(self):
        """
        Returns parameters on the last validation set optimized
        """
        ro.r("fit <- zinbFit(m_val)")
        ro.r("mu <- fit[[1]]")
        ro.r("logitPi <- fit[[2]]")
        ro.r("theta <- fit[[3]]")
        return ro.r("mu"), ro.r("logitPi"), ro.r("theta")