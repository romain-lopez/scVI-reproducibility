
# S4 class file for Shekhar et al., "Comprehensive classification of retinal bipolar cells using single-cell transcriptomics", Cell, 2016

# Required packages
require(Matrix)
require(igraph)
require(gmodels)
require(ggplot2)
require(sva)
require(RANN)
require(reshape)

# Define slots
scDrop <- setClass("scDrop", slots = 
                     c(count.data="data.frame", data="data.frame",batch.data="data.frame", scale.data="matrix", 
                       group="vector", pca.load="data.frame",pca.scores="data.frame",
                       meta="data.frame",tsne.y="data.frame", cell.names="vector"))


# Initializes the S4 object
setGeneric("initialize", function(object,  min.cells=3, min.genes=2500, min.counts=10, scale=TRUE, center=TRUE, maxexpr=5000,...) standardGeneric("initialize"))
setMethod("initialize","scDrop",
          function(object,  min.cells=3, min.genes=2500, min.counts=10, scale=TRUE, center=TRUE, maxexpr=5000,...) {
            print("Initializing S4 object")
            
            # Cell filtering
            num.genes = colSums(object@count.data > 0); names(num.genes) = colnames(object@count.data)
            cells.use = names(num.genes[which(num.genes>min.genes)])
            temp.data=object@count.data[,cells.use]
            
            # Gene filtering
            num.cells= rowSums(temp.data > 0)            
            genes.use=names(num.cells[which(num.cells>min.cells)])
            genes.use1 = names(which(rowSums(object@count.data[,cells.use]) > min.counts))
            genes.use = intersect(genes.use, genes.use1)
            maxgene = apply(object@count.data[, cells.use],1, max)
            genes.use2 = rownames(object@count.data)[maxgene < maxexpr]
            genes.use=intersect(genes.use, genes.use2)
            temp.data=temp.data[genes.use,]
            
            #ADDED FILTERING OF MATRIX
            object@count.data = temp.data
            
            # Normalize each library to the median of the transcript counts across all cells 
            # Then, log transform expression values   
            print("Median normalizing counts and log-transforming")
            col.sums=apply(temp.data,2, sum)
            med_trans = median(col.sums)
            norm_counts = med_trans* scale(temp.data, center=FALSE, scale=col.sums)
            rm(temp.data)
            object@data=as.data.frame(log(norm_counts+ 1))
            
            # Group IDs for each cell
            object@group=factor(unlist(lapply(colnames(object@data),function(x) strsplit(x,"_")[[1]][1] )))
            names(object@group)=colnames(object@data)
            
            print("z-scoring each gene")
            object@cell.names=names(object@group)
            object@scale.data=t(scale(t(object@data),center=center,scale=scale))
            object@scale.data=object@scale.data[complete.cases(object@scale.data),]
            object@meta=data.frame(num.genes[object@cell.names]); colnames(object@meta)[1]="num.genes"
            object@meta[cells.use,"num.trans"] = colSums(object@count.data[,cells.use])
            
            
            object@meta[names(object@group),"sample"]=object@group
            return(object)
          }         
)

setGeneric("doBatchCorrection", function(object,  batch.cov=NULL, max.val=6 ) standardGeneric("doBatchCorrection"))
# Batch Correction using ComBat
# Memory heavy - consider running on a cluster
setMethod("doBatchCorrection","scDrop",
          function(object,   batch.cov=NULL, max.val=6) {
            correct.data = ComBat(object@data,batch.cov, prior.plots=FALSE, par.prior=TRUE)
            correct.data[correct.data > max.val] = max.val
            
            object@batch.data = as.data.frame(correct.data)
            rm(correct.data)
            object@scale.data = t(scale(t(object@batch.data), center=TRUE, scale=TRUE))
            return(object)
          }
)

setGeneric("doPCA", function(object,pcs.store=100) standardGeneric("doPCA"))
setMethod("doPCA", "scDrop", 
          function(object,pcs.store=100) {
            data.use=object@scale.data
            pc.genes = rownames(object@scale.data)
            
            #Remove genes with zero variation
            pc.genes.var = apply(data.use[pc.genes,],1,function(x) var(x))
            genes.use = pc.genes[pc.genes.var>0]
            pc.data = data.use[genes.use,]
            
            pca.obj = fast.prcomp(t(pc.data),center=FALSE, scale=FALSE)
            object@pca.scores=data.frame(pca.obj$x[,1:pcs.store])
            object@pca.load=data.frame(pca.obj$rotation[,1:pcs.store])
            return(object)
          }
)

setGeneric("violinplot", function(object,genes.plot,nCol=NULL,ylab.max=12, ...)  standardGeneric("violinplot"))
setMethod("violinplot","scDrop",
          function(object,genes.plot,nCol=NULL,ylab.max=12, ...) {
            if (is.null(nCol)) {
              nCol=1
              if (length(genes.plot)>6) nCol=3
              if (length(genes.plot)>9) nCol=4
            }
            genes.plot1 =intersect(genes.plot, rownames(object@data))
            plot.data = data.frame()
            if (length(genes.plot1) > 0) {
              plot.data = object@data[genes.plot1,]
            }
            
            genes.plot2 = intersect(genes.plot, c("num.genes", "num.trans"))
            if (length(genes.plot2) > 0) {
              plot.data = rbind(plot.data, t(data.frame(object@meta[,genes.plot2])))
            }
            genes.plot = intersect(genes.plot, rownames(plot.data))
            group.use=object@group
            
            #Make individual violins
            pList=lapply(genes.plot,function(x) makeviolin(x,plot.data[x,],group.use))
            
            multiplotList(pList,cols = nCol)
            rp()
            
          }
)  

makeviolin=function(x,data,group) {
  data$x=as.character(rownames(data))
  data.use=data.frame(data[x,])
  group = factor(group, levels=sort(unique(group)))
  if (length(x)==1) {
    data.melt=data.frame(rep(x,length(group))); colnames(data.melt)[1]="x"
    data.melt$value=as.numeric(data[1,1:length(group)])
    data.melt$id=names(data)[1:length(group)]
  }
  
  if (length(x)>1) data.melt=melt(data.use,id="x")
  data.melt$group=group
  
  noise <- rnorm(length(data.melt$value))/100000
  data.melt$value=as.numeric(as.character(data.melt$value))+noise
  p=ggplot(data.melt,aes(factor(group),value))
  p2=p + geom_violin(scale="width",adjust=1,trim=TRUE,aes(fill=factor(group))) + ylab(x) + xlab("Group/Cluster")
  
  p3=p2+theme(legend.position="top")+guides(fill=guide_legend(title=NULL))+geom_jitter(height=0,size=0.7)
  p4=p3+theme_bw() + ylim(0.8*min(data.melt$value), quantile(data.melt$value,0.99))
  p5=(p4+theme(axis.title.y = element_text(face="bold", colour="black", size=10), axis.text.y  = element_text(angle=90, size=8),
               axis.title.x = element_text(face="bold", colour="black", size=10), axis.text.x  = element_text(size=8)))
  
  return(p5)
}

multiplotList <- function(plots, file, cols=1, layout=NULL) {
  require(grid)
  
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

rp=function() {par(mfrow=c(1,1))}


# Graph clustering
setGeneric("doGraph_clustering", function(object,cells.use=NULL,pcs.use=1:10, num.nn=30, do.jaccard=FALSE, method="Louvain") standardGeneric("doGraph_clustering"))
setMethod("doGraph_clustering", "scDrop", function(object,cells.use=NULL,pcs.use=1:10,num.nn=30, do.jaccard=FALSE, method="Louvain") {
  
  
  if (do.jaccard){
    weights=TRUE;
    method_print = paste0(method,"-","Jaccard")
  } else {
    weights=NULL;
    method_print = method
  }
  
  print(paste0("Performing ", method_print, " clustering. Using ", num.nn, " nearest neighbors, and ", max(pcs.use), " PCs"))
  
  if (is.null(cells.use)){
    data.use=object@pca.scores[,pcs.use]
  } else {
    data.use=object@pca.scores[cells.use,pcs.use]
  } 
  
  Adj = get_edges(data.use,nn=num.nn,do.jaccard=do.jaccard)
  
  
  g=graph.adjacency(Adj, mode = "undirected", weighted=weights)
  if (method=="Louvain") graph.out = cluster_louvain(g)
  if (method=="Infomap") graph.out = cluster_infomap(g)
  
  clust.assign = factor(graph.out$membership, levels=sort(unique(graph.out$membership)))
  names(clust.assign) = graph.out$names
  k=order(table(clust.assign), decreasing = TRUE)
  new.levels = rep(1,length(unique(graph.out$membership)))
  new.levels[k] = 1:length(unique(graph.out$membership))
  levels(clust.assign) = new.levels
  clust.assign = factor(clust.assign, levels=1:length(unique(graph.out$membership)))
  print("Outputting clusters ..")
  object@meta$clust = NULL
  object@meta[names(clust.assign),"clust"]=clust.assign
  object@group=clust.assign; names(object@group)=names(clust.assign);               
  
  
  return(object) 
}
)

# Build a nearest neighbor graph with or without edge weights, and return an adjacency matrix
get_edges=function(X,nn=30,do.jaccard=TRUE) {
  nearest=nn2(X,X,k=nn+1, treetype = "bd", searchtype="priority")
  print("Found nearest neighbors")
  nearest$nn.idx = nearest$nn.idx[,-1]
  nearest$nn.dists = nearest$nn.dists[,-1] #Convert to a similarity score
  nearest$nn.sim = 1*(nearest$nn.dists >= 0 )
  
  
  edges = melt(t(nearest$nn.idx)); colnames(edges) = c("B", "A", "C"); edges = edges[,c("A","B","C")]
  edges$B = edges$C; edges$C=1
  
  #Remove repetitions
  edges = unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))
  
  if (do.jaccard){
    
    NN = nearest$nn.idx
    jaccard_dist = apply(edges, 1, function(x) length(intersect(NN[x[1], ],NN[x[2], ]))/length(union(NN[x[1], ], NN[x[2], ])) )
    
    edges$C = jaccard_dist
    edges = subset(edges, C != 0)
    edges$C = edges$C/max(edges$C)
  }
  
  
  
  
  Adj = matrix(0, nrow=nrow(X), ncol=nrow(X))
  rownames(Adj) = rownames(X); colnames(Adj) = rownames(X)
  Adj[cbind(edges$A,edges$B)] = edges$C
  Adj[cbind(edges$B,edges$A)] = edges$C
  return(Adj)
  
}

#Visualize single cells as clusters in a tSNE plot
plot.tsne=function(object) {
  
  cols=rainbow(length(levels(object@group))); cols[1]="lightgrey"
  group=as.numeric(object@group)
  
  order = sample(c(1:dim(object@tsne.y)[1]), replace=FALSE)
  plot(object@tsne.y[,1],object@tsne.y[,2],col=cols[as.integer(object@group)],pch=16,xlab="tSNE1",ylab="tSNE2",cex=0.3)
  k.centers=t(sapply(levels(object@group),function(x) apply(object@tsne.y[names(object@group[which(object@group %in% x)]),],2,mean)))
  points(k.centers[,1],k.centers[,2],cex=1.3,col="white",pch=16); text(k.centers[,1],k.centers[,2],levels(object@group),cex=1.25)
  
}

# scatter plot coloring genes by their expression levels
setGeneric("gene.expression.scatter", function(object, genes,cells.use=NULL,cols.use=terrain.colors(10),pch.use=16,nCol=NULL, xlim=NULL, ylim=NULL) standardGeneric("gene.expression.scatter"))
setMethod("gene.expression.scatter", "scDrop", 
          function(object, genes,cells.use=NULL,cols.use=terrain.colors( 10),pch.use=16,nCol=NULL, xlim=NULL, ylim=NULL) {
            
            if (is.null(cells.use)){
              cells.use = colnames(object@data)
            } else {
              cells.use = cells.use[cells.use %in% colnames(object@data)]
            }
            
            
            if (is.null(nCol)) {
              nCol=2
              if (length(genes)>6) nCol=3
              if (length(genes)>9) nCol=4
            }         
            num.row=floor(length(genes)/nCol-1e-5)+1
            
            par(mfrow=c(num.row,nCol))
            data.plot=object@tsne.y[cells.use,]
            
            group.id=as.factor(object@group[cells.use])
            data.plot$group=group.id
            x1="tSNE1"; x2="tSNE2"
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
            data.plot$pt.size=1
            for(i in genes) {
              data.gene=as.numeric(object@data[i,cells.use])
              data.cut=as.numeric(as.factor(cut(data.gene,breaks = length(cols.use))))
              data.col=rev(cols.use)[data.cut]
              plot(data.plot$x,data.plot$y,col=data.col,cex=0.7,pch=16,main=i,xlab=x1,ylab=x2, xlim=xlim, ylim=ylim)
            }
            rp()
          }
)

# Binomial test to evaluate differentially expressed genes between two clusters
# if only one cluster is provided, then it will be compared against the rest of the cells
setGeneric("markers.binom", function(object, clust.1,clust.2=NULL,effect.size=log(2), TPM.mat=NULL, Count.mat=NULL) standardGeneric("markers.binom"))
setMethod("markers.binom", "scDrop",
          function(object, clust.1,clust.2=NULL,effect.size=log(2), TPM.mat=NULL, Count.mat=NULL) {
            genes.use=rownames(object@data)
            clust.use=object@group
            cells.1=names(clust.use[which(clust.use%in%clust.1)])
            
            if (is.null(clust.2)) {
              clust.2="rest"
              cells.2=names(clust.use)
              cells.2=cells.2[!(cells.2%in%cells.1)]
            } else {
              cells.2=names(clust.use[which(clust.use%in%clust.2)])
            }
            
            Count.mat = object@count.data
            if (is.null(TPM.mat)) TPM.mat = exp(object@data[, c(cells.1, cells.2)])-1
            if (is.null(Count.mat)) Count.mat = object@count.data[genes.use, c(cells.1, cells.2)]
            result=binomcount.test(object, cells.1, cells.2, effect.size, TPM.mat, Count.mat)
            
            
            posFrac.1 = apply(object@data[rownames(result),cells.1],1,function(x) round(sum(x > 0)/length(x),2))
            posFrac.2 = apply(object@data[rownames(result),cells.2],1,function(x) round(sum(x > 0)/length(x),2))
            
            if (clust.2=="rest"){
              genes.include = posFrac.1 >= 0.1
            } else{
              genes.include = (posFrac.1 >= 0.1) | (posFrac.2 >= 0.1)
            }
            
            result = result[genes.include,]
            result = result[order(abs(result$log.effect), decreasing=TRUE),]
            
            #Mean number of transcripts per cell
            if (!is.null(attr(object,"count.data"))){
              nTrans.1 = apply(object@count.data[rownames(result), cells.1], 1, function(x) round(mean(x),3))
              nTrans.2 = apply(object@count.data[rownames(result), cells.2], 1, function(x) round(mean(x),3))
              result[,paste0("nTrans_", clust.1)] = nTrans.1
              result[, paste0("nTrans_", clust.2)] = nTrans.2
            }
            
            return(result)
          } 
)

setGeneric("binomcount.test", function(object, cells.1,cells.2, effect.size, TPM.mat, Count.mat) standardGeneric("binomcount.test"))
setMethod("binomcount.test", "scDrop",
          function(object, cells.1,cells.2, effect.size, TPM.mat, Count.mat) {
            
            x=TPM.mat
            y=Count.mat
            
            #Test for enrichments in cluster #1
            m = apply(x[, cells.2], 1, function(x) sum(x>0)) #Number of cells expressing marker in cluster #2
            m1 = m; m1[m==0]=1; # Regularization. Add a pseudo count of 1 to unexpressed markers to avoid false positives
            n = apply(x[, cells.1], 1, function(x) sum(x>0)) #Number of cells expressing marker in cluster #1
            #Find the probability of finding n or more cells +ve for marker in cluster 1 given the fraction in cluster 2
            pv1 = pbinom(n, length(cells.1), m1/length(cells.2), lower.tail = FALSE) + dbinom(n, length(cells.1), m1/length(cells.2))
            
            log_fold_express = log(n*length(cells.2)/(m*length(cells.1))) #log proportion of expressing cells
            d1 <- data.frame(log.effect=log_fold_express,pval=pv1)
            d1 <- subset(d1, log.effect >= effect.size)
            d1 <- d1[order(d1$pval,decreasing=FALSE),]
            
            #Enrichments in cells.2
            n1 = n; n1[n==0]=1; # Regularization.  Add a pseudo count of 1 to unexpressed markers to avoid false positives
            #Find the probability of finding m or more cells +ve for marker in cluster 2 given the fraction in cluster 1
            pv2 = pbinom(m, length(cells.2), n1/length(cells.1), lower.tail=FALSE) + dbinom(m, length(cells.2), n1/length(cells.1))
            d2 <- data.frame(log.effect=log_fold_express,pval=pv2)
            d2 <- subset(d2, log.effect <= -effect.size)
            d2 <- d2[order(d2$pval,decreasing=FALSE),]
            
            d = rbind(d1, d2);
            d = d[order(d$pval, decreasing=FALSE),]
            return(d)
          } 
)




setGeneric("merge.clusters.DE", function(object, min.de.genes = 25, effect.size=log(2), pval.cutoff=0.01, pcs.use=1:10,TPM.mat=NULL, Count.mat=NULL) standardGeneric("merge.clusters.DE"))
setMethod("merge.clusters.DE", "scDrop", 
          function(object,min.de.genes=25, effect.size=log(2),pval.cutoff=0.01, pcs.use=1:10, TPM.mat=NULL, Count.mat=NULL) {
            
            genes.use = rownames(object@data)
            clust.test = as.numeric(levels(object@group))
            # if (is.null(tag)){
            #   filename = "CLUSTER_PAIRWISE_MARKERS.txt"
            # } else {
            #   filename = paste0("CLUSTER_PAIRWISE_MARKERS_",tag,".txt")
            # }
            # zz = file(filename,open="wt")
            
            
            num.clust=length(clust.test) 
            print(paste0("Starting with ", num.clust, " clusters"))
            
            pass.thresh=1e6*data.frame(diag(length(levels(object@group)))); 
            
            for (i in setdiff(as.numeric(levels(object@group)), clust.test)){
              pass.thresh[i,]=1e6; pass.thresh[,i]=1e6;
            } 
            
            dist.clust = pass.thresh
            
            #Find the number of differentially expressed genes between every pair of clusters
            for(k in 1:num.clust) {
              i=clust.test[k]
              print(paste0("Testing Cluster ", i))
              for(m in ((k+1):num.clust)) {
                j=clust.test[m]
                #print(j)
                if (m>num.clust) break
                if (pass.thresh[i,j]==0) {
                  marker=markers.binom(object,i,j,effect.size=effect.size, TPM.mat=TPM.mat, Count.mat=Count.mat)
                  P2 = p.adjust(marker$pval, method="fdr")
                  marker$pval = P2
                  
                  marker.pass=subset(marker,pval<pval.cutoff)
                  #print(paste("Test b/w Clusters ",i,j, "-# DE = ",nrow(marker.pass)))
                  #print(head(subset(marker.pass, log.effect > 0),5))
                  #print(head(subset(marker.pass, log.effect < 0),5))
                  
                  num.de.genes = 2*min(nrow(subset(marker.pass, log.effect > 0)), nrow(subset(marker.pass, log.effect < 0)))
                  pass.thresh[i,j]=num.de.genes; pass.thresh[j,i]=pass.thresh[i,j];
                  
                }
                
              }
              
              print(pass.thresh[i,])
              
            }
            
            colnames(pass.thresh) = levels(object@group)
            rownames(pass.thresh) = levels(object@group)
            
            write.table(pass.thresh, file=paste0("DE_genes_matrix_2.txt"), sep="\t", quote=FALSE)
            
            #iteratively merge clusters
            min.val = min(pass.thresh)
            min.val.ind = as.data.frame(which(pass.thresh==min.val, arr.ind = TRUE))
            min.val.ind = min.val.ind[min.val.ind$row < min.val.ind$col,]
            min.val.ind$val = min(pass.thresh)
            rownames(min.val.ind) = 1:nrow(min.val.ind)
            
            merge.ind=-1
            while(min.val <= min.de.genes) {
              merge.ind=merge.ind+1
              
              #In case of ties, merge clusters that are closest in PC space
              clust.dists = ComputeClusterDistances(object, reduction.use="pca", dist.type="centroid", pcs.use=pcs.use) 
              ind.min=which.min(clust.dists[cbind(min.val.ind$row, min.val.ind$col)])
              test.1 = min.val.ind[ind.min,]$row; test.2 = min.val.ind[ind.min,]$col
              
              if (pass.thresh[test.1,test.2]<= min.de.genes) {
                object@group[which(object@group==test.2)]=test.1
                pass.thresh = pass.thresh[-test.2,]; pass.thresh = pass.thresh[,-test.2]
                old.group.levels = as.numeric(levels(object@group))
                old.group.levels = setdiff(old.group.levels, test.2)
                clust.test = setdiff(clust.test, test.2)
                
                object@group = droplevels(object@group)
                levels(object@group) = c(1:length(levels(object@group)))
                object@meta[,"clust"] = object@group
                
                new.group.levels = as.numeric(levels(object@group))
                names(new.group.levels) = as.character(old.group.levels)
                clust.test = new.group.levels[as.character(clust.test)]
                
                
                
                #Recompute pairwise markers for merged cluster
                print(paste0("Recomputing pairwise markers for new clust ", test.1))
                for (i in setdiff(clust.test, test.1)){
                  print(i)
                  marker= markers.binom(object,test.1,i,effect.size=effect.size, TPM.mat=TPM.mat, Count.mat=Count.mat)
                  P2 = p.adjust(marker$pval, method="fdr")
                  marker$pval = P2
                  marker.pass=subset(marker,pval<pval.cutoff)
                  pass.thresh[test.1,i]=2*min(nrow(subset(marker.pass, log.effect>0)),nrow(subset(marker.pass, log.effect<0))); pass.thresh[i,test.1]=pass.thresh[test.1,i];
                  #pass.thresh[test.1,i]=nrow(marker.pass); 
                  pass.thresh[i,test.1]=pass.thresh[test.1,i];
                  
                }
                
              }
              colnames(pass.thresh) = 1:length(levels(object@group))
              rownames(pass.thresh) = colnames(pass.thresh)
              
              min.val = min(pass.thresh)
              min.val.ind = as.data.frame(which(pass.thresh==min.val, arr.ind = TRUE))
              min.val.ind = min.val.ind[min.val.ind$row < min.val.ind$col,]
              min.val.ind$val = min(pass.thresh)
              rownames(min.val.ind) = 1:nrow(min.val.ind)
              
            }
            return(object)
          }
)

setGeneric("ComputeClusterDistances", function(object, reduction.use="pca", dist.type="nn", pcs.use = 1:10, cells.use=NULL, genes.use=NULL) standardGeneric("ComputeClusterDistances"))
setMethod("ComputeClusterDistances", "scDrop", 
          function(object,reduction.use="pca",dist.type="nn", pcs.use = 1:10, cells.use=NULL, genes.use=NULL) {
            cells.use =  colnames(object@data)
            group.use=object@group[cells.use]
            if (reduction.use == "pca"){
              data.use = object@pca.scores[cells.use,pcs.use]
              centroids = ClusterCentroids(object, reduction.use="pca", pcs.use=pcs.use, cells.use=cells.use)
            }
            
            if (dist.type=="centroid"){
              clust.dists = as.matrix(dist(centroids, upper=TRUE))
              diag(clust.dists) = 1e6
            }
            
            num.clust = length(levels(group.use))
            
            
            if (dist.type == "nn"){
              clust.dists = matrix(0, nrow=num.clust, ncol=num.clust)
              diag(clust.dists) = 1e6
              rownames(clust.dists) = levels(group.use)
              colnames(clust.dists) = rownames(clust.dists)
              for (i in 1:nrow(clust.dists)){
                for(j in ((i+1):ncol(clust.dists))){
                  if (j>nrow(clust.dists)) break
                  cells.in.cluster_i = names(object@group)[object@group %in% i]
                  cells.in.cluster_i = cells.in.cluster_i[cells.in.cluster_i %in% cells.use]
                  cells.in.cluster_j = names(object@group)[object@group %in% j]
                  cells.in.cluster_j = cells.in.cluster_j[cells.in.cluster_j %in% cells.use]
                  
                  nnA = nn2(data.use[cells.in.cluster_i,], query = centroids[j,], k=1)
                  nnB = nn2(data.use[cells.in.cluster_j,], query = centroids[i,],k=1)
                  clust.dists[i,j] = min(c(nnA$nn.dists, nnB$nn.dists))
                  
                  clust.dists[j,i] = clust.dists[i,j]
                }
              }
            }
            
            colnames(clust.dists) = c(1:ncol(clust.dists))
            rownames(clust.dists) = colnames(clust.dists)
            return(clust.dists)
          }
          
          
)

setGeneric("ClusterCentroids", function(object,reduction.use="pca", pcs.use = 1:10, cells.use=NULL, genes.use=NULL) standardGeneric("ClusterCentroids"))
setMethod("ClusterCentroids", "scDrop", 
          function(object,reduction.use="pca", pcs.use = 1:10, cells.use=NULL, genes.use=NULL) {
            cells.use = colnames(object@data)
            group.use=object@group[cells.use]
            if (reduction.use == "pca"){
              data.use = object@pca.scores[cells.use,pcs.use]
            }
            
            centroids = c()
            for (i in levels(group.use)){
              cells.in.cluster = names(object@group)[object@group %in% i]
              cells.in.cluster = cells.in.cluster[cells.in.cluster %in% cells.use]
              centroids = rbind(centroids, colMeans(data.use[cells.in.cluster,]))
            }
            centroids = as.data.frame(centroids)
            colnames(centroids) = colnames(data.use)
            rownames(centroids) = as.numeric(levels(object@group))
            
            return(centroids)
          }
          
          
)

setGeneric("dot.plot", function(object,features.use=NULL, group.use=NULL, group.names=NULL, thresh.use=0,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,max.size=10,min.perc=0,...) standardGeneric("dot.plot"))
setMethod("dot.plot", "scDrop", 
          function(object,features.use=NULL,group.use=NULL, group.names = NULL,thresh.use=0,do.transpose=FALSE,max.val.perc=NULL, max.val.exp=NULL,max.size=10,min.perc=0,...) {
            
            
            features.use=features.use[features.use %in% rownames(object@data)]
            if (is.null(group.use)) group.use = levels(object@group)
            if (is.null(group.names)) group.names = group.use
            
            if (length(group.names) != length(group.use)){
              print("Error : group.names must be of the same length as the groups.use/ number of clusters. Using cluster numbers as labels ...")
              group.names = group.use
            }
            
            #Initialize matrix of percent expressing cells
            PercMat = matrix(0, nrow=length(features.use), ncol = 0)
            rownames(PercMat) = features.use; 
            
            #Initialize matrix of average transcript levels
            ExpMat = PercMat;
            
            #Count mat
            Count.mat = object@count.data[features.use, colnames(object@data)]
            
            
            for (i in group.use){
              cells.in.cluster = names(object@group)[which(object@group== i)]
              vec.exp = apply(object@data[features.use, cells.in.cluster], 1, function(x) sum(x>thresh.use)/length(x)) 
              PercMat = cbind(PercMat,vec.exp)
              
              vec.exp = apply(Count.mat[features.use, cells.in.cluster], 1, function(x) if (sum(x>0) > 1){ mean(x[x>0]) } else {sum(x)})
              ExpMat = cbind(ExpMat, vec.exp)
            }
            colnames(ExpMat) = group.names
            colnames(PercMat) = group.names
            
            
            
            rows.use = rownames(PercMat)[apply(PercMat, 1, function(x) max(x) >= min.perc)];
            PercMat = PercMat[rows.use,]
            ExpMat = ExpMat[rows.use,]
            features.use = rows.use
            if (!is.null(max.val.perc)) PercMat[PercMat > max.val.perc] = max.val.perc
            if (!is.null(max.val.exp)) ExpMat[ExpMat > max.val.exp] = max.val.exp
            
            
            ExpVal = melt(ExpMat)
            PercVal = melt(PercMat)
            colnames(ExpVal) = c("gene","cluster","nTrans")
            ExpVal$percExp = PercVal$value*100
            
            if (!do.transpose){
              ExpVal$gene = factor(ExpVal$gene, levels=features.use)
              ExpVal$cluster = factor(ExpVal$cluster, levels= rev(group.names))
              p=ggplot(ExpVal, aes(y = factor(cluster),  x = factor(gene))) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                scale_color_gradient(low ="blue",   high = "red", limits=c( 1, max(ExpVal$nTrans) ))+scale_size(range = c(0, max.size))+   theme_bw() +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
              p = p + ylab("Cluster") + xlab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", angle=45, hjust=1)) + 
                theme(axis.text.y=element_text(size=12, face="italic"))
              print(p)
            } else {
              ExpVal$gene = factor(ExpVal$gene, levels=rev(features.use))
              ExpVal$cluster = factor(ExpVal$cluster, levels= group.names)
              p=ggplot(ExpVal, aes(y = factor(gene),  x = factor(cluster))) + geom_point(aes(colour = nTrans,  size =percExp)) + 
                scale_color_gradient(low ="blue",   high = "red", limits=c(1, max(ExpVal$nTrans) ))+scale_size(range = c(0, max.size))+   theme_bw() +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
              p = p + xlab("Cluster") + ylab("Gene") + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
                theme(axis.text.y=element_text(size=12, face="italic")) 
              print(p)
              
              
            }
            
          }
)
