#binarize a gene expression matrix
binarize = function(x, method = "median") {
  ##todo: check if sparse or not
  #x = as.matrix(x)

  if (method == "median") {
    pi = median(x@x)
  } else if (method == "naive") {
    pi = min(x@x)
  } else {
    #todo: issue an error
  }


  ww = which(x@x >= pi)
  dp = diff(x@p)
  colInd = (rep(seq_along(dp),dp))[ww]
  rowInd = (x@i+1)[ww]
  genes = rownames(x)
  x = Matrix::sparseMatrix(rowInd[1], colInd[1], dims=x@Dim)
  x[cbind(rowInd,colInd)] = 1
  x = as(x, "dgCMatrix")
  rownames(x) = genes

  return(list(mat = x, cutoff = pi))
}



#get prior probability for each gene across the entire data set
getGeneProb = function(x) {

  ##todo: check if sparse or not. If not sparse, issue an error


  ncell = ncol(x)

  probs = as.vector( (Matrix::rowSums(x)) / ncell )
  names(probs) = rownames(x)

  return(probs)
}


#get the prior probability of each cluster across the entire data set
getClassProb = function(x) {

  #x = as.vector(x)

  probs = table(x) / length(x)

  return(probs)
}



#get gene probability conditional on cluster
getGeneConditionalCluster_stats = function(mat, classProb, classLabels, cores = 1) {

  lab = names(classProb)

  if (cores > 1) {
    geneProb = do.call(cbind, mclapply(1:length(lab), function(i) Matrix::rowMeans(mat[,which(classLabels == lab[i])]), mc.cores = cores))
  } else {
    geneProb = do.call(cbind, lapply(1:length(lab), function(i) Matrix::rowMeans(mat[,which(classLabels == lab[i])])))
  }

  geneProb = as(geneProb, "dgCMatrix")
  colnames(geneProb) = lab

  return(geneProb)

}


#get gene probability conditional on cluster
getClusterPostGene = function(condMat, geneProb, classProb) {

  postProb = t(apply(log(condMat), 1, function(x) x + log(classProb)))
  postProb = exp(apply(postProb, 2, function(x) x - log(geneProb)))

  colnames(postProb) = names(classProb)
  postProb = as(postProb, "dgCMatrix")

  return(postProb)

}


#get gene-cluster specificity score
getSpecScore = function(postMat, condMat) {

  specScore = postMat * condMat

  colnames(specScore) = colnames(condMat)
  specScore = as(specScore, "dgCMatrix")
  return(specScore)

}

#get entropy of a probability vector
getEntropy = function(x) {
  ent = 0
  w = which(x > 0)
  if (length(w) > 0) {
    w = x[w]
    ent = -(sum(w * (log(w))))
  }
  return(ent)
}

#get mutual information from two vectors
getMutInfo = function(x, y) {

  ##check lengths of both vectors is equal

  hx = getEntropy(table(x) / length(x))
  hy = getEntropy(table(y) / length(y))
  hxy = getEntropy(table(x,y) / length(y))

  mutInfo = hx + hy - hxy
  return(mutInfo)
}

#rescale a vector so that the highest value is 1 and lowest is 0
score = function(x) {
  x  = ((x-min(x))/(max(x)-min(x)))
  return(x)
}

#perform permutations
getPerma = function(gs, subsetCells = NULL, cores = 1) {

  if (is.null(subsetCells)) {
    shuff = base::sample(gs$inputClass, length(gs$inputClass), replace = FALSE)
    condGeneCluster = getGeneConditionalCluster_stats(gs$binary, gs$classProb, shuff, cores = cores)
    clusterPostGene = getClusterPostGene(condGeneCluster, gs$geneProb, gs$classProb)
    specScore = getSpecScore(clusterPostGene, condGeneCluster)

  } else {

    shuff = base::sample(gs$inputClass[subsetCells], length(gs$inputClass[subsetCells]), replace = FALSE)
    condGeneCluster = getGeneConditionalCluster_stats(gs$binary[,subsetCells], gs$classProb, shuff, cores = cores)
    clusterPostGene = getClusterPostGene(condGeneCluster, gs$geneProb, gs$classProb)
    specScore = getSpecScore(clusterPostGene, condGeneCluster)

  }
  return(specScore)
}

#return area under the curve
getAUC = function(x, y) {

  dX = c(diff(x), 0)
  dY = c(diff(y), 0)
  return(sum(x * dY) + sum(dX * dY)/2)

}
