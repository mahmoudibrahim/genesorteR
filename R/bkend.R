#rescale a vector so that the highest value is 1 and lowest is 0
score = function(x) {
	x  = ((x-min(x))/(max(x)-min(x)))
	return(x)
}

#binarize a gene expression matrix
binarize = function(x, method = "median") {

	if (method == "median") {
		pi = median(x@x)
	} else if (method == "mean") {
		pi = mean(x@x)
	} else if (method == "naive") {
		pi = min(x@x)
	} else if (method == "adaptiveMedian") {
		mm = Mclust(Matrix::rowSums(x), 1:20, modelNames=c("V"), verbose = FALSE)
		if (mm$G == 1) {
			stop("Error: you set binarizeMethod to adaptiveMedian but the optimal number of gene groups based on average expression is 1. Please use a different binarization method or use a different gene expression normalization. Also please consider reporting this error to mmibrahim@pm.me or on https://github.com/mahmoudibrahim/genesorteR/issues (preferred).")
		} else {
			pi = rep(0, mm$G)
			for (i in 1:mm$G) {
				pi[i] = median(x[which(mm$classification == i),]@x)
			}
		}
	} else if ((is.numeric(method)) & (method >= 0)) {
		pi = method
	} else {
		stop("Unrecognized binarization method! genesorteR stopped.")
	}


	if (method == "adaptiveMedian") {
	
		mat = list()
		for (i in 1:mm$G) {
			tx = x[which(mm$classification == i),]
			ww = which(tx@x >= pi[i])
			dp = diff(tx@p)
			colInd = (rep(seq_along(dp),dp))[ww]
			rowInd = (tx@i+1)[ww]
			genes = rownames(tx)
			mat[[i]] = Matrix::sparseMatrix(rowInd[1], colInd[1], dims=tx@Dim)
			mat[[i]][cbind(rowInd,colInd)] = 1
			mat[[i]] = as(mat[[i]], "dgCMatrix")
			rownames(mat[[i]]) = genes
		}
		mat = do.call(rbind, mat)
		x = mat[match(rownames(x),rownames(mat)),]
	
	} else {
	
		ww = which(x@x >= pi)
		dp = diff(x@p)
		colInd = (rep(seq_along(dp),dp))[ww]
		rowInd = (x@i+1)[ww]
		genes = rownames(x)
		x = Matrix::sparseMatrix(rowInd[1], colInd[1], dims=x@Dim)
		x[cbind(rowInd,colInd)] = 1
		x = as(x, "dgCMatrix")
		rownames(x) = genes
	}
	
	return(list(mat = x, cutoff = pi))
}


#get prior probability for each gene across the entire data set
getGeneProb = function(x) {
	ncell = ncol(x)

	probs = as.vector( (Matrix::rowSums(x)) / ncell )
	names(probs) = rownames(x)

	return(probs)
}


#get the prior probability of each cluster across the entire data set
getClassProb = function(x) {
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

	specScore = exp(log(postMat) + log(condMat))

	colnames(specScore) = colnames(condMat)
	specScore = as(specScore, "dgCMatrix")
	return(specScore)

}


#get scaled gene-cluster specificity score
getScaledSpecScore = function(postMat, condMat) {

	specScore = postMat * condMat
	specScore = apply(specScore, 2, score)
	
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


getClusterTFIDF = function(condMat, cores = 1) {
	numDoc = ncol(condMat)
	if (cores > 1) {
		idf = mclapply(1:nrow(condMat), function(x) log(numDoc) - log(length(which(condMat[x,] != 0))), mc.cores = cores)
	} else {
		idf = apply(condMat, 1, function(x) log(numDoc) - log(length(x[x!=0])))
	}
	tf_idf = exp(log(condMat) + log(idf))
	tf_idf = as(tf_idf, "dgCMatrix")
	return(tf_idf)
}

#get mutual information from two vectors
getMutInfo = function(x, y) {

	##check lengths of both vectors is equal
	if (length(y) == length(x)) {
		hx = getEntropy(table(x) / length(x))
		hy = getEntropy(table(y) / length(y))
		hxy = getEntropy(table(x,y) / length(y))

		mutInfo = hx + hy - hxy
	} else {
		stop("Couldn't calculate mutual information!")
	}
	
	return(mutInfo)
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
fastAUC = function(x, y) {
	da = ( (sum( (y - x) )) / length(x) ) * 2
	return(da)
}

#bin means of a vector, thanks to the mysterious user on stackoverflow: https://stackoverflow.com/a/43635971/2435654
binMean = function (vec, every, na.rm = FALSE) {
	n = length(vec)
	x = .colMeans(vec, every, n %/% every, na.rm)
	r = n %% every
	if (r) {
		x = c(x, mean.default(vec[(n - r + 1):n], na.rm = na.rm))
	}
	return(x)
}

#getClassIndeces
getClassIndeces = function(x) {
	uniq = unique(x)
	x = lapply(1:length(uniq), function(a) which(!is.na(match(x, uniq[a]))))
	return(x)
}
