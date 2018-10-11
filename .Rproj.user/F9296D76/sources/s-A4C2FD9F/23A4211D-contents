
run_clustering <- function(dataset, classes, transformations, algorithm, iter=1){

  library(mclust)
  library(cluster)
  library(clValid)
  library(e1071)
  library(cba)
	'%ni%' <- Negate('%in%')
  N = nrow(dataset)
  K=length(unique(classes))
	partitionings = matrix(0, iter, N)

  j=1
	while(j<=iter){

		##MAHA NORMAL
		if(transformations=="RM100"){
			tdataset= mahalanobis(rotacao(dataset))
		}
		if(transformations=="RM20"){
			tdataset= select20_mahalanobis(rotacao(dataset))
		}
		if(transformations=="RM50"){
			tdataset= select50_mahalanobis(rotacao(dataset))
		}
		if(transformations=="RM80"){
			tdataset= select80_mahalanobis(rotacao(dataset))
		}
		##MAHA ROTACAO NORMAL
		if(transformations=="M100"){
			tdataset= mahalanobis(dataset)
		}
		if(transformations=="M20"){
			tdataset= select20_mahalanobis(dataset)
		}
		if(transformations=="M50"){
		 tdataset= select50_mahalanobis(dataset)
		}
		if(transformations=="M80"){
		 tdataset= select80_mahalanobis(dataset)
		}

		##MAHA UNIFORME
		if(transformations=="RUM100"){
		 tdataset= uniform_mahalanobis(rotacao(dataset))
		}
		if(transformations=="RUM20"){
		 tdataset= select20_unimahalanobis(rotacao(dataset))

		}
		if(transformations=="RUM50"){
		 tdataset= select50_unimahalanobis(rotacao(dataset))
		}
		if(transformations=="RUM80"){
		 tdataset= select80_unimahalanobis(rotacao(dataset))
		}

		if(transformations=="UM100"){
		 tdataset= uniform_mahalanobis(dataset)
		}
		if(transformations=="UM20"){
		 tdataset= select20_unimahalanobis(dataset)
		}
		if(transformations=="UM50"){
		 tdataset= select50_unimahalanobis(dataset)
		}
		if(transformations=="UM80"){
		 tdataset= select80_unimahalanobis(dataset)
		}

		##ROTACAO E NATURAL
		if(transformations=="R"){
		 tdataset= rotacao(dataset)
		}
		if(transformations=="N"){
		 tdataset= dataset
		}

		## Densidade
		if(transformations=="DA"){
		 tdataset= tadensidade(dataset,K)
		}
		if(transformations=="DR"){
		 tdataset= trdensidade(dataset,K)
		}

		dists = dist(tdataset)

		if(algorithm=="kmeans"){
    a <- NULL
		th <- 0
  		while(is.null(a) || th>10 ){
  				try(a <- kmeans(tdataset, iter.max=40, centers=K, nstart=1)$cluster)
  				th = th+1
  		}
  		result_partition = organizavetor(a)

  		partitionings[j,] = result_partition
		}

		if(algorithm=="hclust"){
		  a <- NULL
		  th <- 0
		  while(is.null(a) || th>10 ){
		    try(a <- hclust(dists, method='complete'))
		    th = th+1
		  }
		  result_partition = organizavetor(cutree(a,K))

		  partitionings[j,] = result_partition
		}

		if(algorithm=="bclust"){
		  a <- NULL
		  th <- 0
		  while(is.null(a) || th>10 ){
		    try(a <- bclust(tdataset,centers=K,verbose=FALSE)$cluster)
		    th = th+1
		  }
		  result_partition = organizavetor(a)

		  partitionings[j,] = result_partition
		}

		if(algorithm=="cba"){
		  a <- NULL
		  th <- 0
		  partition1 <- kmeans(tdataset, iter.max=40, centers=K, nstart=1)$cluster

		  partition_aux <- hclust(dists, method='complete')
		  partition2 <- organizavetor(cutree(partition_aux,K))
		  while(is.null(a) || th>10 ){
		    try(a <- cba(tdataset,partition1,partition2,K))
		    th = th+1
		  }
		  result_partition = organizavetor(a)

		  partitionings[j,] = result_partition
		}

		if(th>10){
			cat("\n threshold exceeded \n")
		}
		else {
			j=j+1
		}
	}

  ari = rep(0,iter)
  slh = rep(0,iter)
  dun = rep(0,iter)

  accuracy = 1-signif(avaliacao(partitionings, classes, K)/N, digits = 3)
  for(i in 1:iter){

    ari[i] = signif(adjustedRandIndex(partitionings[i,],classes), digits = 3)
		if(length(unique(partitionings[i,]))>1)
			slh[i] = signif(mean(silhouette(partitionings[i,], dists)[,3]), digits = 3)
		else
			slh[i] = 0
		if(length(unique(partitionings[i,]))==1)
			dun[i] = 0
		else{
			dun[i] = signif(dunn(dists,partitionings[i,]), digits = 3)
			if(dun[i]==Inf)
				dun[i] = 1
		}
  }

  avg_results <- NULL
	avg_results <- data.frame("algorithm"=rep(algorithm,4),
	                          "measure"=c("accuracy","ari","dunn","slh"),
	                          "transformation"=rep(transformations,4),
	                          "value"=c(mean(accuracy),mean(ari),mean(dun),mean(slh)))

	return(avg_results)

}




