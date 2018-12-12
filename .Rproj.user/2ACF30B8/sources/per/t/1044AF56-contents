
run_clustering <- function(dataset, classes, transformations, algorithm, iter,K){

  suppressMessages(library(mclust))
  suppressMessages(library(cluster))
  suppressMessages(library(clValid))


	'%ni%' <- Negate('%in%')
  N = nrow(dataset)
  if(is.null(classes)==FALSE){
    K=length(unique(classes))
  }
	partitionings = matrix(0, iter, N)

  j=1
	while(j<=iter){

		tdataset<-transformdata(transformations)

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
		  suppressMessages(library(e1071))
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
		  suppressMessages(library(cba))
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

		if(algorithm=="hkclustering"){
		  suppressMessages(library(hkclustering))
		  a <- NULL
		  th <- 0
		  while(is.null(a)){
		    try(a <- hkclustering(as.data.frame(tdataset),K,50)$cluster_number)
		    th = th+1
		  }
		  result_partition = organizavetor(a)

		  partitionings[j,] = result_partition
		}

		if(algorithm=="dbscan"){
		  suppressMessages(library(dbscan))
		  suppressMessages(library(SiZer))
		  a <- NULL
		  th <- 0
		  while(is.null(a)){
		    knndist<-sort(kNNdist(tdataset,K))
		    pl<-piecewise.linear(1:length(knndist),knndist)
		    eps<-knndist[as.integer(pl$change.point)]
		    try(a <- dbscan(tdataset,eps)$cluster)
		    result_partition = organizavetor(a)

		    partitionings[j,] = result_partition
		  }
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

  if(is.null(classes)==FALSE)
    accuracy = 1-signif(avaliacao(partitionings, classes, K)/N, digits = 3)

  for(i in 1:iter){
    if(is.null(classes)==FALSE)
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

  if(is.null(classes)==FALSE){
    measures = c("accuracy","ari","dunn","slh")
    values = c(mean(accuracy),mean(ari),mean(dun),mean(slh))
  }
  else{
    measures = c("dunn","slh")
    values = c(mean(dun),mean(slh))
  }
  avg_results <- NULL
	avg_results <- data.frame("algorithm"=rep(algorithm,length(measures)),
	                          "measure"=measures,
	                          "transformation"=rep(transformations,length(measures)),
	                          "value"=values)

	return(avg_results)

}




