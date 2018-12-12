# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

ola <- function() {
  print("Hello world")
}


randomclust <- function(dataset, classes=NULL, K=1,iterations=1,transformations=c("all"), algorithms=c("all")) {

  dataset = as.matrix(dataset)
  if(algorithms[1]=="all")
    algorithms <- c("hclust","kmeans","bclust","cba","hkclustering","dbscan")
  if(transformations[1]=="all"){
      if(ncol(dataset)<8)
        transformations = c("N100","RM100","M100","RUM100","UM100","R100","DA","DR")
      else
        transformations = c("N100","N80","N50","N20","RM100","RM20","RM50","RM80","M100","M20","M50","M80","RUM100","RUM20","RUM50","RUM80","UM100","UM20","UM50","UM80","R100","R80","R50","R20","DA","DR")
  }
  result_partial<-data.frame()
  for(i in 1:length(algorithms)){
    for(j in 1:length(transformations)){
      result_partial_aux <- run_clustering(dataset,classes,transformations[j],algorithm = algorithms[i],iterations,K)
      result_partial <- rbind(result_partial,result_partial_aux)
    }
  }

  result_partial=normalizeValues(result_partial,classes)

  result<-matrix(-1,length(transformations),length(algorithms)+1)
  rownames(result)=transformations
  colnames(result)=c(algorithms,"s")
  for(i in 1:length(transformations)){
    for(j in 1:length(algorithms)){
      x1<-transformations[i]
      x2<-algorithms[j]
      result[rownames(result)==x1,
             colnames(result)==x2]=
        sum(subset(result_partial,transformation==x1 & algorithm==x2)[,4])
    }
    result[i,length(algorithms)+1]=sum(result[i,1:length(algorithms)])
  }



  return(result)
}
