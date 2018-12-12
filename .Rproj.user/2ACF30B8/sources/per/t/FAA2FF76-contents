random_rotation_matrix_incl_flip <- function (n){
  QR <- qr(matrix(rnorm(n^2),ncol=n))
  M <- qr.Q(QR) %*% diag(sign(diag(qr.R(QR))))
  return (M)
}

random_rotation_matrix <- function (n){
  M <- random_rotation_matrix_incl_flip(n)
  if(det(M)<0) M[,1] <- -M[,1]
  return (M)
}

mahalanobis <- function(docs) {

  n = nrow(docs)
  m = ncol(docs)
  rdocs = matrix(0,n,m)

  maha<-matrix(rnorm(m^2,mean=0,sd=1),ncol=m)
  maha<- maha**2
  for (i in 1:m){maha[i,i]=1}

  rdocs <- docs%*%maha

  rdocs
}

uniform_mahalanobis <- function(docs) {

  n = nrow(docs)
  m = ncol(docs)
  rdocs = matrix(0,n,m)

  maha<-matrix(runif(m^2),ncol=m)

  for (i in 1:m){maha[i,i]=1}

  rdocs <- docs%*%maha

  rdocs
}

select20_mahalanobis <- function(docs) {

  n = nrow(docs)
  m = ncol(docs)
  rdocs = matrix(0,n,m)

  maha<-matrix(runif(m^2),ncol=m)

  s = sample(1:m, floor(m-m*.2), replace=FALSE)

  for(i in 1:length(s)){maha[,s[i]]=0}
  for (i in 1:m){maha[i,i]=1}

  rdocs <- docs%*%maha

  rdocs
}

select50_mahalanobis <- function(docs) {

  n = nrow(docs)
  m = ncol(docs)
  rdocs = matrix(0,n,m)

  maha<-matrix(runif(m^2),ncol=m)

  s = sample(1:m, floor(m-m*.5), replace=FALSE)

  for(i in 1:length(s)){maha[,s[i]]=0}
  for (i in 1:m){maha[i,i]=1}

  rdocs <- docs%*%maha
  rdocs
}

select80_mahalanobis <- function(docs) {

  n = nrow(docs)
  m = ncol(docs)
  rdocs = matrix(0,n,m)

  maha<-matrix(runif(m^2),ncol=m)

  s = sample(1:m, floor(m-m*.8), replace=FALSE)

  for(i in 1:length(s)){maha[,s[i]]=0}
  for (i in 1:m){maha[i,i]=1}

  rdocs <- docs%*%maha

  rdocs
}

select20_unimahalanobis <- function(docs) {

  n = nrow(docs)
  m = ncol(docs)
  rdocs = matrix(0,n,m)

  maha<-matrix(runif(m^2),ncol=m)

  s = sample(1:m, floor(m-m*.2), replace=FALSE)

  for(i in 1:length(s)){maha[,s[i]]=0}
  for (i in 1:m){maha[i,i]=1}

  rdocs <- docs%*%maha

  rdocs
}

select50_unimahalanobis <- function(docs) {

  n = nrow(docs)
  m = ncol(docs)
  rdocs = matrix(0,n,m)

  maha<-matrix(runif(m^2),ncol=m)

  s = sample(1:m, floor(m-m*.5), replace=FALSE)

  for(i in 1:length(s)){maha[,s[i]]=0}
  for (i in 1:m){maha[i,i]=1}

  rdocs <- docs%*%maha

  rdocs
}

select80_unimahalanobis <- function(docs) {

  n = nrow(docs)
  m = ncol(docs)
  rdocs = matrix(0,n,m)

  maha<-matrix(runif(m^2),ncol=m)

  s = sample(1:m, floor(m-m*.8), replace=FALSE)

  for(i in 1:length(s)){maha[,s[i]]=0}
  for (i in 1:m){maha[i,i]=1}


  rdocs <- docs%*%maha

  rdocs
}

rotacao <- function(docs){

  N = nrow(docs)
  M = ncol(docs)
  rdocs = matrix(0,N,M)

  rmatrix = random_rotation_matrix(M)

  rdocs<-docs%*%rmatrix

  return(rdocs)

}

trdensidade <- function(docs,k){

  N = nrow(docs)
  M = ncol(docs)
  rdocs = matrix(0,N,M)

  centers <<- generate_centers(docs,k)

  alpha = rnorm(n=1, mean=3, sd=1)
  if(alpha<0)
    alpha=alpha*-1

  for(i in 1:N){

    c = maisProximo(docs[i,],centers)
    rdocs[i,]=docs[i,]+(alpha*(docs[i,]-centers[c,]))

  }

  return(rdocs)

}

tadensidade <- function(docs,k){

  N = nrow(docs)
  M = ncol(docs)
  rdocs = matrix(0,N,M)

  centers <<- generate_centers(docs,k) # kmeans++

  alpha = rnorm(n=1, mean=3, sd=1)
  if(alpha<0)
    alpha=alpha*-1

  for(i in 1:N){

    c = maisProximo(docs[i,],centers)
    rdocs[i,]=docs[i,]-(alpha*(docs[i,]-centers[c,]))
    #rdocs[i,]=alpha*(docs[i,]+docs[c,])

  }

  return(rdocs)

}





generate_centers <- function(x, k, iter.max = 10, nstart = 3, ...) {
  n <- nrow(x) # number of data points
  centers <- numeric(k) # IDs of centers
  distances <- matrix(numeric(n * (k - 1)), ncol = k - 1) # distances[i, j]: The distance between x[i,] and x[centers[j],]
  res.best <- list(tot.withinss = Inf) # the best result among <nstart> iterations
  for (rep in 1:nstart) {
    pr <- rep(1, n) # probability for sampling centers
    for (i in 1:(k - 1)) {
      centers[i] <- sample.int(n, 1, prob = pr) # Pick up the ith center
      distances[, i] <- colSums((t(x) - x[centers[i], ])^2) # Compute (the square of) distances to the center
      pr <- distances[cbind(1:n, max.col(-distances[, 1:i, drop = FALSE]))] # Compute probaiblity for the next sampling
    }
    centers[k] <- sample.int(n, 1, prob = pr)
    ## Perform k-means with the obtained centers
    res <- kmeans(x, x[centers, ], iter.max = iter.max, nstart = 1, ...)
    res$inicial.centers <- x[centers, ]
    ## Store the best result
    if (res$tot.withinss < res.best$tot.withinss) {
      res.best <- res
    }
  }
  res.best$centers
}
