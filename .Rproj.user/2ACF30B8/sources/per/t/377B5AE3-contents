avaliacao <- function(G, modelo, K) {


  N = nrow(G)
  C = ncol(G)

  ## vetor para guardar a quantidade de erros de cada agrupamento em G referente ao modelo perfeito
  resultado = c(rep(-1,N))

  for(i in 1:N){

    ## usado uma variavel auxiliar para que se calule as equivalencias de G com o conjunto original
    modelo_aux = modelo

    ## Montar as associacoes entre groupIDs
    A = assossiation(G[i,], modelo_aux,K)

    ## Ajuste dos groupIDs do modelo conforme as associacoes
    modelo_aux = ajustarAssociacoes(modelo, A)
    #modelo_aux = ajustarAssociacoes2(modelo, A,unique(G[i,]))

    ## Calculo no numero de diferencas entre os dois agrupamentos
    if (length(modelo_aux) != length(G[i,]))
      stop()

    resultado[i] = sum(modelo_aux != G[i,])
  }

  return(resultado)
}

assossiation <- function(c1, c2, K){

  A = rep(NA,K)
  J = jaccard(c1,c2,K)

  J[is.nan(J)] = 0.0

  for(i in 1:K){
    p = which(J==max(J), arr.ind=TRUE)[1,1] #linha do max
    q = which(J==max(J), arr.ind=TRUE)[1,2] #coluna do max
    A[p]=q
    J[p,]=-Inf
    J[,q]=-Inf
  }

  return(A)
}

jaccard <-function(G_kmeans, G_hier, K){

  mj = matrix(0, K, K)
  for (i in 0:(K))
    for (j in 0:(K)) {
      juncao = c(which(G_kmeans==i), which(G_hier == j))
      interseccao = length(juncao)-length(unique(juncao))
      uniao = length(unique(juncao))
      mj[i,j] = (interseccao)/(uniao)
    }

  return(mj)

}

ajustarAssociacoes <- function(modelo, associacao) {
  modelo_aux = rep(-1, length(modelo))

  for (i in 1:length(modelo)) {
    modelo_aux[i] = which(associacao==modelo[i])
  }
  return(modelo_aux)
}
maisProximo <- function(data, centroides) {

  distancias = calculaDistanciaDocumentoCentroid(data, centroides)
  pos = which.max(distancias)
  return(pos)
}
calculaDistanciaDocumentoCentroid <- function(doc, centroides){
  N = nrow(centroides)
  for (i in 1:N) {
    norma = sqrt(sum(centroides[i,]**2))
    if(norma>0)
      centroides[i,] = centroides[i,]/norma
  }

  dists=rep(0,N)
  N = nrow(centroides)
  for (i in 1:N)
    dists[i] = sum(doc*centroides[i,])

  return(dists)
}

organizavetor <- function(v){

  uv=unique(v)
  n=length(v)
  nv=rep(-1,n)

  for(i in 1:length(uv)){
    nv[v==uv[i]]=i
  }
  return(nv)
}


heat <- function(data){

  library(gplots)
  library(RColorBrewer)

  s<-rep(0,nrow(data))
  for(i in 1:nrow(data))
    s[i]=sum(data[i,2:ncol(data)])
  data<-cbind(data,s)
  data<-data[order(data[,ncol(data)]),]
  data=data[,1:ncol(data)-1]
  mat_data <- data.matrix(data)
  rownames(mat_data) <- rownames(data)

  mat_data <- mat_data[,names(sort(apply(mat_data, 2, sum)))]


  my_palette <- colorRampPalette(c("#ffffd9", "#41b6c4", "#081d58"))(n = 299) #results1
  col_breaks = c(seq(-1,0,length=8), # for red
                 seq(0,0.8,length=7),  # for yellow
                 seq(0.81,20,length=7)) # for green


  heatmap.2(mat_data,
            main = "Heatmap", # heat map title
            notecol="black",      # change font color of cell labels to black
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(7,5),     # widens margins around plot
            col=my_palette,       # use on color palette defined earlier
            Rowv=FALSE,
            key=TRUE,
            keysize=1,
            key.par=list(mar=c(3.5,0,3,0)),
            Colv=FALSE,
            #lwid=c(0.3,4), lhei=c(1,4),
            dendrogram="none"
  )            # turn off column clustering
}

normalizeValues <- function(df){
  val = subset(df,df$measure=="ari")[,4]
  rn = row.names(subset(df,df$measure=="ari"))
  nval = (val-min(val))/(max(val)-min(val))
  df[rn,4]=nval

  val = subset(df,df$measure=="dunn")[,4]
  rn = row.names(subset(df,df$measure=="dunn"))
  nval = (val-min(val))/(max(val)-min(val))
  df[rn,4]=nval

  val = subset(df,df$measure=="slh")[,4]
  rn = row.names(subset(df,df$measure=="slh"))
  nval = (val-min(val))/(max(val)-min(val))
  df[rn,4]=nval

  return(df)
}
