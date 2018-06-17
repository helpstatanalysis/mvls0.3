#' mvls
#'
#' Takes a longitudinal dataset and impute missing value with a machine learning-based method.
#'
#' @author Lorenzo Querci <lorenzo.querci@studio.unibo.it>
#'
#' @param data A dataset (more than two longitudinal mesurements)
#' @param d Percentage of change between two-sided mesurements to consider it bigger, smaller or the same, It useful to built the var.matrix.
#' @param method It represent the type of machine learning algorithm. 'k' for k-mean and 'h' for hierical.
#' @param cluster It's the number of cluster. Default setting it's 6. It depends on number of longitudinal mesurements. It could be use mvls.print to decide best cluster number.
#' @param nstart It is the nstart setting of function k-mean. Defualt it'20. Not requested for 'h' method.
#' @param pre.imp TRUE/FALSE (default F). It permit to pre-impute data to built the vari.matrix, It could be reduce cluster with only missing value.
#' @param imp.method It's the type of pre-imputation. Defaul it's 'mean', but there is also 'locf' possibility.
#'
#' @return **data** it's the data-set with imputation, **cluster** it's the cluster matrix, **matrix** it's the vari.matrix and **sd.1** contains the sd for each data imputed at single imputation method. Different from sd.2
#'
#' @import "zoo"
#' @import "mice"
#'
#' @export

mvls<-function(data, d=0.1, method='k', cluster=6, nstart=20, pre.imp=F, imp.method='mean'){
  if(pre.imp==T){
    results<-toclusterfunc.imp(data,d,imp.method)
  }else if(pre.imp==F){
    results<-toclusterfunc.noimp(data,d)
  }
  data.f<-results$data
  index.f<-results$index
  vari.matrix.f<-results$vari.matrix
  result<-imputation.kh(data.f, vari.matrix.f, method, cluster, nstart)
  db<-result$data*index.f
  sd.1.j<-result$sd.1.j*index.f
  sd.1.k=result$sd.1.k*index.f
  sd.1<-sd.1.j-sd.1.k
  return(list(data=db, vari.matrix=vari.matrix.f, data.norm=result$data, cluster=result$clu.matrix, matrix=results$matrix, sd.1=sd.1))
}
