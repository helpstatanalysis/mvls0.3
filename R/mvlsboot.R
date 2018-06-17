#' mvlsboot
#'
#' Takes a longitudinal dataset and impute missing value with a machine learning-based method many time.
#'
#' @author Lorenzo Querci <lorenzo.querci@studio.unibo.it>
#'
#' @param data A dataset (more than two longitudinal mesurements)
#' @param d It is the percentage of change between two-sided mesurements to consider it bigger, smaller or the same, It useful to built the var.matrix
#' @param method It represent the type of machine learning algorithm. 'k' for k-mean and 'h' for hierical
#' @param cluster It's the number of cluster. Default setting it's 6. It depends on number of longitudinal mesurements. It could be use mvls.print to decide best cluster number.
#' @param nstart It is the nstart setting of function k-mean. Defualt it'20. Not requested for 'h' method.
#' @param pre.imp TRUE/FALSE (default F). It permit to pre-impute data to built the vari.matrix, It could be reduce cluster with only missing value.
#' @param imp.method It's the type of pre-imputation. Defaul it's 'mean', but there is also 'locf' possibility
#' @param boot Impute many times to reduce error from cluster imputation. 'low' for 5 imputation, 'medium' for 10 imputation and 'high' for 15 imputation.
#'
#' @return **data** it's the data-set with imputation ans **sd.2** It contain the sd for each data imputed at multiple imputation method. Different from sd.1
#'
#' @import "zoo"
#' @import "mice"
#'
#' @export

mvlsboot<-function(data, d=0.1, method='k', cluster=6, nstart=20, pre.imp=F, imp.method='mean', boot='medium'){
  sd.2<-matrix(c(rep(NA,(dim(data)[2]*dim(data)[1]))),ncol=dim(data)[2],nrow=dim(data)[1])
  result.boot<-data
  if(boot=='low'){
    db.boot.1<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.2<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.3<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.4<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.5<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    for (i in 1:dim(result.boot)[2]){
      for (l in 1:dim(result.boot)[1]) {
        sd.2[l,i]<-rnorm(1,mean=0,sd=sd(c(db.boot.1[l,i],db.boot.2[l,i],db.boot.3[l,i],db.boot.4[l,i],db.boot.5[l,i]),na.rm = T))
        result.boot[l,i]<-mean(c(db.boot.1[l,i],db.boot.2[l,i],db.boot.3[l,i],db.boot.4[l,i],db.boot.5[l,i]),na.rm = T)+sd.2[l,i]
      }
    }
  }
  if(boot=='medium'){
    db.boot.1<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.2<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.3<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.4<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.5<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.6<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.7<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.8<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.9<-mvls(data, d, method, cluster, nstart)$data
    db.boot.10<-mvls(data, d, method, cluster, nstart)$data
    for (i in 1:dim(result.boot)[2]){
      for (l in 1:dim(result.boot)[1]) {
        sd.2[l,i]<-rnorm(1,mean=0, sd=sd(c(db.boot.1[l,i],db.boot.2[l,i],db.boot.3[l,i],db.boot.4[l,i],db.boot.5[l,i],db.boot.6[l,i],db.boot.7[l,i],db.boot.8[l,i],db.boot.9[l,i],db.boot.10[l,i]),na.rm = T))
        result.boot[l,i]<-mean(c(db.boot.1[l,i],db.boot.2[l,i],db.boot.3[l,i],db.boot.4[l,i],db.boot.5[l,i],db.boot.6[l,i],db.boot.7[l,i],db.boot.8[l,i],db.boot.9[l,i],db.boot.10[l,i]),na.rm = T)+sd.2[l,i]
      }
    }
  }
  if(boot=='high'){
    db.boot.1<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.2<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.3<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.4<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.5<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.6<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.7<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.8<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.9<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.10<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.11<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.12<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.13<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.14<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    db.boot.15<-mvls(data, d, method, cluster, nstart,pre.imp, imp.method)$data
    for (i in 1:dim(result.boot)[2]){
      for (l in 1:dim(result.boot)[1]) {
        sd.2[l,i]<-rnorm(1,mean=0, sd=sd(c(db.boot.1[l,i],db.boot.2[l,i],db.boot.3[l,i],db.boot.4[l,i],db.boot.5[l,i],db.boot.6[l,i],db.boot.7[l,i],db.boot.8[l,i],db.boot.9[l,i],db.boot.10[l,i],db.boot.11[l,i],db.boot.12[l,i],db.boot.13[l,i],db.boot.14[l,i],db.boot.15[l,i]),na.rm = T))
        result.boot[l,i]<-mean(c(db.boot.1[l,i],db.boot.2[l,i],db.boot.3[l,i],db.boot.4[l,i],db.boot.5[l,i],db.boot.6[l,i],db.boot.7[l,i],db.boot.8[l,i],db.boot.9[l,i],db.boot.10[l,i],db.boot.11[l,i],db.boot.12[l,i],db.boot.13[l,i],db.boot.14[l,i],db.boot.15[l,i]),na.rm = T)+sd.2[l,i]
      }
    }
  }
  return(list(data=result.boot,sd.2=sd.2))
}
