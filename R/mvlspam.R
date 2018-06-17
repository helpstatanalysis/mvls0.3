#' mvls.pam
#'
#' Takes a longitudinal dataset and impute missing value with pam (k by pamk) machine learning-based method.
#'
#' @author Lorenzo Querci <lorenzo.querci@studio.unibo.it>
#'
#' @param data A dataset (more than two longitudinal mesurements).
#' @param imp.method Type od pre-imputation. pam and pamk don't use missing value. Possible "mean" or "locf".
#' @param krange Range for k in pamk (fpc package) function.
#' @param scaling TRUE/FALSE value (default T). Scale data before clustering.
#'
#' @return **data** it's the data-set with imputation, **pam object**, **sd.1** contains the sd for each data imputed at single imputation method. Different from sd.2.
#'
#' @import "fpc"
#' @import "cluster"
#' @import "mice"
#' @import "zoo"
#'
#' @export

mvls.pam<-function(data, imp.method='mean', krange=2:12, scaling=T){
  data<-exclude(data)
  data.pi<-preimputation(data,imp.method)
  k<-pamk(data.pi,krange, scaling)$nc
  pam<-pam(data.pi, k)
  pam$clustering
  clusterCut<-pam$clustering
  clu.matrix<-matrix(rep(clusterCut,dim(data)[2]),ncol=(dim(data)[2]))
  sd.1.j<-matrix(c(rep(NA,(dim(data)[2]*dim(data)[1]))),ncol=dim(data)[2],nrow=dim(data)[1])
  sd.1.k<-matrix(c(rep(NA,(dim(data)[2]*dim(data)[1]))),ncol=dim(data)[2],nrow=dim(data)[1])
  for(i in 1:dim(data)[2]){
    for(l in 1:dim(data)[1]){
      if(is.na(data[l,i])==T){
        if(i<(dim(data)[2]-1)){
          if(is.na(data[l,(i+1)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i+1)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i+1)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i+1)]+(j-k))>1){h<-1
              }else if((data[l,(i+1)]+(j-k))<0){h<-0
              }else(h<-data[l,(i+1)]+(j-k))
              data[l,i]<-h
            }
          }
          if(is.na(data[l,(i+1)])==T && is.na(data[l,(i+2)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i+2)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i+2)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i+2)]+(j-k))>1){h<-1
              }else if((data[l,(i+2)]+(j-k))<0){h<-0
              }else(h<-data[l,(i+2)]+(j-k))
              data[l,i]<-h
            }
          }
        }
        if(i==(dim(data)[2]-1)){
          if(is.na(data[l,(i-1)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i-1)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i-1)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i-1)]-(j-k))>1){h<-1
              }else if((data[l,(i-1)]-(j-k))<0){h<-0
              }else(h<-data[l,(i-1)]-(j-k))
              data[l,i]<-h
            }
          }
          if(is.na(data[l,(i-1)])==T && is.na(data[l,(i-2)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i-2)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i-2)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i-2)]-(j-k))>1){h<-1
              }else if((data[l,(i-2)]-(j-k))<0){h<-0
              }else(h<-data[l,(i-2)]-(j-k))
              data[l,i]<-h
            }
          }
        }
        if(i==(dim(data)[2])){
          if(is.na(data[l,(i-1)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i-1)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i-1)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i-1)]-(j-k))>1){h<-1
              }else if((data[l,(i-1)]-(j-k))<0){h<-0
              }else(h<-data[l,(i-1)]-(j-k))
              data[l,i]<-h
            }
          }
          if(is.na(data[l,(i-1)])==T && is.na(data[l,(i-2)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i-2)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i-2)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i-2)]-(j-k))>1){h<-1
              }else if((data[l,(i-2)]-(j-k))<0){h<-0
              }else(h<-data[l,(i-2)]-(j-k))
              data[l,i]<-h
            }
          }
        }
      }
    }
  }

  for(i in 1:dim(data)[2]){
    for(l in 1:dim(data)[1]){
      if(is.na(data[l,i])==T){
        if(i<(dim(data)[2]-1)){
          if(is.na(data[l,(i+1)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i+1)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i+1)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i+1)]+(j-k))>1){h<-1
              }else if((data[l,(i+1)]+(j-k))<0){h<-0
              }else(h<-data[l,(i+1)]+(j-k))
              data[l,i]<-h
            }
          }
          if(is.na(data[l,(i+1)])==T && is.na(data[l,(i+2)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i+2)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i+2)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i+2)]+(j-k))>1){h<-1
              }else if((data[l,(i+2)]+(j-k))<0){h<-0
              }else(h<-data[l,(i+2)]+(j-k))
              data[l,i]<-h
            }
          }
        }
        if(i==(dim(data)[2]-1)){
          if(is.na(data[l,(i-1)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i-1)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i-1)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i-1)]-(j-k))>1){h<-1
              }else if((data[l,(i-1)]-(j-k))<0){h<-0
              }else(h<-data[l,(i-1)]-(j-k))
              data[l,i]<-h
            }
          }
          if(is.na(data[l,(i-1)])==T && is.na(data[l,(i-2)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i-2)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i-2)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i-2)]-(j-k))>1){h<-1
              }else if((data[l,(i-2)]-(j-k))<0){h<-0
              }else(h<-data[l,(i-2)]-(j-k))
              data[l,i]<-h
            }
          }
        }
        if(i==(dim(data)[2])){
          if(is.na(data[l,(i-1)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i-1)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i-1)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i-1)]-(j-k))>1){h<-1
              }else if((data[l,(i-1)]-(j-k))<0){h<-0
              }else(h<-data[l,(i-1)]-(j-k))
              data[l,i]<-h
            }
          }
          if(is.na(data[l,(i-1)])==T && is.na(data[l,(i-2)])==F){
            val<-as.vector(as.matrix(which(clusterCut==clu.matrix[l,i]))[,1])
            sd.1.j[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),i])))
            sd.1.k[l,i]<-rnorm(1, mean=0, sd=sd(!is.na(data[c(val),(i-2)])))
            j<-mean(data[c(val),i],na.rm = T)+sd.1.j[l,i]
            k<-mean(data[c(val),(i-2)], na.rm = T)+sd.1.k[l,i]
            if(is.na(j)==F && is.na(k)==F){
              if((data[l,(i-2)]-(j-k))>1){h<-1
              }else if((data[l,(i-2)]-(j-k))<0){h<-0
              }else(h<-data[l,(i-2)]-(j-k))
              data[l,i]<-h
            }
          }
        }
      }
    }
  }

  result<-list(data=data, pam=pam, sd.1=sd.1.j-sd.1.k)
}
