imputation.kh<-function(data, vari.matrix, method='k', cluster=6, nstart=20){
  if(method=='k'){
    clusters <- kmeans(vari.matrix, cluster, nstart)
    clusterCut <- clusters$cluster
  }else if(method=='h'){
    clusters <- hclust(dist(vari.matrix))
    clusterCut <- cutree(clusters, cluster)
  }
  clu.matrix<-matrix(rep(clusterCut,dim(data)[2]),ncol =(dim(data)[2]))
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

  result<-list(data=data, clu.matrix=clu.matrix, sd.1.j=sd.1.j, sd.1.k=sd.1.k)
}
