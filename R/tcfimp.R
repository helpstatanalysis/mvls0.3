toclusterfunc.imp<-function(data, d, imp.method='mean'){
  data<-exclude(data)
  matrix<-data
  data.imp<-preimputation(data, imp.method)
  result.imp<-normalize(data)
  db.var<-var.matrix(result.imp$data,d)
  result<-normalize(data)
  data<-result$data
  index<-result$index
  result<-return(list(data=data,index=index,vari.matrix=db.var,matrix=matrix))
}
