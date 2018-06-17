preimputation<-function(data, imp.method='mean'){
  if(imp.method=='mean'){
    data<-mice(data, method = 'mean',printFlag = F)
    datapreimput<-complete(data)
  }
  if(imp.method=='locf'){
    data<-t(db.prov)
    data<-na.locf(data)
    data<-t(data)
    data<-mice(data, method = 'mean', printFlag = F)
    datapreimput<-complete(data)
  }
  return(datapreimput)
}
