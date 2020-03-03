getSents <- function(data,lagV,lagB,rep=32,Amat=diag(1,ncol(data),ncol(data)),Bmat=matrix(as.numeric(diag(NA,ncol(data),ncol(data))),ncol(data),ncol(data)),exodata=NULL,bvec=c(1,1,1,1,1,1,1),qm=4){ 
  
  #Sense check
  if (sum(class(data) %in% "ts")==0){
    cat("Dataset is not a time series (ts). Please provide proper dataset \n")
    break
  }
  
  # Models
  model_var <- vars::VAR(data,p=lagV)
  model_svar <- vars::SVAR(model_var,max.iter=1000,Amat=Amat,Bmat=Bmat,estmethod="direct")
  
  BVAR_sent <- list()
  
  for(i in 1:rep){
    model_bsvar <- szbsvar(data,lagB,exodata,bvec[1],bvec[2],bvec[3],bvec[4],bvec[5],bvec[6],bvec[7],model_svar$A,qm=qm)
    signal <- t(solve(model_bsvar$A0.mode) %*% t(model_bsvar$structural.innovations))
    
    x2<-NULL
    for(ix in 1:(ncol(data)/2)){
      x2<-cbind(x2,
                signal[,(ncol(data)/2)+ix]-apply(as.matrix(signal[,c(1:(ncol(data)/2))[-ix]]),1,sum)
      )
    }
    
    for(j in 2:(ncol(data)/2)){
      if(cor(x2[,1],x2[,j])<0){
        x2[,j] <- (-1)*x2[,j]
      }
    }
    
    BVAR_sent[[i]] <- apply(x2,1,sum)/(ncol(data)/2)
    
    if(i>1){
      if(cor(BVAR_sent[[1]],BVAR_sent[[i]])<0){BVAR_sent[[i]]<-(-1)*BVAR_sent[[i]]}
    }
    
  }
  
  sent <- apply(matrix(unlist(BVAR_sent),ncol=rep),1,median)
 
  return(list(
    sent=ts(sent,end=end(data))
  ))
}
