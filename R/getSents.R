#' getSents function
#' 
#' This function allows you to extract martker sentiments from historical data.
#' @param data  - time series data , ordered ia specific way ($IP_{1}, IP_{2},...,CR_{1},CR_{2}$)
#' @param lagV - number of lags for SVAR model ($A_{0}$ is approximated based on it)
#' @param lagB - number og lags for SZBSVAR model
#' @param rep [optional] - numbre of repetitions of SZSBVAR model. Default is 32
#' @param Amat [optional]- $A_{0}$ matrix for SvAR model. Default is an identity matrix
#' @param Bmat [optional] - $B$ matrix for SVAR model. Default is an identity matrix
#' @param exodata [optional] - matrix with exogeneous data. Defaultis NULL
#' @param bvec [optional] - vector with metaparameters required for Sim-Zha BSVAR model. Default is a vector with 1.
#' @param qm [optional] - scalar representing the frequency of the data. As of now either 4 or 12 are available. Any other number different than 12 is regarding as 4. Default is 4.
#' @author Michal Chojnowski
#' @keywords Bayesian SVAR, sentiments, extraction
#' @example getSents()



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
