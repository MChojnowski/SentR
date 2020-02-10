STSVAR <- function(data,ph,pl,prob){
  
  dataH <- data * prob
  dataL <- data*(1-prob)
  
  varH <- vars::VAR(data,p=ph)
  varL <- vars::VAR(data,p=pl)
  
  dane_all_org <- dane_all <- varH$y
  
  if(ph>pl){
    dane_all_org <- dane_all <- varH$y
  }
  
  for (i in 1:ph){
    dane_all<-cbind(dane_all,lag(dataH,-i))
  }
  for (i in 1:pl){
    dane_all<-cbind(dane_all,lag(dataL,-i))
  }
  
  varh<-1:(varH$K*varH$p)+ncol(dane_all_org)
  varl<-1:(varL$K*varL$p)+max(varh)
  
  r1<-lm(dane_all[,1]~dane_all[,c(
    if(colnames(dane_all_org)[1] %in% colnames(varH$y)){varh}else{NULL}
    ,if(colnames(dane_all_org)[1] %in% colnames(varL$y)){varl}else{NULL})
    ])
  r2<-lm(dane_all[,2]~dane_all[,c(
    if(colnames(dane_all_org)[2] %in% colnames(varH$y)){varh}else{NULL}
    ,if(colnames(dane_all_org)[2] %in% colnames(varL$y)){varl}else{NULL})
    ])
  r3<-lm(dane_all[,3]~dane_all[,c(
    if(colnames(dane_all_org)[3] %in% colnames(varH$y)){varh}else{NULL}
    ,if(colnames(dane_all_org)[3] %in% colnames(varL$y)){varl}else{NULL})
    ])
  r4<-lm(dane_all[,4]~dane_all[,c(
    if(colnames(dane_all_org)[4] %in% colnames(varH$y)){varh}else{NULL}
    ,if(colnames(dane_all_org)[4] %in% colnames(varL$y)){varl}else{NULL})
    ])
  
  r1$coef<-c(r1$coef[1]
             ,if(!(colnames(dane_all_org)[1] %in% colnames(varH$y))){rep(NA,length(varh))}else{NULL}
             ,r1$coef[-1]
             ,if(!(colnames(dane_all_org)[1] %in% colnames(varL$y))){rep(NA,length(varl))}else{NULL}
  )[1:length(c(varh,varl,1))]
  
  r2$coef<-c(r2$coef[1]
             ,if(!(colnames(dane_all_org)[2] %in% colnames(varH$y))){rep(NA,length(varh))}else{NULL}
             ,r2$coef[-1]
             ,if(!(colnames(dane_all_org)[2] %in% colnames(varL$y))){rep(NA,length(varl))}else{NULL}
  )[1:length(c(varh,varl,1))]
  
  r3$coef<-c(r3$coef[1]
             ,if(!(colnames(dane_all_org)[3] %in% colnames(varH$y))){rep(NA,length(varh))}else{NULL}
             ,r3$coef[-1]
             ,if(!(colnames(dane_all_org)[3] %in% colnames(varL$y))){rep(NA,length(varl))}else{NULL}
  )[1:length(c(varh,varl,1))]
  
  r4$coef<-c(r4$coef[1]
             ,if(!(colnames(dane_all_org)[4] %in% colnames(varH$y))){rep(NA,length(varh))}else{NULL}
             ,r4$coef[-1]
             ,if(!(colnames(dane_all_org)[4] %in% colnames(varL$y))){rep(NA,length(varl))}else{NULL}
  )[1:length(c(varh,varl,1))]
  
  comat1<-rbind(r1$coef,r2$coef,r3$coef,r4$coef)[,c(2:(varH$K*varH$p+1))]
  comat2<-rbind(r1$coef,r2$coef,r3$coef,r4$coef)[,c((varH$K*varH$p+2):length(r1$coef))]
  
  for (i in colnames(varH$y)){
    nazwy<-names(varH$varresult[[i]]$coefficients)
    varH$varresult[[i]]$coefficients<-c(comat1[which(colnames(dane_all_org)==i),],0)  
    names(varH$varresult[[i]]$coefficients)<-nazwy
  }
  
  for (i in colnames(varL$y)){  
    nazwy<-names(varL$varresult[[i]]$coefficients)
    varL$varresult[[i]]$coefficients<-c(comat2[which(colnames(dane_all_org)==i),],0)
    names(varL$varresult[[i]]$coefficients)<-nazwy
  }
  
  ### DODAJ CONSTANT
  dane_h<-varH$datamat[,-c(1:varH$K,ncol(varH$datamat))]
  dane_l<-varL$datamat[,-c(1:varL$K,ncol(varL$datamat))]
  
  #przeliczamy reszty
  for (i in 1:varH$K){
    varH$varresult[[i]]$fitted.values<-as.matrix(dane_h)%*%as.matrix(head(varH$varresult[[i]]$coefficients,-1))
    varH$varresult[[i]]$residuals<-varH$datamat[,i]-varH$varresult[[i]]$fitted.values
  }
  for (i in 1:varL$K){
    varL$varresult[[i]]$fitted.values<-as.matrix(dane_l)%*%as.matrix(head(varL$varresult[[i]]$coefficients,-1))
    varL$varresult[[i]]$residuals<-varL$datamat[,i]-varL$varresult[[i]]$fitted.values
  }
  
  #Obliczamy stale
  c1<-residuals(varH)
  const1<-apply(c1,2,mean)
  #const1
  
  c2<-residuals(varL)
  const2<-apply(c2,2,mean)
  #const2
  
  
  ### podmieniamy sta?e
  for (i in 1:varH$K){
    varH$varresult[[i]]$coefficients[varH$K*varH$p+1]<-const1[i]
  }
  for (i in 1:varL$K){
    varL$varresult[[i]]$coefficients[varL$K*varL$p+1]<-const2[i]
  }
  
  #przeliczamy reszty uwzgledniajac stala
  for (i in 1:varH$K){
    varH$varresult[[i]]$fitted.values<-varH$varresult[[i]]$fitted.values+tail(varH$varresult[[i]]$coefficients,1)
    varH$varresult[[i]]$residuals<-varH$varresult[[i]]$residuals-tail(varH$varresult[[i]]$coefficients,1)
  }
  for (i in 1:varL$K){
    varL$varresult[[i]]$fitted.values<-varL$varresult[[i]]$fitted.values+tail(varL$varresult[[i]]$coefficients,1)
    varL$varresult[[i]]$residuals<-varL$varresult[[i]]$residuals-tail(varL$varresult[[i]]$coefficients,1)
  }
  
  return(list(varH=varH,varL=varL))
}