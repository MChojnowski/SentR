
gridSearch <- function(data,ph,pl,sent,mean.sent=mean(sent),sd.sent=sd(sent),min.observ=60,grid.length=50,prob.lim=c(0,1),stable_param=1){

  #Models  
  var_dflt <- vars::VAR(data,p=max(ph,pl))
  var_dflt_h <- vars::VAR(data,p=ph)
  var_dflt_l <- vars::VAR(data,p=pl)
  
  # Parameters
  gridM <- qnorm(seq(0.01,0.99,length.out=grid.length),mean.sent,sd.sent)
  gridS <- seq(0.05,3,length.out=grid.length)

  
  if(min.observ<=1){
    min.observ <- min.observ*nrow(var_dflt_h$y)
  }
  
  # Objects

  
  mat.SSE <- NULL
  mat.AIC <- NULL
  mat.ROOTS <- NULL
  mat.IRF <- list(NULL,NULL)
  
  pb <- txtProgressBar(min = 0, max = length(gridM)*length(gridS), style = 3)
  
  for(FctM in gridM){
    
    sum.SSE <- NULL
    sum.AIC <- NULL
    sum.ROOTS <- NULL
    sum.IRF1 <- NULL
    sum.IRF2 <- NULL
    
    for(FctS in gridS){
      
      setTxtProgressBar(pb, ((which(FctM==gridM)-1)*length(gridS)+which(FctS==gridS)))
      
      prob<-pnorm(sent,FctM,FctS)
      
      if(length(prob[round(prob,2)>prob.lim[1] & round(prob,2)<prob.lim[2]])>min.observ){

        dane_monet1<-var_dflt_h$y*prob
        dane_monet2<-var_dflt_l$y*(1-prob)
        
        colnames(dane_monet1)<-colnames(var_dflt_h$y)
        colnames(dane_monet2)<-colnames(var_dflt_l$y)
        
        dane_all_org<-dane_all<-var_dflt_l$y
        
        if(var_dflt_h$K>var_dflt_l$K){
          dane_all_org<-dane_all<-var_dflt_h$y
        }
        
        for (i in 1:var_dflt_h$p){
          dane_all<-cbind(dane_all,lag(dane_monet1,-i))
        }
        for (i in 1:var_dflt_l$p){
          dane_all<-cbind(dane_all,lag(dane_monet2,-i))
        }
        
        varh<-1:(var_dflt_h$K*var_dflt_h$p)+ncol(dane_all_org)
        varl<-1:(var_dflt_l$K*var_dflt_l$p)+max(varh)
        
        r1<-lm(dane_all[,1]~dane_all[,c(
          if(colnames(dane_all_org)[1] %in% colnames(var_dflt_h$y)){varh}else{NULL}
          ,if(colnames(dane_all_org)[1] %in% colnames(var_dflt_l$y)){varl}else{NULL})
          ])
        r2<-lm(dane_all[,2]~dane_all[,c(
          if(colnames(dane_all_org)[2] %in% colnames(var_dflt_h$y)){varh}else{NULL}
          ,if(colnames(dane_all_org)[2] %in% colnames(var_dflt_l$y)){varl}else{NULL})
          ])
        r3<-lm(dane_all[,3]~dane_all[,c(
          if(colnames(dane_all_org)[3] %in% colnames(var_dflt_h$y)){varh}else{NULL}
          ,if(colnames(dane_all_org)[3] %in% colnames(var_dflt_l$y)){varl}else{NULL})
          ])
        r4<-lm(dane_all[,4]~dane_all[,c(
          if(colnames(dane_all_org)[4] %in% colnames(var_dflt_h$y)){varh}else{NULL}
          ,if(colnames(dane_all_org)[4] %in% colnames(var_dflt_l$y)){varl}else{NULL})
          ])
        
        r1$coef<-c(r1$coef[1]
                   ,if(!(colnames(dane_all_org)[1] %in% colnames(var_dflt_h$y))){rep(NA,length(varh))}else{NULL}
                   ,r1$coef[-1]
                   ,if(!(colnames(dane_all_org)[1] %in% colnames(var_dflt_l$y))){rep(NA,length(varl))}else{NULL}
        )[1:length(c(varh,varl,1))]
        
        r2$coef<-c(r2$coef[1]
                   ,if(!(colnames(dane_all_org)[2] %in% colnames(var_dflt_h$y))){rep(NA,length(varh))}else{NULL}
                   ,r2$coef[-1]
                   ,if(!(colnames(dane_all_org)[2] %in% colnames(var_dflt_l$y))){rep(NA,length(varl))}else{NULL}
        )[1:length(c(varh,varl,1))]
        
        r3$coef<-c(r3$coef[1]
                   ,if(!(colnames(dane_all_org)[3] %in% colnames(var_dflt_h$y))){rep(NA,length(varh))}else{NULL}
                   ,r3$coef[-1]
                   ,if(!(colnames(dane_all_org)[3] %in% colnames(var_dflt_l$y))){rep(NA,length(varl))}else{NULL}
        )[1:length(c(varh,varl,1))]
        
        r4$coef<-c(r4$coef[1]
                   ,if(!(colnames(dane_all_org)[4] %in% colnames(var_dflt_h$y))){rep(NA,length(varh))}else{NULL}
                   ,r4$coef[-1]
                   ,if(!(colnames(dane_all_org)[4] %in% colnames(var_dflt_l$y))){rep(NA,length(varl))}else{NULL}
        )[1:length(c(varh,varl,1))]
        
        comat1<-rbind(r1$coef,r2$coef,r3$coef,r4$coef)[,c(2:(var_dflt_h$K*var_dflt_h$p+1))]
        comat2<-rbind(r1$coef,r2$coef,r3$coef,r4$coef)[,c((var_dflt_h$K*var_dflt_h$p+2):length(r1$coef))]
        
        for (i in colnames(var_dflt_h$y)){
          nazwy<-names(var_dflt_h$varresult[[i]]$coefficients)
          var_dflt_h$varresult[[i]]$coefficients<-c(comat1[which(colnames(dane_all_org)==i),],0)  
          names(var_dflt_h$varresult[[i]]$coefficients)<-nazwy
        }
        
        for (i in colnames(var_dflt_l$y)){  
          nazwy<-names(var_dflt_l$varresult[[i]]$coefficients)
          var_dflt_l$varresult[[i]]$coefficients<-c(comat2[which(colnames(dane_all_org)==i),],0)
          names(var_dflt_l$varresult[[i]]$coefficients)<-nazwy
        }
        
        ### DODAJ CONSTANT
        dane_h<-var_dflt_h$datamat[,-c(1:var_dflt_h$K,ncol(var_dflt_h$datamat))]
        dane_l<-var_dflt_l$datamat[,-c(1:var_dflt_l$K,ncol(var_dflt_l$datamat))]
        
        #przeliczamy reszty
        for (i in 1:var_dflt_h$K){
          var_dflt_h$varresult[[i]]$fitted.values<-as.matrix(dane_h)%*%as.matrix(head(var_dflt_h$varresult[[i]]$coefficients,-1))
          var_dflt_h$varresult[[i]]$residuals<-var_dflt_h$datamat[,i]-var_dflt_h$varresult[[i]]$fitted.values
        }
        for (i in 1:var_dflt_l$K){
          var_dflt_l$varresult[[i]]$fitted.values<-as.matrix(dane_l)%*%as.matrix(head(var_dflt_l$varresult[[i]]$coefficients,-1))
          var_dflt_l$varresult[[i]]$residuals<-var_dflt_l$datamat[,i]-var_dflt_l$varresult[[i]]$fitted.values
        }
        
        #Obliczamy stale
        c1<-residuals(var_dflt_h)
        const1<-apply(c1,2,mean)
        #const1
        
        c2<-residuals(var_dflt_l)
        const2<-apply(c2,2,mean)
        #const2
        
        
        ### podmieniamy sta?e
        for (i in 1:var_dflt_h$K){
          var_dflt_h$varresult[[i]]$coefficients[var_dflt_h$K*var_dflt_h$p+1]<-const1[i]
        }
        for (i in 1:var_dflt_l$K){
          var_dflt_l$varresult[[i]]$coefficients[var_dflt_l$K*var_dflt_l$p+1]<-const2[i]
        }
        
        #przeliczamy reszty uwzgledniajac stala
        for (i in 1:var_dflt_h$K){
          var_dflt_h$varresult[[i]]$fitted.values<-var_dflt_h$varresult[[i]]$fitted.values+tail(var_dflt_h$varresult[[i]]$coefficients,1)
          var_dflt_h$varresult[[i]]$residuals<-var_dflt_h$varresult[[i]]$residuals-tail(var_dflt_h$varresult[[i]]$coefficients,1)
        }
        for (i in 1:var_dflt_l$K){
          var_dflt_l$varresult[[i]]$fitted.values<-var_dflt_l$varresult[[i]]$fitted.values+tail(var_dflt_l$varresult[[i]]$coefficients,1)
          var_dflt_l$varresult[[i]]$residuals<-var_dflt_l$varresult[[i]]$residuals-tail(var_dflt_l$varresult[[i]]$coefficients,1)
        }
        
        checkpoint.reszta<-list(apply(residuals(var_dflt_h),2,mean),apply(residuals(var_dflt_l),2,mean))
        
        if(var_dflt_h$K == var_dflt_l$K){
          VarFitted<-prob*ts(fitted(var_dflt_h),freq=12,end=end(var_dflt_h$y)) +
            (1-prob)*ts(fitted(var_dflt_l),freq=12,end=end(var_dflt_l$y))
        } else {
          VarFitted<-prob*ts(cbind(fitted(var_dflt_h),as.matrix(array(0,dim=c(nrow(fitted(var_dflt_h)),sum(!colnames(dane_all_org) %in% colnames(fitted(var_dflt_h))))))),freq=12,end=end(var_dflt_h$y)) +
            (1-prob)*ts(cbind(fitted(var_dflt_l),as.matrix(array(0,dim=c(nrow(fitted(var_dflt_h)),sum(!colnames(dane_all_org) %in% colnames(fitted(var_dflt_l))))))),freq=12,end=end(var_dflt_l$y))
        }
        
        colnames(VarFitted)<-colnames(var_dflt$y)
        
        resTempVar<-(VarFitted-var_dflt$y)
        
        resTempVar[,1]<-resTempVar[,1]/(apply(dane_monet,2,sd)[1])
        resTempVar[,2]<-resTempVar[,2]/(apply(dane_monet,2,sd)[2])
        resTempVar[,3]<-resTempVar[,3]/(apply(dane_monet,2,sd)[3])
        resTempVar[,4]<-resTempVar[,4]/(apply(dane_monet,2,sd)[4])
        
        sum.SSE<-c(sum.SSE,sqrt(sum(resTempVar^2)))
        sum.ROOTS<-c(sum.ROOTS,max(roots(var_dflt_h),roots(var_dflt_l)))
        sum.AIC<-rbind(sum.AIC,c(AIC(var_dflt_h),AIC(var_dflt_l)))
        sum.IRF1<-c(sum.IRF1,-1)
        sum.IRF2<-c(sum.IRF2,-1)
        
      }else{
        
        sum.SSE<-c(sum.SSE,Inf)
        sum.ROOTS<-c(sum.ROOTS,Inf)
        sum.AIC<-rbind(sum.AIC,c(Inf,Inf))
        sum.IRF1<-rbind(sum.IRF1,c(Inf,Inf))
        sum.IRF2<-rbind(sum.IRF2,c(Inf,Inf))
      }
      
      
    }
    
    mat.SSE<-cbind(mat.SSE,sum.SSE)
    mat.ROOTS<-cbind(mat.ROOTS,sum.ROOTS)
    mat.AIC<-cbind(mat.AIC,apply(sum.AIC,1,sum))
    mat.IRF[[1]]<-cbind(mat.IRF[[1]],sum.IRF1)
    mat.IRF[[2]]<-cbind(mat.IRF[[2]],sum.IRF1)
    
  }
  
  cord_sse <- mat.SSE
  cord_sse[mat.ROOTS >= stable_param] <- Inf
  cord_sse <- which(cord_sse == min(cord_sse), arr.ind = TRUE)
  
  return(list(gridM=gridM
              ,gridS=gridS
              ,sse=mat.SSE
              ,roots=mat.ROOTS
              ,aic=mat.AIC
              ,irf=mat.IRF
              ,sent=sent
              ,optparam=c(gridS[cord_sse[1]],gridM[cord_sse[2]])
              ,stable_param=stable_param
  ))
}
