
##### IRF 3D #####
irf3D <- function(model,impulse,response,Amat,Bmat,h=12,n=100){
  
  Amat <<- Amat
  Bmat <<- Bmat
  
  RESULT <- list()
  
  for(p in seq(0,1,by=(1/n))){
    
    model_var_irf <- model$varH
    if(model$varL$p > model$varH$p){
      model_var_irf <- model$varL
    }
    
    for (i in 1:4){
      vc1 <- c(head(model$varH$varresult[[i]]$coefficients,-1),rep(0,4*(max(model$varH$p,model$varL$p) - model$varH$p)))
      vc2 <- c(head(model$varL$varresult[[i]]$coefficients,-1),rep(0,4*(max(model$varH$p,model$varL$p) - model$varL$p)))
      const <- p * tail(model$varH$varresult[[i]]$coefficients,1) + 
        (1-p) * tail(model$varL$varresult[[i]]$coefficients,1)
      
      nazwy <- names(model_var_irf$varresult[[i]]$coefficients)
      model_var_irf$varresult[[i]]$coefficients <- c(p*vc1+(1-p)*vc2,const)
      names(model_var_irf$varresult[[i]]$coefficients) <- nazwy
      
      model_var_irf$varresult[[i]]$fitted.values<-as.matrix(model_var_irf$datamat[,-c(1:4)])%*%as.matrix(model_var_irf$varresult[[i]]$coefficients)
      model_var_irf$varresult[[i]]$residuals<-model_var_irf$datamat[,i]-model_var_irf$varresult[[i]]$fitted.values
    }
    
    model_svar_irf <- SVAR(model_var_irf,Amat=Amat,Bmat=Bmat,lrtest = FALSE, max.iter=1000)
    irf_plot <- irf(model_svar_irf, impulse = impulse, response=response, n.ahead = h, ortho=TRUE, boot=FALSE)
    
    for(im in impulse){
      for(re in response){
        RESULT[[im]][[re]] <- cbind(RESULT[[im]][[re]],irf_plot$irf[[im]][,re])
      }
    }
  }
  return(RESULT)
}


##### Plot IRF 3D #####
plotIRF3D <- function(IRF,impulse,response,plot.levels=48,smooth.param=50){
  
  IRF_BASE <- IRF[[impulse]][[response]]
  IRF_PLOT <- t(apply(IRF_BASE,2,approxIRF,smooth.param=smooth.param))
  colnames(IRF_PLOT) <- as.character((1:ncol(IRF_PLOT)-1)/smooth.param)
    
    IRF_PLOT <- as.tibble(IRF_PLOT) %>%
                add_column(p = seq(0,1,by=(1/(ncol(IRF_BASE)-1))),.before=1) %>%
                gather(key = "Horizon",value="Value",-1) %>%
                mutate(Horizon=as.numeric(Horizon))
  
  
    ggplot(IRF_PLOT,aes(x = Horizon, y = p)) +
      geom_raster(aes(fill = Value)) +
      scale_fill_gradient2(low="red",high="green",mid="white",midpoint=0) 
}

approxIRF <- function(x,smooth.param=1){
  return(approx(0:(length(x)-1),x,n = length(x) * smooth.param)$y)
}