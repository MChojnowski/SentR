plotIRF <- function(model,model_ref=NULL,Amat=NULL,Bmat=NULL,impulse,response=colnames(model$varH),n.ahead=12,ci=0.95,nboot=1000,main=NULL){
  
  Amat <<- Amat
  Bmat <<- Bmat
  ph <<- model$varH$p
  pl <<- model$varL$p
  
  models <- list()
  models[[1]] <- vars::SVAR(model$varH,Amat=Amat,Bmat=Bmat,lrtest = FALSE, max.iter = 1000)
  models[[2]] <- vars::SVAR(model$varL,Amat=Amat,Bmat=Bmat,lrtest = FALSE, max.iter = 1000)
  
  if (!is.null(model_ref)){
    models[[3]] <- vars::SVAR(model_ref,Amat=Amat,Bmat=Bmat,lrtest = FALSE, max.iter = 1000)
  }
  
  for (res in response){
    for (i in 1:length(models)){
      cat(paste("Plotting response of",res,"by",impulse,"\n"))
      plotIRF.aux(models[[i]],impulse,res,n.ahead,ci,nboot,main[i])
      cat("Done\n")
    }
  }
  
}


##### Plot Aux #####

plotIRF.aux <- function(model,imp,res,h,ci,nboot,main){
  
  var_plot <- irf(model
               ,impulse = imp
               ,response = res	
               ,n.ahead = h
               ,ortho = TRUE
               ,boot = TRUE
               ,runs = nboot
               ,ci = ci
  )
  

  plot(x=c(1:(h+1))
       ,y = unlist(var_plot$Lower)
       ,type = "l"	
       ,lwd = 3
       ,lty = 2
       ,col = "red"
       ,ylab = paste("Impulse of ",imp)
       ,xlab = res
       ,ylim = range(c(unlist(var_plot$Lower)
                     ,unlist(var_plot$Upper)))
       ,main = main
  )
  
  lines(x=c(1:(h+1)),y=unlist(var_plot$Upper),type="l",lwd = 3,lty=2,col="red")
  lines(x=c(1:(h+1)),y=unlist(var_plot$irf),type="l", lwd = 3)
  abline(h = 0)
  
}


##### IRF 3D #####
irf3D <- function(){
  
}
##### Plot IRF 3D #####
ploIRF3D <- function(){
  
}

