#' @export
irf_closed<-function(model_svar,impulse,response,closed_channels,Kmax=24){
	
	if(class(model_svar)!="svarest"){
		print("model in not from class 'svarest'")
		break
	}
	
	model <- model_svar$var
	Nimpulse <- colnames(model$y)[impulse]
	Nresponse <- colnames(model$y)[response]

	IRF <- IRF.closed <-rep(NaN,Kmax)

	F <- matrix(unlist(Acoef(model)),nrow=(model$K))
	F <- rbind(F
          	,cbind(diag(model$p-1) %x% diag((model$K)),matrix(0,nrow=(model$K)*(model$p-1),ncol=(model$K)))
  	)


	F.closed <- F
	F.closed[response,seq(0,(model$K)*(model$p-1),by=model$K)+closed_channels] <- 0

	A <- model_svar$A
	A.closed <- model_svar$A
	A.closed[response,closed_channels] <- 0

	temp <- diag(1,(model$K))
	temp <- solve(A)%*%model_svar$B
	temp <- rbind(cbind(temp,matrix(0,nrow=(model$K),ncol=(model$K)*(model$p-1)))
		,matrix(0,nrow=(model$K)*(model$p-1),ncol=(model$K)*(model$p))
	)

	temp.closed <- temp

	for(k in 1:Kmax){
 		IRF[k] <- temp[response,impulse]
		IRF.closed[k] <- temp.closed[response,impulse]
		temp <- F %*% temp
 		temp.closed <- F.closed%*%temp.closed
	}

	return(list(
		model=irf(model,impulse=Nimpulse,response=Nresponse,n.ahead=Kmax)$irf[[1]]
		,calc=IRF
		,Psi=Psi(model,nstep=(Kmax-1))[impulse,response,]
		,closed=IRF.closed
	))
}

	
