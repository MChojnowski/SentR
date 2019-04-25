plotGrid <- function(grid){
  par(mfrow=c(2,2))
  cord<-which(grid$sse == min(grid$sse), arr.ind = TRUE)
  
  graphics::filled.contour(grid$gridS
                 ,grid$gridM
                 ,as.matrix(grid$sse)
                 ,xlab="Sd"
                 ,ylab="Mean"
                 ,color.palette = matlab.like2
                 ,main="GridSearch"
                 ,nlevels=96
                 ,plot.axes={points(grid$gridS[cord[1]],grid$gridM[cord[2]],pch=19,col="orange");axis(1);axis(2)}
  )
  
  graphics::filled.contour(grid$gridS
                           ,grid$gridM
                           ,as.matrix(grid$roots)
                           ,xlab="Sd"
                           ,ylab="Mean"
                           ,color.palette = matlab.like2
                           ,main="GridSearch"
                           ,nlevels=96
  )
  
  grid$sse[grid$roots>=1]<-Inf
  cord<-which(grid$sse == min(grid$sse), arr.ind = TRUE)
  
  graphics::filled.contour(grid$gridS
                           ,grid$gridM
                           ,as.matrix(grid$sse)
                           ,xlab="Sd"
                           ,ylab="Mean"
                           ,color.palette = matlab.like2
                           ,main="GridSearch"
                           ,nlevels=96
                           ,plot.axes={points(grid$gridS[cord[1]],grid$gridM[cord[2]],pch=19,col="orange");axis(1);axis(2)}
  )
  
  plot(pnorm(grid$sent,grid$gridM[cord[2]],grid$gridS[cord[1]]))
}
