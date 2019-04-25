plotGrid <- function(grid,sent){
  par(mfrow=c(2,2))
  
  graphics::filled.contour(grid$gridS
                 ,grid$gridM
                 ,as.matrix(grid$sse)
                 ,xlab="Sd"
                 ,ylab="Mean"
                 ,color.palette = matlab.like2
                 ,main="GridSearch"
                 ,nlevels=96
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
  
  graphics::filled.contour(grid$gridS
                           ,grid$gridM
                           ,as.matrix(grid$sse)
                           ,xlab="Sd"
                           ,ylab="Mean"
                           ,color.palette = matlab.like2
                           ,main="GridSearch"
                           ,nlevels=96
  )
  
  plot(pnorm(grid$sent,grid$sent$gridM[cord[2]],grid$sent$gridS[cord[1]]))
}