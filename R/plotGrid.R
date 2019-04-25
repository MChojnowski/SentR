plotGrid <- function(grid){
  grid_sse <- grid$sse
  colnames(grid_sse) <- 1:ncol(grid_sse)
  
  grid_sse <- as.tibble(grid_sse) %>%
    add_column(Mean=1:nrow(.),.before=1) %>%
    gather(key="Sd",value="SSE",-1) %>%
    mutate(Sd=as.numeric(Sd))
  
  grid_roots <- grid$roots
  colnames(grid_roots) <- 1:ncol(grid_roots)
  
  grid_roots <- as.tibble(grid_roots) %>%
    add_column(Mean=1:nrow(.),.before=1) %>%
    gather(key="Sd",value="MinRoot",-1) %>%
    mutate(Sd=as.numeric(Sd)) 
  
  plot_data <- left_join(grid_sse,grid_roots,by=c("Mean","Sd"))
  
  cord <- plot_data %>%
    filter(MinRoot < 1) %>%
    filter(SSE==min(SSE))
  
  ggplot(plot_data,aes(x = Mean, y = Sd)) +
    geom_raster(aes(fill = SSE)) +
    scale_fill_gradientn(colors = matlab.like2(20),na.value="white") +
    geom_contour(mapping = aes(z = SSE),col="grey",bins=40) +
    geom_contour(mapping = aes(z = MinRoot),col="white",breaks=1,lty=4,lwd=1) +
    geom_point(data = cord,aes(size=3,color="red"),pch= 12)
  }

