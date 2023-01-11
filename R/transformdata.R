# Transform data according to the legend


transformData <- function(data, stationary_transform)
{
  v = as.matrix(data)
  vv <- matrix(NA, NROW(v), NCOL(v))
  for(jj in 1:NCOL(vv))
  {
   
    transf <- stationary_transform[jj]
    # no change
    if(transf==1){
      vv[,jj] <- v[,jj]
    }
    
    # First diff
    if(transf==2){
      vv[,jj] <- c(NA, diff(v[,jj]))
    }
    
    # Second diff
    if(transf==3){
      vv[,jj] <- c(NA,NA, diff(diff(v[,jj])))
    }
    
    # p. change
    if(transf==4){
      vv[,jj] <- c(NA, (v[2:NROW(v),jj]/v[1:(NROW(v)-1),jj])-1)
    }
    
    # log
    if(transf==5){
      vv[,jj] <- log(v[,jj])
    }
    
    # growth  
    if(transf==6){
      vtemp <- log(v[,jj])
      vv[,jj] <- c(NA, (vtemp[2:NROW(vtemp)]/vtemp[1:(NROW(vtemp)-1)])-1)
    }
    # log diff 
    if(transf==7){
      vv[,jj] <- c(NA, diff(log(v[,jj])))
    }
    
    
  }
  vv = ts(vv, start = start(data), end = end(data), frequency = frequency(data))
  colnames(vv) <- colnames(data)
  return(vv)
}
