## data set with missing values replaced with cubic spline (between observed values) and 
## 1-D digital filter fitted values for missing tails calculated at the variable level.
## This is partly taken from nowcasting R package 


fill_NA <-function(X){
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  k <- 3
  idx.na <- is.na(X)
  
  # missing_row_idx <- (rowSums(idx.na)>N*0.8)
  # missing_row <- which(missing_row_idx)
  # 
  # if(sum(missing_row)!=0){
  #   X<-X[-missing_row,]
  # }
  # idx.na=is.na(X)
  
  for (i in 1:p){  
    
    x = X[,i]
    na_x = is.na(x)
    t1 = min(which(!na_x))
    t2 = max(which(!na_x))
    
    x1<-stats::spline(x[t1:t2],xout = 1:(t2-t1+1))
    xx<-x1$y
    x[t1:t2]<-x1$y
    na_x<-is.na(x)
    x[na_x] <- median(x,na.rm = T)
    
    x_MA<-stats::filter(x = c(rep(x[1],k),x,rep(x[length(x)],k)),filter = rep(1,2*k+1)/(2*k+1),sides = 1)
    
    x_MA=x_MA[(2*k+1):length(x_MA)]
    x[idx.na[,i]]=x_MA[idx.na[,i]]
    X[,i]=x;
  }
  
  
  return(list(X = X, idx.na=idx.na))
  
}
