#' @name summary.sparseDFM
#' @aliases print.sparseDFM
#' @aliases summary.sparseDFM
#' 
#' @title
#' sparseDFM Summary Outputs 
#' 
#' @description 
#' Summary and print outputs for class 'sparseDFM'.
#' 
#' @param x an object of class 'sparseDFM'
#' @param \dots Further \code{print} arguments. 
#' 
#' @returns 
#' Information on the model fitted.
#' 
#' @export

print.sparseDFM <- function(x,...){
  
  X = x$data$X
  A = x$params$A
  n = dim(X)[1]
  p = dim(X)[2]
  r = dim(A)[1]
  
  if(x$data$method=='PCA'){
    typeFM = 'Static'
  }else if(x$data$method=='EM-sparse'){
    typeFM = 'Sparse Dynamic'
  }else{
    typeFM = 'Dynamic'
  }
  
  cat('\n', typeFM, 'Factor Model using', x$data$method, 'with: \n\n n =', n,
      'observations, \n p =', p, 'variables, \n r =', r, 'factors, \n err =',
      x$data$err, 'idiosyncratic errors. \n\n Call summary() for estimation details.')
  
  
}

#' @rdname summary.sparseDFM
#' @param object an object of class 'sparseDFM'
#' @param \dots Further \code{summary} arguments.
#' @returns 
#' Summary information on estimation details. 
#' 
#' @export

summary.sparseDFM <- function(object,...){
    
  X = object$data$X
  A = object$params$A
  n = dim(X)[1]
  p = dim(X)[2]
  r = dim(A)[1]
  
  if(object$data$method=='PCA'){
    typeFM = 'Static'
  }else if(object$data$method=='EM-sparse'){
    typeFM = 'Sparse Dynamic'
  }else{
    typeFM = 'Dynamic'
  }
  
  cat('\nCall: \n\n', paste(deparse(object$data$call)))
  
  cat('\n\n',typeFM, 'Factor Model using', object$data$method, 'with: \n\n n =', n,
      'observations, \n p =', p, 'variables, \n r =', r, 'factors, \n err =',
      object$data$err)
  
  cat('\n\n The r x r factor transition matrix A \n')
  print(object$params$A)
  
  cat('\n\n The r x r factor transition error covariance matrix Sigma_u \n')
  print(object$params$Sigma_u)
  

}


#' @title 
#' sparseDFM Plot Outputs 
#' 
#' @description 
#' Make plots for the output of sparseDFM(). Options include:
#' \itemize{
#' \item \code{factor} - plot factor estimate series on top of the original standardized stationary data
#' \item \code{loading.heatmap} - make a heatmap of the loadings matrix
#' \item \code{loading.lineplot} - make a lineplot of variable loadings for a given factor
#' \item \code{loading.grouplineplot} - separate variable groups into colours for better visualisation 
#' \item \code{residual} - boxplot or scatterplot of residuals 
#' \item \code{lasso.bic} - BIC values for the LASSO tuning parameter
#' }
#' 
#' @param x an object of class 'sparseDFM'.
#' @param type character. The type of plot: \code{"factor"}, \code{"loading.heatmap"}, \code{"loading.lineplot"}, \code{"loading.grouplineplot"} or \code{"residual"}. Default is \code{"factor"}.
#' @param which.factors numeric vector of integers representing which factors should be plotted in \code{"factor"} and \code{"loading.heatmap"}. Default is \code{which.factors}=\code{1:(dim(x$state$factors)[2])}, plotting them all. Accepts a single integer if just one factor required. 
#' @param scale.factors logical. Standardize the factor estimates when plotting in \code{"factor"}. Default is \code{TRUE}.
#' @param which.series numeric vector of integers representing which series should be plotted in \code{"loading.heatmap"}, \code{"loading.lineplot"}, \code{"loading.grouplineplot"} and \code{"residual"}. Default is \code{which.series} = \code{1:(dim(x$params$Lambda)[1])}, plotting them all.
#' @param loading.factor integer. The factor to use in \code{"loading.lineplot"} and \code{"loading.grouplineplot"}. Default is 1. 
#' @param series.col character. The colour of the background series plotted in \code{"factor"}. Default is \code{series.col} = \code{"grey"}.
#' @param factor.col character. The colour of the factor estimate line in \code{"factor"}. Default is \code{factor.col} = \code{"black"}.
#' @param factor.lwd integer. The line width of the factor estimate line in \code{"factor"}. Default is \code{factor.lwd} = 2.
#' @param factor.lab vector of characters to label each factor in \code{"loading.heatmap"}. Default is \code{NULL} for standard labeling. 
#' @param use.series.names logical. Set to TRUE if plot should display series names in the data matrix X. Default is \code{FALSE} for numbered series. 
#' @param series.lab vector of characters to label each data series in \code{"loading.heatmap"}. Default is \code{NULL} for standard labeling.
#' @param series.labpos numeric vector of integers representing which series are labeled by \code{series.lab}. Default is \code{NULL} for standard labeling. 
#' @param colorkey logical. Display the colour key of the heatmap in \code{"loading.heatmap"}. Default is \code{TRUE}.
#' @param col.regions vector of gradually varying colors for \code{"loading.heatmap"}, see levelplot package. Default is \code{NULL} for standard colours. 
#' @param group.names vector of characters of the same dimension as \code{which.series} to represent the name of the group for each series in \code{"loading.grouplineplot"}. 
#' @param group.cols vector of characters of the same dimension as the number of different groups in \code{"loading.grouplineplot"} to represent the colours of the groups.
#' @param group.legend logical. Display the legend. Default is \code{TRUE}.
#' @param residual.type character. The type of residual plot: \code{"boxplot"} or \code{"scatterplot"}. Default is \code{"boxplot"}.
#' @param scatter.series integer. The series to plot when \code{residual.type} = \code{"scatterplot"}. Default is series 1. 
#' @param min.bic.col character. Colour for the best \eqn{\alpha}{\alpha} point. Default is \code{'red'}.
#' @param \dots for \code{plot.sparseDFM}. Further plot arguments. 
#' 
#' @returns 
#' Plots for the output of sparseDFM().
#' 
#' @importFrom Matrix Matrix image 
#' @importFrom ggplot2 ggplot aes geom_segment theme_light scale_color_manual theme element_text element_blank xlab ylab ggtitle geom_boxplot
#' @importFrom graphics par matplot lines box axis mtext boxplot plot abline
#' @importFrom reshape2 melt
#' 
#' @export 


plot.sparseDFM <- function(x, type = 'factor', which.factors = 1:(dim(x$state$factors)[2]), scale.factors = TRUE, 
                           which.series = 1:(dim(x$params$Lambda)[1]), loading.factor = 1, series.col = 'grey',
                           factor.col = 'black', factor.lwd = 2, factor.lab = NULL, use.series.names = FALSE, series.lab = NULL,
                           series.labpos = NULL,colorkey = TRUE, col.regions = NULL, group.names = NULL, group.cols = NULL,
                           group.legend = TRUE, residual.type = 'boxplot', scatter.series = 1, min.bic.col = 'red', ...){
  
  # global variable declaration 
  y = Group = Var2 = value = NULL
  
  ## type = 'factor'
  
  if(type == 'factor'){
    
    par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
    par(mfrow = c(1,1))
    
    dots <- list(...)
    
    if(x$data$standardize){
      data = x$data$X.bal
    }else{
      data = scale(x$data$X.bal)
    }
    
    if(scale.factors){
      factors = scale(x$state$factors)
    }else{
      factors = x$state$factors
    }
    
    r = dim(factors)[2]
    
    which.factors = sort(which.factors)
    
    k = length(which.factors)
    
    
    if(k > 1){
      
      oldpar <- par(mar = c(0, 5, 0, 2), oma = c(6, 0, 5, 0), mfrow = c(k, 1L))
      on.exit(par(oldpar))
      
      for(i in which.factors) {
        
        matplot(data, type = 'l', col = series.col, axes = FALSE, xlab = "", ylab = "")
        lines(factors[,i], type = 'l', col = factor.col, lwd = factor.lwd)
        box()
        axis(2, cex.axis = 1.2)
        if(i == max(which.factors)) axis(1, cex.axis = 1.2)
        mtext(paste("Factor", i), 2, line = 3, cex = 1.2)
        if(i == max(which.factors)) mtext(if (is.null(dots$xlab)) 'Time' else dots$xlab, side = 1, line = 3, cex = 1.2)
        
        
      }
      
      par(mfrow = c(1, 1))
      
      mtext(if(is.null(dots$main)) paste(if(scale.factors) "Standardized", "Factor Estimates and Standardized Data") else dots$main,
            side = 3, line = 2, outer=TRUE, cex = 1.3)
      
      
    }else{
      
      matplot(data, type = 'l', col = series.col, ylab = if (is.null(dots$ylab))'Value' else dots$ylab,
              xlab = if (is.null(dots$xlab)) 'Time' else dots$xlab,
              main = if(is.null(dots$main)) paste(if(scale.factors) "Standardized", "Factor", which.factors, "Estimate and Standardized Data") else dots$main,
              cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.2)
      lines(factors[,which.factors], col = factor.col, lwd = factor.lwd)
      
      
    }
    
    # re-set back to original 
    par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
    par(mfrow = c(1,1))
    
  }
  
  ## type = 'loading.heatmap'
  
  else if(type == 'loading.heatmap'){
    
    if(!is.null(series.lab)){
      if(is.null(series.labpos)){
        stop('Please provide a numeric vector indicating which series correspond to series.lab')
      }
    }
    
    
    # make sure it is of normal display 
    par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
    par(mfrow = c(1,1))
    
    which.factors = sort(which.factors)
    
    dots <- list(...)
    
    if(is.null(factor.lab)){
      
      factor.lab = c()
      
      for(i in which.factors){
        
        factor.lab[i] = paste0('F',i)
        
      }
      
    }
    
    trL = t(x$params$Lambda)
    lw = Matrix::Matrix(trL[which.factors,which.series], sparse = TRUE)
    
    series.names = unlist(dimnames(x$params$Lambda)[1])
    
    if(!use.series.names){
      series.names = NULL
    }
 
    Matrix::image(lw, xlab = if(is.null(dots$xlab))'Series' else dots$xlab,
                  ylab = if(is.null(dots$ylab))'Factor' else dots$ylab, sub = NULL,
                  main = if(is.null(dots$main)) 'Loadings Heatmap' else dots$main,
                  colorkey = colorkey, col.regions = col.regions,
                  scales = list(y=list(at = 1:(length(which.factors)), labels = factor.lab),
                                x=if(is.null(series.lab)){ list(at = 1:(length(which.series)), labels = if(!is.null(series.names)) series.names else which.series, rot = if(!is.null(series.names)) 90 else 0)} else {list(at = series.labpos, labels = series.lab, rot = 90)}))
    
    
  }
  
  ## type = 'loading.lineplot'
  
  else if(type == 'loading.lineplot'){
    
    par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
    par(mfrow = c(1,1))
    
    dots <- list(...)
    
    lw = t(x$params$Lambda)
    
    series.names = unlist(dimnames(x$params$Lambda)[1])
    
    if(!use.series.names){
      series.names = NULL
    }
    
    data <- data.frame(
      x=if(!is.null(series.names)) series.names else which.series,
      y=as.numeric(lw[loading.factor, which.series])
    )
    
    data$Group = rep('group', length(which.series)) 
    
    unique_group_names = 'group'
    
    mycolors = 'black'
    
    if(!is.null(series.names)){
      data$x = factor(data$x, levels = data$x)
    }
    
    
    # plot
    ggplot(data, aes(x=x, y=y)) +
      geom_segment( aes(x, xend=x, y=0, yend=y, color=Group), size=1.2, alpha=0.9) +
      theme_light() +
      scale_color_manual(breaks=unique_group_names,values = mycolors) +
      theme(
        legend.position = "none",
        panel.border = element_blank(),
        text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.text.x = if(!is.null(series.names)) element_text(angle = 90, vjust = 0.5, hjust=1) else element_text(angle = NULL, vjust = NULL, hjust=NULL) 
      ) +
      xlab(if(is.null(dots$xlab))"Series"else dots$xlab) +
      ylab(if(is.null(dots$ylab))"Loading Value"else dots$ylab) + 
      ggtitle(if(is.null(dots$main))paste("Loadings for Factor", loading.factor) else dots$main)
    
  }
  
  ## type = 'loading.grouplineplot'
  
  else if(type == 'loading.grouplineplot'){
    
    if(is.null(group.names)){
      stop("You must specify group.names")
    }
    
    if(is.null(group.cols)){
      stop("You must specify group.cols.")
    }
    
    par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
    par(mfrow = c(1,1))
    
    dots <- list(...)
    
    lw = t(x$params$Lambda)
    
    series.names = unlist(dimnames(x$params$Lambda)[1])
    
    if(!use.series.names){
      series.names = NULL
    }

    data <- data.frame(
      x= if(!is.null(series.names)) series.names else which.series,
      y=as.numeric(lw[loading.factor, which.series])
    )
    
    data$Group = group.names 
    
    unique_group_names = unique(group.names)
    
    mycolors = group.cols
  
    if(!is.null(series.names)){
      data$x = factor(data$x, levels = data$x)
    }
    
    # plot
    ggplot(data, aes(x=x, y=y)) +
      geom_segment( aes(x=x, xend=x, y=0, yend=y, color=Group), size=1.2, alpha=0.9) +
      theme_light() +
      scale_color_manual(breaks=unique_group_names,values = mycolors) +
      theme(
        legend.position = if(group.legend) "right" else "none",
        panel.border = element_blank(),
        text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.text.x = if(!is.null(series.names)) element_text(angle = 90, vjust = 0.5, hjust=1) else element_text(angle = NULL, vjust = NULL, hjust=NULL)
      ) +
      xlab(if(is.null(dots$xlab))"Series"else dots$xlab) +
      ylab(if(is.null(dots$ylab))"Loading Value"else dots$ylab) + 
      ggtitle(if(is.null(dots$main))paste("Loadings for Factor", loading.factor) else dots$main)
    
  }
  
  ## type = 'residual'
  
  else if(type == 'residual'){
    
    if(residual.type == 'scatter' && is.null(scatter.series)) stop("No series chosen in scatter.series")
    
    par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
    par(mfrow = c(1,1))
    
    dots = list(...)
    
    series.names = unlist(dimnames(x$params$Lambda)[1])
    
    if(!use.series.names){
      series.names = NULL
    }
    
    resids = x$data$X.bal - x$data$fitted
    
    if(residual.type=='boxplot'){
    
      data_long = melt(resids)
      
      data_long = data_long[,-1]

      # plot
      ggplot(data_long, aes(x=factor(Var2), y=value)) +
        geom_boxplot() +
        theme_light() +
        theme(
          legend.position = "none",
          panel.border = element_blank(),
          text = element_text(size = 10),
          axis.text = element_text(size = 10),
          axis.text.x = if(!is.null(series.names)) element_text(angle = 90, vjust = 0.5, hjust=1) else element_text(angle = NULL, vjust = NULL, hjust=NULL) 
        ) +
        xlab(if(is.null(dots$xlab))"Series"else dots$xlab) +
        ylab(if(is.null(dots$ylab))"Loading Value"else dots$ylab) + 
        ggtitle(if(is.null(dots$main))paste("Loadings for Factor", loading.factor) else dots$main)
    }else{
      
      s.names = unlist(dimnames(resids)[2])
      scatter.series.name = if(!is.null(s.names)) s.names[scatter.series] else paste('Series',scatter.series)
      
      plot(resids[,scatter.series], xlab = if(is.null(dots$xlab)) 'Time' else dots$xlab,
           ylab = if(is.null(dots$ylab)) 'Value' else dots$ylab, 
           main = if(is.null(dots$main)) paste('Residual Scatter Plot for', scatter.series.name) else dots$main,
           cex.axis = if(is.null(dots$cex.axis)) 1.2 else dots$cex.axis, cex.main = if(is.null(dots$cex.main)) 1.2 else dots$cex.main, cex.lab = if(is.null(dots$cex.lab)) 1.2 else dots$cex.lab, pch = if(is.null(dots$main)) 20 else dots$pch) 
      
      
    }
    
  }
  
  ## type = 'lasso.bic'
  else{
    
    par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
    par(mfrow = c(1,1))
    
    dots = list(...)
    
    
    mycols = rep(if(is.null(dots$col)) 'black' else dots$col, length(x$em$alpha_grid))
    mycols[which.min(x$em$bic)] = min.bic.col
    
    plot(log10(x$em$alpha_grid),x$em$bic, type = if(is.null(dots$type)) 'o' else dots$type, pch = if(is.null(dots$pch)) 20 else dots$pch, col = mycols, xlab =  if(is.null(dots$xlab)) expression(paste("log10(", alpha,')')) else dots$xlab, ylab = if(is.null(dots$ylab)) 'BIC' else dots$ylab, main = if(is.null(dots$main)) 'BIC Values for the LASSO Tuning Parameter' else dots$main)
    
    abline(v=log10(x$em$alpha_opt), col = min.bic.col, lty = 'dashed')
    
  }
  
}


#' @name residuals.sparseDFM
#' @aliases residuals.sparseDFM
#' @aliases resid.sparseDFM
#' @aliases fitted.sparseDFM
#' 
#' @title 
#' sparseDFM Residuals and Fitted Values 
#' 
#' @description 
#' Obtain the residuals or fitted values of the \code{sparseDFM} fit. 
#' 
#' @param object an object of class 'sparseDFM'.
#' @param standardize logical. The residuals and fitted values should be standardized. Default is \code{FALSE}, values returned in the original data \eqn{\bm{X}}{X} scale.
#' @param \dots Further \code{fitted} arguments.
#' 
#' @return Residuals or fitted values of \code{sparseDFM}.
#' 
#' @rdname residuals.sparseDFM
#' @export


fitted.sparseDFM <- function(object, standardize = FALSE,...){
  
  if(standardize){
    return(object$data$fitted)
  }else{
    return(object$data$fitted.unscaled)
  }
  
}


#' @rdname residuals.sparseDFM
#' 
#' @param object an object of class 'sparseDFM'.
#' @param standardize logical. The residuals and fitted values should be standardized. Default is \code{FALSE}, values returned in the original data \eqn{\bm{X}}{X} scale.
#' @param \dots Further \code{residuals} arguments.
#' 
#' @export

residuals.sparseDFM <- function(object, standardize = FALSE,...){
  
  if(standardize){
    res = object$data$X.bal - object$data$fitted
  }else{
    n = dim(object$data$X.bal)[1]
    X.bal_raw = kronecker(t(object$data$X.sd),rep(1,n))*object$data$X.bal + kronecker(t(object$data$X.mean),rep(1,n))
    res = X.bal_raw - object$data$fitted.unscaled
  }
  return(res)
  
}

#' @name predict.sparseDFM
#' @aliases predict.sparseDFM
#' @aliases print.sparseDFM_forecast
#' 
#' @title 
#' Forecasting factor estimates and data series. 
#' 
#' @description 
#' Predict the next h steps ahead for the factor estimates and the data series. Given information up to time \eqn{t}{t}, a h-step ahead forecast is \eqn{\bm{X}_{t+h}=\bm{\Lambda}\bm{A}^{h}\bm{F}_t+\bm{\Phi}^h\bm{\epsilon}_t}{X_{t+h}=\Lambda A^h F_t+\Phi^h \epsilon_t}, where \eqn{\bm{\Phi}=0}{\Phi = 0} for the IID idiosyncratic error case.
#' 
#' @param object an object of class 'sparseDFM'.
#' @param h integer. The number of steps ahead to compute the forecast for. Default is \eqn{h=1}{h=1}.
#' @param standardize logical. Returns data series forecasts in the original data scale if set to \code{FALSE}. Default is \code{FALSE}. 
#' @param \dots Further \code{predict} arguments.
#' 
#' @return X_hat \eqn{h \times p}{h x p} numeric matrix of data series forecasts.
#' @return F_hat \eqn{h \times r}{h x r} numeric matrix of factor forecasts.
#' @return e_hat \eqn{h \times p}{h x p} numeric matrix of AR(1) idiosyncratic error forecasts if \code{err}=\code{AR1} in \code{sparseDFM}.
#' @return h forecasts produced for h steps ahead.
#' @return err the type of idiosyncratic errors used in \code{sparseDFM}.
#' 
#' @export

predict.sparseDFM <- function(object, h = 1, standardize = FALSE,...){
  
  A = object$params$A
  Lambda = object$params$Lambda
  n = dim(object$state$factors)[1]
  r = dim(object$state$factors)[2]
  p = dim(object$data$X.bal)[2]
  F_old = object$state$factors[n,]
  
  if(object$data$err == 'AR1'){
    Phi = object$params$Phi
    e_old = object$state$errors
    e_new = matrix(NA, nrow = h, ncol = p)
  }
  
  F_new = matrix(NA, nrow = h, ncol = r)
  X_new = matrix(NA, nrow = h, ncol = p)
  
  for(i in 1:h){
    
    F_new[i,] = A %*% F_old 
    X_new[i,] = F_new[i,] %*% t(Lambda)
    
    if(object$data$err == 'AR1'){
      
      e_new[i,] = Phi %*% e_old
      X_new[i,] = X_new[i,] + e_new[i,]
      e_old = e_new[i,]
      
    }
    
    F_old = F_new[i,]
    
  }
  
  if(standardize){
    X_new = X_new
  }else{
    X_new = kronecker(t(object$data$X.sd),rep(1,h))*X_new + kronecker(t(object$data$X.mean),rep(1,h))
  }
  
  if(object$data$err == 'AR1'){
    output = list(X_hat = X_new,
                  F_hat = F_new,
                  e_hat = e_new,
                  h = h,
                  err = object$data$err)
  }else{
    output = list(X_hat = X_new,
                  F_hat = F_new,
                  h = h,
                  err = object$data$err)
  }
  
  class(output) <- "sparseDFM_forecast"
  return(output)
  
}


#' @rdname predict.sparseDFM 
#' @param x an object of class 'sparseDFM_forecast' from \code{predict.sparseDFM}.
#' @param \dots Further \code{print} arguments.
#' @returns 
#' Prints out the h-step ahead forecast from \code{predict.sparseDFM}.
#' 
#' @export

print.sparseDFM_forecast <- function(x,...){
  
  h = x$h
  X_hat = x$X_hat
  F_hat = x$F_hat
  
  cat('\n The', h, 'step ahead forecast for the data series are \n')
  print(X_hat)
  
  cat('\n The', h, 'step ahead forecast for the factors are \n')
  print(F_hat)
  
  if(x$err == 'AR1'){
    cat('\n The', h, 'step ahead forecast for the AR(1) idiosyncratic errors are \n')
    print(x$e_hat)
  }
  
}




