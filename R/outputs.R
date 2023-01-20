#' @name summary.SparseDFM
#' @aliases print.SparseDFM
#' @aliases summary.SparseDFM
#' 
#' @title
#' SparseDFM Summary Outputs 
#' 
#' @description 
#' Summary and print outputs for class 'SparseDFM'.
#' 
#' @param x an object of class 'SparseDFM'
#' 
#' @returns 
#' Information on the model fitted.
#' 
#' @export

print.SparseDFM <- function(x){
  
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

#' @rdname summary.SparseDFM
#' @param x an object of class 'SparseDFM'
#' @returns 
#' Summary information on estimation details. 
#' 
#' @export

summary.SparseDFM <- function(x){
    
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
  
  cat('\nCall: \n\n', paste(deparse(x$data$call)))
  
  cat('\n\n',typeFM, 'Factor Model using', x$data$method, 'with: \n\n n =', n,
      'observations, \n p =', p, 'variables, \n r =', r, 'factors, \n err =',
      x$data$err)
  
  cat('\n\n The r x r factor transition matrix A \n')
  print(x$params$A)
  
  cat('\n\n The r x r factor transition error covariance matrix Sigma_u \n')
  print(x$params$Sigma_u)
  

}


#' @title 
#' SparseDFM Plot Outputs 
#' 
#' @description 
#' Make plots for the output of SparseDFM(). Options include:
#' \itemize{
#' \item \code{factor} - plot factor estimate series on top of the original standardized stationary data
#' \item \code{loading.heatmap} - make a heatmap of the loadings matrix
#' \item \code{loading.lineplot} - make a lineplot of variable loadings for a given factor
#' \item \code{loading.grouplineplot} - separate variable groups into colours for better visualisation 
#' \item \code{residual} - boxplot or scatterplot of residuals 
#' }
#' 
#' @param x an object of class 'SparseDFM'.
#' @param type character. The type of plot: \code{"factor"}, \code{"loading.heatmap"}, \code{"loading.lineplot"}, \code{"loading.grouplineplot"} or \code{"residual"}. Default is \code{"factor"}.
#' @param which.factors numeric vector of integers representing which factors should be plotted in \code{"factor"} and \code{"loading.heatmap"}. Default is which.factors = 1:(dim(x$state$factors)[2]), plotting them all. Accepts a single integer if just one factor required. 
#' @param scale.factors logical. Standardize the factor estimates when plotting in \code{"factor"}. Default is \code{TRUE}.
#' @param which.series numeric vector of integers representing which series should be plotted in \code{"loading.heatmap"}, \code{"loading.lineplot"}, \code{"loading.grouplineplot"} and \code{"residual"}. Default is which.series = 1:(dim(x$params$Lambda)[1]), plotting them all.
#' @param loading.factor integer. The factor to use in \code{"loading.lineplot"} and \code{"loading.grouplineplot"}. Default is 1. 
#' @param series.col character. The colour of the background series plotted in \code{"factor"}. Default is \code{series.col} = \code{"grey"}.
#' @param factor.col character. The colour of the factor estimate line in \code{"factor"}. Default is \code{factor.col} = \code{"black"}.
#' @param factor.lwd integer. The line width of the factor estimate line in \code{"factor"}. Default is \code{factor.lwd} = 2.
#' @param factor.lab vector of characters to label each factor in \code{"loading.heatmap"}. Default is \code{NULL} for standard labeling. 
#' @param series.lab vector of characters to label each data series in \code{"loading.heatmap"}. Default is \code{NULL} for standard labeling.
#' @param series.labpos numeric vector of integers representing which series are labeled by \code{series.lab}. Default is \code{NULL} for standard labeling. 
#' @param colorkey logical. Display the colour key of the heatmap in \code{"loading.heatmap"}. Default is \code{TRUE}.
#' @param col.regions vector of gradually varying colors for \code{"loading.heatmap"}, see levelplot package. Default is \code{NULL} for standard colours. 
#' @param group.names vector of characters of the same dimension as \code{which.series} to represent the name of the group for each series in \code{"loading.grouplineplot"}. 
#' @param group.cols vector of characters of the same dimension as the number of different groups in \code{"loading.grouplineplot"} to represent the colours of the groups.
#' @param group.legend logical. Display the legend. Default is \code{TRUE}.
#' @param residual.type character. The type of residual plot: \code{"boxplot"} or \code{"scatterplot"}. Default is \code{"boxplot"}.
#' @param scatter.series integer. The series to plot when \code{residual.type} = \code{"scatterplot"}. Default is series 1. 
#' @param \dots for \code{plot.SparseDFM}. Further plot arguments. 
#' 
#' @returns 
#' Plots for the output of SparseDFM().
#' 
#' @importFrom Matrix Matrix image 
#' @importFrom ggplot2 ggplot aes geom_segment theme_light scale_color_manual theme element_text element_blank xlab ylab ggtitle 
#' @importFrom graphics par matplot lines box axis mtext boxplot plot 
#' 
#' @export 


plot.SparseDFM <- function(x, type = 'factor', which.factors = 1:(dim(x$state$factors)[2]), scale.factors = TRUE, 
                           which.series = 1:(dim(x$params$Lambda)[1]), loading.factor = 1, series.col = 'grey',
                           factor.col = 'black', factor.lwd = 2, factor.lab = NULL, series.lab = NULL, series.labpos = NULL,
                           colorkey = TRUE, col.regions = NULL, group.names = NULL, group.cols = NULL, 
                           group.legend = TRUE, residual.type = 'boxplot', scatter.series = 1, ...){
  
  # Do warning checks 
  
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
        
        factor.lab[i] = paste('F',i)
        
      }
      
    }
    
    trL = t(x$params$Lambda)
    lw = Matrix::Matrix(trL[which.factors,which.series], sparse = TRUE)
    
    
    Matrix::image(lw, xlab = if(is.null(dots$xlab))'Series' else dots$xlab,
                  ylab = if(is.null(dots$ylab))'Factor' else dots$ylab, sub = NULL,
                  main = if(is.null(dots$main)) 'Loadings Heatmap' else dots$main,
                  colorkey = colorkey, col.regions = col.regions,
                  scales = list(y=list(at = 1:(length(which.factors)), labels = factor.lab),
                                x=if(is.null(series.lab)) list(at = 1:(length(which.series)), labels = which.series) else list(at = series.labpos, labels = series.lab)))
    
    
  }
  
  ## type = 'loading.lineplot'
  
  else if(type == 'loading.lineplot'){
    
    par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
    par(mfrow = c(1,1))
    
    dots <- list(...)
    
    trL = t(x$params$Lambda)
    lw = Matrix::Matrix(trL, sparse = TRUE)
    
    
    data <- data.frame(
      x=which.series,
      y=as.numeric(lw[loading.factor, which.series])
    )
    
    data$Group = rep('group', length(which.series)) 
    
    unique_group_names = 'group'
    
    mycolors = 'black'
    
    
    # plot
    ggplot(data, aes(x=x, y=y)) +
      geom_segment( aes(x=x, xend=x, y=0, yend=y, color=Group), size=1.2, alpha=0.9) +
      theme_light() +
      scale_color_manual(breaks=unique_group_names,values = mycolors) +
      theme(
        legend.position = "none",
        panel.border = element_blank(),
        text = element_text(size = 15),
        axis.text = element_text(size = 15)
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
    
    trL = t(x$params$Lambda)
    lw = Matrix::Matrix(trL, sparse = TRUE)
    
    
    data <- data.frame(
      x=which.series,
      y=as.numeric(lw[loading.factor, which.series])
    )
    
    data$Group = group.names 
    
    unique_group_names = unique(group.names)
    
    mycolors = group.cols
    
    
    # plot
    ggplot(data, aes(x=x, y=y)) +
      geom_segment( aes(x=x, xend=x, y=0, yend=y, color=Group), size=1.2, alpha=0.9) +
      theme_light() +
      scale_color_manual(breaks=unique_group_names,values = mycolors) +
      theme(
        legend.position = if(group.legend) "right" else "none",
        panel.border = element_blank(),
        text = element_text(size = 17),
        axis.text = element_text(size = 17)
      ) +
      xlab(if(is.null(dots$xlab))"Series"else dots$xlab) +
      ylab(if(is.null(dots$ylab))"Loading Value"else dots$ylab) + 
      ggtitle(if(is.null(dots$main))paste("Loadings for Factor", loading.factor) else dots$main)
    
  }
  
  ## type = 'residual'
  
  else{
    
    if(residual.type == 'scatter' && is.null(scatter.series)) stop("No series chosen in scatter.series")
    
    par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
    par(mfrow = c(1,1))
    
    dots = list(...)
    
    resids = x$data$X.bal - x$data$fitted.scaled
    
    
    if(residual.type=='boxplot'){
      
      boxplot.series = sort(which.series)
      
      boxplot(resids[,boxplot.series], xlab = if(is.null(dots$xlab)) 'Series' else dots$xlab,
              ylab = if(is.null(dots$ylab)) 'Value' else dots$ylab, 
              main = if(is.null(dots$main)) paste('Boxplot of Residuals', scatter.series),
              cex.axis = 1.2, cex.main = 1.2, cex.lab = 1.2, axes = FALSE)
      box()
      axis(1, at = 1:length(boxplot.series), labels = boxplot.series, xpd = NA)
      axis(2, xpd = NA)
      
    }else{
      
      plot(resids[,scatter.series], xlab = if(is.null(dots$xlab)) 'Time' else dots$xlab,
           ylab = if(is.null(dots$ylab)) 'Value' else dots$ylab, 
           main = if(is.null(dots$main)) paste('Residual Scatter Plot for Series', scatter.series),
           cex.axis = 1.2, cex.main = 1.2, cex.lab = 1.2)
      
    }
    
  }
  
}


#' @name residuals.SparseDFM
#' @aliases residuals.SparseDFM
#' @aliases resid.SparseDFM
#' @aliases fitted.SparseDFM
#' 
#' @title 
#' SparseDFM Residuals and Fitted Values 
#' 
#' @description 
#' Obtain the residuals or fitted values of the \code{SparseDFM} fit. 
#' 
#' @param x an object of class 'SparseDFM'.
#' @param standardize logical. The residuals and fitted values should be standardized. Default is \code{FALSE}, values returned in the original data \eqn{\bm{X}}{X} scale.
#' 
#' @return Residuals or fitted values of \code{SparseDFM}.
#' 
#' @export


fitted.SparseDFM <- function(x, standardize = FALSE){
  
  if(standardize){
    return(x$data$fitted.scaled)
  }else{
    return(x$data$fitted)
  }
  
}


#' @rdname residuals.SparseDFM
#' 
#' @param x an object of class 'SparseDFM'.
#' @param standardize logical. The residuals and fitted values should be standardized. Default is \code{FALSE}, values returned in the original data \eqn{\bm{X}}{X} scale.
#' 
#' @export

residuals.SparseDFM <- function(x, standardize = FALSE){
  
  if(standardize){
    res = x$data$X.bal - x$data$fitted
  }else{
    n = dim(x$data$X.bal)[1]
    X.bal_raw = kronecker(t(x$data$X.sd),rep(1,n))*x$data$X.bal + kronecker(t(x$data$X.mean),rep(1,n))
    res = X.bal_raw - x$data$fitted
  }
  return(res)
  
}









