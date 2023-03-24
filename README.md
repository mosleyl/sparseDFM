# sparseDFM: An R Package to Estimate Dynamic Factor Models with Sparse Loadings

The aim of the package **SparseDFM** is to offer the tools for the R user to implement dynamic factor models with the option to induce sparse loadings. In summary we offer three different popular DFM estimation techniques plus the novel Sparse DFM estimation technique of Mosley et al. (2023a). We call the estimation options:

* `PCA` - principal components analysis (PCA) for static factors seen in Stock and Watson (2002). 
* `2Stage` - the two-stage framework of PCA plus Kalman filter \& smoother seen in Giannone et al. (2008) and Doz et al. (2011).
* `EM` - the quasi-maximum likelihood approach using the EM algorithm to handle arbitrary patterns of missing data seen in Banbura and Modugno (2014).
* `EM-sparse` - the novel sparse EM approach allowing for sparse factor loadings seen in Mosley et al. (2023a).

We allow the user the option of two different idiosyncratic error processes:

* `IID` - errors are IID white noise: $\epsilon_{i,t} \sim N(0, \sigma_i^2)$.
* `AR1` - errors follow an AR(1) process: $\epsilon_{i,t} = \phi_i \epsilon_{i,t-1} + e_{i,t}$ with $e_{i,t} \sim N(0,\sigma_i^2)$.

We also allow two different options for estimating the Kalman filter and smoother equations:

* `multivariate` - classic Kalman filter and smoother equations seen in  Shumway and Stoffer (1982).
* `univariate` - univariate treatment (sequential processing) of the multivariate equations for fast Kalman filter and smoother seen in  Koopman and Durbin (2000).

Alternative software that implement classic DFMs in R include the **MARSS** package of Holmeset al. (2012) which allow for more general state space structures. The nowcasting packages include **nowcasting** by  De Valk et al. (2019) and **nowcastDFM** by  Hopp (2021), that allow mixed-frequency time series nowcasting[*Note*]. Recently, the package **dfms** by  Krantz and Bagdziunas (2022) has been uploaded to CRAN which implement the regular DFM methods with IID errors. They make use of C++ code and therefore is computationally faster than the previous three packages whose implementations are solely in R. Our package *sparseDFM* is novel in multiple ways. We implement the general DFM estimation methods found in the above packages also using C++ code for speed. We implement the univariate treatment of the Kalman filter/smoother equations of Koopman and Durbin (2000) which is computationally much faster than the classic multivariate approach when $p$ is large. We consider AR(1) idiosyncratic errors as well as IID. Finally, we allow the ability to estimate a sparse DFM using an efficient sparsified EM algorithm framework.

## Tutorials on Package Use

Installation:

* Install the package into R using: `install.packages('sparseDFM')`
* Load into the R environment using: `library(sparseDFM)`

We provide two vignettes in the package:

1. Analysing quarterly CPI (consumer price inflation) index data for the UK
2. Nowcasting Trade in Goods for UK exports of the 9 main commodities worldwide

that can be accessed using `browseVignettes('sparseDFM')`.

We provide a detailed software paper on the package:

* See Mosley et al. (2023b)

## References 
* Banbura M, Modugno M (2014). “Maximum Likelihood Estimation of Factor Models on Datasets with Arbitrary Pattern of Missing Data.” Journal of Applied Econometrics, 29(1), 133–160.
* De Valk S, de Mattos D, Ferreira P (2019). “Nowcasting: An R Package for Predicting Economic Variables Using Dynamic Factor Models.” The R Journal, 11(1), 230–244
* Doz C, Giannone D, Reichlin L (2011). “A Two-Step Estimator for Large Approximate Dynamic Factor Models Based on Kalman Filtering.” Journal of Econometrics, 164(1), 188–205.
* Giannone D, Reichlin L, Small D (2008). “Nowcasting: The Real-Time Informational Content of Macroeconomic Data.” Journal of Monetary Economics, 55(4), 665–676
* Holmes EE, Ward EJ, Kellie W (2012). “MARSS: Multivariate Autoregressive State-Space Models for Analyzing Time-Series Data.” The R Journal, 4(1), 11.
* Hopp D (2021). “nowcastDFM.” https://github.com/dhopp1/nowcastDFM
* Koopman SJ, Durbin J (2000). “Fast Filtering and Smoothing for Multivariate State Space Models.” Journal of Time Series Analysis, 21(3), 281–296
* Krantz S, Bagdziunas R (2022). dfms: Dynamic Factor Models. R Package Version 0.1.3, URL https://sebkrantz.github.io/dfms/
* Mosley L, Chan TS, Gibberd A (2023a) The Sparse Dynamic Factor Model: A Regularised Quasi-Maximum Likelihood Approach. arXiv Preprint
* Mosley L, Chan TS, Gibberd A (2023b) sparseDFM: An R Package to Estimate Dynamic Factor Models with Sparse Loadings. arXiv Preprint
* Shumway RH, Stoffer DS (1982). “An Approach to Time Series Smoothing and Forecasting Using the EM Algorithm.” Journal of Time Series Analysis, 3(4), 253–264
* Stock JH, Watson MW (2002). “Forecasting Using Principal Components from a Large Number of Predictors.” Journal of the American Statistical Association, 97(460), 1167–1179
* [*Note*] These nowcasting packages have now been removed from the CRAN repository. However, they are still accessible via GitHub.