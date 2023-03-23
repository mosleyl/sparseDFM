# sparseDFM: An R Package to Estimate Dynamic Factor Models with Sparse Loadings

The aim of the package *SparseDFM* is to offer the tools for the R user to implement dynamic factor models with the option to induce sparse loadings. In summary we offer three different popular DFM estimation techniques plus the novel Sparse DFM estimation technique of Mosley et al. (2023). We call the estimation options:

* `PCA` - principal components analysis (PCA) for static factors seen in Stock and Watson (2002). 
* `2Stage` - the two-stage framework of PCA plus Kalman filter \& smoother seen in Giannone et al. (2008) and Doz et al. (2011).
* `EM` - the quasi-maximum likelihood approach using the EM algorithm to handle arbitrary patterns of missing data seen in Banbura and Modugno (2014).
* `EM-sarse` - the novel sparse EM approach allowing for sparse factor loadings seen in Mosley et al. (2023).

We allow the user the option of two different idiosyncratic error processes:

* `IID` - errors are IID white noise: $\epsilon_{i,t} \sim N(0, \sigma_i^2)$.
* `AR1` - errors follow an AR(1) process: $\epsilon_{i,t} = \phi_i \epsilon_{i,t-1} + e_{i,t}$ with $e_{i,t} \sim N(0,\sigma_i^2)$.

We also allow two different options for estimating the Kalman filter and smoother equations:

* `multivariate` - classic Kalman filter and smoother equations seen in  Shumway and Stoffer (1982).
* `univariate` - univariate treatment (sequential processing) of the multivariate equations for fast Kalman filter and smoother seen in  Koopman and Durbin (2000).

Alternative software that implement classic DFMs in R include the *MARSS* package of Holmeset al. (2012) which allow for more general state space structures. The nowcasting packages include *nowcasting* by  De Valk et al. (2019) and *nowcastDFM* by  Hopp (2021), that allow mixed-frequency time series nowcasting^[These nowcasting packages have now been removed from the CRAN repository. However, they are still accessible via GitHub.]. Recently, the package *dfms* by  Krantz and Bagdziunas (2022) has been uploaded to CRAN which implement the regular DFM methods with IID errors. They make use of C++ code and therefore is computationally faster than the previous three packages whose implementations are solely in R. Our package *sparseDFM* is novel in multiple ways. We implement the general DFM estimation methods found in the above packages also using C++ code for speed. We implement the univariate treatment of the Kalman filter/smoother equations of Koopman and Durbin (2000) which is computationally much faster than the classic multivariate approach when $p$ is large. We consider AR(1) idiosyncratic errors as well as IID. Finally, we allow the ability to estimate a sparse DFM using an efficient sparsified EM algorithm framework.