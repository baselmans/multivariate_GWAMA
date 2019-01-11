# Nweighted GWAMA

Nweighted GWAMA is a R function that performs a multivariate GWAMA of genetically correlated traits while correcting for sample overlap. The details of the method is described in Baselmans et al. (Nature Genetics) http://dx.doi.org/10.1038/s41588-018-0320-8

* The current version is 1_2_3

# Getting Started

You can source the function in R using the following  line of code:
```source("https://github.com/baselmans/multivariate_GWAMA/blob/master/Downloads/N_weighted_GWAMA.function.1_2_3.R?raw=TRUE")```
 
Or alternatively, you can download the function in the folder "Downloads" and load it into R yourself.

# Requirements

To read in the GWAS summary data you might want to use the package ```data.table``` (https://cran.r-project.org/web/packages/data.table/index.html)  

With this package installed, you should be all set to install and use the Nweighted GWAMA function.

# How to test the function before deployment

A couple of examples can be found in the "downloads" directory. The number of SNPs for which we provide summary statistics in these test datasets is 100,000

# Using the Nweighted function

To run the Nweighted function The following files, agruments and data points are required:

1. The included summary statistics needs to be formatted with the following columns in the exact order:
   SNPID,CHR,BP,EA,OA,EAF,N,Z,PVAL
 ```  
   SNPID = RS number
   CHR = chromosome
   BP = base pair
   EA = effect allele
   OA = other allele
   EAF = frequency effect allele
   N = sample size
   Z = z-score
   P = p-value
```   
   Thus, the files need exactly (these) 9 columns

2. Make an empty list with length equal to the number of input files:
   For instance: 
```   
   dat<-vector("list",4)
```   
   Read in the data:
```   
   dat[[1]]<-fread("path_to_files",showProgress=F,data.table=F)
   dat[[2]]<-fread("path_to_files",showProgress=F,data.table=F)
   dat[[3]]<-fread("path_to_files",showProgress=F,data.table=F)
   dat[[4]]<-fread("path_to_files",showProgress=F,data.table=F)
```

3. To correct for sample overlap you need a matrix of cross trait intercepts (CTI), which you can get from from LD Score Regression (Bulik-Sullivan, Nature Genetics 2015).
   
   These Crosstrait intercepts can be found in the log files provided by LD Score Regression (genetic correlation)
   The output looks something like:
```   
   Summary of genetic correlation results						
   p1	p2	rg	se	z	p	h2_obs	h2_obs_se	h2_int	h2_int_se	gcov_int	gcov_int_se
   T1	T2	0.8419	0.3384	2.487	0.01287	0.0691	0.004	        1.0058  0.0094   	0.0045	        0.0048

```
   The value listed under gcov_int is what we call the Cross Trait Intercept in our paper
   
   A CTI matrix could look like:
```  	
         T1	 T2      T3      T4
T1	 1       0.1282	 0.0126	 0.0236
T2	 0.1282	 1       0.2189	 0.1605
T3	 0.0126	 0.2189	 1	 0.3474
T4	 0.0236	 0.1605	 0.3474	 1
```
4.  You need a vector of the SNP heritabilities of the included traits
    For instance: 
 ```
    h2list <- c(0.0498, 0.0441, 0.0708, 0.0294)
 ```   
* Note, The CTI values and the SNP heritabilities need to be in the exact same order as the GWAS summary statistics in your list.

# Run the Nweighted function
```  
  multivariate_GWAMA(x=dat,cov_Z=CTI,h2=h2list,out=".",name="your_name",output_gz=F,check_rows=F)
```  
  The function provides a log file which lists all the different steps performed by the function (e.g. aligning) and possible errors (see   file in the Downloads folder). 

# Model Averaging GWAMA

Model averaging GWAMA is R code that performs a multivariate GWAMA of genetically correlated traits while correcting for sample overlap. The details of the method is described in Baselmans et al. (Nature Genetics) http://dx.doi.org/10.1038/s41588-018-0320-8

Note: LD Score Regression has the assumption that the included test statistics follow a standard normal distribution under the null hypothesis of no effect. In MA GWAMA we can't guarantee that this assumption will be met. Interpreting results from LD Score regression should be done with some reservation. 

# Getting Started

You can source the function in R using the following  line of code:
```source(XXX)```
 
Or alternatively, you can download the function in the folder "Downloads" and load it into R yourself.

# Requirements

  MA GWAMA currently requiers seven R packages to be installed. These are 
```  
  sp (https://cran.r-project.org/web/packages/sp/index.html)
  rater (https://cran.r-project.org/web/packages/raster/index.html)
  RcppArmadillo (https://cran.r-project.org/web/packages/RcppArmadillo/index.html)
  unmarked (https://cran.r-project.org/web/packages/unmarked/index.html)
  VGAM (https://cran.r-project.org/web/packages/VGAM/index.html)
  AICcmodavg (https://cran.r-project.org/web/packages/AICcmodavg/index.html)
  metafor (https://cran.r-project.org/web/packages/metafor/index.html)
```  
 With these packages in place, you should be all set to install and use MA GWAMA
 
# How to test the function before deployment

  A couple of examples can be found in the "downloads" directory. The number of SNPs for which we provide summary statistics in these     test datasets is 1000 (due to more computational power the number of SNPs in the test dataset are limited to 1000)
  
# Using the Nweighted function

  To run the Nweighted function The following files, agruments and data points are required:  

   
   
   
   
