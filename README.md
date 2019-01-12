# Nweighted GWAMA

Nweighted GWAMA is a R function that performs a multivariate GWAMA of genetically correlated traits while correcting for sample overlap. The details of the method is described in Baselmans et al. (Nature Genetics) http://dx.doi.org/10.1038/s41588-018-0320-8

* The current version is 1_2_3

# Getting Started

You can source the function in R using the following  line of code:
```source("https://github.com/baselmans/multivariate_GWAMA/blob/master/Test_Data/N_weighted_GWAMA.function.1_2_3.R?raw=TRUE")```
 
Or alternatively, you can download the function in the folder "Downloads" and load it into R yourself.

# Requirements

To read in the GWAS summary data you might want to use the package ```data.table``` (https://cran.r-project.org/web/packages/data.table/index.html)  

With this package installed, you should be all set to install and use the Nweighted GWAMA function.

# How to test the function before deployment

A couple of examples can be found in the "downloads" directory. The number of SNPs for which we provide summary statistics in these test datasets is 100,000

You can download a test script in the Test_Data folder called: example_script.R

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
dat[[1]]<-fread("https://github.com/baselmans/multivariate_GWAMA/blob/master/Test_Data/LS_100K_no23andMe.txt?raw=TRUE",showProgress=F,data.table=F)
dat[[2]]<-fread("https://github.com/baselmans/multivariate_GWAMA/blob/master/Test_Data/PA_100K_no23andMe.txt?raw=TRUE",showProgress=F,data.table=F)
dat[[3]]<-fread("https://github.com/baselmans/multivariate_GWAMA/blob/master/Test_Data/NEU_100K_no23andMe.txt?raw=TRUE",showProgress=F,data.table=F)
dat[[4]]<-fread("https://github.com/baselmans/multivariate_GWAMA/blob/master/Test_Data/DEP_100K_no23andMe.txt?raw=TRUE",showProgress=F,data.table=F)


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
```
CTI <- as.matrix(read.table("https://github.com/baselmans/multivariate_GWAMA/blob/master/Test_Data/cross_trait_intercept", header =T, row.names = 1))
```

4.  You need a vector of the SNP heritabilities of the included traits
    For instance: 
 ```
    h2list <- c(0.0498, 0.0441, 0.0708, 0.0294)
 ```   
* Note, The CTI values and the SNP heritabilities need to be in the exact same order as the GWAS summary statistics in your list.

# Run the Nweighted function
```  
  multivariate_GWAMA(x=dat,cov_Z=CTI,h2=h2list,out=".",name="NGWAMA_test",output_gz=F,check_rows=F)
```  
  The function provides a log file which lists all the different steps performed by the function (e.g. aligning) and possible errors (see   file in the Downloads folder). 

# Model Averaging GWAMA

Model averaging GWAMA is R code that performs a multivariate GWAMA of genetically correlated traits while correcting for sample overlap. The details of the method is described in Baselmans et al. (Nature Genetics) http://dx.doi.org/10.1038/s41588-018-0320-8

Note: LD Score Regression has the assumption that the included test statistics follow a standard normal distribution under the null hypothesis of no effect. In MA GWAMA we can't guarantee that this assumption will be met. Interpreting results from LD Score regression should be done with some reservation. (Automated function will follow as soon a possible)

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

  A couple of examples can be found in the "downloads" directory. The number of SNPs for which we provide summary statistics in these     test datasets is 1000 (due to more computational power requiered for MA GWAMA, the number of SNPs in the test dataset are limited to     1000)
  
#  Using MA GWAMA

   To run the Nweighted function The following files, agruments and data points are required:  

   1. The included summary statistics needs to be formatted with the following columns in the exact order:
   SNPID,CHR,BP,EA,OA,EAF,N,Z,PVAL
 ```  
   cptid = cptid (chromosome:basepare)
   RS = RS number 
   CHR = chromosome
   BP = base pair
   A1 = effect allele
   A2 = other allele
   EAF = frequency effect allele
   N = sample size
   Z = z-score
   PVAL = p-value
   Beta = regression coefficient
   SE = standard error
```
   2. Reading the data
```   
   T1<-fread("path_to_files",showProgress=F,data.table=F)
   T2<-fread("path_to_files",showProgress=F,data.table=F)
   T3<-fread("path_to_files",showProgress=F,data.table=F)
   T4<-fread("path_to_files",showProgress=F,data.table=F)
```   
 # Reformatting the files
 
   3. select right columns (rs, beta, se)

```
   T1_b1_se1 <- T1[c(2,11:12)]
   T2_b2_se2 <- T2[c(2,11:12)]
   T3_b3_se3 <- T3[c(2,11:12)]
   T4_b4_se4 <- T4[c(2,11:12)]
```

   4. Merge data
```
   M1 <- merge(T1_b1_se1,T2_b2_se2,by=1)
   M1 <- merge(M1,T3_b3_se3,by=1)
   M1 <- merge(M1,T4_b4_se4,by=1)
```   
   Omit missing data
```
    M1 <- na.omit(M1)
```
```
   A1 <- merge(T1,M1, by="SNPID")
   A1 <- A1[row.names(unique(A1["SNPID"])),]

   B <- M1[,c(2,4,6,8)]
   SE <- M1[,c(3,5,7,9)]
```

   5. Read in the Cross Trait Intercept to correct for sample overlap (see Nweighted GWAMA)
```   
   cov_Z <- read.table("your_CTI.txt", header =T, sep = " ")
   cov_Z <- as.matrix(cov_Z) 
```
   Now your files are formatted in the right way and MA GWAMA could be run.
   
   # Running MA GWAMA
```
  The following for loop will perform the MA GWAMA analysis
   
   for(i in 1:nrow(B)){
  
  V <- diag(SE[i,]) %*% cov_Z %*% diag(SE[i,])
  
  
  yi <- as.numeric(B[i,])
  
  
  grid <- expand.grid(c(0,1),c(0,1),c(0,1),c(0,1))
  
  grid <- t(grid)
  
  #h2 = vector sqrt(h2)
  h2 <- as.vector(c(sqrt(0.0498),sqrt(0.0441), sqrt(0.0748), sqrt(0.0305)))
  m31 <- rma.mv(yi~0+h2, V=V,  method="ML")
  
  # 2 effects
  m32 <- rma.mv(yi~0+h2 + I(grid[,2]*h2), V=V, method="ML")
  m33 <- rma.mv(yi~0+h2 + I(grid[,3]*h2), V=V, method="ML")
  m34 <- rma.mv(yi~0+h2 + I(grid[,4]*h2), V=V, method="ML")
  m35 <- rma.mv(yi~0+h2 + I(grid[,5]*h2), V=V, method="ML")
  m36 <- rma.mv(yi~0+h2 + I(grid[,6]*h2), V=V, method="ML")
  m37 <- rma.mv(yi~0+h2 + I(grid[,7]*h2), V=V, method="ML")
  m38 <- rma.mv(yi~0+h2 + I(grid[,8]*h2), V=V, method="ML")
  
  
  
  Mods <- list(m31,m32,m33,m34,m35,m36,m37,m38)
  
  K <- c(rep(1,1),rep(2,7))
  
  LL <- lapply(Mods,logLik)
  LL <- unlist(LL)
  aictabCustom(LL,K,nobs=4)
  
  XX <- sapply(Mods,predict)
  
  est2 <- matrix(unlist(XX[1,1]),4,1)
  est3 <- matrix(unlist(XX[1,2:8]),4,7)
  
  est<- cbind(est2,est3)  
  
  se2 <- matrix(unlist(XX[2,1]),4,1)
  se3 <- matrix(unlist(XX[2,2:8]),4,7)
  
  se<- cbind(se2,se3)
  
  y1 <- modavgCustom(LL,K,nobs=4, estimate=est[1,] ,se=se[1,],second.ord = T)
  y2 <- modavgCustom(LL,K,nobs=4, estimate=est[2,] ,se=se[2,],second.ord = T)
  y3 <- modavgCustom(LL,K,nobs=4, estimate=est[3,] ,se=se[3,],second.ord = T)
  y4 <- modavgCustom(LL,K,nobs=4, estimate=est[4,] ,se=se[4,],second.ord = T)
  
  output_SNP[i, ] <- c(y1$Mod.avg.est,y1$Uncond.SE,y2$Mod.avg.est,y2$Uncond.SE,y3$Mod.avg.est,y3$Uncond.SE,y4$Mod.avg.est,y4$Uncond.SE)
  
  
}   

```

# Generate Output

```
Output <- cbind(A1[,c(1,2,3,4,5,6)],output_SNP,pchisq((output_SNP[,1]/output_SNP[,2])^2,1,lower=F),
                pchisq((output_SNP[,3]/output_SNP[,4])^2,1,lower=F),pchisq((output_SNP[,5]/output_SNP[,6])^2,
                1,lower=F),pchisq((output_SNP[,7]/output_SNP[,8])^2,1,lower=F))

names(Output) <- c("RS","cptid","CHR","BP","A1","A2","Beta_LS","SE_LS","Beta_PA","SE_PA","Beta_NEU","SE_NEU","Beta_DEP","SE_DEP","P_LS","P_PA","P_NEU","P_DEP")

write.table("AVGWAS.txt" ,quote=F,row.names=F,col.names=T)

```


