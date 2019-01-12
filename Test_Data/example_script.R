library(data.table)
source("https://github.com/baselmans/multivariate_GWAMA/blob/master/Test_Data/N_weighted_GWAMA.function.1_2_2.R?raw=TRUE")

## matrix with cross-trait-intercepts on off-diagonal, and LDSC intercepts on the diagonal
CTI <- as.matrix(read.table("cross_trait_intercept", header =T, row.names = 1))

## input for the function must be collected in a list
## make empty list
dat<-vector("list",4) ## change "4" according to the number of input files

## read in data
## columns must have this order: "SNPID","CHR","BP","EA","OA","EAF","N","Z","PVAL"
dat[[1]]<-fread("https://github.com/baselmans/multivariate_GWAMA/blob/master/Test_Data/LS_100K_no23andMe.txt",showProgress=F,data.table=F)
dat[[2]]<-fread("https://github.com/baselmans/multivariate_GWAMA/blob/master/Test_Data/PA_100K_no23andMe.txt",showProgress=F,data.table=F)
dat[[3]]<-fread("https://github.com/baselmans/multivariate_GWAMA/blob/master/Test_Data/NEU_100K_no23andMe.txt",showProgress=F,data.table=F)
dat[[4]]<-fread("https://github.com/baselmans/multivariate_GWAMA/blob/master/Test_Data/DEP_100K_no23andMe.txt",showProgress=F,data.table=F)

#Provide SNP heritabilities
h2list <- c(0.0498, 0.0441, 0.0708, 0.0294)
multivariate_GWAMA(x=dat,cov_Z=CTI,h2=h2list,out=".",name="test",output_gz=F,check_rows=F)
