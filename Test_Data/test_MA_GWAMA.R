library(data.table)
library(sp)
library(raster)
library(RcppArmadillo)
library(unmarked)
library(VGAM)
library(AICcmodavg)
library(metafor)


#read in the data
T1<-fread("https://github.com/baselmans/multivariate_GWAMA/blob/master/Test_Data/LS_1K_no23andMe.txt?raw=TRUE",showProgress=F,data.table=F)
T2<-fread("https://github.com/baselmans/multivariate_GWAMA/blob/master/Test_Data/PA_1K_no23andMe.txt?raw=TRUE",showProgress=F,data.table=F)
T3<-fread("https://github.com/baselmans/multivariate_GWAMA/blob/master/Test_Data/NEU_1K_no23andMe.txt?raw=TRUE",showProgress=F,data.table=F)
T4<-fread("https://github.com/baselmans/multivariate_GWAMA/blob/master/Test_Data/DEP_1K_no23andMe.txt?raw=TRUE",showProgress=F,data.table=F)
### If needed create Beta or SE
#Beta
T1$b1 <- T1$z/sqrt(T1$n*2*T1$eaf*(1-T1$eaf))
T2$b2 <- T2$z/sqrt(T2$n*2*T2$eaf*(1-T2$eaf))
T3$b3 <- T3$z/sqrt(T3$n*2*T3$eaf*(1-T3$eaf))
T4$b4 <- T4$z/sqrt(T4$n*2*T4$eaf*(1-T4$eaf))

#SE
T1$se1 <- ((1/sqrt(T1$n)) * (1/sqrt(2*(T1$eaf*(1-T1$eaf)))))
T2$se2 <- ((1/sqrt(T2$n)) * (1/sqrt(2*(T2$eaf*(1-T2$eaf)))))
T3$se3 <- ((1/sqrt(T3$n)) * (1/sqrt(2*(T3$eaf*(1-T3$eaf)))))
T4$se4 <- ((1/sqrt(T4$n)) * (1/sqrt(2*(T4$eaf*(1-T4$eaf)))))

#select right columns (rs, beta, se)

T1_b1_se1 <- T1[c(2,11:12)]
T2_b2_se2 <- T2[c(2,11:12)]
T3_b3_se3 <- T3[c(2,11:12)]
T4_b4_se4 <- T4[c(2,11:12)]

#Merge data

M1 <- merge(T1_b1_se1,T2_b2_se2,by=1)
M1 <- merge(M1,T3_b3_se3,by=1)
M1 <- merge(M1,T4_b4_se4,by=1)

#remove missings
M1 <- na.omit(M1)
#M1 <- M1[1:1000,]

A1 <- merge(T1,M1, by="rs")
A1 <- A1[row.names(unique(A1["rs"])),]

B <- M1[,c(2,4,6,8)]
SE <- M1[,c(3,5,7,9)]


# cross trait intercept
cov_Z <- as.matrix(read.table("https://github.com/baselmans/multivariate_GWAMA/blob/master/Test_Data/cross_trait_intercept.txt?raw=TRUE", header =T))

# heritabilitie estimates
#h2 = vector sqrt(h2)
h2 <- as.vector(c(sqrt(0.0498),sqrt(0.0441), sqrt(0.0748), sqrt(0.0305)))


output_SNP <- matrix(NA, ncol=8,nrow=nrow(B))

for(i in 1:nrow(B)){
  
  V <- diag(SE[i,]) %*% cov_Z %*% diag(SE[i,])
  
  
  yi <- as.numeric(B[i,])
  
  
  grid <- expand.grid(c(0,1),c(0,1),c(0,1),c(0,1))
  
  grid <- t(grid)
  
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

Output <- cbind(A1[,c(1,2,3,4,5,6)],output_SNP,pchisq((output_SNP[,1]/output_SNP[,2])^2,1,lower=F),
                pchisq((output_SNP[,3]/output_SNP[,4])^2,1,lower=F),pchisq((output_SNP[,5]/output_SNP[,6])^2,
                1,lower=F),pchisq((output_SNP[,7]/output_SNP[,8])^2,1,lower=F))

names(Output) <- c("RS","cptid","CHR","BP","A1","A2","Beta_LS","SE_LS","Beta_PA","SE_PA","Beta_NEU","SE_NEU","Beta_DEP","SE_DEP","P_LS","P_PA","P_NEU","P_DEP")

write.table("AVGWAS.txt" ,quote=F,row.names=F,col.names=T)

                                                                                                                                                                                                                     names(Output) <- c("RS","cptid","CHR","BP","A1","A2","Beta_LS","SE_LS","Beta_PA","SE_PA","Beta_NEU","SE_NEU","Beta_DEP","SE_DEP","P_LS","P_PA","P_NEU","P_DEP")
                                                                                                                                                                                                                   write.table(Output,paste0("AVGWAS_", chr,".txt"),quote=F,row.names=F,col.names=T)
                                                                                                                                                                                                                   
                                                                                                                                                                                                                   



