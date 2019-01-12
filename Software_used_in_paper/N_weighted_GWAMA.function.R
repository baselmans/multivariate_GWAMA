Multivariate_GWAMA<-function(x,cov_Z=NULL,out=getwd(),name=NULL,save.matrix=F,genomeWideSignificant=5e-08,order.output=F,output.gz=T){
	
	if(substr(x=out,start=nchar(out),stop=nchar(out))!="/"){
		out<-paste0(out,"/")
	}
	
	out<-paste0(out,name)
	
	log.file<-file(paste0(out,".N_weighted_GWAMA.log"),"wt")
	#error.file<-file(paste0(out,".N_weighted_GWAMA.error"),"wt")
	sink(log.file,append=F,type="output",split=T)
	#sink(error.file,type="message")
	sink(log.file,append=T,type="message")
	
	begin.time<-Sys.time()
	cat("R function to perform N-weighted GWAMA of correlated traits\n")
	cat("For detailed information see: Baselmans et al. (2017). Multivariate Genome-wide and integrated transcriptome and epigenome-wide analyses of the Well-being spectrum\n\n")
	if(is.null(name)==F){
		cat("The following analysis was run: ",name,"\n",sep="")
	}
	cat("Analysis started at ",format(begin.time,trim="%Y-%b-%d %X",usetz=T),"\n",sep="")	
	
	cat("Performing basic sanity checks:\n")
	size.range<-range(unlist(lapply(x,nrow)))
	cat("  ",length(x)," inputs were supplied in x\n  Number of SNPs: ",format(size.range[1],big.mark=",")," - ",format(size.range[2],big.mark=","),"\n",sep="")
	#check length x vs dims cov_Z, and the symmetry of cov_Z
	if(is.null(cov_Z)==F){
		cat("  User supplied the error correlation between the inputs to cov_Z\n")
		if(is.matrix(cov_Z)==F){
			cat("  WARNING: cov_Z was NOT a matrix\n")
			if(is.numeric(as.matrix(cov_Z))==F){
				#close(log.file)
				stop("cov_Z is not-numerical\n",call.=F)
			}else{
				cat("  Successfully changed cov_Z into matrix")
				cov_z<-as.matrix(cov_Z)
			}
		}
		cat("  Checking for missingness in cov_Z...")
		if(sum(is.na(cov_Z))!=0){
			#cat("ERROR: cov_Z contains missing values")
			#close(log.file)
			stop("cov_Z contains missing values")
		}
		cat("done\n")
		cat("  Checking whether cov_Z is symmetric\n")
		if(nrow(cov_Z)!=ncol(cov_Z)){
			#close(log.file)
			stop(paste0("The number of rows (",nrow(cov_Z),") and number of columns (",ncol(cov_Z),") did not match\n"),call.=F)
		}else{
			cat("    Number of rows: ",nrow(cov_Z),"\n",sep="")
			cat("    Number of columns: ",ncol(cov_Z),"\n",sep="")
		}
		if(length(x)!=nrow(cov_Z)){
			#close(log.file)
			stop(paste0("The number of inputs in x (",length(x),") did not correspond to the dimension of cov_Z (r",nrow(cov_Z)," c",ncol(cov_Z),")\n"),call.=F)
		}
		cat("    The dimensions of cov_Z match with the number of inputs in x\n    Checking the values within cov_Z...")
		if(isSymmetric(cov_Z)==F){
			cat("WARNING: cov_Z was NOT symmetric\n    Checking whether removing the names solves the issue...")
			if(isSymmetric(unname(cov_Z))){
				cat("this did the trick (NOTE: the function continues with the named matrix, which may give funny results)\n")
			}else{
				cat("nope. ")
				#close(log.file)
				stop("The values in cov_Z were non-symmetric\n",call.=F)
			}
		}else{
			cat("done\n  cov_Z is symmetric\n")
		}
		#print(round(x=cov_Z,digits=3))
	}
	cat("Sanity checks finished without problems\n")
	
	cat("Setting unified headers...") 
	x<-lapply(x,setNames,nm=c("SNP","MarkerName","CHR","BP","A1","A2","EAF","N","Z","PVAL"))
	cat("done\n")
	
	cat("Removing duplicate entries...")
	x<-lapply(x,function(x){return(x[duplicated(x[,2])==F,])})
	cat("done\n")
	
	cat("Changing alleles to upper-case...")
	x<-lapply(x,function(x){
			dum<-cbind(x[,1:4],apply(x[,5:6],2,toupper),x[,7:ncol(x)])
			dum$A1<-as.character(dum$A1)
			dum$A2<-as.character(dum$A2)
			return(dum)
		}
	)
	cat("done\n")
	
	cat("Extracting the intersection of SNPs in x...")
	extract_SNPs<-function(x,column){
		SNPlist.all<-unique(unlist(lapply(x,"[[",column)))
		SNPlist.mat<-as.data.frame(matrix(0,nrow=length(SNPlist.all),ncol=length(x),dimnames=list(SNPlist.all,names(x))),stringsAsFactors=F)
		for(i in 1:ncol(SNPlist.mat)){
			SNPlist.mat[rownames(SNPlist.mat)%in%x[[i]][,column],i]<-1
		}
		SNPlist.mat$TOTAL<-rowSums(SNPlist.mat)
		return(SNPlist.mat)
	}
	SNPlist.mat<-extract_SNPs(x,2)
	SNPs<-data.frame(MarkerName=rownames(SNPlist.mat[SNPlist.mat$TOTAL==length(x),]),stringsAsFactors=F)
	cat("done\n",format(nrow(SNPs),big.mark=",")," SNPs were found in all files",sep="")
	#willen we dit ook opslaan (i.e. later komt dit nog een keer voor de SNPs die echt in het model gaan)?
	if(save.matrix){
		cat("\nAn overview is saved at ",out,".N_weighted_GWAMA.SNP_list.matrix",sep="")
		SNPlist.mat.out<-as.data.frame(cbind(rownames(SNPlist.mat),as.data.frame(SNPlist.mat)))
		colnames(SNPlist.mat.out)<-c("SNP",colnames(SNPlist.mat))
		if(output.gz){
			cat(".gz\n")
			gz1<-gzfile(paste0(out,".N_weighted_GWAMA.SNP_list.matrix.gz"),"w")
		}else{
			cat("\n")
			gz1<-paste0(out,".N_weighted_GWAMA.SNP_list.matrix")
		}
		write.table(x=SNPlist.mat.out,file=gz1,quote=F,sep="\t",row.names=F)
		if(is.character(gz1)==F){
			close(gz1)
		}
		rm(SNPlist.mat.out,gz1)
		gc()
	}
	cat("\nKeeping only common SNPs between the different inputs...")
	x_common<-lapply(x,function(x,y){return(x[match(y$MarkerName,x$MarkerName),])},y=SNPs)
	cat("done\n")
	rm(x,SNPs,size.range,SNPlist.mat)
	gc()
	
	#SANITY CHECKS
	cat("Checking whether alleles are aligned between the inputs...")
	A1A2_file1<-subset(x_common[[1]],select=c(A1,A2))
	colnames(A1A2_file1)[1]<-"A1.ref"
	colnames(A1A2_file1)[2]<-"A2.ref"
	cat("done (NOTE: The first input in x was used as reference file)\n")
	
	#Check whether the alleles of all data frames in your list are aligned properly
	identicalYesNo<-lapply(x_common,function(x,y){identical(x$A1,y$A1.ref) & identical(x$A2,y$A2.ref)},y=A1A2_file1)
	print(identicalYesNo)
	
	if(sum(unlist(identicalYesNo))!=length(x_common)){
		cat("Unaligned SNPs were found\n")
		#cbind with the "ref file"-> first data frame in list will be considered the ref allele
		x_reffile<-lapply(x_common,function(x){cbind(x,A1A2_file1)})
		rm(x_common,A1A2_file1)
		gc()
		
		#Check which alleles are aligned between files and which alleles are not aligned
		aligned<-lapply(x_reffile,function(x){subset(x,A1==A1.ref & A2==A2.ref)})
		
		flipped<-lapply(x_reffile,function(x){subset(x,A1==A2.ref & A2==A1.ref)})
		N_flipped<-sapply(flipped,nrow)
		cat("How many markers were flipped?\n")
		print(noquote(format(N_flipped,big.mark=",")))
		
		if(sum(unlist(N_flipped))!=0){
			cat("Flipping Z-scores for SNPs that are not aligned...")
			Zlist<-lapply(flipped,"[","Z")
			Zlist<-lapply(Zlist,function(x){x$Z * -1})
			#make data frame again
			Zlist<-lapply(Zlist,function(x){as.data.frame(x)})
			#cbind newbetas with not aligned file
			flipped_new_Z<-Map(cbind,flipped,Zlist) 
			#format right column names
			flipped_new_Z<-lapply(flipped_new_Z,setNames,nm=c("SNP","MarkerName","CHR","BP","A1","A2","EAF","N","ZnotUSE","PVAL","A1.ref","A2.ref","Z"))
			cat("done\n")

			cat("Flipping EAF for SNPs that are not aligned...")
			EAFlist<-lapply(flipped,"[","EAF")
			EAFlist<-lapply(EAFlist,function(x){1-x$EAF})
			EAFlist<-lapply(EAFlist,function(x){as.data.frame(x)})
			flipped_new_Z_EAF<-Map(cbind,flipped_new_Z,EAFlist) 
			flipped_new_Z_EAF<-lapply(flipped_new_Z_EAF,setNames,nm=c("SNP","MarkerName","CHR","BP","A1","A2","EAFnotUSE","N","ZnotUSE","PVAL","A1.ref","A2.ref","Z","EAF"))
			cat("done\n")

			rm(flipped_new_Z)
			gc()

			#select right columns for rbind with data frames in list that are not complete due to misalignment
			flipped_reorder<-lapply(flipped_new_Z_EAF,function(x){subset(x,select=(c(SNP,MarkerName,CHR,BP,A2,A1,EAF,N,Z,PVAL,A1.ref,A2.ref)))})
			flipped_reorder<-lapply(flipped_reorder,setNames,nm=c("SNP","MarkerName","CHR","BP","A1","A2","EAF","N","Z","PVAL","A1.ref","A2.ref"))

			rm(flipped_new_Z_EAF)
			gc()
			
			cat("Combining the newly aligned SNPs with SNPs that were already combined...")
			x_common_aligned<-Map(rbind,aligned,flipped_reorder)
			cat("done\n")

			rm(aligned,flipped_reorder,Zlist,EAFlist)
			gc()
		}
		rm(flipped,N_flipped)
		gc()
		
		if(exists("x_common_aligned")==F){
			x_common_aligned<-aligned
			rm(aligned)
			gc()
		}
		
		all_aligned<-lapply(x_common_aligned,function(x){subset(x,A1==A1.ref & A2==A2.ref)})
		rm(x_common_aligned)
		gc()
		
		unalignable<-lapply(x_reffile,function(x){subset(x,(A1!=A1.ref & A1!=A2.ref)|(A2!=A1.ref & A2!=A2.ref))})
		N_unalignable<-sapply(unalignable,nrow)
		if(sum(N_unalignable)!=0){
			cat("On rare occasions, SNPs cannot be aligned (e.g. A-C in the reference file vs A-T in the input). This might lead to summary statistics of different length\n")
			cat("How many markers could not be aligned?\n")
			print(noquote(format(N_unalignable,big.mark=",")))
			cat("These were removed from the data frame(s)\nA list was saved at ",out,".not_aligned",sep="")
			unalignable_SNPs<-unique(unlist(lapply(unalignable,"[[","MarkerName")))
			if(output.gz){
				cat(".gz\n")
				gz1<-gzfile(paste0(out,".not_aligned.gz"),"w")
			}else{
				cat("\n")
				gz1<-paste0(out,".not_aligned")
			}
			write.table(x=data.frame(MarkerName=unalignable_SNPs,stringsAsFactors=F),file=gz1,quote=F,sep="\t",row.names=F)
			if(is.character(gz1)==F){
				close(gz1)
			}
			rm(unalignable_SNPs,gz1)
			gc()
		}
		rm(unalignable,N_unalignable,x_reffile)
		gc()
		
		cat("Selecting only SNPs that are common again and extracting them from the data frames...")
		SNPlist.mat<-extract_SNPs(all_aligned,2)
		SNPs<-data.frame(MarkerName=rownames(SNPlist.mat[SNPlist.mat$TOTAL==length(all_aligned),]),stringsAsFactors=F)
		all_aligned<-lapply(all_aligned,function(x,y){x[match(y$MarkerName,x$MarkerName),]},y=SNPs)
		cat("done\n")
		rm(SNPlist.mat,SNPs)
		gc()
	}else{
		cat("All SNPs were aligned\n")
		all_aligned<-x_common
		rm(x_common,A1A2_file1)
		gc()
	}
	rm(identicalYesNo)
	gc()
	
	cat(format(nrow(all_aligned[[1]]),big.mark=",")," SNPs will be used in the N-weighted GWAMA\n",sep="")
	if(save.matrix){
		cat("A list was saved at ",out,".N_weighted_GWAMA.input",sep="")
		if(output.gz){
			cat(".gz\n")
			gz1<-gzfile(paste0(out,".N_weighted_GWAMA.input.gz"),"w")
		}else{
			cat("\n")
			gz1<-paste0(out,".N_weighted_GWAMA.input")
		}		
		write.table(x=all_aligned[[1]][,1],file=gz1,quote=F,sep="\t",row.names=F,col.names=c("MarkerName"))
		if(is.character(gz1)==F){
			close(gz1)
		}
		rm(gz1)
	}

	#####################################################################################
	# By now your files: 																#
	# *all contain common SNPs (present in all sumstats) 								#
	# *A1 and A2 are aligned with the each other using the first sumstats as reffile 	#
	# Now we are ready to perform the steps to perform the multivariate GWAMA 			#
	#####################################################################################
	
	cat("Extracting Z-scores from aligned summary statistics...")
	Zlist<-lapply(all_aligned,"[","Z")
	cat("done\n")
	
	cat("Selecting weights and calculating their square roots...")
	h2list <- list(0.0498, 0.0441, 0.0748, 0.0305) #Full
	Weightlist<-lapply(all_aligned,"[","N")
	sqrtWeightlist<-lapply(Weightlist,function(x){sqrt(x)})
	for(i in 1:length(sqrtWeightlist)){
		sqrtWeightlist[[i]]<-sqrtWeightlist[[i]] * sqrt(h2list[[i]])
	}
	Wnumeric<-lapply(sqrtWeightlist,function(x){as.numeric(x$N)})
	Wnumeric_rbind=do.call(rbind,Wnumeric)
	cat("done\n")
	
	cat("Calculating product Z-scores and sqrt-weights per input, and summing them across the different summary statistics...")
	ZxW<-Map('*',Zlist,sqrtWeightlist)
	ZxW<-do.call(cbind,ZxW)
	ZxWtotal<-as.data.frame(rowSums(ZxW,na.rm=FALSE,dims=1))
	colnames(ZxWtotal)[1]<-"Ztotal"
	cat("done\n")

	Z_all_m<-matrix(NA,nrow=nrow(Weightlist[[1]]),1)
	if(is.null(cov_Z)){
		cat("WARNING: Error correlation was NOT pre-specified using cov_Z\nThe error correlation will be calculated by creating a correlation matrix of the Z-scores across the different inputs (NOTE: this is more conservative)\n")
		#create correlation matrix
		Zlist1<-do.call(cbind,Zlist)
		cov_Z<-cor(Zlist1)
		colnames(cov_Z)<-rownames(cov_Z)<-names(all_aligned)
		#print(round(x=cov_Z,digits=3))
		cat("The correlation matrix was saved at ",out,".corZ",sep="")
		if(output.gz){
			cat(".gz\n")
			gz1<-gzfile(paste0(out,".corZ.gz"),"w")
		}else{
			cat("\n")
			gz1<-paste0(out,".corZ")
		}
		write.table(x=cov_Z,file=gz1,quote=F,sep="\t")
		if(is.character(gz1)==F){
			close(gz1)
		}
	}
	cat("Computing weighted Z-scores...")
	for(i in 1:ncol(Wnumeric_rbind)){
		W<-diag(Wnumeric_rbind[,i])
		# hoe doen we dit met de correlatie erin??
		Z_all_m[i,]<-(ZxWtotal[i,]) / sqrt(sum(W%*%cov_Z%*%W))
	}
	cat("done\n")
	
	cat("Computing P values from weighted Z scores...")
	PVAL<-pchisq((Z_all_m)^2,1,lower.tail=F)
	cat("done\n")
	cat(format(length(PVAL[PVAL<=genomeWideSignificant]),big.mark=",")," SNPs reached genome-wide significance (Minimum P value: ",sprintf("%0.4e",min(PVAL)),")\n",sep="")
	
	cat("Computing the N for the N-weighted GWAMA...")
	N<-do.call(cbind,Weightlist)
	Ntotal<-as.data.frame(rowSums(N,na.rm=FALSE,dims=1))
	colnames(Ntotal)[1]<-"N"
	cat("done\n")

	cat("Computing the combined EAF...")
	#extract EAF and N
	EAFlist<-lapply(all_aligned,"[",c("EAF","N"))
	#multiply EAF by sample size
	EAFxN<-lapply(EAFlist,function(x){x$N * x$EAF})
	# cbind alle element of list
	EAFxN<-do.call(cbind,EAFxN)
	#row sum to obtain EAF weighted by sample size per SNP
	EAFtot<-as.data.frame(rowSums(EAFxN,na.rm=FALSE,dims=1))
	colnames(EAFtot)[1]<-"EAFtot"
	#compute multi EAF by divinding EAFtotal by total sample size
	EAFmulti<-EAFtot / Ntotal
	cat("done\n")
	
	#########################################################################################
	### NOW WE HAVE EVERYTHING THAT WE WILL NEED TO ASSEMBLE THE NEW MULTIVARIATE ###########
	### GWAMA SUMMARY STATISTIC #############################################################
	#########################################################################################

	cat("Constructing N-weighted multivariate summary statistics\n")

	#cat("  Extracting data frames from list...")
	#for(i in seq(all_aligned)){
	#	assign(paste("all_aligned",i,sep=""),all_aligned[[i]])
	#}
	#cat("done\n")
	
	#MULTI_Nweighted<-subset(all_aligned1,select=(c(SNP,MarkerName,CHR,BP,A1,A2)))
	MULTI_Nweighted<-subset(all_aligned[[1]],select=(c(SNP,MarkerName,CHR,BP,A1,A2)))
	#MULTI_Nweighted[,5]=toupper(MULTI_Nweighted[,5])
	#MULTI_Nweighted[,6]=toupper(MULTI_Nweighted[,6])

	MULTI_Nweighted<-cbind(MULTI_Nweighted,Z_all_m,PVAL,Ntotal,EAFmulti)

	#change names header
	colnames(MULTI_Nweighted)[7]<-"Z"
	colnames(MULTI_Nweighted)[10]<-"EAF"

	cat("  Computing BETA and SE...")
	#formula $Z/sqrt($N*2*$eaf*(1-$eaf));
	MULTI_Nweighted$BETA<-MULTI_Nweighted$Z / sqrt(MULTI_Nweighted$N*2*MULTI_Nweighted$EAF*(1-MULTI_Nweighted$EAF))
	#Compute SE
	MULTI_Nweighted$SE<-((1/sqrt(MULTI_Nweighted$N)) * (1 /sqrt(2*(MULTI_Nweighted$EAF*(1-MULTI_Nweighted$EAF)))))
	cat("done\n")
	
	#order the data
	if(order.output){
		cat("  Re-ordering the data on chromosome and base-pair position...")
		MULTI_Nweighted<-MULTI_Nweighted[order(MULTI_Nweighted$SNP),]
		cat("done\n")
	}
	
	#Write file to directory
	cat("Successfully created summary statistics for ",format(nrow(MULTI_Nweighted),big.mark=",")," SNPs\nThe output was saved to ",out,".N_weighted_GWAMA.results",sep="")
	if(output.gz){
		cat(".gz")
		gz1<-gzfile(paste0(out,".N_weighted_GWAMA.results.gz"),"w")
	}else{
		gz1<-paste0(out,".N_weighted_GWAMA.results")
	}
	write.table(MULTI_Nweighted,file=gz1,sep="\t",quote=F,row.names=F)
	if(is.character(gz1)==F){
		close(gz1)
	}
	
	end.time<-Sys.time()
	cat("\nAnalysis successfully finished at ",format(x=end.time,trim="%Y-%b-%d %X",usetz=T),"\n",sep="")
	analysis.duration<-difftime(time1=end.time,time2=begin.time,units="secs")
	analysis.mins<-floor(x=floor(x=(analysis.duration/60)))
	analysis.secs<-analysis.duration-analysis.mins*60
	cat("Analysis took ",analysis.mins," minute(s) and ",analysis.secs," seconds\n",sep="")
	
	closeAllConnections()
}

manhattan1<-function(x,chr="CHR",bp="BP",p="P",snp="SNP",col=c("gray10","gray60"),chrlabs=NULL,suggestiveline=-log10(1e-05),genomewideline=-log10(5e-08),highlight1=NULL,highlight2=NULL,highlight3=NULL,logp=TRUE,...){
	CHR=BP=P=index=NULL
	if(!(chr %in% names(x))){
		stop(paste("Column",chr,"not found!"))
	}
	if(!(bp %in% names(x))){
		stop(paste("Column",bp,"not found!"))
	}
	if(!(p %in% names(x))){
		stop(paste("Column",p,"not found!"))	
	}
	if(!(snp %in% names(x))){
		warning(paste("No SNP column found. OK unless you're trying to highlight."))
	}
	if(!is.numeric(x[[chr]])){
		stop(paste(chr,"column should be numeric. Do you have 'X','Y','MT',etc? If so change to numbers and try again."))	
	}
	if(!is.numeric(x[[bp]])){
		stop(paste(bp,"column should be numeric."))
	}
	if(!is.numeric(x[[p]])){
		stop(paste(p,"column should be numeric."))
	}
	d<-data.frame(CHR=x[[chr]],BP=x[[bp]],P=x[[p]])
	if(!is.null(x[[snp]])){
		d<-transform(d,SNP=x[[snp]])
	}
	d<-subset(d,(is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
	d<-d[order(d$CHR,d$BP),]
	if(logp){
		d$logp<--log10(d$P)
	}else{
		d$logp<-d$P
	}
	d$pos<-NA
	d$index<-NA
	ind=0
	for(i in unique(d$CHR)){
		ind=ind + 1
		d[d$CHR==i,]$index<-ind
	}
	nchr=length(unique(d$CHR))
	if(nchr==1){
		options(scipen=999)
		d$pos<-d$BP/1e+06
		ticks<-floor(length(d$pos))/2 + 1
		xlabel<-paste("Chromosome",unique(d$CHR),"position(Mb)")
		labs<-ticks
	}else{
		lastbase<-0
		ticks<-NULL
		for(i in unique(d$index)){
			if(i==1){
				d[d$index==i,]$pos<-d[d$index==i,]$BP
			}else{
				lastbase<-lastbase + tail(subset(d,index==i - 1)$BP,1)
				d[d$index==i,]$pos<-d[d$index==i,]$BP + lastbase
			}
			ticks<-c(ticks,(min(d[d$CHR==i,]$pos) + max(d[d$CHR==i,]$pos))/2 + 1)
		}
		xlabel="Chromosome"
		labs<-unique(d$CHR)
	}
	xmax<-ceiling(max(d$pos) * 1.03)
	xmin<-floor(max(d$pos) * -0.03)
	def_args<-list(xaxt="n",bty="n",xaxs="i",yaxs="i",las=1,pch=20,xlim=c(xmin,xmax),ylim=c(0,ceiling(max(d$logp))),xlab=xlabel,ylab=expression(-log10))
	dotargs<-list(...)
	do.call("plot",c(NA,dotargs,def_args[!names(def_args) %in% names(dotargs)]))
	if(!is.null(chrlabs)){
		if(is.character(chrlabs)){
			if(length(chrlabs)==length(labs)){
				labs<-chrlabs
			}else{
				warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
			}
		}else{
			warning("If you're trying to specify chromosome labels,chrlabs must be a character vector")
		}
	}
	if(nchr==1){
		axis(1,...)
	}else{
		axis(1,at=ticks,labels=labs,...)
	}
	col<-rep(col,max(d$CHR))
	if(nchr==1){
		with(d,points(pos,logp,pch=20,col=col[1],...))
	}else{
		icol<-1
		for(i in unique(d$index)){
			with(d[d$index==unique(d$index)[i],],points(pos,logp,col=col[icol],pch=20,...))
			icol<-icol + 1
		}
	}
	if(suggestiveline){
		abline(h=suggestiveline,col="blue")
	}
	if(genomewideline){
		abline(h=genomewideline,col="red")
	}
	if(!is.null(highlight1)){
		if(any(!(highlight1 %in% d$SNP))){
			warning("You're trying to highlight1 SNPs that don't exist in your results.")
		}
		d.highlight1<-d[which(d$SNP %in% highlight1),]
		with(d.highlight1,points(pos,logp,col="sienna1",pch=17,cex=1.2,...))
	}
	if(!is.null(highlight2)){
		if(any(!(highlight2 %in% d$SNP))){
			warning("You're trying to highlight2 SNPs that don't exist in your results.")
		}
		d.highlight2<-d[which(d$SNP %in% highlight2),]
		with(d.highlight2,points(pos,logp,col="blue4",pch=20,cex=2.5,...))
	}
	if(!is.null(highlight3)){
		if(any(!(highlight3 %in% d$SNP))){
			warning("You're trying to highlight3 SNPs that don't exist in your results.")
		}
		d.highlight3<-d[which(d$SNP %in% highlight3),]
		with(d.highlight3,points(pos,logp,col="mediumblue",pch=20,cex=2.5,...))
	}
}
