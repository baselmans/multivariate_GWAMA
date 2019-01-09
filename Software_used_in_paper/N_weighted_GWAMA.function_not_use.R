##########################################################################
##																		##
##			R function to run a multivariate N-weighted-GWAMA			##
##			Version 1.2.2 (last modified on 12-OCT-2018)
##      Created by Hill Ip & Bart Baselmans           ##
##																		##
##########################################################################

multivariate_GWAMA<-function(	x,cov_Z=NULL,h2,
								out,name,order_output=F,output_gz=T,
								covars=NULL,add_intercept=F,save_coefficients=F,
								check_rows=T,genome_wide_significant=5e-08,type=NA,
								n_cases=NULL,n_controls=NULL,pop_prev=NULL,rg=NULL,h2_dich=NULL,h2_cont=NULL,lj=124.718,M=5961159
							){
	if(missing(out)){
		warning_out_missing<-T
		out<-paste0(getwd(),"/")
	}else{
		if(!is.character(out)){
			stop(paste0("\"out\" (",out,") was not a character\n"))
		}else if(out=="."){
			warning_out_dot<-T
			out<-paste0(getwd(),"/")
		}else if(substr(out,nchar(out),nchar(out))=="/"){
			if(!dir.exists(out)){
				stop(paste0("out=\"",out,"\" does not exist"))
			}
		}else if(substr(out,nchar(out),nchar(out))!="."){
			if(dir.exists(out)){
				warning_out_missing_slash<-T
				out<-paste0(out,"/")
			}else{
				stop(paste0("invalid value in out=\"",out,"\"\nNOTE, either:\n  1) separately supply the output directory via out (e.g. out=\"/desired/output/directory/\" [note the trailing slash!]) and prefix via name (e.g. name=\"prefix\"); OR\n  2) supply a single argument to out, ending with a dot (e.g. out=\"/desired/output/directory/prefix.\")"))
			}
		}
	}
	out_dir<-out

	if(missing(name)){
		warning_name_missing<-T
		if(substr(out,nchar(out),nchar(out))=="."){
			out_org<-out
			out<-get_outDir_name(out)
			out_dir<-out[[1]]
			name<-out[[2]]
		}else{
			error_name_not_found<-T
			name<-""
		}
	}else if(!exists("warning_out_missing")){
		if(substr(name,nchar(name),nchar(name))=="."){
			name<-gsub("\\.$","",name)
		}
		if(substr(out,nchar(out),nchar(out))=="."){
			warning_out_name_both_supplied<-T
			out_org<-out
			name_supplied<-name
			out<-get_outDir_name(out)
			out_dir<-out[[1]]
			name<-out[[2]]
			if(name==name_supplied){
				rm(warning_out_name_both_supplied)
			}
		}
	}
	out<-paste0(out_dir,name,".N_weighted_GWAMA.")

	if(file.exists(paste0(out,"log"))){
		warning_log_exists<-T
		file.rename(paste0(out,"log"),paste0(out,"OLD.log"))
	}
	log.file<-file(paste0(out,"log"),"wt")
	sink(log.file,,"output",T)
	sink(log.file,T,"message")

	ver<-"1.2.2"
	begin.time<-Sys.time()

	cat("-----\n\nR-function (v",ver,") for performing an N-weighted GWAMA over a set of genetically correlated traits\nFor detailed information, see: Baselmans et al. (2017). Multivariate genome-wide and integrated transcriptome and epigenome-wide analyses of the well-being spectrum\n\n-----\n\nAnalysis started at ",format(begin.time,trim="%y-%b-%d %X",usetz=T),"\n",sep="")
	if(exists("warning_out_missing")){
		cat("WARNING: out was not supplied. Output was saved at the working directory\n")
	}else{
		if(exists("warning_out_dot")){
			cat("WARNING: out was a single dot \".\". Assuming user wants to save output in the working directory\n")
		}else if(exists("warning_out_missing_slash")){
			cat("WARNING: a slash was appended to the end of out\n  If out contains the prefix, please end the argument with a dot (e.g. \"/desired/output/directory/prefix.\")\n")
		}
	}
	if(exists("error_name_not_found")){
		stop("user did not specify the output name\nNOTE, either:\n  1) separately supply the output directory via out (e.g. out=\"/desired/output/directory/\" [note the trailing slash!]) and prefix via name (e.g. name=\"prefix\"); OR\n  2) supply a directory and prefix to out by ending it with a dot (e.g. out=\"/desired/output/directory/prefix.\")",call.=F)
	}
	if(exists("warning_name_missing")){
		cat("WARNING: name was not supplied. Output prefix was guessed via out=\"",out_org,"\"\n",sep="")
	}else if(exists("warning_out_name_both_supplied")){
		cat("WARNING: out contained a prefix (\"",out_org,"\") and name was supplied (\"",name_supplied,"\"). name=\"",name_supplied,"\" was ignored\n",sep="")
	}
	cat("The output was saved as:",gsub(".N_weighted_GWAMA.",".",out),"\n")
	if(exists("warning_log_exists")){
		cat("WARNING: \"",out,"log\" already existed and was renamed to \"",out,"OLD.log\"\n",sep="")
	}
	if(file.exists(paste0(out,"results.txt"))){
		cat("WARNING: \"",out,"results.txt\" already existed and was renamed to \"",out,"OLD.results.txt\"\n",sep="")
		file.rename(paste0(out,"results.txt"),paste0(out,"OLD.results.txt"))
	}else if(file.exists(paste0(out,"results.txt.gz"))){
		cat("WARNING: \"",out,"results.txt.gz\" already existed and was renamed to \"",out,"OLD.results.txt.gz\"\n",sep="")
		file.rename(paste0(out,"results.txt.gz"),paste0(out,"OLD.results.txt.gz"))
	}
	output_gz<-check_logical(output_gz,"output_gz",T)

	cat("\n-----\n\nBasic sanity checks:\n\nChecking the input in x...")
	if(!is.list(x)){
		stop("x was not a list",call.=F)
	}
	if(length(x)<2){
		stop(paste("x had length",length(x)),call.=F)
	}
	cat("done\n  x contained",length(x),"entries\n")

	if(is.null(cov_Z)){
		cat("WARNING: user did NOT supply error correlation via cov_Z\n  The error correlation was calculated by taking the correlation between the Z-scores across the different input files (NOTE: this is a conservative approach)\n")
	}else{
		cat("User supplied the error correlation between the inputs via cov_Z\n  Checking whether cov_Z was a matrix...")
		if(is.matrix(cov_Z)){
			cat("done\n")
		}else{
			cat("WARNING: cov_Z was not a matrix\n    Attempting to turn cov_Z into a matrix...")
			if(!is.numeric(as.matrix(cov_Z))){
				stop("cov_Z contained non-numeric entries",call.=F)
			}
			cov_Z<-as.matrix(cov_Z)
			cat("done\n")
		}
		cat("  Checking whether cov_Z was numeric...")
		if(is.numeric(cov_Z)){
			cat("done\n")
		}else{
			cat("WARNING: cov_Z contained non-numeric values\n    Checking whether as.numeric() solved the issue...")
			dum<-matrix(suppressWarnings(as.numeric(cov_Z)),nrow(cov_Z),ncol(cov_Z))
			if(sum(is.na(dum))==length(cov_Z)){
				stop("cov_Z consisted of non-numeric values",call.=F)
			}else{
				cat("done\n")
				cov_Z<-dum
			}
		}
		cat("  Checking whether cov_Z was symmetric...")
		if(isSymmetric(cov_Z)){
			cat("done\n")
		}else{
			cat("WARNING: cov_Z was not symmetric\n    Checking whether removing the names solves the problem...")
			if(!isSymmetric(unname(cov_Z))){
				stop("the values inside cov_Z were non-symmetric",call.=F)
			}
			cat("done\n      NOTE: the function continued with the named matrix, which may return funny results\n")
		}
		cat("  Checking whether the dimensions of cov_Z matched with the number of inputs in x...")
		if(nrow(cov_Z)!=length(x)){
			stop("the dimensions of cov_Z (",nrow(cov_Z),") did not match with the number of inputs in x (",length(x),")",call.=F)
		}
		cat("done\n  Checking for missingness in cov_Z...")
		if(sum(is.na(cov_Z))!=0){
			stop("cov_Z contained missing values",call.=F)
		}
		cat("done\n")
	}
	cat("Checking h2...")
	if(missing(h2)){
		cat("WARNING: user did not supply a vector of SNP-heritabilities\n  Using a vector of 1's instead\n    NOTE: the multivariate Z-scores will be biased towards the input with the highest SNP-heritability\n")
		h2<-rep(1,length(x))
	}else{
		if(!is.vector(h2)){
			stop("h2 must be a vector",call.=F)
		}else if(length(h2)!=length(x)){
			if(length(h2)==1){
				cat("WARNING: user supplied a single value for h2 (",h2,"); the same value was used for all inputs\n",sep="")
				h2<-rep(h2,length(x))
			}else{
				stop(paste0("length of h2 (",length(h2),") differed from length of x (",length(x),")"),call.=F)
			}
		}else{
			cat("done\n  User supplied the following vector of SNP-heritabilities:\n",paste0("    ",h2),sep="\n")
		}
	}
	if(sum(h2<=0)!=0){
		stop(paste0("SNP-heritabilities out of bounds (range: ",min(h2)," - ",max(h2),")"),call.=F)
	}else if(sum(h2>1)!=0){
		cat("WARNING: the vector of SNP-heritabilities contained values larger than 1 (range: ",min(h2)," - ",max(h2),")",sep="")
	}

	if(!is.na(type)){
		type<-tolower(type)
		if(type!="demontis"){
			stop("only type=\"demontis\" is supported at the moment",call.=F)
		}
		if(length(x)!=2){
			stop(paste0("x must have length 2 for type=\"demontis\" [length(x)=",length(x),"]"),call.=F)
		}
	}

	cat("\nSanity checks finished without problems\n\n-----\n\nProcessing files before the GWAMA:\n  Checking the number of columns...")
	if(sum(unlist(lapply(x,ncol))==9)!=length(x)){
		cat("\n")
		print(lapply(x,ncol))
		stop("input files must contain exactly 9 columns",call.=F)
	}
	cat("done\n  Checking column order...")
	for(i in 1:length(x)){
		foo<-x[[i]]
		dum<-gsub("\\.","_",gsub("-","_",tolower(colnames(foo))))

		if(sum(c(
			"cptid","markername","marker",
			"snp","snps","snpid","snp_id",
			"snpids","snp_ids",
			"rs","rsmid","rs_mid",
			"rsid","rs_id",
			"rsnum","rs_num",
			"rsnums","rs_nums",
			"rsnumber","rs_number",
			"rsnumbers","rs_numbers"
		)%in%dum[1])!=1){
			stop(paste("The first column of file",i,"did not appear to contain SNPIDs"),call.=F)
		}

		if(sum(c(
			"chr","chrom","chroms",
			"chromosome","chromosomes"
		)%in%dum[2])!=1){
			stop(paste("The second column of file",i,"did not appear to contain the chromosomes"),call.=F)
		}

		if(sum(c(
			"bp","basepair","basepairs","base_pair",
			"pos","position","positions"
		)%in%dum[3])!=1){
			stop(paste("The third column of file",i,"did not appear to contain base-pair positions"),call.=F)
		}

		if(sum(c(
			"a0","allele0","allele_0",
			"a1","allele1","allele_1",
			"a2","allele2","allele_2",
			"ea","effectallele","effect_allele",
			"oa","otherallele","other_allele",
			"nea","non_effect_allele",
			"ref","reference_allele",
			"alt","alternative_allele",
			"codedall","coded_all",
			"codedallele","coded_allele",
			"noncodedall","non_codedall","noncoded_all","non_coded_all",
			"noncodedallele","non_codedallele","noncoded_allele","non_coded_allele"
		)%in%dum[4:5])!=2){
			stop(paste("Columns four and/or five of file",i,"did not appear to be alleles"),call.=F)
		}

		if(sum(c(
			"eaf","effect_allele_frequency",
			"oaf","other_allele_frequency",
			"frq","freq","frequency",
			"a0_frq","a0_freq","a0_frequency",
			"a1_frq","a1_freq","a1_frequency",
			"a2_frq","a2_freq","a2_frequency",
			"ref_frq","ref_freq","ref_frequency",
			"reference_frq","reference_freq","reference_frequency",
			"alt_frq","alt_freq","alt_frequency",
			"alternative_frq","alternative_freq","alternative_frequency"
		)%in%dum[6])!=1){
			stop(paste("The sixth column of file",i,"did not appear to contain effect allele frequencies"),call.=F)
		}

		if(sum(c(
			"n","nchrobs","n_total",
			"w","weight","weights"
		)%in%dum[7])!=1){
			stop(paste("The seventh column of file",i,"did not appear to contain sample sizes"),call.=F)
		}

		if(sum(c(
			"z","zscore","z_score",
			"zscores","z_scores"
		)%in%dum[8])!=1){
			stop(paste("The eighth column of file",i,"did not appear to contain Z-scores"),call.=F)
		}

		if(sum(c(
			"p","pval","p_val",
			"pvals","p_vals",
			"pvalue","p_value",
			"pvalues","p_values"
		)%in%dum[9])!=1){
			stop(paste("The nineth column of file",i,"did not appear to contain P-values"),call.=F)
		}

		check_rows<-check_logical(check_rows,"check_rows",T)
		if(check_rows){
			bar<-foo[sample(nrow(foo),max(5000,min(10000,0.001*nrow(foo)))),]

			chr_table<-tolower(names(table(bar[,2])))
			chr_table<-gsub("^chr","",chr_table)
			if(sum(chr_table%in%c(1:25,"x","y","m","mt"))!=length(chr_table)){
				stop(paste("Values in chromosome-column for file",i,"out of bounds"),call.=F)
			}

			if(sum(c("A","C","G","T")%in%bar[,4])!=4 & sum(c("a","c","g","t")%in%bar[,4])!=4){
				stop(paste("one or more of the alleles did not occur in effect-allele-column (column 4) or file",i,"which is very unlikely\n  NOTE: set check_rows to \"FALSE\" if this is known to be true in your dataset"),call.=F)
			}
			if(sum(c("A","C","G","T")%in%bar[,5])!=4 & sum(c("a","c","g","t")%in%bar[,5])!=4){
				stop(paste("one or more of the alleles did not occur in other-allele-column (column 5) or file",i,"which is very unlikely\n  NOTE: set check_rows to \"FALSE\" if this is known to be true in your dataset"),call.=F)
			}

			if(min(foo[,6],na.rm=T)<0|max(foo[,6],na.rm=T)>1){
				stop(paste0("EAF (column 6) out of bounds for file ",i," (",min(foo[,6],na.rm=T)," - ",max(foo[,6],na.rm=T),")"),call.=F)
			}

			if(min(foo[,9],na.rm=T)<0|max(foo[,9],na.rm=T)>1){
				stop(paste0("P-values (column 9) out of bounds for file ",i," (",min(foo[,9],na.rm=T)," - ",max(foo[,9],na.rm=T),")"),call.=F)
			}
		}
		rm(foo,dum)
	}
	cat("done\n")
	cat("  Input contained the following headers:\n\n")
	print(lapply(x,colnames))
	cat("  Setting unified headers...")
	x<-lapply(x,setNames,nm=c("SNPID","CHR","BP","EA","OA","EAF","N","Z","PVAL"))
	cat("done\n")
	cat("  Removing all instances of duplicated SNPIDs...")
	duplicated_snpids<-lapply(x,function(x){
		y<-unique(x[duplicated(x$SNPID),]$SNPID)
		return(y)
	})
	cat("done\n")
	cat("  Removing rows with missing data...")
	x<-lapply(x,function(x){
		y<-x[rowSums(is.na(x))==0,]
		return(y)
	})
	cat("done\n")
	cat("  Changing alleles to upper-case...")
	x<-lapply(x,function(x){
		y<-cbind(x[,1:3],apply(x[,4:5],2,toupper),x[,6:9])
		y$EA<-as.character(y$EA)
		y$OA<-as.character(y$OA)
		return(y)
	})
	cat("done\n  Checking whether alleles are aligned between the input files\n    NOTE: first occurrence of the SNP was used as reference\n    NOTE2: function does NOT take strand flips into account\n")
	reffile<-x[[1]][,c("SNPID","EA","OA")]
	for(i in 2:length(x)){
		dum<-x[[i]][,c("SNPID","EA","OA")]
		dum<-dum[!(dum$SNPID%in%reffile$SNPID),]
		reffile<-rbind(reffile,dum)
		rm(dum)
	}
	colnames(reffile)[2:3]<-c("EA.ref","OA.ref")
	x_common<-lapply(x,function(x,y=reffile){
		dum<-y[match(x$SNPID,y$SNPID),2:3]
		return(cbind(x,dum))
	})
	rm(x,reffile)

	identical_yes_no<-lapply(x_common,function(x){
		return(sum(identical(x$EA,x$EA.ref),identical(x$OA,x$OA.ref)))
	})
	if(sum(unlist(identical_yes_no))==(2*length(x_common))){
		cat("  Alleles were aligned across all input files\n")
		all_aligned<-x_common
		rm(x_common)
	}else{
		cat("  WARNING: unaligned SNPs were found\n    How many of these were flipped?\n")
		aligned<-lapply(x_common,function(x){
			return(x[x$EA==x$EA.ref & x$OA==x$OA.ref,])
		})
		flipped<-lapply(x_common,function(x){
			return(x[x$EA==x$OA.ref & x$OA==x$EA.ref,])
		})
		N_flipped<-sapply(flipped,nrow)
		print(noquote(format(N_flipped,big.mark=",")))
		if(sum(unlist(N_flipped))==0){
			all_aligned<-aligned
		}else{
			cat("    Flipping the Z-scores...")
			flipped<-lapply(flipped,function(x){
				x$Z<--x$Z
				return(x)
			})
			cat("done\n    Flipping the EAFs...")
			flipped<-lapply(flipped,function(x){
				x$EAF<-1-x$EAF
				return(x)
			})
			cat("done\n    Flipping the alleles...")
			flipped<-lapply(flipped,function(x){
				x<-x[,c("SNPID","CHR","BP","OA","EA","EAF","N","Z","PVAL","EA.ref","OA.ref")]
				colnames(x)<-c("SNPID","CHR","BP","EA","OA","EAF","N","Z","PVAL","EA.ref","OA.ref")
				return(x)
			})
			cat("done\n    Combining the newly aligned SNPs with the SNPs that were already aligned...")
			all_aligned<-Map(rbind,aligned,flipped)
			cat("done\n\n")
		}
		rm(aligned,flipped,N_flipped)

		allele_mismatch<-lapply(x_common,function(x){
			return(x[!((x$EA==x$EA.ref & x$OA==x$OA.ref)|(x$EA==x$OA.ref & x$OA==x$EA.ref)),]$SNPID)
		})
		N_allele_mismatch<-sapply(allele_mismatch,length)
		if(sum(unlist(N_allele_mismatch))!=0){
			cat("    On rare occasions, SNPs cannot be aligned (e.g. A-C in the reference file vs A-G in the input)\n    Number of unalignable SNPS:\n")
			print(noquote(format(N_allele_mismatch,big.mark=",")))
			cat("    These were removed from the input files\n    A list of SNPIDs was saved at ",out,"not_aligned.txt[.gz]\n\n",sep="")
			allele_mismatch<-lapply(allele_mismatch,function(x){
				return(data.frame(ID=x,SNPID=x,stringsAsFactors=F))
			})
			allele_mismatch_snps<-merge(allele_mismatch[[1]],allele_mismatch[[2]],by="ID",all=T,sort=F)
			if(length(allele_mismatch)>2){
				for(i in 3:length(allele_mismatch)){
					allele_mismatch_snps<-merge(allele_mismatch_snps,allele_mismatch[[i]],by="ID",all=T,sort=F)
				}
			}
			allele_mismatch_snps<-allele_mismatch_snps[,-1]
			for(i in 1:ncol(allele_mismatch_snps)){
				allele_mismatch_snps[,i]<-sort(allele_msimatch_snps[,i],na.last=T)
			}
			allele_mismatch_snps[is.na(allele_mismatch_snps)]<-""
			write_output(allele_mismatch_snps[,-1],out,"not_aligned.txt",output_gz,F)
			rm(allele_mismatch_snps)
		}
		rm(allele_mismatch,N_allele_mismatch,x_common)
	}

	cat("  Aligning rows across input files...")
	SNP_info<-all_aligned[[1]][,1:5]
	for(i in 2:length(all_aligned)){
		dum<-all_aligned[[i]][,1:5]
		dum<-dum[!(dum$SNPID%in%SNP_info$SNPID),]
		if(nrow(dum)!=0){
			SNP_info<-rbind(SNP_info,dum)
		}
	}
	all_aligned<-lapply(all_aligned,function(x,y=data.frame(SNPID=SNP_info$SNPID,stringsAsFactors=F)){
		dum<-merge(x,y,all.y=T,sort=F)
		return(dum[match(y$SNPID,dum$SNPID),])
	})
	cat("done\n\n",format(nrow(SNP_info),big.mark=",")," SNPs were included in the GWAMA\n",sep="")

	if(!is.na(type)){
		if(type=="demontis"){
			cat("\n-----\n\nA dichotomous and a continuous GWA were meta-analyzed using the method described in Demontis et al. (2017)\n  NOTE: if sample sizes and Z-scores have already been adjusted, run the function without type=\"demontis\"\n\n  Performing basic checks")
			if(names(all_aligned)[1]!="dichotomous"|names(all_aligned)[2]!="continuous"){
				stop("the names of x must be \"dichotomous\" and \"continuous\", respectively",call.=F)
			}

			demontis_objects<-list(n_cases,n_controls,pop_prev,rg,h2_dich,h2_cont,lj,M)
			names(demontis_objects)<-c("n_cases","n_controls","pop_prev","rg","h2_dich","h2_cont","lj","M")
			for(i in 1:8){
				if(is.null(demontis_objects[[i]])|is.na(demontis_objects[[i]])){
					stop(paste(names(demontis_objects)[i],"was empty"),call.=F)
				}
				if(!is.numeric(demontis_objects[[i]])){
					stop(paste(names(demontis_objects)[i],"was non-numeric"),call.=F)
				}
				if(length(demontis_objects[[i]])!=1){
					if(names(demontis_objects)[i]=="lj"){
						cat("WARNING: function cannot handle SNP-specific LD scores yet\n  \"lj\" was set to 124.718\n")
					}else{
						stop(paste(names(demontis_objects)[i],"had length >1"),call.=F)
					}
				}
			}
			for(i in 3:6){
				if(demontis_objects[[i]]<0|demontis_objects[[i]]>1){
					stop(paste0(names(demontis_objects)[i]," out of bounds (range: ",min(demontis_objects[[i]])," - ",max(demontis_objects[[i]]),")"),call.=F)
				}
			}
			cat("  Sanity checks completed\n  Adjusting the sample sizes and Z-scores...")
			N<-sum(n_cases,n_controls)
			p_cases<-n_cases/N
			dich<-all_aligned[[1]]
			dich$N<-N*((p_cases*(1-p_cases)*dnorm(qnorm(pop_prev))^2)/((pop_prev*(1-pop_prev))^2))
			all_aligned[[1]]<-dich
			cont<-all_aligned[[2]]
			cont$Z<-cont$Z/sqrt(1+(1-rg^2)*cont$N*h2_cont*lj/M)
			cont$N<-cont$N*((rg^2*h2_cont/h2_dich)/(1+(1-rg^2)*cont$N*h2_cont*lj/M))
			all_aligned[[2]]<-cont
			rm(dich,cont)
			gc()
			cat("done\n")
		}
	}

	cat("\n-----\n\nCalculating GWAMA statistics:\n  Extracting the EAFs from the aligned summary statistics...")
	EAF<-do.call(cbind,lapply(all_aligned,"[[","EAF"))
	cat("done\n  Extracting the Z-scores from the aligned summary statistics...")
	Z<-do.call(cbind,lapply(all_aligned,"[[","Z"))
	cat("done\n")
	if(is.null(cov_Z)){
		cat("  Calculating cov_Z...")
		cov_Z<-cor(Z,use="pairwise")
		rownames(cov_Z)<-colnames(cov_Z)<-names(all_aligned)
		write.table(cov_Z,paste0(out,"corZ.txt"),quote=F,sep="\t",dec=".")
		if(sum(is.na(cov_Z))!=0){
			stop("cov_Z contained missing values. Are all SNPIDs of the same type (e.g. RS-number or CHR:BP)?",call.=F)
		}
		cat("done\n    cov_Z was saved as ",out,"corZ.txt\n",sep="")
	}
	cat("  Extracting the weights from the aligned summary statistics...")
	N<-as.matrix(do.call(cbind,lapply(all_aligned,"[[","N")))
	rm(all_aligned)
	cat("done\n  Multiplying weights with the SNP-heritabilities...")
	W<-t(t(N) * h2)
	cat("done\n  Taking the square root of the weights...")
	sqrt_W<-sqrt(W)
	cat("done\n  Calculating product between Z-scores and sqrt-W...")
	ZxW<-Z*sqrt_W
	cat("done\n  Summing ZxW over the summary statistics...")
	ZxW_total<-matrix(rowSums(ZxW,na.rm=T),ncol=1)
	cat("done\n  Computing multivariate Z-scores...")
	Z_multi<-apply(X=cbind(ZxW_total,sqrt_W),MARGIN=1,FUN=function(x,y=cov_Z){
		idx<-!is.na(x[2:length(x)])
		W<-diag(x[2:length(x)][idx],sum(idx))
		return(x[1]/sqrt(sum(W%*%y[idx,idx]%*%W)))
	})
	cat("done\n  Computing two-sided P-values from multivariate Z-scores...")
	P<-pchisq(q=Z_multi^2,df=1,lower.tail=F)
	cat("done\n    ",format(x=length(x=P[P<=genome_wide_significant & !is.na(P)]),big.mark=",")," SNPs reached genome-wide significance (Minimum P-value: ",sprintf(fmt="%0.4e",min(P,na.rm=T)),")\n  Computing new N...",sep="")
	N_multi<-rowSums(N,na.rm=T)
	cat("done\n  Computing the number of observations per SNP...")
	M<-as.matrix(apply(Z,1,function(x){
		dum<-sign(x)
		dum[dum<0 & !is.na(dum)]<-"-"
		dum[dum>=0 & !is.na(dum)]<-"+"
		dum[is.na(dum)]<-"?"
		dum<-paste0(dum,collapse="")
		return(dum)
	}))
	cat("done\n  Computing new EAFs and MAFs...")
	EAF_multi<-rowSums(EAF*N,na.rm=T)/N_multi
	MAF<-apply(as.matrix(EAF_multi),1,function(x){
		return(min(x,1-x))
	})
	cat("done\n  Computing beta's...")
	B<-Z_multi/sqrt(N_multi*2*EAF_multi*(1-EAF_multi))
	cat("done\n  Computing SE...")
	SE<-(1/sqrt(N_multi))*(1/sqrt(2*EAF_multi*(1-EAF_multi)))

	cat("done\n\n-----\n\nConstructing N-weighted multivariate GWAMA summary statistics...")
	sumstats<-cbind(NA,SNP_info,EAF_multi,MAF,N_multi,M,B,SE,Z_multi,P)
	colnames(sumstats)<-c("CPTID","SNPID","CHR","BP","EA","OA","EAF","MAF","N","M","BETA","SE","Z","PVAL")
	sumstats$CPTID<-paste(sumstats$CHR,sumstats$BP,sep=":")
	cat("done\n")
	order_output<-check_logical(order_output,"order_output",F)
	if(order_output){
		cat("  Sorting the summary statistics on CHR:BP (might take a while)...")
		sumstats<-sumstats[with(sumstats,order(CPTID)),]
		cat("done\n")
	}
	cat("Successfully created summary statistics for",format(x=nrow(sumstats),big.mark=","),"SNPs\n")
	if(output_gz){
		cat("The output was saved as ",out,"results.txt.gz\n",sep="")
	}else{
		cat("The output was saved as ",out,"results.txt\n",sep="")
	}
	write_output(sumstats,out,"results.txt",output_gz,T)

	end.time<-Sys.time()
	analysis_duration<-difftime(time1=end.time,time2=begin.time,units="secs")
	analysis_mins<-floor(x=floor(x=analysis_duration/60))
	analysis_secs<-analysis_duration-(60*analysis_mins)
	cat("\n-----\n\nAnalysis successfully finished at:",format(x=end.time,trim="%y-%b-%d %X",usetz=T),"\nAnalysis took",analysis_mins,"minutes and",analysis_secs,"seconds\n")
	closeAllConnections()
}

get_outDir_name<-function(x){
	y<-unlist(strsplit(x,"/"))
	name<-gsub("\\.$","",tail(y,n=1))
	out_dir<-paste0(c(y[1:length(y)-1],""),collapse="/")
	if(out_dir=="./"){
		out_dir<-paste0(getwd(),"/")
	}
	return(list(out_dir,name))
}

check_logical<-function(x,arg,default){
	if(!is.logical(x)){
		if(x==1){
			return(T)
		}else if(x==0){
			return(F)
		}else{
			cat("WARNING: invalid value in ",arg,"=\"",x,"\"\nReturning to default (",arg,"=,\"",default,"\")\n",sep="")
			return(default)
		}
	}else{
		return(x)
	}
}

write_output<-function(x,out,extension,output_gz,col.names){
	out<-paste0(out,extension)
	if(output_gz){
		gz1<-gzfile(paste0(out,".gz"),"wt")
	}else{
		gz1<-out
	}
	write.table(x,gz1,quote=F,sep="\t",dec=".",row.names=F,col.names=col.names)
	if(output_gz){
		close(gz1)
	}
}
