################################################################################
################################ AMM GWAS ######################################
################################################################################

# Last modified: 13 Sept 2018.
# Author: Envel Kerdaffrec.
# Function to perform GWAS scans genome wide or locally, depending on the provided genotypic data.
# Allows the inclusion of SNPs as cofactors.

# source("/Users/envel.kerdaffrec/GMI/projects/scripts/gwasAmmCof.r")


# Y: a data.frame with at least two columns, Y$id and Y$pheno.name.
# pheno.name: the name of the column containing the phenotype data
# X: the snp matrix to use. Rows are accessions, Cols are SNPs coded with 0 and 1.
# K: kinship matrix.
# mac.maf: MAC (if >= 1) or MAF (if > 0 & < 1) threshold.
# cof: SNPs to be included as cofactors. Should be a vector c('5- 18590327', '5- 18590501') or a data frame with an id column and colums corresponding to SNPs coded with 0 and 1.


# gwasAmmCof(Y=GR21, pheno.name='gr21', X=DOG1_SNPS, K=K, macf=2)


gwasAmmCof <- function(Y, pheno.name, X, K, mac.maf=2, cof){

	#source("/Users/envel.kerdaffrec/GMI/projects/rscripts/gwasRelatedFunctions.r")
	#source("/Users/envel.kerdaffrec/GMI/projects/rscripts/gwasEmma.r")
	

# prepare Y

	Y <- Y[order(Y$id), ]
	pheno_tmp <- matrix(Y[,pheno.name])
	pheno_tmp <- matrix(pheno_tmp[!is.na(pheno_tmp), ])
	rownames(pheno_tmp) <- Y$id
	colnames(pheno_tmp) <- pheno.name
	
	
	
# match X and Y
	
	X <- X[rownames(X) %in% rownames(pheno_tmp), ]
	X <- X[order(as.numeric(rownames(X))),]
	pheno <- as.matrix(pheno_tmp[rownames(pheno_tmp) %in% rownames(X), ])
	cat(dim(X), "dim X", "\n")
	cat(dim(pheno), "dim Y", "\n")
	#print(head(pheno))
	#print(X[1:10,1:10])
	
	

# compute MAF and MAC and do some filtering
	
	snp_info <- generateSnpInfo(X)
	mac_maf <- computeMacMaf(X)
	mac_maf <- merge(snp_info, mac_maf, by='SNP')
	X_ok <- filterXMacMaf(X, mac_maf, mac.maf)	
	
	

# normalisation of K
	
	K_ <- K[rownames(K) %in% rownames(pheno), ]
	K_ok <- K_[, colnames(K_) %in% rownames(pheno)]
	n <- nrow(pheno)
	K_norm <- (n-1)/sum((diag(n)-matrix(1,n,n)/n)*K_ok)*K_ok
	cat("dim K", dim(K_norm), "\n")



# get cofactors

	if (!missing(cof)) {
		if (class(cof)=='character') { 
			cof_ok <- as.matrix(X_ok[,colnames(X_ok) %in% cof])
			#print(dim(cof_ok))
			colnames(cof_ok) <- cof
			rownames(cof_ok) <- rownames(X_ok)
			X_ok <- X_ok[,!colnames(X_ok) %in% cof]
		}
		if (class(cof)=='matrix' | class(cof)=='data.frame') {
			cof_ok <- as.matrix(cof[cof$id %in% rownames(X_ok),-1])
			colnames(cof_ok) <- colnames(cof)[-1] 
			rownames(cof_ok) <- rownames(X_ok)
			class(cof_ok) <- 'numeric'
			X_ok <- X_ok[,!colnames(X_ok) %in% colnames(cof_ok)]
		}
	}		



# set the intercept

	if (!missing(cof)) {
		Xo <- rep(1,nrow(X_ok))
		ex <- cbind(Xo, cof_ok)
	} else {
		ex <- matrix(rep(1, nrow(X_ok)))
	}
		


# REML
	
	null <- emma.REMLE(pheno, ex, K_norm) #run null model, without SNPs
	herit <- null$vg / (null$vg + null$ve)
	cat("pseudo-heritability =", herit, "\n")
	
	M <- solve(chol(null$vg * K_norm + null$ve * diag(dim(K_norm)[1])))
	Y_t <- crossprod(M, pheno)
	ex_t <- crossprod(M, ex)



# EMMAX scan

	cat("EMMAX is running", "\n")
	RSS_H0 <- rep(sum(lsfit(ex_t, Y_t, intercept=FALSE)$residuals^2), ncol(X_ok))
	RSS_H1 <- apply(X_ok, 2, function(x){sum(lsfit(cbind(ex_t, crossprod(M,x)), Y_t, intercept=FALSE)$residuals^2)})

	#models1 <- apply(X_ok,2,function(x){summary(lm(Y_t~0+ex_t+crossprod(M,x)))$coeff[2,]})
	#models2 <- summary(lm(Y_t~0+ex_t+crossprod(M,snp1)+crossprod(M,snp2)))$coeff[2:3,] 

	m <- nrow(Y)
	p_H0 <- ncol(ex_t)
	p_H1 <- ncol(ex_t) + 1
	df1 <- p_H1 - p_H0
	df2 <- m - p_H1

	F <- ((RSS_H0 - RSS_H1) / df1) / (RSS_H1 / df2)
	pval <- pf(F, df1, df2, lower.tail=FALSE)



# generate output

	cat("Creating output", "\n")
	pval_ok <- data.frame(SNP=as.character(names(pval)), Pval=pval)

	output <- merge(mac_maf, pval_ok, by='SNP')
	output <- output[,-c(4,5)]
	output$Score <- -log10(output$Pval)
	output <- output[order(output$Chr, output$Pos),]

	#rm(list=ls(pattern="emma"))

	cat("Done", "\n")

	output

}




############## LOG ################


# 13 Sept 2018: get cofactors code, when cofactors are provided as matrix or data frame. Added "class(cof_ok) <- 'numeric'". 









