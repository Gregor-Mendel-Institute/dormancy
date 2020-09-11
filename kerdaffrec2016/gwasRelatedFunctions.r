
######## GWAS related functions #########
######## Set of functions to prepare data for running gwas and visualize results #########

# source("/Users/envel.kerdaffrec/GMI/projects/rscripts/gwasRelatedFunctions.r")
# source("/Users/envel.kerdaffrec/GMI/projects/rscripts/customPlots.r")




### make SNP info ###
# x <- SNP matrix used for GWAS

generateSnpInfo <- function(x) {
	snp_info <- data.frame(colnames(x), matrix(unlist(strsplit(colnames(x), "- ")), ncol=2, byrow=T))
	snp_info[,2:3] <- apply(snp_info[,2:3], 2, function(x) as.numeric(as.character(x)))
	colnames(snp_info) <- c('SNP', 'Chr', 'Pos')
	snp_info	
}



### Compute MAC and MAF ###
# x <- SNP matrix used for GWAS

computeMacMaf <- function(x) {
	mac_maf <- data.frame(SNP=colnames(x), AC_1=apply(x, 2, sum))	
	mac_maf$AC_0 <- nrow(x) - mac_maf$AC_1
	mac_maf$MAC <- apply(mac_maf[,2:3], 1, min)
	mac_maf$MAF <- mac_maf$MAC/nrow(x)
	mac_maf
}

### Compute MAC and MAF ###
# x <- SNP matrix used for GWAS
# mac_maf <- output of computeMacMaf

filterXMacMaf <- function(x, mac_maf, mac.maf=0) {
	cat("Number of SNPs before filtering:", ncol(x), "\n")
	if (mac.maf == 0) {
		X_ok <- x
		cat("No SNPs were filtered out based on MAC / MAF", "\n")
	} else if (mac.maf >= 1) {
		to_prune <- subset(mac_maf, MAC < mac.maf)
		X_ok <- x[,!colnames(x) %in% to_prune$SNP]
		cat("Number of SNPs after MAC filtering:", ncol(X_ok), "\n")
	} else if (mac.maf < 1) {
		to_prune <- subset(mac_maf, MAF < mac.maf)
		X_ok <- x[,!colnames(x) %in% to_prune$SNP]
		cat("Number of SNPs after MAF filtering:", ncol(X_ok), "\n")		
	}
	X_ok
}



### Filter gwas output based on coordinates ###
# x <- gwas output, with colnames "Chr" and "Pos"

filterRegion <- function(x, chr, pos1, pos2) {	
	cat("Number of SNPs before filtering:", nrow(x), "\n")
	output_chr <- subset(x, Chr == chr)
	output_pos <- subset(output_chr, Pos >= pos1 & Pos <= pos2)
	cat("Number of SNPs after filtering based on coordinates:", nrow(output_pos), "\n")
	output_pos
}


### Filter gwas output based on MAC or MAF ###
# x <- gwas output, with colnames "MAC" and / or "MAF"
# mac.maf.cutoff <- 0 = no filtering 
#					greater than 0 to smaller than 1 = maf filtering
#					1 or greater = mac filtering

filterMacMaf <- function(x, mac.maf=0) {
	if (!exists("mac.maf")) { stop }
	if (mac.maf == 0) {
		output_macf <- x
		cat("No SNPs were filtered out based on MAC / MAF", "\n")
	} else if (mac.maf >= 1) {
		output_macf <- subset(x, MAC >= mac.maf)
		cat("Number of SNPs after MAC filtering:", nrow(output_macf), "\n")
	} else if (mac.maf < 1) {
		output_macf <- subset(x, MAF >= mac.maf)
		cat("Number of SNPs after MAF filtering:", nrow(output_macf), "\n")		
	}
	output_macf
}


### Filter gwas output based on P-values ###
# x <- gwas output, with colname "Pval"
# pval.filt <- -log10 cut off P-values

filterPvalues <- function(x, pval.filt=0) {
	
	if (pval.filt == 0) {
		output_filt <- x
		cat("No SNPs were filtered out based on P-values", "\n")
	} else {
		output_filt <- subset(x, -log10(x$Pval) >= pval.filt)
		cat("Number of SNPs after P-values filtering:", nrow(output_filt), "\n")
	}
	output_filt
}



### Compute LD between SNP of interest and the rest of the markers ###
### Assigns colors to SNPs based on LD ###
# x <- gwas output, with colname "SNP"
# ld.snp <- "5- 18570773"
# X <- SNP matrix
# acc.list <- exact list of accessions used for GWAS to filter out X
# plot_col <- plotColLD(output_filt, "5- 18570773", X, rownames(Y))

plotColLD <- function(x, ld.snp, X, acc.list) {
		
	if (missing("ld.snp")) { stop }
	if (missing("X")) { stop }
	if (missing("acc.list")) { stop }

	cat("Computing LD relative to", ld.snp, "\n")
	ld <- X[, colnames(X)%in%x$SNP]
	ld <- ld[rownames(ld)%in%acc.list,]
	snp_mean <- apply(ld, 2, mean)
	snp_sd <- apply(ld, 2, sd)
	mat_mean <- matrix(nrow=nrow(ld), ncol=ncol(ld), data=snp_mean, byrow=T)
	mat_sd <- matrix(nrow=nrow(ld), ncol=ncol(ld), data=snp_sd, byrow=T)
	ld_stand <- (ld-mat_mean)/mat_sd
	r2 <- (crossprod(ld_stand, ld_stand)/nrow(ld))^2
	r2 <- round(r2/r2[1,1], 8)
	ld_snp <- r2[rownames(r2)%in%ld.snp, ]

	# set colors based on LD
		
	#source("/Users/envel.kerdaffrec/GMI/projects/rscripts/customPlots.r")
	plot_col <- customColors(ld_snp, hues=c("black","#0000FF","#4876FF", "#FFFF00", "#FFA500","#FF0000", "#8B0000"), n=20)
	plot_col
}


### Calculate and add Bonferroni threshold to plots ###
# x <- number of snps

addBonf <- function(x) {
	bonf_thresh <- -log10(0.05 / x)
	cat(paste("Bonferroni treshold =", bonf_thresh), "\n")
	abline(h=bonf_thresh, col=1, lty=3)
}	
	

### Plot P-values genome wide ###

manGenWide <- function(x, d, ax, ax.line, ax.ann, ax.ann.line, ax.ann.cex, ax.lab, ax.lab.line, ax.lab.cex, ...) {

	plot(x$xpos, -log10(x$Pval), ...)

	# x axis
	
	if (ax[1] == 1) {	
		axis(1, lwd=1, tcl=-0.2, labels=NA, at=d, line=ax.line[1])
	}	
	if (ax.ann[1] == 1) {
		axis(1, lwd=0, line=ax.ann.line[1], cex.axis=ax.ann.cex[1], at=d, labels=c(1:max(as.character(unique(x$Chr)))))
	}		
	if (ax.lab[1] == 1) {	
		mtext("Chromosome", side=1, line=ax.lab.line[1], cex=ax.lab.cex[1])
	}

	# y axis
			
	if (ax[2] == 1) {	
		axis(2, lwd=1, tcl=-0.2, labels=NA, line=ax.line[2])
	}			
	if (ax.ann[2] == 1) {		
		axis(2, lwd=0, line=ax.ann.line[2], cex.axis=ax.ann.cex[2], las=1)
	}	
	if (ax.lab[2] == 1) {
		
		lab <- expression(-log[10](italic(P)~value))
		mtext(lab, side=2, line=ax.lab.line[2], cex=ax.lab.cex[2], las=3)
		#mtext("-log10(P-value)", side=2, line=ax.lab.line[2], cex=ax.lab.cex[2], las=3)
	}

}


### Plot P-values locally ###

manLocal <- function(x, chr, pos1, pos2, ax, ax.line, ax.ann, ax.ann.line, ax.ann.cex, ax.lab, ax.lab.line, ax.lab.cex, x.offset, xlim, ...) {
	
	plot(x$Pos, -log10(x$Pval), xlim=c(xlim[1]-(xlim[2]-xlim[1])/x.offset, xlim[2]+(xlim[2]-xlim[1])/x.offset), ...)

	# x axis
	
	if (ax[1] == 1) {	
		axis(1, lwd=1, tcl=-0.2, labels=NA, line=ax.line[1])
	}	
	if (ax.ann[1] == 1) {
		axis(1, lwd=0, line=ax.ann.line[1], cex.axis=ax.ann.cex[1])
	}		
	if (ax.lab[1] == 1) {	
		mtext(paste("Position on Chr. ", chr, " (bp)", sep=""), side=1, line=ax.lab.line[1], cex=ax.lab.cex[1])
	}


	# y axis
			
	if (ax[2] == 1) {	
		axis(2, lwd=1, tcl=-0.2, labels=NA, line=ax.line[2])
	}			
	if (ax.ann[2] == 1) {		
		axis(2, lwd=0, line=ax.ann.line[2], cex.axis=ax.ann.cex[2], las=1)
	}	
	if (ax.lab[2] == 1) {
		lab <- expression(-log[10](italic(P)~value))
		mtext(lab, side=2, line=ax.lab.line[2], cex=ax.lab.cex[2], las=3)
		#mtext("-log10(P-value)", side=2, line=ax.lab.line[2], cex=ax.lab.cex[2], las=3)
	}
}	



### Plot RSS ###

plotRSS <- function(x, ax.ann=c(1,1), ax.ann.line=c(-0.4,-0.4), ax.ann.cex=c(1,1), ax.lab=c(1,1), ax.lab.line=c(1.6,1.8), ax.lab.cex=c(0.6,0.6), ...) {
	
	RSS <- x$RSSout[-(nrow(x$RSSout)/2),]
	
	print(RSS)
	
	plot(0, type="n", xlim=c(0, nrow(RSS)-1), ylim=c(0,1), axes=F, ann=F, xaxs="i", yaxs="i", ...)
	
	axis(1, tcl=-0.2, labels=NA, line=0)
	if (ax.ann[1] == 1) {
		axis(1, lwd=0, line=ax.ann.line[1], cex.axis=ax.ann.cex[1])
	}
	if (ax.lab[1] == 1) {
		mtext("Step number", side=1, line=ax.lab.line[1], cex=ax.lab.cex[1])
	}
		
	axis(2, tcl=-0.2, labels=NA, line=0)
	if (ax.ann[2] == 1) {
		axis(2, lwd=0, line=ax.ann.line[2], cex.axis=ax.ann.cex[2], las=1)
	}
	if (ax.lab[2] == 1) {
		mtext("Variance", side=2, line=ax.lab.line[2], cex=ax.lab.cex[2], las=3)
	}
		
	polygon(c(0:(nrow(RSS)-1), nrow(RSS), 0), c(RSS[,3],0,0), col="SteelBlue3", border=0)
	polygon(c(0:(nrow(RSS)-1), nrow(RSS), 0), c(RSS[,2],0,0), col="SteelBlue1", border=0)
	polygon(c(0:(nrow(RSS)-1), nrow(RSS), 0), c(RSS[,1],0,0), col="Orange2", border=0)
	
	box()

	abline(v=(nrow(RSS)-1)/2, lty=2)

}





### QQ plot ###

gwasQQPlot <- function(x, mac.maf=0.05, pval.filt=2, xlab='Expected', ylab='Observed', ...){
	
	x <- filterMacMaf(x, mac.maf)
	
	e <- -log10(ppoints(nrow(x)))
	o <- -log10(sort(x$Pval))
	
	e <- e[o >= pval.filt]
	o <- o[o >= pval.filt]
	
	lims <- c(pval.filt, ceiling(max(c(e,o))))
	
	customPlot(e, o, xlim=lims, ylim=lims, xlab=xlab, ylab=ylab, ...)
	abline(0,1,col="dark grey")
	box()
	
}


### Prep data for MLMM ###

gwasPrepMlmm <- function(Y, pheno.name, X, K, mac.maf=2){
	
	Y <- Y[order(Y$id), ]
	Y <- Y[!is.na(Y[,pheno.name]), ]
	pheno_tmp <- as.numeric(Y[,pheno.name])
	names(pheno_tmp) <- Y$id
		
	# match X and Y
	
	X <- X[rownames(X) %in% names(pheno_tmp), ]
	pheno <- pheno_tmp[names(pheno_tmp) %in% rownames(X)]
	cat(dim(X), "dim X", "\n")
	cat(length(pheno), "length Y", "\n")

	# compute MAF and MAC and do some filtering
	
	snp_info <- generateSnpInfo(X)
	mac_maf <- computeMacMaf(X)
	mac_maf <- merge(snp_info, mac_maf, by='SNP')
	X_ok <- filterXMacMaf(X, mac_maf, mac.maf)	
	
	# subset K
	
	K_ <- K[rownames(K) %in% names(pheno), ]
	K_ok <- K_[, colnames(K_) %in% names(pheno)]

	# out

	out <- list(pheno, X_ok, K_ok, mac_maf[,-c(4,5)])
	names(out) <- c('Y', 'X', 'K', 'SNP_INFO')
	out
}


gwasPrepMlmmLocal <- function(Y, pheno.name, X, K, mac.maf=2, chr, pos1, pos2){
	
	Y <- Y[order(Y$id), ]
	Y <- Y[!is.na(Y[,pheno.name]), ]
	pheno_tmp <- as.numeric(Y[,pheno.name])
	names(pheno_tmp) <- Y$id
		
	# match X and Y
	
	X <- X[rownames(X) %in% names(pheno_tmp), ]
	pheno <- pheno_tmp[names(pheno_tmp) %in% rownames(X)]
	cat(dim(X), "dim X", "\n")
	cat(length(pheno), "length Y", "\n")

	# compute MAF and MAC and do some filtering
	
	snp_info <- generateSnpInfo(X)	
	snp_info <- subset(snp_info, Chr==5 & Pos >= pos1 & Pos <= pos2)
	X <- X[,colnames(X) %in% snp_info$SNP]
	mac_maf <- computeMacMaf(X)
	mac_maf <- merge(snp_info, mac_maf, by='SNP')
	X_ok <- filterXMacMaf(X, mac_maf, mac.maf)	
	
	# subset K
	
	K_ <- K[rownames(K) %in% names(pheno), ]
	K_ok <- K_[, colnames(K_) %in% names(pheno)]

	# out

	out <- list(pheno, X_ok, K_ok, mac_maf[,-c(4,5)])
	names(out) <- c('Y', 'X', 'K', 'SNP_INFO')
	out
}




### Prep MLMM data for plotting ###


gwasPrepPlotMlmm <- function(out, snp.info, main.col, cof.col) {
	
	cofs <- lapply(out[[2]], function(x) x[[2]])

	pvals <- lapply(out[[2]], function(x) { 
		x <- x[[1]]
		rownames(x) <- NULL
		x <- merge(snp.info, x, by="SNP")
		x <- data.frame(x, Score=-log10(x$pval))
		x <- data.frame(x, Snp_col=main.col, Plot_order=0)
		x$Snp_col <- as.character(x$Snp_col)
		colnames(x)[6] <- 'Pval'
		return(x)
	})

	# assign different color for cofactors and make sure they are plotted last
	
	for (i in (1:length(pvals))) { 
		pvals[[i]]$Snp_col[pvals[[i]]$SNP %in% cofs[[i]]] <- cof.col
		pvals[[i]]$Plot_order[pvals[[i]]$SNP %in% cofs[[i]]] <- 1
		pvals[[i]] <- pvals[[i]][order(pvals[[i]]$Plot_order),]

	}
	
	pvals
	
}


# convert pvals csv files from GWAS browser to .rda file

gwasBrowserCsvToRda <- function(inpdir, inpfile){
	
	input <- read.csv(paste(inpdir, inpfile, sep="/"))
	
	colnames(input) <- sub('chr', 'Chr', colnames(input))
	colnames(input) <- sub('pos', 'Pos', colnames(input))
	colnames(input) <- sub('maf', 'MAF', colnames(input))
	colnames(input) <- sub('mac', 'MAC', colnames(input))
	colnames(input) <- sub('score', 'Score', colnames(input))
	
	input$Pval <- 1/10^input$Score
	input$SNP <- paste(input$Chr, '- ', input$Pos, sep='')
	
	output <- input[,c('SNP', 'Chr', 'Pos', 'MAF', 'MAC', 'Pval', 'Score')]	
	
	output <- output[order(output$Chr, output$Pos),]
	
	obj <- sub(".csv", "", inpfile)
	assign(obj[1], output)
	
	save(list=obj[1], file=paste(inpdir, sub(".csv", ".rda", inpfile), sep="/"))	
}


# convert R GWAS output to .csv file readable by the GWAS browser

gwasRToBrowser <- function(output, fileout){
	
	colnames(output) <- sub('Chr', 'chr', colnames(output))
	colnames(output) <- sub('Pos', 'pos', colnames(output))
	colnames(output) <- sub('MAF', 'maf', colnames(output))
	colnames(output) <- sub('MAC', 'mac', colnames(output))
	colnames(output) <- sub('Pval', 'pvalue', colnames(output))

	output$GVE <- NA
	
	output <- output[,c('chr','pos','pvalue','maf','mac','GVE')]
		
	write.csv(output, fileout, row.names=F)
	
}

# convert limix .csv output to .csv file readable by the GWAS browser

gwasLimixToBrowser <- function(filein, fileout){
	
	output <- read.csv(filein)
	
	colnames(output) <- sub('Chr', 'chr', colnames(output))
	colnames(output) <- sub('Pos', 'pos', colnames(output))
	colnames(output) <- sub('MAF', 'maf', colnames(output))
	colnames(output) <- sub('MAC', 'mac', colnames(output))
	colnames(output) <- sub('Pval', 'pvalue', colnames(output))

	output$GVE <- NA
	
	output <- output[,c('chr','pos','pvalue','maf','mac','GVE')]
		
	write.csv(output, fileout, row.names=F)
	
}




# convert R GWAS MLMM output to .csv file readable by the GWAS browser

gwasMlmmToBrowser <- function(mlmm_pvals, pref, nsteps=5){
	
	x <- mlmm_pvals	
	
	for (step in (1:nsteps)){
		
		xx <- x[[step]]
		
		out_name_browser <- paste(pref, '_step_', step-1, '_browser.csv', sep='')
			
		print(out_name_browser)
				
		xx <- xx[order(xx$Chr, xx$Pos),]
		colnames(xx) <- c('SNP', 'chr', 'pos', 'mac', 'maf', 'pvalue', 'score', 'snp.col', 'cofactor')
		
		xx <- xx[,c(2,3,6,5,4)]
		xx <- cbind(xx, GVE=NA)
		write.csv(xx, file=out_name_browser, row.names=F)
		
	}
}


# filter out gwas output

gwasFilterOutput <- function(output, mac.maf=0.05, chr, pos1, pos2, ord.pval=T){
	
	if (!missing(chr)){
		output <- subset(output, Chr==chr)
		if (!missing(pos1) & !missing(pos1)){
			output <- subset(output, Pos >= pos1 & Pos <= pos2)
		}
	} 
	
	output <- filterMacMaf(output, mac.maf)
	
	if (ord.pval){
		output <- output[order(output$Pval),]
	}
	
	output 
}


# get pvalues of genes from GWAS results


getGenesPvalues <- function(genes, gwas.out, win=15000){
		
	genes_chr <- as.integer(unique(genes$Chr))
	gwas_chr <- as.integer(unique(gwas.out$Chr))
	
	if (identical(genes_chr, gwas_chr)){
		
		out <- list()
		
		for (chr in genes_chr){
			
			cat(chr, '\n')
			
			genes_sub <- subset(genes, Chr==chr)
			gwas_sub <- subset(gwas.out, Chr==chr)
			gwas_sub <- gwas_sub[order(gwas_sub$Pos),]
			
			intervals <- list()
						
			for (g in (1:nrow(genes_sub))){
				
				#print(g)
				intervals[[g]] <- gwas_sub[gwas_sub$Pos >= genes_sub[g,4] - win  & gwas_sub$Pos <= genes_sub[g,5] + win,]
				
			}
			
			snps <- lapply(intervals, function(x) x[which(x$Pval==min(x$Pval)),][1,])
			
			out[[chr]] <- do.call("rbind", snps)
		
		}

		out <- do.call("rbind", out)	
		out <- cbind(genes, out[,c(2:ncol(out))])
		out
	}
}




# filter gwas results for horizontal plotting

prepGWASPlot <- function(gwas.out, mac.maf=10, pval.filt=5, plot.col=c('#63B8FF', '#36648B'), spacer=10^6){
	
	x <- filterMacMaf(gwas.out, mac.maf)
	x <- filterPvalues(x, pval.filt)
	
	if (nrow(x)!=0){
		
		x <- x[order(x$Chr, x$Pos),]
		x$Score <- -log10(x$Pval)

		x$Col <- plot.col[1]
		#x$Col[x$Score >= 5.5] <- plot.col[1]
		x$Col[x$Score >= 6.0] <- plot.col[1]
		x$Col[x$Score >= 7.0] <- plot.col[2]
		x$Col[x$Score >= 7.5] <- plot.col[3]
		
		size <- aggregate(x$Pos, list(x$Chr), length)$x
	
		chrs <- c(30427700, 19698300, 23460000, 18585000, 26975600)
		cumpos <- cumsum(chrs + spacer)
		cumpos[length(chrs)] <- cumpos[length(chrs)] - spacer
		cumpos <- c(0, cumpos)
	
		a <- 0
		ch <- unique(x$Chr)
		for (i in 1:length(ch)){
 			a <- c(a, rep(cumpos[ch[i]], size[i]))
 		}
 		a <- a[-1]
	
		x$xPos <- x$Pos + a
		x <- x[order(x$Score),]	
		x
	
	} else {
		
		x <- matrix(NA, nrow=1, ncol=8, byrow=T)
		colnames (x) <- c('Chr', 'Pos', 'MAC', 'MAF', 'Pval', 'Score', 'Col', 'xPos')
		x <- as.data.frame(x)
		
	}	
	x		
}	


# get the overlap between the GWAS results

getGWASOverlap <- function(filt, cond.names, chr, pos){
	b <- lapply(filt, function(x) subset(x, Chr==chr & Pos==pos))
	for (i in (1:length(b))){
		if (nrow(b[[i]])==0){
			b[[i]][1,] <- NA
		}
	}
	b <- do.call('rbind', b)
	rownames(b) <- cond.names
	b
}

# add adjusted x pos to GWAS output, or whatever output as long as 'Chr' and 'Pos' are provided

gwasAddXPos <- function(x, chrs.length=c(30427700, 19698300, 23460000, 18585000, 26975600), spacer=10^6){

		x$Chr <- as.numeric(x$Chr)		
		x$Pos <- as.numeric(x$Pos)
		
		size <- aggregate(x$Pos, list(x$Chr), length)$x
		cumpos <- cumsum(chrs.length + spacer)
		cumpos[length(chrs.length)] <- cumpos[length(chrs.length)] - spacer
		cumpos <- c(0, cumpos)
	
		a <- 0
		ch <- unique(x$Chr)
		for (i in 1:length(ch)){
 			a <- c(a, rep(cumpos[ch[i]], size[i]))
 		}
 		a <- a[-1]
	
		x$xPos <- x$Pos + a
		x
}	



# Fonction to parse the HDF5 data used by PyGWAS and teh GWAS browser. Creates entire or chunked SNP matrices that are accepted by my R GWAS scripts. Saves them as .rda files. To reduce the weight of the matrices MAC / MAF and accession filtering can be applied.

# hd5f.file = path and name of input file.
# out.file = path and name of output files.
# acc.to.keep = vector containing accessions to keep. Default is none, all accessions are kept.
# macmaf = threshold for filtering. MAF, value < 1. MAC, value >= 1.
# chunk.size = number of snps in each chunk.

##ACTHUNG: for some reason h5read crashes when loading the snp matrix, both locally and on the cluster... so useless function until the problem is fixed.

# hd5f.file="/Users/envel.kerdaffrec/GMI/projects/data/genotypes/1001_all_chromosomes_binary_92Swedes.hdf5"; out.file="/Users/envel.kerdaffrec/GMI/projects/data/genotypes/1001_all_chromosomes_binary_92Swedes"; acc.to.keep=none; macmaf=0; chunk.size=500000

gwasHDF5toRda <- function(hd5f.file, out.file, acc.to.keep, macmaf=0, chunk.size=500000){
    
    require(rhdf5)
    
    acc <- h5read(hd5f.file, "/accessions")
    pos <- h5read(hd5f.file, "/positions")
    snp <- h5read(hd5f.file, "/snps")
    
    chrs <- unlist(h5readAttributes(hd5f.file, "/positions")[2])
    chr_regions <- matrix(unlist(h5readAttributes(hd5f.file, "/positions")[1]), ncol=2, byrow=T)
    chr <- rep(chrs, chr_regions[,2]-chr_regions[,1])
    snp_coords <- paste(chr, "- ", pos)    
    colnames(snp) <- snp_coords
    rownames(snp) <- acc
    print(snp[1:5,1:5])
    
    if (ncol(snp) > chunk.size){
        chunk_start <- seq(1,ncol(snp), chunk.size)
        chunk_end <- c(chunk_start[-1]-1, ncol(snp))
        chunks <- cbind(chunk_start, chunk_end)
        nrow(chunks)
    
        chunk_list <- list()
    
        for (i in (1:nrow(chunks))){
            chunk <- snp[,chunks[i,1]:chunks[i,2]]
            
            if (acc.to.keep[1]!="none"){
            	chunk_acc_filt <- chunk[rownames(chunk)%in%acc.to.keep,]
            	cat("Number of accessions after filtering: ", nrow(chunk_acc_filt), "\n")
            
            } else {
            	chunk_acc_filt <- chunk
            }
            	
            macmaf <- computeMacMaf(chunk_acc_filt)
            X <- filterXMacMaf(chunk, macmaf, 2) 
            out_chunk <- paste(out.file, "_chunk", i, ".rda", sep='')
            save(X, file=out_chunk)
            chunk_list[[i]] <- X
            rm(X)     
        }   
   
    } else {
        snp_acc_filt <- snp[rownames(snp)%in%acc_to_keep,]
        X <- snp
        out <- paste(out.file, ".rda", sep='')
        save(X, file=out_chunk)
    }
}




# Fonction to convert snp csv file (tipically the output from slice_hdf5_to_csv.py) to rda snp file that suits my GWAS R functions.

# snpsCsvToRda(file.dir = "/Users/envel.kerdaffrec/GMI/projects/data/genotypes/CM78", file.pattern ="CM78_binary_mac0_180acc_100k_chunk")


snpsCsvToRda <- function(file.dir, file.pattern){
   
    file_path <- list.files(file.dir, pattern=file.pattern, full.names=T)
    file_path_csv <- file_path[grep("csv", file_path)]
    file_path_out <- sub("csv", "rda", file_path_csv)
    
    for (f in (1:length(file_path_csv))){
        chunk <- read.table(file_path_csv[f], sep=",", header=T, check.names=F)
        chunk[1:5,1:5]
        rownames(chunk) <- chunk$id
        chunk <- chunk[,-1]
        save(chunk, file=file_path_out[f])  
    }
}





# Fonction to convert snp INFO csv file (tipically the output from slice_hdf5_to_csv.py) to rda snp INFO file that suits my GWAS R functions.

# snpInfoCsvToRda(file.dir = "/Users/envel.kerdaffrec/GMI/projects/data/genotypes/CM78", file.pattern ="info", combine = TRUE, file.out = "CM78_binary_mac0_180acc_SNP_INFO.rda")



snpInfoCsvToRda <- function(file.dir, file.pattern, combine = TRUE, file.out){
    
    file_path <- list.files(file.dir, pattern=file.pattern, full.names=T)
    
    snp_info <- list()    
    for (f in (1:length(file_path))){
        snp_info[[f]] <- read.table(file_path[f], sep=",", header=T, check.names=F)
    }
    
    if(combine[1]==TRUE){
        snp_info <- do.call("rbind", snp_info)
        snp_info <- snp_info[order(snp_info$Chr & snp_info$Pos),]
        file_path_out <- paste(file.dir, file.out, sep="/")
        save(snp_info, file=file_path_out)
    } else {
        file_path_out <- sub("csv", "rda", file_path)
        for (i in (1:length(snp_info))){
            save(snp_info[[i]], file=file_path_out[[i]])
        }
    }
}




