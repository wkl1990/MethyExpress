#' Correlation analysis of methylation and expression
#'
#' This function perform correlation analysis of DNA methylation and gene expression between gene and cpg pairs (GCPs).
#'
#' @param methydata Mehtylation data matrix
#' @param expressdata Expression data matrix

#' @param pairsdata The gene and cpg pairs for correlation analysis, column V1 is the cpg probe and column V2 is the gene name
#' @param method The method used to calculate correlation Defaults to \pkg{spearman}
#' @author WKL
#' @keywords GCPs, correlation
#' @return GCPs correlation statistcis
#' @examples
#' MethyExpressCor(methydata, expressdata, pairsdata)
#' @export

MethyExpressCor <- function(methydata, expressdata, pairsdata, method="spearman"){
	methy_overlap_express_sample <- intersect(colnames(methydata),colnames(expressdata))
	methydata_overlap <- methydata[,match(methy_overlap_express_sample,colnames(methydata))]
	expressdata_overlap <- expressdata[,match(methy_overlap_express_sample,colnames(expressdata))]
	#identical(colnames(methydata_overlap),methy_overlap_express_sample)
	#identical(colnames(expressdata_overlap),methy_overlap_express_sample)
	print(paste("Whether methydata and expressdata samples are identical: ", identical(colnames(expressdata_overlap),colnames(methydata_overlap))))

	methyprobe <- as.data.frame(pairsdata[,1])
	methyprobe$id <- 1:nrow(methyprobe)
	Phmethy <- merge(methyprobe,methydata_overlap,by.x="pairsdata[, 1]",by.y="row.names")
	Phmethy <- Phmethy[order(Phmethy$id), ]

	expressprobe <- as.data.frame(pairsdata[,2])
	expressprobe$id <- 1:nrow(expressprobe)
	Pexpress <- merge(expressprobe,expressdata_overlap,by.x="pairsdata[, 2]",by.y="row.names")
	Pexpress <- Pexpress[order(Pexpress$id), ]

	pair_overlap <- as.data.frame(intersect(Phmethy$id, Pexpress$id))
	Phmethy_overlap <- merge(pair_overlap,Phmethy,by.x="intersect(Phmethy$id, Pexpress$id)",by.y="id")
	Pexpress_overlap <- merge(pair_overlap,Pexpress,by.x="intersect(Phmethy$id, Pexpress$id)",by.y="id")

	#calculate the cor
	expdata<-Pexpress_overlap[,-c(1,2)]
	methdata<-Phmethy_overlap[,-c(1,2)]
	expdata<-as.matrix(expdata)
	methdata<-as.matrix(methdata)
	corr <- rep(0,nrow(methdata))
	corrpval <- rep(0,nrow(methdata))
	for(i in 1:nrow(methdata)) {
	corr[i] <- cor.test(expdata[i,],methdata[i,],method=method)$estimate
	corrpval[i] <- cor.test(expdata[i,],methdata[i,],method=method)$p.value
	}
	Ename <- as.vector(Pexpress_overlap[,2])
	Mname <- as.vector(Phmethy_overlap[,2])
	qvalue <- p.adjust(corrpval, method = "fdr")
	correl <- cbind(Ename,Mname,corr,corrpval,qvalue)
	methy_correl <- as.data.frame(correl)
	methy_correl$corr <- as.numeric(as.character(methy_correl$corr))
	methy_correl$corrpval <- as.numeric(as.character(methy_correl$corrpval))
	methy_correl$qvalue <- as.numeric(as.character(methy_correl$qvalue))
	# write.csv(methy_correl,file="methy.corr.spearman.csv")
	methy_correl
}

