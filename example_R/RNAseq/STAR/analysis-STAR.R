### Script for anlaysis of STAR mapped reads ##########

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("edgeR")

require("edgeR")  # for kegga
require("DESeq2"); packageVersion("DESeq2")
require("tximport")
require("readr")
require("ReportingTools")
require("AnnotationDbi")
require("ensembldb")
require("AnnotationHub")
require("RColorBrewer")
require("Biostrings")
require("pvclust")
require("dplyr")
require("ggplot2")
require("pheatmap")
require("PoiClaClu")
require("GOstats")
require("genefilter")
require("Category")
require("ashr")
require("BiocParallel")
require("cowplot")
require("Rsamtools")
require("GenomicFeatures")
require("GenomicAlignments")
require("topGO")

register(MulticoreParam(50))


#library("hexbin")
#BiocManager::install("vidger")
#library("vidger") # Visualization of Differential Gene Expression Results using R

#library("hexbin")
#BiocManager::install("vidger")
#library("vidger") # Visualization of Differential Gene Expression Results using R

## R code for processing BAMS and turning them into a count matrix for deseq2:
## From 
## http://bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#recommended-transcript-abundances-and-the-tximport-pipeline

# R
files = dir(path='../02-STAR/', pattern="*.bam", full.names=T)
bamfiles = BamFileList(files, yieldSize=4e06)

# Defining gene models:
# library("GenomicFeatures")
gfffile = "../ref/ITAG4.0_gene_models.gff"
txdb <- makeTxDbFromGFF(gfffile, format = "gff", circ_seqs = character(), dbxrefTag="gene")
txdb

## According to https://solgenomics.net/organism/Solanum_lycopersicum/genome there should be 34075 protein coding genes in ITAG4.0
ebg <- exonsBy(txdb, by="gene")
ebg

## Read in annotation table, note that ebg uses the "gene" ID without the .1 version, so we trim these in the row names:
txann = read.table("../ref/gene_id-lookup.tsv", header=T, as.is=T, sep='\t', quote='')
rownames(txann) = gsub('..$', '', txann$Name)

## Load GO information:
geneID2GO = readMappings(file="../ref/geneid2go.map")
names(geneID2GO) = gsub('..$', '', names(geneID2GO))

###### Read counting step  ##########
# library("GenomicAlignments")
# library("BiocParallel")

## Count how many reads overlaps each exon in each BAM file, this gets run on ~50 cores, assuming the I/O can keep up.
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=TRUE,
                        ignore.strand=FALSE)

dim(se)

head(assay(se), 3)

# make a little barplot showing how many reads were mapped per sample:
colSums(assay(se))
tbl = colSums(assay(se))
names(tbl) = gsub("Aligned.sortedByCoord.out.bam", '', names(tbl))
png(w=1000, h=800, file="barplot_mapped_reads.png")
par(las=2)
    barplot(tbl, main="Num mapped reads")
    abline(h=2e06, col="red")
dev.off()

## Read metadata and add it to the RangedSummarizedExperiment object
samples = read.table("../SampleMetaData.csv", header=T, as.is=T, sep='\t')
#samples$Treatment = paste(samples$Lipid, samples$Time, samples$Bacteria, sep='.')
#tmp1 = gsub("Aligned.sortedByCoord.out.bam", '', rownames(colData(se)))  

mdata = na.omit(samples[match(rownames(colData(se)), samples$BAMFile), ])
#mdata = rbind(mdata, c("Undetermined", "Mix", "All", "Ap","Mix.All.Ap"))
rownames(mdata) = mdata$SampleID  # Note, that if the rownames aren't set, the assay(se) colnames will get the original names

colData(se) = DataFrame(mdata)

### Ribosomal RNA reads, output counts and samples: ###
se.ann = txann[match(rownames(se), rownames(txann)), ]
rRNA.idx = 1:nrow(se) %in% intersect(grep("rRNA", se.ann$Note), grep("protein", se.ann$Note, invert=T))
rRNA.se = se[rRNA.idx, ]

pdf(file="rRNA_reads_per_sample.pdf", width=10, height=6)
	par(las=2)
	barplot(assay(rRNA.se), main="Total rRNA reads per sample")
	par(mar=c(10,5,5,5))
	boxplot(colSums(assay(rRNA.se)) ~ rRNA.se$Treatment, main="Total rRNA reads by treatment type", xlab='', ylab='reads')
	boxplot(colSums(100*assay(se)[rRNA.idx,])/colSums(assay(se)) ~ se$Treatment, main="Percent rRNA reads by treatment", ylab='% Reads', xlab='')
dev.off()

### Remove rRNA reads before normalization ###
se2 = se[!rRNA.idx, ]

# write.table(file='rRNA-filtered-counts.tsv', row.names=F, 
# 											 data.frame(Sample=colnames(se),
# 											 mapped_to_annotated=colSums(assay(se)), 
# 											 filtered_for_rRNA=colSums(assay(se2))), sep='\t', quote=F)



## Convert to DESeqDataSet for Deseq2 analysis:
#library("DESeq2")
dds = DESeqDataSet(se2, design= ~ Treatment)
#colnames(counts(dds)) = colData(dds)$SampleID


##### log-based transformations for plotting, clustering, etc #####
#  regularized logarithm or rlog, which incorporates a prior on the sample differences (Love, Huber, and Anders 2014)
rld = rlog(dds, blind=T)  # blind=T for QA
head(assay(rld), 3)

# variance stabilizing transformations (VST) (Tibshirani 1988; Huber et al. 2003; Anders and Huber 2010)
vsd = vst(dds, blind=T)
head(assay(vsd), 3)


############# Make some QA/QC plots ###################
pdf(file="QA-QC-plots.pdf", width=16, height=8)

	#Barplots and Boxplots:
	mycolors= colorRampPalette(brewer.pal(11, "Spectral"))(length(levels(dds$Treatment)))
	names(mycolors) = levels(dds$Treatment)

	barplot1 = ggplot2::ggplot(data.frame(sample=colnames(dds), reads=colSums(counts(dds)), Treatment=dds$Treatment), 
					aes(x=sample, y=reads, fill=Treatment)) 
	barplot1 = barplot1	+ geom_col()
	barplot1 = barplot1 + theme(axis.text.x=element_text(angle=50, size=10, vjust=0.5))
	barplot1 = barplot1 + ggtitle("Reads mapped per sample")
	barplot1 = barplot1 + theme(plot.title = element_text(size=20, face="bold", margin=margin(10,0,10,0)))

	df.m = reshape2::melt(log2(counts(dds)+1), id.vars=NULL)
	vplot1 = ggplot2::ggplot(df.m, aes(x=Var2, y=value)) + geom_violin() 
	vplot1 = vplot1 + ggtitle("Log2 mapped read count per feature")
	vplot1 = vplot1 + theme(plot.title = element_text(size=20, face="bold", margin=margin(10,0,10,0)))
	vplot1 = vplot1 + labs(x='', y="log2(read count)")
	vplot1 = vplot1 + theme(axis.text.x=element_text(angle=50, size=10, vjust=0.5))

	plot_grid(barplot1, vplot1, labels="AUTO")


	# par(las=2)
	# par(mfrow=c(1,2))
	# par(mar=c(7,4,4,4))
	# barplot(colSums(counts(dds)), main="Mapped reads", col=mycolors[dds$Treatment])
	# legend('topright', col=mycolors, legend=names(mycolors), pch=20, pt.cex=2)
	# boxplot(log2(counts(dds) + 1 ), col=mycolors[dds$Treatment],
	# 	main="Log2 estimated read count distributions")

	#MDS plots:
	d = dist(t(counts(dds)))
	fit = cmdscale(d, eig=T, k=3)
	par(mfrow=c(1,2))
	plot(x=fit$points[,1], y=fit$points[,2], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type='n',
			xlim=range(fit$points[,1]) * 1.1, ylim=range(fit$points[,2]) * 1.1)
	text(x=fit$points[,1], y=fit$points[,2], labels=row.names(fit$points), col = mycolors[dds$Treatment], font=2)
	legend('bottomright', col=mycolors, legend=names(mycolors), pch=20, pt.cex=2)
	plot(x=fit$points[,1], y=fit$points[,3], xlab="Coordinate 1", ylab="Coordinate 3", main="Metric MDS", type='n',
			xlim=range(fit$points[,1]) * 1.1, ylim=range(fit$points[,3]) * 1.1)
	text(fit$points[,1], fit$points[,3], labels=row.names(fit$points), col = mycolors[dds$Treatment], font=2)

	#Hclust:
	par(mfrow=c(1,2))
	hc = hclust(dist(t(counts(dds))), method='ward.D2')
	plot(hc, main="Hclust with Euclidean distance and Ward.D2 clustering")
	hc = hclust(as.dist(1-abs(cor(counts(dds)))), method='ward.D2')
	plot(hc, main="Hclust with 1-abs(correlation) distance and Ward.D2 clustering")


	# These take forever...
	## See http://research.stowers.org/mcm/efg/R/Visualization/cor-cluster/index.htm for correctional clustering:
	# library(pvclust)
	# cluster.bootstrap = pvclust(counts(dds), nboot=100, method.dist="abscor")
	# plot(cluster.bootstrap)
	# pvrect(cluster.bootstrap)

	#pvclust provides two types of p-values: AU (Approximately Unbiased) p-value and BP 
	# (Bootstrap Probability) value. AU p-value, which is computed by multiscale bootstrap 
	# resampling, is a better approximation to unbiased p-value than BP value computed by normal bootstrap resampling.

	## Plots from THE tutorial 
	# https://www.bioconductor.org/help/workflows/rnaseqGene/#exploratory-analysis-and-visualization

	#Boxplots:
	# par(las=2)
	# boxplot(assay(rld), col=mycolors[as.character(rld$Treatment)],
	# 	main="RLog transformed read count distributions")

	#MDS plots:
	par(mfrow=c(2,1))
	d = dist(t(assay(rld)))
	fit = cmdscale(d, eig=T, k=3)
	par(mfrow=c(1,2))
	plot(x=fit$points[,1], y=fit$points[,2], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS, RLog transformed counts", type='n',
			xlim=range(fit$points[,1]) * 1.1, ylim=range(fit$points[,2]) * 1.1)
	text(x=fit$points[,1], y=fit$points[,2], labels=row.names(fit$points), col = mycolors[rld$Treatment], font=2)
	legend('bottomright', col=mycolors, legend=names(mycolors), pch=20, pt.cex=2)
	plot(x=fit$points[,1], y=fit$points[,3], xlab="Coordinate 1", ylab="Coordinate 3", main="Metric MDS, RLog transformed counts", type='n',
			xlim=range(fit$points[,1]) * 1.1, ylim=range(fit$points[,3]) * 1.1)
	text(fit$points[,1], fit$points[,3], labels=row.names(fit$points), col = mycolors[rld$Treatment], font=2)

	library("dplyr")
	library("ggplot2")

	dds <- estimateSizeFactors(dds)

	df <- bind_rows(
	as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
			mutate(transformation = "log2(x + 1)"),
	as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
	as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))
	
	colnames(df)[1:2] <- c("x", "y")  

	# These plots are supposed to show that the log2(x+1) transformation amplifies
	# differences when the values are close to 0.
	ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
		ggtitle("Low expression genes are more variable in log space, RLD and VST are supposed to fix this") + 
	coord_fixed() + facet_grid( . ~ transformation)  


	# Euclidean distances between samples and pretty heatmaps:
	sampleDists <- dist(t(assay(rld)))
	library("pheatmap")
	library("RColorBrewer")

	sampleDistMatrix <- as.matrix( sampleDists )
	rownames(sampleDistMatrix) <- rld$Treatment
	colnames(sampleDistMatrix) <- NULL
	colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
	pheatmap(sampleDistMatrix,
			clustering_distance_rows = sampleDists,
			clustering_distance_cols = sampleDists,
			col = colors,
			main='Heatmap of sample-to-sample distances using the rlog-transformed values.')

	# Poisson Distance
	library("PoiClaClu")
	poisd <- PoissonDistance(t(counts(dds)))

	samplePoisDistMatrix <- as.matrix( poisd$dd )
	rownames(samplePoisDistMatrix) <- rld$Treatment #paste( rld$dex, rld$cell, sep=" - " )
	colnames(samplePoisDistMatrix) <- NULL
	pheatmap(samplePoisDistMatrix,
			clustering_distance_rows = poisd$dd,
			clustering_distance_cols = poisd$dd,
			col = colors, 
			main='Heatmap of sample-to-sample distances using the Poisson Distance.')


	## PCA plot with ggplot:
	pcaData = plotPCA(rld, intgroup = c("Treatment","Species"), returnData=T)
	percentVar <- round(100 * attr(pcaData, "percentVar"))
	ggplot(pcaData, aes(x = PC1, y = PC2, color = Treatment, shape = Species)) +
	geom_point(size =4) + ggtitle("PCA plot using the rlog-transformed values") +
	xlab(paste0("PC1: ", percentVar[1], "% variance")) +
	ylab(paste0("PC2: ", percentVar[2], "% variance")) +
	coord_fixed()

	## MDS ggplot:
	mds <- as.data.frame(colData(rld))  %>%
			cbind(cmdscale(sampleDistMatrix))
	ggplot(mds, aes(x = `1`, y = `2`, color = Treatment, shape = Species)) +
		ggtitle("MDS plot using rlog-transformed values.") + 
	geom_point(size = 4) + coord_fixed()

	## MDS using PoissonDistance values:
	mdsPois <- as.data.frame(colData(dds)) %>%
	cbind(cmdscale(samplePoisDistMatrix))
	ggplot(mdsPois, aes(x = `1`, y = `2`, color = Treatment, shape = Species)) +
		ggtitle("MDS plot using PoissonDistance values.") + 
	geom_point(size = 4) + coord_fixed()


	### Discard potato and Undetermined samples for a few more plots:
	rld2 = rld[,rld$Species == 'tomato']

	## PCA plot:
	pcaData = plotPCA(rld2, intgroup = c("Treatment","Species"), returnData=T)
	percentVar <- round(100 * attr(pcaData, "percentVar"))
	ggplot(pcaData, aes(x = PC1, y = PC2, color = Treatment, shape = Species)) +
	geom_point(size =4) + ggtitle("PCA plot using the rlog-transformed values (tomato only)") +
	xlab(paste0("PC1: ", percentVar[1], "% variance")) +
	ylab(paste0("PC2: ", percentVar[2], "% variance")) +
	coord_fixed()

	## MDS ggplot:
	## MDS using PoissonDistance values:
	sampleDists <- dist(t(assay(rld2)))

	sampleDistMatrix <- as.matrix( sampleDists )
	rownames(sampleDistMatrix) <- rld2$Treatment
	colnames(sampleDistMatrix) <- NULL

	mds <- as.data.frame(colData(rld2))  %>%
			cbind(cmdscale(sampleDistMatrix))
	ggplot(mds, aes(x = `1`, y = `2`, color = Treatment, shape = Species)) +
		ggtitle("MDS plot using rlog-transformed values (tomato only).") + 
	geom_point(size = 4) + coord_fixed()


dev.off()

#### Assessment of the QA/QC plots  ######
"
-

"

### 

########### Differential Expression with DESeq2 #######


##### Do some data cleanup based on QA/QC and etc #######
## What are the top 10 expressed genes?
#library("AnnotationHub")
#Counts
# counts(dds)[order(rowSums(counts(dds)), decreasing=T)[1:30],]

#Percent
# 100*t(t(counts(dds)[order(rowSums(counts(dds)), decreasing=T)[1:30],])/colSums(counts(dds)))

pct = 100*t(t(counts(dds))/colSums(counts(dds)))
pct2 = pct[order(rowSums(pct), decreasing=T),]
pct2 = data.frame(pct2, txann[match(rownames(pct2), txann$gene_id), ])

#Write out some tables:
write.table(file = "top_30_highest_expressed_percent.tsv", data.frame(gene=rownames(pct), pct2)[1:30,], sep=',', row.names=F)
# write.table(file = "Percent_viral_reads_uninfected.tsv", data.frame(ViralTranscript=vids, pct2[vids,dds$SampleID[dds$Virus == 'None']]),
# 			sep=',', row.names=F)
# write.table(file = "Number_viral_reads_uninfected.tsv", data.frame(ViralTranscript=vids, counts(dds)[vids,dds$SampleID[dds$Virus == 'None']]),
# 			sep=',', row.names=F)

#Nothing removed:
#this is how to remove something: dds = dds[rownames(counts(dds)) != 'ENSDARG00000089382', ]

###### DESeq2 ######
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#standard-workflow

# Filter out low-count genes where at least 3 samples have to have > 10 count:
#keep = rowSums(counts(dds)) >= 10
keep = rowSums(counts(dds) > 40) >= 4
table(keep)
# keep
# FALSE  TRUE 
# 21021 13039 

dds2 = dds[keep,] 

# Plot some post-filter clusters:
pdf(file="QA-QC-filtered-plots.pdf", width=16, height=8)

	#MDS plots:
	d = dist(t(counts(dds2)))
	fit = cmdscale(d, eig=T, k=3)
	par(mfrow=c(1,2))
	plot(x=fit$points[,1], y=fit$points[,2], xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS", type='n',
			xlim=range(fit$points[,1]) * 1.1, ylim=range(fit$points[,2]) * 1.1)
	text(x=fit$points[,1], y=fit$points[,2], labels=row.names(fit$points), col = mycolors[dds2$Treatment], font=2)
	legend('bottomright', col=mycolors, legend=names(mycolors), pch=20, pt.cex=2)
	plot(x=fit$points[,1], y=fit$points[,3], xlab="Coordinate 1", ylab="Coordinate 3", main="Metric MDS", type='n',
			xlim=range(fit$points[,1]) * 1.1, ylim=range(fit$points[,3]) * 1.1)
	text(fit$points[,1], fit$points[,3], labels=row.names(fit$points), col = mycolors[dds2$Treatment], font=2)

	#Hclust:
	par(mfrow=c(1,2))
	hc = hclust(dist(t(counts(dds2))), method='ward.D2')
	plot(hc, main="Hclust with Euclidean distance and Ward.D2 clustering")
	hc = hclust(as.dist(1-abs(cor(counts(dds2)))), method='ward.D2')
	plot(hc, main="Hclust with 1-abs(correlation) distance and Ward.D2 clustering")
dev.off()

# Drop the "Undetermined" sample and reestimate parameters:
dds2 = dds2[,!(dds2$SampleID %in% "Undetermined")]


## The PSINA2.4 sample looks quite different from others in the group:
# This PSINA2-4 sample is sort of strange looking compared to the rest, run the analysis both with and w/o it.
#dds2 = dds2[,!(dds2$SampleID %in% "PSINA2-4")]



## Make sure that the "control" condition is set as such:
## This means that everything is + or - with respect to the None.0H.None condition:
dds2$Treatment = relevel(droplevels(dds2$Treatment), ref="AC")

## Start Differentialy expression:
dds2 = DESeq(dds2)  # estimate parameters

#system('mkdir -p 02-Plots')
all_contrasts = list(SINA1_vs_AC = c("Treatment", "SINA1","AC"),
					 SINA3_vs_AC = c("Treatment", "SINA3","AC"),
					 SINA4_vs_AC = c("Treatment", "SINA4","AC"),
					 SINA5_vs_AC = c("Treatment", "SINA5","AC"),
					 SINA6_vs_AC = c("Treatment", "SINA6","AC"),
					 CRCFDPho_vs_AC = c("Treatment", "CRCFDPho","AC"),
					 CRCFDPhe_vs_AC = c("Treatment", "CRCFDPhe","AC"),
					 PSINA2_vs_PtoR = c("Treatment", "PSINA2", "PtoR"),
					 PSINA3_vs_PtoR = c("Treatment", "PSINA3", "PtoR"),
					 PNAC1_vs_PtoR = c("Treatment", "PNAC1", "PtoR"),
					 BURP_vs_ACb = c("Treatment", "BURP", "ACb"),
					 SINA10E_vs_ACwt = c("Treatment", "SINA10E", "ACwt"),
					 BSDKO_vs_ACwt = c("Treatment", "BSDKO", "ACwt"),
					 RHA1B_vs_Desiree = c("Treatment", "RHA1B", "Desiree")					 
			)

all_results=c()

process_contrast = function(contrast, all_contrasts, dds2){
	print(contrast)
	system(paste('mkdir -p', contrast))
	outp = paste('./', contrast, '/', sep='')
	res = results(dds2, contrast=all_contrasts[[contrast]], alpha=.2, lfcThreshold=0)  # , lfcThreshold=1
	summary(res)
	capture.output(summary(res), file=paste("./", contrast, '/deg_results.txt', sep=''))

	#all_results = c(all_results, res)

	#This moderates the log2 fold change using Emperical Bayes:
	#  2019-03-20 there are now two ways to handle log2 fold change.
	#		https://support.bioconductor.org/p/110307/
	#		One approach (which is now listed in https://f1000research.com/articles/4-1070/v2) appears to be to estimate
	# 			abs(lfc) and threshold results with lfcThreshold=1 (which means threshold on genes which double in expression)
	#       The other approach (the old way) is to threshold on adjusted p-value and then shrink the lfc using a new "apeglm" strategy.
	#			The (better) estimates of lfc could then be used for thresholding as before.
	# 2019-03-21 and see this for s-values and a discussion of methods of shrinking lfc.
	#	http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
	#	https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#moreshrink especially this.
	#   It turns out that you can't use the recommended apeglm shrinkage strategy unless the comparison you are doing is in the list
	#	displayed by resultsNames(dds2). If you use special contrasts, then it won't work, and you have to change your model using
	#	relevel() so that the contrast you need IS listed. Only "normal" and "ashr" support arbitrary contrasts. It appears that 
	#	"ashr" and "apeglm" have roughly (?) the same quality of results, so "ashr" is used.

	print("resShrink...")
	#resShrink = lfcShrink(dds2, contrast=all_contrasts[[contrast]], type="apeglm")
	#coef = paste(all_contrasts[[contrast]][1], all_contrasts[[contrast]][3],	"vs", all_contrasts[[contrast]][2], sep='_')
	#resShrink = lfcShrink(dds2, coef= coef, type="apeglm")
	resShrink = lfcShrink(dds2, contrast=all_contrasts[[contrast]], type="ashr")

	# CITE: Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2. https://doi.org/10.1093/biostatistics/kxw041

	print("Plotting...")
	pdf(file=paste(outp, contrast, '.pdf', sep=''), width=11, height=11)
	#pdf(file=paste("./02-Plots/", contrast, '.pdf', sep=''), width=11, height=11)
		par(mfrow=c(1,2))
		plotMA(res, main="MA Plot, standard Log2 fold change")
		plotMA(resShrink, main="MA Plot, moderated Log2 fold change")

		## Gene clustering, top differentially expressed (by p-value)
		#topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 200)
		l = all_contrasts[[contrast]][c(2,3)]
		topDEGenes = rownames(res)[head(order(res$pvalue), 150)]
		mat = log2(counts(dds2, normalized=T)[topDEGenes, dds2$Treatment %in% l] + 1)
		#mat  <- mat - rowMeans(mat) # center
		anno <- as.data.frame(colData(dds2)[dds2$Treatment %in% l, c("Treatment", "Species")])
		pheatmap(mat, annotation_col = anno, 
				 main=paste(contrast, ': Top 150 DEGs by adjusted p-value, Log2 normalized count values', sep=''),
				 cluster_rows=F, fontsize_row=6)

		## Gene clustering, top differentially expressed (by log2FoldChange)
		topLFCGenes = rownames(resShrink)[head(order(abs(resShrink$log2FoldChange), decreasing=T), 150)]
		mat = log2(counts(dds2, normalized=T)[topLFCGenes, dds2$Treatment %in% l] + 1)
		#mat  <- mat - rowMeans(mat) # center
		anno <- as.data.frame(colData(dds2)[dds2$Treatment %in% l, c("Treatment", "Species")])
		pheatmap(mat, annotation_col = anno, 
			main=paste(contrast, ": Top 150 DEGs by moderated Log2FoldChange, Log2 normalized count values", sep=''),
			cluster_rows=F, fontsize_row=6)

		## Heat map of all DE genes at padj < .1:
		DEGenes = na.omit(rownames(res)[res$padj < .1])
		if(length(DEGenes) > 1){
			mat = log2(counts(dds2, normalized=T)[DEGenes, dds2$Treatment %in% l] + 1)
			#mat  <- mat - rowMeans(mat) # center
			anno <- as.data.frame(colData(dds2)[dds2$Treatment %in% l, c("Treatment", "Species")])
			pheatmap(mat, annotation_col = anno, main='All DEGs, padj < 0.1, Log2 normalized count values',
						show_rownames=F)
		}

	dev.off()
	# Order results:
	resOrdered = res[order(res$pvalue), ]

	# Now the results can be written to a report using the DESeqDataSet object.
	#source("https://bioconductor.org/biocLite.R")
	#biocLite("org.Hs.eg.db")
	print("Annotating")
	#columns(org.Hs.eg.db)

	de.idx = resOrdered$padj < .1 & !(is.na(resOrdered$padj))
	out = resOrdered[de.idx, ]
	out = data.frame(out, moderatedL2FC=resShrink[rownames(out), "log2FoldChange"])
	out = out[,c("baseMean","log2FoldChange","moderatedL2FC","lfcSE","stat","pvalue","padj")]
	out = data.frame(ITAG40=rownames(out), out)
	out = data.frame(out, counts(dds2, normalized=T)[rownames(out), dds2$Treatment %in% l])
	out = data.frame(out, txann[rownames(out), ])

	#what columns are available?
	#columns(annotationdb)

	#ann = AnnotationDbi::select(org.Hs.eg.db, keys=rownames(out), columns=c("GENENAME","SYMBOL","ENTREZID"), keytype='ENSEMBL')
	#ann = ann[!duplicated(ann$ENSEMBL), ]

	#out = data.frame(out, ann[match(rownames(out), ann$ENSEMBL), c("ENSEMBL", "GENENAME","SYMBOL","ENTREZID")])
	#out = data.frame(out, ann[match(rownames(out), ann$ENSEMBL), c("GENENAME","SYMBOL","ENTREZID")])

	#idx = which(rownames(out) %in% names(vids.desc))
	#out$GENENAME[idx] = vids.desc[rownames(out)[idx]]

	write.table(out, file=paste(outp, contrast, "_DEGs.tsv", sep=''), row.names=F, sep='\t')


	##### topGO ######
	# See https://www.bioconductor.org/packages/devel/bioc/vignettes/topGO/inst/doc/topGO.pdf
	# https://ucdavis-bioinformatics-training.github.io/2019_August_UCD_mRNAseq_Workshop/differential_expression/enrichment
	# Using the elim + weight algorithm from Improved scoring of functional groups from gene expression data by decorrelating GO graph structure
	# Alexa et al. 2006.
	degs = rownames(out)
	geneList = factor(as.integer(names(geneID2GO) %in% degs))
	names(geneList) = names(geneID2GO)
	
	GOdata = new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO=geneID2GO)
	resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
	tab <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 200)
	write.table(tab, file=paste(outp, contrast, "_GO_MolecularFunction.tsv", sep=''), sep='\t', row.names=F)

	GOdata = new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO=geneID2GO)
	resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
	tab <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 200)
	write.table(tab, file=paste(outp, contrast, "_GO_BiologicalProcess.tsv", sep=''), sep='\t', row.names=F)

	GOdata = new("topGOdata", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO=geneID2GO)
	resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
	tab <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 200)
	write.table(tab, file=paste(outp, contrast, "_GO_CellularComponent.tsv", sep=''), sep='\t', row.names=F)

}

mclapply(names(all_contrasts), process_contrast, all_contrasts, dds2, mc.cores=10)

# for(contrast in names(all_contrasts)){
# 	process_contrast(contrast, all_contrasts, dds2)
# }



writeLines(capture.output(sessionInfo()), "sessionInfo.txt")












# ## WGCNA ###
# # cluster the normalized read counts stored in "nc"
# # https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-auto.pdf

# library(WGCNA)
# options(stringsAsFactors = F)
# enableWGCNAThreads()

# #Note, WGCNA expects row to be a sample
# # Get the top 5000 most variable genes:
# m = counts(dds2, normalized=T)
# v = apply(m, 1, var)
# o = order(v, decreasing=T)
# #nc.subset = anntable[keep.idx, c("GFPneg1", "GFPneg2", "GFPneg2_2", "GFPneg5_2", "GFPpos1", "GFPpos2", "GFPpos2_2", "GFPpos5_2")]

# m.subset = m[o,][1:5000,]

# m.t = t(m.subset) # WGCNA wants rows as samples

# powers = c(seq(from=1, to=30, by=4), seq(from=30, to=50, by=2))

# # call the network topology analysis function
# sft = pickSoftThreshold(m.t, powerVector = powers, verbose = 5)


# png(file='./reports/pick_WGCNA_params.png', width=1000, height=600)
# #sizeGrWindow(9, 5)
# par(mfrow = c(1,2));
# cex1 = 0.9;
# # Scale-free topology fit index as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.90,col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# dev.off()


# ## A power of 9 seems like a reasonable value based on mean connectivity.
# # Construct gene network and identify modules in one function call:
# net = blockwiseModules(m.t, power = 9,
#                        TOMType = "signed", minModuleSize = 15,
#                        reassignThreshold = 0, mergeCutHeight = 0.15,
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#                        saveTOMs = TRUE,
#                        saveTOMFileBase = "zebrafishROD_TOM",
#                        verbose = 3)
# table(net$colors)


# png(file="./reports/WGCNA_cluster_dendrogram.png", width=2000, height=1800)
# # open a graphics window
# #sizeGrWindow(12, 9)
# # Convert labels to colors for plotting
# mergedColors = labels2colors(net$colors)
# # Plot the dendrogram and the module colors underneath
# plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
# "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
# dev.off()



# ## Make per-cluster centered and scaled expression plots:
# ucolors = unique(net$colors)
# pdf(file='./reports/WGCNA_per_cluster_expression_plots.pdf', width=10, height=6)
# for(color in ucolors){
#   par(mar=c(6,4,4,2))
#   matplot(scale(m.t[,net$colors == color]), type='l', xaxt = 'n', ylab="Normalized read counts, centered and scaled.",
#     main=paste("Cluster", labels2colors(color)))
#   par(las=2)
#   axis(side=1, at=1:8, lab=rownames(m.t))
# }
# dev.off()
