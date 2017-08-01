# RnBeads implementation Skeleton
require(RnBeads)
require(hexbin)
require(wordcloud)

# Global options
rnb.options("analysis.name"="Garovic450k", "assembly"="hg19",
            "disk.dump.big.matrices"=FALSE)
            #,"region.types"=c("genes", "promoters"))  # "tiling"     "genes"      "promoters"  "cpgislands"

# Data Import 
setwd('~/Downloads/Breakefield/raw.data')
load('initialized.RData')
#breakefield450k <- read.table("FinalReport2.txt", header=TRUE, sep="\t", quote="")
#data450k <- breakefield450k[,c(1, as.numeric(grep(".*.Signal_A", colnames(breakefield450k))), 
#                               as.numeric(grep(".*.Signal_B", colnames(breakefield450k))),
#                               as.numeric(grep(".*.Detection.Pval", colnames(breakefield450k))))]
#rm(breakefield450k)
#rownames(data450k) <- data450k[,1]
#data450k <- data450k[,-1]
#colnames(data450k) <- gsub("zappulli", "", colnames(data450k))
#data450k <- data450k[,order(colnames(data450k))]
#data450k[1,]
#betas <- data450k[,seq(3,dim(data450k)[2],3)]/(data450k[,seq(2,dim(data450k)[2],3)]+data450k[,seq(3,dim(data450k)[2],3)])
#rm(data450k)

# Uncomment below if you need to make your own sample_annotation.txt file. You only need to execute this code once.
# Human brain microvascular endothelial cells (HMVECs)
# Supernatant is fluid that remains after centrifugation
#sampleNames <- colnames(betas)                              
#sampleNames[grep("EBM.*", colnames(betas), perl=TRUE)] <- "Endothelial.Basal.Medium"                # HBMVEC + EBS   
#sampleNames[grep("EGM.*", colnames(betas), perl=TRUE)] <- "Endothelial.Growth.Medium"               # HBMVEC + EBS + GFs
#sampleNames[grep("EVpellet.*", colnames(betas), perl=TRUE)] <- "Extracellular.Vesicle"              # HBMVEC + EBS + EV pellet
#sampleNames[grep("EVsuper.*", colnames(betas), perl=TRUE)] <- "Extracellular.Vesicle.Supernatant"   # HBMVEC + EBS + EV pellet + surnatant
#sampleNames[grep("GBM8EBM.*", colnames(betas), perl=TRUE)] <- "GBM8.in.EBM"                         # HBMVEC + GBM8 + supernatant
#sampleNames[grep("GBM8NB.*", colnames(betas), perl=TRUE)] <- "GBM8.in.NB"                           # HBMVEC + GBM8  
#sampleNames[grep("UCMpellet.*", colnames(betas), perl=TRUE)] <- "Unconditioned.Medium"              # GBM8 + normal medium
#sampleNames[grep("UCMsuper.*", colnames(betas), perl=TRUE)] <- "Unconditioned.Medium.Supernatant"   # GBM8 + without growth factors
#cell.lane <- gsub("^.*_MGH1_", "", colnames(betas))
#cell.lane <- gsub(".Signal_B", "", cell.lane)
#temp <- cbind(colnames(betas), sampleNames)
#colnames(temp) <- c("Sample_ID", "Group")
#write.table(temp, "sample_annotation.txt", quote=FALSE, eol="\r\n", row.names=FALSE, col.names=TRUE, sep="\t")
#rm(temp, sampleNames)

#home <- "~/Downloads/Breakefield/"
#home <- "C:/Users/Lillian/OneDrive/Lab/Breakefield"
#dirs <- list(data=paste(home,"raw.data/",sep=""),
#             annotation=paste(home,"raw.data/sample_annotation.txt",sep=""),
#             analysis=paste(home,"results/",sep=""),
#             reports=paste(home,"results/reports",sep=""))
#rnb.getOption("import.table.separator")
#rnb.options(filtering.sex.chromosomes.removal=TRUE,
#            import.table.separator="\t",  identifiers.column="Sample_ID")
#pheno <- read.table(dirs$annotation, sep="\t", header=TRUE)
#data.source <- RnBeadSet(pheno=pheno, betas=as.matrix(betas), probes = rownames(betas), useff=FALSE) # region.types=c("genes", "promoters"),
#rnb.initialize.reports(dirs$reports)
#rnb.set <- data.source
#rm(data.source, betas, pheno, home)
#save.image("initialized.RData")


# Quality Control - Skip if you have no p-values (usually if no p-values are provided, QC has already been done.)
rnb.options("qc"=TRUE)
rnb.options("qc.boxplots"=TRUE, "qc.barplots"=TRUE, "qc.negative.boxplot"=FALSE,
            "qc.snp.distances"=TRUE, "qc.snp.boxplot"=FALSE, "qc.snp.barplot"=FALSE,
            "qc.sample.batch.size"=50, "qc.coverage.plots"=FALSE, "qc.coverage.threshold.plot"=1:10,
            "qc.coverage.histograms"=FALSE, "qc.coverage.violins"=FALSE)
if (rnb.getOption("qc")) {
  rnb.run.qc(rnb.set, dirs$reports)
}

# Preprocessing - Essential module.
rnb.options("preprocessing"=TRUE)
rnb.options("filtering.whitelist"=NULL, "filtering.blacklist"=NULL, 
            "filtering.context.removal"=c("CC", "CAG", "CAH", "CTG", "CTH", "Other"),
            "filtering.snp"="any", # Choose from any of following options: c("no", "3", "5", "any", "yes")
            "filtering.cross.reactive" = TRUE,
            "filtering.greedycut" = FALSE, "filtering.greedycut.pvalue.threshold"=0.05,
            "filtering.sex.chromosomes.removal"=FALSE, "filtering.missing.value.quantile"=1,
            "filtering.coverage.threshold"=5, "filtering.low.coverage.masking"=FALSE, 
            "filtering.high.coverage.outliers"=FALSE, "filtering.deviation.threshold"=0
            )
if (rnb.getOption("preprocessing")) {
  rnb.set <- rnb.run.preprocessing(rnb.set, dirs$reports)$rnb.set
}
save.image("preprocessed.RData")


# Tracks and Tables - This is useful if you are interfacing RnBeads with, say, the Genome Browser and need different file formats of results
#                     like bigwig, bed, etc.
#rnb.run.tnt(rnb.set, dirs$reports)


# Covariate Inference - They call this "optional" but I call it ESSENTIAL! :-)
rnb.options("inference"=TRUE)
rnb.options("inference.targets.sva"=character(), "inference.reference.methylome.column"=character(),
            "inference.max.cell.type.markers"=10000, "inference.top.cell.type.markers"=500,
            "inference.sva.num.method"="leek")
if (rnb.getOption("inference")) {
  rnb.set <- rnb.run.inference(rnb.set, dirs$reports)$rnb.set
}
save.image("inference.RData")



# Exploratory Analysis  -- I usually skip this module, because it takes awhile, but it's certainly useful for discovery!
rnb.options("exploratory"=TRUE)
rnb.options("exploratory.columns"=NULL, "exploratory.top.dimensions"=0, 
            "exploratory.principal.components"=8, "exploratory.correlation.pvalue.threshold"=0.01,
            "exploratory.correlation.permutations"=10000, "exploratory.correlation.qc"=TRUE,
            "exploratory.beta.distribution"=TRUE, "exploratory.intersample"=TRUE, 
            "exploratory.deviation.plots"=NULL, "exploratory.clustering"="all",
            "exploratory.clustering.top.sites"=1000, "exploratory.clustering.heatmaps.pdf"=FALSE,
            "exploratory.region.profiles"=NULL, "exploratory.region.profiles"=NULL,
            "exploratory.gene.symbols"=NULL, "exploratory.custom.loci.bed"=NULL)
if (rnb.getOption("exploratory")) {
  rnb.set <- rnb.run.exploratory(rnb.set, dirs$reports)$rnb.set
}

# Differential Methylation
rnb.options("differential"=TRUE)
rnb.options("differential.site.test.method"="limma",   # alternatively, the "ttest" method
            "differential.permutations"=0,
            "differential.comparison.columns"="Status", 
            "differential.comparison.columns.all.pairwise"="Status",
            "covariate.adjustment.columns"=NULL, "columns.pairing"=NULL, 
            "differential.adjustment.sva"=FALSE,  # this is a popular correction/detection method for batch effects
            "differential.adjustment.celltype"=FALSE,   # set to TRUE if you are doing epigenetic deconvolution
            "differential.enrichment"=TRUE)
if (rnb.getOption("differential")) {
  rnb.set <- rnb.run.differential(rnb.set, dirs$reports)
}
save.image("diffMeth.RData")


