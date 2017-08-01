#This Rscript is written to perform the differential methylation analysis of the
#Breakefield 450K methylation data sets. These experiments are designed so that we
#can discern GBM8-EV treatment-specific methylation effects on HBMVECs in vitro

#The most important output for this script is the set of region-based diffmeth.csv tables
#These can be processed later for example GO enrichment analysis can be performed

source("https://bioconductor.org/biocLite.R")
library(RnBeads)
library(hexbin)
library(wordcloud)

################################################################################
#Set run-specific parameters
raw.data.dir = "../../data/raw/Breakefield_HBMVEC_450k/"
analysis.dir = "./RnBeads/analysis"
report.dir = file.path(analysis.dir, "GetDMRs_9comps_July2016")
##sample.annotation = file.path(raw.data.dir,"sample_annotation.txt")

rnb.options("analysis.name"="Breakefield450k",
            "assembly"="hg19",
            "disk.dump.big.matrices"=FALSE,
            filtering.sex.chromosomes.removal=TRUE,
            import.table.separator="\t",
            identifiers.column="Sample_ID")

################################################################################
#Import any custom annotation regions over which to compute differntial methyaltion:
#Here we manually import bed files from 'http://enhancer.binf.ku.dk/presets/'
#Download using the web utility [do not copy text from webpage]
#
#'Negative control' cell type enhancer regions for our model:
#(We expect fewer DMRs here with GBM8-EV treated HBMVECs)
#   skin fibroblast
#   T-Cell
#   skeletal muscle cell
#'Query' cell type enhancer regions for our model:
#(We suspect some overlap in DMRs over these regions with GBM8-EV treated HBMVECs)
#   astrocyte
#   blood vessel endothelial cell
###'CL-' indicates cell-type active enhancer set
my.enhancer.profile.files <- list.files(path = './input', pattern = 'CL-')

my.profile.names <- gsub(pattern = '_differentially_expressed_enhancers.bed',
                         replacement = '',x = my.enhancer.profile.files)
my.profile.names <- gsub(pattern = '^...........',
                         replacement = '',x = my.profile.names)

my.enhancer.profiles = vector("list", length(my.profile.names))
names(my.enhancer.profiles) = paste(my.profile.names,".enh",sep = '')
for(ct in my.enhancer.profile.files){
    profile <- suppressWarnings(read.table( paste('./input/',ct ,sep = ''),
                                            sep = '\t', header = F, skip = 1))
    profile <- profile[,1:3]
    colnames(profile) = c("chromosome", "start", "end")
    my.enhancer.profiles[[which(my.enhancer.profile.files==ct)]] = profile
}

# Create RnBeads annotations by providing data frames for each set of regions
for(enh.prof in names(my.enhancer.profiles)){
    rnb.set.annotation(enh.prof,
                       my.enhancer.profiles[[enh.prof]],
                       assembly="hg19")
}


# Set the options to include the enhancer annotations in region.types
rnb.options(region.types=c("genes", "promoters", names(my.enhancer.profiles)))
rnb.getOption('region.types')

################################################################################
#Methjyaltion Data loading and pre-formatting
#I created the 'bkfld450k.rda' Rdata object for faster loading times
breakefield450K = readRDS(paste(raw.data.dir,'bkfld450k.rda',sep = '')) 
colnames(breakefield450K) = gsub("_zappulliMGH1_[ABCDEF][1234]", "", colnames(breakefield450K))
colnames(breakefield450K) = gsub("..", ".", colnames(breakefield450K),fixed = T)
breakefield450K = breakefield450K[,order(colnames(breakefield450K))]

betas = breakefield450K[,c(as.numeric(grep("AVG_Beta", colnames(breakefield450K))))]
rownames(betas) = breakefield450K$TargetID
p.vals = breakefield450K[,c(as.numeric(grep("Detection.Pval", colnames(breakefield450K))))]

################################################################################
#Import the comparisons table
my.pheno <- read.csv(file = './input/my.pheno.txt',header = T,sep = ',',quote = '',na.strings = 'NA')

#set each object to the proper type for inclusion in rnbeadset object#
my.pheno = as.data.frame(my.pheno)
betas = as.matrix(betas)
p.vals = as.matrix(p.vals)
dim(my.pheno)
#Select the desired subset of samples for this analysis
#In this case we want to make all comparisons
my.range = c(1:24)
my.pheno = my.pheno[my.range,]
betas = betas[,my.range]
p.vals = p.vals[,my.range]

my.diff.comps <- colnames(my.pheno)[5:13]
#############################################################################################################
#############################################################################################################

#############################################################################################################
#############################################################################################################
# Global options
rnb.initialize.reports(report.dir)
#Turn the data into an RnBeadSet object:
my.rnb.set <- RnBeadSet(pheno= my.pheno,
                         betas= betas,
                         probes = rownames(betas),
                         useff=FALSE,
                         p.values = p.vals, 
                         platform = '450k')

#Remove   probes   with   a   p-value   >   0.01   in   any   sample
any.bad.p.val = apply(dpval(my.rnb.set)>0.01, 1, any)
my.rnb.set = remove.sites(my.rnb.set, any.bad.p.val)

#Make sure the data appears appropriately to RnBeads
rnb.sample.groups(annotations = my.pheno, columns = NULL, columns.pairs = NULL,
                  min.group.size = rnb.getOption("min.group.size"),
                  max.group.count = rnb.getOption("max.group.count"))

rnb.sample.replicates(my.rnb.set, "Replicate")

# Quality Control - Skip if you have no p-values (usually if no p-values are provided, QC has already been done.)
rnb.options("qc"=TRUE)
rnb.options("qc.boxplots"=TRUE, "qc.barplots"=TRUE, "qc.negative.boxplot"=FALSE,
            "qc.snp.distances"=TRUE, "qc.snp.boxplot"=T, "qc.snp.barplot"=FALSE,
            "qc.sample.batch.size"=50, "qc.coverage.plots"=FALSE,
            "qc.coverage.threshold.plot"=1:10,"qc.coverage.histograms"=FALSE,
            "qc.coverage.violins"=T)

if (rnb.getOption("qc")) {
  rnb.run.qc(my.rnb.set, report.dir)
}

# Preprocessing - Essential module.
rnb.options("preprocessing"=TRUE)
rnb.options("filtering.whitelist"=NULL, "filtering.blacklist"=NULL, 
            "filtering.context.removal"=c("CC","CAG","CAH","CTG","CTH","Other"),
            "filtering.snp"="any", #Choose from: c("no", "3", "5", "any", "yes") 
            "filtering.cross.reactive" = TRUE,
            "filtering.greedycut" = FALSE,
            "filtering.greedycut.pvalue.threshold"=0.05,
            "filtering.sex.chromosomes.removal"=FALSE,
            "filtering.missing.value.quantile"=1,
            "filtering.coverage.threshold"=5,
            "filtering.low.coverage.masking"=FALSE, 
            "filtering.high.coverage.outliers"=FALSE,
            "filtering.deviation.threshold"=0
            )

if(rnb.getOption("preprocessing")){
  my.rnb.set <- rnb.run.preprocessing(rnb.set = my.rnb.set, report.dir)$rnb.set
  saveRDS(my.rnb.set, file = paste(report.dir,'/my.rnb.set.pre',sep = ''))
}

#Covariate Inference - They call this "optional" but I call it ESSENTIAL! :-)
rnb.options("inference" = TRUE)
rnb.options("inference.targets.sva" = c('Group','Replicate','Celltype'),
            "inference.sva.num.method" = "leek",
            "inference.age.prediction" = FALSE )

if(rnb.getOption("inference")){
  my.rnb.set <- rnb.run.inference(my.rnb.set,report.dir)$rnb.set
  saveRDS(my.rnb.set, file = paste(report.dir,'/my.rnb.set.inf',sep = ''))
}

# Exploratory Analysis  -- I usually skip this module, because it takes a while, but it's certainly useful for discovery!
rnb.options("exploratory"= F)
rnb.options("exploratory.columns"=NULL, "exploratory.top.dimensions"=0, 
            "exploratory.principal.components"= 8, "exploratory.correlation.pvalue.threshold"=0.01,
            "exploratory.correlation.permutations"=1000, "exploratory.correlation.qc"=TRUE,
            "exploratory.beta.distribution"=TRUE, "exploratory.intersample"=TRUE, 
            "exploratory.deviation.plots"=NULL, "exploratory.clustering"="all",
            "exploratory.clustering.top.sites"= 1000, "exploratory.clustering.heatmaps.pdf"=FALSE,
            "exploratory.region.profiles"=NULL, "exploratory.region.profiles"=NULL,
            "exploratory.gene.symbols"=NULL, "exploratory.custom.loci.bed"=NULL)

if(rnb.getOption("exploratory")){
  my.rnb.set <- rnb.run.exploratory(my.rnb.set,report.dir)$rnb.set
  saveRDS(my.rnb.set, file = paste(report.dir,'/my.rnb.set.expl',sep = ''))
}

nPerms = 0
# Differential Methylation
rnb.options("differential"= T)
rnb.options("differential.site.test.method" = "limma",   # alternatively, the "ttest" method
            "differential.permutations" = nPerms,
            "analyze.sites" = FALSE,
            "differential.comparison.columns" = my.diff.comps, 
            "differential.comparison.columns.all.pairwise" = NULL,
            "covariate.adjustment.columns" = c('Replicate','Group'),
            "columns.pairing" = NULL, 
            "differential.adjustment.sva" = TRUE,  # this is a popular correction/detection method for batch effects
            "differential.adjustment.celltype" = FALSE,   # set to TRUE if you are doing epigenetic deconvolution
            "differential.enrichment" = FALSE, #This is done in a separate script
            "export.to.bed" = FALSE, "export.types" = NULL)

my.dmr.list <- vector('list', length(rnb.getOption("differential.comparison.columns")))
names(my.dmr.list) <- rnb.getOption("differential.comparison.columns")
if(rnb.getOption("differential")){
  for(cur.comp in names(my.dmr.list)){
      print(cur.comp)
      my.dmrs <- rnb.execute.computeDiffMeth(x = my.rnb.set,
                                             pheno.cols = cur.comp,
                                             region.types = rnb.region.types.for.analysis(my.rnb.set),
                                             n.perm = nPerms)
      my.dmr.list[[cur.comp]] <- my.dmrs
#      my.dmrs <- rnb.run.differential(my.rnb.set, report.dir)
  }
  saveRDS(my.dmr.list, file = paste(report.dir,'/my.dmrs.final',sep = ''))
}
  
####END####
