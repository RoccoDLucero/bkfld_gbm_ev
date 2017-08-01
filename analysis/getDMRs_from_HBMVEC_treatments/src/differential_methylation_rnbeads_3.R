#This Rscript is written to perform the differential methylation analysis of the
#Breakefield 450K methylation data sets. These experiments are designed so that we
#can discern GBM8-EV treatment-specific methylation effects on HBMVECs in vitro

#This script calls MapPRESSTOenh-prom_sym.R to prepare the custom annotation
#regions of enhancers mapped to the promoters based on Andersson et. al.
# Analysis of FANTOM5 CAGE data

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
report.dir = file.path(analysis.dir, "GetDMRs_9comps_SVAbe_July19_2016")
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
source('./MapPRESSTOenh_prom_sym.R') #creates object my.mapped.ehnancer.annots

# Create RnBeads annotations by providing data frames for each set of regions
for(enh.prof in names(my.mapped.ehnancer.annots)){
    rnb.set.annotation(type = enh.prof, description = enh.prof,
                       regions = my.mapped.ehnancer.annots[[enh.prof]][c(1:3,6)],
                       assembly="hg19")
}

#rnb.get.annotation(type = "T_cell.enh")

# Set the options to include the enhancer annotations in region.types
rnb.options(region.types=c("genes", "promoters", names(my.mapped.ehnancer.annots)))
rnb.getOption('region.types')

################################################################################
#Methylation Data loading and pre-formatting
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

################################################################################
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
##
##Global options
rnb.initialize.reports(report.dir)
#Turn the data into an RnBeadSet object:
my.rnb.set <- RnBeadSet(pheno= my.pheno,
                         betas= betas,
                         probes = rownames(betas),
                         useff=FALSE,
                         p.values = p.vals, 
                         platform = '450k')


#Remove probes with a p-value > 0.01 in any sample
any.bad.p.val = apply(dpval(my.rnb.set)>0.01, 1, any)
my.rnb.set = remove.sites(my.rnb.set, any.bad.p.val)

#Make sure the data appears appropriately to RnBeads
rnb.sample.groups(annotations = my.pheno, columns = NULL, columns.pairs = NULL,
                  min.group.size = rnb.getOption("min.group.size"),
                  max.group.count = rnb.getOption("max.group.count"))

rnb.sample.replicates(my.rnb.set, "Replicate")

#Quality Control
rnb.options("qc"=TRUE)
rnb.options("qc.boxplots"=TRUE, "qc.barplots"=TRUE, "qc.negative.boxplot"=FALSE,
            "qc.snp.distances"=TRUE, "qc.snp.boxplot"=T, "qc.snp.barplot"=FALSE,
            "qc.sample.batch.size"=50, "qc.coverage.plots"=FALSE,
            "qc.coverage.threshold.plot"=1:10,"qc.coverage.histograms"=FALSE,
            "qc.coverage.violins"=T)

if (rnb.getOption("qc")) {
   rnb.run.qc(my.rnb.set, report.dir)
}

#Preprocessing - Essential module.
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

#Covariate Inference
rnb.options("inference" = TRUE)
rnb.options("inference.targets.sva" = c('Group','Replicate','Celltype'),
            "inference.sva.num.method" = "be",
            "inference.age.prediction" = FALSE )

if(rnb.getOption("inference")){
  my.rnb.set <- rnb.run.inference(my.rnb.set,report.dir)$rnb.set
  saveRDS(my.rnb.set, file = paste(report.dir,'/my.rnb.set.inf',sep = ''))
}

#Exploratory Analysis
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
rnb.options("differential"= TRUE)
rnb.options("differential.site.test.method" = "limma",   # alternatively, the "ttest" method
            "differential.permutations" = nPerms,
            "analyze.sites" = TRUE,
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
full.run = F
if(rnb.getOption("differential")){
  if(!full.run){    
      for(cur.comp in names(my.dmr.list)){
          print(cur.comp)
          #Produce the diffmeth object per comparison for all region types
          my.dmrs <- rnb.execute.computeDiffMeth(x = my.rnb.set,
                                                 pheno.cols = cur.comp,
                                                 
                                                 region.types = rnb.region.types.for.analysis(my.rnb.set),
                                                 n.perm = nPerms)
          #For each diffmeth object
          
          my.dmr.list[[cur.comp]] <- my.dmrs
      }
  } else {     
    #my.dmrs <- rnb.run.differential(my.rnb.set, report.dir,show.report = T)
  }
  saveRDS(my.dmr.list, file = paste(report.dir,'/my.dmrs.noPerms.final',sep = ''))
}

#This script will get the coordinates of annotation regions
#covered on the 450K array
source(file = './FinalAnnotDMRtables.R', verbose = T)


####END####
