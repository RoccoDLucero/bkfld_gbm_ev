library(RnBeads)
library(biomaRt)
library(org.Hs.eg.db)
library(annotate)
source(file = './processPRESSTOmapping.R')

#This needs to be called in the environment of the GETDMRsJuly_2016 (or later) script##
#make sure my.enhancer profiles list objedct is available
#Rnbeads hg19 annotations for genes and promoters
#import the list of DMRs over enhancer regions
if(!exists(x = 'my.dmr.list')){
    my.dmr.list <- readRDS(file = paste(report.dir,'/my.dmrs.final',sep = ''))
}

colnames(my.dmr.list$EBM.EGM@regions$T_cell.enh$`EBM vs. EGM (based on EBM.EGM)`)
quantile(my.dmr.list$EBM.EGM@regions$genes$`EBM vs. EGM (based on EBM.EGM)`[,"diffmeth.p.adj.fdr"])

################################################################
################################################################
#GENE and PROMOTER regions can be processed readily 
#To perform GO analysis on ENHANCER regions we must choose some
# way to assign gene symbols to enhancers. The Pressto data 
# offers one such mapping based on cage tag correlations
#
###############################################################
###############################################################

################################################################################
#Import enhancer annotation regions over which to compute differntial methyaltion:
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
my.enhancer.profile.files <- list.files(path = './input', pattern = 'CL-')
my.profile.names <- gsub(pattern = '_differentially_expressed_enhancers.bed',
                         replacement = '',x = my.enhancer.profile.files)
my.profile.names <- gsub(pattern = '^...........',
                         replacement = '',x = my.profile.names)

my.enhancer.profiles = vector("list", length(my.profile.names))
names(my.enhancer.profiles) = paste(my.profile.names,".enh",sep = '')
for(ct in my.enhancer.profile.files){
    profile <- suppressWarnings(read.table( paste('./input/',ct ,sep = ''), sep = '\t', header = F, skip = 1))
    profile <- profile[,1:3]
    colnames(profile) = c("chromosome", "start", "end")
    my.enhancer.profiles[[which(my.enhancer.profile.files==ct)]] = profile
}

################################################################################
#Import the PRESSTO enhancer promoter mappings, convert to GRanges 
enh.PresstoProm.assoc <- My.Prepare.Enh.Prom.Assoc.Data() #Function from processPRESSTOmapping.R
GR.enh.prom.assoc <- makeGRangesFromDataFrame(enh.PresstoProm.assoc,keep.extra.columns = T)

#For cases in which our universe of terms in the enrichment analyses
#will be all of the genes under promoters associating with an enhancer in the PRESSTO data
#select the gene symbols from all of the annotated enhancers

#Make a function that does something like the following, but keep
#thinking about the most informative way to do this analysis
#
#make all cell/tissue specific enhancer annotation sets into GRANGES objects
#subset the enhancer promoter map according to these regions
#perform GO on the genes associated with these regions 
#to look for overall function of that cell-type specific set
#using all gene symbols in the map as the universe
#THEN...
#Determine if cell-type specific active enhancer associated genes reflect an obvious
#theme associated with that cell type

#NEXT...
#Incorporate the DMR data for a given comparison into the analysis, by subsetting
#further based on regions that gain/lose methylation
#perform GO on these subsets with all cell-type associated gene symbols
#as the universe and finally all symbols as the universe.
#Then..
#Determine whether the genes associated with differential enhancer methylation reflect
#some functional themes
#
#Are the DMRs more prevalent in a given cell-type specific enhancer set?
#I hypothsize that DMRs will prevail in endothelial and Astrocyte  active enhancers
#
#Produce a report that breaks down DMR-gene associations
#by listing all genes significantly affected
#by listing all biological processes significantly affected
#Focus on the Positive control and GBM8EV treated groups
#focus on the enhancer sets astrocyte and blood-vessel endothelial cell
#
GR.enh.prom.assoc.univ <- GR.enh.prom.assoc[,'symbol']$symbol
GR.bv.enhancers = makeGRangesFromDataFrame(bv.endo.enhancers)

assoc.bv.enh.prom = subsetByOverlaps(GR.enh.refseqProm.assoc, GR.bv.enhancers)
assoc.bv.enh.prom.query = assoc.bv.enh.prom[,'symbol']$symbol




my.get.GO.gene.sets <- function(query, univ){
    mart<- useDataset(dataset = "hsapiens_gene_ensembl" , useMart("ensembl"))

    my.query <- getBM(filters= "hgnc_symbol",
                      attributes= c("ensembl_gene_id", "entrezgene", "description"),
                      values= query, mart= mart)

    my.univ <- getBM(filters= "hgnc_symbol",
                     attributes= c("ensembl_gene_id", "entrezgene", "description"),
                     values = univ, mart= mart)
    
    return(list(query = my.query, univ = my.univ))
}


az = performGOenrichment.diffMeth.entrez(gids = na.omit(ax$entrezgene), uids = na.omit(ay$entrezgene) ,ontology = "BP",assembly = 'hg19')
uz = performGOenrichment.diffMeth.entrez(gids = na.omit(ux$entrezgene), uids = na.omit(uy$entrezgene) ,ontology = "BP",assembly = 'hg19')
wz = performGOenrichment.diffMeth.entrez(gids = na.omit(wx$entrezgene), uids = na.omit(wy$entrezgene) ,ontology = "BP",assembly = 'hg19')

tt = summary(az)
ttt = summary(uz)
summary(wz)
setdiff(tt$Term,ttt$Term)
intersect(tt$Term,ttt$Term)
length(ttt$Term)
na.omit(ax)



##############################################
#
############################################################################
###########################################################################
############################################################################


my.dat.dir <- './RnBeads/analysis/noGBM8_noNegCtrl_w_astro_endo_enhancers/differential_methylation_data/'
my.diff.tables <- list.files(my.dat.dir)
my.diff.tables <- my.diff.tables[grep("diffMeth",my.diff.tables)]
#This may also be accessed via the RnBeads diff meth list my.dmrs
minRank <- 100

i <- 6
is.enhancer = T 
cur.diff.table <- read.csv(paste(my.dat.dir,my.diff.tables[i],sep = ''))


#Display the current comparison and region information for the 
#Differential methylation table
meanCols =  grep('mean.mean.',colnames(cur.diff.table))[c(1,2)]
smpNames <- c( gsub(pattern = 'mean.mean.', replacement = '', x = colnames(cur.diff.table[meanCols[1]])),
                   gsub(pattern = 'mean.mean.', replacement = '', x = colnames(cur.diff.table[meanCols[2]])) )
regNames <- (gsub('\\.csv','',gsub(pattern = 'diffMethTable_region_cmp._', replacement = '', x = my.diff.tables[i])))
writeLines(paste(regNames,':: ',paste(smpNames[1],'vs',smpNames[2]),sep = ''))

#Separate the table into tables of hypo and hypermethylated regions
cur.hypo.table <- cur.diff.table[cur.diff.table$mean.mean.diff < 0,]
cur.hyper.table <- cur.diff.table[cur.diff.table$mean.mean.diff > 0,]

#Rank the regions according to 'combined rank'
cur.hypo.table <- cur.hypo.table[order(x = cur.hypo.table$combinedRank,decreasing = F,na.last = T),]
cur.hyper.table <- cur.hyper.table[order(x = cur.hyper.table$combinedRank,decreasing = F,na.last = T),]

#Select just the best n ranked regions
cur.hypo.table <- cur.hypo.table[1:minRank,]
cur.hyper.table <- cur.hyper.table[1:minRank,]

#Print summary statistics for probe coverage over region
table(cur.hypo.table$num.sites)
table(cur.hyper.table$num.sites)

cur.dmr.table = list(hypo = cur.hypo.table, hyper = cur.hyper.table)

methDirection <- c("gained methylation","lost methylation")
j <- 1
ensembl_query <- cur.dmr.table[[j]]$id
ensembl_univ <- cur.diff.table$id

#This is for when we use enhancers. We need to get the association with
#promoters based on the overlap of regions listed in the diffmeth output of RnBeads
#with the PRESSTO enhancer-promoter association map.
#the universe will be created from the overlap with all regions in the raw differential methylayion table
#the query will be created from the overlap with the DMRs
mart<- useDataset(dataset = "hsapiens_gene_ensembl" , useMart("ensembl"))

if(is.enhancer){
    GR.diff.table = makeGRangesFromDataFrame(cur.diff.table,keep.extra.columns = T)
    GR.dmr.table  = makeGRangesFromDataFrame(cur.dmr.table[[j]],keep.extra.columns = T)
    assoc.enh.prom.univ = subsetByOverlaps(GR.enh.refseqProm.assoc, GR.diff.table)
    assoc.enh.prom.query = subsetByOverlaps(GR.enh.refseqProm.assoc, GR.dmr.table)
    ensembl_query <- assoc.enh.prom.query[,'symbol']$symbol
    ensembl_univ <- assoc.enh.prom.univ[,'symbol']$symbol
    
    xx <- getBM(
        filters= "hgnc_symbol", 
        attributes= c("ensembl_gene_id", "entrezgene", "description"),
        values= ensembl_query,
        mart= mart)
    
    yy <- getBM(
        filters= "hgnc_symbol", 
        attributes= c("ensembl_gene_id", "entrezgene", "description"),
        values = ensembl_univ,
        mart= mart)
    
    head(xx)
    head(yy)
    zz = performGOenrichment.diffMeth.entrez(gids = na.omit(xx$entrezgene), uids = na.omit(yy$entrezgene) ,ontology = "BP",assembly = 'hg19')
    
    summary(zz)
    na.omit(xx)
    writeLines(paste(regNames,'::',methDirection[j],'::',paste(smpNames[1],'->',smpNames[2]),sep = ''))
}

else{

xx <- getBM(
    filters= "ensembl_gene_id", 
    attributes= c("ensembl_gene_id", "entrezgene", "description"),
    values= ensembl_query,
    mart= mart)

yy <- getBM(
    filters= "ensembl_gene_id", 
    attributes= c("ensembl_gene_id", "entrezgene", "description"),
    values = ensembl_univ,
    mart= mart)

zz = performGOenrichment.diffMeth.entrez(gids = na.omit(xx$entrezgene), uids = na.omit(yy$entrezgene) ,ontology = "BP",assembly = 'hg19')

summary(zz)
na.omit(xx)
writeLines(paste(regNames,'::',methDirection[j],'::',paste(smpNames[1],'->',smpNames[2]),sep = ''))
}

####################################################################
###########################################################################
#############################################################################
############################################################################
methDirection <- c("hypo","hyper")
j <- 2
ensembl_query <- "MIER1" #as a test query to lookup entrezID by gene symbol .. or refseqID
ensembl_univ <- na.omit(assoc.astro.enh.prom['symbol'])
foof <- mget(ensembl_univ, org.Hs.egSYMBOL2EG, ifnotfound=NA)
mart<- useDataset(dataset = "hsapiens_gene_ensembl" , useMart("ensembl"))

xx <- getBM(
    filters= "ensembl_gene_id", 
    attributes= c("ensembl_gene_id", "entrezgene", "description",'refseq_mrna'),
    values= ensembl_query,
    mart= mart)

yy <- getBM(
    filters= "ensembl_gene_id", 
    attributes= c("ensembl_gene_id", "entrezgene", "description"),
    values = ensembl_univ,
    mart= mart)

zz = performGOenrichment.diffMeth.entrez(gids = na.omit(xx$entrezgene), uids = na.omit(yy$entrezgene) ,ontology = "BP",assembly = 'hg19')

summary(zz)
writeLines(paste(regNames,'::',methDirection[j],': ',paste(smpNames[1],'->',smpNames[2]),sep = ''))


