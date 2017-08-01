library(RnBeads)
library(biomaRt)
library(org.Hs.eg.db)
library(annotate)
#source(file = './processPRESSTOmapping.R')

#This needs to be called in the environment of the GETDMRsJuly_2016 (or later) script##
#make sure my.final.annotslist objedct is available
#Rnbeads hg19 annotations

analysis.dir = "./RnBeads/analysis"
report.dir = file.path(analysis.dir, "GetDMRs_9comps_July19_2016")
if(!exists(x = 'my.final.annots')){
    my.dmr.list <- readRDS(file = paste(report.dir,'/my.annots.final',sep = ''))
    print("loading annotations")
} else{if(exists('my.final.annots')){my.dmr.list <- my.final.annots}}

{
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



#For cases in which our universe of terms in the enrichment analyses
#will be all of the genes under promoters associating with an enhancer in the PRESSTO data
#select the gene symbols from all of the annotated enhancers in the given set

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
}
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

comparisons <- names(my.dmr.list)
comparison <- comparisons[2]
regions <- names(my.dmr.list[[comparison]]) 
region <- regions[1]
geneset <- my.dmr.list[[comparison]][[region]]
cat(comparison)
cat(region)
minRank <- 100
#Separate the table into tables of hypo and hypermethylated regions
diff.col <- colnames(geneset)[grep('mean.diff',colnames(geneset))]
cur.hypo.table <- geneset[geneset[,diff.col] < 0,]
cur.hyper.table <- geneset[geneset[,diff.col] > 0,]

#Rank the regions according to 'combined rank'
cur.hypo.table <- cur.hypo.table[order(x = cur.hypo.table$combinedRank,decreasing = F,na.last = T),]
cur.hyper.table <- cur.hyper.table[order(x = cur.hyper.table$combinedRank,decreasing = F,na.last = T),]

#Select just the best n ranked regions
cur.hypo.table <- cur.hypo.table[1:minRank,]
cur.hyper.table <- cur.hyper.table[1:minRank,]

#Print summary statistics for probe coverage over region
table(cur.hypo.table$probes.covered)
table(cur.hypo.table$num.sites)
table(cur.hyper.table$probes.covered)
dim(cur.hypo.table)
dim(cur.hyper.table)
#######
assigned.promoter.col <- colnames(geneset)[grep('symbol',colnames(geneset))]
ensembl_univ <- geneset[,assigned.promoter.col]
ensembl_query <- cur.hypo.table[,assigned.promoter.col]
ensembl_univ <- ensembl_univ[complete.cases(ensembl_univ)]
ensembl_query <- ensembl_query[complete.cases(ensembl_query)]


mart<- useDataset(dataset = "hsapiens_gene_ensembl" , useMart("ensembl"))


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
    
    cat(comparison)
    cat(region)
    summary(zz)
    na.omit(xx)
    writeLines(paste(regNames,'::',methDirection[j],'::',paste(smpNames[1],'->',smpNames[2]),sep = ''))




############################################################################
cur.dmr.table = list(hypo = cur.hypo.table, hyper = cur.hyper.table)

methDirection <- c("gained methylation","lost methylation")
j <- 1
ensembl_query <- cur.dmr.table[[j]]$id
ensembl_univ <- cur.diff.table$id


#############################################################################
############################################################################
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


