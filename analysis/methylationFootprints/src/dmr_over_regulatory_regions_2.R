#Rocco Lucero June 29 2016
#Show that the DMRs over regulatory regions in the EGF group occur
#over genes and pathways implicated in the angiogenesis phenotype.
#To do this we will look for the genes listed in the MSIGDB GSEA 'Angiogenesis'
#Gene set

source('./My.get_DMR_genes_in_genesets.R')
AngioGenes <- read.delim('./input/Angiogenesis_geneset.txt')

DMR.data.src <- '../2016_04_13_getDMRs_from_HBMVEC_treatments/RnBeads/analysis/noGBM8_noNegCtrl_w_astro_endo_enhancers/differential_methylation_data/'
DMRs <- list.files(DMR.data.src)

EGF.DMRs <- DMRs[grep(DMRs,pattern = '_region_cmp1_')]    
EGF.Prom.DMR <- read.csv(paste(DMR.data.src,EGF.DMRs[4],sep = ''),header = T,
                        sep = ',',quote = "")

EV.DMRs <- DMRs[grep(DMRs,pattern = '_region_cmp2_')]    
EV.Prom.DMR <- read.csv(paste(DMR.data.src,EV.DMRs[4],sep = ''),header = T,
                        sep = ',',quote = "")

SUP.DMRs <- DMRs[grep(DMRs,pattern = '_region_cmp3_')]    
SUP.Prom.DMR <- read.csv(paste(DMR.data.src,SUP.DMRs[4],sep = ''),header = T,
                         sep = ',',quote = "")

a <- myDMRGenesFromGeneset(DMRs = EGF.Prom.DMR,geneSet = AngioGenes,maxRank = 1000)
length(a[[3]]$description)
length(a[[4]]$description)
head(a[[1]]$symbol)


b <- myDMRGenesFromGeneset(DMRs = EV.Prom.DMR,geneSet = AngioGenes,maxRank = 1000)
length(b[[3]]$description)
length(b[[4]]$description)
head(b[[1]]$symbol)


c <- myDMRGenesFromGeneset(DMRs = SUP.Prom.DMR,geneSet = AngioGenes,maxRank = 1000)
length(c[[3]]$description)
length(c[[4]]$description)
head(c[[1]]$symbol)
