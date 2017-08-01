#Rocco Lucero June 29 2016
#Show that the DMRs over regulatory regions in the EGF group occur
#over genes and pathways implicated in the angiogenesis phenotype.
#To do this we will look for the genes listed in the MSIGDB GSEA 'Angiogenesis'
#Gene set


bkfld450k <- readRDS(file = '../../data/raw/Breakefield_HBMVEC_450k/bkfld450k.rda')
AngioGenes <- read.delim('./input/Angiogenesis_geneset.txt')

pathRoot <- '../2016_04_13_getDMRs_from_HBMVEC_treatments/RnBeads/analysis/'
analysisFolders = c('/noGBM8_noNegCtrl_w_astro_endo_enhancers')
datFolder = '/differential_methylation_data/'
    
EGF.DMR.data.src <- paste(pathRoot,analysisFolders[1],datFolder,sep = '')
EGF.DMR.dat <- list.files(EGF.DMR.data.src)

#The relevant comparisons must be determined by looking at the file contents
EGF.DMRFiles <- EGF.DMR.dat[grep(EGF.DMR.dat,pattern = '_region_cmp1_')]
EGF.DMPFiles <- EGF.DMR.dat[grep(EGF.DMR.dat,pattern = 'site_cmp1')]

EGF.Prom.DMR <- read.csv(paste(EGF.DMR.data.src,EGF.DMRFiles[3],sep = ''),header = T,
                         sep = ',',quote = "")


EGF.DMPs <- read.csv(paste(EGF.DMR.data.src,EGF.DMPFiles[1],sep = ''),header = T,
                         sep = ',',quote = "")

EGF.DMPs <- EGF.DMPs[order(EGF.DMPs$combinedRank),][1:1000,]

EGF.DMP.gain <- EGF.DMPs[EGF.DMPs$mean.diff > 0.02,]
EGF.DMP.loss <- EGF.DMPs[EGF.DMPs$mean.diff < -0.02,]

EGF.DMP.loss <- bkfld450k[bkfld450k$TargetID %in% EGF.DMP.loss$cgid,]
EGF.DMP.loss <- EGF.DMP.loss[complete.cases(EGF.DMP.loss$UCSC_REFGENE_NAME),]
out.dmps = (unique(EGF.DMP.loss$UCSC_REFGENE_NAME))
write.table(out.dmps,'./output/egfDMPloss.txt',
            quote = F,col.names = F,row.names = F,sep = '\t')




library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id",
                attributes= c("ensembl_gene_id","description"),
                values = EGF.Prom.DMR$id,
                mart= mart)

EGF.Prom.DMR <- EGF.Prom.DMR[EGF.Prom.DMR$id %in% G_list$ensembl_gene_id,]
colnames(G_list) = c('id', 'description')
colnames(EGF.Prom.DMR)
EGF.Prom.DMR <- merge(EGF.Prom.DMR,G_list)
EGF.Prom.DMR1 <- EGF.Prom.DMR[,c(1,5,9,11:17)]

EGF.Prom.DMR1 = EGF.Prom.DMR1[complete.cases(EGF.Prom.DMR1),]


EGF.Prom.DMR1.gain = EGF.Prom.DMR1[EGF.Prom.DMR1$mean.mean.diff>0,]
EGF.Prom.DMR1.loss = EGF.Prom.DMR1[EGF.Prom.DMR1$mean.mean.diff<0,]

EGF.Prom.DMR1.gain = EGF.Prom.DMR1.gain[order(EGF.Prom.DMR1.gain$combinedRank),]
EGF.Prom.DMR1.loss = EGF.Prom.DMR1.loss[order(EGF.Prom.DMR1.loss$combinedRank),]

head(EGF.Prom.DMR1.gain)
head(EGF.Prom.DMR1.loss)

EGF.Prom.DMR1.gain$description[1:100]
EGF.Prom.DMR1.gain$symbol[1:100]

EGF.Prom.DMR1.gain[EGF.Prom.DMR1.gain$symbol %in% AngioGenes[,1],]
EGF.Prom.DMR1.loss[EGF.Prom.DMR1.loss$symbol %in% AngioGenes[,1],]


################################################################################

EGF.Gene.DMR <- read.csv(paste(EGF.DMR.data.src,EGF.DMRs[3],sep = ''),header = T,
                         sep = ',',quote = "")

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id",
                attributes= c("ensembl_gene_id","description"),
                values = EGF.Gene.DMR$id,
                mart= mart)

EGF.Gene.DMR <- EGF.Gene.DMR[EGF.Gene.DMR$id %in% G_list$ensembl_gene_id,]
colnames(G_list) = c('id', 'description')
EGF.Gene.DMR <- merge(EGF.Gene.DMR,G_list)
EGF.Gene.DMR1 <- EGF.Gene.DMR[,c(1,5,9,11:17)]

EGF.Gene.DMR1 = EGF.Gene.DMR1[complete.cases(EGF.Gene.DMR1),]


EGF.Gene.DMR1.gain = EGF.Gene.DMR1[EGF.Gene.DMR1$mean.mean.diff>0,]
EGF.Gene.DMR1.loss = EGF.Prom.DMR1[EGF.Gene.DMR1$mean.mean.diff<0,]

EGF.Gene.DMR1.gain = EGF.Gene.DMR1.gain[order(EGF.Gene.DMR1.gain$combinedRank),]
EGF.Gene.DMR1.loss = EGF.Gene.DMR1.loss[order(EGF.Gene.DMR1.loss$combinedRank),]

head(EGF.Gene.DMR1.gain)
head(EGF.Gene.DMR1.loss)


a = EGF.Gene.DMR1.gain[EGF.Gene.DMR1.gain$symbol %in% AngioGenes[,1],]
b = EGF.Gene.DMR1.loss[EGF.Gene.DMR1.loss$symbol %in% AngioGenes[,1],]
a$mean.mean.diff
b$mean.mean.diff


