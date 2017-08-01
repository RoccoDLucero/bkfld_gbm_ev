

mirTARBase <- read.delim('./input/hsa_MTI.txt',header = T,stringsAsFactors = F)
mirTAR.Strong <- mirTARBase[mirTARBase$Support.Type == unique(mirTARBase$Support.Type)[1],]

#smallRNAFolders <- list.dirs("./input/krichevsky_EV_smallRNA/")
s1 <- read.table(file = "./input/krichevsky_EV_smallRNA/sample_GBM1_exosome_clean_fq/readCounts_miRNAmature_sense.txt",header = T)
s2 <- read.table(file = "./input/krichevsky_EV_smallRNA/sample_GBM2_exosome_clean_fq/readCounts_miRNAmature_sense.txt",header = T)
s3 <- read.table(file = "./input/krichevsky_EV_smallRNA/sample_GBM3_exosome_clean_fq/readCounts_miRNAmature_sense.txt",header = T)
s4 <- read.table(file = "./input/krichevsky_EV_smallRNA/sample_GBM4_exosome_clean_fq/readCounts_miRNAmature_sense.txt",header = T)

depth <- 25
s1 <- s1[order(s1$uniqueReadCount,decreasing = T),][s1$uniqueReadCount >= quantile(s1$uniqueReadCount)[4],]
s2 <- s2[order(s2$uniqueReadCount,decreasing = T),][s2$uniqueReadCount >= quantile(s2$uniqueReadCount)[4],]
s3 <- s3[order(s3$uniqueReadCount,decreasing = T),][s3$uniqueReadCount >= quantile(s3$uniqueReadCount)[4],]
s4 <- s4[order(s4$uniqueReadCount,decreasing = T),][s4$uniqueReadCount >= quantile(s4$uniqueReadCount)[4],]

a <- intersect(s1$ReferenceID,s2$ReferenceID)
b <- intersect(s3$ReferenceID,s4$ReferenceID)
c1 <- intersect(a,b)



##########################################################
pth_ct1 <- "./input/KJENS_healthy_Controls/sample_SAMPLE_0438_CONTROL_CSF_fastq/readCounts_miRNAmature_sense.txt"
pth_ct2 <- "./input/KJENS_healthy_Controls/sample_SAMPLE_0758_CONTROL_CSF_fastq/readCounts_miRNAmature_sense.txt"
pth_ct3 <- "./input/KJENS_healthy_Controls/sample_SAMPLE_0935_CONTROL_CSF_fastq/readCounts_miRNAmature_sense.txt"
pth_ct4 <- "./input/KJENS_healthy_Controls/sample_SAMPLE_9635_CONTROL_CSF_fastq/readCounts_miRNAmature_sense.txt"

ct1 <- read.table(file = pth_ct1,header = T)
ct2 <- read.table(file = pth_ct2,header = T)
ct3 <- read.table(file = pth_ct3,header = T)
ct4 <- read.table(file = pth_ct4,header = T)

depth <- 20
ct1 <- ct1[order(ct1$uniqueReadCount,decreasing = T),][ct1$uniqueReadCount >= quantile(ct1$uniqueReadCount)[4],]
ct2 <- ct2[order(ct2$uniqueReadCount,decreasing = T),][ct2$uniqueReadCount >= quantile(ct2$uniqueReadCount)[4],]
ct3 <- ct3[order(ct3$uniqueReadCount,decreasing = T),][ct3$uniqueReadCount >= quantile(ct3$uniqueReadCount)[4],]
ct4 <- ct4[order(ct4$uniqueReadCount,decreasing = T),][ct4$uniqueReadCount >= quantile(ct4$uniqueReadCount)[4],]

a <- intersect(ct1$ReferenceID,ct2$ReferenceID)
b <- intersect(ct3$ReferenceID,ct4$ReferenceID)
c2 <- intersect(a,b)

#####################################################
gbm <- matrix(unlist(strsplit(c1,':|\\|')),ncol = 5,byrow = T)[,1]
ctr <- matrix(unlist(strsplit(c2,':|\\|')),ncol = 5,byrow = T)[,1]
#gbm <- matrix(unlist(strsplit(gbm.unq.mirna,':')),ncol = 5,byrow = T)[,1]
#ctr <- matrix(unlist(strsplit(ctl.unq.mirna,':')),ncol = 5,byrow = T)[,1]


#take only unique top mirna
gbm <- setdiff(gbm,ctr)
ctr <- setdiff(ctr,gbm)

gbmTargs.df <- (mirTAR.Strong[mirTAR.Strong$miRNA %in% gbm,])
gbmGenes <- (unique(gbmTargs.df$Target.Gene))

ctTargs.df <- (mirTAR.Strong[mirTAR.Strong$miRNA %in% ctr,])
ctGenes <- (unique(ctTargs.df$Target.Gene))

length(intersect(gbmGenes,ctGenes))
length(gbmGenes)
length(ctGenes)
###
write.table(gbmGenes,"./output/gbm.TOPgenes.txt",sep = "\t",quote = F,row.names = F,col.names = F)
write.table(ctGenes,"./output/controlTOPgenes.txt",sep = "\t",quote = F,row.names = F,col.names = F)

write.table(gbm,"./output/gbm.unq.mirna.txt",sep = "\t",quote = F,row.names = F,col.names = F)
write.table(ctr,"./output/control.unq.mirna.txt",sep = "\t",quote = F,row.names = F,col.names = F)

write.table(intersect(gbmGenes,ctGenes),'./output/gbm.crt_int.txt',sep = "\t",quote = F,row.names = F,col.names = F)
write.table(setdiff(gbmGenes,ctGenes),'./output/gbm.jnk.txt',sep = "\t",quote = F,row.names = F,col.names = F)
write.table(setdiff(ctGenes,gbmGenes),'./output/ctr.jnk.txt',sep = "\t",quote = F,row.names = F,col.names = F)



####################################################################################################
q <- read.table(file = "./input/krichevsky_EV_smallRNA/sample_GBM1_exosome_clean_fq/readCounts_gencode_sense_geneLevel.txt",header = T)
qn <- read.table(file = "./input/KJENS_healthy_Controls/sample_SAMPLE_0438_CONTROL_CSF_fastq/readCounts_gencode_sense_geneLevel.txt",header = T)

my.tx <- 'mgmt' # "^nf1"
q[grep(my.tx,q$ReferenceID,ignore.case = T),]
qn[grep(my.tx,qn$ReferenceID,ignore.case = T),]

quantile(q$uniqueReadCount)
quantile(qn$uniqueReadCount)


Endothelial_enriched_Stanford <- read.delim("F:/Dropbox/BRL/proj/bkfld_gbm_ev/analysis/miRNA_targets_GBM8_EVs/input/Endothelial_enriched_Stanford.txt", 
                                          sep = "\t",header = F,stringsAsFactors = F)
Endothelial_enriched_Stanford <- as.data.frame(matrix(Endothelial_enriched_Stanford$V1,ncol = 4,byrow = T))
Endothelial_enriched_Stanford$V1 <- toupper(Endothelial_enriched_Stanford$V1)
intersect(gbmGenes,Endothelial_enriched_Stanford$V1)
