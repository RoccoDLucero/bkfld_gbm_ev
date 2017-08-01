


library(GenomicRanges)

##################################################################################
#Data Import for CpG Methylation of Samples and Sample selection
##################################################################################
#Read in the 450K data for the 9 conditions by three replicates
#This code is specific to the Breakefield Data set
my.450K.file <- "../../data/raw/Breakefield_HBMVEC_450k/FinalReport2.txt"
my.450K.file <- read.delim(my.450K.file, header = T, stringsAsFactors=F)
colnames(my.450K.file) <- gsub("_zappulliMGH1_..", "", colnames(my.450K.file))
colnames(my.450K.file) <- gsub("..", ".", colnames(my.450K.file),fixed = T)
#Subset and modify the data to get genomic coordiantes and Beta Values that can
#be used to build GRANGES objects
methy <- cbind(my.450K.file[,c('CHR','MAPINFO','GENOME_BUILD')],
               my.450K.file[,grep('AVG_Beta',colnames(my.450K.file))])

#In this test we take the EBM, EGM and EV-treated Samples
my.sample.columns = c(4:6,12:14,20:22)
my.sample.names = colnames(methy[,my.sample.columns])
methy <- methy[,c(1:3,my.sample.columns)]
methy <- methy[order(paste(methy$CHR,methy$MAPINFO),as.numeric(methy$MAPINFO),na.last = T),]
methy <- methy[complete.cases(methy),]
methy$CHR <- paste('chr',methy$CHR,sep = '')

#Make GRANGES objects that we will need:
gr450Kprobes <- GRanges(seqnames = methy$CHR,
                       ranges = IRanges(start = methy$MAPINFO, end = methy$MAPINFO),
                       mcols = methy[,my.sample.names])

##################################################################################
#Testing Here: We will only look at REG2MAP regulatory regions of HUVECS for now
#Get open chromatin enhancers and promoters for HUVEC primary cells

##Later this should go directly to the site to download files specified in
##a "REG2MAP.cell-types.csv" file
regions_enh_E122 <- read.delim("../../data/raw/reg2map/regions_enh_E122.bed", header=FALSE)
regions_prom_E122 <- read.delim("../../data/raw/reg2map/regions_prom_E122.bed", header=FALSE)

#Make REG2MAP regions into GRANGES objects 
grEnh = GRanges(seqnames = regions_enh_E122$V1,
                      ranges = IRanges(start = regions_enh_E122$V2,end = regions_enh_E122$V3))
grProm = GRanges(seqnames = regions_prom_E122$V1,
                       ranges = IRanges(start = regions_prom_E122$V2,end = regions_prom_E122$V3))

##All of the above in this section needs to be made more generic

####################################################################################
#Get all of the 450K probes that fall within chosen Regualtory Regions
####################################################################################
#These are the most improtant probes for footprinting
enhProbes = subsetByOverlaps(gr450Kprobes,grEnh)
promProbes = subsetByOverlaps(gr450Kprobes,grProm)

#Get statistics about the number of elements with at least one probe and
#average numbers of probes per element
table(countOverlaps(grEnh,enhProbes)) # indicates that 92.5% of HUVEC enhancers ARE NOT on 450K
    #~2500 enhancers have at least 5 probes
table(countOverlaps(grProm,promProbes)) #indicates that in HUVEC 36% of promoters are not covered
summary(grProm)
199162/215301
8084/22298
#Get all of the 450K probes that are not in endothelial cell enhancers or promoters
#These may be likened to a negative control for footprinting
non.enhProbes = setdiff(gr450Kprobes,grEnh) 
non.enhProbes = subsetByOverlaps(gr450Kprobes,non.enhProbes)

non.promProbes = setdiff(gr450Kprobes,grProm)
non.promProbes = subsetByOverlaps(gr450Kprobes,non.promProbes)
#############################################################################
################################################################################
#Incorporate the DMRs from RnBEADS output
################################################################################
my.input.files <- list.files(path = './input/')
current.dmr.file <- paste('./input/',my.input.files[1],sep = '')
my.current.DMRs <- read.csv(file = current.dmr.file, header = T)
DMRcols <- colnames(my.current.DMRs)[c(2,8:10,18,19)]
grDMRprobes <- GRanges(seqnames = my.current.DMRs$Chromosome,
                       ranges = IRanges(start = my.current.DMRs$Start, end = my.current.DMRs$Start),
                       mcols = my.current.DMRs[,DMRcols],strand = my.current.DMRs$Strand)
DMR.REG1 <- subsetByOverlaps(grDMRprobes,enhProbes)

current.dmr.file <- paste('./input/',my.input.files[2],sep = '')
my.current.DMRs <- read.csv(file = current.dmr.file, header = T)
DMRcols <- colnames(my.current.DMRs)[c(2,8:10,18,19)]
grDMRprobes <- GRanges(seqnames = my.current.DMRs$Chromosome,
                       ranges = IRanges(start = my.current.DMRs$Start, end = my.current.DMRs$Start),
                       mcols = my.current.DMRs[,DMRcols],strand = my.current.DMRs$Strand)
DMR.REG2 <- subsetByOverlaps(grDMRprobes,enhProbes)

plot(ecdf(DMR.REG1@elementMetadata@listData$mcols.mean.diff),col = 'red1', cex = .2)
lines(ecdf(DMR.REG2@elementMetadata@listData$mcols.mean.diff),col = 'blue1', cex = .2)


#############################################################################
##Incorporate the ENCODE ChIP-seq data 
#################################################################################
#Here we provide a folder of ChIP-seq peak BED files
my.filenames <- list.files(path = '../../data/raw/tfchip/')
TFnames <- gsub(pattern = ".ENCODE.sorted.bed",replacement ='',x = my.filenames)

#we need to loop over this list importing one factor at a time
#running the analysis and generating the output for each factor
for(i in 1:101){
my.current.file = my.filenames[i]
my.current.TF <- read.delim(paste("../../data/raw/tfchip/",
                                  my.current.file,sep = ''),
                            header=FALSE)

#E2F1 is the first TF we will test
#Make the E2F1 TFBS bed fie into a granges object
grTF.ChIPseq = GRanges(seqnames = my.current.TF$V1,
                        ranges = IRanges(start = my.current.TF$V2, end =my.current.TF$V3))

REG.probes = enhProbes
#Get probes where ChIP-seq peaks overlap Regualtory Regions
TF.REG = subsetByOverlaps(REG.probes,grTF.ChIPseq)
#All probes not overlapping ChIP-seq peaks in Regulatory region
non.TF.REG = setdiff(REG.probes,grTF.ChIPseq)
non.TF.REG = subsetByOverlaps(REG.probes,non.TF.REG)

#GET all probes overlapping ChIP-seq peak but not in regulatory region
#Needs to be coded later
#non.REG.TF = setdiff(grTF.ChIPseq,REG.probes)
#non.REG.TF = subsetByOverlaps(non.REG.TF,REG.probes)

summary(TF.REG)
summary(non.TF.REG)
#summary(non.REG.TF)

#Do methylation values differ where TF binds?
if(F){
plot(ecdf(TF.REG@elementMetadata@listData$mcols.EGM1.AVG_Beta),col = 'red1', cex = .2)
lines(ecdf(TF.REG@elementMetadata@listData$mcols.EGM2.AVG_Beta),col = 'red3', cex = .2)
lines(ecdf(TF.REG@elementMetadata@listData$mcols.EGM3.AVG_Beta),col = 'red4', cex = .2)
lines(ecdf(non.TF.REG@elementMetadata@listData$mcols.EGM1.AVG_Beta), col = 'blue1', cex =.2)
lines(ecdf(non.TF.REG@elementMetadata@listData$mcols.EGM2.AVG_Beta), col = 'blue3', cex =.2)
lines(ecdf(non.TF.REG@elementMetadata@listData$mcols.EGM3.AVG_Beta), col = 'blue4', cex =.2)

plot(ecdf(TF.REG@elementMetadata@listData$mcols.EVpellet1.AVG_Beta),col = 'red1', cex = .2)
lines(ecdf(TF.REG@elementMetadata@listData$mcols.EVpellet2.AVG_Beta),col = 'red3', cex = .2)
lines(ecdf(TF.REG@elementMetadata@listData$mcols.EVpellet3.AVG_Beta),col = 'red4', cex = .2)
lines(ecdf(non.TF.REG@elementMetadata@listData$mcols.EVpellet1.AVG_Beta), col = 'blue1', cex =.2)
lines(ecdf(non.TF.REG@elementMetadata@listData$mcols.EVpellet2.AVG_Beta), col = 'blue3', cex =.2)
lines(ecdf(non.TF.REG@elementMetadata@listData$mcols.EVpellet3.AVG_Beta), col = 'blue4', cex =.2)

plot(ecdf(REG.probes@elementMetadata@listData$mcols.EBM1.AVG_Beta),col = 'red1', cex = .2)
lines(ecdf(REG.probes@elementMetadata@listData$mcols.EBM2.AVG_Beta),col = 'red3', cex = .2)
lines(ecdf(REG.probes@elementMetadata@listData$mcols.EBM3.AVG_Beta),col = 'red4', cex = .2)
lines(ecdf(REG.probes@elementMetadata@listData$mcols.EGM1.AVG_Beta), col = 'blue1', cex =.2)
lines(ecdf(REG.probes@elementMetadata@listData$mcols.EGM2.AVG_Beta), col = 'blue3', cex =.2)
lines(ecdf(REG.probes@elementMetadata@listData$mcols.EGM3.AVG_Beta), col = 'blue4', cex =.2)
lines(ecdf(REG.probes@elementMetadata@listData$mcols.EVpellet1.AVG_Beta),col = 'green1', cex = .2)
lines(ecdf(REG.probes@elementMetadata@listData$mcols.EVpellet2.AVG_Beta),col = 'green3', cex = .2)
lines(ecdf(REG.probes@elementMetadata@listData$mcols.EVpellet3.AVG_Beta),col = 'green4', cex = .2)
}
#YES
if(length(TF.REG)!=0){
if(F){
plot(ecdf(TF.REG@elementMetadata@listData$mcols.EBM1.AVG_Beta),col = 'red1', cex = .2)
lines(ecdf(TF.REG@elementMetadata@listData$mcols.EBM2.AVG_Beta),col = 'red3', cex = .2)
lines(ecdf(TF.REG@elementMetadata@listData$mcols.EBM3.AVG_Beta),col = 'red4', cex = .2)
lines(ecdf(TF.REG@elementMetadata@listData$mcols.EVpellet1.AVG_Beta),col = 'blue1', cex = .2)
lines(ecdf(TF.REG@elementMetadata@listData$mcols.EVpellet2.AVG_Beta),col = 'blue3', cex = .2)
lines(ecdf(TF.REG@elementMetadata@listData$mcols.EVpellet3.AVG_Beta),col = 'blue4', cex = .2)
lines(ecdf(TF.REG@elementMetadata@listData$mcols.EGM1.AVG_Beta),col = 'green1', cex = .2)
lines(ecdf(TF.REG@elementMetadata@listData$mcols.EGM2.AVG_Beta),col = 'green3', cex = .2)
lines(ecdf(TF.REG@elementMetadata@listData$mcols.EGM3.AVG_Beta),col = 'green4', cex = .2)
}
DMR.TF.REG1 <- subsetByOverlaps(DMR.REG1,grTF.ChIPseq)
DMR.TF.REG2 <- subsetByOverlaps(DMR.REG2,grTF.ChIPseq)

if(length(DMR.TF.REG1)!=0 & length(DMR.TF.REG2)!=0){
plot(ecdf(DMR.TF.REG1@elementMetadata@listData$mcols.mean.diff),col = 'red1', cex = .2)
lines(ecdf(DMR.TF.REG2@elementMetadata@listData$mcols.mean.diff),col = 'blue1', cex = .2)
 }


}
####################################################################################
#6/22/2016 There is some kind of systematic bias in methylation values, especially for the EBM condition
# That is contributing a lot of variance within the enhancer and promoter regions and
#probably across the entire probe set.





