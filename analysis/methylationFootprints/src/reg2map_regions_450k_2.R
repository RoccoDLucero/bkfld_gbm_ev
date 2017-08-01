library(GenomicRanges)

raw.data.dir = "../../data/raw/Breakefield_HBMVEC_450k/"
#Data Import
#Here we have ENCOD ChIP-seq data for E2F1
E2F1.ENCODE <- read.delim("H:/Dropbox/BRL/proj/bkfld_gbm_ev/data/raw/tfchip/E2f1.ENCODE.sorted.bed", header=FALSE)

#Using REG2MAP we get open chromatin enhancers and promoters for HUVEC primary cells
regions_enh_E122 <- read.delim("H:/Dropbox/BRL/proj/bkfld_gbm_ev/data/raw/reg2map/regions_enh_E122.bed", header=FALSE)
regions_prom_E122 <- read.delim("H:/Dropbox/BRL/proj/bkfld_gbm_ev/data/raw/reg2map/regions_prom_E122.bed", header=FALSE)

#Read in the 450K data for the 9 conditions by three replicates
breakefield450K = read.delim(paste(raw.data.dir,"FinalReport2.txt",sep = ''),header = T, stringsAsFactors=F)
colnames(breakefield450K) = gsub("_zappulliMGH1_..", "", colnames(breakefield450K))
colnames(breakefield450K) = gsub("..", ".", colnames(breakefield450K),fixed = T)
methy = cbind(breakefield450K[,c('CHR','MAPINFO','GENOME_BUILD')],breakefield450K[,grep('AVG_Beta',colnames(breakefield450K))])
methy = methy[,c(1:3,4:6,12:14,20:22)]
methy = methy[order(paste(methy$CHR,methy$MAPINFO),as.numeric(methy$MAPINFO),na.last = T),]
methy = methy[complete.cases(methy),]
methy$CHR = paste('chr',methy$CHR,sep = '')

#################################################################################
#Make GRANGES objects that we will need:
#Testing here... Take the EBM samples and the Growth Factor treated samples
#Make them into a granges object so we can quickly test enrichment, etc.
gr450Kprobes = GRanges(seqnames = methy$CHR,
                       ranges = IRanges(start = methy$MAPINFO, end = methy$MAPINFO),
                       mcols = methy[,4:12])

#E2F1 is the first TF we will test
#Make the E2F1 TFBS bed fie into a granges object
grE2F1.ENCODE = GRanges(seqnames = E2F1.ENCODE$V1,
                        ranges = IRanges(start = E2F1.ENCODE$V2, end = E2F1.ENCODE$V3))

#Make REG2MAP regions into GRANGES objects
grHUVEC.enh = GRanges(seqnames = regions_enh_E122$V1,
                           ranges = IRanges(start = regions_enh_E122$V2,end = regions_enh_E122$V3))
grHUVEC.prom = GRanges(seqnames = regions_prom_E122$V1,
                           ranges = IRanges(start = regions_prom_E122$V2,end = regions_prom_E122$V3))

####################################################################################
#Get all of the 450K probes that fall within an endothelial cell
#enhancer region or promoter region
#These are the most improtant probes for footprinting

enhProbes = subsetByOverlaps(gr450Kprobes,grHUVEC.enh)
promProbes = subsetByOverlaps(gr450Kprobes,grHUVEC.prom)

#Get all of the 450K probes that are not in endothelial cell enhancers or promoters
#These may be likened to a negative control for footprinting
non.enhProbes = setdiff(gr450Kprobes,grHUVEC.enh) 
non.enhProbes = subsetByOverlaps(gr450Kprobes,non.enhProbes)

non.promProbes = setdiff(gr450Kprobes,grHUVEC.prom)
non.promProbes = subsetByOverlaps(gr450Kprobes,non.promProbes)

#####################################################################################
######################################################################################
#As a test for the methylation values at enhancers vs all other elements in basal condition:
#plot cumulative distribution function for beta values of probes wihin and
#outside of enhancers for basal media
plot(ecdf(enhProbes@elementMetadata@listData$mcols.EBM1.AVG_Beta),col = 'red')
lines(ecdf(enhProbes@elementMetadata@listData$mcols.EBM2.AVG_Beta),col='blue')
lines(ecdf(enhProbes@elementMetadata@listData$mcols.EBM3.AVG_Beta),col = 'green')
lines(ecdf(non.enhProbes@elementMetadata@listData$mcols.EBM1.AVG_Beta))
lines(ecdf(non.enhProbes@elementMetadata@listData$mcols.EBM2.AVG_Beta))
lines(ecdf(non.enhProbes@elementMetadata@listData$mcols.EBM3.AVG_Beta))
ks.test(enhProbes@elementMetadata@listData$mcols.EBM1.AVG_Beta,
        non.enhProbes@elementMetadata@listData$mcols.EBM1.AVG_Beta)

#As a test for the methylation values at enhancers vs other elements after GF treatment:
#plot cumulative distribution function for beta values of probes wihin and
#outside of enhancers for basal media
plot(ecdf(enhProbes@elementMetadata@listData$mcols.EGM1.AVG_Beta),col = 'red')
lines(ecdf(enhProbes@elementMetadata@listData$mcols.EGM2.AVG_Beta),col='blue')
lines(ecdf(enhProbes@elementMetadata@listData$mcols.EGM3.AVG_Beta),col = 'green')
lines(ecdf(non.enhProbes@elementMetadata@listData$mcols.EGM1.AVG_Beta))
lines(ecdf(non.enhProbes@elementMetadata@listData$mcols.EGM2.AVG_Beta))
lines(ecdf(non.enhProbes@elementMetadata@listData$mcols.EGM3.AVG_Beta))

#Does enhancer methylation change globally from EBM to EGM? --> Not so for batch uncorrected data 
plot(ecdf(enhProbes@elementMetadata@listData$mcols.EBM1.AVG_Beta),col = 'red')
lines(ecdf(enhProbes@elementMetadata@listData$mcols.EBM2.AVG_Beta),col='red')
lines(ecdf(enhProbes@elementMetadata@listData$mcols.EBM3.AVG_Beta),col = 'red')
lines(ecdf(enhProbes@elementMetadata@listData$mcols.EGM1.AVG_Beta),col = 'purple')
lines(ecdf(enhProbes@elementMetadata@listData$mcols.EGM2.AVG_Beta),col='purple')
lines(ecdf(enhProbes@elementMetadata@listData$mcols.EGM3.AVG_Beta),col = 'purple')
lines(ecdf(enhProbes@elementMetadata@listData$mcols.EVpellet1.AVG_Beta),col = 'green')
lines(ecdf(enhProbes@elementMetadata@listData$mcols.EVpellet2.AVG_Beta),col='green')
lines(ecdf(enhProbes@elementMetadata@listData$mcols.EVpellet3.AVG_Beta),col = 'green')

#Does promoter methylation change globally from EBM to EGM? --> Not so for batch uncorrected data 
plot(ecdf(promProbes@elementMetadata@listData$mcols.EBM1.AVG_Beta),col = 'black')
lines(ecdf(promProbes@elementMetadata@listData$mcols.EBM2.AVG_Beta),col='black')
lines(ecdf(promProbes@elementMetadata@listData$mcols.EBM3.AVG_Beta),col = 'black')
lines(ecdf(promProbes@elementMetadata@listData$mcols.EGM1.AVG_Beta),col = 'purple')
lines(ecdf(promProbes@elementMetadata@listData$mcols.EGM2.AVG_Beta),col='purple')
lines(ecdf(promProbes@elementMetadata@listData$mcols.EGM3.AVG_Beta),col = 'purple')
lines(ecdf(promProbes@elementMetadata@listData$mcols.EVpellet1.AVG_Beta),col = 'green')
lines(ecdf(promProbes@elementMetadata@listData$mcols.EVpellet2.AVG_Beta),col='green')
lines(ecdf(promProbes@elementMetadata@listData$mcols.EVpellet3.AVG_Beta),col = 'green')
#####################################################################################
#####################################################################################
##Incorporate the ENCODE ChIP-seq data for E2F1...

#Get all the 450K probes that overalp HUVEC enhancers
#Try to get stats on the number of probes in each enhancer to see if this is representative
# or whether a few enhancers are biasing the result
grHUVEC.enh.probes = subsetByOverlaps(gr450Kprobes,grHUVEC.enh)

#All probes overlapping E2F1 peaks within endothelial cell active enhancers
grE2F1.HUVEC.enh.probes = subsetByOverlaps(grHUVEC.enh.probes,grE2F1.ENCODE)

#All probes not overlapping E2F1 peaks in endothelial cell active enhancers
grNonE2F1.HUVEC.enh.probes =  setdiff(grHUVEC.enh.probes,grE2F1.HUVEC.enh.probes)
grNonE2F1.HUVEC.enh.probes = subsetByOverlaps(grHUVEC.enh.probes,grNonE2F1.HUVEC.enh.probes)

#Comparing E2F1 sites likely Active to non-active enhancers in endothelial cells we see that the
#activ
summary(grHUVEC.enh.probes)
summary(grE2F1.HUVEC.enh.probes)
summary(grNonE2F1.HUVEC.enh.probes)
head(grHUVEC.enh.probes)
head(grE2F1.HUVEC.enh.probes)
head(grNonE2F1.HUVEC.enh.probes)


plot(ecdf(grNonE2F1.HUVEC.enh.probes@elementMetadata@listData$mcols.EBM2.AVG_Beta),col = 'red')
#lines(ecdf(grNonE2F1.HUVEC.enh.probes@elementMetadata@listData$mcols.EGM3.AVG_Beta),col = 'orange')
lines(ecdf(grNonE2F1.HUVEC.enh.probes@elementMetadata@listData$mcols.EVpellet2.AVG_Beta),col = 'black')
lines(ecdf(grE2F1.HUVEC.enh.probes@elementMetadata@listData$mcols.EBM2.AVG_Beta),col = 'blue',cex = .2)
#lines(ecdf(grE2F1.HUVEC.enh.probes@elementMetadata@listData$mcols.EGM3.AVG_Beta),col = 'pink',cex = .2)
lines(ecdf(grE2F1.HUVEC.enh.probes@elementMetadata@listData$mcols.EVpellet2.AVG_Beta),col = 'purple',cex = .2)


ks.test(endoEnh.E2F1.probes@elementMetadata@listData$mcols.EBM1.AVG_Beta,
        non.endoEnh.E2F1.probes@elementMetadata@listData$mcols.EBM1.AVG_Beta)
median(endoEnh.E2F1.probes@elementMetadata@listData$mcols.EBM3.AVG_Beta)
median(non.endoEnh.E2F1.probes@elementMetadata@listData$mcols.EBM3.AVG_Beta)

##################################################################################################
##################################################################################################
##Incorporate the ENCODE ChIP-seq data for E2F1...

#Get all the 450K probes that overalp HUVEC promoters
#Try to get stats on the number of probes in each promoter to see if this is representative
# or whether a few enhancers are biasing the result
grHUVEC.prom.probes = subsetByOverlaps(gr450Kprobes,grHUVEC.prom)

#All probes overlapping E2F1 peaks within endothelial cell active enhancers
grE2F1.HUVEC.prom.probes = subsetByOverlaps(grHUVEC.prom.probes,grE2F1.ENCODE)

#All probes not overlapping E2F1 peaks in endothelial cell active enhancers
grNonE2F1.HUVEC.prom.probes =  setdiff(grHUVEC.prom.probes,grE2F1.HUVEC.prom.probes)
grNonE2F1.HUVEC.prom.probes = subsetByOverlaps(grHUVEC.prom.probes,grNonE2F1.HUVEC.prom.probes)

#Comparing E2F1 sites likely Active to non-active enhancers in endothelial cells we see that the
#activ
summary(grHUVEC.prom.probes)
summary(grE2F1.HUVEC.prom.probes)
summary(grNonE2F1.HUVEC.prom.probes)
head(grHUVEC.enh.probes)
head(grE2F1.HUVEC.prom.probes)
head(grNonE2F1.HUVEC.prom.probes)

plot(ecdf(grNonE2F1.HUVEC.prom.probes@elementMetadata@listData$mcols.EBM3.AVG_Beta),col = 'red')
lines(ecdf(grE2F1.HUVEC.prom.probes@elementMetadata@listData$mcols.EBM3.AVG_Beta),col = 'blue')
lines(ecdf(grHUVEC.prom.probes@elementMetadata@listData$mcols.EBM3.AVG_Beta),col = 'black')
lines(ecdf(grE2F1.HUVEC.prom.probes@elementMetadata@listData$mcols.EGM3.AVG_Beta),col = 'pink')
lines(ecdf(grNonE2F1.HUVEC.prom.probes@elementMetadata@listData$mcols.EGM3.AVG_Beta),col = 'orange')


plot(ecdf(grE2F1.HUVEC.prom.probes@elementMetadata@listData$mcols.EBM2.AVG_Beta),col = 'blue')
lines(ecdf(grE2F1.HUVEC.prom.probes@elementMetadata@listData$mcols.EGM2.AVG_Beta),col = 'pink')
ks.test(grE2F1.HUVEC.prom.probes@elementMetadata@listData$mcols.EBM3.AVG_Beta,
        grE2F1.HUVEC.prom.probes@elementMetadata@listData$mcols.EGM3.AVG_Beta,alternative = 'tw')

###################################################################################
#Incorporate the ENCODE ChIP-seq data for E2F1...
#All probes overlapping E2F1 peaks within endothelial cell active enhancers
endoEnh.E2F1.X = subsetByOverlaps(enhProbes,grE2F1.ENCODE)
#All probes not overlapping E2F1 peaks in endothelial cell active enhancers
non.endoEnh.E2F1.X = setdiff(enhProbes,grE2F1.ENCODE)
non.endoEnh.E2F1.X = subsetByOverlaps(gr450Kprobes,non.endoEnh.E2F1.X)

plot(ecdf(endoEnh.E2F1.X@elementMetadata@listData$mcols.EBM1.AVG_Beta))
lines(ecdf(non.endoEnh.E2F1.X@elementMetadata@listData$mcols.EGM1.AVG_Beta))


####################################################################################
