#In this script I will use the HoneyBadger2 database at
#http://www.broadinstitute.org/~meuleman/reg2map/HoneyBadger2_release/
#to generate the sets of regions of open chromatin and closed chromatin
#for mesodermal and endothelial celltypes.
# will look specifically at the enhancers and promoters ignoring dyadic at firstreg


library(curl)
##Make the download happen automatically:
downloadDir = '../../data/raw/reg2map/'
#First get the sample info to identify relevant cell-type epigenomes:
indexLink ='http://www.broadinstitute.org/~meuleman/reg2map/sample_info/sample_info.RData'
download.file(url = indexLink,destfile = downloadDir)
load('../../data/raw/reg2map/sample_info.RData')
grep('Endoth',sample_info$name)  #E120
grep('Mesoderm',sample_info$name) #E13
#I should use the epigenomes 13, and 120 which are mesoderm and endothelial in origin
head(sample_info[c(13,120),])
#EID       mnemonic                                           name   color group_name position
#13  E013 ESDR.CD56.MESO     hESC Derived CD56+ Mesoderm Cultured Cells #4178AE   ES-deriv       18
#120 E122      VAS.HUVEC HUVEC Umbilical Vein Endothelial Primary Cells #000000 ENCODE2012      120

#Now 
DownloadFromHost = 'http://www.broadinstitute.org/~meuleman/reg2map/HoneyBadger2_release/DNase/p10/'
#download.file()
#'enh'
#'prom'
#'BED_files_per_sample/'
#paste('regions_enh_E', sample ,'.bed',sep = '')

regions_enh_E120 <- read.delim("H:/Dropbox/BRL/proj/bkfld_gbm_ev/data/raw/reg2map/regions_enh_E120.bed", header=FALSE)
head(regions_enh_E120)
mean(regions_enh_E120$V3 - regions_enh_E120$V2) #~440bp ranges for open enhancers in endothelial ~huvec

