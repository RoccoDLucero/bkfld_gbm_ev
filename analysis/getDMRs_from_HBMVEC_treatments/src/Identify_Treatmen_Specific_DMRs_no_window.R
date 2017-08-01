



#Import the tables of differentially methylated probes from the
#RnBeads differential methylation analysis.
DMR_output_Dir = "../2016_04_13_getDMRs_from_HBMVEC_treatments/RnBeads/analysis/reports_pairwise/differential_methylation_data/"
siteLevel = T
geneLevel = F
promoterLevel = F
if(siteLevel){
HBMVEC_EBM_GBM8.EVs_site = read.csv(
    paste(DMR_output_Dir,"diffMethTable_site_cmp1.csv",sep = ''),
    stringsAsFactors=FALSE
    )

HBMVEC_EGM_GBM8.EVs_site = read.csv(
    paste(DMR_output_Dir,"diffMethTable_site_cmp2.csv",sep = ''),
    stringsAsFactors=FALSE
    )

HBMVEC_EBM_EGM_site = read.csv(
    paste(DMR_output_Dir,"diffMethTable_site_cmp3.csv",sep = ''),
    stringsAsFactors=FALSE
    )

}

if(geneLevel){
    HBMVEC_EBM_GBM8.EVs_gene = read.csv(
        paste(DMR_output_Dir,"diffMethTable_region_cmp1_genes.csv",sep = ''),
        stringsAsFactors=FALSE
    )
    
    HBMVEC_EGM_GBM8.EVs_gene = read.csv(
        paste(DMR_output_Dir,"diffMethTable_region_cmp2_genes.csv",sep = ''),
        stringsAsFactors=FALSE
    )
    
    HBMVEC_EBM_EGM_gene = read.csv(
        paste(DMR_output_Dir,"diffMethTable_region_cmp3_genes.csv",sep = ''),
        stringsAsFactors=FALSE
    )
    
}

if(promoterLevel){
    HBMVEC_EBM_GBM8.EVs_prom = read.csv(
        paste(DMR_output_Dir,"diffMethTable_region_cmp1_promoters.csv",sep = ''),
        stringsAsFactors=FALSE
    )
    
    HBMVEC_EGM_GBM8.EVs_prom = read.csv(
        paste(DMR_output_Dir,"diffMethTable_region_cmp2_promoters.csv",sep = ''),
        stringsAsFactors=FALSE
    )
    
    HBMVEC_EBM_EGM_prom = read.csv(
        paste(DMR_output_Dir,"diffMethTable_region_cmp3_promoters.csv",sep = ''),
        stringsAsFactors=FALSE
    )
    
}
###############################################################################
#With the data loaded look for differences in the DMRs between 'treatment' and 'control'
#First Identify shared DMRs of notable rank, eg top 500 genes. If the overlap is large,
#Remove those that are in the same direction in both groups relative to reference methylation
#The new list can be used to investigate specific effects of a given treatment 

#This pulls in a function that returns only the desired number of best ranked
#featues from the diffmeth table output of RnBeads
#Complete the sorting in unix with the Version sort -V option to make the
#file compatible with BEDtools
#source('../../Rdat/my.getTopDMRs1.R')
source('../../Rdat/my.getTopDMRs2.5.R')

outName = paste(deparse(substitute(HBMVEC_EBM_GBM8.EVs_site)),"noFlank",sep = '')
Top400HMVECplusEVs = my.getTopRanked(HBMVEC_EBM_GBM8.EVs_site,
                                feature.type = 'site',rankCut = 400,
                                flankSize = 1,
                                fileOut = outName
)

outName = paste(deparse(substitute(HBMVEC_EBM_EGM_site)),"noFlank",sep ='')
Top400HMVECinEGM = my.getTopRanked(HBMVEC_EBM_EGM_site,
                                   feature.type = 'site',rankCut = 400,
                                   flankSize = 1,
                                   fileOut = outName
)

#We can also get the treatment specific differentially methylated CpGs
#By removing those that change states under other treatments

#Set this up later...

#commonGenes = intersect(EBM_EV_top$symbol,EBM_EGM_top$symbol)
#length(commonGenes)
#commonGenes


shared.DMRs.EV_EGM <- read.delim("H:/Dropbox/BRL/proj/bkfld_gbm_ev/analysis/2016_04_13_TF_at_DMR/shared.DMRs.EV_EGM.bed", header=FALSE)

library(data.table)
keys <- c("Chromosome", "Start")
merge(shared.DMRs.EV_EGM, HBMVEC_EBM_EGM_site, by=keys)
tData <- data.table(testData, key=keys)
tBounce <- data.table(testBounce, key=keys)





