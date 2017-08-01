



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

HBMVEC_EBM_EGM_gene = na.exclude(HBMVEC_EBM_EGM_gene[order(HBMVEC_EBM_EGM_gene$combinedRank),])
HBMVEC_EBM_GBM8.EVs_gene = na.exclude(HBMVEC_EBM_GBM8.EVs_gene[order(HBMVEC_EBM_GBM8.EVs_gene$combinedRank),])

top = 250
my.range = 10000
EBM_EGM_top = HBMVEC_EBM_EGM_gene[1:top,2:4]
EBM_EGM_top$Start = EBM_EGM_top$Start - my.range
EBM_EGM_top$End = EBM_EGM_top$End + my.range
EBM_EGM_top = EBM_EGM_top[order(EBM_EGM_top$Chromosome,EBM_EGM_top$Start),]
#Complete the sorting in unix with the Version sort -V option
write.table(x = EBM_EGM_top,file = './Top250gene.CmbRnk.ebm_egm.dmr.bed',
            quote = F,row.names = F,sep = '\t')

EBM_EV_top = HBMVEC_EBM_GBM8.EVs_gene[1:top,2:4]
EBM_EV_top$Start = EBM_EV_top$Start - my.range
EBM_EV_top$End = EBM_EV_top$End + my.range
EBM_EV_top = EBM_EV_top[order(EBM_EV_top$Chromosome),]
#Complete the sorting in unix with the Version sort -V option
write.table(x = EBM_EV_top,file = './Top250gene.CmbRnk.ebm_ev.dmr.bed',
            quote = F,row.names = F,sep = '\t')

commonGenes = intersect(EBM_EV_top$symbol,EBM_EGM_top$symbol)
length(commonGenes)
commonGenes
