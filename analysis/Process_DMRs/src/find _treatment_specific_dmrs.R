cmp1_promoters <- read.csv("H:/Dropbox/BRL/proj/bkfld_gbm_ev/analysis/2016_04_13_getDMRs_from_HBMVEC_treatments/RnBeads/analysis/noGBM8_noNegCtrl_w_astro_endo_enhancers/differential_methylation_data/diffMethTable_region_cmp1_promoters.csv")
cmp2_promoters <- read.csv("H:/Dropbox/BRL/proj/bkfld_gbm_ev/analysis/2016_04_13_getDMRs_from_HBMVEC_treatments/RnBeads/analysis/noGBM8_noNegCtrl_w_astro_endo_enhancers/differential_methylation_data/diffMethTable_region_cmp2_promoters.csv")
cmp3_promoters <- read.csv("H:/Dropbox/BRL/proj/bkfld_gbm_ev/analysis/2016_04_13_getDMRs_from_HBMVEC_treatments/RnBeads/analysis/noGBM8_noNegCtrl_w_astro_endo_enhancers/differential_methylation_data/diffMethTable_region_cmp3_promoters.csv")

cmp1_promoters = cmp1_promoters[(order(cmp1_promoters$combinedRank)),]
cmp2_promoters = cmp2_promoters[(order(cmp2_promoters$combinedRank)),]
cmp3_promoters = cmp3_promoters[(order(cmp3_promoters$combinedRank)),]

cmp1_promoters.100 = cmp1_promoters[1:100,]
cmp2_promoters.100 = cmp2_promoters[1:100,]
cmp3_promoters.100 = cmp3_promoters[1:100,]
colnames(cmp1_promoters.100)


a =  intersect(cmp1_promoters.100$entrezID,cmp2_promoters.100$entrezID)
b = intersect(cmp1_promoters.100$entrezID,cmp3_promoters.100$entrezID) 
c = intersect(cmp3_promoters.100$entrezID,cmp2_promoters.100$entrezID)




cmp1_genes <- read.csv("H:/Dropbox/BRL/proj/bkfld_gbm_ev/analysis/2016_04_13_getDMRs_from_HBMVEC_treatments/RnBeads/analysis/noGBM8_noNegCtrl_w_astro_endo_enhancers/differential_methylation_data/diffMethTable_region_cmp1_genes.csv")
cmp2_genes <- read.csv("H:/Dropbox/BRL/proj/bkfld_gbm_ev/analysis/2016_04_13_getDMRs_from_HBMVEC_treatments/RnBeads/analysis/noGBM8_noNegCtrl_w_astro_endo_enhancers/differential_methylation_data/diffMethTable_region_cmp2_genes.csv")
cmp3_genes <- read.csv("H:/Dropbox/BRL/proj/bkfld_gbm_ev/analysis/2016_04_13_getDMRs_from_HBMVEC_treatments/RnBeads/analysis/noGBM8_noNegCtrl_w_astro_endo_enhancers/differential_methylation_data/diffMethTable_region_cmp3_genes.csv")

cmp1_genes = cmp1_genes[(order(cmp1_genes$combinedRank)),]
cmp2_genes = cmp2_genes[(order(cmp2_genes$combinedRank)),]
cmp3_genes = cmp3_genes[(order(cmp3_genes$combinedRank)),]

cmp1_genes.100 = cmp1_genes[1:100,]
cmp2_genes.100 = cmp2_genes[1:100,]
cmp3_genes.100 = cmp3_genes[1:100,]



aa = intersect(cmp1_genes.100$entrezID,cmp2_genes.100$entrezID)
bb = intersect(cmp1_genes.100$entrezID,cmp3_genes.100$entrezID) 
cc = intersect(cmp3_genes.100$entrezID,cmp2_genes.100$entrezID)

intersect(cmp1_genes.100$symbol,cmp2_genes.100$symbol)
intersect(cmp1_genes.100$symbol,cmp3_genes.100$symbol)
intersect(cmp2_genes.100$symbol,cmp3_genes.100$symbol)



intersect(a,aa)
intersect(b,bb)
intersect(c,cc)

intersect(b,c)
intersect(aa,bb)
intersect(aa,cc)


