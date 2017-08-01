#As a quick check for possible batch effects we will run unsupervised
# Heirarchical clustering on the methylation beta values for all of the 24
# 450k samples. We expect to see that the samples cluster by cell type,
# and then by treatment.

#The script randomly selects twenty five thousand probe beta values and shuffles
#the Column order of the beta value matrix before unsupervised heirarchical 
#clustering

##RESULT# We saw that replicate number may have a stronger effect than treatment
##Indicating a possible batch effect.

require(gplots)

#bkfld450k = read.delim(
#    file = '../../data/raw/Breakefield_HBMVEC_450k/FinalReport2.txt',
#    header = T, sep = '\t'
#    )
#saveRDS(bkfld450k, './bkfld450k.rda')
rm(bkfld450k)
bkfld450k = readRDS('./bkfld450k.rda')
#file.link(from = '../../data/raw/Breakefield_HBMVEC_450k/FinalReport2.txt',
#         to = './input/bkfld450k.txt'
#          )

betas = bkfld450k[,grep('AVG_Beta',colnames(bkfld450k))]
colnames(betas) <- gsub(pattern = '_zappulliMGH1_[ABCDEF][1234].AVG_Beta',
                       replacement = '',
                       x = colnames(betas)
                       )
detection.p = bkfld450k[,grep('Detection.Pval',colnames(bkfld450k))]
metas = bkfld450k[,c(1:3,221:254)]

betas = betas[complete.cases(betas),]
detection.p = detection.p[which(complete.cases(betas)),]
metas = metas[which(complete.cases(betas)),]
head(betas)
set.seed(20160520) #date of first analysis
pdf(file = './output/2d_cluster_25000probesx5.pdf')
for(i in 1:5){
    heatmap.2(x = as.matrix(betas[sample(rownames(betas),25000),sample(1:ncol(betas),ncol(betas))]),dendrogram = 'column',
          trace = 'none',labRow = F,margins = c(8,5)
    )
}

dev.off()
gc()

