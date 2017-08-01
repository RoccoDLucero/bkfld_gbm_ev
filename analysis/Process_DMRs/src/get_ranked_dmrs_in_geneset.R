#Rocco Lucero June 29 2016
#Function to take as input:
#   RnBeads differential methylation output file
#   A geneSet
#Function returns#
#   For Genes with a gene symbol in biomart  
#   Regions gaining methylation in combined rank order and associated gene
#   Regions losing methylation in combined rank order and associated gene
#          
#
#
#
################################################################################

myDMRGenesFromGeneset = function(DMRs, geneSet, maxRank){
    
    library('biomaRt')
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    G_list <- getBM(filters= "ensembl_gene_id",
                    attributes= c("ensembl_gene_id","description"),
                    values = DMRs$id,
                    mart= mart)
    
    DMRs <- DMRs[DMRs$id %in% G_list$ensembl_gene_id,]
    colnames(G_list) = c('id', 'description')
    
    DMRs  <- merge(DMRs,G_list)
    #DMRs <- DMRs[,c(1,5,9,11:17)]
    
    DMRs <- DMRs[complete.cases(DMRs),]
    
    
    DMRs.gain <- DMRs[DMRs$mean.mean.diff>0,]
    DMRs.loss <- DMRs[DMRs$mean.mean.diff<0,]
    
    DMRs.gain <- DMRs.gain[order(DMRs.gain$combinedRank),]
    DMRs.loss <- DMRs.loss[order(DMRs.loss$combinedRank),]

    DMRs.gain <- DMRs.gain[1:maxRank,]
    DMRs.loss <- DMRs.loss[1:maxRank,]
    
    DMRs.gain.in.geneset <- DMRs.gain[DMRs.gain$symbol %in% geneSet[,1],]
    DMRs.loss.in.geneset <- DMRs.loss[DMRs.loss$symbol %in% geneSet[,1],]
    
    list(DMRs.gain,DMRs.loss,
         DMRs.gain.in.geneset,DMRs.loss.in.geneset)
    
}
