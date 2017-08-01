#This function takes in a differential methylation data table (RnBeads output)
#and outputs the top sites or regions based on the combindeRank column.
my.getTopRanked = function(difMethTab, feature.type = 'region' ,rankCut = 250, flankSize = 10000, fileOut = F){
    require(gtools)
    difMethTab = difMethTab[order(difMethTab$combinedRank),]
    topDM.roi = difMethTab[1:rankCut,]
    
    topDM.roi$Start = topDM.roi$Start - flankSize
    if(feature.type == 'region'){topDM.roi$End = (topDM.roi$End + flankSize)}
    if(feature.type == 'site'){topDM.roi$End = (topDM.roi$Start + flankSize)}
    
    #Partially Sort by chromosome
    topDM.roi = topDM.roi[mixedorder(paste(topDM.roi$Chromosome,topDM.roi$Start)),]
    topDM.roi = topDM.roi[,c('Chromosome','Start','End')]
    #Write table to file
    if(fileOut!=F){
        write.table(x = topDM.roi,file = fileOut,
                    quote = F,row.names = F,sep = '\t')
    }
    #Return the narrowed table
    topDM.roi
}