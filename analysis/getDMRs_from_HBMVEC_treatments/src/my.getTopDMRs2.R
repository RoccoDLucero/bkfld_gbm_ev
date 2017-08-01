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
    #topDM.roi = topDM.roi[mixedorder(paste(topDM.roi$Chromosome,topDM.roi$Start)),]
    topDM.roi = topDM.roi[,c('Chromosome','Start','End')]
    
    #Generate the Null Set as regions that did not count as DMRs
    NON.DM.roi = difMethTab[((rankCut+1):nrow(difMethTab)),]
    NON.DM.roi$Start = NON.DM.roi$Start - flankSize
    if(feature.type == 'region'){NON.DM.roi$End = (NON.DM.roi$End + flankSize)}
    if(feature.type == 'site'){NON.DM.roi$End = (NON.DM.roi$Start + flankSize)}
    #Partially Sort by chromosome
    #NON.DM.roi = NON.DM.roi[mixedorder(paste(NON.DM.roi$Chromosome,NON.DM.roi$Start)),]
    NON.DM.roi = NON.DM.roi[,c('Chromosome','Start','End')]
    
    #Write table to file
    if(fileOut!=F){
        #Out = deparse(substitute(difMethTab))
        write.table(x = topDM.roi,file = paste('./top',rankCut,fileOut,'.bed',sep = ''),
                    quote = F,row.names = F,sep = '\t')
        write.table(x = NON.DM.roi,file = paste('./top',rankCut,fileOut,'_NON_DMR.bed',sep = ''),
                    quote = F,row.names = F,sep = '\t')
    }
    #Return the narrowed table
    list(topDM.roi,NON.DM.roi)
}