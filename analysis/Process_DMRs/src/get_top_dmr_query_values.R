my.get.top.DMR.query.values <- function(my.dmr.list, comparison, region, minRank = 100){
    #The parameter 'uniq' sets the script to ensure that the DMRs analyzed
    # are differntially methylated in a given comparison relative to
    # all others (default when parings is set to NULL)
    # or a user defined set of comparisons based on an input matrix
    # of dimension comparison x comparison
    ###############################################
    geneset <- my.dmr.list[[comparison]][[region]]
    #Separate the table into tables of hypo and hypermethylated regions
    diff.col <- colnames(geneset)[grep('mean.diff',colnames(geneset))]
    cur.hypo.table <- geneset[geneset[,diff.col] < 0,]
    cur.hyper.table <- geneset[geneset[,diff.col] > 0,]
    
    #Rank the regions according to 'combined rank'
    cur.hypo.table <- cur.hypo.table[order(x = cur.hypo.table$combinedRank,decreasing = F,na.last = T),]
    cur.hyper.table <- cur.hyper.table[order(x = cur.hyper.table$combinedRank,decreasing = F,na.last = T),]
    
    #make sure the regions are unique
    cur.hypo.table <- cur.hypo.table[!duplicated(cur.hypo.table$combinedRank),]
    cur.hyper.table <- cur.hyper.table[!duplicated(cur.hyper.table$combinedRank),]
    
    #Select just the best n ranked regions
    cur.hypo.table.top <- cur.hypo.table[1:minRank,]
    cur.hyper.table.top <- cur.hyper.table[1:minRank,]
    
    return(list(hypoTopUniqRank = cur.hypo.table.top,
                hyperTopUniqRank = cur.hyper.table.top,
                hypoAllUniqRank = cur.hypo.table,
                hyperAllUniqRank = cur.hyper.table))
}
