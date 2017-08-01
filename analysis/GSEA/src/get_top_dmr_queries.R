my.get.top.DMR.queries.2 <- function(dmr.list, comparison1, comparison2 = NULL, region,
                                     minRank = 100, diff.cutoff = NULL){

    library(biomaRt)
    library(org.Hs.eg.db)
    library(annotate)
    listEnsembl()
    listEnsembl(GRCh=37)
    grch37 = useEnsembl(biomart="ensembl",GRCh=37)
    listDatasets(grch37)[31:35,]
    ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    head(listFilters(ensembl))
    
#####PLACE this outside of the function    
    symbolMapping <- read.delim("H:/Dropbox/BRL/proj/bkfld_gbm_ev/analysis/GSEA/input/symbolMapping.txt")
    symbolMapping <- symbolMapping[,c("Approved.Symbol","Ensembl.Gene.ID")]
    symbolMapping <- symbolMapping[symbolMapping$Ensembl.Gene.ID != "" ,]
#####    
    ###############################################
    
    get.dmr.table <- function(dmr.list,comparison,region){
        if(is.null(comparison)){return("Comparison Null")}
        geneset <- dmr.list[[comparison]][[region]]
        geneset$LocID <- paste(geneset$Chromosome,geneset$Start,geneset$End, sep = ':')
    
        #Separate the table into tables of hypo and hypermethylated regions
        diff.col <- colnames(geneset)[grep('mean.diff',colnames(geneset))]
        cur.hypo.table  <- geneset[geneset[,diff.col] < 0,]
        cur.hyper.table <- geneset[geneset[,diff.col] > 0,]
        
        #Rank the regions according to 'combined rank'
        cur.hypo.table <- cur.hypo.table[order(x = cur.hypo.table$combinedRank,decreasing = F,na.last = T),]
        cur.hyper.table <- cur.hyper.table[order(x = cur.hyper.table$combinedRank,decreasing = F,na.last = T),]
        
        dmrTables <-list(cur.hyper.table,cur.hypo.table)
        names(dmrTables) <- c('lost_methylayion','gained_methylation')
        return(dmrTables)
    }
    
    tableList1 <- get.dmr.table(my.dmr.list,comparison1,region)
    tableList2 <-get.dmr.table(my.dmr.list,comparison2,region)
        
    if(is.null(comparison2)){
        #Select just the best n ranked regions
        lost <- tableList1$lost_methylayion[1:minRank,]
        gained <- tableList1$gained_methylation[1:minRank,]
        
        
        lost$updateSymbol <- lost[(rownames(lost) %in% symbolMapping$Ensembl.Gene.ID),] 
        
        
        assigned.reg.col <- colnames(dmr.list[[comparison1]][[region]])
        assigned.reg.col <- assigned.reg.col[grep('symbol',assigned.reg.col)]
        
        ensembl_univ <- geneset1[,assigned.reg.col]
        ensembl_query.gained <- gained[,assigned.reg.col]
        ensembl_query.lost <- lost[,assigned.reg.col]
        #return
       xx <- (list(gained.me = ensembl_query.gained,
                    lost.me = ensembl_query.lost,
                    univ = ensembl_univ ))
    }
    
    if(!is.null(comparison2)){
        #Get the unique DMR locations for geneset1
        cur.hypo.rows <- setdiff(cur.hypo.table1$LocID,cur.hypo.table2$LocID)
        cur.hyper.rows <- setdiff(cur.hyper.table1$LocID,cur.hyper.table2$LocID)
        
        cur.hypo.table1 <- cur.hypo.table1[cur.hypo.rows,]
        cur.hyper.table1 <- cur.hyper.table1[cur.hyper.rows,]
        
        #This parameter can be used to set a minimum methylation difference
        if(!is.null(diff.cutoff)){
            cur.hypo.table1 <- cur.hypo.table1[abs(cur.hypo.table1$mean.mean.diff) > diff.cutoff,]
            cur.hyper.table1 <- cur.hyper.table1[abs(cur.hyper.table1$mean.mean.diff) > diff.cutoff,]
        }
            
        #For the remaining DMRS take as many as the minRank allows, or if
        # fewer are available take everything leftover
        cur.hypo.table1 <- cur.hypo.table1[1:min(nrow(cur.hypo.table1),minRank),]
        cur.hyper.table1 <- cur.hyper.table1[1:min(nrow(cur.hyper.table1),minRank),]
        
        assigned.reg.col <- colnames(geneset1)[grep('symbol',colnames(geneset1))]
        ensembl_univ <- geneset1[,assigned.reg.col]
        ensembl_query.hypo <- cur.hypo.table1[,assigned.reg.col]
        ensembl_query.hyper <- cur.hyper.table1[,assigned.reg.col]
        return(list(hypo = ensembl_query.hypo,
                    hyper = ensembl_query.hyper,
                    univ = ensembl_univ ))
        
        
    }
    
    
    
}
