library(RnBeads)
library(biomaRt)
library(org.Hs.eg.db)
library(annotate)
#source(file = './processPRESSTOmapping.R')

#This needs to be called in the environment of the GETDMRsJuly_2016 (or later) script##
#make sure my.final.annotslist objedct is available
#Rnbeads hg19 annotations

analysis.dir = "./RnBeads/analysis"
report.dir = file.path(analysis.dir, "GetDMRs_9comps_July19_2016")
if(!exists(x = 'my.dmr.list')){
    my.dmr.list <- readRDS(file = paste(report.dir,'/my.annots.final',sep = ''))
    print("loading annotations")
} else{if(exists('my.final.annots')){my.dmr.list <- my.final.annots}}


my.get.top.DMR.queries <- function(my.dmr.list, comparison, region, minRank = 100){
    geneset <- my.dmr.list[[comparison]][[region]]
    #Separate the table into tables of hypo and hypermethylated regions
    diff.col <- colnames(geneset)[grep('mean.diff',colnames(geneset))]
    cur.hypo.table <- geneset[geneset[,diff.col] < 0,]
    cur.hyper.table <- geneset[geneset[,diff.col] > 0,]
    
    #Rank the regions according to 'combined rank'
    cur.hypo.table <- cur.hypo.table[order(x = cur.hypo.table$combinedRank,decreasing = F,na.last = T),]
    cur.hyper.table <- cur.hyper.table[order(x = cur.hyper.table$combinedRank,decreasing = F,na.last = T),]
    
    #Select just the best n ranked regions
    cur.hypo.table <- cur.hypo.table[1:minRank,]
    cur.hyper.table <- cur.hyper.table[1:minRank,]
    
    assigned.promoter.col <- colnames(geneset)[grep('symbol',colnames(geneset))]
    ensembl_univ <- geneset[,assigned.promoter.col]
    ensembl_query.hypo <- cur.hypo.table[,assigned.promoter.col]
    ensembl_query.hyper <- cur.hyper.table[,assigned.promoter.col]
    return(list(hypo = ensembl_query.hypo,
                hyper = ensembl_query.hyper,
                univ = ensembl_univ ))
}

my.get.GO.gene.sets <- function(query, univ, ontology = 'BP', assembly = 'hg19'){
    mart<- useDataset(dataset = "hsapiens_gene_ensembl" , useMart("ensembl"))
    
    my.query <- getBM(filters= "hgnc_symbol",
                      attributes= c("ensembl_gene_id", "entrezgene", "description"),
                      values= query, mart= mart)
    
    my.univ <- getBM(filters= "hgnc_symbol",
                     attributes= c("ensembl_gene_id", "entrezgene", "description"),
                     values = univ, mart= mart)
    
    terms <- performGOenrichment.diffMeth.entrez(gids = na.omit(my.query$entrezgene),
                                                 uids = na.omit(my.univ$entrezgene),
                                                 ontology = ontology,assembly = assembly)
    
    
    return(list(query = my.query, univ = my.univ, GO_terms = terms))
}


minRank <- 100
comparisons <- names(my.dmr.list)
AllGoResults <- vector('list',length(comparisons))
names(AllGoResults) <- comparisons

for(cmp in comparisons){
    
    regions <- names(my.dmr.list[[cmp]])
    cmp.GO.Result <- vector(mode = 'list',length = length(regions))
    names(cmp.GO.Result) <-regions
    
    print(paste('current comparison is',cmp))
    
        for(rg in regions){
        print('Getting Query Gene Symbols')
        my.ensembl.queries <- my.get.top.DMR.queries(my.dmr.list = my.dmr.list,
                                                     comparison = cmp,
                                                     region = rg,
                                                     minRank = minRank)
        
        print('Performing GO Enrichment Analysis for region:')
        print(rg)
        hypoSet <- my.get.GO.gene.sets(query = my.ensembl.queries[['hypo']],
                                          univ = my.ensembl.queries[['univ']])
        
        hyperSet <- my.get.GO.gene.sets(query = my.ensembl.queries[['hyper']],
                                           univ = my.ensembl.queries[['univ']])
        
        print('Completed Enrichment Analysis for region:')
        print(rg)
        summaryHypo  <- summary(hypoSet$GO_terms)
        summaryHyper <- summary(hyperSet$GO_terms)
        print("generated GO BP term summaries")
        
        regionalResult <- list(summaryHypo,summaryHyper,hyperSet,hypoSet)
        names(regionalResult) <- c('summaryHypo','summaryHyper',
                                   'hyperSet','hypoSet')
        
        print('compiled regional result:')
        print(paste(cmp,rg))
        cmp.GO.Result[[rg]] <-regionalResult
        print('Added Regional Result to Comparison level list object')
    }
    
    AllGoResults[[cmp]] <- cmp.GO.Result
    print("Added Comparison level GO result to AllResults object")
    
}   

saveRDS(AllGoResults, paste(report.dir,'/AllGO.Results1.RDS',sep = ''))

#methDirection <- c("gained methylation","lost methylation")
############################################################################
if(F){
#my.GO.results <- readRDS(paste(report.dir,'/AllGO.Results.RDS',sep = ''))    
for(cmp in comparisons){
    
    regions <- names(AllGoResults[[cmp]])
    cmp.GO.Result <- vector(mode = 'list',length = length(regions))
    names(cmp.GO.Result) <-regions
    
    print(paste('current comparison is',cmp))
    write(x = paste('current comparison is',cmp), 
          file = './output/test_tables.txt',
          append = T, sep = '\t')
    
    for(rg in regions){
        print(paste(cmp,rg))
        write(x = paste(cmp,rg), 
              file = './output/test_tables.txt',
              append = T, sep = '\t')
        
        write.table(AllGoResults[[cmp]][[rg]]$summaryHypo,
                    file = './output/test_tables.txt',row.names = F,
                    append = T, sep = '\t')
        
        write.table(AllGoResults[[cmp]][[rg]]$summaryHyper,
                    file = './output/test_tables.txt', row.names = F,
                    append = T, sep = '\t')
    }
}
}

if(F){
    #my.GO.results <- readRDS(paste(report.dir,'/AllGO.Results.RDS',sep = ''))
    my.file = './output/querySets.txt'
    for(cmp in comparisons){
        
        regions <- names(AllGoResults[[cmp]])
        cmp.GO.Result <- vector(mode = 'list',length = length(regions))
        names(cmp.GO.Result) <-regions
        
        print(paste('current comparison is',cmp))
        write(x = paste('current comparison is',cmp), 
              file = my.file,
              append = T, sep = '\t')
        
        for(rg in regions){
            print(paste(cmp,rg))
            write(x = paste(cmp,rg), 
                  file = my.file,
                  append = T, sep = '\t')
            
            write.table(AllGoResults[[cmp]][[rg]]$hypoSet$query,
                        file = my.file,row.names = F,
                        append = T, sep = '\t')
            
            write.table(AllGoResults[[cmp]][[rg]]$hyperSet$query,
                        file = my.file,row.names = F,
                        append = T, sep = '\t')
        }
    }
}
getTerms = F
if(getTerms){
    hypoAll <- vector("list",length = length(comparisons))
    hyperAll <- vector("list",length = length(comparisons))
    names(hypoAll) <- comparisons
    names(hyperAll) <- comparisons
    for(cmp in comparisons){
        
        regions <- names(AllGoResults[[cmp]])
        cmp.GO.Result <- vector(mode = 'list',length = length(regions))
        names(cmp.GO.Result) <-regions
        
        print(paste('current comparison is',cmp))
        
        hypoTerms <- c()
        hyperTerms <- c()
        for(rg in regions){
            
            hypoTerms <- c(hypoTerms, AllGoResults[[cmp]][[rg]]$summaryHypo$Term)
            hyperTerms <- c(hyperTerms, AllGoResults[[cmp]][[rg]]$summaryHyper$Term)
                        
        }
        hypoAll[[cmp]] <- hypoTerms
        hyperAll[[cmp]] <- hyperTerms
    }
    
    
    excluded <- c("HBMVECinEBM.GBM8inEBM","HBMVECinEBM.GBM8inNB",
                  "EGM.GBM8Evpel","GBM8EVpel.GBM8EVsup")
    comparisons.restricted <- setdiff(comparisons,excluded)
    for(cmp in comparisons.restricted){
        subj <- comparisons.restricted[grep(cmp,comparisons.restricted)]
        rest <- setdiff(comparisons.restricted, subj)
        
        bckg <- c()
        for(i in rest){
            rest.hypo <- hypoAll[[i]]
            rest.hyper <- hyperAll[[i]]
            bckg <- c(bckg,rest.hypo, rest.hyper)
        }
        
        focus <- c(hypoAll[[subj]], hyperAll[[subj]])
        uniqTerms <- setdiff(focus, bckg)
        #print(length(uniqTerms))
        write.table(uniqTerms,file = paste('./output/TermsUniqToCmp',subj,
                                           '.txt',sep = ''),
                    sep = '\t',quote = FALSE,row.names = FALSE)
    }
}

getGenes = T
if(getGenes){
    hypoAll <- vector("list",length = length(comparisons))
    hyperAll <- vector("list",length = length(comparisons))
    names(hypoAll) <- comparisons
    names(hyperAll) <- comparisons
    for(cmp in comparisons){
        
        regions <- names(AllGoResults[[cmp]])
        cmp.GO.Result <- vector(mode = 'list',length = length(regions))
        names(cmp.GO.Result) <-regions
        
        print(paste('current comparison is',cmp))
        
        hypoTerms <- c()
        hyperTerms <- c()
        for(rg in regions){
            
            hypoTerms <- c(hypoTerms, AllGoResults[[cmp]][[rg]]$hypoSet$query)
            hyperTerms <- c(hyperTerms, AllGoResults[[cmp]][[rg]]$hyperSet$query)
            
        }
        hypoAll[[cmp]] <- hypoTerms
        hyperAll[[cmp]] <- hyperTerms
    }
    
    
    excluded <- c("HBMVECinEBM.GBM8inEBM","HBMVECinEBM.GBM8inNB",
                  "EGM.GBM8Evpel","GBM8EVpel.GBM8EVsup")
    comparisons.restricted <- setdiff(comparisons,excluded)
    for(cmp in comparisons.restricted){
        subj <- comparisons.restricted[grep(cmp,comparisons.restricted)]
        rest <- setdiff(comparisons.restricted, subj)
        
        bckg <- c()
        for(i in rest){
            rest.hypo <- hypoAll[[i]]$description
            rest.hyper <- hyperAll[[i]]$description
            bckg <- c(bckg,rest.hypo, rest.hyper)
        }
        
        focus <- c(hypoAll[[subj]]$description, hyperAll[[subj]]$description)
        uniqTerms <- setdiff(focus, bckg)
        #print(length(uniqTerms))
        write.table(uniqTerms,file = paste('./output/GenesUniqToCmp',subj,
                                           '.txt',sep = ''),
                    sep = '\t',quote = FALSE,row.names = FALSE)
    }
}
