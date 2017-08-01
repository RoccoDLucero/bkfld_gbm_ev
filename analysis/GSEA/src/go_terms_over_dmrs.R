##Rocco Lucero August 2 2016
##This script should use RNBEADS DMR output object from e.g. script "GetDMRs_9comps_July19_2016"
## Here we will produce the lists of GO terms and genes over DMRs at various annotation regions 
## 
library(RnBeads)
library(biomaRt)
library(org.Hs.eg.db)
library(annotate)

source('my.get.top.DMR.queries.R')
source('my.get.top.DMR.queries.2.R')
#source('my.get.top.DMR.queries.any.R')
source('my.get.GO.gene.sets.R')

#make sure my.final.annotslist object is available
#Using: Rnbeads hg19 annotations
if(!exists(x = 'my.dmr.list')){
    #From 'GetDMRs_9comps_July19_2016.R' ->... -> 'FinalAnnotDMRTables.R'
    #Only the regions covered on the 450K array are included
    my.dmr.list <- readRDS(file = './input/my.annots.final') 
    print("loading DMRs and annotations")
}

################################################################################
#Here we get the gene symbols for both unique non-unique DMRs
if(F){
minRank <- 100
comparisons <- names(my.dmr.list)[1]
AllGoResults <- vector('list',length(comparisons))
names(AllGoResults) <- comparisons

for(cmp in comparisons){
    
    regions <- names(my.dmr.list[[cmp]])[1]
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
}

################################################################################
################################################################################
## Here we use my.get.top.DMR.queries.2 to get the GO Terms over regions
## that are uniquely differentially methylated in one comparison group versus all other
## comparison groups in the analysis
if(F){ #All of this code needs to be tuned for use with the 'my.get.top...2' function
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
}