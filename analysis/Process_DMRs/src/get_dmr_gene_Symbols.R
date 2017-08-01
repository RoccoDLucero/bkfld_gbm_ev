library(RnBeads)
library(biomaRt)
library(org.Hs.eg.db)
library(annotate)

source('my.get.top.DMR.queries.R')
#source('my.get.top.DMR.queries.2.R')
#source('my.get.top.DMR.queries.any.R')
source('my.get.GO.gene.sets.R')
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

################################################################################
#Here we get the gene symbols for both unique non-unique DMRs
if(F){
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

################################################################################
################################################################################
if(F){
    #AllGoResults <- readRDS(paste(report.dir,'/AllGO.Results1.RDS',sep = ''))    
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


################################################################################
################################################################################
## Here we generate a file that contains the gene descriptions unique to each
## comparison by region by gain/loss of methylation 
##
if(F){
    #AllGoResults <- readRDS(paste(report.dir,'/AllGO.Results1.RDS',sep = ''))
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
            
            write(x = 'Methylation Gained', 
                  file = my.file,
                  append = T, sep = '\n')
            
            write.table((unique(AllGoResults[[cmp]][[rg]]$hypoSet$query$description)),
                        file = my.file,row.names = F,
                        append = T, sep = '\t')
            
            write(x = 'Methylation Lost', 
                  file = my.file,
                  append = T, sep = '\n')
            
            write.table((unique(AllGoResults[[cmp]][[rg]]$hyperSet$query$description)),
                        file = my.file,row.names = F,
                        append = T, sep = '\t')
        }
    }
}

################################################################################
################################################################################
## Here we write a file that contains the lists of unique GO terms returned for
## each comparison for each region for methyaltion gain/loss versus
## all others comparisons

#AllGoResults <- readRDS(paste(report.dir,'/AllGO.Results1.RDS',sep = ''))
getTerms = T
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
            #hypoTerms <- c(hypoTerms, AllGoResults[[cmp]][[rg]]$summaryHypo)
            #hyperTerms <- c(hyperTerms, AllGoResults[[cmp]][[rg]]$summaryHyper)
                        
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
        names(focus) <- c('Gained Methylation','Lost Methylation')
        uniqTerms <- setdiff(focus, bckg)
        print(length(uniqTerms))
        write.table(uniqTerms,file = paste('./output/TermsUniqToCmp_2',subj,
                                           '.txt',sep = ''),
                    sep = '\t',quote = FALSE,row.names = FALSE)
    }
}
################################################################################
################################################################################
##
##
##

AllGoResults <- readRDS(paste(report.dir,'/AllGO.Results1.RDS',sep = ''))
getTerms = T
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
        #hypoTerms <- c(hypoTerms, AllGoResults[[cmp]][[rg]]$summaryHypo)
        #hyperTerms <- c(hyperTerms, AllGoResults[[cmp]][[rg]]$summaryHyper)
        
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
    names(focus) <- c('Gained Methylation','Lost Methylation')
    uniqTerms <- setdiff(focus, bckg)
    print(length(uniqTerms))
    write.table(uniqTerms,file = paste('./output/TermsUniqToCmp_2',subj,
                                       '.txt',sep = ''),
                sep = '\t',quote = FALSE,row.names = FALSE)
}
}
################################################################################
################################################################################
## Here we write a file that contains the lists of unique GENE descriptions returned for
## each comparison for each region for methyaltion gain/loss versus
## all other comparisons
getDiffGenes = T
if(getDiffGenes){
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
        #print(length(unique(focus)))
        #print(length(uniqTerms))
        #print(length(uniqTerms)/length(unique(focus)))
        write.table(uniqTerms,file = paste('./output/GenesUniqToCmp',subj,
                                           '.txt',sep = ''),
                    sep = '\t',quote = FALSE,row.names = FALSE)
    }
}
################################################################################
################################################################################
## Here we look for overlap in the annotation over DMRs between a defined pair of
## comparisons. 
##
##
getSharedGenes = F
subj1 <- "EBM.EGM"
subj2 <- "EBM.GBM8EVpel"

if(getSharedGenes){
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
            #hypoTerms <- c(hypoTerms, AllGoResults[[cmp]][[rg]]$hypoSet$univ)
            #hyperTerms <- c(hyperTerms, AllGoResults[[cmp]][[rg]]$hyperSet$univ)
            
            
        }
        hypoAll[[cmp]] <- hypoTerms
        hyperAll[[cmp]] <- hyperTerms
    }
    
    
    subj1genes <- c(hypoAll[[subj1]]$description, hyperAll[[subj1]]$description)
    subj2genes <- c(hypoAll[[subj2]]$description, hyperAll[[subj2]]$description)
        
    shared <- intersect(subj1genes,subj2genes)   
        
    write.table(shared,file = paste('./output/GenesSharedIn',subj1, subj2,
                                           '.txt',sep = ''),
                    sep = '\t',quote = FALSE,row.names = FALSE)
}
################################################################################
################################################################################
