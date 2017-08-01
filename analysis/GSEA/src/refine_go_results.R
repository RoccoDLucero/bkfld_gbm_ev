#refine.GO.results.R
#This script is used to refine the results of the GO analysis


##The code here should be functionalized and 
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
