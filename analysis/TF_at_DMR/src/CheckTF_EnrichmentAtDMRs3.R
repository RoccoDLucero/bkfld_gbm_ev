##Rocco Lucero August 4 2016
##
##Here we will look for Transcription Factor Binding Sites
## or MOTIFS That enrich near the Annotation sets and our DMRS.
library(grid)
library(parallel)
source("https://bioconductor.org/biocLite.R")
dep.pckgs <- c('BiocGenerics', 'Biostrings','MotifDb','testthat',
           'gtools','BiocStyle', 'knitr')
biocLite(dep.pckgs,suppressUpdates = T)
pckgs <- c('PWMEnrich','PWMEnrich.Hsapiens.background')
biocLite(pckgs,suppressUpdates = T)
library(PWMEnrich)
library(PWMEnrich.Hsapiens.background)
data("PWMCutoff4.hg19.MotifDb.Hsap")
data("PWMLogn.hg19.MotifDb.Hsap")

if(!exists(x = 'my.dmrSeq.list')){
    #From 'GetDMRs_9comps_July19_2016.R' ->... -> 'FinalAnnotDMRTables.R'
    #Only the regions covered on the 450K array are included
    my.dmrSeq.list <- readRDS(file = './output/All.DMR.Seqs.RDS') 
    print("loading DMRs, sequences, and annotations")
}


#Perform Motif enrichment for the hyper and hypo methylated DMRs for each comparison
enhancer.all.bg.seqs <- c(my.dmrSeq.list$EBM.EGM$blood_vessel_endothelial_cell.enh$univ$AnnotationSequence,
                          my.dmrSeq.list$EBM.EGM$T_cell.enh$univ$AnnotationSequence,
                          my.dmrSeq.list$EBM.EGM$astrocyte.enh$univ$AnnotationSequence,
                          my.dmrSeq.list$EBM.EGM$skeletal_muscle_cell.enh$univ$AnnotationSequence,
                          my.dmrSeq.list$EBM.EGM$skin_fibroblast.enh$univ$AnnotationSequence)
enhancer.all.bg.seqs <- unique(enhancer.all.bg.seqs)

my.enhancers.logn.bg1 <- makeBackground(motifs = PWMLogn.hg19.MotifDb.Hsap$pwms,
                                       bg.seq = enhancer.all.bg.seqs,
                                       algorithm = 'human')

if(T){
    comparisons <- names(my.dmrSeq.list)
    my.cmp.motif.enr <- vector('list',length(comparisons))
    names(my.cmp.motif.enr) <- comparisons
    
    for(cmp in comparisons){
        regions <- names(my.dmrSeq.list[[cmp]])
        print(paste('current comparison is',cmp))

        seqsHypo <- c()
        seqsHyper <- c()
        reg.motif.enr <- vector('list', 2)
        names(reg.motif.enr) <- c('Hypo','Hyper')
        print('Concatenating sequences by DMR directionality')
        for(rg in regions[1:2]){  # For enhancers only use[3:length(regions)]
            print(rg)
            #hypo <- (my.dmrSeq.list[[cmp]][[rg]]$hypo$AnnotationSequence)
            #hyper <- (my.dmrSeq.list[[cmp]][[rg]]$hyper$AnnotationSequence)
            hypo <- unique(my.dmrSeq.list[[cmp]][[rg]]$hypo$AnnotationSequence)
            hyper <- unique(my.dmrSeq.list[[cmp]][[rg]]$hyper$AnnotationSequence)
            
            seqsHypo <- c(seqsHypo, hypo )
            seqsHyper <- c(seqsHyper, hyper)
            #print('done')
        }
        seqsHypo <- do.call(c,seqsHypo)
        seqsHyper <- do.call(c,seqsHyper)
        #Motif Enrichment for HypoMethylated Set
        reg.motif.enr[['Hypo']]  <-  motifEnrichment(sequences = seqsHypo,
                                                     pwms = PWMLogn.hg19.MotifDb.Hsap    )#,
                                                     #bg = my.enhancers.logn.bg1@pwms)
        #Motif Enrichment for HyperMethylated Set
        reg.motif.enr[['Hyper']] <-  motifEnrichment(sequences = seqsHyper,
                                                     pwms = PWMLogn.hg19.MotifDb.Hsap    )#,
                                                     #bg = my.enhancers.logn.bg1@pwms)
        my.cmp.motif.enr[[cmp]] <- reg.motif.enr
    }   
#saveRDS(my.cmp.motif.enr, './output/All.cmp.genesAndProm.motifEnr.test.RDS')
my.cmp.motif.enr<-  readRDS('./output/All.cmp.ehnancerOnly.motifEnr.RDS')
}    

comparisons <- names(my.cmp.motif.enr)
my.enrichment.reports <- vector('list', length(comparisons))
names(my.enrichment.reports) <- comparisons
for(cmp in comparisons){
    cmp.reports <- vector('list', 2)
    names(cmp.reports) <- c('Hypo','Hyper')
    
    cmp.reports[['Hypo']] <- groupReport(my.cmp.motif.enr[[cmp]][['Hypo']])
    cmp.reports[['Hyper']] <- groupReport(my.cmp.motif.enr[[cmp]][['Hyper']])
    
    my.enrichment.reports[[cmp]] <- cmp.reports
}

maxRank = 500
allmotifs <- c()
for(cmp in comparisons){
    if(F){
        print(cmp)
        print('Hypo')
        print(my.enrichment.reports[[cmp]]$Hypo$target[1:25])
        print(max(my.enrichment.reports[[cmp]]$Hypo$p.value[1:25]))
        print(cmp)
        print('Hyper')
        print(my.enrichment.reports[[cmp]]$Hyper$target[1:25])
        print(max(my.enrichment.reports[[cmp]]$Hyper$p.value[1:25]))
    }
    #a <-(setdiff(my.enrichment.reports[[cmp]]$Hypo$target[1:25],my.enrichment.reports[[cmp]]$Hyper$target[1:25]))
    b <-c(my.enrichment.reports[[cmp]]$Hyper$target[1:maxRank],my.enrichment.reports[[cmp]]$Hypo$target[1:maxRank])
    allmotifs <- c(allmotifs,b)
    #allmotifs <- unique(allmotifs)
    aa <- allmotifs
    aa <- aa[!(duplicated(aa)|duplicated(aa,fromLast = T))]
    table(aa)
}

aa <- allmotifs
aa <- aa[!(duplicated(aa)|duplicated(aa,fromLast = T))]
table(aa)



for(cmp in comparisons){
    print(cmp)
    print(intersect(aa,my.enrichment.reports[[cmp]]$Hypo$target[1:maxRank]))
    print(intersect(aa,my.enrichment.reports[[cmp]]$Hyper$target[1:maxRank]))
}



#
mdb = query(MotifDb, 'UW.Motif.0156')
mdb
cc <- values(mdb) 
cc$geneSymbol
#I may need to find a new way to get gene symbols from motifs...



# perform motif enrichment!
res = motifEnrichment(sequence, PWMLogn.hg19.MotifDb.Hsap)
report <- sequenceReport(res,1)
report
pdf(file = './output/testplt.pdf',width = 15,height = 10)
plot(report[1:10],fontsize=5)
dev.off()



##SCRATCH
dd <- c(my.dmrSeq.list[[cmp]][[1]]$hypo$AnnotationSequence,
        my.dmrSeq.list[[cmp]][[2]]$hypo$AnnotationSequence,
        my.dmrSeq.list[[cmp]][[3]]$hypo$AnnotationSequence)
dd

dd <-as.character(sequenceReport(my.cmp.motif.enr$EBM.GBM8EVpel$Hyper,1)$target)
intersect(dd,miRNA_targ_toptargets4.5.2016$V1)

groupReport(my.cmp.motif.enr$EBM.GBM8EVpel$Hyper)$target[1:25]
