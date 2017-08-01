####To get the coordinates of annotation regions covered on the 450K array use...
library(RnBeads)
my.probes450 <- rnb.get.annotation(type = 'probes450')
my.reg.types <- rnb.region.types.for.analysis(my.rnb.set)
my.final.annots <- vector("list", length(my.dmr.list))
names(my.final.annots) <- names(my.dmr.list)

for(cur.comp.dmrs in names(my.dmr.list)){
    my.comp.annots <-  vector("list", length(my.reg.types))
    names(my.comp.annots) <- my.reg.types
    comp.dmrs <- my.dmr.list[[cur.comp.dmrs]] #This is a diffmeth object
 for(reg.type in my.reg.types){
    #dim(my.enhancer.profiles[[reg.type]])
    annot.tab <- cbind(annotation(my.rnb.set, type = reg.type,add.names = F),
          get.table(comp.dmrs,comp.dmrs@comparisons,reg.type))
    
    #use GRanges tools to compute probe overlap per region
    ann <- annotation(my.rnb.set, type = reg.type,add.names = F)
    ann <- makeGRangesFromDataFrame(ann)
    probes.covered <- (countOverlaps(ann,my.probes450[[1]]))
    
    annot.tab <- cbind(annot.tab,probes.covered)
    my.comp.annots[[reg.type]] <- annot.tab
 }
    my.final.annots[[cur.comp.dmrs]] <-my.comp.annots    
}

saveRDS(my.final.annots, file = paste(report.dir,'/my.annots.final',sep = ''))
