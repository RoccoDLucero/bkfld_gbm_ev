#This script should find unque DMRs for each comparison group versus
# user defined set of other comparison groups
#this requires setting a threshold on the combined rank 
#The output of this script will be a reduced set of DMRs per treatment,
# that should be input to the GOBP analysis

analysis.dir = "./RnBeads/analysis"
report.dir = file.path(analysis.dir, "GetDMRs_9comps_July19_2016")
if(!exists(x = 'my.dmr.list')){
    my.dmr.list <- readRDS(file = paste(report.dir,'/my.annots.final',sep = ''))
    print("loading annotations")
} else{if(exists('my.final.annots')){my.dmr.list <- my.final.annots}}

#first concatenate all the DMRs for a given comparison and get the unique
#coordinates in a single vector. To to this concatenate the coordinates into
#single values

rankThreshold = 1000
comps <- names(my.dmr.list)
regions <- names(my.dmr.list$EBM.EGM)

for(cmp in comps){
    for(rg in regions){
        my.dmr.list[[cmp]][[rg]]$Chromosome
        my.dmr.list[[cmp]][[rg]]$Start
        my.dmr.list[[cmp]][[rg]]$End
    }
}
