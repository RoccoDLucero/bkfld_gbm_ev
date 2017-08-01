#This script will be used to 
# Get hg19 sequence data in the annotation regions of interest
# based on annotation coordinates

# Load the BSgenome package and download the hg19 reference sequence  
# For documentation see http://www.bioconductor.org/packages/release/bioc/html/BSgenome.html 
source(file = '../GSEA/my.get.GO.gene.sets.R')
source(file =  'my.get.top.DMR.sequences.R')
source("http://www.bioconductor.org/biocLite.R") 
biocLite("BSgenome",suppressUpdates = T) 
biocLite("BSgenome.Hsapiens.UCSC.hg19",suppressUpdates = T) #installs the human genome (~850 MB download). 
library('BSgenome.Hsapiens.UCSC.hg19') 
library('MotifDb')
#use the final annots object produced by 
#processing the RNBEads output
#This should exclude all regions not covered on the 450K array
if(!exists(x = 'my.dmr.list')){
    #From 'GetDMRs_9comps_July19_2016.R' ->... -> 'FinalAnnotDMRTables.R'
    #Only the regions covered on the 450K array are included
    my.dmr.list <- readRDS(file = './input/my.annots.final') 
    print("loading DMRs and annotations")
}

#For each comparison in each region
#for each annotation
#pull the genomic sequence
#Turn the sequence into a DNAStringSet" object from(attr(,"package") [1] "Biostrings")
#and add these to the annots object.


if(T){
    minRank <- 100
    comparisons <- names(my.dmr.list)
   
    for(cmp in comparisons){
        
        regions <- names(my.dmr.list[[cmp]])
        print(paste('current comparison is',cmp))
       
        for(rg in regions){
            print('Getting Query Gene Sequences')
            #add the sequences before subsetting any further:
            #make sure the sequence object comes out as class DNAStringset
            
            #Add the sequence information to each entry of the DMR list
            my.dmr.list[[cmp]][[rg]]$AnnotationSequence  <- getSeq(Hsapiens, my.dmr.list[[cmp]][[rg]]$Chromosome,
                                                                   my.dmr.list[[cmp]][[rg]]$Start,
                                                                   my.dmr.list[[cmp]][[rg]]$End) 
            
            #Get the best ranked regions for hyper/hypo methylation
            # and replace the entry at this index
            my.dmr.list[[cmp]][[rg]] <- my.get.top.DMR.sequences(my.dmr.list = my.dmr.list,
                                                         comparison = cmp,
                                                         region = rg,
                                                         minRank = minRank)
        }
        
    }   
    
    saveRDS(my.dmr.list, './output/All.DMR.Seqs.RDS')
    my.dmrSeq.list <- my.dmr.list #Rename the object
    rm(my.dmr.list)
}    
