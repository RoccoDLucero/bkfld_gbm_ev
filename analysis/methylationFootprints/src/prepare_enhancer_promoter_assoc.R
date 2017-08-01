#This function imports and prepares the PRESSTO enhancer promoter maps
#described in 'An atlas of active enhancers across human cell types and tissues'
#A data frame describing the genomic span of the interaction, the enhancer location,
#and identifiers for the promoter is returned

My.Prepare.Enh.Prom.Assoc.Data <- function(){
    #these enhancer promoter maps are from GRCH37/HG19 (Pressto)
    
    enh.refseqProm.assoc = read.table('./input/enhancer_tss_associations.bed',
                                      sep = '\t', header = T, comment.char = '!' ) 
    head(enh.refseqProm.assoc)
    #split the 'name' field into fields {'Enh.coord','Ref.seq.id','symbol','express.corrR','FDR'}
    expected.fields = 5
    assoc.info = enh.refseqProm.assoc$name
    assoc.info = as.character(assoc.info)
    assoc.info = strsplit(assoc.info,';')
    assoc.fields = sapply(X = assoc.info,FUN = length)
    #select records for which name field has all 5 entries
    assoc.info = assoc.info[assoc.fields==expected.fields]
    enh.refseqProm.assoc = enh.refseqProm.assoc[assoc.fields==expected.fields,]
    
    #merge the split enhancer promoter association data back into the original data frame
    assoc.info = as.data.frame(matrix(data = unlist(assoc.info),ncol = expected.fields,byrow = T))
    enh.refseqProm.assoc = cbind(enh.refseqProm.assoc[,-4],assoc.info)
    colnames(enh.refseqProm.assoc) = c( c("chromosome", "start", "end"),
                                        colnames(enh.refseqProm.assoc)[4:11],
                                        c('Enh.coord','Ref.seq.id','symbol','express.corrR','FDR'))
    
    return(enh.refseqProm.assoc)
}
