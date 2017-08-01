my.get.GO.gene.sets <- function(query, univ, ontology = 'BP', assembly = 'hg19'){
    library(RnBeads)
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
