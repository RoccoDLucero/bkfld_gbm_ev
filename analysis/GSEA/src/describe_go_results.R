#describe.results.R
#This script should generate descriptive statistics of the DMR/GO analyis#

comps <- names(AllGoResults)[c(1:3,7:8)] #Make sure to disclude comparisons involving GBM8 Cells
regions <- names(AllGoResults[[1]])

ups <- c()
downs <- c()
nms <- c()
downQ <- c()
upQ <- c()
for(cmp in comps){
    for(rg in regions){
        nms <- c(nms,paste(cmp,rg)) 
        downs <- c(downs,
                   (length(unique(AllGoResults[[cmp]][[rg]]$hyperSet$query$ensembl_gene_id))))
        ups <- c(ups,(length(unique(AllGoResults[[cmp]][[rg]]$hypoSet$query$ensembl_gene_id))))
        downQ <- c(downQ,
                   (length(unique(AllGoResults[[cmp]][[rg]]$hyperSet$univ$ensembl_gene_id))))
        upQ <- c(upQ,(length(unique(AllGoResults[[cmp]][[rg]]$hypoSet$univ$ensembl_gene_id))))
    }
}

hist(downs, breaks = seq(0,(max(downs)+10),10), col = 'red')
hist(ups, breaks = seq(0,(max(ups)+10),10), col = 'blue')
dn <- matrix(data = downs/downQ,nrow = 5,byrow = T)
up <- matrix(data = ups/upQ, nrow = 5, byrow = T)
dn <- matrix(data = downs,nrow = 5,byrow = T)
up <- matrix(data = ups, nrow = 5, byrow = T)
colnames(up) <- regions
rownames(up) <- comps
up

colnames(dn) <- regions
rownames(dn) <- comps
dn

hist(dn, breaks = seq(0,(max(dn)+.01),.01), col = 'red')
hist(up, breaks = seq(0,(max(up)+.01),.01), col = 'blue')


numDMRs <- rbind(ups,downs)
colnames(numDMRs) <- nms
numDMRs
