# Script for extracting highly ranked loci from RnBeads pipeline, identifying genes affected by gene / promoter regions
#setwd('~/Downloads/Breakefield/cmp1/results_.5/reports/differential_methylation_data')

#Compare targeted top 100 loci to untargeted top 100 loci
targeted.loci <- NULL
#sites
data <- as.matrix(read.csv("diffMethTable_site_cmp1_targeted.csv", header=TRUE, sep=","))
targeted.loci <-data[,c("cgid", "Chromosome", "Start", "Strand", "mean.diff", "mean.quot.log2", "diffmeth.p.val","combinedRank")]
#targeted.loci <- targeted.loci[intersect(which(as.numeric(targeted.loci[,"diffmeth.p.val"])<0.05), 
#                                                which(abs(as.numeric(targeted.loci[,"mean.diff"]))>0.3)),]
targeted.loci[,"mean.diff"] <- round(as.numeric(targeted.loci[,"mean.diff"]),2)
targeted.loci[,"mean.quot.log2"] <- round(as.numeric(targeted.loci[,"mean.quot.log2"]),2)
targeted.loci[,"diffmeth.p.val"] <- format(as.numeric(targeted.loci[,"diffmeth.p.val"]),width=4, digits=2)
targeted.loci <- targeted.loci[order(as.numeric(targeted.loci[,"combinedRank"])),]

#setwd('~/Downloads/Breakefield/cmp1/results/reports/differential_methylation_data')
untargeted.loci <- NULL
#sites
data <- as.matrix(read.csv("diffMethTable_site_cmp1_untargeted.csv", header=TRUE, sep=","))
untargeted.loci <-data[,c("cgid", "Chromosome", "Start", "Strand", "mean.diff", "mean.quot.log2", "diffmeth.p.val","combinedRank")]
#untargeted.loci <- untargeted.loci[intersect(which(as.numeric(untargeted.loci[,"diffmeth.p.val"])<0.05), 
#                                         which(abs(as.numeric(untargeted.loci[,"mean.diff"]))>0.3)),]
untargeted.loci[,"mean.diff"] <- round(as.numeric(untargeted.loci[,"mean.diff"]),2)
untargeted.loci[,"mean.quot.log2"] <- round(as.numeric(untargeted.loci[,"mean.quot.log2"]),2)
untargeted.loci[,"diffmeth.p.val"] <- format(as.numeric(untargeted.loci[,"diffmeth.p.val"]),width=4, digits=2)
untargeted.loci <- untargeted.loci[order(as.numeric(untargeted.loci[,"combinedRank"])),]

df <- data.frame(targeted=as.matrix(targeted.loci[1:100,1]),
                 untargeted=as.matrix(untargeted.loci[1:100,1]))
setwd('~/Downloads/Breakefield/cmp1/')
write.table(df, "EBM-EVs.txt", quote=FALSE, row.names=FALSE)

load("~/Downloads/MethylationChipAnnotation.rda")
MethylAnno[which(MethylAnno[,1]%in% targeted.loci[1:100,1]), c("GeneSymbol", "UCSC_RefGene_Accession")]
write.table(unique(na.omit(MethylAnno[which(MethylAnno[,1]%in% targeted.loci[which(as.numeric(targeted.loci[,"combinedRank"])<500),1]), "UCSC_RefGene_Accession"])), "targeted.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(unique(na.omit(MethylAnno[which(MethylAnno[,1]%in% untargeted.loci[which(as.numeric(untargeted.loci[,"combinedRank"])<500),1]), "UCSC_RefGene_Accession"])), "untargeted.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

write.table(unique(na.omit(MethylAnno[which(MethylAnno[,"IlmnID"] %in% data[,"cgid"]),
                                      "UCSC_RefGene_Accession"])), "background.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)




# Set Operation -----------------------------------------------------------

#[(A-C) intersect (A-E)] \ [(B-C) U (D-C) U (E-C) U (F-C) U (B-F)]
# cmp2 intersect cmp17 \ cmp3 U cmp1 U cmp6 U cmp7 U cmp22
#setdiff(intersect(significant.genes.cmp2, significant.genes.cmp17), 
#         unique(c(significant.genes.cmp3, significant.genes.cmp1, 
#               significant.genes.cmp6, significant.genes.cmp7,significant.genes.cmp22)))