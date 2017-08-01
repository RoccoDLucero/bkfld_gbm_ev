#breakefield450k <- read.table("FinalReport2.txt", header=TRUE, sep="\t", quote="")
#data450k <- breakefield450k[,c(1, as.numeric(grep(".*.Signal_A", colnames(breakefield450k))), 
#                               as.numeric(grep(".*.Signal_B", colnames(breakefield450k))),
#                               as.numeric(grep(".*.Detection.Pval", colnames(breakefield450k))))]
#rm(breakefield450k)
#rownames(data450k) <- data450k[,1]
#data450k <- data450k[,-1]
#colnames(data450k) <- gsub("zappulli", "", colnames(data450k))
#data450k <- data450k[,order(colnames(data450k))]
#data450k[1,]
#betas <- data450k[,seq(3,dim(data450k)[2],3)]/(data450k[,seq(2,dim(data450k)[2],3)]+data450k[,seq(3,dim(data450k)[2],3)])
#rm(data450k)

# Uncomment below if you need to make your own sample_annotation.txt file. You only need to execute this code once.
# Human brain microvascular endothelial cells (HMVECs)
# Supernatant is fluid that remains after centrifugation
#sampleNames <- colnames(betas)                              
#sampleNames[grep("EBM.*", colnames(betas), perl=TRUE)] <- "Endothelial.Basal.Medium"                # HBMVEC + EBS   
#sampleNames[grep("EGM.*", colnames(betas), perl=TRUE)] <- "Endothelial.Growth.Medium"               # HBMVEC + EBS + GFs
#sampleNames[grep("EVpellet.*", colnames(betas), perl=TRUE)] <- "Extracellular.Vesicle"              # HBMVEC + EBS + EV pellet
#sampleNames[grep("EVsuper.*", colnames(betas), perl=TRUE)] <- "Extracellular.Vesicle.Supernatant"   # HBMVEC + EBS + EV pellet + surnatant
#sampleNames[grep("GBM8EBM.*", colnames(betas), perl=TRUE)] <- "GBM8.in.EBM"                         # HBMVEC + GBM8 + supernatant
#sampleNames[grep("GBM8NB.*", colnames(betas), perl=TRUE)] <- "GBM8.in.NB"                           # HBMVEC + GBM8  
#sampleNames[grep("UCMpellet.*", colnames(betas), perl=TRUE)] <- "Unconditioned.Medium"              # GBM8 + normal medium
#sampleNames[grep("UCMsuper.*", colnames(betas), perl=TRUE)] <- "Unconditioned.Medium.Supernatant"   # GBM8 + without growth factors
#cell.lane <- gsub("^.*_MGH1_", "", colnames(betas))
#cell.lane <- gsub(".Signal_B", "", cell.lane)
#temp <- cbind(colnames(betas), sampleNames)
#colnames(temp) <- c("Sample_ID", "Group")
#write.table(temp, "sample_annotation.txt", quote=FALSE, eol="\r\n", row.names=FALSE, col.names=TRUE, sep="\t")
#rm(temp, sampleNames)

#home <- "~/Downloads/Breakefield/"
#home <- "C:/Users/Lillian/OneDrive/Lab/Breakefield"
#dirs <- list(data=paste(home,"raw.data/",sep=""),
#             annotation=paste(home,"raw.data/sample_annotation.txt",sep=""),
#             analysis=paste(home,"results/",sep=""),
#             reports=paste(home,"results/reports",sep=""))
#rnb.getOption("import.table.separator")
#rnb.options(filtering.sex.chromosomes.removal=TRUE,
#            import.table.separator="\t",  identifiers.column="Sample_ID")
#pheno <- read.table(dirs$annotation, sep="\t", header=TRUE)
#data.source <- RnBeadSet(pheno=pheno, betas=as.matrix(betas), probes = rownames(betas), useff=FALSE) # region.types=c("genes", "promoters"),
#rnb.initialize.reports(dirs$reports)
#rnb.set <- data.source
#rm(data.source, betas, pheno, home)
#save.image("initialized.RData")