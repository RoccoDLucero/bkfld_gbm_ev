library(RnBeads)
library(RnBeads.hg19)

#This script should take a set of bed files (e.g TFBSs, or Enhancers, or their intersection)
# and convert that set into RnBeads annotation regions for use in 
# downstream analysis with the rest of the pipeline e.g. as is done for
# promoters, and genebodies


#Group D - EBM-2 was added with the following growth supplements:
#    human Epidermal Growth Factor (hEGF), 0.1%;
#    Vascular Endothelial Growth Factor (VEGF), 0.1%;
#    R3-Insulin-like Growth Factor-1 (R3-IGF-1), 0.1%;
#    Ascorbic Acid, 0.1%;
#    Hydrocortisone, 0.04%;
#    human Fibroblast Growth Factor-Beta (hFGF-β), 0.04%;
#    Gentamicin/Amphotericin-B (GA), 0.1%; P/S 1% grep('Endo',sample_info$name)

#Use bedtools to create files that contain E2f1 sites in open chromatin
#and E2f1 sites in closed chromatin in HUVECS or MesodermalCells.
#Get each of the files as an annotation set that can be compared
#to the DMRs of GrowthFactor Treated vs. Untreated HBMVECs

###################################################################################
###################################################################################
E2F1.open = rnb.set.annotation(type = "E2f1.controls.endo" ,regions = './E2F1.open.endothel.bed',assembly = "hg19")

#Hence you can use the annotation command in order to annotate obtained methylation or differential methylation values:
    aa <- annotation(rnb.set, type="promoters") annotated.dmrs <- data.frame(aa, dmrs) head(annotated.dmrs)

# Retrieve the chromHMM state segmentation from UCSC
    library(rtracklayer) mySession <- browserSession()
    genome(mySession) <- "hg19"
    tab.chromHMM.h1 <- getTable(ucscTableQuery(mySession, track="wgEncodeBroadHmm", table="wgEncodeBroadHmmH1hescHMM"))
# Filter for enhancer states
    tab.enhancers <- tab.chromHMM.h1[grep("Enhancer", tab.chromHMM.h1 $name), ]
# Select the interesting parts of the table and rename columns
    tab.enhancers <- tab.enhancers[, c("chrom", "chromStart", "chromEnd", "name")]
    colnames(tab.enhancers) <- c("chromosome", "start", "end", "name")
# Create RnBeads annotation by providing a data.frame rnb.set.annotation("enhancersH1EscChromHMM", tab.enhancers,assembly="hg19")
# Set the options to include the enhancer annotation
    rnb.options(region.types=c(rnb.getOption("region.types"),"enhancersH1EscChromHMM"))
    # Parse the input again, this time with the enhancer annotation added
    rnb.set.enh <rnb.execute.import(data.source=data.source, data.type="idat.dir") rnb.set.enh
# Annotation and methylation levels of enhancer regions in this object
    annot.enh <- annotation(rnb.set.enh, "enhancersH1EscChromHMM")
    head(annot.enh)
    meth.enh <- meth(rnb.set.enh, "enhancersH1EscChromHMM")
    head(meth.enh)
#Note that the included genomic regions remain available to RnBeads
#in the current R session only.
#If you later want to reuse custom annotation data, use the rnb.save.annotation() and rnb.load.annotation() functions:
    annotation.filename <- file.path(analysis.dir, "annotation_hg19_enhancersH1EscChromHMM.RData")
# Save the enhancer annotation to disk
    rnb.save.annotation(annotation.filename, "enhancersH1EscChromHMM", assembly="hg19")
    # Load the enhancer annotation as a duplicate
    rnb.load.annotation(annotation.filename, "enhancersH1EscChromHMM_duplicate")
    # Check that the annotation has been successfully loaded
    rnb.region.types()
