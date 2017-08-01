
################################################################################
#Set these variables to access Differential Methylation Data files
# in the subfolders created by RnBeads pipeline
#file pathe will be: 'pathroot/analysisFolder/datFolder/comps' 
#
pathRoot <- './RnBeads/analysis/'

analysisFolder <- c('noGBM8_noNegCtrl_w_astro_endo_enhancers/',
                   'reports_Batch_Corrected/',
                   'reports_pairwise3/')

datFolder <- 'differential_methylation_data/'

comps <- c('diffMethTable_site_cmp1.csv',
          'diffMethTable_site_cmp2.csv',
          'diffMethTable_site_cmp3.csv')

################################################################################
#This funcrion
DMR.HIST.GROB <- function(diffDat,xlab,ylab,xlim = c(-1,1),ylim = c(0,(4*10^5))){
    library(ggplot2)
    library(grid)
    library(gridExtra)
    my.diffMethTable <- read.csv(diffDat,header = T) 
    
    my.histo.grob <- qplot(my.diffMethTable$mean.diff,
                           geom = 'histogram',
                           ylim = ylim, xlim = xlim,
                           xlab = xlab)   
}
################################################################################
##Make ggplot GROB Histograms and produce a grid plot output
##
##
my.savePath <- './output/DRM_HistoGrid.pdf'

title <- 'Basal vs. GF-treated HBMVEC'
diffDat <- paste(pathRoot,analysisFolder[1],datFolder,comps[1],sep = '')
h1 <- DMR.HIST.GROB(diffDat,xlab = title)

title <- 'Basal vs. EV-treated HBMVEC'
diffDat <- paste(pathRoot,analysisFolder[1],datFolder,comps[2],sep = '')
h2 <- DMR.HIST.GROB(diffDat,xlab = title)

title <- 'BASAL vs. EV.sup treated HBMVEC'
diffDat <- paste(pathRoot,analysisFolder[1],datFolder,comps[3],sep = '')
h3 <- DMR.HIST.GROB(diffDat,xlab = title)

title <- 'HBMVEC vs. GBM8'
diffDat <- paste(pathRoot,analysisFolder[3],datFolder,comps[3],sep = '')
h4 <- DMR.HIST.GROB(diffDat,xlab = title)

my.grobOrder <- list(h1,h2,h3,h4)

grid.arrange(h1,h2,h3,h4)
g = arrangeGrob(grobs = my.grobOrder)
ggsave(my.savePath,g,device = 'pdf')

