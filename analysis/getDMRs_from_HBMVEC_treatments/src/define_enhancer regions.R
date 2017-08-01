library(rtracklayer)
library(vioplot)

if(F){
mySession <- browserSession()
genome(mySession) <- "hg19"
tab.chromHMM.h1 <- getTable(ucscTableQuery(mySession,
                                           track="wgEncodeBroadHmm", table="wgEncodeBroadHmmH1hescHMM"))
# Filter for enhancer states
tab.enhancers <- tab.chromHMM.h1[grep("Enhancer", tab.chromHMM.h1$name), ]
class(tab.enhancers)
head(tab.enhancers)

# Select the interesting parts of the table and rename columns
tab.enhancers <- tab.enhancers[, c("chrom", "chromStart", "chromEnd", "name")]
colnames(tab.enhancers) <- c("chromosome", "start", "end", "name")
# Create RnBeads annotation by providing a data.frame
rnb.set.annotation("enhancersH1EscChromHMM", tab.enhancers, assembly="hg19")
# Set the options to include the enhancer annotation
rnb.options(region.types=c(rnb.getOption("region.types"),"enhancersH1EscChromHMM"))

# Parse the input again, this time with the enhancer annotation added
rnb.set.enh <-
    rnb.execute.import(data.source=data.source, data.type="idat.dir")
rnb.set.enh
# Annotation and methylation levels of enhancer regions in this object
annot.enh <- annotation(rnb.set.enh, "enhancersH1EscChromHMM")
head(annot.enh)
meth.enh <- meth(rnb.set.enh, "enhancersH1EscChromHMM")


head(meth.enh)
hist(my.rnb.set@meth.regions$bv.endoEnhancersPRESSTO[,1:3],col = rgb(.1,.1,.1,.5), ylim = c(0,1200))
hist(my.rnb.set@meth.regions$bv.endoEnhancersPRESSTO[,4:6],add = T, col = rgb(1,.8,.8,.5))
hist(my.rnb.set@meth.regions$bv.endoEnhancersPRESSTO[,7:9],add = T, col = rgb(1,0,1,.5))

hist(my.rnb.set@meth.regions$genes)

hist(my.rnb.set@meth.regions$promoters[,1:3],col = rgb(.1,.1,.7,.5), ylim = c(0,20000))
hist(my.rnb.set@meth.regions$promoters[,7:9],add = T, col = rgb(1,.8,.2,.3))

hist(my.rnb.set@meth.regions$astroEnhancersPRESSTO[,1:3],col = rgb(.1,.1,.7,.5), ylim = c(0,700))
hist(my.rnb.set@meth.regions$astroEnhancersPRESSTO[,7:9],add = T, col = rgb(1,.8,.2,.3))

a = my.rnb.set@meth.regions$astroEnhancersPRESSTO[,1:3] - my.rnb.set@meth.regions$astroEnhancersPRESSTO[,7:9]
b = apply(X = a,MARGIN = 1,FUN = min)
min(b,na.rm = T)
b[b>.05]
class(my.rnb.set)
}


#############################
############################
library(ggplot2)
for(choop in 1:7){
    z <- names(my.dmr.list[[1]])[choop]
    z
    a = my.dmr.list$EBM.EGM[[z]][order(my.dmr.list$EBM.EGM[[z]]$combinedRank),]$mean.mean.diff[1:100]
    b = my.dmr.list$EBM.GBM8EVpel[[z]][order(my.dmr.list$EBM.GBM8EVpel[[z]]$combinedRank),]$mean.mean.diff[1:100]
    c = my.dmr.list$EBM.GBM8sup[[z]][order(my.dmr.list$EBM.GBM8sup[[z]]$combinedRank),]$mean.mean.diff[1:100]
    d = my.dmr.list$EBM.UCMpel[[z]][order(my.dmr.list$EBM.UCMpel[[z]]$combinedRank),]$mean.mean.diff[1:100]
    e = my.dmr.list$EBM.UCMsup[[z]][order(my.dmr.list$EBM.UCMsup[[z]]$combinedRank),]$mean.mean.diff[1:100]
    cols <- list(rgb(.8,.1,.1,.3),rgb(.1,.5,.1,.3),
                 rgb(.1,.1,.5,.3),rgb(.5,.5,.1,.3),rgb(.1,.5,.5,.3))
    cols <- unlist(cols)
    #boxplot(a[abs(a)>.02],b[abs(b)>.02],c[abs(c)>.02],d[abs(d)>.02],e[abs(e)>.02],
    #        col = cols, ylim = c(-.3,.3),main =z)
    pdf(file = paste('./output/junkplots',choop,'.pdf',sep=''))
    boxplot(-a,-b,-c,-d,-e,col = cols, ylim = c(-.2,.2),main=z )
    dev.off()
}

pdf(file = './output/junkplots.pdf')
hist(-a[abs(a)>.02], col = cols[[1]], ylim = c(0,800), breaks = seq(-.4,.4,.01))
hist(-b[abs(b)>.02],add = T,col= cols[[2]])
hist(-c[abs(c)>.02],add = T, col = cols[[3]])
hist(-d[abs(d)>.02],add = T,col= cols[[4]])
hist(-e[abs(e)>.02],add = T, col = cols[[5]])
legend(legend = names(my.dmr.list)[c(1,2,3,6,7)], "topright", c(a,b,c,d,e), col= cols , lwd=10)



hist(c[abs(c)>.02],add = T, col = cols[[3]])
hist(d[abs(d)>.02],add = T,col= cols[[4]])
hist(e[abs(e)>.02],add = T, col = cols[[5]])

dev.off()
length(a[abs(a)>.02])
length(b[abs(b)>.02])
length(c[abs(c)>.02])
length(d[abs(d)>.02])
length(e[abs(e)>.02])


#vioplot(a,b,c,ylim = c(-1,1))

#mplot <- ggplot(cbind(a,b), aes(x = ) + geom_violin()
