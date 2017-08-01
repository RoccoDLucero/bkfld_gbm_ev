#Cataloging genes and pathways that impact exRNA communication •
#Cross collaborations to further investigate the genes and pathways that 
#impact exRNA communication.
#The focus is on exRNAs that change the epigenetic state of cells and thereby 
#their transcriptome. It involves a collaborative interaction between Aleks, 
#Anna K, Louise Laurent and myself. Anna has already done her part in sequencing 
#exRNA fractions from glioblastoma cells and has sent that data to Aleks. 
#We have also done the methylome of the endothelial cells. 
#All that is left is the RNA sequencing of the endothelial cells with 
#and without exRNA exposure and the bioinformatics analysis

#I am wondering if any of Rocco's data yet suggest a pathway for cell permeability
#or another parameter that would increase penetrability of the endothelial cells.
#Can you check on that? or should I contact Rocco directly


library(reshape2)
library(ggplot2)


head(mirTARBase[,1:3])

Atlas <- readRDS("../../../exrna_atlas/atlas.all.counts.and.biogps.meta")
Atlas <- Atlas$miRNA

#Get GBM and control samples for the ExRNA Atlas miRNA
#This will serve as our test group and control data
#We aim to test whether certain mRNA are differentially targeted by extracellular miRNA

my.get.miRNA.Levels.of.target <- function(target.mRNA, test.EV.data,condition,control.EV.data, mirTARBase){
    tar.M <- target.mRNA
    tar.M <- paste("^",tar.M,"$",sep = "")
    tst <- test.EV.data
    ctrl <- control.EV.data
    
    miRNA.rows <- grep(tar.M,mirTARBase$Target.Gene,ignore.case = T) #For given target pull all miRNA from MirTARBASE
    miRNA <- mirTARBase[miRNA.rows,"miRNA"]
    
    #apply(AKRIC_GBM_mirnaRPM[,which(colnames(AKRIC_GBM_mirnaRPM) %in% zo1mirna)],2,median)
    #apply(healthy_control_mirnaRPM[,which(colnames(healthy_control_mirnaRPM) %in% zo1mirna)],2,median)
    
    a1 <-(tst[,which(colnames(tst) %in% miRNA)])
    a2 <-(ctrl[,which(colnames(ctrl) %in% miRNA)])
    a1$EXP <- condition
    a2$EXP <- "HEALTHY"
    a <- rbind(a1,a2)
    
    #my.is.diff <- ks.test(a1[,1],a2[,1],alternative = "less")
    
    df_long <- melt(a[,1:ncol(a)], id.vars= 'EXP')
    p2 <- ggplot(df_long, aes(x=factor(EXP),y=value,fill=factor(EXP)))+
        geom_boxplot() + labs(title= paste("miRNA Targeting", target.mRNA) )+
        facet_wrap(~variable,scales = "free")
    
    print(p2)
    
}

my.check.targ <- function(targ, mirTARBase){
    ptt <- paste("^",targ,"$",sep = "")
    targeted.mirna <- grep(ptt, mirTARBase$Target.Gene,ignore.case = T)
    print(mirTARBase[targeted.mirna, c("miRNA","Target.Gene","Experiments")])
    
}

################################################################################
#For some RNA of interest check to see if it is targeted by miRNA
new.targ <- "ccl20"
my.check.targ(targ = new.targ,mirTARBase = mirTARBase)


################################################################################
#Look for differential miRNA expression: GBM-cell exrna vs healthy conrol exrna
AKRIC_GBM_mirnaRPM <- Atlas[grep("AKRIC.*-10-17",Atlas$Study.x),]
healthy_control_AllBFmirnaRPM <- Atlas[grep("Healthy",Atlas$condition),]
healthy_control_CSFmirnaRPM <- Atlas[which(Atlas$condition == "Healthy" | Atlas$biofluid_name == "Cerebrospinal fluid"),]
healthy_control_NONCSFmirnaRPM <- Atlas[which(Atlas$condition == "Healthy" | Atlas$biofluid_name != "Cerebrospinal fluid"),]

levels(Atlas$biofluid_name)

tested.trgs <- NULL
trgs.vec <- c("tjp1","cdh5","tnfaip1","tnfaip3","timp1","cxcl1","cxcl3","cx3cl1",
              "cxcl10","ccl20","vegfa","icam1")

tst.dat <- AKRIC_GBM_mirnaRPM
ctrl.dat <- healthy_control_CSFmirnaRPM
cnd <- "GBM8 cell culture"
for(trg in trgs.vec[1:3]){
    if(trg %in% tested.trgs){
        print(paste(trg,"already tested."))
        next
    }else{
        print(paste("Testing:",trg))
        my.get.miRNA.Levels.of.target(target.mRNA = trg,
                              test.EV.data =  tst.dat,
                              condition = cnd,
                              control.EV.data = ctrl.dat,
                              mirTARBase)
        tested.trgs <- c(tested.trgs,trg)
        
    }
}
    
################################################################################
#Look for differential miRNA expression: aSAH exrna vs healthy conrol exrna
levels(Atlas$Study.x)
tst.dat <- Atlas[grep("-aSAH_",Atlas$Study.x),]
ctrl.dat <- healthy_control_AllBFmirnaRPM
cnd <- "aSAH"

tested.trgs <- NULL

for(trg in trgs.vec[1:3]){
    if(trg %in% tested.trgs){
        print(paste(trg,"already tested."))
        next
    }else{
        print(paste("Testing:",trg))
        my.get.miRNA.Levels.of.target(target.mRNA = trg,
                                      test.EV.data =  tst.dat,
                                      condition = cnd,
                                      control.EV.data = ctrl.dat,
                                      mirTARBase)
        tested.trgs <- c(tested.trgs,trg)
        
    }
}

################################################################################

test.dat <- healthy_control_CSFmirnaRPM
ctrl.dat <- healthy_control_NONCSFmirnaRPM    
cnd <- "Healthy CSF"
tested.trgs <- NULL

for(trg in trgs.vec[1:3]){
    if(trg %in% tested.trgs){
        print(paste(trg,"already tested."))
        next
    }else{
        print(paste("Testing:",trg))
        my.get.miRNA.Levels.of.target(target.mRNA = trg,
                                      test.EV.data =  tst.dat,
                                      condition = cnd,
                                      control.EV.data = ctrl.dat,
                                      mirTARBase)
        tested.trgs <- c(tested.trgs,trg)
        
    }
}

################################################################################
apply(AKRIC_GBM_mirnaRPM[,which(colnames(AKRIC_GBM_mirnaRPM) %in% zo1mirna)],2,median)
apply(healthy_control_mirnaRPM[,which(colnames(healthy_control_mirnaRPM) %in% zo1mirna)],2,median)


apply(AKRIC_GBM_mirnaRPM[,which(colnames(AKRIC_GBM_mirnaRPM) %in% VECadherinmirna)],2,median)
apply(healthy_control_mirnaRPM[,which(colnames(healthy_control_mirnaRPM) %in% VECadherinmirna)],2,median)


a1 <-(AKRIC_GBM_mirnaRPM[,which(colnames(AKRIC_GBM_mirnaRPM) %in% zo1mirna)])
a2 <-(healthy_control_mirnaRPM[,which(colnames(healthy_control_mirnaRPM) %in% zo1mirna)])
b1 <-(AKRIC_GBM_mirnaRPM[,which(colnames(AKRIC_GBM_mirnaRPM) %in% VECadherinmirna)])
b2 <-(healthy_control_mirnaRPM[,which(colnames(healthy_control_mirnaRPM) %in% VECadherinmirna)])

a1$EXP <- "GBM"
a2$EXP <- "HEALTHY"
b1$EXP <- "GBM"
b2$EXP <- "HEALTHY"


a <- rbind(a1,a2)
b <- rbind(b1,b2)

df_long <- melt(a[,1:ncol(a)], id.vars= 'EXP')
    
#then plot
p2 <- ggplot(df_long, aes(x=factor(EXP),y=value,fill=factor(EXP)))+
    geom_boxplot() + labs(title="miRNA Targeting ZO1") +facet_wrap(~variable,scales = "free")

print(p2)




df_long <- melt(b[,c(1:ncol(b))], id.vars= 'EXP')

#then plot
p2 <- ggplot(df_long, aes(x=factor(EXP),y=value,fill=factor(EXP)))+
    geom_boxplot() +
    labs(title="Exosomal miRNA Targeting VE Cadherin", y="Reads per Million")+
    facet_wrap(~variable, scales ="free_y")

print(p2)


dim(b)
