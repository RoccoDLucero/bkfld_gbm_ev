#"CheckTF_EnrichmentAtDMRs.R"
#Rocco Lucero April 21 2016
#This script is meant to be a temporary solution, while the
#linux ouput files are ugly and alternate lines of TFBS filenames,
# and intersect counts. This can only handle 2 files at a time to compute
# monte carlo simulated p-values that are not corrected
# for multiple hypothesis testing. 

my.comp.TF.enr = function(inputDir = "./", nameElems = c('*count*','*.bed'),
                          DMRcount = 400, nonDMRcount = 485177, sim.p.val = T)
    {
    temp = list.files(path = inputDir)
    temp1 = unlist(temp)
    for(elem in nameElems){
        temp1 = temp1[grep(elem,temp1)]
    }
    myfiles = lapply(temp1, read.delim,header = F,stringsAsFactors = F)
    names(myfiles) = temp1


    print(names(myfiles))
    my.intersects = cbind(myfiles[[1]],myfiles[[2]])
    TF.list = my.intersects[seq(from = 1,nrow(my.intersects),by = 2),1]
    TF.list = gsub('../../data/raw/tfchip/','',TF.list)
    TF.list = gsub('.sorted.bed','',TF.list)
    my.intersects = my.intersects[seq(from = 2,nrow(my.intersects),by = 2),]
    colnames(my.intersects) = c('DMR','non.DMR')
    my.intersects = cbind(TF.list, my.intersects)
    
    my.intersects$DMR = as.numeric(my.intersects$DMR)
    my.intersects$non.DMR = as.numeric(my.intersects$non.DMR)

    inputVals = c(DMRcount,nonDMRcount)

    chisq.p.vals = c()
    enr.ratios = c()
    for(row in 1:nrow(my.intersects)){
        tab = as.matrix(rbind(my.intersects[row,c('DMR','non.DMR')],inputVals))
        tmp = (chisq.test(tab,simulate.p.value = sim.p.val))
        chisq.p.vals = c(chisq.p.vals,tmp$p.value)
    
        enr.tmp = (tab[1,1]/tab[2,1])/(tab[1,2]/tab[2,2])
        enr.ratios =  c(enr.ratios,enr.tmp)
    }
    my.intersects =  cbind(my.intersects,chisq.p.vals,enr.ratios)
    colnames(my.intersects) = c('TF','DMRcount','non.DMRcount','chisq.p.vals','enr.ratios')
    my.intersects
}
