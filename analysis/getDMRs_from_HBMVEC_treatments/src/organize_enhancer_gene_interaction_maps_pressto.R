#This Script organizes the enhancer gene mapping files from
# 'http://enhancer.binf.ku.dk/presets/' so that they can be used as
#custom annotations in the RNBBEADS based DMR analysis

#library(RnBeads)
library(biomaRt)
library(org.Hs.eg.db)
library(annotate)
#source(file = './processPRESSTOmapping.R') #contains My.Prepare.Enh.Prom.Assoc.Data()
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

#First import the PRESSTO cell-type active enhancer regions
#Create a list containing these objects
my.enhancer.profile.files <- list.files(path = './input', pattern = 'CL-')
my.profile.names <- gsub(pattern = '_differentially_expressed_enhancers.bed',
                         replacement = '',x = my.enhancer.profile.files)
my.profile.names <- gsub(pattern = '^...........',
                         replacement = '',x = my.profile.names)

my.enhancer.profiles = vector("list", length(my.profile.names))
names(my.enhancer.profiles) = paste(my.profile.names,".enh",sep = '')

for(ct in my.enhancer.profile.files){
    profile <- suppressWarnings(read.table( paste('./input/',ct ,sep = ''), sep = '\t', header = F, skip = 1))
    profile <- profile[,1:3]
    colnames(profile) = c("chromosome", "start", "end")
    my.enhancer.profiles[[which(my.enhancer.profile.files==ct)]] = profile
}

#The variable my.enhancer.profiles contains a list of several cell-type-active
#PRESSTO enhancer regions.
#These enhnacers need to be mapped to gene symbols or ensembl IDs

################################################################################
#First Define a functions to get the enhancer coordinate out of the mapping file:
#This function prepares the data for conversion to a GRANGES object
coord.convert.a <- function(factorColumn){
    coord.convert.b <- function(coordinate){
        c <- coordinate
        c <- unlist(strsplit(c,split = ':|-'))
    }
    
    fac <- as.character(factorColumn)
    out.coords <- sapply(X = fac,FUN = coord.convert.b)
    out.coords <- matrix(out.coords,ncol = 3, byrow = T)
    ch <- as.character(out.coords[,1])
    st <- as.numeric(out.coords[,2])
    en <- as.numeric(out.coords[,3])
    out.coords <- cbind.data.frame(ch,st,en)
    colnames(out.coords) <- c("chromosome", "start", "end")
    return(out.coords)

}
##################################################################
#Import the PRESSTO enhancer promoter mappings
#Select out the enhancer coordinates and mapped interacting promoter (symbol)
enh.PresstoProm.assoc <- My.Prepare.Enh.Prom.Assoc.Data() #Function from processPRESSTOmapping.R
enh.PresstoProm.assoc <- cbind.data.frame(coord.convert.a(enh.PresstoProm.assoc$Enh.coord),
                                "interaction_symbol" =  enh.PresstoProm.assoc$symbol)
GR.enh.prom.assoc <- makeGRangesFromDataFrame(enh.PresstoProm.assoc,keep.extra.columns = T)

#Convert the global enhancer map and cell-type specific maps to GRANGES objects
#Then use the mapping to assign gene symbols to the ct-active enhancers
#produce a list of data frames containing each cell-type assignment 

my.mapped.ehnancer.annots <- vector("list", length(names(my.enhancer.profiles)))
names(my.mapped.ehnancer.annots) <- names(my.enhancer.profiles)
asDF = T
for(ct.enh in names(my.enhancer.profiles)){
    GR.ct.enhancers = makeGRangesFromDataFrame(my.enhancer.profiles[[ct.enh]])
    assoc.enh.prom = subsetByOverlaps(GR.enh.prom.assoc,GR.ct.enhancers)
    if(asDF){
        df <- as.data.frame(assoc.enh.prom)
        colnames(df)[1] = "chromosome"
        my.mapped.ehnancer.annots[[ct.enh]] <- df
    } else {
    my.mapped.ehnancer.annots[[ct.enh]] <- assoc.enh.prom
    }
}

#We now have a list of enhancer regions that should be intersected with
#DRMs from our experiments.
