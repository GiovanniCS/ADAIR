library(here)
library(lilikoi)
setwd(paste(here(),"/Data",sep=""))

samples.to.keep = 6
genes.to.keep = 50 #24554
filter.out.genes.not.in.pathway = TRUE
filter.out.metabolites.not.in.pathway = FALSE
cases = as.integer(samples.to.keep/2)
controls = as.integer(samples.to.keep/2)

dt = lilikoi.Loaddata(file=system.file("extdata", "plasma_breast_cancer.csv", package = "lilikoi"))

tot.sample = dim(dt[["Metadata"]])[1]
metab.matrix = dt[["Metadata"]]
metab.matrix = metab.matrix[c(1:cases,tot.sample:(tot.sample-controls+1)),]
dataSet = dt$dataSet

#metabolite_names = read.csv("Tosca/metabolite_names.csv", sep=";")
#dataSet[["cmpd"]] = metabolite_names$cohort1

convertResults = lilikoi.MetaTOpathway('name')

metab_pathway_mapping = convertResults[["table"]]
if(filter.out.metabolites.not.in.pathway){
  metab_pathway_mapping = metab_pathway_mapping[!is.na(metab_pathway_mapping$pathway),]
}
metab_pathway_mapping = as.data.frame(metab_pathway_mapping)
metab_pathway_mapping$pathway = as.character(metab_pathway_mapping$pathway)
for(i in 1:nrow(metab_pathway_mapping)){
  pathways = metab_pathway_mapping$pathway[[i]]
  if(!is.na(pathways)){
    type = strsplit(pathways,split = " ")[[1]]
    type = trimws(type[length(type)])
    pathways = strsplit(pathways,split = "and")[[1]]
    pathways = trimws(unlist(strsplit(pathways,",")))
    temp = strsplit(pathways[length(pathways)]," ")[[1]]
    temp = temp[-length(temp)]
    pathways[length(pathways)] = paste(temp,collapse=" ")
    pathways = sapply(pathways, function(p){paste(p,type,sep=" ")})
    pathways = trimws(pathways)
  }
  else{
    pathways = "None"
    names(a) = "None"
  }
  metab_pathway_mapping$pathway[i] = list(pathways)
}

pathways_temp = unique(unlist(metab_pathway_mapping$pathway))
metabolites.pathways = list()
for(i in 1:length(pathways_temp)){
  metabolites.pathways[[pathways_temp[i]]] = list()
}

for(i in 1:nrow(metab_pathway_mapping)){
  pws = metab_pathway_mapping[i,]$pathway[[1]]
  for(j in 1:length(pws)){
    p = pws[j]
    l = metabolites.pathways[[p]]
    l[[length(l)+1]] = metab_pathway_mapping[i,]$Query
    metabolites.pathways[[p]] = l
  }
}

metabolites.pathways
metab.matrix = t(metab.matrix)
metab.matrix

#SNPs mapping to genes
library(vcfR)
library(biomaRt)
library(KEGGREST)
library(tidyverse)

mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl", mart)
attributes <- c("entrezgene_id","hgnc_symbol","start_position","end_position","chromosome_name")
filters <- c("chromosome_name")
values <- list(chromosome=c(1:23))
genes.names <- getBM(attributes=attributes, filters=filters, values=values, mart=mart)
genes.names <- genes.names[!is.na(genes.names$entrezgene_id),]

hsa_pathways <- keggList("pathway", "hsa") %>% 
    tibble(pathway = names(.), description = .)
hsa_path_eg  <- keggLink("pathway", "hsa") %>% 
    tibble(pathway = ., egid = sub("hsa:", "", names(.))) %>% 
    group_by(pathway) %>% 
    summarize(genes = list(egid))
hsa <- left_join(hsa_path_eg, hsa_pathways)
snps.pathways = as.list(hsa$genes)
names(snps.pathways) = hsa$pathway

if(filter.out.genes.not.in.pathway){
  genes.names <- genes.names[genes.names$entrezgene_id %in% unlist(hsa_path_eg$genes),]
}else {
  genes.without.pathway <- genes.names[!(genes.names$entrezgene_id %in% unlist(hsa_path_eg$genes)),]
  snps.pathways$`path::none` = as.character(genes.without.pathway$entrezgene_id)
}

genes.names <- genes.names[1:genes.to.keep,]
genes.names

system("wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz")
system("mv ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz All.vcf.gz")
system(paste("vcftools --gzvcf All.vcf.gz",
  " --from-bp 1 --to-bp 1 --chr 1 --recode --out SampleIDs",
  sep="")
)
vcf <- read.vcfR("SampleIDs.recode.vcf", verbose = FALSE )
samples = colnames(vcf@gt)[-1]
samples = samples[1:samples.to.keep]
write.table(samples,file="names.txt",sep="\n",quote = FALSE,col.names=FALSE,row.names=FALSE)
system("vcftools --gzvcf All.vcf.gz --keep names.txt --recode --out filtered_vcf")

cat("Sample ids retrived and vcf filtered\n")

SNPs.to.genes <- matrix(0,nrow = nrow(genes.names),ncol = length(samples))
rownames(SNPs.to.genes) <- genes.names$entrezgene_id
colnames(SNPs.to.genes) <- samples
for(i in seq_len(nrow(genes.names))){
  cat(paste("Gene ",genes.names$ensembl_gene_id[i],"\n",sep=""))
  start <- genes.names$start_position[i]
  end <- genes.names$end_position[i]
  chromosome <- genes.names$chromosome_name[i]
  system(paste("vcftools --gzvcf filtered_vcf.recode.vcf",
    " --from-bp ",start," --to-bp ",end," --chr ",chromosome,
    " --recode --out current_gene",sep="")
  )
  vcf <- read.vcfR( "current_gene.recode.vcf", verbose = FALSE )
  variants = dim(vcf)[1]
  if(variants > 0){
    variants <- vcf@gt[,-1]
    is.single.vector = is.null(dim(variants))
    if(is.single.vector){
      for(j in 1:length(variants)){
        snp = variants[j]
        if(!is.na(snp) && snp!= "0/0" && snp!= "0|0"){
          SNPs.to.genes[i,j] = SNPs.to.genes[i,j] + 1
        }
      }
    }
    else{
      rows = nrow(variants)
      cols = ncol(variants)
      for(j in 1:rows){
        for(k in 1:cols){
          snp = variants[j,k]
          if(!is.na(snp) && snp!= "0/0" && snp!= "0|0"){
            SNPs.to.genes[i,k] = SNPs.to.genes[i,k] + 1
          }
        }
      }
    }
  }
}
rownames(SNPs.to.genes) = genes.names$entrezgene_id
SNPs.to.genes = SNPs.to.genes[rowSums(SNPs.to.genes) > 0,]
colnames(SNPs.to.genes) = colnames(metab.matrix)


#netDx integration
library(SummarizedExperiment)
library(MultiAssayExperiment)
suppressWarnings(suppressMessages(require(netDx)))
class = metab.matrix[1,]
metabolites = metab.matrix[-1,]
se <- SummarizedExperiment(assays=list(metabolites = metabolites,snps=SNPs.to.genes))
exp <- list(metabolites = se)
patient.data <- data.frame(ID=colnames(metabolites),STATUS=class)
rownames(patient.data) <- patient.data$ID
multiomics <- MultiAssayExperiment(experiments=exp,colData = patient.data)
# create grouping rules
groupList <- list()
groupList[["metabolites"]] <- metabolites.pathways
# clinical data is not grouped; each variable is its own feature
groupList[["snps"]] <- snps.pathways

makeNets <- function(dataList, groupList, netDir,...) {
    netList <- c() # initialize before is.null() check
    # correlation-based similarity for mRNA, RPPA and methylation data
    # (Pearson correlation)
    for (nm in names(groupList)) {
       # NOTE: the check for is.null() is important!
        if (!is.null(groupList[[nm]])) {
        netList <- makePSN_NamedMatrix(dataList[[nm]],
                                      rownames(dataList[[nm]]),
                                      groupList[[nm]],
                                      netDir,
                                      verbose=FALSE,
                                      writeProfiles=TRUE,...) 
        }
    }
    return(netList)
}

set.seed(42) # make results reproducible

# location for intermediate work
# set keepAllData to TRUE to not delete at the end of the 
# predictor run.
# This can be useful for debugging.
outDir <- paste(tempdir(),"pred_output",sep=getFileSep()) # use absolute path
system(paste("mkdir -p ",outDir,sep=""))
numSplits <- 2L
out <- suppressMessages(
   buildPredictor(dataList=multiomics,groupList=groupList,
      makeNetFunc=makeNets,
      outDir=outDir, ## netDx requires absolute path
      numSplits=numSplits, featScoreMax=2L, featSelCutoff=1L,
       numCores=1L)
)
