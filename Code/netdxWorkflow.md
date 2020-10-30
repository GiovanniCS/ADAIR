# Workflow for setup and execution of netDx R package

````bash
# Open a bash terminal and follow the instructions below.
# Check if Docker is installed in your system, otherwise go to https://docs.docker.com/get-docker/
docker -v
# Pull this Docker image with netDx (and all its dependencies) ready to be executed
# This image weight about 2GB, it may take some time to download
docker pull giovannics/netdx:latest
# Create a container from the image and enter in the container
docker run -it giovannics/netdx:latest /bin/bash
# You have now access to a bash terminal running inside the container.
# From here you can start an R session
R
````
The following script is an example of how netDx can be used for build patient similarity networks (PSNs) from transcriptomics and clinical data.   
Reference: https://www.embopress.org/doi/full/10.15252/msb.20188497  
Once the R session has started you can follow the script copying the commands in sequential order.
````R
suppressWarnings(suppressMessages(require(netDx)))
suppressWarnings(suppressMessages(library(curatedTCGAData)))

# fetch data remotely
brca <- suppressMessages(curatedTCGAData("BRCA",c("mRNAArray"),FALSE))

# process input variables
staget <- sub("[abcd]","",sub("t","",colData(brca)$pathology_T_stage))
staget <- suppressWarnings(as.integer(staget))
colData(brca)$STAGE <- staget

pam50 <- colData(brca)$PAM50.mRNA
pam50[which(!pam50 %in% "Luminal A")] <- "notLumA"
pam50[which(pam50 %in% "Luminal A")] <- "LumA"
colData(brca)$pam_mod <- pam50

tmp <- colData(brca)$PAM50.mRNA
idx <- union(which(tmp %in% c("Normal-like","Luminal B","HER2-enriched")),
                    which(is.na(staget)))
pID <- colData(brca)$patientID
tokeep <- setdiff(pID, pID[idx])
brca <- brca[,tokeep,]

smp <- sampleMap(brca)
samps <- smp[which(smp$assay=="BRCA_mRNAArray-20160128"),]
# remove duplicate assays mapped to the same sample
notdup <- samps[which(!duplicated(samps$primary)),"colname"]
brca[[1]] <- suppressMessages(brca[[1]][,notdup])

# colData must have ID and STATUS columns
pID <- colData(brca)$patientID
colData(brca)$ID <- pID
colData(brca)$STATUS <- colData(brca)$pam_mod

# create grouping rules
groupList <- list()
# genes in mRNA data are grouped by pathways
pathList <- readPathways(fetchPathwayDefinitions("January",2018))
groupList[["BRCA_mRNAArray-20160128"]] <- pathList[1:3]
# clinical data is not grouped; each variable is its own feature
groupList[["clinical"]] <- list(
      age="patient.age_at_initial_pathologic_diagnosis",
       stage="STAGE"
)

# create function to tell netDx how to build features (PSN) from your data
makeNets <- function(dataList, groupList, netDir,...) {
    netList <- c() # initialize before is.null() check
    # make RNA nets (NOTE: the check for is.null() is important!)
    # (Pearson correlation)
    if (!is.null(groupList[["BRCA_mRNAArray-20160128"]])) { 
    netList <- makePSN_NamedMatrix(dataList[["BRCA_mRNAArray-20160128"]],
                rownames(dataList[["BRCA_mRNAArray-20160128"]]),
                groupList[["BRCA_mRNAArray-20160128"]],
                netDir,verbose=FALSE, 
                writeProfiles=TRUE,...) 
    }

    # make clinical nets (normalized difference)
    netList2 <- c()
    if (!is.null(groupList[["clinical"]])) {
    netList2 <- makePSN_NamedMatrix(dataList$clinical, 
        rownames(dataList$clinical),
        groupList[["clinical"]],netDir,
        simMetric="custom",customFunc=normDiff, # custom function
        writeProfiles=FALSE,
        sparsify=TRUE,verbose=TRUE,...)
    }
    netList <- c(unlist(netList),unlist(netList2))
    return(netList)
}

# run predictor 
set.seed(42) # make results reproducible
outDir <- paste(tempdir(),randAlphanumString(),
    "pred_output",sep=getFileSep())
# To see all messages, remove suppressMessages() and set logging="default".
# To keep all intermediate data, set keepAllData=TRUE
out <- buildPredictor(
      dataList=brca,groupList=groupList,
      makeNetFunc=makeNets,
      outDir=outDir, ## netDx requires absolute path
      numSplits=2L,featScoreMax=2L, featSelCutoff=1L,
      numCores=1L,logging="none",
      keepAllData=FALSE,debugMode=TRUE
   )

# collect results
numSplits <- 2
st <- unique(colData(brca)$STATUS)
acc <- c()         # accuracy
predList <- list() # prediction tables
featScores <- list() # feature scores per class
for (cur in unique(st)) featScores[[cur]] <- list()

for (k in 1:numSplits) { 
    pred <- out[[sprintf("Split%i",k)]][["predictions"]];
    # predictions table
    tmp <- pred[,c("ID","STATUS","TT_STATUS","PRED_CLASS",
                     sprintf("%s_SCORE",st))]
    predList[[k]] <- tmp 
    # accuracy
    acc <- c(acc, sum(tmp$PRED==tmp$STATUS)/nrow(tmp))
    # feature scores
    for (cur in unique(st)) {
       tmp <- out[[sprintf("Split%i",k)]][["featureScores"]][[cur]]
       colnames(tmp) <- c("PATHWAY_NAME","SCORE")
       featScores[[cur]][[sprintf("Split%i",k)]] <- tmp
    }
}
````