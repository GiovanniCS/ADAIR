library(MultiAssayExperiment)

# load metabolites names we have received from Rotterdam
script.directory <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script.directory)
load("../Data/Tosca/metabolite_names.Rdata")

# metabolites dummy dataset has normally distributed features.
# mu and sigma bounds are extracted from dummy data
# 245 (and 231) features comes from the metabolites names
mu <- runif(245, 0.0027, 1.2)
sigma <- runif(245, 2.6e-05, 0.09)
# original dummy dataset sample size
sample.size <- 2001

RSI.metabolites <- mapply(function(x,y){rnorm(x,y,n=sample.size)},x=mu,y=sigma)
RSI.metabolites <- as.data.frame(RSI.metabolites)
colnames(RSI.metabolites) = metabolite_names[metabolite_names[,1] != "",1]

mu <- runif(231, 0.0027, 1.2)
sigma <- runif(231, 2.6e-05, 0.09)

RSII.metabolites <- mapply(function(x,y){rnorm(x,y,n=sample.size)},x=mu,y=sigma)
RSII.metabolites <- as.data.frame(RSII.metabolites)
colnames(RSII.metabolites) = metabolite_names[metabolite_names[,2] != "",2]

# get common metabolites between cohorts
common.metab <- intersect(colnames(RSI.metabolites),colnames(RSII.metabolites))
metab.to.keep.RSI <- common.metab %in% colnames(RSI.metabolites)
RSI.metabolites <- RSI.metabolites[,metab.to.keep.RSI]

metab.to.keep.RSII <- common.metab %in% colnames(RSII.metabolites)
RSII.metabolites <- RSI.metabolites[,metab.to.keep.RSII]

# order features and merge cohort samples
RSI.metabolites <- RSI.metabolites[,sort(colnames(RSI.metabolites))]
RSII.metabolites <- RSII.metabolites[,sort(colnames(RSII.metabolites))]
Rotterdam.metabolites <- rbind(RSI.metabolites,RSII.metabolites)

# add patient id
rownames(Rotterdam.metabolites) <- sample(as.character(c(10000:20000)),4002)

# adapt to netDX input requirements
Rotterdam.metabolites <- t(Rotterdam.metabolites)

exp = list(metabolites = Rotterdam.metabolites)
#IMPORTANT: dummy dataset don't have case/control label
patient.data <- data.frame(ID=colnames(Rotterdam.metabolites),STATUS=sample(c("case","control"),4002,replace = TRUE))
rownames(patient.data) = patient.data$ID
Rotterdam.metabolites <- MultiAssayExperiment(experiments=exp,colData = patient.data)

# save dataset (e.g. to import into a docker container)
save(Rotterdam.metabolites,file="Rotterdam.metabolites.Rdata")


# INSIDE NETDX ENVIRONMENT ...
groupList <- list()
metab.groups <- list(rownames(Rotterdam.metabolites))
names(metab.groups[[1]]) <- rownames(Rotterdam.metabolites)
groupList[["Metabolites"]] <- metab.groups[[1]]

outDir <- paste(tempdir(),randAlphanumString(), "pred_output",sep=getFileSep())
system(paste("mkdir -p ",outDir,sep=""))
netList <- makePSN_NamedMatrix(Rotterdam.metabolites,rownames(Rotterdam.metabolites),groupList[["Metabolites"]],outDir,verbose=FALSE,writeProfiles=TRUE)
# Pearson similarity chosen - enforcing min. 5 patients per net.
# Error in { : task 1 failed - "'match' requires vector arguments"

