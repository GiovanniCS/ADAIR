library(caret)
library(edgeR)
library(DaMiRseq)

load("/Users/giovanni/Desktop/ADAIR/Luca/op1_prepared_data.rda")
pd_iso_veh = data_raw_l[["infoc"]][["conds"]] == "VEH"
data_veh = data_raw_l
data_veh$counts = data_veh$counts[,pd_iso_veh]
data_veh$infoc = data_veh$infoc[pd_iso_veh,]

# Removing genes that have equal expression level in at least 5 out of 7 samples 
# (this includes alzo genes that are not expressed in at least 5 out of 7 samples)
data_veh[["counts"]] = t(data_veh[["counts"]])
nzv = preProcess(data_veh[["counts"]],method="nzv",uniqueCut = 30)
data_veh[["counts"]] = predict(nzv,data_veh[["counts"]])
data_veh[["counts"]] = t(data_veh[["counts"]])

# TEMP FIX FOR NaN ISSUE
# temp = data_veh[["counts"]]
# colnames(temp) = paste(colnames(temp),"bis",sep="")
# data_veh[["counts"]] = cbind(data_veh[["counts"]],temp)

# Normalization: Each sample is divided by total count and multiply by 10^6
data_veh[["norm"]] = cpm(data_veh[["counts"]])

#Select most variable genes
compute_cv <- function(x) sd(x) / mean(x)
cv <- apply(data_veh[["norm"]], 1, compute_cv)
cutoff <- 0.25
x <- data_veh[["norm"]][rank(cv) / length(cv) > 1 - cutoff, ]

#prepare for model training
x = t(x)
x = cbind(data.frame(x),y = factor(data_veh[["infoc"]][["comb_names"]]))

set.seed(17)
mtry=as.integer(sqrt(dim(x)[2]))
model_caret <- train(y~., 
                     data = x,
                     method = "ranger",
                     trControl = trainControl(method = "LOOCV"),
                     tuneGrid = data.frame(mtry=mtry,
                                           min.node.size = 1,
                                           splitrule="gini"),
                     #importance = "permutation"
                     importance = "impurity"
)
# What impurity and permutation mean?
# https://alexisperrier.com/datascience/2015/08/27/feature-importance-random-forests-gini-accuracy.html
# https://github.com/imbs-hl/ranger/issues/237
plot(varImp(model_caret),top=10)





# DaMiRseq workflow
data_veh = data_raw_l
data_veh$counts = data_veh$counts[,pd_iso_veh]
data_veh$infoc = data_veh$infoc[pd_iso_veh,]

x = data_veh[["counts"]]
x = lapply(x, function(x) as.integer(x))
x = data.frame(x)
y = data.frame(class = rep(data_veh[["infoc"]][["comb_names"]],2))
rownames(y) = colnames(x)
SE = DaMiR.makeSE(x, y)

data_norm <- DaMiR.normalization(SE, minCounts=10, fSample=0.7, th.cv=3, type="rlog")
set.seed(12345)
data_clean<-DaMiR.transpose(assay(data_norm))
df<-colData(data_norm)
# data_reduced <- DaMiR.FSelect(data_clean, df, th.corr=0.0001)
# data_reduced <- DaMiR.FReduct(data_reduced$data)
# There are not PCs with a correlation greater than 0.0001 with class. Please decrease 'th.corr'
data_reduced <- DaMiR.FReduct(data_clean)
# 13701 Highly correlated features have been discarded for classification. 5 Features remained. 
# Too few samples..
DaMiR.MDSplot(data_reduced, df)
df.importance <- DaMiR.FSort(data_reduced, df)
