library(Xeva)

# Define your desired directory
output_dir <- "/Users/guanqiaofeng/Documents/BHK/Xeva/Xevacurate_KRAS/results/explore"

# Check if the directory exists; if not, create it
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Set working directory
setwd(output_dir)

# read in xevaset
x.set <- readRDS("../xevaset/UHN_Tsao_Lung_2022_v1.rds")
# compute resposne
#x.set <- setResponse(x.set, res.measure = c("mRECIST", "slope", "AUC", "angle", "abc", "TGI", "lmm"), verbose=FALSE)
x.set <- setResponse(x.set, res.measure = "slope") # no error
x.set <- setResponse(x.set, res.measure = "AUC") # no error
x.set <- setResponse(x.set, res.measure = "angle") # no error
x.set <- setResponse(x.set, res.measure = "abc") # no error
x.set <- setResponse(x.set, res.measure = "mRECIST") # >=50 warnings
x.set <- setResponse(x.set, res.measure = "TGI") # error: Error in x[[jj]][iseq] <- vjj : replacement has length zero 
x.set <- setResponse(x.set, res.measure = "lmm") # error: Error in MEEM(object, conLin, control$niterEM) : NA/NaN/Inf in foreign function call (arg 1)

saveRDS(x.set, file = "../xevaset/UHN_Tsao_Lung_2022_DrugResponse_v1.rds")

# xset.mr <- summarizeResponse(x.set, response.measure = "mRECIST")
# write.table(xset.mr, file = "xset_mr.tsv", sep = "\t", quote = FALSE)

#### start here ####

x.set <- readRDS("../xevaset/UHN_TNBC_ALL_202507_DrugResponse_added.rds")

# access model data
x.mod <- modelInfo(x.set)
dim(x.mod)
x.mod[1:4, ]

# access experiment for a model
model.data <- getExperiment(x.set, model.id='REF004.ABT-263.m1')
head(model.data)

# access batches
batch.name <- batchInfo(x.set)
batch.name[1:4]
batchInfo(x.set, batch="REF004.ABT-263")

# visualize PDX growth curve
plotPDX(x.set, batch = "01-1104.CX-5461")
# normalize volume, change color of lines
plotPDX(x.set, batch = "01-1104.CX-5461", 
        vol.normal = TRUE, 
        control.col = "#a6611a", 
        treatment.col = "#018571",
        major.line.size = 1,
        max.time = 100)

# change color of lines, error visualization with errorbar
plotPDX(x.set, batch = "01-1104.CX-5461", 
        control.col = "#a6611a", 
        treatment.col = "#018571",
        major.line.size = 1,
        max.time = 100,
        SE.plot='errorbar')

# normalize volume, change color of lines, error visualization with ribbon
plotPDX(x.set, batch = "01-1104.CX-5461", 
        vol.normal = TRUE, 
        control.col = "#a6611a", 
        treatment.col = "#018571",
        major.line.size = 1,
        max.time = 100,
        SE.plot='ribbon')

# NOT normalize volume, change color of lines, error visualization with ribbon
plotPDX(x.set, batch = "01-1104.CX-5461", 
        vol.normal = FALSE, 
        control.col = "#a6611a", 
        treatment.col = "#018571",
        major.line.size = 1,
        max.time = 100,
        SE.plot='ribbon')



# in patient level
plotPDX(x.set, patient.id = "73230", drug = "AZD-5305", control.name = "H2O") # this function is working
# is it possible to plot three treatments (AZD-8201, ADC-Control, H2O), no for current xeva package

# PDX model drug response
x.mr <- summarizeResponse(x.set,  response.measure = "mRECIST", group.by = "patient.id")
x.mr <- summarizeResponse(x.set, response.measure = "AUC", model.id = "69693.AZD-5305.m1") ### not working/working for earlier TNBC xevaset
waterfall(x.set, drug="Carboplatin", res.measure="best.average.response")

### summarizeMolecularProfiles not working properly
mut <- summarizeMolecularProfiles(x.set,drug = "BMN-673", mDataType="RNASeq") ### not working
mut <- summarizeMolecularProfiles(x.set,drug = "BMN-673", mDataType="mutation") ### not working
### Only the H2O and Carboplatin work, but also have issues in the output
mut <- summarizeMolecularProfiles(x.set,drug = "Carboplatin", mDataType="mutation")
mut <- summarizeMolecularProfiles(x.set,drug = "H2O", mDataType="mutation")

# angle between treatment and control
response(x.set, batch="16720_P6_5305_RES.BMN-673", res.measure="angle")
# linear mixed-effects model (lmm)
response(x.set, batch="16720_P6_5305_RES.BMN-673", res.measure="lmm")

# work for no duplication case
data(brca)
brca.mr <- summarizeResponse(brca, response.measure = "mRECIST")
brca.mr[1:5, 1:4]
waterfall(brca, drug="binimetinib", res.measure="best.average.response")
mut <- summarizeMolecularProfiles(brca,drug = "BYL719", mDataType="mutation")
model.type <- Biobase::exprs(mut)["CDK13", ]
model.type[grepl("Mut", model.type)] <- "mutation"
model.type[model.type!="mutation"] <- "wild type"
model.color <- list("mutation"="#b2182b", "wild type"="#878787")
waterfall(brca, drug="BYL719", res.measure="best.average.response",
          model.id=names(model.type), model.type= model.type,
          type.color = model.color)

# exploring repdx dataset 
# does not work for case where there are duplications
data(repdx)
repdx <- setResponse(repdx, res.measure = c("mRECIST", "slope", "AUC", "angle", "abc"), verbose=FALSE)
repdx <- setResponse(repdx, res.measure = "slope") # no error
plotPDX(repdx, batch = "P3", SE.plot = "errorbar") # works
plotPDX(repdx, batch = "P4", vol.normal = TRUE, SE.plot = "ribbon") # works
repdx.mr <- summarizeResponse(repdx, response.measure = "slope") # error
waterfall(repdx, drug="binimetinib", res.measure="best.average.response") # error
response(repdx, batch="P1", res.measure="angle") # works

### Gene drug association not working, related with summarizeMolecularProfiles
drugSensitivitySig(object=x.set, drug="H2O", mDataType="RNASeq",
                   features=c(1,2), sensitivity.measure="slope", fit="lm")


# brca dataset works
drugSensitivitySig(object=brca, drug="tamoxifen", mDataType="RNASeq",
                   features=c(1,2,3,4,5),
                   sensitivity.measure="best.average.response", fit = "lm")


