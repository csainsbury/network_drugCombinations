source("~/R/_workingDirectory/_perAdmissionRewriteDataTableFunctions.R")
library(gtools)
library(igraph)
library(data.table)

id_per_location <- function(ID) {
  return(length(unique(ID)))
}

flagMove <- function(ID, charL) {
  
  charLreport <- charL
  charLnumeric <- as.numeric(factor(charL))
  
  testFrame <- data.frame(charLreport, charLnumeric)
  
  testFrame$flagMove <- 0
  testFrame$flagMove[1:nrow(testFrame)-1] <- diff(testFrame$charLnumeric)
  testFrame$nextL <- c("spacer")
  testFrame$nextL[1:(nrow(testFrame)-1)] <- charLreport[2:length(charLreport)]
  
  testFrame$charLreport <- as.character(factor(charL))
  
  outputList <- list(testFrame$charLreport, testFrame$nextL, testFrame$flagMove)
  
  return(outputList)
  
}

returnUnixDateTime<-function(date) {
  returnVal<-as.numeric(as.POSIXct(date, format="%Y-%m-%d", tz="GMT"))
  return(returnVal)
}

findSimilarDrugs <- function(inputFrame) {
  
  # inputFrame <- interestSet
  # inputFrame <- inputFrame[1:10000,]
  
  inputFrame$DrugName.original <- inputFrame$DrugName
  inputFrame$DrugNameNew <- inputFrame$DrugName
  
  inputFrame$DrugNameNew[grep("Metformin", inputFrame$DrugName, ignore.case = TRUE)] <- "Metformin"
  inputFrame$DrugNameNew[grep("Glucophage", inputFrame$DrugName, ignore.case = TRUE)] <- "Metformin"
  
  inputFrame$DrugNameNew[grep("Gliclazide", inputFrame$DrugName, ignore.case = TRUE)] <- "Gliclazide"
  inputFrame$DrugNameNew[grep("Diamicron", inputFrame$DrugName, ignore.case = TRUE)] <- "Gliclazide"
  
  inputFrame$DrugNameNew[grep("Rosiglitazone", inputFrame$DrugName, ignore.case = TRUE)] <- "Rosiglitazone"
  inputFrame$DrugNameNew[grep("Avandia", inputFrame$DrugName, ignore.case = TRUE)] <- "Rosiglitazone"
  
  inputFrame$DrugNameNew[grep("Victoza", inputFrame$DrugName, ignore.case = TRUE)] <- "Liraglutide"
  
  
  inputFrame$DrugNameNew[grep("Pioglitazone", inputFrame$DrugName, ignore.case = TRUE)] <- "Pioglitazone"
  
  inputFrame$DrugNameNew[grep("Sitagliptin", inputFrame$DrugName, ignore.case = TRUE)] <- "Sitagliptin"
  inputFrame$DrugNameNew[grep("Januvia", inputFrame$DrugName, ignore.case = TRUE)] <- "Sitagliptin"
  
  
  inputFrame$DrugNameNew[grep("Humalog Mix25", inputFrame$DrugName, ignore.case = TRUE)] <- "Humalog Mix 25"
  
  inputFrame$DrugNameNew[grep("Lantus", inputFrame$DrugName, ignore.case = TRUE)] <- "Insulin Glargine"
  inputFrame$DrugNameNew[grep("Levemir", inputFrame$DrugName, ignore.case = TRUE)] <- "Insulin Detemir"
  
  inputFrame$DrugNameNew[grep("Insulatard", inputFrame$DrugName, ignore.case = TRUE)] <- "Insulatard"
  
  inputFrame$DrugNameNew[grep("Actrapid", inputFrame$DrugName, ignore.case = TRUE)] <- "Actrapid"
  
  inputFrame$DrugNameNew[grep("Novorapid", inputFrame$DrugName, ignore.case = TRUE)] <- "Novorapid"
  
  inputFrame$DrugNameNew[grep("Novomix 30", inputFrame$DrugName, ignore.case = TRUE)] <- "Novomix 30"
  
  inputFrame$DrugNameNew[grep("Mixtard 30", inputFrame$DrugName, ignore.case = TRUE)] <- "Mixtard 30"
  inputFrame$DrugNameNew[grep("Mixtard 20", inputFrame$DrugName, ignore.case = TRUE)] <- "Mixtard 20"
  
  inputFrame$DrugNameNew[grep("Humulin M3", inputFrame$DrugName, ignore.case = TRUE)] <- "Humulin M3"
  
  inputFrame$DrugNameNew[grep("strip", inputFrame$DrugName, ignore.case = TRUE)] <- "Test Strips"
  
  outputFrame <- inputFrame
  
  outputFrame$DrugName.original <- NULL
  outputFrame$DrugName <- outputFrame$DrugNameNew
  outputFrame$DrugNameNew <- NULL
  
  return(outputFrame)
}

# generate node and link files
drugDataSet <- read.csv("~/R/GlCoSy/SDsource/Export_all_diabetes_drugs.txt",header=TRUE,row.names=NULL)

# drugDataSet <- read.csv("~/R/GlCoSy/SDsource/test_drug_out_second100kIDs_allTime.txt",header=TRUE,row.names=NULL)
 drugDataSet$BNFCode <- as.character(drugDataSet$BNFCode)
 drugDataSet$DrugName <- as.character(drugDataSet$DrugName)
 drugDataSet$LinkId <- as.numeric(levels(drugDataSet$LinkId))[drugDataSet$LinkId]
# drugDataSet$LinkId <- as.numeric(drugDataSet$LinkId)
# drugDataSet <- read.csv("./test_drug_out_second100kIDs_allTime.txt",header=TRUE,row.names=NULL)

# restrict to diabetes drugs
interestSet <- subset(drugDataSet, substr(drugDataSet$BNFCode,1,3) == "6.1" | substr(drugDataSet$BNFCode,1,4) == "0601")
interestSet <- findSimilarDrugs(interestSet)

# generate a top-100 etc list for merging back
# meeds a bit of data cleaning - merging synonymous drugs etc
n = 50
topNdrugs_DrugNames <- as.data.frame(table(interestSet$DrugName))
topNdrugs_DrugNames <- topNdrugs_DrugNames[order(topNdrugs_DrugNames$Freq), ]

topNdrugs <- tail(topNdrugs_DrugNames, n)

# merge top drugs back with interestSet to generate working data frame:
interestSet_topN_merge <- merge(interestSet, topNdrugs, by.x="DrugName", by.y="Var1")

## start drug data manipulation

drugsetDT <- data.table(interestSet_topN_merge)
drugsetDT$prescription_dateplustime1 <- returnUnixDateTime(drugsetDT$PrescriptionDateTime)
drugsetDT_original <-drugsetDT # preserve an original full dataset incase needed
# drugsetDT$LinkId <- as.numeric(levels(drugsetDT$LinkId))[drugsetDT$LinkId]

# limit analysis to 10y period 1/1/2006 - 1/1/2017
drugsetDT <- drugsetDT[prescription_dateplustime1 > returnUnixDateTime('2006-01-01') &
                 prescription_dateplustime1 < returnUnixDateTime('2017-01-01')]

# scale time to 0 to 1 range
drugsetDT$prescription_dateplustime1.original <- drugsetDT$prescription_dateplustime1
drugsetDT$prescription_dateplustime1 <- (drugsetDT$prescription_dateplustime1 - min(drugsetDT$prescription_dateplustime1)) / (max(drugsetDT$prescription_dateplustime1) - min(drugsetDT$prescription_dateplustime1))

    # master list of top n drugs to form framework for plot - later use
    masterLocTable <- as.data.frame(table(drugsetDT$DrugName)); colnames(masterLocTable) <- c("D", "f")
    masterLocTable$f <- 0

returnEdgeSet <- function(drugNames) {
  # drugNames <- drugsetDT[LinkId == "2147644003"]$DrugName
  drugNamesVector <- as.data.frame(table(drugNames), stringsAsFactors = F)$drugNames
  
  edgesOutput <- as.data.frame(combinations(length(drugNamesVector), 2, drugNamesVector))
  colnames(edgesOutput) <- c("from", "to")
  
  edgesOutput$weight <- 1
  edgesOutput$type <- c("drugPrescription")
  
  return(edgesOutput)

}

returnEdgeSet_withDateInfo <- function(drugNames, prescription_dateplustime1) {
  
  # drugNames <- drugsetDT[LinkId == "2147644003"]$DrugName
  # prescription_dateplustime1 <- drugsetDT[LinkId == "2147644003"]$prescription_dateplustime1
  
  interestTable <- data.table(drugNames, prescription_dateplustime1)
  interestTable <- unique(interestTable)
  drugNamesVector <- as.data.frame(table(drugNames), stringsAsFactors = F)$drugNames
  
  edgesOutput <- as.data.frame(combinations(length(drugNamesVector), 2, drugNamesVector))
  colnames(edgesOutput) <- c("from", "to")
  edgesOutput$medianCombinationPrescription <- 0
  
    for (jj in seq(1, nrow(edgesOutput), 1)) {
      drug1 <- edgesOutput$from[jj]
      drug2 <- edgesOutput$to[jj]
      
      allPrescriptions <- merge(interestTable[drugNames == drug1], interestTable[drugNames == drug2], by = "prescription_dateplustime1")
      
      edgesOutput$medianCombinationPrescription[jj] <-quantile(allPrescriptions$prescription_dateplustime1)[3]
      
    }
  
  # test approach
  edgesOutput$numbers <- 0
  
  
  
  edgesOutput$weight <- 1
  edgesOutput$type <- c("drugPrescription")
  
  return(edgesOutput)
  
}

edgesOutputStart <- as.data.frame(matrix(nrow=0, ncol=4))
colnames(edgesOutputStart) <- c("from", "to", "weight", "type")

IDframe <- as.data.frame(table(drugsetDT$LinkId))
IDframe$Var1 <- as.numeric(levels(IDframe$Var1))[IDframe$Var1]

##
plotfilename <- paste("./plots/single_patient_n50_","multiPlot_91kpatientSet",".pdf",sep="")
pdf(plotfilename, width=100, height=100)
par(mfrow=c(10,10))

# define colors for time plotting
maxColorValue <- 1000
palette <- colorRampPalette(c("red","blue"))(maxColorValue)

for (j in seq(1000, 4000-1, 1)) {
# for (j in seq(1, nrow(IDframe), 1)) {
  
  if (j%%100 == 0) {print(j)}
  id_subset <- drugsetDT[LinkId == IDframe$Var1[j]]
  
  if (nrow(id_subset)>1) {
    if (length(unique(id_subset$DrugName))>1) {
      
      # time bin here - ie time slices with generation of edgesets and on to images for each bin
      
      id_sub_edgeset <- returnEdgeSet_withDateInfo(id_subset$DrugName, id_subset$prescription_dateplustime1)
      
    }
  }
# }

  # remove combinations not prescribed together
  id_sub_edgeset$medianCombinationPrescription[is.na(id_sub_edgeset$medianCombinationPrescription)] <- 0
  id_sub_edgeset <- subset(id_sub_edgeset, medianCombinationPrescription > 0)

## 
links_output <- id_sub_edgeset

##

locTable_id <- as.data.frame(table(id_subset$DrugName)); colnames(locTable_id) <- c("D", "f")

locTable <- merge(masterLocTable, locTable_id, by.x="D", by.y="D", all.x = TRUE)
locTable$f.x <- NULL
locTable$f.y[is.na(locTable$f.y)] <- 0
colnames(locTable) <- c("D", "f")


nodes_output <- as.data.frame(matrix(ncol = 5, nrow = nrow(locTable)))
colnames(nodes_output) <- c("id", "media", "media.type", "type.label", "n.prescription")
nodes_output$id <- locTable$D
nodes_output$n.prescription <- locTable$f

nodes_output$media.type <- "drug"
nodes_output$type.label <- 1

nodes_output$medianPrescriptionDate <- 0

  # add median prescription data per drug
  for (k in seq(1, nrow(nodes_output), 1)) {
    
    if(nodes_output$n.prescription[k] > 0) {
      nodes_output$medianPrescriptionDate[k] <- quantile(id_subset[DrugName == nodes_output$id[k]]$prescription_dateplustime1)[3]
    }
    
  }

nodes <- nodes_output
links <- links_output

# # aggregate links - shouldn't be needed due to the links being all combinations
# links <- aggregate(links[,3], links[,-3], sum)
# links <- links[order(links$from, links$to),]
# colnames(links)[4] <- "weight"
# rownames(links) <- NULL



# Converting the data to an igraph object:
net <- graph.data.frame(links, nodes, directed=F) 

# Removing loops from the graph:
net <- simplify(net, remove.multiple = F, remove.loops = T) 

# Let's and reduce the arrow size and remove the labels:
# plot(net, edge.arrow.size=.4,vertex.label=NA)


# plotfilename <- paste("./plots/single_patient_n50_",j,".pdf",sep="")
# pdf(plotfilename, width=10, height=10)

set.seed(1234)

l <- layout_in_circle(net)

plot(net, layout=l,
     edge.curved=-0.1,
      edge.color=rgb(200, 100, 100, 60, names = NULL, maxColorValue = 255),
     # edge.color = palette[cut(E(net)$medianCombinationPrescription, maxColorValue)],
    # edge.color=rgb(200, 100, 100, E(net)$medianCombinationPrescription * 255, names = NULL, maxColorValue = 255),
     edge.width = ((E(net)$weight)/max(E(net)$weight))*15,
     vertex.size=sqrt(V(net)$n.prescription)*2,
     vertex.label.cex=0.6
)


}

dev.off()

# combinations(length(namesFrame), 2, namesFrame)
