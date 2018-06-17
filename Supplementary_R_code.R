#Supplementary R code

#R Code for the construction of Boolean Network from the RNASeq data

#pre-processing and normalization of RNASeq data. library(limma)

library(BoolNet)

#reading the data file (.csv) FileName<-GSE75748_bulk_time_course_ec.csv ReadCountMatrix<- read.csv(FileName, rownames = 1)

NormalisedData<-voom(ReadCountMatrixFiltered, normalize.method=”quantile”)

#extracting the normalized Expression value

Data<-NormalisedData$E genes<-
c("NANOG","POU5F1","SOX2","PDGFRA","FOXA2","LHX1","ZIC3","HESX1","PRDM14","H

AND1","LIFR","STAT3","PAX6","GATA3","KAT6A","ZFHX3","TBX3","BMP4","FGF2","CDC2

5C","CDK6","SALL4","L1TD1","PBX1","MYC","CTNNB1","GSK3B","TCF3","KLF4","LIF","
NEUROG1","HNF4A","UTF1","MYF5","INHBA","LEFTY1","SMAD1","FGFR2","REST","NO
G",”FRIZZLED”,”DISHEVELED”,”WNT5A”)

#filtering the 45 genes in the network FinalData<-Data[rownames(Data) %in% genes, ] #binarizing the expression data bin<-binarizeTimeSeries(FinalData)

#constructing the reconstructed network using “bestfit” method, K is max no.

# Of input interactions for a gene net<-reconstructNetwork(bin$binarizedMeasurements,method = "bestfit",maxK = 10)

#storing output in file since it is too large to be displayed fully sink("data1.txt")

print(net)

# It resulted in pethora of boolen functions for each node which was picked manually #generating a list that carry values substituted for don't care condition if any

lis<-list(c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(0),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(), c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c())

#creating index for each of the boolean function selected manually idx<-

c(3,4,2,11,11,4,16,11,190,2,122,11,4,9,10,11,11,1,7,1,234,1,1,1,40,4,11,1,11,1,1,11,10,1,11,11,1,1,1

0,1)

#constructing the final network net2<-chooseNetwork(net,idx,lis,readableFunctions = "canonical") #saving the network

toSBML(net2,"final_network.sbml") saveNetwork(net2,"final_network2")

The reconstructed network topology was studied and novel interactions were identified.

R code for the simulation of the boolean network

library(BoolNet)

#the network file was constructed containg boolean functions( given in Table S2,supplementary) #loading the Network

network1<-loadNetwork("integratednet.txt")

#simulation with different initial conditions #when only nanog is active OnlyNanog<-c(1,rep(0,44))

#when only PPOU5F1 is active OnlyPOU5F1<-c(0,1,rep(0,43))

#similarly activating genes with 1 based on their index #getting the stable states(attracters)

getAttractors(network1, startStates=OnlyNanog,type = "synchronous") #similarly done for the others

#simulation with 100 random initial states sink(“output.txt”)

getAttractors(network1, method="random", startStates=100,type = "synchronous") #from the many attractors the one with maximum basin size was choosen

R code for the perturbations

library(BoolNet)

network1<-loadNetwork("integratednet.txt") #creating a network where NANOG is knockeddown netnan<-fixGenes(net,"NANOG",0)

#viewing attractors

getAttractors(netnan, method="random", startStates=100,type = "synchronous") #choosing the attractor with maximum basin size netpou<-fixGenes(net,"POU5F1",0)

getAttractors(netpou, method="random", startStates=100,type = "synchronous") netsox<-fixGenes(net,"SOX2",0)

getAttractors(netsox, method="random", startStates=100,type = "synchronous") netposox<-fixGenes(net,c(“POU5F1”,”SOX2”),0)

getAttractors(netpousox, method="random", startStates=100,type = "synchronous") netpounan<-fixGenes(net,c(“POU5F1”,”NANOG”,0)

getAttractors(netpunan, method="random", startStates=100,type = "synchronous") netpounanl1<-fixGenes(net,c(“POU5F1”,”NANOG”,”L1TD1”),0) getAttractors(netpounanl1, method="random", startStates=100,type = "synchronous") net1<-fixGenes(net,"LITD1",0)

getAttractors(net1, method="random", startStates=100,type = "synchronous") netbmp<-fixGenes(net,"BMP4",0)

getAttractors(netbmp, method="random", startStates=100,type = "synchronous") netfg<-fixGenes(net,"FGF2",0)

getAttractors(netfg, method="random", startStates=100,type = "synchronous") netwn<-fixGenes(net,"WNT5A",0)

getAttractors(netwn, method="random", startStates=100,type = "synchronous")

#final result stored as dataframe dd were analzed with a function perturbs<-function(dd){

for(j in 1:15){

cat("\n","perturbation results:",colnames(dd)[j]); for(i in 1:44){


if(dd[i,1]=="0" & dd[i,j]=="0"){

} else if(dd[i,1]=="1" &dd[i,j]=="1"){

} else if(dd[i,1]=="0" & dd[i,j]=="1"){ cols<-rownames(dd)[i] cat("\n","upregulated:",cols)

}else if(dd[i,1]=="1" & dd[i,j]=="0"){ colu<-rownames(dd)[i] cat("\n","downregulated:",colu);
}
}}
}


R code for k-means clustering for binarization of expression data

library(limma)
library(BoolNet)
library(Binarize)

#reading the data file (.csv) FileName<-GSE75748_bulk_cell_type_ec.csv ReadCountMatrix<- read.csv(FileName, rownames = 1)

NormalisedData<-voom(ReadCountMatrixFiltered, normalize.method=”quantile”)

#extracting the normalized Expression value

Data<-NormalisedData$E

genes<-c("NANOG","POU5F1","SOX2","PDGFRA","FOXA2","LHX1","ZIC3","HESX1","PRDM14","H AND1","LIFR","STAT3","PAX6","GATA3","KAT6A","ZFHX3","TBX3","BMP4","FGF2","CDC2 5C","CDK6","SALL4","L1TD1","PBX1","MYC","CTNNB1","GSK3B","TCF3","KLF4","LIF"," NEUROG1","HNF4A","UTF1","MYF5","INHBA","LEFTY1","SMAD1","FGFR2","REST","NO
G",”FRIZZLED”,”DISHEVELED”,”WNT5A”)

#filtering the 45 genes in the network FinalData<-Data[rownames(Data) %in% genes, ]

#filtering the hESC data and geting average of the data AvgFinalData<-(FinalData$H9_rep1 + FinalData$H9_rep2 + FinalData$H9_rep3)/3

#applying k-means to the data

Bin<-binarize.kmeans(AvgFinalData) bins<-Bin@binarizedMeasurments

#saving the data org<-Bin@originalMeasurments gene.names<-rownames(AvgFinalData)

Result<-data.frame(gene.names,org,bins)

#plot for the clustering library(ggplot2) colnames
ggplot(Result,aes(x=binary, y=avg, group=Category))+

geom_point(aes(shape=Category, color=Category, size=Category))+scale_x_continuous(0,1,1)+ scale_shape_manual(values = c(17,16,10,18))+scale_size_manual(values=c(3,3,3,3))

#barplot for the results bardata<-Result[order(ourdata1$org,decreasing = TRUE),]

barplot(bardata$org,names.arg = bardata$gene.names,las=2)