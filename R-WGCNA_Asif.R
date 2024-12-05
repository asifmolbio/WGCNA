# Load necessary Bioconductor installer and install specific packages for functional enrichment and preprocessing
source("http://bioconductor.org/biocLite.R")
biocLite(c("GO.db", "preprocessCore", "impute"))  # Install functional annotation and preprocessing packages

# Install required R packages for various statistical and WGCNA-related functions
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "reshape", 
                   "fastcluster", "dynamicTreeCut", "survival"))  # Utilities for data processing and clustering
install.packages("WGCNA")  # Main package for Weighted Gene Co-Expression Network Analysis
install.packages("scatterplot3d")  # For 3D scatter plots

# Load the WGCNA library
library(WGCNA)

# Set the working directory to where input data is located (update with your path)
setwd("C:\\Users\\AsifMolBio\\Desktop\\WGCNA-R")  # Adjust this to your own data folder

# Load expression data from a file
inputdata1="tcga.txt"  # Update with your own dataset filename
data0=read.table(inputdata1, sep="\t", row.names=1, header=TRUE, check.names=FALSE, quote="!")  # Read expression data
datSummary=rownames(data0)  # Extract row names (gene names) for reference
datExpr = t(data0)  # Transpose the data for WGCNA (samples as rows, genes as columns)
no.samples = dim(datExpr)[[1]]  # Get the number of samples
dim(datExpr)  # Check the dimensions of the dataset

# Selecting the soft-thresholding power
powers1=c(seq(1,10,by=1), seq(12,20,by=2))  # Range of power values to test
RpowerTable=pickSoftThreshold(datExpr, powerVector=powers1)[[2]]  # Identify suitable power for scale-free topology
cex1=1
par(mfrow=c(1,2))  # Create a 2-panel plot layout
pdf("beta.pdf")  # Save plot to PDF
plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2], xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2", type="n")  # Plot scale-free topology model fit
text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2], labels=powers1, cex=cex1, col="red")
abline(h=0.85, col="red")  # Add a threshold line for good fit (e.g., R^2 > 0.85)
plot(RpowerTable[,1], RpowerTable[,5], xlab="Soft Threshold (power)", 
     ylab="Mean Connectivity", type="n")  # Plot mean connectivity
text(RpowerTable[,1], RpowerTable[,5], labels=powers1, cex=cex1, col="red")
dev.off()

# Constructing adjacency matrix
beta1=6  # Choose a soft-thresholding power (modify based on the above plot)
Connectivity=softConnectivity(datExpr, power=beta1)  # Compute network connectivity using selected power
pdf("scalefree.pdf", 15, 10)  # Save scale-free topology plot to a PDF
par(mfrow=c(1,1))
scaleFreePlot(Connectivity, main=paste("soft threshold, power=", beta1), truncated=TRUE)  # Assess scale-free fit
dev.off()

# Module detection
ConnectivityCut = 1000  # Filter to top connected genes (adjust as needed)
ConnectivityRank = rank(-Connectivity)  # Rank genes by connectivity
restConnectivity = ConnectivityRank <= ConnectivityCut  # Subset to top-ranked genes
ADJrest = adjacency(datExpr[,restConnectivity], power=beta1)  # Compute adjacency matrix
dissTOM = TOMdist(ADJrest)  # Calculate Topological Overlap Matrix (TOM) distance
hierTOM = hclust(as.dist(dissTOM), method="average")  # Hierarchical clustering of TOM
colorh1 = cutreeStaticColor(hierTOM, cutHeight=0.8, minSize=3)  # Assign module colors
pdf("module.pdf")
par(mfrow=c(2,1), mar=c(2,4,1,1))
plot(hierTOM, main="Cluster Dendrogram", labels=FALSE, xlab="", sub="")  # Dendrogram
plotColorUnderTree(hierTOM, colors=data.frame(module=colorh1))  # Plot module colors
title("Module (branch) color")
dev.off()

# TOM plot
pdf("TOM.pdf")
TOMplot(dissTOM, hierTOM, colorh1, terrainColors=TRUE)  # Visualize TOM with modules
dev.off()

# Multidimensional scaling plot (MDS)
pdf("cmd.pdf")
cmd1=cmdscale(as.dist(dissTOM), 3)  # Perform classical MDS
pairs(cmd1, col=as.character(colorh1), main="MDS plot")  # Pairwise scatterplots
dev.off()

# 3D scatterplot
library(scatterplot3d)
pdf("3d.pdf")
par(mfrow=c(1,1), mar=c(4,3,2,3)+0.1)
scatterplot3d(cmd1, color=colorh1, angle=250, 
              xlab="Scaling Axis 1", ylab="Scaling Axis 2", zlab="Scaling Axis 3")
dev.off()

# Module eigengenes analysis
datME = moduleEigengenes(datExpr[,restConnectivity], colorh1)[[1]]  # Calculate module eigengenes
dissimME = 1 - (t(cor(datME, method="p"))) / 2  # Dissimilarity matrix of module eigengenes
hclustdatME = hclust(dist(dissimME), method="average")  # Cluster module eigengenes
pdf("modul_cluster.pdf")
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based on the module eigengenes of modules")  # Dendrogram of module eigengenes
dev.off()

# Module correlation plot
pdf("modul_cor.pdf")
pairs(datME)  # Pairwise scatterplots of module eigengenes
dev.off()

modul <- signif(cor(datME, use="p"), 2)  # Correlation matrix for module eigengenes
write.table(modul, "modul_cor.txt", sep="\t", quote=FALSE)

# Export gene module memberships
datME = moduleEigengenes(datExpr, colorh1)[[1]]  # Recalculate module eigengenes
color1 = rep("grey", dim(datExpr)[[2]])  # Assign module colors
color1 = as.character(colorh1)
datKME = signedKME(datExpr, datME)  # Compute module membership
datout = data.frame(datSummary, colorNEW=color1, datKME)  # Combine with gene information
write.table(datout, "gene_module.xls", sep="\t", row.names=FALSE, quote=FALSE)

# Export network to Cytoscape
exportNetworkToCytoscape(ADJrest, edgeFile="edge.txt", nodeFile="node.txt", threshold=0.5)  # Export for visualization

# Relating modules to clinical traits
inputclinial = "clinical_data.txt"  # Clinical data file
dataclinial = read.table(inputclinial, sep="\t", row.names=1, header=TRUE, check.names=FALSE, quote="!")  # Load clinical data
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr, colorh1)$eigengenes  # Module eigengenes
MEsFemale = orderMEs(MEs0)  # Order eigengenes
modul_clinical_cor = cor(MEsFemale, dataclinial, use="p")  # Correlation with clinical traits
write.table(modul_clinical_cor, "module-clinial-cor.xls", sep="\t", quote=FALSE)
modul_clinical_p = corPvalueStudent(modul_clinical_cor, nSamples)  # P-values
write.table(modul_clinical_p, "modul-clinical-p.xls", sep="\t", quote=FALSE)

# Heatmap of module-trait relationships
textMatrix = paste(signif(modul_clinical_cor, 2), " (", signif(modul_clinical_p, 1), ")", sep="")
dim(textMatrix) = dim(modul_clinical_cor)

pdf("modul-clinical.pdf")
par(mar=c(6, 8.5, 3, 3))  # Set margins for the heatmap plot
labeledHeatmap(Matrix=modul_clinical_cor,  # Module-trait correlation matrix
               xLabels=names(dataclinial),  # Clinical trait names
               yLabels=names(MEsFemale),  # Module eigengene names
               ySymbols=names(MEsFemale),  # Symbols for module eigengenes
               colorLabels=FALSE,  # Use color labels (FALSE for text-based labels)
               colors=greenWhiteRed(50),  # Heatmap color gradient
               textMatrix=textMatrix,  # Annotated correlation matrix with p-values
               setStdMargins=FALSE, 
               cex.text=1,  # Font size for text annotations
               zlim=c(-1, 1),  # Color scale limits
               main="Module-trait relationships")  # Plot title
dev.off()
