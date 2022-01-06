#Code derived from Prof.Steve Hovarth's brilliant tutorials found here - https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/index.html#

rm(list = ls(all = TRUE))

getwd()
R.Version()
.libPaths()
install.packages("dynamicTreeCut")
install.packages("fastcluster")
install.packages("WGCNA")
BiocManager::install("preprocessCore")
BiocManager::install("impute")

library(WGCNA)

#WGCNA package setting#
##1a##
#Removing auxiliary data and transposing the expression data for further analysis#
getwd()
rna_norm_counts <- read.csv("~/Bioinformatics MSc Resources/Module5_ComputationalBioforComplexSystems/rna_norm_counts.csv")
View(rna_norm_counts)
options(stringsAsFactors = FALSE);

datExpr0 <- as.data.frame(t(rna_norm_counts[, -c(1)]));
names(datExpr0) <- rna_norm_counts$X;
rownames(datExpr0) <- names(rna_norm_counts)[-c(1)];

##1b##
#Checking data for excessive missing values and ID of outlier microarray samples#

BiocManager::install("GO.db")
library(WGCNA)

#Check for genes and samples with too many missing values#
getAnywhere(goodSamplesGenes)
gsg <- goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
#The values returned is TRUE, therefore all genes have passed the cut - none need to be removed#

#Clustering of samples to detect any obvious outliers#

sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(24,20)
par(cex = 0.6);
par(mar = c(0,5,2.4,2))
plot(sampleTree, main = "Clustering dendrogram of samples based on their Euclidean distance", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
abline(h = 0.5e+06, col = "red");

clust = cutreeStatic(sampleTree, cutHeight = 0.5e+06, minSize = 0)
table(clust)
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

head(datExpr)
##1c##
#Loading in our sample data sheet to try and relate that to gene expression#

sample_sheet <- read.csv("~/Bioinformatics MSc Resources/Module5_ComputationalBioforComplexSystems/sample_sheet.csv")

str(sample_sheet)

sample_sheet = sample_sheet[c(1:2)]
str(sample_sheet)

RiverSamples <- rownames(datExpr);
traitRows <- match(RiverSamples, sample_sheet$SampleID);
head(traitRows)
datsamples <- sample_sheet[traitRows,-1];
datsamples <- as.matrix(datsamples)
rownames(datsamples) <- sample_sheet[traitRows, 1];


collectGarbage();

#The expression data is in the 'datexpr' variable, and the corresponding sample traits in the variable 'datSamples'.#
#Before continuing to network construction and module detection, we visualise how the sample features relate to the sample dendrogram#

#Re-cluster samples#

#editing datsamples matrix to remove all characters so it's fully numeric# 
write.csv(datsamples,"datsamples.csv")


datsamples <- read.csv("~/Bioinformatics MSc Resources/Module5_ComputationalBioforComplexSystems/datsamples.csv")

sampleTree2 = hclust(dist(datExpr), method = "average")
typeof(datsamples)
traitColors = numbers2colors(datsamples, signed = FALSE);

plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datsamples),
                    main = "Sample dendrogram and trait heatmap")

save(datExpr, datsamples, file="Part1Jan.Rdata")

##2 - Automatic construction of the gene network and ID of modules#
#2.a.1 Choosing the soft-thresholding power: analysis of network topology#

options(stringsAsFactors = FALSE);
lnames = load(file = "Part1Jan.RData");
lnames

print(datsamples)
#Choose a set of soft-thresholding powers#
powers = c(c(1:10), seq(from = 12, to=20, by=2))

#Call the network topology analysis function#
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

#Plot the results#
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")

#Mean connectivity as a function of the soft-thresholding power#
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


##2.a.2##
#One-step network construction and module detection#
library(WGCNA)
net = blockwiseModules(datExpr, power = 4,maxBlockSize = 10000, 
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "RNA_normcounts_Jan",
                       verbose = 3)

table(net$colors)

net$dendrograms
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "Part2Jan.RData")

##3 - Relating modules to external clinical traits##
#3.a Quantifying module-trait associations#

#Here the aim is to ID modules that have significant association with different traits, i.e. river location.#

load(file = "Part2Jan.RData")
# Defining the numbers of genes and samples#
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

# Recalculating MEs with colour labels#
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datsamples
, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
#Constructing a figure to assist in interpreting the large number of modules and traits#
#Each association is colour coded with a correlation value#

# Display correlations and their p-values#
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values via a heatmap plot#
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datsamples),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#3.b Gene relationship to important modules: Gene significance and Module membership#

#We use this step to see if there's correlation between different sample extraction sites and module membership#

#In the following code, associations of individual genes with location are define by Gene significance#
#Also, for each module, the quantitative measure of module membership MM is the correlation of the module eigengene and the gene expression profile.#

#All this allows us to quantify the similarity of all genes to every module.#

# Define variable site containing the site column of datsamples
site = as.data.frame(datsamples$Site);
names(site) = "site"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, site, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(site), sep="");
names(GSPvalue) = paste("p.GS.", names(site), sep="");

#3.c Intramodular analysis: identifying genes with high GS and MM#

module = "pink"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for extraction location(site(",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#correlation is very high at 0.63!Genes are significantly associated with this trait#

##4. Visualising networks in R##

##4.a Visualising the gene network##

#One way to visualize a weighted network is to plot its heatmap.#
#Each row and column of the heatmap correspond to a single gene.#
#The heatmap shows adjacencies or topological overlaps, with light colors denoting low adjacency (overlap) and darker colors higher adjacency (overlap).#
# The gene dendrograms and module colors are plotted along the top and left side of the heatmap.#

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

##4.b Visualising the eigengene network##

#Visualising the eigengene network representing relationships among modules and between site modules.#
#The following dendrogram and heatmap ID groups of correlated eigengenes - meta-modules.#

# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
site = as.data.frame(datsamples$site);
names(site) = "site"
# Add the site to existing module eigengenes
MET = orderMEs(cbind(MEs, site))

#Plot the relationships among the eigengenes and the trait#
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)

#We split the dendrogram and heatmap plots using the following code#
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

##5 Export network data to Cytoscape##

#Cytoscape allows one to input an edge file and a node file#

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 6);
# Select modules
modules = c("brown", "red");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = rna_norm_counts$X[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);