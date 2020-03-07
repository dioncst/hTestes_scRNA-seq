# Load libraries
library(Seurat)
library(cowplot)
library(magrittr)
library(ggplot2)
library(dplyr)
options(future.globals.maxSize = 10000 * 1024^2)

# healthy
## Read in data
healthyData <- read.table(file = "./GSE112013_Combined_UMI_table.txt.gz", 
                          header = T, 
                          row.names=1, 
                          sep="", as.is=T)
## Standard Workflow Dataset preprocessing                       
healthy <- CreateSeuratObject(counts = healthyData, 
                              project = "Healthy", 
                              min.cells = 3, 
                              min.features = 200)
healthy[["percent.mt"]] <- PercentageFeatureSet(healthy, pattern = "^MT-")
healthy <- subset(healthy, subset = nFeature_RNA > 200 & percent.mt < 20)
healthy <- NormalizeData(healthy, verbose = FALSE)
healthy <- FindVariableFeatures(healthy, selection.method = "vst", nfeatures = 2000)

# OAa
## Read in data
OAaData <- Read10X("./A1/")
## Standard Workflow Dataset preprocessing
OAa <- CreateSeuratObject(counts = OAaData, 
                          project = "OAa",
                          min.cells = 3, 
                          min.features = 200)
OAa[["percent.mt"]] <- PercentageFeatureSet(OAa, pattern = "^MT-")
OAa <- subset(OAa, subset = nFeature_RNA > 200 & percent.mt < 20)
OAa <- NormalizeData(OAa, verbose = FALSE)
OAa <- FindVariableFeatures(OAa, selection.method = "vst", nfeatures = 2000)

# OAb
## Read in data
OAbData <- Read10X("./A2/")
## Standard Workflow Dataset preprocessing
OAb <- CreateSeuratObject(counts = OAbData, 
                          project = "OAb",
                          min.cells = 3, 
                          min.features = 200)
OAb[["percent.mt"]] <- PercentageFeatureSet(OAb, pattern = "^MT-")
OAb <- subset(OAb, subset = nFeature_RNA > 200 & percent.mt < 20)
OAb <- NormalizeData(OAb, verbose = FALSE)
OAb <- FindVariableFeatures(OAb, selection.method = "vst", nfeatures = 2000)

# Integration of 3 pancreatic islet cell datasets
tList = list(OAa, OAb, healthy)
tAnchors <- FindIntegrationAnchors(object.list = tList, dims = 1:30)
tCombined <- IntegrateData(anchorset = tAnchors, dims = 1:30)
# set during IntegrateData
DefaultAssay(tCombined) <- "integrated"

# Run the standard workflow for visualization and clustering
tCombined <- ScaleData(tCombined, verbose = FALSE)
tCombined <- RunPCA(tCombined, npcs = 30, verbose = FALSE)
tCombined <- RunUMAP(tCombined, reduction = "pca", dims = 1:30)
tCombined <- FindNeighbors(tCombined, reduction = "pca", dims = 1:20)
tCombined <- FindClusters(tCombined, resolution = 0.6)

#assign celltype
DotPlot(spcNew, features = c("ID4", "UTF1",  # SSCs
                             "DMRT1", "KIT",  # Sgonia
                             "DMC1", "SYCP3",  # Scytes 
                             "SUN5", "BRDT", "CAMK4",  # RS
                             "CCDC179", "TNP1",  # ES
                             "PRM2",  # Stids
                             "SOX9", "AMH",  # Sertoli
                             "DLK1", "INHBA",  # Leydig
                             "MYH11","ACTA2",  # Myoid
                             "VWF",  # Endothelial
                             "CD14"  # Macrophage  
                              ))
tCombined <- RenameIdents(tCombined, "4" = "SSCs", "1" = "SSCs", "5" = "Sgonia", 
                                        "13" = "Scytes", "10" = "Scytes","12" = "Scytes",
                                        "8" = "RS", "7" = "ES", "14" = "ES",
                                        "6" = "Stids", "3" = "Endothelial", "11" = "Macrophage",
                                        "0" = "Leydig", "2" = "Myoid", 
                                        "9" = "Myoid", "15" = "Myoid")

spcNew <- subset(tCombined, ident = c("SSCs", "Sgonia", "Scytes", "RS", "ES", "Stids", 
                                      "Leydig", "Myoid", "Endothelial", "Macrophage"))
levels(spcNew) <- c("SSCs", "Sgonia", "Scytes", "RS", "ES", "Stids", 
                    "Leydig", "Myoid", "Endothelial", "Macrophage")

spcNew$celltype <- paste(Idents(spcNew), spcNew$orig.ident, sep = "_")
spcNew$celltype.orig <- Idents(spcNew) ##保存原本的Idents为celleype
Idents(spcNew) <- "celltype"

save(spcNew,file = "./spcNew20200306.Rdata")

