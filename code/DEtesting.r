library(Seurat)
library(cowplot)
library(magrittr)
library(ggplot2)
library(dplyr)
load("./spcNew20200306.Rdata")
Idents(spcNew) <- "celltype"

# healthy_vs_OAa
## 1.SSCs
SSCs_OAa <- FindMarkers(spcNew, ident.2 = "SSCs_Healthy", ident.1 = "SSCs_OAa",
                                min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
SSCs_OAa <- SSCs_OAa %>% mutate(SYMBOL = rownames(SSCs_OAa), Celltype = "SSCs", Sample = "OAa")
                     %>% filter(p_val_adj < 0.05)  

## 2.Sgonia
Sgonia_OAa <- FindMarkers(spcNew, ident.2 = "Sgonia_Healthy" , ident.1 = "Sgonia_OAa", 
                                  min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
Sgonia_OAa <- Sgonia_OAa %>% mutate(SYMBOL = rownames(Sgonia_OAa), Celltype = "Sgonia", Sample = "OAa") 
                         %>% filter(p_val_adj < 0.05) 


## 3. Scytes
Scytes_OAa <- FindMarkers(spcNew, ident.2 = "Scytes_Healthy" , ident.1 = "Scytes_OAa", 
                                  min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
Scytes_OAa <- Scytes_OAa %>% mutate(SYMBOL = rownames(Scytes_OAa), Celltype = "Scytes", Sample = "OAa") 
                         %>% filter(p_val_adj < 0.05) 


# somatic cells
## 4.Leydig
Leydig_OAa <- FindMarkers(spcNew, ident.2 = "Leydig_Healthy" , ident.1 = "Leydig_OAa",
                                  min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
Leydig_OAa <- Leydig_OAa %>% mutate(SYMBOL = rownames(Leydig_OAa), Celltype = "Leydig", Sample = "OAa") 
                         %>% filter(p_val_adj < 0.05)


## 5.Macrophage
Macrophage_OAa <- FindMarkers(spcNew, ident.2 = "Macrophage_Healthy" , ident.1 = "Macrophage_OAa",min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
Macrophage_OAa <- Macrophage_OAa %>% mutate(SYMBOL = rownames(Macrophage_OAa), Celltype = "Macrophage", Sample = "OAa") %>% filter(p_val_adj < 0.05) 


## 6.Myoid
Myoid_OAa <- FindMarkers(spcNew, ident.2 = "Myoid_Healthy" , ident.1 = "Myoid_OAa",
                                 min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
Myoid_OAa <- Myoid_OAa %>% mutate(SYMBOL = rownames(Myoid_OAa), Celltype = "Myoid", Sample = "OAa") 
                       %>% filter(p_val_adj < 0.05) 


## 7. Endothelial
Endothelial_OAa <- FindMarkers(spcNew, ident.2 = "Endothelial_Healthy" , ident.1 = "Endothelial_OAa",
                                       min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
Endothelial_OAa <- Endothelial_OAa %>% mutate(SYMBOL = rownames(Endothelial_OAa), Celltype = "Endothelial", Sample = "OAa") 
                                   %>% filter(p_val_adj < 0.05) 



# healthy_vs_OAb
## 1.SSCs
SSCs_OAb <- FindMarkers(spcNew, ident.2 = "SSCs_Healthy", ident.1 = "SSCs_OAb", 
                                min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
SSCs_OAb <- SSCs_OAb %>% mutate(SYMBOL = rownames(SSCs_OAb), Celltype = "SSCs", Sample = "OAb") 
                     %>% filter(p_val_adj < 0.05)  

## 2.Sgonia
Sgonia_OAb <- FindMarkers(spcNew, ident.2 = "Sgonia_Healthy" , ident.1 = "Sgonia_OAb", 
                                  min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
Sgonia_OAb <- Sgonia_OAb %>% mutate(SYMBOL = rownames(Sgonia_OAb),Celltype = "Sgonia", Sample = "OAb") 
                         %>% filter(p_val_adj < 0.05) 


## 3.Scytes
Scytes_OAb <- FindMarkers(spcNew, ident.2 = "Scytes_Healthy" , ident.1 = "Scytes_OAb", 
                                  min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
Scytes_OAb <- Scytes_OAb %>% mutate(SYMBOL = rownames(Scytes_OAb), Celltype = "Scytes", Sample = "OAb") 
                         %>% filter(p_val_adj < 0.05) 


#somatic cells
## 4.Leydig
Leydig_OAb <- FindMarkers(spcNew, ident.2 = "Leydig_Healthy" , ident.1 = "Leydig_OAb",
                                  min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
Leydig_OAb <- Leydig_OAb %>% mutate(SYMBOL = rownames(Leydig_OAb), Celltype = "Leydig", Sample = "OAb") 
                         %>% filter(p_val_adj < 0.05)


## 5.Macrophage
Macrophage_OAb <- FindMarkers(spcNew, ident.2 = "Macrophage_Healthy" , ident.1 = "Macrophage_OAb",
                                      min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
Macrophage_OAb <- Macrophage_OAb %>% mutate(SYMBOL = rownames(Macrophage_OAb), Celltype = "Macrophage", Sample = "OAa") 
                                 %>% filter(p_val_adj < 0.05) 


## 6.Myoid
Myoid_OAb <- FindMarkers(spcNew, ident.2 = "Myoid_Healthy" , ident.1 = "Myoid_OAb",
                                 min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
Myoid_OAb <- Myoid_OAb %>% mutate(SYMBOL = rownames(Myoid_OAb), Celltype = "Myoid", Sample = "OAb") 
                       %>% filter(p_val_adj < 0.05) 


## 7. Endothelial
Endothelial_OAb <- FindMarkers(spcNew, ident.2 = "Endothelial_Healthy" , ident.1 = "Endothelial_OAb",
                                       min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
Endothelial_OAb <- Endothelial_OAb %>% mutate(SYMBOL = rownames(Endothelial_OAb), Celltype = "Endothelial", Sample = "OAb") 
                                   %>% filter(p_val_adj < 0.05) 



