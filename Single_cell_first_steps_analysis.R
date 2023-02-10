
## ---------------------------
## The script process single cell RNA-seq data from Brouiller et al..
##
## The script takes cellranger output directory, creates a Seurat object and performs 
## downstream analysis steps (normalization, PCA, clustering, differential expression analysis).
## 
## Author: Francesca Nadalin
##
## ---------------------------

# paths
PATIENT_32_DATA <- "../../../data/2007566/count_180817_P32_exvivoDC__Homo_sapiens/filtered_gene_bc_matrices/GRCh38/"
OBJ_AXLDC <- "object_AxlDC_merged.Robj"
OBJ_CDC2 <- "object_cDC2_merged.Robj"
FUNCTIONS <- "/data/users/ext-fnadalin/curie/scripts/R/seurat_scripts/functions.R"
ISG <- "../data/gene_signatures/List_Human_ISG_curated.txt"
TCELL <- "../data/gene_signatures/180907_T_cell_help.txt"
CL2 <- "../data/cluster_signatures/Signature_cluster2_top100_DEGs.txt"
CL3 <- "../data/cluster_signatures/Signature_cluster3_top100_DEGs.txt"
CL4 <- "../data/cluster_signatures/Signature_cluster4_top100_DEGs.txt"
OUTDIR <- "out"

dir.create(OUTDIR, showWarnings = FALSE)

# libs
source(FUNCTIONS)

# cluster P32 cells
data <- Read10X(PATIENT_32_DATA)
object <- CreateSeuratObject(raw.data = data)
object <- NormalizeData(object)
object <- ScaleData(object)
object <- FindVariableGenes(object)
object <- RunPCA(object)
object <- FindClusters(object, dims.use = 1:15, resolution = 0.3)
object@meta.data$orig.ident <- object@meta.data$res.0.3
object <- SetIdent(object, ident.use = object@meta.data$orig.ident)

allMarkers <- FindAllMarkers(object) 
write.table(allMarkers, file = file.path(OUTDIR, "allMarkers.tsv"), sep = "\t", quote = FALSE)

saveRDS(object, file = file.path(OUTDIR, "object_P32_all.Rds"))

# remove contaminants and assign cell IDs

cells <- object@cell.names[object@ident %in% 0:3]
object <- SubsetData(object, cells.use = cells, subset.raw = TRUE)
object <- RenameClusters(object = object, renaming = "cDC2 HIV1+,cDC1 HIV1+,pDC HIV1+,Axl+DC HIV1+", id = "res.0.3")
object@meta.data$orig.ident <- object@meta.data$res.0.3

saveRDS(object, file = file.path(OUTDIR, "object_P32.Rds"))

# load Axl+DC and cDC2 from healthy donor

obj_axldc <- LoadObject(OBJ_AXLDC)
cells <- colnames(obj_axldc@data)[obj_axldc@meta.data$orig.ident == "MOCK-24h"]
obj_axldc <- SubsetData(obj_axldc, cells.use = cells, subset.raw = TRUE, do.clean = TRUE)
obj_axldc@meta.data$orig.ident <- rep("Axl+DC HIV1-", length(cells))
obj_axldc <- SetIdent(obj_axldc, ident.use = obj_axldc@meta.data$orig.ident)
new.names <- gsub("^MOCK-24h", "AxlDC", obj_axldc@cell.names)
colnames(obj_axldc@raw.data) <- rownames(obj_axldc@meta.data) <- obj_axldc@cell.names <- colnames(obj_axldc@data) <- new.names

obj_cdc2 <- LoadObject(OBJ_CDC2)
cells <- colnames(obj_cdc2@data)[obj_cdc2@meta.data$orig.ident == "MOCK-24h"]
obj_cdc2 <- SubsetData(obj_cdc2, cells.use = cells, subset.raw = TRUE, do.clean = TRUE)
obj_cdc2@meta.data$orig.ident <- rep("cDC2 HIV1-", length(cells))
obj_cdc2 <- SetIdent(obj_cdc2, ident.use = obj_cdc2@meta.data$orig.ident)
new.names <- gsub("^MOCK-24h", "cDC2", obj_cdc2@cell.names)
colnames(obj_cdc2@raw.data) <- rownames(obj_cdc2@meta.data) <- obj_cdc2@cell.names <- colnames(obj_cdc2@data) <- new.names

# merge

object_merge <- MergeSeurat(object1 = object, object2 = obj_axldc, do.normalize = FALSE, do.scale = FALSE, do.center = FALSE)
object_merge <- MergeSeurat(object1 = object_merge, object2 = obj_cdc2, do.normalize = TRUE, do.scale = TRUE, do.center = TRUE)

# compute ISG and T-cell help signatures

isg <- read.table(ISG, sep = "\t")[,1]
tcell <- read.table(TCELL, sep = "\t")[,1]
gene_lists <- list("List Human ISG" = isg, "180907 Gene T_cell_help" = tcell) 

object_merge <- AddModuleScore(object_merge, genes.list = gene_lists)
idx <- which(colnames(object_merge@meta.data) == "Cluster1")
colnames(object_merge@meta.data)[idx] <- names(gene_lists)[1]
idx <- which(colnames(object_merge@meta.data) == "Cluster2")
colnames(object_merge@meta.data)[idx] <- names(gene_lists)[2]

object_merge <- SetIdent(object_merge, ident.use = object_merge@meta.data$orig.ident)
saveRDS(object_merge, file = file.path(OUTDIR, "object_patient_Axldc_cdc2.Rds"))

# differential expression analysis

cdc2 <- file.path(OUTDIR, "cDC2_DEGs.tsv")
axldc <- file.path(OUTDIR, "AxlDC_DEGs.tsv")

GeneMarkersTable(object = object_merge, out.name = cdc2, ident.1 = "cDC2 HIV1+", ident.2 = "cDC2 HIV1-", test.use = "MAST")
GeneMarkersTable(object = object_merge, out.name = axldc, ident.1 = "Axl+DC HIV1+", ident.2 = "Axl+DC HIV1-", test.use = "MAST")

# subset on Axl+DC HIV+ cells and compute cluster module score

cells <- object_merge@cell.names[object_merge@ident == "Axl+DC HIV1+"]
object_axldc_patient <- SubsetData(object_merge, cells.use = cells, subset.raw = TRUE)

cl2 <- read.table(CL2, sep = "\t")[,1]
cl3 <- read.table(CL3, sep = "\t")[,1]
cl4 <- read.table(CL4, sep = "\t")[,1]
gene_lists <- list("cluster 2" = cl2, "cluster 3" = cl3, "cluster 4" = cl4)

object_axldc_patient <- AddModuleScore(object_axldc_patient, genes.list = gene_lists)
for (i in 1:3) {
    idx <- which(colnames(object_axldc_patient@meta.data) == paste0("Cluster",i))
    colnames(object_axldc_patient@meta.data)[idx] <- names(gene_lists)[i]
}

saveRDS(object_axldc_patient, file = file.path(OUTDIR, "object_Axldc_patient.Rds"))


q()


