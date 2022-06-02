
## ---------------------------
## The script produces the heatmap with annotated genes produced 
## The script takes 3 input : - the matrix table containing TPM counts for each gene
##                            - the position of genes in the matrix that have to be annotated
##                            - the labels of genes in the matrix that have to be annotated
##
## Author: Ouardia Ait-Mohamed
##
## ---------------------------

## Load packages ##
# R version 4.1.0
library("ComplexHeatmap") # ComplexHeatmap_2.8.0 
library("circlize") # circlize_0.4.14 


## Arguments ##

if (sys.nframe() == 0){

  rm(list=ls())
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 3) {
    print("Error : Not enough arguments provided")
    print("Usage : Heatmap_with_annotated_genes.R ")
    q()
  } else if (length(args) > 3) {
    print("Error : Too many arguments provided")
    print("Usage : Heatmap_with_annotated_genes.R")  
    q()
  }
  
}


## files that should be provided ##
df <- args[1]
list1 <- args[2]
list2 <- args[3]

## Matrix and vectors transformation ##

mtx <- read.table(df, header=TRUE, sep="\t")

scaledata <- t(scale(t(as.matrix(mtx[,-1]))))
rownames(scaledata) <- mtx$SYMBOL
genpos <- read.delim(list1, header=F, sep = ",")

genpos<-as.vector(unlist(t(as.matrix(genpos))))
labels <- read.delim(list2, header=F, sep = ",")

labels <-as.vector(unlist(t(as.matrix(labels))))

## Heatmap format features ##

pdf(file="/output_path/heatmap_with_annotated_genes.pdf", onefile=FALSE)

col_fun = colorRamp2(c(-4,0,4), c("dodgerblue", "white", "firebrick1"))
col_fun(seq(-4, 4))

ha1 = HeatmapAnnotation(Cell=c(rep("Axl+DC", 12), rep("cDC2", 12)),
Condition=rep(c(rep("mock", 3), rep("AD8", 3), rep("AD8vpx", 3), rep("AD8vpx_AZTNVP", 3)),2),
col = list(Cell=c("Axl+DC"="firebrick4", "cDC2"="gold"),
Condition=c("mock"="gray70", 
            "AD8"="darkseagreen3", 
            "AD8vpx"="green4", 
            "AD8vpx_AZTNVP"="darkorange")),
show_annotation_name = FALSE, simple_anno_size = unit(0.4, "cm"))



ha1@anno_list$Condition@color_mapping@levels <- c("mock", "AD8", "AD8vpx", "AD8vpx_AZTNVP")
ha1@anno_list$Cell@color_mapping@levels <- c("Axl+DC", "cDC2")



anno = anno_mark(at =genpos,
link_width = unit(7, "mm"), padding = 0.4,
 labels = labels, which = "row",
 labels_gp = gpar(fontsize = 9, col=(rep("black",83))))


p <- Heatmap(as.matrix(scaledata), column_title="",
  name = "TPM \n Z-scores",
 col=col_fun,
 row_title = NULL,
 show_column_names = FALSE,
 show_row_dend = FALSE,
  row_split =5,
 row_dend_reorder=TRUE,
 cluster_columns = FALSE,
 row_title_gp=gpar(fontsize =12),
 top_annotation = ha1) + rowAnnotation(mark = anno)

draw(p, padding = unit(c(40,10, 10, 2), "mm"))


for (i in 5:1){
  decorate_heatmap_body("TPM \n Z-scores", {
    grid.lines(c(12/24, 12/24), c(0, 1), 
    gp = gpar(col = "white", lty = 1, lwd=4))
  }, slice=(i))
}


dev.off()





