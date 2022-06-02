
## ---------------------------
## The script produces the PCA plot for CD11c pos and neg samples from Axl+DC 
## The script takesthe counting files, creates a DESeq2 object et returns the PCA plot 
##
## Author: Ouardia Ait-Mohamed
##
## ---------------------------


## Load packages ##
# R version 4.1.0
library("DESeq2") # DESeq2_1.32.0 
library("ggplot2") # ggplot2_3.3.5 


##1. Read the featureCounts counts files 

# Path to the folder that contains the counting files

basedir <- "/path_to/preDC"

cntdir <- paste(basedir, "count_files", sep="/")
# the counting files are with .txt extension
pat <- ".txt"
featureCounts.all <- list.files(path = cntdir,
                         pattern = pat,
                         all.files = TRUE,
                         recursive = FALSE,
                         ignore.case = FALSE,
                         include.dirs = FALSE)

# the folder must contain all counting files that should be used 
myfiles <- featureCounts.all
DT <- list()

# read each file as array element of DT and rename the last 2 cols
# we created a list of single sample tables
for (i in 1:length(myfiles) ) {
  infile = paste(cntdir, myfiles[i], sep = "/")
  DT[[myfiles[i]]] <- read.table(infile, header = F, stringsAsFactors = FALSE)
  cnts <- gsub("(.*)_all_counts.txt", "\\1", myfiles[i])
  colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
}

# merge all elements based on protein ID columns
data <- DT[[myfiles[1]]]

# check the data 
head(data)
for (i in 2:length(myfiles)) {
  y <- DT[[myfiles[i]]]
  z <- merge(data, y, by = c("ID"))
  data <- z
}

# change the ID column to become rownames
rownames(data) <- data$ID
data <- data[,-1]

# add total counts per sample
data <- rbind(data, tot.counts=colSums(data))


# create a table like organized in the colData file
rawCountTable <- data

names(rawCountTable) = gsub(pattern = ".txt", replacement = "", x = names(rawCountTable))

# change the table to matrix
cts <- as.matrix(rawCountTable)
cts <- cts[-nrow(cts),] 


##2. Create a DESeq object

# Read the design file
# the design file is a table with two variables : 
#Donor and treatment, whether cells are infected or MOCK

coldata <- read.csv("/path_to/design_file_preDC_CD11C.csv")
coldata <- as.data.frame(coldata)
coldata$Donor <- factor(coldata$Donor)
coldata$Treatment <- factor(coldata$Treatment)



dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = ~ Treatment)



dds$Treatment <- relevel(dds$Treatment, ref = "MOCK")

f_dds <- DESeq(dds)


##3. PCA plot

pca_plot <- function(dds){
  # use the rlog transformation from DESeq2
  rld <- rlog(dds)
  pcaData <- plotPCA(rld, intgroup = c("Treatment", "Donor"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar")) 
  p <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Treatment, shape = Donor)) +
    geom_point(size = 5)+
    scale_color_manual(values = c("CD11cpos_AD8vpx" = "firebrick1", "CD11cneg_AD8vpx" = "blue", "CD11cneg_MOCK" = "gray80", "CD11cpos_MOCK" = "gray30"))+
    scale_shape_manual(values = c("D723" = 19, "D724" = 15, "D730" = 17))+
    geom_hline(yintercept=0, size=0.2)+
    geom_vline(xintercept=0, size=0.2)+
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ylim(-30, 30)+
    xlim(-50,+50)+
    theme(legend.text=element_text(size=14))+
    theme(legend.title=element_text(size=16, face = "bold"))+
    theme(legend.key=element_blank())+
    theme(panel.background = element_rect(fill = "white",
                                          colour = "black",
                                          size = 0.5, linetype = "solid"))
  return(p)
}



























