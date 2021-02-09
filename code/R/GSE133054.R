# Written: 7/11/2019 (BVD)

# Load applicable libraries
library(GEOquery)
library(BiocManager)
library(xlsx)
library(tidyr)
library(dplyr)

library(tximport)
library(biomaRt)
library(DESeq2)
library(vsn)
library(rhdf5)

library(ggplot2)
library(tibble)

#Source files from bioconductor
#BiocManager::install("tximport") 
#BiocManager::install("biomaRt")
#BiocManager::install("DESeq2")
#BiocManager::install('vsn')
#BiocManager::install('rhdf5)


# Create file names to get abundance files from Kallisto 
dir <- "F:/Transcriptomics/kallisto/outputs/GSE133054"

# Check for all files in the directory to ensure they are present 
list.files(dir)

# Read in meta data table
metadata <- read.xlsx(file.path("F:/Transcriptomics/kallisto/outputs/GSE133054/metadata.xlsx"), header = TRUE, sheetIndex =1) %>% 
  filter(!is.na(Folder.for.kallisto))

# Read in each file based off the format, and label sample numbers
files <- file.path(dir,metadata$Folder.for.kallisto,"abundance.tsv")
names(files) <- paste0("sample",1:dim(metadata)[1])

# Use biomArt to gather transcript and gene IDs
mart  <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
t2g6 <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name", "entrezgene_id", "transcript_version","description","phenotype_description"), mart = mart)
t2g6 <- within(t2g6,target_id <- paste(ensembl_transcript_id, transcript_version,sep = "."))
t2g6 <- dplyr::rename(t2g6, transcript_id = ensembl_transcript_id,ens_gene = ensembl_gene_id, ext_gene = external_gene_name, etz_gene = entrezgene_id, descript = description, phenotype = phenotype_description)
tx2gene <- t2g6[,c(8,4)] %>% filter(!is.na(etz_gene)) %>% unique()
dds2model <- t2g6[,c(2,3,4)]
dds2model <- dds2model[!duplicated(dds2model),]

###### Use TxImport to summarize transcript changes to the gene level #####################################
txi <- tximport::tximport(files,type = "kallisto", txOut = FALSE, tx2gene = tx2gene)

# Label all conditions in the matrix in order of appearance in the each of the files variables
condition = metadata$condition

# Change condition vector from character vector to factor vector
condition = factor(condition)

# Convert conditions to data frame for maniupaltion 
sampleTable <- data.frame(condition)

# Create row names based on the genes in the same order
rownames(sampleTable) <- colnames(txi$counts)

###################################### Use DESeq2 to analyze differential gene expression ########################################
# Generate a DESeq2 dataset
# Import TxImport Results into DESeq2 structure
dds <- DESeq2::DESeqDataSetFromTximport(txi,sampleTable,~condition)

# Pre-filter data to remove low counts
nrow(dds)
# determine the number of samples
metadata %>% group_by(condition) %>% dplyr::count()

keep <- rowSums(counts(dds) >= 5) >= 7
dds <- dds[keep,]
nrow(dds)

# vst transformation is recommended for n > 30 (n = 65)
# blind = FALSE means that the experimental design is not used directly in the transformation
# but is used in estimating the global amount of variability in the counts
vsd <- vst(dds, blind = FALSE)


# Making comparisons between different normalizations for heteroskedasticity
# account for sequencing depth normalization, already done with vst and rld
dds <- estimateSizeFactors(dds)

# PCA analysis
# Use the vsd data rather than original counts data
PCAdata <- plotPCA(vsd, returnData = TRUE)
percentVar <- round(100 * attr(PCAdata, "percentVar"))

# The results seen in the heatmap are also shown here
ggplot(PCAdata, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")

#### Running differential expression analysis with DESeq2
dds <- DESeq(dds)

# Create results data frames
failure_DEGs <- data.frame(DESeq2::results(dds, 
                                           contrast = c("condition", "heart_failure", "healthy"), 
                                           alpha = 0.1, 
                                           pAdjustMethod = "BH", 
                                           independentFiltering = TRUE)) %>% 
  rownames_to_column("Entrez_ID") %>% 
  mutate(padj = ifelse(is.na(padj), 1, padj))

failure_DEGs %>% group_by(padj < 0.1) %>% dplyr::count()

# Write results to files for TIDEs analysis
all.results = failure_DEGs %>% 
  mutate(failure_logFC = ifelse(padj < 0.1, log2FoldChange, 0), 
         failure = ifelse(padj < 0.1, 
                          ifelse(log2FoldChange < 0, -1, 1), 0)) %>% 
  dplyr::select(Entrez_ID, failure, failure_logFC)

write.table(all.results, file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEdata/GSE133054_failure_DEGs.csv", 
            col.names = TRUE, row.names = FALSE, sep = ",")

# Determine how many DEGs map to the model
genes_in_model <- openxlsx::read.xlsx("C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/20190823 -- genes present in heart model.xlsx", sheet = 1) %>% 
  filter(gene_present == 1) %>% 
  dplyr::select(gene)
genes_in_model <- data.frame(sapply(genes_in_model, function(x) gsub("\'", "", x)))
colnames(genes_in_model) <- c('GENE.NAME')
genes_in_model$GENE.NAME <- (genes_in_model$GENE.NAME)
model.genes.data <- genes_in_model %>% 
  inner_join(all.results, by = c("GENE.NAME"="Entrez_ID"))

model.genes.data %>% group_by(failure) %>% dplyr::count()

write.table(model.genes.data, file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEdata/GSE133054_failure_model_DEGs.csv", 
            col.names = TRUE, row.names = FALSE, sep = ",")
