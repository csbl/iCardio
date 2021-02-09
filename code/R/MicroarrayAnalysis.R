# Written by Bonnie Dougherty, 2019-10-02


library(dplyr)
library(tidyr)
library(VennDiagram)
library(xlsx)

#source("http://www.bioconductor.org/biocLite.R")
#biocLite("GEOquery")
#biocLite("ALL")
#biocLite("Biobase")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("genefilter")

library(GEOquery)
library(limma)
library(ggplot2)


############################## GSE5406 #########################################################################
# Download all data from HF study (210 samples)
# Returns a list of expression sets
gse_data <- getGEO("GSE5406", GSEMatrix = TRUE)

# Returns list of expression sets, only 1 entry
gse_data <- gse_data[[1]]

# Write raw probe data to file for GSEA
data.to.save <- data.frame(NAME = row.names(exprs(gse_data)), DESCRIPTION = "na", exprs(gse_data))
write.table(x = data.to.save, 
            file = "data/GSE5406.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t")

# Boxplot of values to make sure things are normalized correctly
#boxplot(exprs(gse_data))

# Samples are separated by type of HF, fit only as control, ischemic, idiopathic
disease_state <- data.frame(gse_data$title) %>% 
  mutate(Control = grepl("NonFailing", gse_data.title), 
         Ischemic = grepl("Ischemic", gse_data.title), 
         Idiopathic = grepl("Idiopathic",gse_data.title)) %>% 
  mutate(Summary = ifelse(Control == TRUE, "Control", ifelse(Ischemic == TRUE, "Ischemic", "Idiopathic")))

design <- model.matrix(~0 + factor(disease_state$Summary))
colnames(design) <- c("Control","Idiopathic","Ischemic")
contrast.matrix <- makeContrasts(Ischemic - Control, Idiopathic - Control, levels = design)

# Create a cls file based on number of samples
# total samples, total classes, 1
# user visible names for classes, must be in the same order that data appears
# sample labels
fileConn <- file("data/GSE5406_phenotype.cls")
writeLines(c("210 3 1",
             "# Ischemic Idiopathic Control",
             paste(unlist(disease_state$Summary), collapse = " ")), 
           fileConn)
close(fileConn)


# Previous research has indicated that filtering data to top 50% most variable genes increases power
library(genefilter)
gse_data_filter <- varFilter(gse_data)

# Pull out the gene annotation data
# Data only has gene symbol and Entrez Gene ID
gene.data <- fData(gse_data_filter) %>% 
  dplyr::select(ID, ENTREZ_GENE_ID)

# Each probe can have multiple Entrez Gene IDs
s <- strsplit(gene.data$ENTREZ_GENE_ID, split = " /// ")
gene.data <- data.frame(ID = rep(gene.data$ID, sapply(s, length)), ENTREZ_GENE_ID = unlist(s))


# Run limma analysis
gse.data.process <- getEAWP(gse_data_filter)
lmfit <- lmFit(exprs(gse_data_filter), design)
# Run contrasts
lmfit.cont <- contrasts.fit(lmfit, contrast.matrix)
# Apply empirical Bayes
lmfit.cont.ebayes <- eBayes(lmfit.cont)
lmfit.cont.ebayes$genes <- gse.data.process$probes$ID

# Collect logFCs for saving to file
testResults <- topTable(lmfit.cont.ebayes, number = nrow(lmfit.cont.ebayes))

# Adjust for multiple tests, save logFCs
# This data is still with transcripts not grouped to the gene level yet
lmfit.results <- data.frame(ID = gse.data.process$probes$ID, 
                            decideTests(lmfit.cont.ebayes, adjust.method = "fdr", p.value = 0.1))
colnames(lmfit.results) <- c('ID','Ischemic', 'Idiopathic')
# Removing Ischemic.Idiopathic from datasets; onlye 6 decreased, 7 increased probes
lmfit.results <- left_join(lmfit.results, testResults, by = c("ID" = "ProbeID")) %>% 
  dplyr::select('ID', 'Ischemic','Idiopathic', 'Ischemic...Control', 'Idiopathic...Control')
colnames(lmfit.results) <- c('ID', 'Ischemic', 'Idiopathic', 'Ischemic_logFC', 'Idiopathic_logFC')

# Original length
all.results <- full_join(gene.data, lmfit.results)
all.results <- all.results %>% 
  filter(!(is.na(ENTREZ_GENE_ID)))

# Group results by ENTREZ_GENE_ID - remaining is 8242 genes
all.results <- all.results %>% 
  group_by(ENTREZ_GENE_ID) %>% 
  summarize(Ischemic = sum(Ischemic), 
            Ischemic_logFC = ifelse(Ischemic > 0, max(Ischemic_logFC), ifelse(Ischemic < 0, min(Ischemic_logFC), 0)), 
            Idiopathic = sum(Idiopathic), 
            Idiopathic_logFC = ifelse(Idiopathic > 0, max(Idiopathic_logFC), ifelse(Idiopathic < 0, min(Idiopathic_logFC), 0))) %>% 
  mutate(Ischemic = ifelse(Ischemic > 0, 1, ifelse(Ischemic < 0, -1, 0)), 
         Idiopathic = ifelse(Idiopathic > 0, 1, ifelse(Idiopathic < 0, -1, 0)))

all.data <- all.results
colnames(all.data)[1] <- 'GENE.NAME'

all.data %>% group_by(Ischemic) %>% dplyr::count()
all.data %>% group_by(Idiopathic) %>% dplyr::count()

write.table(all.data, file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEdata/GSE5406_DEGs.csv", 
            col.names = TRUE, row.names = FALSE, sep = ",")

# Determine how many DEGs map to the model
genes_in_model <- openxlsx::read.xlsx("C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/20190823 -- genes present in heart model.xlsx", sheet = 1) %>% 
  filter(gene_present == 1) %>% 
  dplyr::select(gene)
genes_in_model <- data.frame(sapply(genes_in_model, function(x) gsub("\'", "", x)))
colnames(genes_in_model) <- c('GENE.NAME')
genes_in_model$GENE.NAME <- (genes_in_model$GENE.NAME)
model.genes.data <- genes_in_model %>% 
  inner_join(all.data, by = c("GENE.NAME"))

model.genes.data %>% group_by(Ischemic == 0) %>% dplyr::count()
model.genes.data %>% group_by(Idiopathic == 0) %>% dplyr::count()

############################################ GSE57345 ##########################################################
# Returns a list of expression sets
gse_data <- getGEO("GSE57345", GSEMatrix = TRUE)

# Returns list of expression sets, only 1 entry
#GPL11532
gse_data <- gse_data[[1]]

# Write raw probe data to file for GSEA
data.to.save <- data.frame(NAME = row.names(exprs(gse_data)), DESCRIPTION = "na", exprs(gse_data))
write.table(x = data.to.save, 
            file = "data/GSE57345.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t")

# Check to see what normalization method was
gse_data@phenoData@data[["data_processing"]][1]

# Boxplot of values to make sure things are normalized correctly
#boxplot(exprs(gse_data))

# Samples are separated by type of HF, fit only as control, ischemic, idiopathic
disease_state <- data.frame(gse_data$title) %>% 
  mutate(Control = grepl("Non-failing", gse_data.title), 
         Ischemic = grepl("Ischemic", gse_data.title), 
         Idiopathic = grepl("Idiopathic", gse_data.title)) %>% 
  mutate(Summary = ifelse(Control == TRUE, "Control", ifelse(Ischemic == TRUE, "Ischemic", "Idiopathic")))

design <- model.matrix(~0 + factor(disease_state$Summary))
colnames(design) <- c("Control","Idiopathic","Ischemic")
contrast.matrix <- makeContrasts(Ischemic - Control, Idiopathic - Control, Ischemic - Idiopathic, levels = design)

# Create a cls file based on number of samples
# total samples, total classes, 1
# user visible names for classes, must be in the same order that data appears
# sample labels
fileConn <- file("data/GSE57345_phenotype.cls")
writeLines(c("313 3 1",
             "# Idiopathic Ischemic Control",
             paste(unlist(disease_state$Summary), collapse = " ")), 
           fileConn)
close(fileConn)



# Previous research has indicated that filtering data to top 50% most variable genes increases power
library(genefilter)
gse_data_filter <- varFilter(gse_data)

library(annotate)
library(hugene10sttranscriptcluster.db)
annodb <- "hugene10sttranscriptcluster.db"
ID <- featureNames(gse_data_filter)
Symbol <- as.character(lookUp(ID,annodb,"SYMBOL"))
Name <- as.character(lookUp(ID,annodb,"GENENAME"))
Entrez <- as.character(lookUp(ID,annodb, "ENTREZID"))

# Can annotate approximatly 55% of the probe sets
length(which(Name == "NA"))
length(ID)

gene.data <- data.frame(ID = as.numeric(ID), ENTREZ_GENE_ID <- Entrez)
colnames(gene.data) <- c('ID','ENTREZ_GENE_ID')
gene.data <- gene.data %>% 
  filter(ENTREZ_GENE_ID != "NA")


# Run limma analysis
gse.data.process <- getEAWP(gse_data_filter)
lmfit <- lmFit(exprs(gse_data_filter), design)
# Run contrasts
lmfit.cont <- contrasts.fit(lmfit, contrast.matrix)
# Apply empirical Bayes
lmfit.cont.ebayes <- eBayes(lmfit.cont)
lmfit.cont.ebayes$genes <- gse.data.process$probes$ID

# Collect logFCs for saving to file
testResults <- topTable(lmfit.cont.ebayes, number = nrow(lmfit.cont.ebayes))

# Adjust for multiple tests, save logFCs
# This data is still with transcripts not grouped to the gene level yet
lmfit.results <- data.frame(ID = gse.data.process$probes$ID, 
                            decideTests(lmfit.cont.ebayes, adjust.method = "fdr", p.value = 0.1))
colnames(lmfit.results) <- c('ID','Ischemic', 'Idiopathic', 'Ischemic.Idiopathic')
# Removed Ischemic...Idiopathic comparison from analysis: only 15 decreased and 30 increased genes across probes
lmfit.results <- left_join(lmfit.results, testResults, by = c("ID" = "ProbeID")) %>% 
  dplyr::select('ID', 'Ischemic','Idiopathic', 'Ischemic...Control', 'Idiopathic...Control')
colnames(lmfit.results) <- c('ID', 'Ischemic', 'Idiopathic', 'Ischemic_logFC', 'Idiopathic_logFC')

# Original length
all.results <- full_join(gene.data, lmfit.results)
all.results <- all.results %>% 
  filter(!(is.na(ENTREZ_GENE_ID)))

# Group results by ENTREZ_GENE_ID - remaining is 8242 genes
all.results <- all.results %>% 
  group_by(ENTREZ_GENE_ID) %>% 
  summarize(Ischemic = sum(Ischemic), 
            Ischemic_logFC = ifelse(Ischemic > 0, max(Ischemic_logFC), ifelse(Ischemic < 0, min(Ischemic_logFC), 0)), 
            Idiopathic = sum(Idiopathic), 
            Idiopathic_logFC = ifelse(Idiopathic > 0, max(Idiopathic_logFC), ifelse(Idiopathic < 0, min(Idiopathic_logFC), 0))) %>% 
  mutate(Ischemic = ifelse(Ischemic > 0, 1, ifelse(Ischemic < 0, -1, 0)), 
         Idiopathic = ifelse(Idiopathic > 0, 1, ifelse(Idiopathic < 0, -1, 0)))

all.data <- all.results
colnames(all.data)[1] <- 'GENE.NAME'

all.data %>% group_by(Ischemic) %>% dplyr::count()
all.data %>% group_by(Idiopathic) %>% dplyr::count()

write.table(all.data, file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEdata/GSE57345_DEGs.csv", 
            col.names = TRUE, row.names = FALSE, sep = ",")

# Determine how many DEGs map to the model
genes_in_model <- openxlsx::read.xlsx("C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/20190823 -- genes present in heart model.xlsx", sheet = 1) %>% 
  filter(gene_present == 1) %>% 
  dplyr::select(gene)
genes_in_model <- data.frame(sapply(genes_in_model, function(x) gsub("\'", "", x)))
colnames(genes_in_model) <- c('GENE.NAME')
genes_in_model$GENE.NAME <- (genes_in_model$GENE.NAME)
model.genes.data <- genes_in_model %>% 
  inner_join(all.data, by = c("GENE.NAME"))

model.genes.data %>% group_by(Ischemic == 0) %>% dplyr::count()
model.genes.data %>% group_by(Idiopathic == 0) %>% dplyr::count()

############################################ GSE1869 ##########################################################
# Returns a list of expression sets
gse_data <- getGEO("GSE1869", GSEMatrix = TRUE)

# Returns list of expression sets
gse_data <- gse_data[[1]]

# Write raw probe data to file for GSEA
data.to.save <- data.frame(NAME = row.names(exprs(gse_data)), DESCRIPTION = "na", exprs(gse_data))
write.table(x = data.to.save, 
            file = "data/GSE1869.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t")

# Check to see what normalization method was
gse_data@phenoData@data[["data_processing"]]

# Boxplot of values to make sure things are normalized correctly
#boxplot((exprs(gse_data)))

# Samples are separated by type of HF, fit only as control, ischemic, idiopathic
disease_state <- data.frame(gse_data$description) %>% 
  mutate(Control = grepl("donor heart", gse_data.description), 
         Ischemic = grepl("Pre-LVAD, ischemic", gse_data.description), 
         Idiopathic = grepl("Pre-LVAD, nonischemic",gse_data.description)) %>% 
  mutate(Summary = ifelse(Control == TRUE, "Control", ifelse(Ischemic == TRUE, "Ischemic", 
                                                             ifelse(Idiopathic == TRUE, "Idiopathic", "Excluded"))))

design <- model.matrix(~0 + factor(disease_state$Summary))
colnames(design) <- c("Control","Excluded","Idiopathic","Ischemic")
contrast.matrix <- makeContrasts(Ischemic - Control, Idiopathic - Control, Ischemic - Idiopathic, levels = design)

# Create a cls file based on number of samples
# total samples, total classes, 1
# user visible names for classes, must be in the same order that data appears
# sample labels
fileConn <- file("data/GSE1869_phenotype.cls")
writeLines(c("37 4 1",
             "# Excluded Ischemic Idiopathic Control",
             paste(unlist(disease_state$Summary), collapse = " ")), 
           fileConn)
close(fileConn)

# Previous research has indicated that filtering data to top 50% most variable genes increases power
library(genefilter)
gse_data_filter <- varFilter(gse_data)

# Pull out the gene annotation data
# Data only has gene symbol and Entrez Gene ID
gene.data <- fData(gse_data_filter) %>% 
  dplyr::select(ID, ENTREZ_GENE_ID)

# Each probe can have multiple Entrez Gene IDs
s <- strsplit(gene.data$ENTREZ_GENE_ID, split = " /// ")
gene.data <- data.frame(ID = rep(gene.data$ID, sapply(s, length)), ENTREZ_GENE_ID = unlist(s))


# Run limma analysis
gse.data.process <- getEAWP(gse_data_filter)
lmfit <- lmFit(exprs(gse_data_filter), design)
# Run contrasts
lmfit.cont <- contrasts.fit(lmfit, contrast.matrix)
# Apply empirical Bayes
lmfit.cont.ebayes <- eBayes(lmfit.cont)
lmfit.cont.ebayes$genes <- gse.data.process$probes$ID

# Collect logFCs for saving to file
testResults <- topTable(lmfit.cont.ebayes, number = nrow(lmfit.cont.ebayes))

# Adjust for multiple tests, save logFCs
# This data is still with transcripts not grouped to the gene level yet
lmfit.results <- data.frame(ID = gse.data.process$probes$ID, 
                            decideTests(lmfit.cont.ebayes, adjust.method = "fdr", p.value = 0.1))
colnames(lmfit.results) <- c('ID','Ischemic', 'Idiopathic', 'Ischemic.Idiopathic')
lmfit.results <- left_join(lmfit.results, testResults, by = c("ID" = "ProbeID")) %>% 
  dplyr::select('ID', 'Ischemic','Idiopathic', 'Ischemic...Control', 'Idiopathic...Control')
colnames(lmfit.results) <- c('ID', 'Ischemic', 'Idiopathic', 'Ischemic_logFC', 'Idiopathic_logFC')

# Original length
all.results <- full_join(gene.data, lmfit.results)
all.results <- all.results %>% 
  filter(!(is.na(ENTREZ_GENE_ID)))

# Group results by ENTREZ_GENE_ID - remaining is 8242 genes
all.results <- all.results %>% 
  group_by(ENTREZ_GENE_ID) %>% 
  summarize(Ischemic = sum(Ischemic), 
            Ischemic_logFC = ifelse(Ischemic > 0, max(Ischemic_logFC), ifelse(Ischemic < 0, min(Ischemic_logFC), 0)), 
            Idiopathic = sum(Idiopathic), 
            Idiopathic_logFC = ifelse(Idiopathic > 0, max(Idiopathic_logFC), ifelse(Idiopathic < 0, min(Idiopathic_logFC), 0))) %>% 
  mutate(Ischemic = ifelse(Ischemic > 0, 1, ifelse(Ischemic < 0, -1, 0)), 
         Idiopathic = ifelse(Idiopathic > 0, 1, ifelse(Idiopathic < 0, -1, 0)))

all.results %>% group_by(Ischemic) %>% dplyr::count()
all.results %>% group_by(Idiopathic) %>% dplyr::count()

write.table(all.results, file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEdata/GSE1869_DEGs.csv", 
            col.names = TRUE, row.names = FALSE, sep = ",")

# Determine how many DEGs map to the model
genes_in_model <- openxlsx::read.xlsx("C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/20190823 -- genes present in heart model.xlsx", sheet = 1) %>% 
  filter(gene_present == 1) %>% 
  dplyr::select(gene)
genes_in_model <- data.frame(sapply(genes_in_model, function(x) gsub("\'", "", x)))
colnames(genes_in_model) <- c('GENE.NAME')
genes_in_model$GENE.NAME <- (genes_in_model$GENE.NAME)
model.genes.data <- genes_in_model %>% 
  inner_join(all.results, by = c("GENE.NAME"="ENTREZ_GENE_ID"))

model.genes.data %>% group_by(Ischemic == 0) %>% dplyr::count()
model.genes.data %>% group_by(Idiopathic == 0) %>% dplyr::count()

#### Additional dataset included for validation ####
#### GSE26887 ####
# Returns a list of expression sets
gse_data <- getGEO("GSE26887", GSEMatrix = TRUE)

# Returns list of expression sets, only 1 entry
gse_data <- gse_data[[1]]

# Write raw probe data to file for GSEA
data.to.save <- data.frame(NAME = row.names(exprs(gse_data)), DESCRIPTION = "na", exprs(gse_data))
write.table(x = data.to.save, 
            file = "C:/Users/bvd5nq/Documents/GSEA/GEOdata/GSE26887.txt",
            row.names = FALSE, col.names = TRUE, sep = "\t")

# Check to see what normalization method was
gse_data@phenoData@data[["data_processing"]][1]

# Boxplot of values to make sure things are normalized correctly
boxplot((exprs(gse_data)))

# Samples are separated by type of HF, fit only as control, ischemic, idiopathic
disease_state <- data.frame(gse_data$title) %>% 
  mutate(Control = grepl("control", gse_data.title), 
         Failure = grepl("non diabetic", gse_data.title), 
         Idiopathic = grepl("_diabetic",gse_data.title)) %>% 
  mutate(Summary = ifelse(Control == TRUE, "Control", ifelse(Failure == TRUE, "Failure", "Diabetic")))

design <- model.matrix(~0 + factor(disease_state$Summary))
colnames(design) <- c("Control","Failure","Diabetic")
contrast.matrix <- makeContrasts(Failure - Control, levels = design)

# Create a cls file based on number of samples
# total samples, total classes, 1
# user visible names for classes, must be in the same order that data appears
# sample labels
fileConn <- file("C:/Users/bvd5nq/Documents/GSEA/GEOdata/GSE26887_phenotype.cls")
writeLines(c("24 3 1",
             "# Control Diabetic Failure",
             paste(unlist(disease_state$Summary), collapse = " ")), 
           fileConn)
close(fileConn)


# Previous research has indicated that filtering data to top 50% most variable genes increases power
library(genefilter)
gse_data_filter <- varFilter(gse_data)

library(annotate)
library(hugene10sttranscriptcluster.db)
annodb <- "hugene10sttranscriptcluster.db"
ID <- featureNames(gse_data_filter)
Symbol <- as.character(lookUp(ID,annodb,"SYMBOL"))
Name <- as.character(lookUp(ID,annodb,"GENENAME"))
Entrez <- as.character(lookUp(ID,annodb, "ENTREZID"))

# Can annotate approximatly 55% of the probe sets
length(which(Name == "NA"))
length(ID)

gene.data <- data.frame(ID = as.numeric(ID), ENTREZ_GENE_ID <- Entrez)
colnames(gene.data) <- c('ID','ENTREZ_GENE_ID')
gene.data <- gene.data %>% 
  filter(ENTREZ_GENE_ID != "NA")

# Run limma analysis
gse.data.process <- getEAWP(gse_data_filter)
lmfit <- lmFit(exprs(gse_data_filter), design)
# Run contrasts
lmfit.cont <- contrasts.fit(lmfit, contrast.matrix)
# Apply empirical Bayes
lmfit.cont.ebayes <- eBayes(lmfit.cont)
lmfit.cont.ebayes$genes <- gse.data.process$probes$ID

# Collect logFCs for saving to file
testResults <- topTable(lmfit.cont.ebayes, number = nrow(lmfit.cont.ebayes)) %>% 
  filter(!is.na(ID))

# Adjust for multiple tests, save logFCs
# This data is still with transcripts not grouped to the gene level yet
lmfit.results <- data.frame(ID = gse.data.process$probes$ID, decideTests(lmfit.cont.ebayes, adjust.method = "fdr", p.value = 0.1))
colnames(lmfit.results) <- c('ID','Failing')
lmfit.results <- left_join(lmfit.results, testResults, by = c("ID")) %>% 
  dplyr::select('ID', 'Failing','logFC')
colnames(lmfit.results) <- c('ID', 'Failing', 'logFC')

# Original length
all.results <- full_join(gene.data, lmfit.results)
colnames(all.results) <- c("ID","ENTREZ_GENE_ID","Failing",'logFC')
all.results <- all.results %>% 
  filter(!(is.na(ENTREZ_GENE_ID)))

# Group results by ENTREZ_GENE_ID - remaining is 8242 genes
all.results <- all.results %>% 
  group_by(ENTREZ_GENE_ID) %>% 
  summarize(Failing = sum(Failing), 
            logFC = ifelse(Failing > 0, max(logFC), ifelse(Failing < 0, min(logFC), 0))) %>% 
  mutate(Failing = ifelse(Failing > 0, 1, ifelse(Failing < 0, -1, 0)))

all.results %>% group_by(Failing == 0) %>% dplyr::count()

# Filter out the genes that are in the human reconstruction
library(openxlsx)

# new iCardio model has 2173 genes
genes_in_model <- openxlsx::read.xlsx("C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/20190823 -- genes present in heart model.xlsx", sheet = 1) %>% 
  filter(gene_present == 1) %>% 
  dplyr::select(gene)
genes_in_model <- data.frame(sapply(genes_in_model, function(x) gsub("\'", "", x)))
colnames(genes_in_model) <- c('GENE.NAME')
genes_in_model$GENE.NAME <- (genes_in_model$GENE.NAME)
model.genes.data <- genes_in_model %>% 
  inner_join(all.results, by = c("GENE.NAME"="ENTREZ_GENE_ID"))

model.genes.data %>% group_by(Failing == 0) %>% dplyr::count()

all.data <- all.results
colnames(all.data)[1] <- 'GENE.NAME'



write.table(all.data, file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEdata/GSE26887_DEGs.csv", 
            col.names = TRUE, row.names = FALSE, sep = ",")

write.table(model.genes.data, file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEdata/GSE26887_model_DEGs.csv", 
            col.names = TRUE, row.names = FALSE, sep = ",")





