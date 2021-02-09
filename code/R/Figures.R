### FIGURES FOR CARDIOMYOCYTE MANUSCRIPT ###
options(java.parameters = "-Xmx1024m")

library(dplyr)
library(ggplot2)
library(xlsx)
library(tidyr)
library(scales)
library(ggplot2)
library(randomForest)
library(readxl)

library(grid)
library(gridExtra)
library(lattice)
library(ggpubr)
library(cowplot)
library(forcats)

library(ggnewscale)


#### Figure 2 ####

# Figure 2B - comparing ATP yields between iCardio and MitoCore
ATP.yields <- read.xlsx(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/2020-04-10 -- ATP yields.xlsx", 
                        sheetIndex = 1)

Figure2 <- data.frame(ATP.yields) %>% dplyr::filter(MitoCore.ATP.yield != 0) %>% 
  ggplot(aes(x = MitoCore.ATP.yield, y = iCardio.ATP.yield)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, color = "#39568CFF") + 
  xlab("MitoCore ATP yield per unit carbon source") + 
  ylab("iCardio ATP yield per unit carbon source") +
  theme_bw() + 
  theme(strip.background = element_blank())

ggsave('C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/figures/Figure2B.png',Figure2, 
       dpi = 600, width = 17, height = 6, units = 'in')


#### Figure 4 ####
# Figure 4 - examples of differentially expressed tasks
ischemic.random <- readxl::read_xlsx(path= "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/TIDEs/GSE5406_Ischemic_random.xlsx", 
                                   sheet = 1, 
                                   col_names = FALSE)
colnames(ischemic.random)[1:1000] <- paste0('X',c(1:1000))
  
# colnames are 1:1000
# Need to read in some example task data for task names
GSE5406.Ischemic = read.xlsx(file = "data/TIDEs/GSE5406_allGenes_Ischemic.xlsx", sheetIndex = 1, header = 0)
colnames(GSE5406.Ischemic) = c('ID','description','taskScore','significance')

ischemic.random <- cbind(GSE5406.Ischemic %>% select(ID, description), ischemic.random)

# Pick out the TIDEs from GSE5406
GSE5406.Ischemic <- GSE5406.Ischemic %>% 
  mutate(significance = ifelse(abs(significance) >= 0.025, 0, ifelse(significance < 0, -1, 1))) %>% 
  select(ID, description, taskScore, significance)


# Select only a few tasks 
Figure4.legend.data <- ischemic.random %>% 
  filter(ID == "323" | description == "B-hydroxybutyrate uptake") %>% 
  gather(key = "random", value = "value", X1:X1000) %>% 
  select(-random) %>% 
  inner_join(GSE5406.Ischemic, by = c("ID","description")) %>% 
  mutate(description = ifelse(ID == "323", "de novo synthesis of Neu5Ac", "B-hydroxybutyrate uptake"))

Figure4.legend <- ggplot(Figure4.legend.data) + 
  geom_rect(aes(fill = significance), xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf, show.legend = FALSE) +
  scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") + 
  geom_histogram(aes(x = value), color = "black", fill = "black") + 
  geom_vline(aes(xintercept = taskScore), color = "red", linetype = "dashed", size = 1.2) + 
  facet_wrap(~description, scales = "free", nrow = 1) +
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())


# Figure 4 - all metabolic tasks 
Figure4.data <- ischemic.random %>% 
  gather(key = "random", value = "value", X1:X1000) %>% 
  select(-random) %>% 
  inner_join(GSE5406.Ischemic, by = c("ID", "description"))


# Inner join with task categories
task.categories <- read.xlsx(file = "data/TIDEs/TIDEs_categories.xlsx", sheetIndex = 1, header = 1)
task.categories$category <- factor(task.categories$category, 
                                   levels = levels(factor(task.categories$category)))

Figure4.data <- Figure4.data %>% 
  inner_join(task.categories, by = c("description")) %>% 
  mutate(change = ifelse(significance == 0, FALSE, TRUE))

Figure4.data$description <- factor(Figure4.data$description, 
                                   levels = task.categories$description) 

# Creating individual plots for large layout
# generate a function that outputs a plot
generate_plot <- function(data, task_category, ncol){
  figure <- data %>% 
    filter(category == task_category) %>% 
    ggplot() +
    geom_rect(aes(fill = significance), xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf, show.legend = FALSE) +
    scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") + 
    new_scale("fill") +
    geom_histogram(aes(x = value, fill = change)) + 
    scale_fill_manual(values = c("grey","black")) +
    geom_vline(aes(xintercept = taskScore, color = change), linetype = "dashed", size = 0.75) + 
    scale_color_manual(values = c("grey100","red")) +
    facet_wrap(~description, scales = "free", ncol = ncol) +
    theme_bw() + 
    ggtitle(task_category)+
    theme(strip.background = element_blank(), 
          strip.text.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          panel.grid = element_blank(), 
          axis.ticks = element_blank(), 
          plot.title = element_text(hjust = 0.5, face = "bold", size = 15), 
          legend.position = "none")
  return(figure)
}

# Generate individual plots
Figure4A <- generate_plot(Figure4.data, "Amino acid degradation", 10)    #2
Figure4B <- generate_plot(Figure4.data, "Amino acid synthesis", 5)      #3 
Figure4C <- generate_plot(Figure4.data, "Central carbon metabolism", 10) #4
Figure4D <- generate_plot(Figure4.data, "Cofactor synthesis", 5)        #5
Figure4E <- generate_plot(Figure4.data, "Fatty acid degradation", 10)    #6
Figure4F <- generate_plot(Figure4.data, "Fatty acid synthesis", 10)      #7
Figure4G <- generate_plot(Figure4.data, "Ketone body metabolism", 10)    #8
Figure4H <- generate_plot(Figure4.data, "Lipid synthesis", 10)   #9
Figure4I <- generate_plot(Figure4.data, "Nucleotide metabolism", 5)      #10
Figure4J <- generate_plot(Figure4.data, "Signaling metabolism", 10)       #11
Figure4K <- generate_plot(Figure4.data, "Transport", 5)

gridlayout <- read.xlsx(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/figures/Figure4 gridlayout.xlsx", 
                        sheetIndex = 1,
                        header = FALSE) %>% 
  as.matrix(dimnames = NULL)

col_width = 0.06
spacer_width = (1 - 16*col_width)/4
relative_widths = c(col_width, col_width, col_width, col_width, spacer_width, 
                    col_width, col_width, col_width, col_width, spacer_width,
                    col_width, col_width, col_width, col_width, spacer_width,
                    col_width, col_width, spacer_width, col_width, col_width)

col_hgt <- 0.07
spacer_hgt = (1-14*col_hgt)
relative_heights = c(col_hgt, col_hgt,col_hgt,col_hgt,col_hgt,col_hgt,
                     spacer_hgt, 
                     col_hgt,col_hgt,col_hgt,col_hgt,col_hgt,col_hgt,col_hgt,col_hgt)

# Generate a figure of all tasks separated by category
Figure4 <- arrangeGrob(Figure4.legend, Figure4A, Figure4B, Figure4C, 
                       Figure4D, Figure4E, Figure4F, 
                       Figure4G, Figure4H, Figure4I, 
                       Figure4J, Figure4K,
                       ncol = 20, 
                       nrow = 15, 
                       heights = relative_heights,
                       widths = relative_widths,
                       layout_matrix = gridlayout)

ggsave('C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/figures/Figure4.png',Figure4, 
       dpi = 600, width = 17, height = 6, units = 'in')

# Save the name of each task, the description, and the significance (-1,0,1) to Excel spreadsheet for supplemental information
SupplementalFile3 = Figure4.data %>% 
  select(ID, description, category, significance) %>% 
  unique()
SupplementalFile3 = with(SupplementalFile3, SupplementalFile3[order(category, description),])
xlsx::write.xlsx(SupplementalFile3, 
           file = "C:/Users/bvd5nq/Desktop/From Laptop/Papin Lab - Graduate Work/Paper Drafts/Heart Reconstruction/Supplemental Files/Supplemental File 3.xlsx", 
           row.names = FALSE)


#### Figure 5 ####
# Figure 5A - aggregating TIDE results across studies
# Read in all the files than merge into large data table
GSE1869.Ischemic = read.xlsx(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/TIDEs/GSE1869_allGenes_Ischemic.xlsx", sheetIndex = 1, header = 0)
colnames(GSE1869.Ischemic) = c('ID','description','GSE1869.Ischemic.taskScore','GSE1869.Ischemic')

GSE1869.Idiopathic = read.xlsx(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/TIDEs/GSE1869_allGenes_Dilated.xlsx", sheetIndex = 1, header = 0)
colnames(GSE1869.Idiopathic) = c('ID','description','GSE1869.Idiopathic.taskScore','GSE1869.Idiopathic')

GSE5406.Idiopathic = read.xlsx(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/TIDEs/GSE5406_allGenes_Idiopathic.xlsx", sheetIndex = 1, header = 0)
colnames(GSE5406.Idiopathic) = c('ID','description','GSE5406.Idiopathic.taskScore','GSE5406.Idiopathic')

GSE5406.Ischemic = read.xlsx(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/TIDEs/GSE5406_allGenes_Ischemic.xlsx", sheetIndex = 1, header = 0)
colnames(GSE5406.Ischemic) = c('ID','description','GSE5406.Ischemic.taskScore','GSE5406.Ischemic')

GSE57345.Ischemic = read.xlsx(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/TIDEs/GSE57345_allGenes_Ischemic.xlsx", sheetIndex = 1, header = 0)
colnames(GSE57345.Ischemic) = c('ID','description','GSE57345.Ischemic.taskScore','GSE57345.Ischemic')

GSE57345.Idiopathic = read.xlsx(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/TIDEs/GSE57345_allGenes_Idiopathic.xlsx", sheetIndex = 1, header = 0)
colnames(GSE57345.Idiopathic) = c('ID','description','GSE57345.Idiopathic.taskScore','GSE57345.Idiopathic')


all.data <- full_join(GSE1869.Idiopathic %>% dplyr::select(ID, description, GSE1869.Idiopathic), 
                      GSE1869.Ischemic %>% dplyr::select(ID, description, GSE1869.Ischemic), 
                      by = c("ID","description")) %>% 
  full_join(GSE5406.Ischemic %>% dplyr::select(ID, description, GSE5406.Ischemic)) %>% 
  full_join(GSE5406.Idiopathic %>% dplyr::select(ID, description, GSE5406.Idiopathic)) %>%
  full_join(GSE57345.Ischemic %>% dplyr::select(ID, description, GSE57345.Ischemic)) %>% 
  full_join(GSE57345.Idiopathic %>% dplyr::select(ID, description, GSE57345.Idiopathic)) 

colnames(all.data) <- c('ID','description', 
                        'GSE1869 Idiopathic', 'GSE1869 Ischemic',
                        'GSE5406 Ischemic','GSE5406 Idiopathic',
                        'GSE57345 Ischemic','GSE57345 Idiopathic')

all.data <- all.data %>% 
  gather(key = "dataset", value = "significance", 'GSE1869 Idiopathic':'GSE57345 Idiopathic') %>% 
  mutate(significance = ifelse(abs(significance) <= 0.025, ifelse(significance > 0, 1, -1), 0))

# Find the most common metabolic functions associated with changes in transcription
all.data %>% filter(significance != 0) %>% group_by(description) %>% dplyr::count() %>% View()

orig.common.tasks = all.data %>% filter(significance != 0) %>% group_by(ID, description) %>% dplyr::count()

# Save data for supplemental file
save.data <- all.data %>% 
  pivot_wider(names_from = dataset, values_from = significance)
xlsx::write.xlsx(data.frame(save.data), 
           file = "C:/Users/bvd5nq/Desktop/From Laptop/Papin Lab - Graduate Work/Paper Drafts/Heart Reconstruction/Supplemental Files/Supplemental File 4.xlsx", 
           sheetName = "Ischemic vs Idiopathic",
           row.names = FALSE)

# Supplemental figure with all of the TIDEs
# Load reaction subcategories
task.categories <- read.xlsx(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/TIDEs/TIDEs_categories.xlsx", 
                             sheetIndex = 1, 
                             header = 1)
task.categories$category <- factor(task.categories$category, 
                                   levels = levels(factor(task.categories$category)))

# Count the number of tasks from each dataset that are differentially expressed
all.data %>% filter(!grepl('S', ID)) %>%  group_by(dataset) %>% dplyr::summarize(n = sum(significance != 0))

# Change the all.data ID from integer to character
all.data$ID <- as.character(all.data$ID)

# Look to see what TIDEs are shared across all data sets
common.TIDEs <- all.data %>%
  # filter out the subsystem analysis
  filter(!grepl('S', ID)) %>%
  group_by(ID,description) %>% 
  dplyr::summarize(count = sum(significance))

# Pick out the most commonly changed metabolic tasks across studies 
Figure5A.data = all.data %>% 
  filter(ID %in% c(323, 122:149, 260, 'C50', 'C78', 'C85'))

Figure5A.data$description <- factor(Figure5A.data$description, 
                                    levels = task.categories$description) 

Figure5A.data$dataset <- factor(Figure5A.data$dataset,
                                levels = c("GSE5406 Idiopathic","GSE5406 Ischemic",
                                           "GSE1869 Idiopathic","GSE1869 Ischemic",
                                           "GSE57345 Idiopathic","GSE57345 Ischemic"))

# Read in reaction plot names for plotting
# Make sure to have the correct order for plotting
reaction.metadata = read.xlsx(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/reaction_metadata.xlsx", 
                              header = TRUE, sheetIndex = 1)
reaction.metadata$plot.name = factor(reaction.metadata$plot.name, 
                                     levels = reaction.metadata$plot.name)

Figure5A.data = left_join(Figure5A.data, reaction.metadata, 
                                  by = c("ID","description"))

dataset.metadata = data.frame(dataset = c('GSE1869 Ischemic', 'GSE1869 Idiopathic',
                                          'GSE5406 Ischemic','GSE5406 Idiopathic',
                                          'GSE57345 Ischemic','GSE57345 Idiopathic'),
                              GEO = c('GSE1869','GSE1869','GSE5406','GSE5406','GSE57345','GSE57345'), 
                              HF = c('Idiopathic','Ischemic','Ischemic','Idiopathic','Ischemic','Idiopathic'))
dataset.metadata$GEO = factor(dataset.metadata$GEO, levels = c("GSE5406","GSE1869","GSE57345"))

Figure5A.data = Figure5A.data %>% left_join(dataset.metadata, by = c("dataset"))

Figure5A = ggplot(Figure5A.data, aes(x = HF, y = fct_rev(plot.name), fill = significance)) + 
  geom_tile(height = 0.9, width = 0.9) + 
  scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") + 
  theme_bw() + 
  facet_wrap(~GEO) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        #axis.text.y = element_blank(), 
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + 
  scale_y_discrete(position = "right")
Figure5A


# Pick out metabolism that is commonly noted to change in late stage heart failure
# ketone body metabolism, BCAA metabolism, glucose metabolism/glycolysis
Figure5B.data = all.data %>% 
  filter(ID %in% c('C43','C42','C5','C7',
                   56, 'C11', 54, 'C9', 68, 'C10', 
                   309, 310, 
                   181:200))

reaction.metadata = read.xlsx(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/reaction_metadata.xlsx", 
                              header = TRUE, sheetIndex = 2)
reaction.metadata$plot.name = factor(reaction.metadata$plot.name, 
                                     levels = reaction.metadata$plot.name)

Figure5B.data = Figure5B.data %>% left_join(reaction.metadata, by = c("ID","description"))
Figure5B.data$description = factor(Figure5B.data$description, levels = reaction.metadata$description)

dataset.metadata = data.frame(dataset = c('GSE1869 Ischemic', 'GSE1869 Idiopathic',
                                          'GSE5406 Ischemic','GSE5406 Idiopathic',
                                          'GSE57345 Ischemic','GSE57345 Idiopathic'),
                              GEO = c('GSE1869','GSE1869','GSE5406','GSE5406','GSE57345','GSE57345'), 
                              HF = c('Idiopathic','Ischemic','Ischemic','Idiopathic','Ischemic','Idiopathic'))

Figure5B.data = Figure5B.data %>% 
  left_join(dataset.metadata, by = c("dataset"))

Figure5B.data$GEO = factor(Figure5B.data$GEO, 
                                   levels = c("GSE5406","GSE1869","GSE57345"))

Figure5B = ggplot(Figure5B.data, aes(x = HF, y = fct_rev(plot.name), fill = significance)) + 
  geom_tile(height = 0.9, width = 0.9) + 
  scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") + 
  theme_bw() + 
  facet_wrap(~GEO) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        #axis.text.y = element_blank(), 
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + 
  scale_y_discrete(position = "right")
Figure5B

# Save Figure 5
gridlayout = matrix(c(NA,1,NA,2,NA)) %>% t()

space = 0.05
Figure5 = arrangeGrob(Figure5A, Figure5B,
                      ncol = 5, 
                      nrow = 1, 
                      widths = c(0.025, 0.37, 0.115, 0.375, 0.125), 
                      layout_matrix = gridlayout)

# Add text for panel labels
# Add A, B, C labels to the plot
Figure5 <- as_ggplot(Figure5) + 
  draw_plot_label(label = c('A','B'), 
                  size = 16, 
                  x = c(0.0, 0.48),
                  y = c(1,1)) + 
  draw_text(text = c("Fatty acid\nsynthesis",
                     "Ketone body\nmetabolism",
                     "Branched chain\namino acid metabolism",
                     "Glucose metabolism",
                     "Fatty acid oxidation"), 
            size = 10,
            hjust = 0,
            x = c(0.39, 0.88, 0.88, 0.88, 0.88), 
            y = c(0.45, 0.89, 0.76, 0.65, 0.37))

ggsave('C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/figures/Figure5.png', Figure5, 
       dpi = 600, width = 12, height = 6, units = 'in')

########## Supplemental Figure 1 ###############
# Clustering of all DEGs from microarray studies

# Read in DEGs from studies
GSE5406 = read.table(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEdata/GSE5406_DEGs.csv", 
                     header = TRUE, sep = ",") %>% 
  select(GENE.NAME, Ischemic_logFC, Idiopathic_logFC)
colnames(GSE5406) = c('GENE.NAME','GSE5406.Ischemic','GSE5406.Idiopathic')
GSE1869 = read.table(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEdata/GSE1869_DEGs.csv", 
                     header = TRUE, sep = ",") %>% 
  select(ENTREZ_GENE_ID, Ischemic_logFC, Idiopathic_logFC)
colnames(GSE1869) = c('GENE.NAME','GSE1869.Ischemic','GSE1869.Idiopathic')
GSE57345 = read.table(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEdata/GSE57345_DEGs.csv", 
                      header = TRUE, sep = ",") %>% 
  select(GENE.NAME, Ischemic_logFC,  Idiopathic_logFC)
colnames(GSE57345) = c('GENE.NAME','GSE57345.Ischemic','GSE57345.Idiopathic')

all.DEGs = inner_join(GSE5406, GSE57345) %>% 
  inner_join(GSE1869)


# Another approach for heatmaps
library(gplots)
mycol <- colorpanel(1000,"blue","white","red")
colnames(all.DEGs) = c('GENE.NAME','GSE5406 \nIschemic','GSE5406 \nIdiopathic',
                       'GSE57345 \nIschemic','GSE57345 \nIdiopathic',
                       'GSE1869 \nIschemic','GSE1869 \nIdiopathic')

png(filename = 'C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/figures/SuppFigure1A.png')
heatmap.2(as.matrix(all.DEGs[,-1]), 
          #scale="row",
          labRow='', 
          col=mycol, trace="none", density.info="none", 
          #margin=c(8,6), lhei=c(2,10), 
          dendrogram = "both", 
          cexCol = 0.9) 
dev.off()

all.TIDEs = all.data %>% 
  pivot_wider(names_from = "dataset", values_from = "significance")
colnames(all.TIDEs)[3:8] = c('GSE1869 \nIschemic','GSE1869 \nIdiopathic',
                             'GSE5406 \nIschemic','GSE5406 \nIdiopathic',
                             'GSE57345 \nIschemic','GSE57345 \nIdiopathic')

png(filename = 'C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/figures/SuppFigure1B.png')
heatmap.2(as.matrix(all.TIDEs[,c(-1,-2)]), 
          #scale="row",
          labRow='', 
          col=mycol, trace="none", density.info="none", 
          #margin=c(8,6), lhei=c(2,10), 
          dendrogram = "both", 
          cexCol = 0.9) 
dev.off()


#### Supplemental Figure 2A - all TIDEs categorized based on general metabolic categories ####
SuppFigure2A.data <- all.data %>%
  filter(!grepl('S', ID)) %>% 
  left_join(task.categories, by = c("description")) %>% 
  inner_join(dataset.metadata, by = c("dataset"))

SuppFigure2A.data$dataset <- factor(SuppFigure2A.data$dataset,
                                    levels = c("GSE5406 Idiopathic","GSE5406 Ischemic",
                                               "GSE1869 Idiopathic","GSE1869 Ischemic",
                                               "GSE57345 Idiopathic","GSE57345 Ischemic"))

levels(SuppFigure2A.data$dataset) <- gsub(" ", "\n", levels(SuppFigure2A.data$dataset))

# Generate individual panels for groups of tasks
generate_plot <- function(data, task_category){
  figure <- data %>%  
    filter(category == task_category) %>% 
    ggplot(aes(y = HF, x = description, fill = significance)) + 
    geom_tile() + 
    facet_wrap(~GEO) + 
    scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") +
    #scale_fill_gradientn(colours = c("white","blue","red","white")) +
    #scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    theme_bw() + 
    coord_flip() +
    theme(axis.text.x = element_blank(), 
          axis.text.y = element_blank(), 
          legend.position = "none", 
          axis.ticks = element_blank(),
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          strip.text.x = element_blank(), 
          strip.background.x = element_blank(),
          plot.margin = unit(c(0.05,0.05,0.05,0.05), "pt")) 
  return(figure)
}

# Generate individual plots
SuppFigure2A <- SuppFigure2A.data %>% filter(category == 'Amino acid degradation') %>% 
  ggplot(aes(y = HF, x = description, fill = significance)) + 
  geom_tile() + 
  facet_wrap(~GEO) + 
  scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") +
  #scale_fill_gradientn(colours = c("white","blue","red","white")) +
  #scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_bw() + 
  coord_flip() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none", 
        axis.ticks = element_blank(),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.05,0.05,0.05,0.05), "pt")) 
SuppFigure2B <- generate_plot(SuppFigure2A.data, "Amino acid synthesis")      #3 
SuppFigure2C <- generate_plot(SuppFigure2A.data, "Central carbon metabolism") #4
SuppFigure2D <- generate_plot(SuppFigure2A.data, "Cofactor synthesis")        #5
SuppFigure2E <- generate_plot(SuppFigure2A.data, "Fatty acid degradation")    #6
SuppFigure2F <- generate_plot(SuppFigure2A.data, "Fatty acid synthesis")      #7
SuppFigure2G <- generate_plot(SuppFigure2A.data, "Ketone body metabolism")    #8
SuppFigure2H <- generate_plot(SuppFigure2A.data, "Lipid synthesis")           #9
SuppFigure2I <- generate_plot(SuppFigure2A.data, "Nucleotide metabolism")      #10
SuppFigure2J <- generate_plot(SuppFigure2A.data, "Signaling metabolism")       #11
SuppFigure2K <- SuppFigure2A.data %>%  
  filter(category == "Transport") %>% 
  ggplot(aes(y = HF, x = description, fill = significance)) + 
  geom_tile() + 
  facet_wrap(~GEO) + 
  scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") +
  #scale_fill_gradientn(colours = c("white","blue","red","white")) +
  #scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_bw() + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        axis.text.y = element_blank(), 
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(),
        strip.background.x = element_blank(), 
        strip.text.x = element_blank(), 
        plot.margin = unit(c(0.05,0.05,0.05,0.05), "pt")) # order for plot margin: top, right, bottom, left

#### Supplemental Figure 2B - GSEA data ####

# All GSEA data is saved in the GSEA data folder
GSE5406_Idiopathic <- read.delim(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEA data/GSE5406_Idiopathic_Healthy.txt") %>% 
  rbind(read.delim(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEA data/GSE5406_Idiopathic_Idiopathic.txt"))%>% 
  mutate(NOM.p.val = ifelse(NES < 0, -NOM.p.val, NOM.p.val)) %>% 
  select(NAME, NOM.p.val)
colnames(GSE5406_Idiopathic)[2] <- 'GSE5406 Idiopathic'

GSE5406_Ischemic <- read.delim(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEA data/GSE5406_Ischemic_Healthy.txt") %>% 
  rbind(read.delim(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEA data/GSE5406_Ischemic_Ischemic.txt"))%>% 
  mutate(NOM.p.val = ifelse(NES < 0, -NOM.p.val, NOM.p.val)) %>% 
  select(NAME, NOM.p.val)
colnames(GSE5406_Ischemic)[2] <- 'GSE5406 Ischemic'

GSE57345_Idiopathic <- read.delim(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEA data/GSE57345_Idiopathic_Healthy.txt") %>% 
  rbind(read.delim(file  = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEA data/GSE57345_Idiopathic_Idiopathic.txt"))%>% 
  mutate(NOM.p.val = ifelse(NES < 0, -NOM.p.val, NOM.p.val)) %>% 
  select(NAME, NOM.p.val)
colnames(GSE57345_Idiopathic)[2] <- 'GSE57345 Idiopathic'

GSE57345_Ischemic <- read.delim(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEA data/GSE57345_Ischemic_Healthy.txt") %>% 
  rbind(read.delim(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEA data/GSE57345_Ischemic_Ischemic.txt"))%>% 
  mutate(NOM.p.val = ifelse(NES < 0, -NOM.p.val, NOM.p.val)) %>% 
  select(NAME, NOM.p.val)
colnames(GSE57345_Ischemic)[2] <- 'GSE57345 Ischemic'

GSE1869_Idiopathic <- read.delim(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEA data/GSE1869_Idiopathic_Control.txt") %>% 
  rbind(read.delim(file  = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEA data/GSE1869_Idiopathic_Idiopathic.txt"))%>% 
  mutate(NOM.p.val = ifelse(NES < 0, -NOM.p.val, NOM.p.val)) %>% 
  select(NAME, NOM.p.val)
colnames(GSE1869_Idiopathic)[2] <- 'GSE1869 Idiopathic'

GSE1869_Ischemic <- read.delim(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEA data/GSE1869_Ischemic_Control.txt") %>% 
  rbind(read.delim(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEA data/GSE1869_Ischemic_Ischemic.txt"))%>% 
  mutate(NOM.p.val = ifelse(NES < 0, -NOM.p.val, NOM.p.val)) %>% 
  select(NAME, NOM.p.val)
colnames(GSE1869_Ischemic)[2] <- 'GSE1869 Ischemic'


SuppFigure2GSEA.data <- full_join(GSE5406_Idiopathic, GSE5406_Ischemic, by = c("NAME")) %>% 
  full_join(GSE57345_Idiopathic, by = c("NAME")) %>% 
  full_join(GSE57345_Ischemic, by = c("NAME")) %>% 
  full_join(GSE1869_Idiopathic, by = c("NAME")) %>% 
  full_join(GSE1869_Ischemic, by = c("NAME")) %>% 
  gather(key = dataset, value = NOM.p.val, 'GSE5406 Idiopathic':'GSE1869 Ischemic') %>% 
  mutate(NOM.p.val = ifelse(abs(NOM.p.val) >= 0.05, 0, ifelse(NOM.p.val < 0, -1, 1))) %>% 
  mutate(NOM.p.val = ifelse(is.na(NOM.p.val), 0, NOM.p.val))

common.sets <- SuppFigure2GSEA.data %>% group_by(NAME) %>% summarize(counts = sum(NOM.p.val))

SuppFigure2GSEA.data$dataset <- factor(SuppFigure2GSEA.data$dataset,
                                   levels = c("GSE5406 Idiopathic","GSE5406 Ischemic",
                                              "GSE1869 Idiopathic","GSE1869 Ischemic",
                                              "GSE57345 Idiopathic","GSE57345 Ischemic"))
levels(SuppFigure2GSEA.data$dataset) <- gsub(" ", "\n", levels(SuppFigure2GSEA.data$dataset))
SuppFigure2GSEA.data$NAME <- as.factor(SuppFigure2GSEA.data$NAME)
SuppFigure2GSEA.data$NAME <- fct_rev(SuppFigure2GSEA.data$NAME)

no_change <- SuppFigure2GSEA.data %>% 
  group_by(NAME) %>% 
  summarize(n = sum(abs(NOM.p.val))) %>% 
  filter(n != 0)


# Plot the results of all KEGG pathways
SuppFigure4GSEA <- SuppFigure2GSEA.data %>%
  #inner_join(no_change) %>%
  ggplot(aes(y = dataset, x = NAME, fill = NOM.p.val)) +
  geom_tile() +
  scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") +
  #scale_fill_gradient(colours = c("white","blue","red","white")) +
  #scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_bw() +
  #coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        axis.text.y = element_text(hjust = 1, size = 8),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

ggsave('C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/figures/SuppFigure4.png', SuppFigure4GSEA, 
       dpi = 600, width = 15, height = 6, units = 'in')


# Filter to only include metabolic KEGG pathways for Supplemental Figure 2B
KEGG.metabolic.pathways <- read.table(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEA data/KEGG_metabolic_pathways.csv", 
                                      header = TRUE, sep = ",")

# Remove the spaces from the SuppFigure3 dataset variable
levels(SuppFigure2GSEA.data$dataset) <- gsub("\n", " ", levels(SuppFigure2GSEA.data$dataset))
SuppFigure2GSEA.data$dataset <- factor(SuppFigure2GSEA.data$dataset,
                                   levels = c("GSE5406 Idiopathic","GSE5406 Ischemic",
                                              "GSE1869 Idiopathic","GSE1869 Ischemic",
                                              "GSE57345 Idiopathic","GSE57345 Ischemic"))


SuppFigure2GSEA.data <- inner_join(SuppFigure2GSEA.data, KEGG.metabolic.pathways, by = c("NAME" = "KEGG.pathway")) %>% 
  select(-ID) %>% 
  left_join(dataset.metadata)
SuppFigure2GSEA.data$NAME <- factor(SuppFigure2GSEA.data$NAME, levels = KEGG.metabolic.pathways$KEGG.pathway)
SuppFigure2GSEA.data$NAME <- fct_rev(SuppFigure2GSEA.data$NAME)
SuppFigure2GSEA.data$plot.name = factor(SuppFigure2GSEA.data$plot.name, levels = KEGG.metabolic.pathways$plot.name)

# Plot the results of metabolic KEGG pathways
SuppFigure2B.GSEA <- SuppFigure2GSEA.data %>% 
  ggplot(aes(x = HF, y = fct_rev(plot.name), fill = NOM.p.val)) + 
  geom_tile(height = 0.9, width = 0.9) + 
  scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") +
  facet_wrap(~GEO) + 
  theme_bw() + 
  theme(axis.text.y = element_text(size = 7), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + 
  scale_y_discrete(position = "right")


gridlayout <- read.xlsx(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/figures/SuppFigure2A gridlayout.xlsx", 
                        sheetIndex = 1,
                        header = FALSE) %>% 
  as.matrix(dimnames = NULL)

# Create base layout
SuppFigure2 <- arrangeGrob(SuppFigure2A, SuppFigure2B, SuppFigure2C,
                           SuppFigure2D, SuppFigure2E, SuppFigure2F, 
                           SuppFigure2G, SuppFigure2H, SuppFigure2I, 
                           SuppFigure2J, SuppFigure2K, SuppFigure2B.GSEA,
                           ncol = 5, 
                           nrow = 345, 
                           widths = c(0.03, 0.26, 0.2, 0.5, 0.01),
                           layout_matrix = cbind(gridlayout[,2],gridlayout,gridlayout[,2]))


adjust = 0
offset = 0.3
nrow = 345
category_location = c(22, 50, 79, 105, 143, 206.5, 238, 245, 262, 277, 294) + 10
category_location = (nrow - category_location)/nrow
# Next, convert to a ggplot object and add labels for Figure 5A
SuppFigure2 <- as_ggplot(SuppFigure2) +
  draw_text(text = c("Amino acid degradation",
                     "Amino acid synthesis",
                     "Central carbon metabolism",
                     "Cofactor synthesis",
                     "Fatty acid degradation", 
                     "Fatty acid synthesis", 
                     "Ketone body metabolism", 
                     "Lipid membrane synthesis", 
                     "Nucleotide metabolism", 
                     "Signaling metabolism",
                     "Transport"), 
            size = 8,
            hjust = 0,
            x = c(offset+adjust, offset+adjust, offset+adjust, offset+adjust, offset+adjust, offset+adjust, offset+adjust, offset+adjust, offset+adjust, 
                  offset+adjust, offset+adjust), 
            y = category_location) + 
  draw_plot_label(label = c('A','B'), 
                  size = 16, 
                  x = c(0.0, 0.46),
                  y = c(1,1))

ggsave('C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/figures/SuppFigure2.png', SuppFigure2, 
       dpi = 600, width = 9, height = 6, units = 'in')


############### Supplemental Figure 3 - additional heart failure datasets  #####################
# Load in GSE76701, GSE26887, GSE71613, GSE133054 results from TIDEs

# file location: C:\Users\bvd5nq\Documents\R scripts\Cardiomyocyte Model\data\TIDEs\
GSE26887 = read.xlsx(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/TIDEs/GSE26887_allGenes_Failure.xlsx", sheetIndex = 1, header = 0)
colnames(GSE26887) = c('ID','description','GSE26887.taskScore','GSE26887')

GSE71613 = read.xlsx(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/TIDEs/GSE71613_allGenes_Failure.xlsx", sheetIndex = 1, header = 0)
colnames(GSE71613) = c('ID','description','GSE71613.taskScore','GSE71613')

GSE133054 = read.xlsx(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/TIDEs/GSE133054_allGenes_Failure.xlsx", sheetIndex = 1, header = 0)
colnames(GSE133054) = c('ID','description','GSE133054.taskScore','GSE133054 Failure')

GSE141910 = read.xlsx(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/TIDEs/GSE141910_allGenes_Failure.xlsx", sheetIndex = 1, header = 0)
colnames(GSE141910) = c('ID','description','GSE141910.taskScore','GSE141910')

all.failure.results = inner_join(GSE71613 %>% dplyr::select(ID, description, GSE71613),
                                 GSE26887 %>% dplyr::select(ID, description, GSE26887)) %>% 
  inner_join(GSE133054 %>% dplyr::select(ID, description, `GSE133054 Failure`)) %>% 
  inner_join(GSE141910 %>% dplyr::select(ID, description, GSE141910)) %>% 
  pivot_longer(GSE71613:`GSE141910`, 
               names_to = "dataset", 
               values_to = "significance") %>%  
  mutate(label = ifelse(abs(significance) < 0.025, '**', 
                        ifelse(abs(significance) < 0.05, '*', ''))) %>% 
  mutate(significance = ifelse(abs(significance) < 0.05, 
                               ifelse(significance < 0, -1, 1), 0)) %>% 
  # filter out TIDEs that are subsystems
  filter(!grepl('S',ID))

# Number of TIDEs per dataset
all.failure.results %>% filter(significance != 0) %>% group_by(dataset) %>% dplyr::count()

# Write results to Supplemental File 4
# Save data for supplemental file
save.data <- all.failure.results %>% 
  select(-label) %>% 
  pivot_wider(names_from = dataset, values_from = significance)
xlsx::write.xlsx(data.frame(save.data), 
                 file = "C:/Users/bvd5nq/Desktop/From Laptop/Papin Lab - Graduate Work/Paper Drafts/Heart Reconstruction/Supplemental Files/Supplemental File 5.xlsx", 
                 sheetName = "Heart Failure",
                 row.names = FALSE)

common.failure.TIDEs = all.failure.results %>% 
  filter(significance != 0) %>% 
  group_by(ID, description) %>% 
  dplyr::count()
colnames(common.failure.TIDEs)[3] = "failure"
colnames(orig.common.tasks)[3] = "original"

all.common.TIDEs = inner_join(common.failure.TIDEs, orig.common.tasks)


# Generate Supplemental Figure 3
# Panel A is the commonly changed metabolic functions across studies
# Pick out the most commonly changed metabolic tasks across studies 
SuppFigure3A.data = all.failure.results %>% 
  filter(ID %in% c(323, 122:149, 260, 'C50', 'C78', 'C85'))

SuppFigure3A.data$dataset = factor(SuppFigure3A.data$dataset, 
                                   levels = c("GSE26887","GSE71613","GSE133054 Hypertrophic","GSE133054 Failure","GSE141910"))

# Read in reaction plot names for plotting
# Make sure to have the correct order for plotting
reaction.metadata = read.xlsx(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/reaction_metadata.xlsx", 
                              header = TRUE, sheetIndex = 1)
reaction.metadata$plot.name = factor(reaction.metadata$plot.name, 
                                     levels = reaction.metadata$plot.name)

SuppFigure3A.data = left_join(SuppFigure3A.data, reaction.metadata, 
                              by = c("ID","description"))

SuppFigure3A = ggplot(SuppFigure3A.data, aes(x = dataset, y = fct_rev(plot.name), fill = significance)) + 
  geom_tile(height = 0.9, width = 0.9) + 
  scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") + 
  theme_bw() + 
  geom_text(aes(label=label)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        #axis.text.y = element_blank(), 
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + 
  scale_y_discrete(position = "right")
SuppFigure3A


# Panel B is the same metabolic functions from Figure 5B
SuppFigure3B.data = all.failure.results %>% 
  filter(ID %in% c('C43','C42','C5','C7',
                   56, 'C11', 54, 'C9', 68, 'C10', 
                   309, 310, 
                   181:200))

SuppFigure3B.data$dataset = factor(SuppFigure3B.data$dataset, 
                                   levels = c("GSE26887","GSE71613","GSE133054 Hypertrophic","GSE133054 Failure","GSE141910"))

# Read in reaction plot names for plotting
# Make sure to have the correct order for plotting
reaction.metadata = read.xlsx(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/data/reaction_metadata.xlsx", 
                              header = TRUE, sheetIndex = 2)
reaction.metadata$plot.name = factor(reaction.metadata$plot.name, 
                                     levels = reaction.metadata$plot.name)

SuppFigure3B.data = left_join(SuppFigure3B.data, reaction.metadata, 
                              by = c("ID","description"))

SuppFigure3B = ggplot(SuppFigure3B.data, aes(x = dataset, y = fct_rev(plot.name), fill = significance)) + 
  geom_tile(height = 0.9, width = 0.9) + 
  scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") + 
  theme_bw() + 
  geom_text(aes(label=label)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        #axis.text.y = element_blank(), 
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank()) + 
  scale_y_discrete(position = "right")
SuppFigure3B

# Save Figure 5
gridlayout = matrix(c(NA,1,NA,2,NA)) %>% t()

SuppFigure3 = arrangeGrob(SuppFigure3A, SuppFigure3B,
                          ncol = 5, 
                          nrow = 1, 
                          widths = c(0.025, 0.37, 0.115, 0.375, 0.125), 
                          layout_matrix = gridlayout)

# Add text for panel labels
# Add A, B, C labels to the plot
SuppFigure3 <- as_ggplot(SuppFigure3) + 
  draw_plot_label(label = c('A','B'), 
                  size = 16, 
                  x = c(0.0, 0.48),
                  y = c(1,1)) + 
  draw_text(text = c("Fatty acid\nsynthesis",
                     "Ketone body\nmetabolism",
                     "Branched chain\namino acid metabolism",
                     "Glucose metabolism",
                     "Fatty acid oxidation"), 
            size = 10,
            hjust = 0,
            x = c(0.4, 0.88, 0.88, 0.88, 0.88), 
            y = c(0.55, 0.94, 0.81, 0.71, 0.42))

ggsave('C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/figures/SuppFigure3.png', SuppFigure3, 
       dpi = 600, width = 12, height = 6, units = 'in')


# Read in all DEGs for the additional heart failure studies
# Read in DEGs from studies
GSE26887 = read.table(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEdata/GSE26887_DEGs.csv", 
                      header = TRUE, sep = ",") %>% 
  select(GENE.NAME, Failing, logFC)
colnames(GSE26887) = c('GENE.NAME','GSE26887.Failing','GSE26887.logFC')
GSE71613 = read.table(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEdata/GSE71613_DEGs.csv", 
                      header = TRUE, sep = ",") 
colnames(GSE71613) = c('GENE.NAME','GSE71613.Failure','GSE71613.logFC')
GSE133054 = read.table(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEdata/GSE133054_DEGs.csv", 
                       header = TRUE, sep = ",") %>% 
  select(Entrez_ID, failure, failure_logFC)
colnames(GSE133054) = c('GENE.NAME','GSE133054.Failure','GSE133054.logFC')
GSE141910 = read.table(file = "C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/GSEdata/GSE141910_failure_DEGs.csv", 
                       header = TRUE, sep = ",") 
colnames(GSE141910) = c('GENE.NAME','GSE141910.Failure','GSE141910.logFC')

all.DEGs.failure = inner_join(GSE26887, GSE71613) %>% 
  inner_join(GSE141910) %>% 
  inner_join(GSE133054)



#### Supplemental Figure 5 ####
# Analyze TIDE data by reaction subsystem
SuppFigure5.data <- all.data %>%
  filter(grepl('S', ID)) 
SuppFigure5.data$dataset <- factor(SuppFigure5.data$dataset,
                                levels = c("GSE5406 Idiopathic","GSE5406 Ischemic",
                                           "GSE1869 Idiopathic","GSE1869 Ischemic",
                                           "GSE57345 Idiopathic","GSE57345 Ischemic"))
SuppFigure5.data = inner_join(SuppFigure5.data, dataset.metadata)
SuppFigure5.data$description <- as.factor(SuppFigure5.data$description)
SuppFigure5.data$description <- fct_rev(SuppFigure5.data$description)

SuppFigure5 <- SuppFigure5.data %>% 
  ggplot(aes(y = HF, x = description, fill = significance)) + 
  geom_tile() + 
  facet_wrap(~GEO) +
  scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") +
  #scale_fill_gradientn(colours = c("white","blue","red","white")) +
  #scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_bw() + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        axis.text.y = element_text(), 
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        plot.margin = unit(c(0.05,0.05,0.05,0.05), "pt"))

ggsave('C:/Users/bvd5nq/Documents/R scripts/Cardiomyocyte Model/figures/SuppFigure5.png', SuppFigure5, 
       dpi = 600, width = 4.8, height = 6)

