### FIGURES FOR iCardio MANUSCRIPT ###
# Written by Bonnie Dougherty, 2020-04-08

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


#### Figure 2 ####

# Figure 2B - comparing ATP yields between iCardio and MitoCore
ATP.yields <- read.xlsx(file = "data/ATP yields.xlsx", 
                        sheetIndex = 1)

Figure.2 <- ATP.yields %>% filter(MitoCore.ATP.yield != 0) %>% 
  ggplot(aes(x = MitoCore.ATP.yield, y = iCardio.ATP.yield)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE, color = "#39568CFF") + 
  xlab("MitoCore ATP yield per unit carbon source") + 
  ylab("iCardio ATP yield per unit carbon source") +
  theme_bw() + 
  theme(strip.background = element_blank())




#### Figure 4 ####
# Figure 4 - examples of differentially expressed tasks
ischemic.random <- readxl::read_xlsx(path= "data/GSE5406_Ischemic_random.xlsx", 
                                   sheet = 1, 
                                   col_names = FALSE)
colnames(ischemic.random)[1:1000] <- paste0('X',c(1:1000))
  
# colnames are 1:1000
# Need to read in some example task data for task names
GSE5406.Ischemic = read.xlsx(file = "data/GSE5406_allGenes_Ischemic.xlsx", sheetIndex = 1, header = 0)
colnames(GSE5406.Ischemic) = c('ID','description','taskScore','significance')

ischemic.random <- cbind(GSE5406.Ischemic %>% select(ID, description), ischemic.random)

# Pick out the TIDEs from GSE5406
GSE5406.Ischemic <- GSE5406.Ischemic %>% 
  mutate(significance = ifelse(abs(significance) >= 0.025, 0, ifelse(significance < 0, -1, 1))) %>% 
  select(ID, description, taskScore, significance)


# Select only a few tasks 
Figure4.legend.data <- ischemic.random %>% 
  filter(ID == "104" | description == "Acetoacetate uptake") %>% 
  gather(key = "random", value = "value", X1:X1000) %>% 
  select(-random) %>% 
  inner_join(GSE5406.Ischemic, by = c("ID","description")) %>% 
  mutate(description = ifelse(ID == "104", "PC de novo synthesis", "Acetoacetate uptake"))

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
task.categories <- read.xlsx(file = "data/TIDEs_categories.xlsx", sheetIndex = 1, header = 1)
task.categories$category <- factor(task.categories$category, 
                                   levels = levels(factor(task.categories$category)))

Figure4.data <- Figure4.data %>% 
  inner_join(task.categories, by = c("description"))

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
    geom_histogram(aes(x = value), color = "black", fill = "black") + 
    geom_vline(aes(xintercept = taskScore), color = "red", linetype = "dashed", size = 0.9) + 
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
          plot.title = element_text(hjust = 0.5, face = "bold", size = 15))
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
Figure4H <- generate_plot(Figure4.data, "Lipid membrane synthesis", 10)   #9
Figure4I <- generate_plot(Figure4.data, "Nucleotide metabolism", 5)      #10
Figure4J <- generate_plot(Figure4.data, "Signaling metabolism", 10)       #11
Figure4K <- generate_plot(Figure4.data, "Transport", 5)

gridlayout <- read.xlsx(file = "data/Figure4 gridlayout.xlsx", 
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

ggsave('figures/Figure4.png',Figure4, 
       dpi = 600, width = 17, height = 6, units = 'in')

# Count the number of TIDEs
Figure4.data %>% select(description, significance) %>% unique() %>% group_by(significance) %>% dplyr::count() %>% View()
# total of 306 tasks tested
# 94 TIDEs



#### Figure 5 ####
# Figure 5A - aggregating TIDE results across studies
# Read in all the files than merge into large data table
GSE1869.Ischemic = read.xlsx(file = "data/GSE1869_allGenes_Ischemic.xlsx", sheetIndex = 1, header = 0)
colnames(GSE1869.Ischemic) = c('ID','description','GSE1869.Ischemic.taskScore','GSE1869.Ischemic')

GSE1869.Idiopathic = read.xlsx(file = "data/GSE1869_allGenes_Dilated.xlsx", sheetIndex = 1, header = 0)
colnames(GSE1869.Idiopathic) = c('ID','description','GSE1869.Idiopathic.taskScore','GSE1869.Idiopathic')

GSE5406.Idiopathic = read.xlsx(file = "data/GSE5406_allGenes_Idiopathic.xlsx", sheetIndex = 1, header = 0)
colnames(GSE5406.Idiopathic) = c('ID','description','GSE5406.Idiopathic.taskScore','GSE5406.Idiopathic')

GSE5406.Ischemic = read.xlsx(file = "data/GSE5406_allGenes_Ischemic.xlsx", sheetIndex = 1, header = 0)
colnames(GSE5406.Ischemic) = c('ID','description','GSE5406.Ischemic.taskScore','GSE5406.Ischemic')

GSE57345.Ischemic = read.xlsx(file = "data/GSE57345_allGenes_Ischemic.xlsx", sheetIndex = 1, header = 0)
colnames(GSE57345.Ischemic) = c('ID','description','GSE57345.Ischemic.taskScore','GSE57345.Ischemic')

GSE57345.Idiopathic = read.xlsx(file = "data/GSE57345_allGenes_Idiopathic.xlsx", sheetIndex = 1, header = 0)
colnames(GSE57345.Idiopathic) = c('ID','description','GSE57345.Idiopathic.taskScore','GSE57345.Idiopathic')


all.data <- full_join(GSE1869.Idiopathic %>% select(ID, description, GSE1869.Idiopathic), 
                      GSE1869.Ischemic %>% select(ID, description, GSE1869.Ischemic), 
                      by = c("ID","description")) %>% 
  full_join(GSE5406.Ischemic %>% select(ID, description, GSE5406.Ischemic)) %>% 
  full_join(GSE5406.Idiopathic %>% select(ID, description, GSE5406.Idiopathic)) %>%
  full_join(GSE57345.Ischemic %>% select(ID, description, GSE57345.Ischemic)) %>% 
  full_join(GSE57345.Idiopathic %>% select(ID, description, GSE57345.Idiopathic)) 

colnames(all.data) <- c('ID','description', 'GSE1869 Idiopathic', 'GSE1869 Ischemic',
                        'GSE5406 Ischemic','GSE5406 Idiopathic',
                        'GSE57345 Ischemic','GSE57345 Idiopathic')

all.data <- all.data %>% 
  gather(key = "dataset", value = "significance", 'GSE1869 Idiopathic':'GSE57345 Idiopathic') %>% 
  mutate(significance = ifelse(abs(significance) >= 0.025, 0, ifelse(significance > 0, 1, -1)))

# Save data for supplemental file
save.data <- all.data %>% 
  pivot_wider(names_from = dataset, values_from = significance)
write.xlsx(save.data, 
           file = "data/Supplemental File 3.xlsx", 
           row.names = FALSE)

# Count the number of tasks from each dataset that are differentially expressed
all.data %>% group_by(dataset) %>% dplyr::summarize(n = sum(significance != 0))

# Change the all.data ID from integer to character
all.data$ID <- as.character(all.data$ID)

# Load reaction subcategories
task.categories <- read.xlsx(file = "data/TIDEs_categories.xlsx", 
                             sheetIndex = 1, 
                             header = 1)
task.categories$category <- factor(task.categories$category, 
                                   levels = levels(factor(task.categories$category)))
# Produce Figure 5
Figure5A.data <- all.data %>%
  filter(!grepl('S', ID)) %>% 
  left_join(task.categories, by = c("description"))

Figure5A.data$description <- factor(Figure5A.data$description, 
                                   levels = task.categories$description) 

Figure5A.data$dataset <- factor(Figure5A.data$dataset,
                               levels = c("GSE5406 Idiopathic","GSE5406 Ischemic",
                                          "GSE1869 Idiopathic","GSE1869 Ischemic",
                                          "GSE57345 Idiopathic","GSE57345 Ischemic"))
levels(Figure5A.data$dataset) <- gsub(" ", "\n", levels(Figure5A.data$dataset))

# Count the number of tasks without subsystems
Figure5A.data %>% group_by(dataset) %>% dplyr::summarize(n = sum(significance != 0))


# 275 x 600
Figure5A <- Figure5A.data %>%  
  ggplot(aes(y = dataset, x = description, fill = significance)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") +
  #scale_fill_gradientn(colours = c("white","blue","red","white")) +
  #scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_bw() + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text.y = element_blank(), 
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

# Generate individual panels for groups of tasks
generate_plot <- function(data, task_category){
  figure <- data %>%  
    filter(category == task_category) %>% 
    ggplot(aes(y = dataset, x = description, fill = significance)) + 
    geom_tile() + 
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
  return(figure)
}

# Generate individual plots
Figure5A <- generate_plot(Figure5A.data, "Amino acid degradation")    #2
Figure5B <- generate_plot(Figure5A.data, "Amino acid synthesis")      #3 
Figure5C <- generate_plot(Figure5A.data, "Central carbon metabolism") #4
Figure5D <- generate_plot(Figure5A.data, "Cofactor synthesis")        #5
Figure5E <- generate_plot(Figure5A.data, "Fatty acid degradation")    #6
Figure5F <- generate_plot(Figure5A.data, "Fatty acid synthesis")      #7
Figure5G <- generate_plot(Figure5A.data, "Ketone body metabolism")    #8
Figure5H <- generate_plot(Figure5A.data, "Lipid membrane synthesis")   #9
Figure5I <- generate_plot(Figure5A.data, "Nucleotide metabolism")      #10
Figure5J <- generate_plot(Figure5A.data, "Signaling metabolism")       #11
Figure5K <- Figure5A.data %>%  
  filter(category == "Transport") %>% 
  ggplot(aes(y = dataset, x = description, fill = significance)) + 
  geom_tile() + 
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
        plot.margin = unit(c(0.05,0.05,0.05,0.05), "pt")) # order for plot margin: top, right, bottom, left

# Look to see what TIDEs are shared across all data sets
common.TIDEs <- Figure5A.data %>% 
  group_by(description) %>% 
  dplyr::summarize(count = sum(significance))


# Analyze TIDE data by reaction subsystem
SuppFigure2.data <- all.data %>%
  filter(grepl('S', ID)) 
SuppFigure2.data$dataset <- factor(SuppFigure2.data$dataset,
                                levels = c("GSE5406 Idiopathic","GSE5406 Ischemic",
                                           "GSE1869 Idiopathic","GSE1869 Ischemic",
                                           "GSE57345 Idiopathic","GSE57345 Ischemic"))
SuppFigure2.data$description <- as.factor(SuppFigure2.data$description)
SuppFigure2.data$description <- fct_rev(SuppFigure2.data$description)
levels(SuppFigure2.data$dataset) <- gsub(" ", "\n", levels(SuppFigure2.data$dataset))

SuppFigure2 <- SuppFigure2.data %>% 
  ggplot(aes(y = dataset, x = description, fill = significance)) + 
  geom_tile() + 
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

ggsave('figures/SuppFigure2.png', SuppFigure2, 
       dpi = 600, width = 4.8, height = 6)

### Figure 5B - GSEA data ####
# All GSEA data is saved in the GSEA data folder
GSE5406_Idiopathic <- read.delim(file = "data/GSE5406_Idiopathic_Healthy.txt") %>% 
  rbind(read.delim(file = "data/GSE5406_Idiopathic_Idiopathic.txt"))%>% 
  mutate(NOM.p.val = ifelse(NES < 0, -NOM.p.val, NOM.p.val)) %>% 
  select(NAME, NOM.p.val)
colnames(GSE5406_Idiopathic)[2] <- 'GSE5406 Idiopathic'

GSE5406_Ischemic <- read.delim(file = "data/GSE5406_Ischemic_Healthy.txt") %>% 
  rbind(read.delim(file = "data/GSE5406_Ischemic_Ischemic.txt"))%>% 
  mutate(NOM.p.val = ifelse(NES < 0, -NOM.p.val, NOM.p.val)) %>% 
  select(NAME, NOM.p.val)
colnames(GSE5406_Ischemic)[2] <- 'GSE5406 Ischemic'

GSE57345_Idiopathic <- read.delim(file = "data/GSE57345_Idiopathic_Healthy.txt") %>% 
  rbind(read.delim(file  = "data/GSE57345_Idiopathic_Idiopathic.txt"))%>% 
  mutate(NOM.p.val = ifelse(NES < 0, -NOM.p.val, NOM.p.val)) %>% 
  select(NAME, NOM.p.val)
colnames(GSE57345_Idiopathic)[2] <- 'GSE57345 Idiopathic'

GSE57345_Ischemic <- read.delim(file = "data/GSE57345_Ischemic_Healthy.txt") %>% 
  rbind(read.delim(file = "data/GSE57345_Ischemic_Ischemic.txt"))%>% 
  mutate(NOM.p.val = ifelse(NES < 0, -NOM.p.val, NOM.p.val)) %>% 
  select(NAME, NOM.p.val)
colnames(GSE57345_Ischemic)[2] <- 'GSE57345 Ischemic'

GSE1869_Idiopathic <- read.delim(file = "data/GSE1869_Idiopathic_Control.txt") %>% 
  rbind(read.delim(file  = "data/GSE1869_Idiopathic_Idiopathic.txt"))%>% 
  mutate(NOM.p.val = ifelse(NES < 0, -NOM.p.val, NOM.p.val)) %>% 
  select(NAME, NOM.p.val)
colnames(GSE1869_Idiopathic)[2] <- 'GSE1869 Idiopathic'

GSE1869_Ischemic <- read.delim(file = "data/GSE1869_Ischemic_Control.txt") %>% 
  rbind(read.delim(file = "data/GSE1869_Ischemic_Ischemic.txt"))%>% 
  mutate(NOM.p.val = ifelse(NES < 0, -NOM.p.val, NOM.p.val)) %>% 
  select(NAME, NOM.p.val)
colnames(GSE1869_Ischemic)[2] <- 'GSE1869 Ischemic'


SuppFigure1.data <- full_join(GSE5406_Idiopathic, GSE5406_Ischemic, by = c("NAME")) %>% 
  full_join(GSE57345_Idiopathic, by = c("NAME")) %>% 
  full_join(GSE57345_Ischemic, by = c("NAME")) %>% 
  full_join(GSE1869_Idiopathic, by = c("NAME")) %>% 
  full_join(GSE1869_Ischemic, by = c("NAME")) %>% 
  gather(key = dataset, value = NOM.p.val, 'GSE5406 Idiopathic':'GSE1869 Ischemic') %>% 
  mutate(NOM.p.val = ifelse(abs(NOM.p.val) >= 0.05, 0, ifelse(NOM.p.val < 0, -1, 1))) %>% 
  mutate(NOM.p.val = ifelse(is.na(NOM.p.val), 0, NOM.p.val))

common.sets <- SuppFigure1.data %>% group_by(NAME) %>% summarize(counts = sum(NOM.p.val))

SuppFigure1.data$dataset <- factor(SuppFigure1.data$dataset,
                                   levels = c("GSE5406 Idiopathic","GSE5406 Ischemic",
                                              "GSE1869 Idiopathic","GSE1869 Ischemic",
                                              "GSE57345 Idiopathic","GSE57345 Ischemic"))
levels(SuppFigure1.data$dataset) <- gsub(" ", "\n", levels(SuppFigure1.data$dataset))
SuppFigure1.data$NAME <- as.factor(SuppFigure1.data$NAME)
SuppFigure1.data$NAME <- fct_rev(SuppFigure1.data$NAME)

no_change <- SuppFigure1.data %>% 
  group_by(NAME) %>% 
  summarize(n = sum(abs(NOM.p.val))) %>% 
  filter(n != 0)
  

# Plot the results of all KEGG pathways
SuppFigure1 <- SuppFigure1.data %>% 
  inner_join(no_change) %>% 
  ggplot(aes(y = dataset, x = NAME, fill = NOM.p.val)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") +
  #scale_fill_gradientn(colours = c("white","blue","red","white")) +
  #scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  theme_bw() + 
  #coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6), 
        axis.text.y = element_text(hjust = 1, size = 8), 
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

# Save Supplemental Figure 1 - all KEGG pathways
ggsave('figures/SuppFigure1.png', SuppFigure1, 
       dpi = 600, width = 14, height = 6, units = 'in')


common.pathways <- SuppFigure1.data %>% 
  group_by(NAME) %>% 
  summarize(n = sum(abs(NOM.p.val))) %>% 
  filter(n >=4) %>% 
  inner_join(SuppFigure1.data)


# Filter to only include metabolic KEGG pathways
KEGG.metabolic.pathways <- read.table(file = "data/KEGG_metabolic_pathways.csv", 
                                     header = TRUE, sep = ",")

Figure5B.data <- inner_join(SuppFigure1.data, KEGG.metabolic.pathways, by = c("NAME"="KEGG.pathway")) %>% 
  select(-ID)
Figure5B.data$NAME <- factor(Figure5B.data$NAME, levels = KEGG.metabolic.pathways$KEGG.pathway)
Figure5B.data$NAME <- fct_rev(Figure5B.data$NAME)

# Plot the results of metabolic KEGG pathways
Figure5_GSEA <- Figure5B.data %>% 
  ggplot(aes(y = dataset, x = NAME, fill = NOM.p.val)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "skyblue2", mid = "white", high = "#FF9999") +
  #scale_fill_gradientn(colours = c("white","blue","red","white")) +
  #scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  #xlim(levels(Figure5B.data$NAME)) + 
  theme_bw() + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8), 
        axis.text.y = element_blank(), 
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        plot.margin = unit(c(0.05,0.05,0.05,0.05), "pt"))

gridlayout <- read.xlsx(file = "data/Figure5 gridlayout.xlsx", 
                        sheetIndex = 1,
                        header = FALSE) %>% 
  as.matrix(dimnames = NULL)

# Create base layout
Figure5 <- arrangeGrob(Figure5A, Figure5B, Figure5C,
                       Figure5D, Figure5E, Figure5F, 
                       Figure5G, Figure5H, Figure5I, 
                       Figure5J, Figure5K, 
                       Figure5_GSEA, 
                       ncol = 6, 
                       nrow = 345, 
                       widths = c(0.5/15, 3/15, 2.5/15, 3/15, 3/15, 3/15),
                       layout_matrix = cbind(gridlayout[,2],gridlayout,gridlayout[,2]))


adjust = 0.005
offset = 3.5/15
nrow = 345
category_location = c(21, 50, 79, 105, 143, 206.5, 238, 245, 262, 277, 294)
category_location = (nrow - category_location)/nrow
# Next, convert to a ggplot object and add labels for Figure 5A
Figure5 <- as_ggplot(Figure5) +
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
            y = category_location)

# Convert KEGG pathway names for plots
KEGG.pathways <- read.xlsx(file = "data/Convert KEGG pathway.xlsx", 
                           header = TRUE, sheetIndex = 1) %>% 
  right_join(KEGG.metabolic.pathways, by = c('NAME'='KEGG.pathway'))
#labels <- rev(as.character(KEGG.pathways$plot_name))
labels <- as.character(KEGG.pathways$plot_name)
Figure5 <- Figure5 + 
  draw_text(text = labels, 
            size = 8, 
            hjust = 0, 
            lineheight = 1, 
            x = rep(9/15+0.005, length(labels)), 
            y = (51:-1:1)/51*308/345 + 33/345)


# Add in labels for categories for KEGG pathways
labels <- levels(KEGG.pathways$TIDE.category)
nrow = 345
category_location = c(40, 101, 140, 172, 214, 246, 276, 302)
category_location = (nrow - category_location)/nrow
Figure5 <- Figure5 + 
  draw_text(text = labels, 
            size = 8, 
            hjust = 0, 
            lineheight = 1, 
            x = rep(12/15+0.05, length(labels)), 
            y = category_location)



# Add A, B, C labels to the plot
Figure5 <- Figure5 + 
  draw_plot_label(label = c('A','B'), 
                  size = 16, 
                  x = c(0.1,4.5)/12,
                  y = c(1, 1))

ggsave('figures/Figure5.png', Figure5, 
       dpi = 600, width = 12, height = 6, units = 'in')



# Clustering results from Figure5A.data
cluster.data <- Figure5A.data %>% 
  select(description, dataset, significance) %>% 
  spread(key = dataset, value = significance) %>% 
  select(-description)
cluster.data <- as.matrix(cluster.data)

# create heatmap and don't reorder columns
heatmap(cluster.data)

