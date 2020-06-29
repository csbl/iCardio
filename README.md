#Identifying common shifts in metabolic functions from omics heart failure data using a new, cardiomyocyte-specific genome scale model of metabolism

### Abstract
The heart is a metabolic omnivore, known to consume many different carbon substrates in order to maintain function. In diseased states, the heart’s metabolism can shift between different carbon substrates; however, there is some disagreement in the field as to the metabolic shifts seen in end-stage heart failure and whether all heart failure converges to a common metabolic phenotype. Here, we present a new, validated cardiomyocyte-specific GEnome-scale metabolic Network REconstruction (GENRE), iCardio, and use the model to identify common shifts in metabolic functions across heart failure omics datasets. iCardio was built by integrating protein expression data from the Human Protein Atlas (HPA) with iHsa, a human GENRE of metabolism, followed by manual curation using pre-defined metabolic tasks. We demonstrate the utility of iCardio in interpreting heart failure gene expression data by identifying Tasks Inferred from Differential Expression (TIDEs) which represent metabolic functions associated with changes in gene expression for a diseased state. In contrast to outputs from other gene set enrichment approaches, TIDEs represent a reaction-centric approach that accounts for the complex relationship between gene expression and metabolic function that is captured through the gene-protein-reaction relationships in metabolic models. We can connect individual TIDEs back to changes in specific reactions and metabolic genes, thereby discovering potential metabolic contributors to, drug targets for, or biomarkers of heart failure. Here, we identify decreased NO and Neu5Ac synthesis as common metabolic markers of heart failure across datasets. Further, we highlight the differences in metabolic functions seen across studies, both in the TIDEs analysis and in GSEA analysis, further highlighting the complexity of heart failure, not only between patients but through the progression of the disease. 


### Overview

	project
	|- README          		# the top level description of content
  |
  |- code/              # code used for the presented analysis
  | |- MATLAB/          # code for reproducing building draft iCardio models, model curation, and TIDEs analysis
  | | |- functions/     # functions necessary for reproducing results
  | |- R/               # code for microarray DEG analysis
  |
  |- results/           # code to reproduce figures and supplementary tables
  | |- figures/
  | |- tables/
