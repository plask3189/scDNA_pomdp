# shortcuts: ctrl I is indent
library(scDNA)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(magick)
library(viridis)
library(readxl)
library(stringr)
library(HDF5Array)
source("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/scDNA copy/R/new_dev/visualize_clonal_evolution_kp.R")
source("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/GitHub/scDNA copy/R/clonograph.R")

sample_file<-"BRAF/A5330braf.dna+protein.h5"
#sample_file<- "BRAF/4629_braf.dna+protein.h5"
#sample_file<- "BRAF/M1912braf.dna+protein.h5"
#sample_file<- "BRAF/A0634braf.dna+protein.h5"
file_name <- basename(sample_file) # Extract the file name without the path
sample_name <- sub("\\..*", "", file_name) # Use sub() to remove "dna+protein.h5" ... everything after the first period (.)
data_kp = paste0("data_kp")

if (!dir.exists(data_kp)) { # Create the directory if it doesn't already exist
  dir.create(data_kp)
}
save_sample_data_path = paste0(data_kp,"/",sample_name)
if (!dir.exists(save_sample_data_path)) {# Create the directory if it doesn't already exist
  dir.create(save_sample_data_path)
}

variant_output<-variant_ID(file=sample_file,
                           panel="MSK_RL", # "UCSC" can be used for other panels
                           GT_cutoff=0,  # mimimum percent of cells where a successful genotyping call was made
                           VAF_cutoff=0) # mimimum variant allele frequency 
genes_of_interest <- c("IDH2","NRAS","NPM1","TET2","FLT3","IDH1")

variants_of_interest <- variant_output %>%
  dplyr::filter(Class != 'Intronic' & !is.na(Class)) %>%
  dplyr::filter(VAF > 0.01) %>%
  dplyr::filter(genotyping_rate > 85) %>%
  dplyr::filter(!is.na(CONSEQUENCE) & CONSEQUENCE != 'synonymous') %>%
  dplyr::filter(SYMBOL %in% genes_of_interest) %>%
  dplyr::filter(WT != 0) %>% # because if it is then that's fishy
  dplyr::arrange(desc(VAF)) %>%
  dplyr::slice(1:2)


sce <- tapestri_h5_to_sce(file = sample_file, variant_set = variants_of_interest)
sce@metadata[["sample_name"]]<-sample_name 
sce <- enumerate_clones(sce)
sce <- compute_clone_statistics(sce, skip_ploidy = FALSE)
clono<-clonograph(sce, complete_only = TRUE, num_bars_to_keep = 5, title=sample_name)
clono
sce <- trajectory_analysis(sce, use_ADO = FALSE)




final_vis<-visualize_tree(sce, variants_of_interest, remove_low_reward_edges = TRUE)
final_vis

clonograph_save_file_name <- file.path(save_sample_data_path, paste0(sample_name, ".png"))
network_name <- file.path(save_sample_data_path, paste0(sample_name, ".html"))
# Save the plot and network visualization
ggsave(clonograph_save_file_name)  # Save the plot
visSave(final_vis, file = network_name)  # Save the network visualization


