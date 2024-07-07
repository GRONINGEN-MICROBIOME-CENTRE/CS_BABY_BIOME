############### CREATING FIGURE 1 B, C, D, E, F, G #######################
# Title: "CS_BABY_BIOME_FIGURE_1_COMBINE"
# Author: "Trishla Sinha"
# Date: "19/10/2023""
# Last update: "01/04/2024"

setwd("/Users/trishlasinha/Desktop/CS_Baby_Biome/submission/2024_submission/figures/figures_individual/figure_1/")
load("figure_1_D_E_F_G.RData")
load("figure_1_B_C.RData")

pdf('figure_1_B_C_D_E_F_G.pdf', width=10, height=12)

Figure_1 <- ggarrange(shannon_AB_CS_Baby_Biome, richness_AB_CS_Baby_Biome,overall, merged_plot, overall_AB_extra, tcam_extra_plot_AB,
                                     labels = c("B", "C", "D", "E", "F", "G"), 
                                     ncol = 2, nrow = 3)  
Figure_1


dev.off()
