# Analysis of Goyder - CLLMM revegetation project (Soil physchem data)
# RJH 2026

# Set working directory --------------------------------------------------------
setwd("/PATH/TO/DATA")

# libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# create output directory
outdir <- "/PATH/TO/DATA/R_output_soilpc"
# dir.create(file.path(outdir))
outdir <- file.path(outdir)

# Read data
all_soil_data <- read_excel("csbp-goyder-formatted.xlsx")
all_soil_metadata <- read_excel("Goyder-soil_analysis-metadata.xlsx")
all_soil_metadata$Site_ID <- paste0(all_soil_metadata$Site_ID, "P")

# Format correctly (convert to numeric variables where incorrectly characters)
# Also conver values like < 2 to 2 (limit of detection becomes the presumed value)
Texture_Categories <- data.frame(
  Texture <- c(1, 1.5, 2, 2.5, 3, 3.5),
  Texture_Description <- c("Sand", "Sand/Loam", "Loam", "Loam/Clay", "Clay", "Heavy Clay")
)
colnames(Texture_Categories) <- c("Texture", "Texture_Description")

all_soil_data <- all_soil_data %>%
  select( -c(`Lab Number\r\n`, `Date Received`, Depth, Colour, Gravel)) %>%
  rename(Site_ID = `Sample Name 2`) %>%
  mutate(across(
    where(is.character),
    ~ str_remove_all(.x, "< ")
  )) %>%
  mutate(across(-(1), as.numeric)) %>% 
  left_join(Texture_Categories , by = "Texture") %>% 
  select(-Texture) %>%
  merge(all_soil_metadata, by = "Site_ID", all = TRUE)

all_soil_data$Latitude <- ifelse(all_soil_data$Latitude > 0, all_soil_data$Latitude*-1, all_soil_data$Latitude)

# Comparisons ------------------------------------------------------------------
# View(all_soil_data)
ggplot(all_soil_data, aes(x = Longitude, y = Latitude, colour = Texture_Description))+
  geom_point(size=2.5)+
  theme_test()

colnames(all_soil_data)

soil_long <- all_soil_data %>%
  pivot_longer(
    cols = 2:20,
    names_to = "variable",
    values_to = "value"
  )

unique(soil_long$variable)
units_df <- tibble::tibble(
  variable = c(
    "Ammonium Nitrogen", "Nitrate Nitrogen",   "Phosphorus Colwell", "Potassium Colwell",  "Sulfur"       ,    
    "Organic Carbon"   , "Conductivity"    ,   "pH Level (CaCl2)"  , "pH Level (H2O)"   ,  "DTPA Copper"  ,    
    "DTPA Iron"        , "DTPA Manganese"  ,   "DTPA Zinc"         , "Exc. Aluminium"   ,  "Exc. Calcium" ,    
    "Exc. Magnesium"   , "Exc. Potassium"  ,   "Exc. Sodium"       , "Boron Hot CaCl2"  
  ),
  units = c(
    "mg kg-¹"          , "mg kg-¹"         ,   "mg kg-¹"           , "mg kg-¹"           , "mg kg-¹"      ,
    "%"                , "dS m-¹"          ,   "pH"                , "pH"                , "mg kg-¹"      ,
    "mg kg-¹"          , "mg kg-¹"         ,   "mg kg-¹"           , "meq 100g-¹"        , "meq 100g-¹"   ,
    "meq 100g-¹"       , "meq 100g-¹"      ,   "meq 100g-¹"        , "meq 100g-¹"
  )
)

# colours
soil_long$Treatment <- factor(soil_long$Treatment, levels = c("Degraded", "Revegetation", "Remnant"))

TreatCols <- c(
  "Degraded" = "#FFDB58", 
  "Revegetation" = "#008000", 
  "Remnant" = "#546B54"
)

soil_long <- soil_long %>%
  left_join(units_df, by = "variable")

## Treatment comparisons -------------------------------------------------------
ggplot(soil_long, aes(x = Treatment, y = value, fill=Treatment)) +
  geom_violin()+
  geom_boxplot(outlier.shape = NA, width = 0.1, fill = "white") +
  geom_jitter(width = 0.1, height = 0, alpha = 0.4, size = 1) +
  facet_wrap(.~ variable, scales = "free_y",
             labeller = labeller(
               variable = function(x) {
                 paste0(x, " (", soil_long$units[match(x, soil_long$variable)], ")")
                 }
               ))+
  scale_fill_manual(values = TreatCols)+
  labs(x = "Treatment", y = NULL) +  
  theme_test()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(filename = "Soil-boxplot-comparisons.pdf", path = outdir, width = 12, height = 10)

# log transform some variables
ggplot(soil_long, aes(x = Treatment, y = log(value), fill=Treatment)) +
  geom_violin()+
  geom_boxplot(outlier.shape = NA, width = 0.1, fill = "white") +
  geom_jitter(width = 0.1, height = 0, alpha = 0.4, size = 1) +
  facet_wrap(.~ variable, scales = "free_y",
             labeller = labeller(
               variable = function(x) {
                 paste0(x, " (", soil_long$units[match(x, soil_long$variable)], ")")
               }
             ))+
  scale_fill_manual(values = TreatCols)+
  labs(x = "Treatment", y = NULL) +
  # ggtitle("Log transformed variable")+
  theme_test()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(filename = "Soil-boxplot-comparisons-LOG.pdf", path = outdir, width = 12, height = 10)

## Management comparisons ------------------------------------------------------
ggplot(soil_long, aes(x = Latitude, y = value, fill=Treatment)) +
  geom_smooth(method = "lm")+
  geom_point() +
  facet_wrap(.~variable, scales = "free",
             labeller = labeller(
               variable = function(x) {
                 paste0(x, " (", soil_long$units[match(x, soil_long$variable)], ")")
               }
             ))+
  # scale_fill_manual(values = TreatCols)+
  theme_bw()
# ggsave(filename = "Soil-smooth-PC_Latitude-comparisons.pdf", path = outdir, width = 15, height = 15)

ggplot(soil_long, aes(x = Longitude, y = value, fill=Treatment)) +
  geom_smooth(method = "lm")+
  geom_point() +
  facet_wrap(.~variable, scales = "free",
             labeller = labeller(
               variable = function(x) {
                 paste0(x, " (", soil_long$units[match(x, soil_long$variable)], ")")
               }
             ))+
  # scale_fill_manual(values = TreatCols)+
  theme_bw()
# ggsave(filename = "Soil-smooth-PC_Longitude-comparisons.pdf", path = outdir, width = 15, height = 15)


# PCAs -------------------------------------------------------------------------
rownames(all_soil_data) <- all_soil_data$Site_ID
  
all_soil_data_pca <- all_soil_data[, sapply(all_soil_data, is.numeric)]

all_soil_data_pca <- all_soil_data_pca %>%
  select(-c(iEcosystem, iWptID, iTreatID, Latitude, Longitude, `pH Level (H2O)`))
metadata_pca <- all_soil_data[, sapply(all_soil_data, is.character)]


library(vegan)
all_soil_data_pca_std <- decostand(all_soil_data_pca, method = "standardize")
# View(all_soil_data_pca_std)

# Vegan
soil_pca <- rda(all_soil_data_pca_std)
summary(soil_pca)
plot(soil_pca)
biplot(soil_pca)
class(soil_pca)

# prep for ggplot2
# uscores is where ieach sample PC coordinates are
uscores <- data.frame(soil_pca$CA$u)
uscores$Site_ID <- rownames(uscores)

# sanity check
plot(uscores$PC1, uscores$PC2)
dim(all_soil_data_pca_std)
dim(uscores)

pca_data <- left_join(metadata_pca, uscores, by = "Site_ID")

# To get the variable coordinates across each principal component
vscores <- data.frame(soil_pca$CA$v)
vscores$variable <- rownames(vscores)

pca_data_Treatment_Centroids <- pca_data %>%
  select(Site_ID, Treatment, PC1, PC2) %>%
  group_by(Treatment) %>%
  reframe(
    Treatment_centroid_pc1 = mean(PC1),
    Treatment_centroid_pc2 = mean(PC2)
  )

pca_data_bigger <- left_join(
  pca_data, pca_data_Treatment_Centroids, by = "Treatment"
)

pca_data_bigger$Treatment <- factor(pca_data_bigger$Treatment, levels = c("Degraded", "Revegetation", "Remnant"))

## Plot PCA TREATMENT --------------------------------------------------------------------
# HULLs
data_for_hulls <- pca_data_bigger %>%
  select(Site_ID, PC1, PC2, Treatment, tEcoName, Texture_Description)

find_hull <- function(data_for_hulls) data_for_hulls[chull(data_for_hulls$PC1, data_for_hulls$PC2), ]
hulls <- plyr::ddply(data_for_hulls, "Treatment", find_hull)
hulls_tEcoName <- plyr::ddply(data_for_hulls, "tEcoName", find_hull)
hulls_Texture<- plyr::ddply(data_for_hulls, "Texture_Description", find_hull)

# Variance expained (Way longer than it should be!!!!!)
va_imp <- as.data.frame(summary(soil_pca)$cont$importance)
va_pc1 <- va_imp["Proportion Explained", "PC1"]
va_pc2 <- va_imp["Proportion Explained", "PC2"]

xlims <- range(range(pca_data_bigger$PC1), range(vscores$PC1))
ylims <- range(range(pca_data_bigger$PC2), range(vscores$PC2))
xrange <- diff(xlims)
yrange <- diff(ylims)
max_range <- max(xrange, yrange)

xmid <- mean(xlims)
ymid <- mean(ylims)
xlims_equal <- c(xmid - max_range/2, xmid + max_range/2)
ylims_equal <- c(ymid - max_range/2, ymid + max_range/2)

# ggplot
pca_treatment <- ggplot()+
  # 0 intercept lines
  geom_vline(xintercept = 0, linewidth = 0.1, colour = "black", linetype = "dotted")+
  geom_hline(yintercept = 0, linewidth = 0.1, colour = "black", linetype = "dotted")+
  # Hulls for treatments
  geom_polygon(data=hulls, aes(x = PC1, y = PC2,
    group = Treatment, fill = Treatment), alpha=0.2, colour = "black", linewidth=0.2)+
  # Fitted vectors (physchem variabeles)
  geom_segment(data = vscores, colour="blue", linewidth = 0.2, alpha = 0.5,
               aes(x = 0, y = 0, xend = PC1, yend = PC2))+
  geom_text(data = vscores, colour="blue", hjust = 1, 
            aes(x = PC1, y = PC2, label = variable))+
  # sample points (physchem variabeles). fitted twice to get nice black border
  geom_point(data = pca_data_bigger, size = 2.5, colour = "black",
             aes(x = PC1, y = PC2, shape = Treatment))+
  geom_point(data = pca_data_bigger, size = 2, 
             aes(x = PC1, y = PC2, colour = Treatment, shape = Treatment))+
  # # Centroids for Treatments
  # geom_point(data = pca_data_bigger, size = 5, colour = "black",
  #            aes(x = Treatment_centroid_pc1, y = Treatment_centroid_pc2,
  #                shape = Treatment))+
  # geom_point(data = pca_data_bigger, size = 4,
  #            aes(x = Treatment_centroid_pc1, y = Treatment_centroid_pc2,
  #                shape = Treatment, colour = Treatment))+
  # geom_segment(data = pca_data_bigger, alpha = 0.5,
  #              aes(x = Treatment_centroid_pc1,
  #                  y = Treatment_centroid_pc2,
  #                  xend = PC1, yend = PC2,
  #                  colour = Treatment))+
  # Fine tuning aesthetics
  scale_colour_manual(values = TreatCols)+
  scale_fill_manual(values = TreatCols)+
  coord_equal(ratio = 1 , xlim=xlims_equal, ylim = ylims_equal)+
  labs(x = paste0("PC1 (", round(va_pc1, 3)*100, "%)"), 
       y = paste0("PC2 (", round(va_pc2, 3)*100, "%)"))+
  theme_test()
pca_treatment
# ggsave(filename = "Soil-PCA-treatment.pdf", path = outdir, width = 8, height = 6)

# Soil physchem treatment stats
# all_soil_data_pca_std <- decostand(all_soil_data_pca, method = "standardize")

PhysChem_dist_mat <- dist(all_soil_data_pca_std)  
rownames(pca_data) <- pca_data$Site_ID

adonis2(PhysChem_dist_mat ~ Treatment, data = pca_data)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# adonis2(formula = PhysChem_dist_mat ~ Treatment, data = pca_data)
#           Df SumOfSqs      R2      F Pr(>F)  
# Treatment  2     91.8 0.09445 2.7117  0.014 *
# Residual  52    880.2 0.90555                

B_disp_model <- betadisper(PhysChem_dist_mat, pca_data$Treatment)
permutest(B_disp_model, permutations = 999)

# tEcoName
unique(pca_data_bigger$tEcoName)

EcoColours <- c(
  "NA"            = "darkgrey",
  "E.fasciculosa" = "#1B9E77",
  "C.gracilis"    = "#D95F02",
  "E.diversifolia"= "#7570B3",
  "A.verticillata"= "#E7298A",
  "E.porosa"      = "#66A61E",
  "E.leucoxylon"  = "#E6AB02",
  "E.odorata"     = "#A6761D",
  "E.incrassata"  = "#666666"
  )

pca_eco <- ggplot()+
  # 0 intercept lines
  geom_vline(xintercept = 0, linewidth = 0.1, colour = "black", linetype = "dotted")+
  geom_hline(yintercept = 0, linewidth = 0.1, colour = "black", linetype = "dotted")+
  # Hulls for treatments
  geom_polygon(data=hulls_tEcoName, aes(x = PC1, y = PC2,
                               group = tEcoName, fill = tEcoName), alpha=0.2, colour = "black", linewidth=0.2)+
  # Fitted vectors (physchem variabeles)
  geom_segment(data = vscores, colour="blue", linewidth = 0.2, alpha = 0.5,
               aes(x = 0, y = 0, xend = PC1, yend = PC2))+
  geom_text(data = vscores, colour="blue", hjust = 1, 
            aes(x = PC1, y = PC2, label = variable))+
  # sample points (physchem variabeles). fitted twice to get nice black border
  geom_point(data = pca_data_bigger, size = 2.5, colour = "black",
             aes(x = PC1, y = PC2))+
  geom_point(data = pca_data_bigger, size = 2, 
             aes(x = PC1, y = PC2, colour = tEcoName))+
  # Fine tuning aesthetics
  scale_colour_manual(values = EcoColours, name= "Ecosystem")+
  scale_fill_manual(values = EcoColours, name= "Ecosystem")+
  coord_equal(ratio = 1 , xlim=xlims_equal, ylim = ylims_equal)+
  labs(x = paste0("PC1 (", round(va_pc1, 3)*100, "%)"), 
       y = paste0("PC2 (", round(va_pc2, 3)*100, "%)"))+
  theme_test()
pca_eco
# ggsave(filename = "Soil-PCA-ecosystem.pdf", path = outdir, width = 8, height = 6)

TexColours <- c(
  "Loam"      = "#4DA3FF", 
  "Loam/Clay" = "#FF6F61", 
  "Sand"      = "#C77CFF",   
  "Sand/Loam" = "#7FDBB6"
  )

pca_texture <- ggplot()+
  # 0 intercept lines
  geom_vline(xintercept = 0, linewidth = 0.1, colour = "black", linetype = "dotted")+
  geom_hline(yintercept = 0, linewidth = 0.1, colour = "black", linetype = "dotted")+
  # Hulls for treatments
  geom_polygon(data=hulls_Texture, aes(x = PC1, y = PC2,
                                        group = Texture_Description, fill = Texture_Description), alpha=0.2, colour = "black", linewidth=0.2)+
  # Fitted vectors (physchem variabeles)
  geom_segment(data = vscores, colour="blue", linewidth = 0.2, alpha = 0.5,
               aes(x = 0, y = 0, xend = PC1, yend = PC2))+
  geom_text(data = vscores, colour="blue", hjust = 1, 
            aes(x = PC1, y = PC2, label = variable))+
  # sample points (physchem variabeles). fitted twice to get nice black border
  geom_point(data = pca_data_bigger, size = 2.5, colour = "black",
             aes(x = PC1, y = PC2))+
  geom_point(data = pca_data_bigger, size = 2, 
             aes(x = PC1, y = PC2, colour = Texture_Description))+
  # Fine tuning aesthetics
  scale_colour_manual(values = TexColours, name= "Texture")+
  scale_fill_manual(values = TexColours, name= "Texture")+
  coord_equal(ratio = 1 , xlim=xlims_equal, ylim = ylims_equal)+
  labs(x = paste0("PC1 (", round(va_pc1, 3)*100, "%)"), 
       y = paste0("PC2 (", round(va_pc2, 3)*100, "%)"))+
  theme_test()
pca_texture
# ggsave(filename = "Soil-PCA-texture.pdf", path = outdir, width = 8, height = 6)

library(ggpubr)
ggarrange(pca_treatment, pca_eco, pca_texture, align = "hv",
          labels = c("a", "b", "c"), ncol = 1)
# ggsave(filename = "Soil-PCA-combined.pdf", path = outdir, width = 6, height = 14)

# Correllation matrix and plots ----------------------------------------------

library(corrplot)
cor_pc <-cor(all_soil_data_pca)
cor_pc_res <- cor.mtest(all_soil_data_pca, conf.level = .95)
cor_pc_res <- ifelse()
corplot_obj <- corrplot(cor_pc, type="upper", order="hclust", "ellipse",
         p.mat = cor_pc_res$p)




