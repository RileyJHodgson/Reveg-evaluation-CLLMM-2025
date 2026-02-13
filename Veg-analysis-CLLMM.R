# Analysis of Goyder - CLLMM revegetation project (veg data)
# RJH 2026

# Set working directory -------------------
setwd("/PATH/TO/DATA/data_2025")

# libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)

# create output directory
outdir <- "/PATH/TO/DATA/data_2025/R_output_veg_updated"
# dir.create(file.path(outdir))
outdir <- file.path(outdir)

# Read data
# metadata_s <- read_excel("Goyder-veg_analysis-metadata.xlsx", sheet = "Site level")
metadata_v <- read_excel("Goyder-veg_analysis-metadata.xlsx", sheet = "Visit level")
all_veg_data <- read_excel("Joint_Vegetation_data_210126.xlsx", sheet= "Final_Data")
colnames(metadata_v)

dim(metadata_v)
dim(all_veg_data)

unique(all_veg_data$tEcoName...20)
unique(all_veg_data$tEcoName)

# rename iLump, more accurately
all_veg_data$iLumpEco

ilumpDF <- data.frame(
  "iLumpEco" = c(1, 4, 5, 6, 10),
  "tEcoName" = c("E. fasciculosa", "E. diversifolia", "A. verticillata", "Mixed mallee woodland", "C. gracilis"),
  "tEcosystem_Desc" = c("Pink Gum (Eucalyptus fasciculosa) Low Open Grassy Woodland of the Mount Lofty Ranges", "Eucalyptus diversifolia Mallee Communities of the South East", "Sheoak low woodland with shrubby understorey", "Mallee Grassy Woodland", "Non eucalypt grassy woodland")
)
all_veg_data <- left_join(all_veg_data, ilumpDF, by = "iLumpEco")
all_veg_data %>%
  select(tEcoName, tEcoName...20, iLumpEco) %>%
  distinct()

# Hindmarsh Island Landcare Group # not sure if 619 is included (Mundoo Island)
HILG_sites <- c(586, 669, 668, 580, 653, 654, 580)

all_veg_data <- all_veg_data %>%
  mutate(HILG = if_else(WptID %in% HILG_sites, "HILG", NA))

# add addditional metadata on sampling sites
all_veg_data <- left_join(all_veg_data, metadata_v, by = c("WptID"="iWptID"))
all_veg_data <- all_veg_data %>% select(-ogr_Survey)
# write.csv(all_veg_data, file ="Joint_Vegetation_data_210126_COMBINED.csv")

all_veg_data$Survey_Year <- as.factor(ifelse(all_veg_data$ProjectID == 3, "2015", "2025"))

# Adjust variable for easier downstream analysis
all_veg_data$Site_ID <- paste0("G", all_veg_data$WptID,"_", all_veg_data$Survey_Year)

# Extract cover codes for relative cover proportion (proxy for abundance)
cover_codes <- read_excel("Joint_Vegetation_data_210126.xlsx", sheet= "Cover")
cover_codes$midpoints_perc <- c(NA, 0.5, 1, 2.5, 15, 37.5, 62.5, 87.5)

midpoints_perc_vec <- c("0"=0, "N"=0.5, "T"=1, "1"=2.5, "2"=15, "3"=37.5, "4"=62.5, "5"=87.5)

# Add new relative cover values to maind dataset
all_veg_data <- all_veg_data %>% mutate(cover_midpoint = midpoints_perc_vec[COVCODE])

# get average species abundances per visit

all_veg_data_v <- all_veg_data
all_veg_data_v$cover_midpoint <- unname(all_veg_data_v$cover_midpoint)

all_veg_data_v$VisitID <- paste0("G", all_veg_data_v$WptID,"_", all_veg_data_v$VisitID, "_y", all_veg_data_v$Survey_Year)
all_veg_data_v <- as.data.frame(all_veg_data_v)

all_veg_data_v_com <- all_veg_data_v %>%
  distinct() %>% #make sure we don't have duplicates
  select(VisitID, `SCIENTIFIC NAME`, cover_midpoint) %>%
  group_by(VisitID, `SCIENTIFIC NAME`) %>%
  summarise(total_cover = sum(cover_midpoint, na.rm = TRUE)) %>%
  as.data.frame()
nrow(all_veg_data_v)
nrow(all_veg_data_v_com)

veg_species_wide <- all_veg_data_v_com %>%
  pivot_wider(
    id_cols = VisitID,             
    names_from = `SCIENTIFIC NAME`,
    values_from = total_cover,
    values_fill = 0) %>% 
  as.data.frame()

rownames(veg_species_wide) <- veg_species_wide$VisitID
veg_species_wide <- veg_species_wide %>% select(-VisitID)

# Add management data
management_data <- read_excel("Veg_Survey_data_2025_v3.xlsx", sheet = "ManagementIssues")
colnames(management_data)

management_data <- management_data %>%
  mutate(grazing_level = case_when(
    str_detect(`Domestic grazing`, regex("Low", ignore_case = TRUE))    ~ "Low",
    str_detect(`Domestic grazing`, regex("Moderate", ignore_case = TRUE)) ~ "Moderate",
    str_detect(`Domestic grazing`, regex("High", ignore_case = TRUE))   ~ "High",
    TRUE                                                       ~ NA_character_
  ))
# View(management_data)

grazing_info_only <- management_data %>%
  select(WptID, grazing_level, `Domestic grazing`)

grazing_info_only_tab <- grazing_info_only
grazing_info_only_tab$grazing_level <- ifelse(is.na(grazing_info_only_tab$grazing_level) == TRUE, "No grazing observed", grazing_info_only_tab$grazing_level)
table(grazing_info_only_tab$grazing_level)
table(grazing_info_only$grazing_level)


# Alpha Diversity --------------------------------------------------------------
library(vegan)
veg_species_wide_R <- specnumber(veg_species_wide)
veg_species_wide_Sh <- diversity(veg_species_wide, index = "shannon")
veg_species_wide_Si <- diversity(veg_species_wide, index = "simpson")
veg_species_wide_ISi <- diversity(veg_species_wide, index = "invsimpson")
diversity_data <- data.frame(VisitID = names(veg_species_wide_R),
                             Richness = veg_species_wide_R,
                             Shannons = veg_species_wide_Sh,
                             Simpsons = veg_species_wide_Si,
                             InvSimpsons = veg_species_wide_ISi)

diversity_data$EffSpp <- exp(diversity_data$Shannons)

colnames(all_veg_data_v)
all_veg_site <- all_veg_data_v %>% 
  select(Survey_Year, WptID, iPlantYear, Treatment, VisitID, iLumpEco, tEcoName, Site_ID, HILG) %>%
  distinct()
head(all_veg_site)

colnames(all_veg_site)
colnames(all_veg_site)

diversity_data2 <- left_join(diversity_data, all_veg_site, by = "VisitID")
colnames(diversity_data2)
diversity_data2$RemRev <- ifelse(diversity_data2$iPlantYear == "Remnant", "Remnant", "Revegetation")
diversity_data2$RemRev <- factor(diversity_data2$RemRev, levels = c("Revegetation", "Remnant"))

unique(diversity_data2$iLumpEco)
unique(diversity_data2$tEcoName)

ggplot(diversity_data2, aes(x= iPlantYear, y= EffSpp))+
  geom_violin(fill = "brown")+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(Survey_Year ~ tEcoName, scales = "free_x", space = "free_x")+
  theme_test()
# ggsave(filename = "Veg-EffectiveSpecies-PlantYear.pdf", path = outdir, width = 14, height = 7)

remrev_cols <- c("Remnant" = "brown", "Revegetation" = "orange")
ggplot(diversity_data2, aes(x= RemRev, y= EffSpp))+
  geom_violin(aes(fill = RemRev))+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(Survey_Year ~ tEcoName, scales = "free_x", space = "free_x")+
  scale_fill_manual(values = remrev_cols)+
  labs(x = "Treatment", y = "Effective number of species")+
  theme_test()
# ggsave(filename = "Veg-EffectiveSpecies-RemRev.pdf", path = outdir, width = 14, height = 7)

diversity_data3 <- diversity_data2 %>%
  group_by(Survey_Year, WptID, iPlantYear, Treatment, iLumpEco, tEcoName, RemRev, HILG) %>%
  reframe(mean_EffSpp =  mean(EffSpp),
          mean_Richness = mean(Richness),
          mean_Shannons = mean(Shannons),
          mean_InvSimps = mean(InvSimpsons))

ggplot(diversity_data3, aes(x= RemRev, y= mean_EffSpp))+
  geom_violin(aes(fill = RemRev))+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(.~Survey_Year, scales = "free_x", space = "free_x")+
  scale_fill_manual(values = remrev_cols)+
  theme_test()
# ggsave(filename = "Veg-EffectiveSpecies-REMREV-fac_survey_year.pdf", path = outdir, width = 7, height = 4)

ggplot(diversity_data3, aes(x= tEcoName, y= mean_EffSpp))+
  geom_violin(aes(fill = tEcoName))+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(RemRev~Survey_Year, scales = "free_x", space = "free_x")+
  scale_fill_brewer(palette = "Set1")+
  theme_test()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# ggsave(filename = "Veg-EffectiveSpecies-Ecosystem-bySurveyRemRev.pdf", path = outdir, width = 10, height = 7)

# QUICK HILG check
ggplot(diversity_data3, aes(x= tEcoName, y= mean_EffSpp))+
  geom_violin(aes(fill = tEcoName))+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 1, aes(colour = HILG)) +
  facet_grid(RemRev~Survey_Year, scales = "free_x", space = "free_x")+
  scale_fill_brewer(palette = "Set1")+
  scale_colour_manual(values = c("HILG"= "red"), na.value = "black")+
  theme_test()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

all_alpha_plot_Rich <- ggplot(diversity_data3, aes(x= RemRev, y= mean_Richness))+
  geom_violin(aes(fill = RemRev))+ geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) + 
  facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
  scale_fill_brewer(palette = "Set2")+ theme_test()+ theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(y= "Richness", x = "") 
all_alpha_plot_Inv <- ggplot(diversity_data3, aes(x= RemRev, y= mean_InvSimps))+
  geom_violin(aes(fill = RemRev))+ geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
  scale_fill_brewer(palette = "Set2")+ theme_test()+ theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(y= "Inverse Simpson's diversity", x = "")
all_alpha_plot_sha <- ggplot(diversity_data3, aes(x= RemRev, y= mean_Shannons))+
  geom_violin(aes(fill = RemRev))+ geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
  scale_fill_brewer(palette = "Set2")+ theme_test()+ theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(y= "Shannon's diversity index", x = "")
all_alpha_plot_EFS <- ggplot(diversity_data3, aes(x= RemRev, y= mean_EffSpp))+
  geom_violin(aes(fill = RemRev))+ geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
  scale_fill_brewer(palette = "Set2")+ theme_test()+ theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(y= "Effective number of species", x = "")


all_alpha_plots <- ggpubr::ggarrange(all_alpha_plot_Rich, all_alpha_plot_Inv, all_alpha_plot_sha, all_alpha_plot_EFS, common.legend = TRUE, align = "hv")

all_alpha_plots
# ggsave(plot = all_alpha_plots, filename = "Veg-all_alpha_diversity.pdf", path = outdir, width = 10, height = 8)

all_alpha_plot_Rich <- ggplot(diversity_data3, aes(x= RemRev, y= mean_Richness))+
  geom_violin(aes(fill = RemRev))+ geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 1, aes(colour = HILG)) + 
  facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
  scale_fill_brewer(palette = "Set2")+ theme_test()+ theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(y= "Richness", x = "") +
  scale_colour_manual(values = c("HILG"= "red"), na.value = "black")
all_alpha_plot_Inv <- ggplot(diversity_data3, aes(x= RemRev, y= mean_InvSimps))+
  geom_violin(aes(fill = RemRev))+ geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 1, aes(colour = HILG)) + 
  facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
  scale_fill_brewer(palette = "Set2")+ theme_test()+ theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(y= "Inverse Simpson's diversity", x = "")+
  scale_colour_manual(values = c("HILG"= "red"), na.value = "black")
all_alpha_plot_sha <- ggplot(diversity_data3, aes(x= RemRev, y= mean_Shannons))+
  geom_violin(aes(fill = RemRev))+ geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 1, aes(colour = HILG)) + 
  facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
  scale_fill_brewer(palette = "Set2")+ theme_test()+ theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(y= "Shannon's diversity index", x = "")+
  scale_colour_manual(values = c("HILG"= "red"), na.value = "black")
all_alpha_plot_EFS <- ggplot(diversity_data3, aes(x= RemRev, y= mean_EffSpp))+
  geom_violin(aes(fill = RemRev))+ geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 1, aes(colour = HILG)) + 
  facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
  scale_fill_brewer(palette = "Set2")+ theme_test()+ theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(y= "Effective number of species", x = "")+
  scale_colour_manual(values = c("HILG"= "red"), na.value = "black")


all_alpha_plots_HILG <- ggpubr::ggarrange(all_alpha_plot_Rich, all_alpha_plot_Inv, all_alpha_plot_sha, all_alpha_plot_EFS, common.legend = TRUE, align = "hv")

## Alpha diversity stats ------------------------------------------------------
library(lme4)
getwd() # "PATH/TO/DATA/data_2025"
source("Permute-LMEM-Toolkit.R")

# LMEM effective no species by treatment and sample year
LMEM_eff_RRSY_int <- LMEM_permute_anova( data = diversity_data2,
                                         response = "EffSpp", fixed_effects = c("RemRev", "Survey_Year", "RemRev:Survey_Year"),
                                         random_effects = c("WptID", "tEcoName"),
                                         nreps = 999)

LMEM_eff_RRSY_int$formula
# EffSpp ~ RemRev + Survey_Year + RemRev:Survey_Year + (1 | WptID) + (1 | tEcoName)

LMEM_eff_RRSY_int$Anova
# Fixed effects:
#                      Chisq Df Pr(>Chisq)    
# RemRev             116.953  1  < 2.2e-16 ***
# Survey_Year         14.711  1  0.0001253 ***
# RemRev:Survey_Year  22.725  1  1.869e-06 ***

LMEM_eff_RRSY_int$permutation_p_values
# RemRev        Survey_Year RemRev:Survey_Year 
#      0                  0                  0 

# LMEM effective no species.
LMEM_eff_EcSY_int <- LMEM_permute_anova(data = diversity_data2,
                                        response = "EffSpp", fixed_effects = c("tEcoName", "Survey_Year", "tEcoName:Survey_Year"),
                                        random_effects = c("WptID", "RemRev"),
                                        nreps = 999)
LMEM_eff_EcSY_int$formula
# EffSpp ~ tEcoName + Survey_Year + tEcoName:Survey_Year + (1 | WptID) + (1 | RemRev)

LMEM_eff_EcSY_int$Anova
#                        Chisq Df Pr(>Chisq)    
# tEcoName              3.0518  4  0.5492022    
# Survey_Year          14.6354  1  0.0001304 ***
# tEcoName:Survey_Year 17.1937  4  0.0017724 ** 

LMEM_eff_EcSY_int$permutation_p_values
#    tEcoName          Survey_Year tEcoName:Survey_Year
# 0.529529530          0.000000000          0.001001001

### Native plant diversity only ------------------------------------------------
Native_veg_data <- subset(all_veg_data_v, is.na(INTRODUCED) == TRUE)

Native_veg_com <- Native_veg_data %>%
  distinct() %>% #make sure we don't have duplicates
  select(VisitID, `SCIENTIFIC NAME`, cover_midpoint) %>%
  group_by(VisitID, `SCIENTIFIC NAME`) %>%
  summarise(total_cover = sum(cover_midpoint, na.rm = TRUE)) %>%
  as.data.frame()
nrow(Native_veg_data) # 14635
nrow(Native_veg_com)  # 14599

Native_veg_wide <- Native_veg_com %>%
  pivot_wider(
    id_cols = VisitID,             
    names_from = `SCIENTIFIC NAME`,
    values_from = total_cover,
    values_fill = 0) %>% 
  as.data.frame()

rownames(Native_veg_wide) <- Native_veg_wide$VisitID
Native_veg_wide <- Native_veg_wide %>% select(-VisitID)

# Get diversity data
Native_veg_wide_R <- specnumber(Native_veg_wide)
Native_veg_wide_Sh <- diversity(Native_veg_wide, index = "shannon")
Native_veg_wide_Si <- diversity(Native_veg_wide, index = "simpson")
Native_veg_wide_ISi <- diversity(Native_veg_wide, index = "invsimpson")

Native_diversity_data <- data.frame(VisitID = names(Native_veg_wide_R),
                                    Richness = Native_veg_wide_R,
                                    Shannons = Native_veg_wide_Sh,
                                    Simpsons = Native_veg_wide_Si,
                                    InvSimpsons = Native_veg_wide_ISi)
Native_diversity_data$EffSpp <- exp(Native_diversity_data$Shannons)

all_veg_site_native <- all_veg_data_v %>% 
  select(Survey_Year, WptID, iPlantYear, Treatment, VisitID, iLumpEco, tEcoName, Site_ID, VisitID, HILG) %>%
  distinct()
head(all_veg_site_native)

colnames(Native_diversity_data)
colnames(all_veg_site_native)

Native_diversity_data2 <- left_join(Native_diversity_data, all_veg_site_native, by = "VisitID")
colnames(Native_diversity_data2)
Native_diversity_data2$RemRev <- ifelse(Native_diversity_data2$iPlantYear == "Remnant", "Remnant", "Revegetation")
Native_diversity_data2$RemRev <- factor(Native_diversity_data2$RemRev, levels = c("Revegetation", "Remnant"))

unique(Native_diversity_data2$tEcoName)
# Native_diversity_data2$tEcoName <- factor(Native_diversity_data2$tEcoName, levels = c("A.verticillata", "C.gracilis", "E.diversifolia", "E.fasciculosa", "E.porosa", "E.incrassata", "E.leucoxylon", "E.odorata"))

ggplot(Native_diversity_data2, aes(x= iPlantYear, y= EffSpp))+
  geom_violin(fill = "darkgreen")+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(Survey_Year ~ tEcoName, scales = "free_x", space = "free_x")+
  theme_test()
# ggsave(filename = "Veg-EffectiveSpecies-PlantYear-NATIVEONLY.pdf", path = outdir, width = 14, height = 7)

remrev_cols <- c("Remnant" = "darkgreen", "Revegetation" = "#9DC183")
ggplot(Native_diversity_data2, aes(x= RemRev, y= EffSpp))+
  geom_violin(aes(fill = RemRev))+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(Survey_Year ~ tEcoName, scales = "free_x", space = "free_x")+
  scale_fill_manual(values = remrev_cols)+
  theme_test()
# ggsave(filename = "Veg-EffectiveSpecies-RemRev-NATIVEONLY.pdf", path = outdir, width = 14, height = 7)

Native_diversity_data3 <- Native_diversity_data2 %>%
  group_by(Survey_Year, WptID, iPlantYear, Treatment, iLumpEco, tEcoName, RemRev, HILG) %>%
  reframe(mean_EffSpp =  mean(EffSpp))

ggplot(Native_diversity_data3, aes(x= RemRev, y= mean_EffSpp))+
  geom_violin(aes(fill = RemRev))+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 1, aes(colour = HILG)) +
  facet_grid(.~Survey_Year, scales = "free_x", space = "free_x")+
  scale_fill_manual(values = remrev_cols)+
  theme_test()
# ggsave(filename = "Veg-EffectiveSpecies-PlantYear-NATIVEONLY-survey_year.pdf", path = outdir, width = 7, height = 4)

ggplot(Native_diversity_data3, aes(x= tEcoName, y= mean_EffSpp))+
  geom_violin(aes(fill = tEcoName))+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(RemRev~Survey_Year, scales = "free_x", space = "free_x")+
  scale_fill_brewer(palette = "Set2")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  theme_test()
# ggsave(filename = "Veg-EffectiveSpecies-Ecosystem-NATIVEONLY-bySurveyRemRev.pdf", path = outdir, width = 7, height = 4)

# LUMPED

ggplot(Native_diversity_data2, aes(x= iPlantYear, y= EffSpp))+
  geom_violin(fill = "darkgreen")+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(Survey_Year ~ iLumpEco, scales = "free_x", space = "free_x")+
  theme_test()
# ggsave(filename = "Veg-EffectiveSpecies-PlantYear-NATIVEONLY-LumpedEco.pdf", path = outdir, width = 14, height = 7)

remrev_cols <- c("Remnant" = "darkgreen", "Revegetation" = "#9DC183")
ggplot(Native_diversity_data2, aes(x= RemRev, y= EffSpp))+
  geom_violin(aes(fill = RemRev))+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(Survey_Year ~ iLumpEco, scales = "free_x", space = "free_x")+
  scale_fill_manual(values = remrev_cols)+
  theme_test()
# ggsave(filename = "Veg-EffectiveSpecies-RemRev-NATIVEONLY-LumpedEco.pdf", path = outdir, width = 14, height = 7)

Native_diversity_data3 <- Native_diversity_data2 %>%
  group_by(Survey_Year, WptID, iPlantYear, Treatment, iLumpEco, tEcoName, RemRev) %>%
  reframe(mean_EffSpp =  mean(EffSpp))

ggplot(Native_diversity_data3, aes(x= RemRev, y= mean_EffSpp))+
  geom_violin(aes(fill = RemRev))+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(. ~ Survey_Year, scales = "free_x", space = "free_x")+
  scale_fill_manual(values = remrev_cols)+
  theme_test()
# ggsave(filename = "Veg-EffectiveSpecies-REMREV-NATIVEONLY-fac_survey_year.pdf", path = outdir, width = 7, height = 4)

ggplot(Native_diversity_data3, aes(x= tEcoName, y= mean_EffSpp))+
  geom_violin(aes(fill = tEcoName))+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(RemRev ~ Survey_Year, scales = "free_x", space = "free_x")+
  scale_fill_brewer(palette = "Set1")+
  theme_test()
# ggsave(filename = "Veg-EffectiveSpecies-Ecosystem-NATIVEONLY-fac_survey_year.pdf", path = outdir, width = 7, height = 4)

# stats
model1_native <- LMEM_permute_anova(data = Native_diversity_data2,
                                    response = "EffSpp", fixed_effects = c("RemRev", "Survey_Year", "RemRev:Survey_Year"),
                                    random_effects = c("WptID", "tEcoName"),
                                    nreps = 999)

model1_native$formula
# EffSpp ~ RemRev + Survey_Year + RemRev:Survey_Year + (1 | WptID) + (1 | tEcoName)

model1_native$Anova
#                     Chisq Df Pr(>Chisq)    
# RemRev             89.716  1  < 2.2e-16 ***
# Survey_Year        57.041  1  4.268e-14 ***
# RemRev:Survey_Year 71.960  1  < 2.2e-16 ***

model2_native <- LMEM_permute_anova(data = Native_diversity_data2,
                                    response = "EffSpp", fixed_effects = c("tEcoName", "Survey_Year", "tEcoName:Survey_Year"),
                                    random_effects = c("WptID", "RemRev"),
                                    nreps = 999)

model2_native$formula
# EffSpp ~ tEcoName + Survey_Year + tEcoName:Survey_Year + (1 | WptID) + (1 | RemRev)
model2_native$Anova
#                        Chisq Df Pr(>Chisq)    
# tEcoName              3.9797  4     0.4088    
# Survey_Year          53.9399  1  2.067e-13 ***
# tEcoName:Survey_Year  4.8644  4     0.3015    

### Invasive plant diversity only ----------------------------------------------
Invasive_veg_data <- subset(all_veg_data_v, INTRODUCED == "*")

Invasive_veg_com <- Invasive_veg_data %>%
  distinct() %>% #make sure we don't have duplicates
  select(VisitID, `SCIENTIFIC NAME`, cover_midpoint) %>%
  group_by(VisitID, `SCIENTIFIC NAME`) %>%
  summarise(total_cover = sum(cover_midpoint, na.rm = TRUE)) %>%
  as.data.frame()
nrow(Invasive_veg_data) # 99108
nrow(Invasive_veg_com)  # 11010

Invasive_veg_wide <- Invasive_veg_com %>%
  pivot_wider(
    id_cols = VisitID,             
    names_from = `SCIENTIFIC NAME`,
    values_from = total_cover,
    values_fill = 0) %>% 
  as.data.frame()

rownames(Invasive_veg_wide) <- Invasive_veg_wide$VisitID
Invasive_veg_wide <- Invasive_veg_wide %>% select(-VisitID)

# Get diversity data
Invasive_veg_wide_R <- specnumber(Invasive_veg_wide)
Invasive_veg_wide_Sh <- diversity(Invasive_veg_wide, index = "shannon")
Invasive_veg_wide_Si <- diversity(Invasive_veg_wide, index = "simpson")
Invasive_veg_wide_ISi <- diversity(Invasive_veg_wide, index = "invsimpson")

Invasive_diversity_data <- data.frame(VisitID = names(Invasive_veg_wide_R),
                                      Richness = Invasive_veg_wide_R,
                                      Shannons = Invasive_veg_wide_Sh,
                                      Simpsons = Invasive_veg_wide_Si,
                                      InvSimpsons = Invasive_veg_wide_ISi)
Invasive_diversity_data$EffSpp <- exp(Invasive_diversity_data$Shannons)

# Recombine with metadata
all_veg_site_Invasive <- all_veg_data_v %>% 
  select(Survey_Year, WptID, iPlantYear, Treatment, VisitID, iLumpEco, tEcoName, Site_ID) %>%
  distinct()
head(all_veg_site_Invasive)

colnames(Invasive_diversity_data)
colnames(all_veg_site_Invasive)

Invasive_diversity_data2 <- left_join(Invasive_diversity_data, all_veg_site_Invasive, by = "VisitID")
colnames(Invasive_diversity_data2)
Invasive_diversity_data2$RemRev <- ifelse(Invasive_diversity_data2$iPlantYear == "Remnant", "Remnant", "Revegetation")
Invasive_diversity_data2$RemRev <- factor(Invasive_diversity_data2$RemRev, levels = c("Revegetation", "Remnant"))

unique(Invasive_diversity_data2$tEcoName)
# Invasive_diversity_data2$tEcoName <- factor(Invasive_diversity_data2$tEcoName, levels = c("A.verticillata", "C.gracilis", "E.diversifolia", "E.fasciculosa", "E.porosa", "E.incrassata", "E.leucoxylon", "E.odorata"))

ggplot(Invasive_diversity_data2, aes(x= iPlantYear, y= EffSpp))+
  geom_violin(fill = "#D1B000")+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(Survey_Year ~ tEcoName, scales = "free_x", space = "free_x")+
  theme_test()
# ggsave(filename = "Veg-EffectiveSpecies-PlantYear-INVASIVEONLY.pdf", path = outdir, width = 14, height = 7)

# remrev_cols <- c("Remnant" = "darkgreen", "Revegetation" = "#9DC183")
remrev_cols_yellow <- c("Remnant" = "#B38600", "Revegetation" = "#FFF176")

ggplot(Invasive_diversity_data2, aes(x= RemRev, y= EffSpp))+
  geom_violin(aes(fill = RemRev))+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(Survey_Year ~ tEcoName, scales = "free_x", space = "free_x")+
  scale_fill_manual(values = remrev_cols_yellow)+
  theme_test()
# ggsave(filename = "Veg-EffectiveSpecies-RemRev-INVASIVEONLY.pdf", path = outdir, width = 14, height = 7)

Invasive_diversity_data3 <- Invasive_diversity_data2 %>%
  group_by(Survey_Year, WptID, iPlantYear, Treatment, tEcoName, RemRev) %>%
  reframe(mean_EffSpp =  mean(EffSpp))

ggplot(Invasive_diversity_data3, aes(x= RemRev, y= mean_EffSpp))+
  geom_violin(aes(fill = RemRev))+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(.~Survey_Year, scales = "free_x", space = "free_x")+
  scale_fill_manual(values = remrev_cols_yellow)+
  theme_test()
# ggsave(filename = "Veg-EffectiveSpecies-PlantYear-INVASIVEONLY-survey_year.pdf", path = outdir, width = 7, height = 4)

# LUMPED

ggplot(Invasive_diversity_data2, aes(x= iPlantYear, y= EffSpp))+
  geom_violin(fill = "#D1B000")+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(Survey_Year ~ iLumpEco, scales = "free_x", space = "free_x")+
  theme_test()
# ggsave(filename = "Veg-EffectiveSpecies-PlantYear-INVASIVEONLY-LumpEco.pdf", path = outdir, width = 14, height = 7)

# remrev_cols <- c("Remnant" = "darkgreen", "Revegetation" = "#9DC183")
remrev_cols_yellow <- c("Remnant" = "#B38600", "Revegetation" = "#FFF176")

ggplot(Invasive_diversity_data2, aes(x= RemRev, y= EffSpp))+
  geom_violin(aes(fill = RemRev))+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(Survey_Year ~ iLumpEco, scales = "free_x", space = "free_x")+
  scale_fill_manual(values = remrev_cols_yellow)+
  theme_test()
# ggsave(filename = "Veg-EffectiveSpecies-RemRev-INVASIVEONLY-LumpEco.pdf", path = outdir, width = 14, height = 7)

Invasive_diversity_data3 <- Invasive_diversity_data2 %>%
  group_by(Survey_Year, WptID, iPlantYear, Treatment, tEcoName, RemRev) %>%
  reframe(mean_EffSpp =  mean(EffSpp))

ggplot(Invasive_diversity_data3, aes(x= tEcoName, y= mean_EffSpp))+
  geom_violin(aes(fill = tEcoName))+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(RemRev~Survey_Year, scales = "free_x", space = "free_x")+
  scale_fill_brewer(palette = "Set1")+
  theme_test()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
# ggsave(filename = "Veg-EffectiveSpecies-Ecosystem-INVASIVEONLY-bySurveyRemRev.pdf", path = outdir, width = 7, height = 4)

ggplot(Invasive_diversity_data3, aes(x= tEcoName, y= mean_EffSpp))+
  geom_violin(aes(fill = tEcoName))+
  geom_boxplot(outlier.shape = NA, width = 0.1)+
  geom_jitter(height = 0, width = 0.2, alpha = 0.25) +
  facet_grid(RemRev ~ Survey_Year, scales = "free_x", space = "free_x")+
  scale_fill_brewer(palette = "Set1")+
  theme_test()
# ggsave(filename = "Veg-EffectiveSpecies-Ecosystem-INVASIVEONLY-fac_survey_year.pdf", path = outdir, width = 7, height = 4)

# stats
model1_invasive <- LMEM_permute_anova(data = Invasive_diversity_data2,
                                      response = "EffSpp", fixed_effects = c("RemRev", "Survey_Year", "RemRev:Survey_Year"),
                                      random_effects = c("WptID", "tEcoName"),
                                      nreps = 999)
model1_invasive$formula
# EffSpp ~ RemRev + Survey_Year + RemRev:Survey_Year + (1 | WptID) + (1 | tEcoName)

model1_invasive$Anova
#                      Chisq Df Pr(>Chisq)    
# RemRev             16.8451  1  4.056e-05 ***
# Survey_Year        39.1936  1  3.838e-10 ***
# RemRev:Survey_Year  8.8367  1   0.002952 ** 

model1_invasive$permutation_p_values
# RemRev        Survey_Year RemRev:Survey_Year 
#      0                  0                  0 

model2_invasive <- LMEM_permute_anova(data = Invasive_diversity_data2,
                                      response = "EffSpp", fixed_effects = c("tEcoName", "Survey_Year", "tEcoName:Survey_Year"),
                                      random_effects = c("WptID", "RemRev"),
                                      nreps = 999)
model2_invasive$formula
# EffSpp ~ tEcoName + Survey_Year + tEcoName:Survey_Year + (1 | WptID) + (1 | RemRev)

model2_invasive$Anova
#                        Chisq Df Pr(>Chisq)    
# tEcoName              1.4021  4  0.8438287    
# Survey_Year          39.4967  1  3.286e-10 ***
# tEcoName:Survey_Year 20.3942  4  0.0004174 ***

model2_invasive$permutation_p_values
#  tEcoName          Survey_Year tEcoName:Survey_Year 
# 0.8188188            0.0000000            0.0000000 

# Stacked bar plots ------------------------------------------------------------
# Stacked bar plots to visualise vegetation community differences
all_veg_data_v_dataonly <- all_veg_data_v %>% 
  select(VisitID, Survey_Year, iPlantYear, Treatment, NSXCODE, tLF, INTRODUCED, tEcoName, `SCIENTIFIC NAME`, WptID, VisitID)

all_veg_data_v_dataonly$species_site <- paste0(all_veg_data_v_dataonly$VisitID, "_", all_veg_data_v_dataonly$`SCIENTIFIC NAME`)
all_veg_data_v_dataonly <- all_veg_data_v_dataonly %>% select(-`SCIENTIFIC NAME`, -VisitID)
all_veg_data_v_com$species_site <- paste0(all_veg_data_v_com$VisitID, "_", all_veg_data_v_com$`SCIENTIFIC NAME`)

all_veg_data_v_com_wData <- merge(all_veg_data_v_com, all_veg_data_v_dataonly, by = "species_site")

all_veg_data_v_com_wData$INTRODUCED <- ifelse(is.na(all_veg_data_v_com_wData$INTRODUCED)==TRUE, 0, 1)

veg_data_long <- all_veg_data_v_com_wData
veg_data_long$Wpt_yr <- paste0("G", veg_data_long$WptID, "_y", veg_data_long$Survey_Year)
head(veg_data_long)

# reduce resolution to functional group + introduced vs native species
Lifeform_types <- as.data.frame(read_excel("Joint_Vegetation_data_210126.xlsx", sheet= "Life_Form"))###
nrow(Lifeform_types)
Lifeform_types$Functional_gp <- c("Tree", "Tree", "Tree", "Tree",
                                  "Tree", "Tree",
                                  "Shrub", "Shrub", "Shrub", "Shrub","Shrub", "Shrub",
                                  "Grass", "Grass", "Grass",
                                  "Herbaceous",
                                  "Sedges",
                                  "Sedges",
                                  "Other (Vines, Mosses, Lichens etc.)",
                                  "Other (Vines, Mosses, Lichens etc.)",
                                  "Other (Vines, Mosses, Lichens etc.)",
                                  "Other (Vines, Mosses, Lichens etc.)",
                                  "Other (Vines, Mosses, Lichens etc.)")

veg_data_long2 <- veg_data_long %>%
  left_join(
    Lifeform_types %>%
      select(tLF, Functional_gp),
    by = "tLF"
  )

# Summarise at vegetation level 
# (total cover summarised across taxa within visits)
veg_data_long_LF <- veg_data_long2 %>% 
  select(VisitID, INTRODUCED, Functional_gp, total_cover, Wpt_yr) %>%
  group_by(VisitID, Wpt_yr, INTRODUCED, Functional_gp) %>%
  summarise(total_cover_LF = sum(total_cover, na.rm = FALSE)) %>%
  as.data.frame() 
head(veg_data_long_LF)

colnames(veg_data_long_LF)

# Find averages by functional group (shrub, forb etc,)
# Averages are across visits, within sites
veg_data_long_LF2 <- veg_data_long_LF %>%
  group_by(Wpt_yr, Functional_gp, INTRODUCED) %>%
  summarise(Avg_Cover_LF = mean(total_cover_LF, na.rm = FALSE)) %>%
  as.data.frame()

# A) native or invasive functional group
head(veg_data_long_LF2)
veg_data_long_LF2$tLF_NvI <- paste0(veg_data_long_LF2$Functional_gp, "_", veg_data_long_LF2$INTRODUCED)
veg_data_long_LF2$Native_Invasive_LFs <- ifelse(veg_data_long_LF2$INTRODUCED == 1, 
                                                paste0("Invasive_", veg_data_long_LF2$Functional_gp),
                                                paste0("Native_", veg_data_long_LF2$Functional_gp))

unique(veg_data_long_LF2$tLF_NvI)
unique(veg_data_long_LF2$Native_Invasive_LFs)

# B) Native or invasive visit-level funcitonal groups
head(veg_data_long_LF)
veg_data_long_LF$tLF_NvI <- paste0(veg_data_long_LF$Functional_gp, "_", veg_data_long_LF$INTRODUCED)
veg_data_long_LF$Native_Invasive_LFs <- ifelse(veg_data_long_LF$INTRODUCED == 1, 
                                               paste0("Invasive_", veg_data_long_LF$Functional_gp),
                                               paste0("Native_", veg_data_long_LF$Functional_gp))

unique(veg_data_long_LF$tLF_NvI)
unique(veg_data_long_LF$Native_Invasive_LFs)

# now add remnant vs reveg ecosystems....
colnames(veg_data_long_LF2)
colnames(veg_data_long_LF)
colnames(veg_data_long)

library(stringr)
veg_data_long_LF2$Survey_year <- substr(veg_data_long_LF2$Wpt_yr,
                                        nchar(veg_data_long_LF2$Wpt_yr) - 3,
                                        nchar(veg_data_long_LF2$Wpt_yr))
veg_data_long_LF$Survey_year <- substr(veg_data_long_LF$Wpt_yr,
                                       nchar(veg_data_long_LF2$Wpt_yr) - 3,
                                       nchar(veg_data_long_LF2$Wpt_yr))


wpt_year_lookup <- veg_data_long %>%
  distinct(Wpt_yr, iPlantYear, tEcoName)
wpt_year_lookup$RevRem <- ifelse(wpt_year_lookup$iPlantYear == "Remnant", "Remnant", "Revegetation")

veg_data_long_LF3 <- veg_data_long_LF2 %>%
  left_join(
    wpt_year_lookup,
    by = "Wpt_yr"
  ) %>%
  rename(Treatment = iPlantYear)
veg_data_long_LF1.2 <- veg_data_long_LF %>%
  left_join(
    wpt_year_lookup,
    by = "Wpt_yr"
  ) %>%
  rename(Treatment = iPlantYear)


# Get average per functional group across samples/site
veg_data_long_LF4 <- veg_data_long_LF3 %>%
  group_by(Survey_year, Treatment, tEcoName, Native_Invasive_LFs) %>% 
  summarise(Average_cover = mean(Avg_Cover_LF))

unique(veg_data_long_LF4$Native_Invasive_LFs)

### colour gradient manual ---------
Lifeform_types$Functional_gp
unique(veg_data_long_LF4$Native_Invasive_LFs)

big_palette_8 <- c(
  "Native_Other (Vines, Mosses, Lichens etc.)" = "#e5f5e0",  "Invasive_Other (Vines, Mosses, Lichens etc.)"= "#fee391",
  "Native_Sedges" = "#c7e9c0",                               "Invasive_Sedges"         = "#fec44f",
  "Native_Grass" = "#a1d99b",                                "Invasive_Grass"          = "#fe9929",
  "Native_Herbaceous" = "#74c476",                           "Invasive_Herbaceous"     = "#ec7014",
  "Native_Shrub" = "#41ab5d",                                "Invasive_Shrub"          = "#cc4c02",
  "Native_Tree" = "#006d2c",                                 "Invasive_Tree"           = "#8c2d04")

veg_data_long_LF1.2$Native_Invasive_LFs <- factor(veg_data_long_LF1.2$Native_Invasive_LFs, 
                                                  levels = c("Invasive_Other (Vines, Mosses, Lichens etc.)", "Invasive_Sedges", "Invasive_Grass", "Invasive_Herbaceous", "Invasive_Shrub", "Invasive_Tree",
                                                             "Native_Other (Vines, Mosses, Lichens etc.)", "Native_Sedges", "Native_Grass", "Native_Herbaceous", "Native_Shrub", "Native_Tree"))

veg_data_long_LF4$Native_Invasive_LFs <- factor(veg_data_long_LF4$Native_Invasive_LFs, 
                                                levels = c("Invasive_Other (Vines, Mosses, Lichens etc.)", "Invasive_Sedges", "Invasive_Grass", "Invasive_Herbaceous", "Invasive_Shrub", "Invasive_Tree",
                                                           "Native_Other (Vines, Mosses, Lichens etc.)", "Native_Sedges", "Native_Grass", "Native_Herbaceous", "Native_Shrub", "Native_Tree"))

veg_data_long_LF4$Treatment <- factor(veg_data_long_LF4$Treatment, levels = c("2015", "2014", "2013", "2012", "Remnant"))

## PlantYear Plot graphs ------------------------------------------------------------------

colnames(veg_data_long_LF4) <- c("Survey_year", "Treatment", "Ecosystem", "Native_Invasive_LFs", "Average_cover")
table(veg_data_long_LF4$Ecosystem)
unique(veg_data_long_LF4$Ecosystem)

absol_stacked <- ggplot(veg_data_long_LF4, aes(x = Treatment, y = Average_cover, fill = Native_Invasive_LFs)) + 
  geom_col(colour="black")+
  scale_fill_manual(values = big_palette_8) +
  facet_grid(Survey_year~Ecosystem, scales = "free_x", space = "free_x")+
  labs(y = "Average cover (%)", x = "Planting year")+
  theme_test()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
absol_stacked
# ggsave(filename = "Veg-Functional_stackplot-survey-Ecosystem-LumpEco.pdf", path = outdir, width = 14, height = 7)

head(veg_data_long_LF4)
veg_data_long_RA <- veg_data_long_LF4 %>%
  group_by(Survey_year, Treatment, Ecosystem) %>%
  mutate(
    Relative_cover = Average_cover / sum(Average_cover, na.rm = TRUE)
  ) %>%
  ungroup()
veg_data_long_RA$Treatment <- factor(veg_data_long_RA$Treatment, levels = c("2015", "2014", "2013", "2012", "Remnant"))

relabun_stacked <- ggplot(veg_data_long_RA, aes(x = Treatment, y = Relative_cover, fill = Native_Invasive_LFs)) + 
  geom_col(colour="black")+
  scale_fill_manual(values = big_palette_8) +
  facet_grid(Survey_year ~ Ecosystem, scales = "free_x", space = "free_x")+
  labs(y = "Average cover (%)", x = "Planting year")+
  theme_test()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
relabun_stacked
# ggsave(filename = "Veg-Functional_stackplot-RELABUN-survey-Ecosystem-LumpEco.pdf", path = outdir, width = 14, height = 7)

## Reveg vs Remnant onl
# View(veg_data_long_LF3)

veg_data_long_LF3$RevRem <- ifelse(veg_data_long_LF3$Treatment == "Remnant", "Remnant", "Revegetation")


# Get average per functional group in across samples (change treatment)
veg_data_long_LF4_RevRem <- veg_data_long_LF3 %>%
  group_by(Survey_year, RevRem, tEcoName, Native_Invasive_LFs) %>% 
  summarise(Average_cover = mean(Avg_Cover_LF))

veg_data_long_LF4_RevRem$Native_Invasive_LFs <- factor(veg_data_long_LF4_RevRem$Native_Invasive_LFs, 
                                                       levels = c("Invasive_Other (Vines, Mosses, Lichens etc.)",
                                                                  "Invasive_Sedges",
                                                                  "Invasive_Grass",
                                                                  "Invasive_Herbaceous",
                                                                  "Invasive_Shrub",
                                                                  "Invasive_Tree",
                                                                  "Native_Other (Vines, Mosses, Lichens etc.)",
                                                                  "Native_Sedges",
                                                                  "Native_Grass",
                                                                  "Native_Herbaceous",
                                                                  "Native_Shrub",
                                                                  "Native_Tree"))

## RemRev Plot graphs ------------------------------------------------------------------
veg_data_long_LF4_RevRem$RevRem <- factor(veg_data_long_LF4_RevRem$RevRem, levels = c("Revegetation", "Remnant"))

colnames(veg_data_long_LF4_RevRem) <- c("Survey_year", "RevRem", "Ecosystem", "Native_Invasive_LFs", "Average_cover")
absol_stacked_RR <- ggplot(veg_data_long_LF4_RevRem, aes(x = RevRem, y = Average_cover, fill = Native_Invasive_LFs)) + 
  geom_col(colour="black")+
  scale_fill_manual(values = big_palette_8) +
  facet_grid(Survey_year ~ Ecosystem, scales = "free_x", space = "free_x")+
  labs(y = "Average cover (%)", x = "Treatment")+
  theme_test()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
absol_stacked_RR
# ggsave(filename = "Veg-Functional_stackplot-survey-Ecosystem_RemRev.pdf", path = outdir, width = 11, height = 7)

absol_stacked_RR_swap <- ggplot(veg_data_long_LF4_RevRem, aes(x = Survey_year, y = Average_cover, fill = Native_Invasive_LFs)) + 
  geom_col(colour="black")+
  scale_fill_manual(values = big_palette_8) +
  facet_grid(RevRem ~ Ecosystem, scales = "free_x", space = "free_x")+
  labs(y = "Average cover (%)", x = "Survey Year")+
  theme_test()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
absol_stacked_RR_swap
# ggsave(filename = "Veg-Functional_stackplot-RevRems-Ecosystem_SurveyYear-lumpEco.pdf", path = outdir, width = 11, height = 7)


head(veg_data_long_LF4_RevRem)
veg_data_long_LF4_RevRem_RA <- veg_data_long_LF4_RevRem %>%
  group_by(Survey_year, RevRem, Ecosystem) %>%
  mutate(
    Relative_cover = Average_cover / sum(Average_cover, na.rm = TRUE)
  ) %>%
  ungroup()
veg_data_long_LF4_RevRem_RA$RevRem <- factor(veg_data_long_LF4_RevRem_RA$RevRem , levels = c("Revegetation", "Remnant"))
# veg_data_long_LF4_RevRem_RA$Ecosystem <- factor(veg_data_long_LF4_RevRem_RA$Ecosystem, 
#                                                 levels = c("A.verticillata", "C.gracilis", "E.diversifolia", "E.fasciculosa", 
#                                                            "E.porosa", "E.incrassata", "E.leucoxylon", "E.odorata"))

relabun_stacked_RevRem <- ggplot(veg_data_long_LF4_RevRem_RA, aes(x = RevRem, y = Relative_cover, fill = Native_Invasive_LFs)) + 
  geom_col(colour="black")+
  scale_fill_manual(values = big_palette_8) +
  facet_grid(Survey_year ~ Ecosystem, scales = "free_x", space = "free_x")+
  labs(y = "Average cover (%)", x = "Treatment")+
  theme_test()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
relabun_stacked_RevRem
# ggsave(filename = "Veg-Functional_stackplot-RELABUN-survey-Ecosystem_RemRev.pdf", path = outdir, width = 11, height = 7)

relabun_stacked_RevRem_swap <- ggplot(veg_data_long_LF4_RevRem_RA, aes(x = Survey_year, y = Relative_cover, fill = Native_Invasive_LFs)) + 
  geom_col(colour="black")+
  scale_fill_manual(values = big_palette_8) +
  facet_grid(RevRem ~ Ecosystem, scales = "free_x", space = "free_x")+
  labs(y = "Average cover (%)", x = "Survey year")+
  theme_test()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
relabun_stacked_RevRem_swap
# ggsave(filename = "Veg-Functional_stackplot-RELABUN-RevRems-Ecosystem_SurveyYear.pdf", path = outdir, width = 11, height = 7)

# heatmaps treatlevel
tile_plot <- ggplot(veg_data_long_LF4_RevRem, aes(x = RevRem, y = Native_Invasive_LFs, fill = Average_cover )) + 
  geom_tile()+
  scale_fill_gradient2(low = "#EBCC2A", mid = "#FDFDFD", high = "#3B9AB2") + 
  facet_grid(Survey_year ~ Ecosystem, scales = "free_x", space = "free_x")+
  labs(y = "Native vs Invasive functional group", x = "Treatment")+
  theme_test()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
tile_plot

# Log fold change in functional groups...
remnant_baseline <- veg_data_long_LF4_RevRem %>%
  filter(RevRem == "Remnant") %>%
  select(
    Survey_year,
    Ecosystem,
    Native_Invasive_LFs,
    Remnant_cover = Average_cover)

veg_logratio_av <- veg_data_long_LF4_RevRem %>%
  left_join(
    remnant_baseline,
    by = c("Survey_year", "Ecosystem", "Native_Invasive_LFs")) %>%
  mutate(log_ratio_change = log((Average_cover + 1) / (Remnant_cover + 1)))

unique(veg_logratio_plot$RevRem.x)
ggplot(veg_logratio_av, aes(x = Survey_year, y = Native_Invasive_LFs, fill = log_ratio_change)) +
  geom_tile(colour = "black") +
  facet_grid( ~ Ecosystem, scales="free_x") +
  scale_fill_gradient2(high = "#2166ac", mid = "beige", low = "#b2182b", midpoint = 0, 
                       name = "Log-ratio\ndifference to remnant", na.value="white") +
  labs(y = "Functional group") +
  theme_test() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(veg_logratio_av, aes(x = log_ratio_change, y = Native_Invasive_LFs, fill = )) +
  geom_col(colour = "black", width = 0.8) +
  facet_grid(Survey_year ~ Ecosystem, scales="free_x") +
  # scale_fill_gradient2(high = "#2166ac", mid = "beige", low = "#b2182b", midpoint = 0, 
  #                      name = "Log-ratio\ndifference to remnant", na.value="white") +
  labs(y = "Functional group") +
  theme_test()

## Log fold change differences -------------------------------------------------
veg_data_long_LF1.2$RevRem <- factor(veg_data_long_LF1.2$RevRem, levels = c("Revegetation", "Remnant"))
tile_plot_vis <- ggplot(veg_data_long_LF1.2, aes(x = RevRem, y = Native_Invasive_LFs, fill = total_cover_LF )) + 
  geom_tile()+
  scale_fill_gradient2(low = "#EBCC2A", mid = "#FDFDFD", high = "#3B9AB2") + 
  facet_grid(Survey_year ~ tEcoName, scales = "free_x", space = "free_x")+
  labs(y = "Native vs Invasive functional group", x = "Treatment")+
  theme_test()+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
tile_plot_vis

### Compare Reveg and Remnant within across survey years   ----
colnames(veg_data_long_LF1.2)

results <- list()

for (Year in unique(veg_data_long_LF1.2$Survey_year)) {
  for (Eco in unique(veg_data_long_LF1.2$tEcoName)) {
    for (fun_group in unique(veg_data_long_LF1.2$Native_Invasive_LFs)) {
      
      df_sub <- subset(veg_data_long_LF1.2, Survey_year == Year & tEcoName == Eco & Native_Invasive_LFs == fun_group)
      
      if (nrow(df_sub) == 0) next
      
      df_remnant <- subset(df_sub, RevRem == "Remnant")
      df_reveg   <- subset(df_sub, RevRem == "Revegetation")
      
      if (nrow(df_remnant) == 0 | nrow(df_reveg) == 0) next
      
      obs_lr <- log(
        (mean(df_reveg$total_cover_LF) + 1) /
          (mean(df_remnant$total_cover_LF) + 1)
      )
      
      results[[paste(Year, Eco, fun_group, sep = "_")]] <- data.frame(
        Survey_year = Year,
        Ecosystem = Eco,
        Functional_group = fun_group,
        log_ratio_obs = obs_lr
      )
    }
  }
}

results_df <- do.call(rbind, results)
# results_df

nreps <- 999
results_perm <- list()

for (Year in unique(veg_data_long_LF1.2$Survey_year)) {
  for (Eco in unique(veg_data_long_LF1.2$tEcoName)) {
    for (fun_group in unique(veg_data_long_LF1.2$Native_Invasive_LFs)) {
      
      df_sub <- subset(veg_data_long_LF1.2, Survey_year == Year & tEcoName == Eco & Native_Invasive_LFs == fun_group)
      
      if (nrow(df_sub) == 0) next
      if (length(unique(df_sub$RevRem)) < 2) next
      
      # observed
      obs_lr <- log((mean(df_sub$total_cover_LF[df_sub$RevRem == "Revegetation"]) + 1) /(mean(df_sub$total_cover_LF[df_sub$RevRem == "Remnant"]) + 1))
      
      # permutation null
      perm_lr <- numeric(nreps)
      
      for (i in seq_len(nreps)) {
        perm_labels <- sample(df_sub$RevRem)
        perm_lr[i] <- log((mean(df_sub$total_cover_LF[perm_labels == "Revegetation"]) + 1) /(mean(df_sub$total_cover_LF[perm_labels == "Remnant"]) + 1))
      }
      
      # two-sided p-value
      p_perm <- mean(abs(perm_lr) >= abs(obs_lr))
      
      results_perm[[paste(Year, Eco, fun_group, sep = "_")]] <- data.frame(
        Survey_year = Year,
        Ecosystem = Eco,
        Functional_group = fun_group,
        log_ratio_obs = obs_lr,
        p_perm = p_perm
      )
    }
  }
}


perm_results_df <- do.call(rbind, results_perm)

perm_results_df <- perm_results_df %>%
  mutate(
    sig = case_when(
      p_perm <= 0.001 ~ "***",
      p_perm <= 0.01  ~ "**",
      p_perm <= 0.05  ~ "*",
      TRUE            ~ ""
    )
  )
head(perm_results_df)
perm_results_df$Functional_group <- factor(perm_results_df$Functional_group, 
                                           levels = c("Invasive_Other (Vines, Mosses, Lichens etc.)",
                                                      "Invasive_Sedges",
                                                      "Invasive_Grass",
                                                      "Invasive_Herbaceous",
                                                      "Invasive_Shrub",
                                                      "Invasive_Tree",
                                                      "Native_Other (Vines, Mosses, Lichens etc.)",
                                                      "Native_Sedges",
                                                      "Native_Grass",
                                                      "Native_Herbaceous",
                                                      "Native_Shrub",
                                                      "Native_Tree"))
head(perm_results_df)

# perm_results_df$NatInv <- str_(perm_results_df$Functional_group))
ggplot(perm_results_df, aes(x = Survey_year, y = Functional_group, fill = log_ratio_obs))+
  geom_tile()+
  facet_grid(~Ecosystem)+
  scale_fill_gradient2(high = "#2166ac", mid = "beige", low = "#b2182b", midpoint = 0,
                       name = "Log-ratio\ndifference to remnant", na.value="white") +
  geom_text(aes(label = sig), colour = "black", size = 5, fontface = "bold")+
  theme_test()+
  labs(x = "Survey year", y = "Functional group")
# ggsave(filename = "Veg-Functional_heatmap-RELABUN-Ecosystem_RemRev.pdf", path = outdir, width = 12, height = 7)


##### Combined ecosystem comparisons acros Rem REV ---------
results_RR2 <- list()

for (Year in unique(veg_data_long_LF1.2$Survey_year)) {
  for (fun_group in unique(veg_data_long_LF1.2$Native_Invasive_LFs)) {
    
    df_sub <- subset(veg_data_long_LF1.2, Survey_year == Year & Native_Invasive_LFs == fun_group)
    
    if (nrow(df_sub) == 0) next
    
    df_reveg <- subset(df_sub, RevRem == "Revegetation")
    df_remnant <- subset(df_sub, RevRem == "Remnant")
    
    if (nrow(df_reveg) == 0 | nrow(df_remnant) == 0) next
    
    obs_lr <- log(
      (mean(df_reveg$total_cover_LF) + 1) /
        (mean(df_remnant$total_cover_LF) + 1)
    )
    
    results_RR2[[paste(Year, fun_group, sep = "_")]] <- data.frame(
      Survey_year = Year,
      Functional_group = fun_group,
      log_ratio_obs = obs_lr
    )
  }
}

results_RR2_df <- do.call(rbind, results_RR2)

nreps <- 999
results_RR2_perm <- list()

for (Year in unique(veg_data_long_LF1.2$Survey_year)) {
  for (fun_group in unique(veg_data_long_LF1.2$Native_Invasive_LFs)) {
    
    df_sub <- subset(veg_data_long_LF1.2, Survey_year == Year  & Native_Invasive_LFs == fun_group)
    
    if (nrow(df_sub) == 0) next
    # if (length(unique(df_sub$RemRev)) < 2) next

    # observed
    obs_lr <- log((mean(df_sub$total_cover_LF[df_sub$RevRem == "Revegetation"]) + 1) /(mean(df_sub$total_cover_LF[df_sub$RevRem == "Remnant"]) + 1))
    
    # permutation null
    perm_lr <- numeric(nreps)
    
    for (i in seq_len(nreps)) {
      perm_labels <- sample(df_sub$RevRem)
      perm_lr[i] <- log((mean(df_sub$total_cover_LF[perm_labels == "Revegetation"]) + 1) /(mean(df_sub$total_cover_LF[perm_labels == "Remnant"]) + 1))
    }
    
    # two-sided p-value
    p_perm <- mean(abs(perm_lr) >= abs(obs_lr))
    
    results_RR2_perm[[paste(Year, fun_group, sep = "_")]] <- data.frame(
      Survey_year = Year,
      Functional_group = fun_group,
      log_ratio_obs = obs_lr,
      p_perm = p_perm
    )
  }
}

results_RR2_perm_df <- do.call(rbind, results_RR2_perm)
head(results_RR2_perm_df)

results_RR2_perm_df <- results_RR2_perm_df %>%
  mutate(
    sig = case_when(
      p_perm <= 0.001 ~ "***",
      p_perm <= 0.01  ~ "**",
      p_perm <= 0.05  ~ "*",
      TRUE            ~ ""
    )
  )
head(results_RR2_perm_df)
results_RR2_perm_df$Functional_group <- factor(results_RR2_perm_df$Functional_group, 
                                               levels = c("Invasive_Other (Vines, Mosses, Lichens etc.)", "Invasive_Sedges", "Invasive_Grass", "Invasive_Herbaceous", "Invasive_Shrub", "Invasive_Tree",
                                                          "Native_Other (Vines, Mosses, Lichens etc.)", "Native_Sedges", "Native_Grass", "Native_Herbaceous", "Native_Shrub", "Native_Tree"))
head(results_RR2_perm_df)

results_RR2_perm_df$Survey_year <- factor(results_RR2_perm_df$Survey_year, levels = c("2015", "2025"))

ggplot(results_RR2_perm_df, aes(x = Survey_year, y = Functional_group, fill = log_ratio_obs))+
  geom_tile()+
  scale_fill_gradient2(high = "#2166ac", mid = "beige", low = "#b2182b", midpoint = 0,
                       name = "Log-ratio\ndifference to Remnant", na.value="white") +
  geom_text(aes(label = sig), colour = "black", size = 5, fontface = "bold")+
  theme_test()+
  labs(x = "Survey year", y = "Functional group")

### Compare surveys years within each remnant and reveg treatment   ----
colnames(veg_data_long_LF1.2)

results_sy <- list()

for (RR in unique(veg_data_long_LF1.2$RevRem)) {
  for (Eco in unique(veg_data_long_LF1.2$tEcoName)) {
    for (fun_group in unique(veg_data_long_LF1.2$Native_Invasive_LFs)) {
      
      df_sub <- subset(veg_data_long_LF1.2, RevRem == RR & tEcoName == Eco & Native_Invasive_LFs == fun_group)
      
      if (nrow(df_sub) == 0) next
      
      df_2015 <- subset(df_sub, Survey_year == "2015")
      df_2025   <- subset(df_sub, Survey_year == "2025")
      
      if (nrow(df_2015) == 0 | nrow(df_2025) == 0) next
      
      obs_lr <- log(
        (mean(df_2025$total_cover_LF) + 1) /
          (mean(df_2015$total_cover_LF) + 1)
      )
      
      results_sy[[paste(RR, Eco, fun_group, sep = "_")]] <- data.frame(
        RevRem = RR,
        Ecosystem = Eco,
        Functional_group = fun_group,
        log_ratio_obs = obs_lr
      )
    }
  }
}

results_sy_df <- do.call(rbind, results_sy)
# results_df

nreps <- 999
results_sy_perm <- list()

for (RR in unique(veg_data_long_LF1.2$RevRem)) {
  for (Eco in unique(veg_data_long_LF1.2$tEcoName)) {
    for (fun_group in unique(veg_data_long_LF1.2$Native_Invasive_LFs)) {
      
      df_sub <- subset(veg_data_long_LF1.2, RevRem == RR & tEcoName == Eco & Native_Invasive_LFs == fun_group)
      
      if (nrow(df_sub) == 0) next
      if (length(unique(df_sub$Survey_year)) < 2) next
      
      # observed
      obs_lr <- log((mean(df_sub$total_cover_LF[df_sub$Survey_year == "2025"]) + 1) /(mean(df_sub$total_cover_LF[df_sub$Survey_year == "2015"]) + 1))
      
      # permutation null
      perm_lr <- numeric(nreps)
      
      for (i in seq_len(nreps)) {
        perm_labels <- sample(df_sub$Survey_year)
        perm_lr[i] <- log((mean(df_sub$total_cover_LF[perm_labels == "2025"]) + 1) /(mean(df_sub$total_cover_LF[perm_labels == "2015"]) + 1))
      }
      
      # two-sided p-value
      p_perm <- mean(abs(perm_lr) >= abs(obs_lr))
      
      results_sy_perm[[paste(RR, Eco, fun_group, sep = "_")]] <- data.frame(
        RevRem = RR,
        Ecosystem = Eco,
        Functional_group = fun_group,
        log_ratio_obs = obs_lr,
        p_perm = p_perm
      )
    }
  }
}


results_sy_perm_df <- do.call(rbind, results_sy_perm)

results_sy_perm_df <- results_sy_perm_df %>%
  mutate(
    sig = case_when(
      p_perm <= 0.001 ~ "***",
      p_perm <= 0.01  ~ "**",
      p_perm <= 0.05  ~ "*",
      TRUE            ~ ""
    )
  )
head(results_sy_perm_df)
results_sy_perm_df$Functional_group <- factor(results_sy_perm_df$Functional_group, 
                                              levels = c("Invasive_Other (Vines, Mosses, Lichens etc.)",
                                                         "Invasive_Sedges",
                                                         "Invasive_Grass",
                                                         "Invasive_Herbaceous",
                                                         "Invasive_Shrub",
                                                         "Invasive_Tree",
                                                         "Native_Other (Vines, Mosses, Lichens etc.)",
                                                         "Native_Sedges",
                                                         "Native_Grass",
                                                         "Native_Herbaceous",
                                                         "Native_Shrub",
                                                         "Native_Tree"))
head(results_sy_perm_df)

results_sy_perm_df$RevRem <- factor(results_sy_perm_df$RevRem, levels = c("Revegetation", "Remnant"))
ggplot(results_sy_perm_df, aes(x = RevRem, y = Functional_group, fill = log_ratio_obs))+
  geom_tile()+
  facet_grid(~Ecosystem)+
  scale_fill_gradient2(high = "#2166ac", mid = "beige", low = "#b2182b", midpoint = 0,
                       name = "Log-ratio\ndifference to 2015", na.value="white") +
  geom_text(aes(label = sig), colour = "black", size = 5, fontface = "bold")+
  theme_test()+
  labs(x = "Treatment", y = "Functional group")
# ggsave(filename = "Veg-Functional_heatmap-RELABUN-Ecosystem_SurveyYear.pdf", path = outdir, width = 12, height = 7)

##### Combined ecosystem comparisons acros years ---------
results_sy2 <- list()

for (RR in unique(veg_data_long_LF1.2$RevRem)) {
  for (fun_group in unique(veg_data_long_LF1.2$Native_Invasive_LFs)) {
    
    df_sub <- subset(veg_data_long_LF1.2, RevRem == RR & Native_Invasive_LFs == fun_group)
    
    if (nrow(df_sub) == 0) next
    
    df_2015 <- subset(df_sub, Survey_year == "2015")
    df_2025 <- subset(df_sub, Survey_year == "2025")
    
    if (nrow(df_2015) == 0 | nrow(df_2025) == 0) next
    
    obs_lr <- log(
      (mean(df_2025$total_cover_LF) + 1) /
        (mean(df_2015$total_cover_LF) + 1)
    )
    
    results_sy2[[paste(RR, fun_group, sep = "_")]] <- data.frame(
      RevRem = RR,
      Functional_group = fun_group,
      log_ratio_obs = obs_lr
    )
  }
  # }
}

results_sy2_df <- do.call(rbind, results_sy2)
# results_sy2

nreps <- 999
results_sy2_perm <- list()

for (RR in unique(veg_data_long_LF1.2$RevRem)) {
  for (fun_group in unique(veg_data_long_LF1.2$Native_Invasive_LFs)) {
    
    df_sub <- subset(veg_data_long_LF1.2, RevRem == RR & Native_Invasive_LFs == fun_group)
    
    if (nrow(df_sub) == 0) next
    if (length(unique(df_sub$Survey_year)) < 2) next
    
    # observed
    obs_lr <- log((mean(df_sub$total_cover_LF[df_sub$Survey_year == "2025"]) + 1) /(mean(df_sub$total_cover_LF[df_sub$Survey_year == "2015"]) + 1))
    
    # permutation null
    perm_lr <- numeric(nreps)
    
    for (i in seq_len(nreps)) {
      perm_labels <- sample(df_sub$Survey_year)
      perm_lr[i] <- log((mean(df_sub$total_cover_LF[perm_labels == "2025"]) + 1) /(mean(df_sub$total_cover_LF[perm_labels == "2015"]) + 1))
    }
    
    # two-sided p-value
    p_perm <- mean(abs(perm_lr) >= abs(obs_lr))
    
    results_sy2_perm[[paste(RR, fun_group, sep = "_")]] <- data.frame(
      RevRem = RR,
      Functional_group = fun_group,
      log_ratio_obs = obs_lr,
      p_perm = p_perm
    )
  }
}


results_sy2_perm_df <- do.call(rbind, results_sy2_perm)

results_sy2_perm_df <- results_sy2_perm_df %>%
  mutate(
    sig = case_when(
      p_perm <= 0.001 ~ "***",
      p_perm <= 0.01  ~ "**",
      p_perm <= 0.05  ~ "*",
      TRUE            ~ ""
    )
  )
head(results_sy2_perm_df)
results_sy2_perm_df$Functional_group <- factor(results_sy2_perm_df$Functional_group, 
                                               levels = c("Invasive_Other (Vines, Mosses, Lichens etc.)",
                                                          "Invasive_Sedges",
                                                          "Invasive_Grass",
                                                          "Invasive_Herbaceous",
                                                          "Invasive_Shrub",
                                                          "Invasive_Tree",
                                                          "Native_Other (Vines, Mosses, Lichens etc.)",
                                                          "Native_Sedges",
                                                          "Native_Grass",
                                                          "Native_Herbaceous",
                                                          "Native_Shrub",
                                                          "Native_Tree"))
head(results_sy2_perm_df)

results_sy2_perm_df$RevRem <- factor(results_sy2_perm_df$RevRem, levels = c("Revegetation", "Remnant"))
ggplot(results_sy2_perm_df, aes(x = RevRem, y = Functional_group, fill = log_ratio_obs))+
  geom_tile()+
  scale_fill_gradient2(high = "#2166ac", mid = "beige", low = "#b2182b", midpoint = 0,
                       name = "Log-ratio\ndifference to 2015", na.value="white") +
  geom_text(aes(label = sig), colour = "black", size = 5, fontface = "bold")+
  theme_test()+
  labs(x = "Treatment", y = "Functional group")

### Review plots -----
ggpubr::ggarrange(ggplot(perm_results_df, aes(x = Survey_year, y = Functional_group, fill = log_ratio_obs))+
                    geom_tile()+
                    facet_grid(~Ecosystem)+
                    scale_fill_gradient2(high = "#2166ac", mid = "beige", low = "#b2182b", midpoint = 0,
                                         name = "Log-ratio\ndifference to remnant", na.value="white") +
                    geom_text(aes(label = sig), colour = "black", size = 5, fontface = "bold")+
                    theme_test()+
                    labs(x = "Survey year", y = "Functional group"),
                  ggplot(results_sy_perm_df, aes(x = RevRem, y = Functional_group, fill = log_ratio_obs))+
                    geom_tile()+
                    facet_grid(~Ecosystem)+
                    scale_fill_gradient2(high = "#2166ac", mid = "beige", low = "#b2182b", midpoint = 0,
                                         name = "Log-ratio\ndifference to 2015", na.value="white") +
                    geom_text(aes(label = sig), colour = "black", size = 5, fontface = "bold")+
                    theme_test()+
                    labs(x = "Treatment", y = "Functional group"),
                  align= "hv")

results_sy2_perm_df$Ecosystem <- "Combined vegetation communities"
min(results_sy2_perm_df$log_ratio_obs);max(results_sy2_perm_df$log_ratio_obs)
min(results_sy_perm_df$log_ratio_obs);max(results_sy_perm_df$log_ratio_obs)
# replace Invasive with "non-native"
results_sy_perm_df$Functional_group <- gsub("^Invasive", "Non-native", results_sy_perm_df$Functional_group)
results_sy_perm_df$Functional_group <- gsub("_", " ", results_sy_perm_df$Functional_group)
results_sy2_perm_df$Functional_group <- gsub("^Invasive", "Non-native", results_sy2_perm_df$Functional_group)
results_sy2_perm_df$Functional_group <- gsub("_", " ", results_sy2_perm_df$Functional_group)

results_sy2_perm_df$Functional_group <- factor(results_sy2_perm_df$Functional_group, levels = c("Non-native Other (Vines, Mosses, Lichens etc.)", "Non-native Sedges", "Non-native Grass", "Non-native Herbaceous","Non-native Shrub", "Non-native Tree",
                                                                                                "Native Other (Vines, Mosses, Lichens etc.)", "Native Sedges", "Native Grass", "Native Herbaceous", "Native Shrub", "Native Tree"))
results_sy_perm_df$Functional_group <- factor(results_sy_perm_df$Functional_group, levels = c("Non-native Other (Vines, Mosses, Lichens etc.)", "Non-native Sedges", "Non-native Grass", "Non-native Herbaceous", "Non-native Shrub", "Non-native Tree",
                                                                                                "Native Other (Vines, Mosses, Lichens etc.)", "Native Sedges", "Native Grass", "Native Herbaceous", "Native Shrub", "Native Tree"))


results_sy_perm_df$RevRem <- factor(results_sy_perm_df$RevRem, levels = c("Revegetation", "Remnant"))
results_sy2_perm_df$RevRem <- factor(results_sy2_perm_df$RevRem, levels = c("Revegetation", "Remnant"))

ggpubr::ggarrange(ggplot(results_sy_perm_df, aes(x = RevRem, y = Functional_group, fill = log_ratio_obs))+
                    geom_tile()+
                    facet_grid(~Ecosystem)+
                    scale_fill_gradient2(high = "#2166ac", mid = "beige", low = "#b2182b", midpoint = 0,
                                         name = "Log-ratio\ndifference to 2015", na.value="white",
                                         breaks = c(-3, -2, -1, 0, 1, 2, 3),
                                         limits = c(-3, 3)) +
                    geom_text(aes(label = sig), colour = "black", size = 5, fontface = "bold")+
                    theme_test()+
                    labs(x = "", y = ""),
                  ggplot(results_sy2_perm_df, aes(x = RevRem, y = Functional_group, fill = log_ratio_obs))+
                    geom_tile()+
                    facet_grid(~Ecosystem)+
                    scale_fill_gradient2(high = "#2166ac", mid = "beige", low = "#b2182b", midpoint = 0,
                                         name = "Log-ratio\ndifference to 2015", na.value="white",
                    breaks = c(-3, -2, -1, 0, 1, 2, 3), 
                    limits = c(-3, 3)) +
                    geom_text(aes(label = sig), colour = "black", size = 5, fontface = "bold")+
                    theme_test() +
                    theme(panel.border = element_rect(colour = "black", fill = NA, size = 2),  # thick border
                          axis.text.y = element_blank()  
                    )+
                    labs(x = "", y = ""),
                  align= "hv", common.legend = TRUE, widths = c(5, 1), legend = "bottom")
# ggsave(filename = "Veg-Functional_heatmap-ALL_Ecosystem_SurveyYear.pdf", path = outdir, width = 15, height = 7)

results_RR2_perm_df$Ecosystem <- "Combined vegetation communities"
min(perm_results_df$log_ratio_obs); max(perm_results_df$log_ratio_obs)
min(results_RR2_perm_df$log_ratio_obs); max(results_RR2_perm_df$log_ratio_obs)

# replace Invasive with "non-native"
perm_results_df$Functional_group <- gsub("^Invasive", "Non-native", perm_results_df$Functional_group)
perm_results_df$Functional_group <- gsub("_", " ", perm_results_df$Functional_group)
results_RR2_perm_df$Functional_group <- gsub("^Invasive", "Non-native", results_RR2_perm_df$Functional_group)
results_RR2_perm_df$Functional_group <- gsub("_", " ", results_RR2_perm_df$Functional_group)

perm_results_df$Functional_group <- factor(perm_results_df$Functional_group, levels = c("Non-native Other (Vines, Mosses, Lichens etc.)", "Non-native Sedges", "Non-native Grass", "Non-native Herbaceous", "Non-native Shrub", "Non-native Tree",
                                                                                              "Native Other (Vines, Mosses, Lichens etc.)", "Native Sedges", "Native Grass", "Native Herbaceous", "Native Shrub", "Native Tree"))
results_RR2_perm_df$Functional_group <- factor(results_RR2_perm_df$Functional_group, levels = c("Non-native Other (Vines, Mosses, Lichens etc.)", "Non-native Sedges", "Non-native Grass", "Non-native Herbaceous","Non-native Shrub", "Non-native Tree",
                                                                                                "Native Other (Vines, Mosses, Lichens etc.)", "Native Sedges", "Native Grass", "Native Herbaceous", "Native Shrub", "Native Tree"))

perm_results_df$Survey_year <- factor(perm_results_df$Survey_year , levels = c("2015", "2025"))
results_RR2_perm_df$Survey_year <- factor(results_RR2_perm_df$Survey_year, levels = c("2015", "2025"))

ggpubr::ggarrange(ggplot(perm_results_df, aes(x = Survey_year, y = Functional_group, fill = log_ratio_obs))+
                    geom_tile()+
                    facet_grid(~Ecosystem)+
                    scale_fill_gradient2(high = "#2166ac", mid = "beige", low = "#b2182b", midpoint = 0,
                                         name = "Log-ratio\ndifference to Remnant", na.value="white",
                                         breaks = c(-3, -2, -1, 0, 1, 2, 3),
                                         limits = c(-3, 3)) +
                    geom_text(aes(label = sig), colour = "black", size = 5, fontface = "bold")+
                    theme_test()+
                    labs(x = "", y = ""),
                  ggplot(results_RR2_perm_df, aes(x = Survey_year, y = Functional_group, fill = log_ratio_obs))+
                    geom_tile()+
                    facet_grid(~Ecosystem)+
                    scale_fill_gradient2(high = "#2166ac", mid = "beige", low = "#b2182b", midpoint = 0,
                                         name = "Log-ratio\ndifference to Remnant", na.value="white",
                                         breaks = c(-3, -2, -1, 0, 1, 2, 3), 
                                         limits = c(-3, 3)) +
                    geom_text(aes(label = sig), colour = "black", size = 5, fontface = "bold")+
                    theme_test() +
                    theme(panel.border = element_rect(colour = "black", fill = NA, size = 2),  # thick border
                          axis.text.y = element_blank()  
                    )+
                    labs(x = "", y = ""),
                  align= "hv", common.legend = TRUE, widths = c(5, 1), legend = "bottom")
# ggsave(filename = "Veg-Functional_heatmap-ALL_Ecosystem_Treatment_RR.pdf", path = outdir, width = 15, height = 7)

### Proportion plants by management --------------------------------------------
head(veg_data_long_LF)
head(veg_data_long_LF2)
head(veg_data_long_LF3)
head(veg_data_long_LF4)
head(management_data)

# Seperate by year
veg_data_long_LF_2025 <- subset(veg_data_long_LF, Survey_year == "2025")
# get WPTID from datafames as a column
veg_data_long_LF_2025 <- veg_data_long_LF_2025 %>%
  separate(col = Wpt_yr,
           into = c("WptID", "Year"),
           sep = "_y")

management_data$WptID <- paste0("G", management_data$WptID)
# Join with management_data
LF_M <- left_join(veg_data_long_LF_2025, management_data, by = "WptID")
unique(LF_M$grazing_level)

LF_M <- LF_M %>%
  separate(col = Native_Invasive_LFs,
           into = c("Status", "Functional Group"),
           sep = "_")

LF_M$grazing_level <- ifelse(is.na(LF_M$grazing_level) == TRUE, "Not observed", LF_M$grazing_level)
LF_M$grazing_level <- factor(LF_M$grazing_level, levels = c("Not observed", "Low", "Moderate", "High"))
LF_M$Status <- ifelse(LF_M$Status == "Invasive", "Non-native", "Native")

LF_M$`Functional Group` <- factor(LF_M$`Functional Group`, levels = c("Other (Vines, Mosses, Lichens etc.)", "Sedges", "Herbaceous", "Grass", "Shrub", "Tree"))

# Plot functional cover by grazing level
ggplot(LF_M, aes(x = grazing_level, y = total_cover_LF, fill= grazing_level))+
  geom_violin()+
  geom_boxplot(outlier.shape = NA, width = 0.1, fill= "white")+
  # geom_jitter(width = 0.1, alpha = 0.5)+
  labs(x ="Observed evidence of grazing", y = "Species agglomerated cover")+
  facet_grid(Status~`Functional Group`)+
  # facet_grid(~`Functional Group`)+
  theme_test()

ggplot(LF_M, aes(x = grazing_level, y = `Functional Group`, fill= total_cover_LF))+
  geom_tile()+
  scale_fill_gradient(low = "orange", high = "red")+
  facet_grid(~Status)+
  theme_test()

# Beta Diversity ---------------------------------------------------------------
# Here i need to get average cover/abundance for each site, calcualted across all 
# visitIDs

library(vegan)

# veg_species_wide
# Native_veg_wide 
# Invasive_veg_wide

all_veg_data_v_com_beta <- all_veg_data_v %>%
  distinct() %>% #make sure we don't have duplicates
  select(VisitID, Site_ID, WptID, `SCIENTIFIC NAME`, cover_midpoint) %>%
  group_by(VisitID, Site_ID, WptID, `SCIENTIFIC NAME`) %>%
  summarise(total_cover = sum(cover_midpoint, na.rm = TRUE)) %>%
  as.data.frame() %>%
  select(Site_ID, WptID, `SCIENTIFIC NAME`, total_cover) %>%
  group_by(Site_ID, WptID, `SCIENTIFIC NAME`) %>%
  summarise(mean_cover = mean(total_cover, na.rm = TRUE)) %>%
  as.data.frame()

nrow(all_veg_data_v_com_beta)
head(all_veg_data_v_com_beta)

veg_species_wide_beta <- all_veg_data_v_com_beta %>%
  pivot_wider(
    id_cols = Site_ID,             
    names_from = `SCIENTIFIC NAME`,
    values_from = mean_cover,
    values_fill = 0) %>% 
  as.data.frame()

rownames(veg_species_wide_beta) <- veg_species_wide_beta$Site_ID
veg_species_wide_beta <- veg_species_wide_beta %>% select(-Site_ID)
dim(veg_species_wide_beta) # 149 686

# View(veg_species_wide_beta)
# save Veg wide site data_frame
# saveRDS(object = veg_species_wide_beta, file = "veg_species_wide_beta_nw.RDS")

set.seed(123)
NMDS_bray <- metaMDS(veg_species_wide_beta, distance = "bray")

NMDS_bray$stress # 0.1944129
plot(NMDS_bray)

points_bray <- NMDS_bray$points

points_bray <- as.data.frame(points_bray)
points_bray$Site_ID <- rownames(points_bray)

# Attach metadata
all_veg_site_beta <- all_veg_data_v %>% 
  select(Site_ID, Survey_Year, WptID, iPlantYear, Treatment, iLumpEco, HILG) %>%
  distinct()
head(all_veg_site_beta)
head(points_bray)

metadata_s <- left_join(all_veg_site_beta, points_bray, by = "Site_ID")
metadata_s$Survey_Year <- as.factor(metadata_s$Survey_Year)

colnames(metadata_s)
# "Site_ID"      "Survey_Year"  "iWptID"       "iPlantYear"   "iTreatID"     "iEcosystemID"
# "tEcoName"     "MDS1"         "MDS2"

# HULLs
metadata_s$RemRev <- ifelse(metadata_s$iPlantYear == "Remnant", "Remnant", "Revegetation")

data_for_hulls <- metadata_s %>%
  select(WptID, MDS1, MDS2, iLumpEco, RemRev, iPlantYear, Survey_Year)

data_for_hulls$hull_combs <- paste0(data_for_hulls$RemRev, "_", data_for_hulls$Survey_Year)
data_for_hulls$hull_lumps <- paste0(data_for_hulls$iLumpEco, "_", data_for_hulls$Survey_Year)
metadata_s$hull_combs <- paste0(metadata_s$RemRev, "_", metadata_s$Survey_Year)
metadata_s$hull_lumps <- paste0(metadata_s$iLumpEco, "_", metadata_s$Survey_Year)

find_hull <- function(data_for_hulls) data_for_hulls[chull(data_for_hulls$MDS1, data_for_hulls$MDS2), ]
hulls_RemRev <- plyr::ddply(data_for_hulls, "hull_combs", find_hull)
hulls_iLumpEco <- plyr::ddply(data_for_hulls, "hull_lumps", find_hull)
hulls_iPlantYear <- plyr::ddply(data_for_hulls, "iPlantYear", find_hull)

# HILG sites
ggplot(metadata_s, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =HILG))+
  geom_path(aes(group=WptID), alpha=0.75, colour = "grey")+
  geom_point(size= 2.5, colour = "black")+
  geom_point(size= 3.5, colour = "black")+
  geom_point(size= 3)+
  labs(x = "NMDS1", y = "NMDS2")+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-SurveyYear_HILG-path.pdf", path = outdir, width = 7, height = 6)


# management plot
colnames(grazing_info_only)
metadata_s <- left_join(metadata_s, grazing_info_only, by = "WptID")

ggplot(metadata_s, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =grazing_level))+
  geom_path(aes(group=WptID), alpha=0.75, colour = "grey")+
  geom_point(size= 2.5, colour = "black")+
  geom_point(size= 3.5, colour = "black")+
  geom_point(size= 3)+
  labs(x = "NMDS1", y = "NMDS2")+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()

ggplot(metadata_s[metadata_s$RemRev == "Revegetation" & metadata_s$Survey_Year == "2025",], aes(x = MDS1, y = MDS2, colour =grazing_level))+
  geom_path(aes(group=WptID), alpha=0.75, colour = "grey")+
  geom_point(size= 2.5, colour = "black")+
  geom_point(size= 3.5, colour = "black")+
  geom_point(size= 3)+
  labs(x = "NMDS1", y = "NMDS2", title = "Grazing detected (2025 - revegetation)")+
  # scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  scale_colour_manual(values =  c("High" = "#d73027", "Moderate" = "#ffffbf", "Low" = "#4575b4", "Not observed" = "lightgrey"))+
  theme_test()
ggsave(filename = "NMDS-reveg_only-grazing intensity.pdf", path = outdir)

# Remnant reveg survey year 
unique(hulls_RemRev$hull_combs)
ggplot(metadata_s, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =hull_combs))+
  geom_path(aes(group=WptID), alpha=0.75, colour = "grey")+
  geom_polygon(data=hulls_RemRev, aes(x = MDS1, y = MDS2, group = hull_combs, fill = hull_combs), alpha=0.2, colour = "black", linewidth=0.2)+  
  scale_colour_manual(values =  c("Remnant_2015" = "#a1d99b", "Revegetation_2015" = "#fe9929", "Remnant_2025" = "#006d2c", "Revegetation_2025" = "#cc4c02"))+ #"#006d2c" "#8c2d04" "#cc4c02",
  scale_fill_manual(values =  c("Remnant_2015" = "#a1d99b", "Revegetation_2015" = "#fe9929", "Remnant_2025" = "#006d2c", "Revegetation_2025" = "#cc4c02"))+ #"#006d2c" "#8c2d04" "#cc4c02",
  geom_point(size= 2.5, colour = "black")+
  geom_point(size= 3.5, colour = "black")+
  geom_point(size= 3)+
  labs(x = "NMDS1", y = "NMDS2")+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-SurveyYear_RemRev_HULLs.pdf", path = outdir, width = 7, height = 6)

ggplot(metadata_s, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =RemRev))+
  geom_path(aes(group=WptID), alpha=0.75, colour = "grey")+
  scale_colour_manual(values =  c("Remnant" = "#006d2c", "Revegetation" = "#cc4c02"))+
  scale_fill_manual(values =  c("Remnant" = "#006d2c", "Revegetation" = "#cc4c02"))+
  geom_point(size= 2.5, colour = "black")+
  geom_point(size= 3.5, colour = "black")+
  geom_point(size= 3)+
  labs(x = "NMDS1", y = "NMDS2")+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-SurveyYear_RemRev-connected.pdf", path = outdir, width = 7, height = 6)

# reattach ecosystem names
metadata_s <- left_join(metadata_s, ilumpDF, by = "iLumpEco")
colnames(metadata_s)
ggplot(metadata_s, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =tEcoName))+
  geom_path(aes(group=WptID), alpha=0.75, colour = "grey")+
  geom_point(size= 2.5, colour = "black")+
  geom_point(size= 3.5, colour = "black")+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  scale_colour_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  labs(x = "NMDS1", y = "NMDS2")+
  geom_point(size= 3)+
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-SurveyYear_Ecosystem-connected.pdf", path = outdir, width = 7, height = 6)

# unique(hulls_iLumpEco$hull_lumps)
# RColorBrewer::brewer.pal(n = 10, palette="Paired")
# ggplot(metadata_s, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =as.factor(hull_lumps)))+
#   geom_polygon(data=hulls_iLumpEco, aes(x = MDS1, y = MDS2, group = hull_lumps, fill = hull_lumps), alpha=0.2, colour = "black", linewidth=0.2)+  
#   geom_point(size= 2.5, colour = "black")+
#   geom_point(size= 3.5, colour = "black")+
#   
#   # scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
#   # scale_fill_manual(values = c("1_2015", "10_2015", "4_2015", "5_2015", "6_2015",
#   #                              "1_2025", "10_2025", "4_2025", "5_2025", "6_2025" ))
#   scale_colour_brewer(palette = "Paired")+
#   scale_fill_brewer(palette = "Set1")+
#   geom_point(size= 3)+
#   theme_test()
# # ggsave(filename = "Veg-BETA-NMDS-SurveyYear_Ecosystem-LumpEco.pdf", path = outdir, width = 7, height = 6)

ggplot(metadata_s, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =iPlantYear))+
  geom_path(aes(group=WptID), alpha=0.75, colour = "grey")+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  scale_colour_brewer(palette = "Set1")+
  geom_point(size= 2.5, colour = "black")+
  geom_point(size= 3.5, colour = "black")+
  labs(x = "NMDS1", y = "NMDS2")+
  geom_point(size= 3)+
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-SurveyYear_PlantYear.pdf", path = outdir, width = 7, height = 6)

# Jaccard
veg_species_PA <- ifelse(veg_species_wide_beta == 0, 0, 1) 

set.seed(123)
NMDS_Jacc <- metaMDS(veg_species_PA, distance = "jaccard", binary = TRUE)
NMDS_Jacc$stress # 0.2035536

plot(NMDS_Jacc)
plot(NMDS_bray)

points_Jacc <- NMDS_Jacc$points

points_Jacc <- as.data.frame(points_Jacc)
points_Jacc$Site_ID <- rownames(points_Jacc)

# Attach metadata
all_veg_site_beta <- all_veg_data_v %>% 
  select(Site_ID, Survey_Year, WptID, iPlantYear, Treatment, tEcoName) %>%
  distinct()
head(all_veg_site_beta)
head(points_Jacc)

metadata_s_Jacc <- left_join(all_veg_site_beta, points_Jacc, by = "Site_ID")
metadata_s_Jacc$Survey_Year <- as.factor(metadata_s_Jacc$Survey_Year)

colnames(metadata_s_Jacc)
# "Site_ID"      "Survey_Year"  "iWptID"       "iPlantYear"   "iTreatID"     "iEcosystemID"
# "tEcoName"     "MDS1"         "MDS2"

metadata_s_Jacc$RemRev <- ifelse(metadata_s_Jacc$iPlantYear == "Remnant", "Remnant", "Revegetation")


ggplot(metadata_s_Jacc, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =RemRev))+
  geom_path(aes(group=WptID), alpha=0.75, colour = "grey")+
  scale_colour_manual(values =  c("Remnant" = "#006d2c", "Revegetation" = "#cc4c02"))+
  scale_fill_manual(values =  c("Remnant" = "#006d2c", "Revegetation" = "#cc4c02"))+
  geom_point(size= 2.5, colour = "black")+
  geom_point(size= 3.5, colour = "black")+
  geom_point(size= 3)+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-Jaccard-SurveyYear_RemRev-connected.pdf", path = outdir, width = 7, height = 6)

ggplot(metadata_s_Jacc, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =iPlantYear))+
  geom_path(aes(group=WptID), alpha=0.75, colour = "grey")+
  scale_colour_brewer(palette = "Set1")+
  geom_point(size= 2.5, colour = "black")+
  geom_point(size= 3.5, colour = "black")+
  geom_point(size= 3)+  
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-Jaccard-SurveyYear_PlantYear.pdf", path = outdir, width = 7, height = 6)

ggplot(metadata_s_Jacc, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =tEcoName))+
  geom_path(aes(group=WptID), alpha=0.75, colour = "grey")+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  scale_colour_brewer(palette = "Set1")+
  geom_point(size= 2.5, colour = "black")+
  geom_point(size= 3.5, colour = "black")+
  geom_point(size= 3)+  
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-Jaccard-SurveyYear_Ecosystem-Connected.pdf", path = outdir, width = 7, height = 6)

# stats
dist_mat_bray_all <- vegdist(veg_species_wide_beta, method = "bray")
# labels(dist_mat_bray_all)
# dist_mat_bray_all_revonly <- vegdist(veg_species_wide_beta[, method = "bray")
# labels(dist_mat_bray_all)

# I need to control for two variables in my data: ecosystem type (tEcoName), and repeated sites (iWptID)
metadata_s$strata_ecosys_wpt <- interaction(
  metadata_s$iLumpEco,
  metadata_s$WptID,
  drop = TRUE
)

set.seed(123)
adonis2(dist_mat_bray_all ~ Survey_Year * RemRev, data = metadata_s,
        strata = metadata_s$strata_ecosys_wpt)
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Blocks:  strata 
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist_mat_bray_all ~ Survey_Year * RemRev, data = metadata_s, strata = metadata_s$strata_ecosys_wpt)
#                     Df SumOfSqs      R2       F Pr(>F)    
# Survey_Year          1    1.789 0.03406 5.2741  0.001 ***
# RemRev               1    0.785 0.01494 2.3138  0.003 ** 
# Survey_Year:RemRev   1    0.771 0.01468 2.2736  0.022 *  
# Residual           145   49.192 0.93632                  
# Total              148   52.537 1.00000    

# Interactions are tricky for beta disperson
metadata_s$Year_RemRev_interaction <- interaction(
  metadata_s$Survey_Year,
  metadata_s$RemRev,
  drop = TRUE
)

bd_int.Y_RR <- betadisper(dist_mat_bray_all, group = metadata_s$Year_RemRev_interaction)
permutest(bd_int.Y_RR, permutations = 999)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# Response: Distances
#            Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups      3 0.05020 0.0167343 6.7304    999  0.001 ***
# Residuals 145 0.36053 0.0024864        

set.seed(123)
adonis2(dist_mat_bray_all ~ Survey_Year * iPlantYear, data = metadata_s,
        strata = metadata_s$strata_ecosys_wpt)
#                         Df SumOfSqs      R2      F Pr(>F)    
# Survey_Year              1    1.789 0.03406 5.2765  0.001 ***
# iPlantYear               4    1.860 0.03540 1.3712  0.001 ***
# Survey_Year:iPlantYear   4    1.753 0.03337 1.2924  0.085 .  
# Residual               139   47.135 0.89717                  
# Total                  148   52.537 1.00000     

# Interactions are tricky for beta disperson
metadata_s$Year_PlantY_interaction <- interaction(
  metadata_s$Survey_Year,
  metadata_s$iPlantYear,
  drop = TRUE
)
bd_int.Y_PY <- betadisper(dist_mat_bray_all, group = metadata_s$Year_PlantY_interaction)
set.seed(123)
permutest(bd_int.Y_PY, permutations = 999)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# Response: Distances
#            Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups      9 0.03957 0.0043968 1.3234    999  0.219
# Residuals 139 0.46181 0.0033223       

# How much variation is due to ecosystem differences?
# AND How much is uniquely attributable to restoration or year effects after accounting for ecosystem?
set.seed(123)
adonis2(dist_mat_bray_all ~ iLumpEco + Survey_Year * RemRev,
        data = metadata_s,
        by   = "terms")
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = dist_mat_bray_all ~ tEcoName + Survey_Year * RemRev, data = metadata_s, by = "terms")
#                     Df SumOfSqs      R2       F Pr(>F)    
# tEcoName        7    2.808 0.05345 1.1912  0.024 *  
# Survey_Year          1    1.799 0.03425 5.3432  0.001 ***
# RemRev               1    0.690 0.01313 2.0479  0.004 **
# Survey_Year:RemRev   1    0.769 0.01463 2.2827  0.002 ** 
# Residual           138   46.472 0.88455                  
# Total              148   52.537 1.00000       

# betadisperson checks
bd_int.Y_PY <- betadisper(dist_mat_bray_all, group = metadata_s$iLumpEco)
permutest(bd_int.Y_PY, permutations = 999)
#            Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)    
# Groups      7 0.14380 0.0205434 7.0869    999  0.001 ***
# Residuals 141 0.40873 0.0028988         

# As calculated prior:
bd_int.Y_RR
permutest(bd_int.Y_RR, permutations = 999)
#            Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups      3 0.05020 0.0167343 6.7304    999  0.001 ***
# Residuals 145 0.36053 0.0024864          

##### distance to 2015 plot -----
dim(dist_mat_bray_all)
labels(dist_mat_bray_all)
dist_mat <- as.matrix(dist_mat_bray_all)

bray_2015_2025 <- metadata_s %>%
  select(Site_ID, WptID, Survey_Year) %>%
  # mutate(Survey_Year = Survey_Year) %>%
  arrange(WptID, Survey_Year) %>%
  group_by(WptID) %>%
  filter(n() == 2) %>%   # keep only sites with both years
  summarise(
    bray_dist = dist_mat[Site_ID[1], Site_ID[2]],
    .groups = "drop"
  )
summary(bray_2015_2025$bray_dist)
range(bray_2015_2025$bray_dist)

bray_2015_2025 <- bray_2015_2025 %>%
  left_join(
    metadata_s %>% distinct(WptID, iLumpEco, RemRev),
    by = "WptID"
  )
ggplot(bray_2015_2025, aes(x = RemRev, y = bray_dist, fill = as.factor(iLumpEco))) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width=0.1, fill = "white") +
  geom_jitter(width = 0.1, alpha = 0.5) +
  theme_test() +
  labs(y = "2015 - 2025 site distances (Bray–Curtis)",
       x = "Ecosystem")+
  facet_grid(.~as.factor(iLumpEco))

# Management (reveg only)

library(vegan)

# veg_species_wide
# Native_veg_wide 
# Invasive_veg_wide

all_veg_data_v_com_beta_reveg <- all_veg_data_v[all_veg_data_v$iPlantYear != "Remnant" & all_veg_data_v$Survey_Year == "2025" ,] %>%
  distinct() %>% #make sure we don't have duplicates
  select(VisitID, Site_ID, WptID, `SCIENTIFIC NAME`, cover_midpoint) %>%
  group_by(VisitID, Site_ID, WptID, `SCIENTIFIC NAME`) %>%
  summarise(total_cover = sum(cover_midpoint, na.rm = TRUE)) %>%
  as.data.frame() %>%
  select(Site_ID, WptID, `SCIENTIFIC NAME`, total_cover) %>%
  group_by(Site_ID, WptID, `SCIENTIFIC NAME`) %>%
  summarise(mean_cover = mean(total_cover, na.rm = TRUE)) %>%
  as.data.frame()

nrow(all_veg_data_v_com_beta_reveg)
head(all_veg_data_v_com_beta_reveg)

veg_species_wide_beta_reveg <- all_veg_data_v_com_beta_reveg %>%
  pivot_wider(
    id_cols = Site_ID,             
    names_from = `SCIENTIFIC NAME`,
    values_from = mean_cover,
    values_fill = 0) %>% 
  as.data.frame()

rownames(veg_species_wide_beta_reveg) <- veg_species_wide_beta_reveg$Site_ID
veg_species_wide_beta_reveg <- veg_species_wide_beta_reveg %>% select(-Site_ID)
dim(veg_species_wide_beta_reveg) # 57 348

# View(veg_species_wide_beta)
# save Veg wide site data_frame

dist_mat_bray_all_rev <- vegdist(veg_species_wide_beta_reveg, method = "bray")

metadata_s$grazing_level <- ifelse(is.na(metadata_s$grazing_level) == TRUE, "Not observed", metadata_s$grazing_level)

s <- metadata_s[metadata_s$RemRev == "Revegetation" & metadata_s$Survey_Year == "2025" ,]

set.seed(123)
adonis2(dist_mat_bray_all_rev ~ grazing_level, data = s,
        strata = s$strata_ecosys_wpt)
# adonis2(formula = dist_mat_bray_all_rev ~ grazing_level, data = s, strata = s$strata_ecosys_wpt)
#               Df SumOfSqs      R2      F Pr(>F)
# grazing_level  3   1.0744 0.06138 1.1552  0.136
# Residual      53  16.4296 0.93862              
# Total         56  17.5039 1.00000       

## Beta diversity (Natives only) -----------------------------------------------
Native_veg_data <- subset(all_veg_data_v, is.na(INTRODUCED) == TRUE)

native_veg_data_v_com_beta <- Native_veg_data %>%
  distinct() %>% #make sure we don't have duplicates
  select(VisitID, Site_ID, WptID, `SCIENTIFIC NAME`, cover_midpoint) %>%
  group_by(VisitID, Site_ID, WptID, `SCIENTIFIC NAME`) %>%
  summarise(total_cover = sum(cover_midpoint, na.rm = TRUE)) %>%
  as.data.frame() %>%
  select(Site_ID, WptID, `SCIENTIFIC NAME`, total_cover) %>%
  group_by(Site_ID, WptID, `SCIENTIFIC NAME`) %>%
  summarise(mean_cover = mean(total_cover, na.rm = TRUE)) %>%
  as.data.frame()

nrow(native_veg_data_v_com_beta)
head(native_veg_data_v_com_beta)

veg_native_wide_beta <- native_veg_data_v_com_beta %>%
  pivot_wider(
    id_cols = Site_ID,             
    names_from = `SCIENTIFIC NAME`,
    values_from = mean_cover,
    values_fill = 0) %>% 
  as.data.frame()

rownames(veg_native_wide_beta) <- veg_native_wide_beta$Site_ID
veg_native_wide_beta <- veg_native_wide_beta %>% select(-Site_ID)
dim(veg_native_wide_beta) # 149 495

set.seed(123)
NMDS_bray_native <- metaMDS(veg_native_wide_beta, distance = "bray")

NMDS_bray_native$stress # 0.2288311
plot(NMDS_bray_native)

points_bray_native <- NMDS_bray_native$points

points_bray_native <- as.data.frame(points_bray_native)
points_bray_native$Site_ID <- rownames(points_bray_native)

# Attach metadata
native_veg_site_beta <- Native_veg_data %>% 
  select(Site_ID, Survey_Year, WptID, iPlantYear, Treatment, tEcoName) %>%
  distinct()
head(native_veg_site_beta)
head(points_bray_native)

metadata_s_native <- left_join(native_veg_site_beta, points_bray_native, by = "Site_ID")
metadata_s_native$Survey_Year <- as.factor(metadata_s_native$Survey_Year)

colnames(metadata_s_native)
metadata_s_native$RemRev <- ifelse(metadata_s_native$iPlantYear == "Remnant", "Remnant", "Revegetation")

ggplot(metadata_s_native, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour = RemRev))+
  geom_path(aes(group=WptID), alpha=0.75, colour = "grey")+
  scale_colour_manual(values =  c("Remnant" = "#006d2c", "Revegetation" = "#cc4c02"))+
  scale_fill_manual(values =  c("Remnant" = "#006d2c", "Revegetation" = "#cc4c02"))+
  geom_point(size= 2.5, colour = "black")+
  geom_point(size= 3.5, colour = "black")+
  geom_point(size= 3)+  
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-Bray-SurveyYear_RemRev-NATIVEONLY-Connection.pdf", path = outdir, width = 7, height = 6)

ggplot(metadata_s_native, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =iPlantYear))+
  geom_path(aes(group=WptID), alpha=0.75, colour = "grey")+
  geom_point(size= 2.5, colour = "black")+
  geom_point(size= 3.5, colour = "black")+
  geom_point(size= 3)+  
  scale_colour_brewer(palette = "Set1")+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-Bray-SurveyYear_PlantYear-NATIVEONLY-Connected.pdf", path = outdir, width = 7, height = 6)

metadata_s_native <- left_join(metadata_s_native, ilumpDF, by = "iLumpEco")
ggplot(metadata_s_native, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =tEcoName))+
  scale_colour_brewer(palette = "Set1")+
  geom_path(aes(group=WptID), alpha=0.75, colour = "grey")+
  geom_point(size= 2.5, colour = "black")+
  geom_point(size= 3.5, colour = "black")+
  geom_point(size= 3)+  
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-Bray-SurveyYear_Ecosystem-NATIVEONLY-Connected.pdf", path = outdir, width = 7, height = 6)

# Jaccard
veg_native_PA <- ifelse(veg_native_wide_beta == 0, 0, 1) 

set.seed(123)
NMDS_Jacc_native <- metaMDS(veg_native_PA, distance = "jaccard", binary = TRUE)
NMDS_Jacc_native$stress # 0.2280751

plot(NMDS_Jacc_native)
plot(NMDS_bray_native)

points_Jacc_native <- NMDS_Jacc_native$points

points_Jacc_native <- as.data.frame(points_Jacc_native)
points_Jacc_native$Site_ID <- rownames(points_Jacc_native)

# Attach metadata
native_veg_site_beta <- Native_veg_data %>% 
  select(Site_ID, Survey_Year, WptID, iPlantYear, Treatment, iLumpEco) %>%
  distinct()
head(native_veg_site_beta)
head(points_Jacc_native)

metadata_s_native_Jacc <- left_join(native_veg_site_beta, points_Jacc_native, by = "Site_ID")
metadata_s_native_Jacc$Survey_Year <- as.factor(metadata_s_native_Jacc$Survey_Year)

colnames(metadata_s_native_Jacc)
# "Site_ID"      "Survey_Year"  "iWptID"       "iPlantYear"   "iTreatID"     "iEcosystemID"
# "tEcoName"     "MDS1"         "MDS2"

metadata_s_native_Jacc$RemRev <- ifelse(metadata_s_native_Jacc$iPlantYear == "Remnant", "Remnant", "Revegetation")

ggplot(metadata_s_native_Jacc, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =RemRev))+
  geom_path(aes(group=WptID), alpha=0.75, colour = "grey")+
  scale_colour_manual(values =  c("Remnant" = "#006d2c", "Revegetation" = "#cc4c02"))+
  scale_fill_manual(values =  c("Remnant" = "#006d2c", "Revegetation" = "#cc4c02"))+
  geom_point(size= 2.5, colour = "black")+
  geom_point(size= 3.5, colour = "black")+
  geom_point(size= 3)+  
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-Jaccard-SurveyYear_RemRev-NATIVEONLY.pdf", path = outdir, width = 7, height = 6)

ggplot(metadata_s_native_Jacc, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =iPlantYear))+
  geom_point(size= 2)+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-Jaccard-SurveyYear_PlantYear-NATIVEONLY.pdf", path = outdir, width = 7, height = 6)

ggplot(metadata_s_native_Jacc, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =as.factor(iLumpEco)))+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  geom_point(size= 2)+
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-Jaccard-SurveyYear_Ecosystem-NATIVEONLY.pdf", path = outdir, width = 7, height = 6)


## Beta diversity (Invasives only) -----------------------------------------------
Invasive_veg_data <- subset(all_veg_data_v, INTRODUCED == "*")

invasive_veg_data_v_com_beta <- Invasive_veg_data %>%
  distinct() %>% #make sure we don't have duplicates
  select(VisitID, Site_ID, WptID, `SCIENTIFIC NAME`, cover_midpoint) %>%
  group_by(VisitID, Site_ID, WptID, `SCIENTIFIC NAME`) %>%
  summarise(total_cover = sum(cover_midpoint, na.rm = TRUE)) %>%
  as.data.frame() %>%
  select(Site_ID, WptID, `SCIENTIFIC NAME`, total_cover) %>%
  group_by(Site_ID, WptID, `SCIENTIFIC NAME`) %>%
  summarise(mean_cover = mean(total_cover, na.rm = TRUE)) %>%
  as.data.frame()

nrow(invasive_veg_data_v_com_beta)
head(invasive_veg_data_v_com_beta)

veg_invasive_wide_beta <- invasive_veg_data_v_com_beta %>%
  pivot_wider(
    id_cols = Site_ID,             
    names_from = `SCIENTIFIC NAME`,
    values_from = mean_cover,
    values_fill = 0) %>% 
  as.data.frame()

rownames(veg_invasive_wide_beta) <- veg_invasive_wide_beta$Site_ID
veg_invasive_wide_beta <- veg_invasive_wide_beta %>% select(-Site_ID)
dim(veg_invasive_wide_beta) # 149 191

set.seed(123)
NMDS_bray_invasive <- metaMDS(veg_invasive_wide_beta, distance = "bray")

NMDS_bray_invasive$stress # 0.2400563
plot(NMDS_bray_invasive)

points_bray_invasive <- NMDS_bray_invasive$points

points_bray_invasive <- as.data.frame(points_bray_invasive)
points_bray_invasive$Site_ID <- rownames(points_bray_invasive)

# Attach metadata
invasive_veg_site_beta <- Invasive_veg_data %>% 
  select(Site_ID, Survey_Year, WptID, iPlantYear, Treatment, tEcoName) %>%
  distinct()
head(invasive_veg_site_beta)
head(points_bray_invasive)

metadata_s_invasive <- left_join(invasive_veg_site_beta, points_bray_invasive, by = "Site_ID")
metadata_s_invasive$Survey_Year <- as.factor(metadata_s_invasive$Survey_Year)

colnames(metadata_s_invasive)
# "Site_ID"      "Survey_Year"  "iWptID"       "iPlantYear"   "iTreatID"     "iEcosystemID"
# "tEcoName"     "MDS1"         "MDS2"

metadata_s_invasive$RemRev <- ifelse(metadata_s_invasive$iPlantYear == "Remnant", "Remnant", "Revegetation")

ggplot(metadata_s_invasive, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour = RemRev))+
  geom_path(aes(group=WptID), alpha=0.75, colour = "grey")+
  scale_colour_manual(values =  c("Remnant" = "#006d2c", "Revegetation" = "#cc4c02"))+
  scale_fill_manual(values =  c("Remnant" = "#006d2c", "Revegetation" = "#cc4c02"))+
  geom_point(size= 2.5, colour = "black")+
  geom_point(size= 3.5, colour = "black")+
  geom_point(size= 3)+  
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-Bray-SurveyYear_RemRev-INVASIVEONLY.pdf", path = outdir, width = 7, height = 6)

ggplot(metadata_s_invasive, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =iPlantYear))+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  geom_point(size= 2)+
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-Bray-SurveyYear_PlantYear-INVASIVEONLY.pdf", path = outdir, width = 7, height = 6)

ggplot(metadata_s_invasive, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =tEcoName))+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  geom_point(size= 2)+
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-Bray-SurveyYear_Ecosystem-INVASIVEONLY.pdf", path = outdir, width = 7, height = 6)

# Jaccard
veg_invasive_PA <- ifelse(veg_invasive_wide_beta == 0, 0, 1) 

set.seed(123)
NMDS_Jacc_invasive <- metaMDS(veg_invasive_PA, distance = "jaccard", binary = TRUE)
NMDS_Jacc_invasive$stress # 0.2388271

plot(NMDS_Jacc_invasive)
plot(NMDS_bray_invasive)

points_Jacc_invasive <- NMDS_Jacc_invasive$points

points_Jacc_invasive <- as.data.frame(points_Jacc_invasive)
points_Jacc_invasive$Site_ID <- rownames(points_Jacc_invasive)

# Attach metadata
invasive_veg_site_beta <- Invasive_veg_data %>% 
  select(Site_ID, Survey_Year, WptID, iPlantYear, Treatment, tEcoName) %>%
  distinct()
head(invasive_veg_site_beta)
head(points_Jacc_invasive)

metadata_s_invasive_Jacc <- left_join(invasive_veg_site_beta, points_Jacc_invasive, by = "Site_ID")
metadata_s_invasive_Jacc$Survey_Year <- as.factor(metadata_s_invasive_Jacc$Survey_Year)

colnames(metadata_s_invasive_Jacc)
# "Site_ID"      "Survey_Year"  "iWptID"       "iPlantYear"   "iTreatID"     "iEcosystemID"
# "tEcoName"     "MDS1"         "MDS2"

metadata_s_invasive_Jacc$RemRev <- ifelse(metadata_s_invasive_Jacc$iPlantYear == "Remnant", "Remnant", "Revegetation")

ggplot(metadata_s_invasive_Jacc, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =RemRev))+
  geom_path(aes(group=WptID), alpha=0.75, colour = "grey")+
  scale_colour_manual(values =  c("Remnant" = "#006d2c", "Revegetation" = "#cc4c02"))+
  scale_fill_manual(values =  c("Remnant" = "#006d2c", "Revegetation" = "#cc4c02"))+
  geom_point(size= 2.5, colour = "black")+
  geom_point(size= 3.5, colour = "black")+
  geom_point(size= 3)+  
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-Jaccard-SurveyYear_RemRev-INVASIVEONLY.pdf", path = outdir, width = 7, height = 6)

ggplot(metadata_s_invasive_Jacc, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =iPlantYear))+
  geom_point(size= 2)+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-Jaccard-SurveyYear_PlantYear-INVASIVEONLY.pdf", path = outdir, width = 7, height = 6)

ggplot(metadata_s_invasive_Jacc, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =tEcoName))+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  geom_point(size= 2)+
  theme_test()
# ggsave(filename = "Veg-BETA-NMDS-Jaccard-SurveyYear_Ecosystem-INVASIVEONLY.pdf", path = outdir, width = 7, height = 6)


# Percentage exotic plants -----------------------------------------------------
colnames(all_veg_data_v)

# To get the proportion of exotic plants at the Site level (mean of each visit)
prop_exotic_df <- all_veg_data_v %>%
  distinct() %>%
  select(Survey_Year, WptID, VisitID, iPlantYear, tEcoName, `SCIENTIFIC NAME`, INTRODUCED, cover_midpoint, HILG) %>%
  group_by(Survey_Year, WptID, `SCIENTIFIC NAME`, INTRODUCED) %>%
  summarise(total_cover = mean(cover_midpoint, na.rm = TRUE), .groups = "drop") %>%
  group_by(Survey_Year, WptID, `SCIENTIFIC NAME`, INTRODUCED) %>%
  summarise(average_cover = mean(total_cover))

prop_exotic_df$INTRODUCED <- ifelse(is.na(prop_exotic_df$INTRODUCED) ==TRUE, 0, 1)
head(prop_exotic_df)

# Get averages across visits
prop_exotic_df2 <- prop_exotic_df %>% 
  group_by(Survey_Year, WptID) %>% 
  mutate(average_cover_p = average_cover / sum(average_cover))

head(prop_exotic_df2)
# # sanity check
# prop_exotic_df2%>%
#   group_by(Survey_Year, iWptID) %>%
#   reframe(sum(average_cover_p))

# Get percentage exotic output
prop_exotic_df3 <- prop_exotic_df2 %>% 
  group_by(Survey_Year, WptID) %>% 
  summarise(
    prop_exotic = sum(average_cover_p[INTRODUCED == 1], na.rm = TRUE),
    .groups = "drop"
  )
head(prop_exotic_df3)
# # A tibble: 6 × 3
#   Survey_Year iWptID prop_exotic
#  <dbl>  <dbl>       <dbl>
# 1        2015    466       0.372
# 2        2015    470       0.124
# 3        2015    476       0.335
# 4        2015    479       0.203
# 5        2015    480       0.192
# 6        2015    482       0.260

# To get the proportion of exotic plants within each Site and Visit 
# Get percentage exotic output

prop_exotic_df_visits <- all_veg_data_v %>%
  distinct() %>%
  select(Survey_Year, WptID, VisitID, iPlantYear, tEcoName, `SCIENTIFIC NAME`, INTRODUCED, cover_midpoint) %>%
  group_by(Survey_Year, WptID, VisitID, `SCIENTIFIC NAME`, INTRODUCED) %>%
  summarise(total_cover = mean(cover_midpoint, na.rm = TRUE), .groups = "drop")

prop_exotic_df_visits$INTRODUCED <- ifelse(is.na(prop_exotic_df_visits$INTRODUCED) ==TRUE, 0, 1)
head(prop_exotic_df_visits)

prop_exotic_df_visits2 <- prop_exotic_df_visits %>% 
  group_by(Survey_Year, WptID, VisitID) %>% 
  mutate(average_cover_p = total_cover / sum(total_cover))
head(prop_exotic_df_visits2)

# Get averages across visits
prop_exotic_df_visits3 <- prop_exotic_df_visits2 %>% 
  group_by(Survey_Year, WptID, VisitID) %>% 
  summarise(
    prop_exotic = sum(average_cover_p[INTRODUCED == 1], na.rm = TRUE),
    .groups = "drop"
  )
head(prop_exotic_df_visits3)

# Put site and visit level data back into a more full dataframe with associated site info
head(all_veg_data_v)
all_veg_data_visit_infoonly <- all_veg_data_v %>%
  select(Survey_Year, WptID, VisitID, iPlantYear, tEcoName, HILG) %>%
  distinct()
head(all_veg_data_visit_infoonly)

all_veg_data_site_infoonly <- all_veg_data_v %>%
  select(Survey_Year, WptID, HILG, iPlantYear, tEcoName) %>%
  distinct()
head(all_veg_data_site_infoonly)


prop_exotic_df_visits4 <- left_join(prop_exotic_df_visits3, all_veg_data_visit_infoonly,
                                    by=c("Survey_Year", "WptID", "VisitID"))
prop_exotic_df4 <- left_join(prop_exotic_df3, all_veg_data_site_infoonly,
                             by=c("Survey_Year", "WptID"))

# Plot and analyse
prop_exotic_df_visits4$RevRem <- ifelse(prop_exotic_df_visits4$iPlantYear == "Remnant", "Remnant", "Revegetation")
prop_exotic_df_visits4$RevRem <- factor(prop_exotic_df_visits4$RevRem, levels = c("Revegetation", "Remnant"))

head(prop_exotic_df4)
prop_exotic_df4$Survey_Year <- as.factor(prop_exotic_df4$Survey_Year)
prop_exotic_df4$iPlantYear <- factor(prop_exotic_df4$iPlantYear, levels = c("2015","2014","2013","2012","Remnant"))
# prop_exotic_df4$tEcoName <- factor(prop_exotic_df4$tEcoName, levels = c("A.verticillata","C.gracilis","E.diversifolia","E.fasciculosa","E.porosa","E.incrassata","E.leucoxylon","E.odorata"))

prop_exotic_df4$RevRem <- ifelse(prop_exotic_df4$iPlantYear == "Remnant", "Remnant", "Revegetation")
prop_exotic_df4$RevRem <- factor(prop_exotic_df4$RevRem, levels = c("Revegetation", "Remnant"))

colnames(prop_exotic_df4) <- c("Survey_Year", "WptID", "prop_exotic", "HILG", "PlantYear", "Ecosystem", "RevRem" )
ggplot(prop_exotic_df4, aes(x = PlantYear, y = prop_exotic, fill=PlantYear))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(height = 0, width = 0.1)+
  facet_grid(Survey_Year~Ecosystem, scales = "free_x", space = "free_x")+
  theme_test()
# ggsave(filename = "Veg-Proportion_Exotic-SurveyYear_PlantYear.pdf", path = outdir, width = 14, height = 7)

ggplot(prop_exotic_df4, aes(x = RevRem, y = prop_exotic, fill=RevRem))+
  geom_violin()+
  geom_boxplot(outlier.shape = NA, width = 0.1, fill= "white")+
  geom_jitter(height = 0, width = 0.1, aes(colour = HILG))+
  facet_grid(Survey_Year~Ecosystem, scales = "free_x", space = "free_x")+
  theme_test()+
  scale_color_manual(values = ("HILG" = "red"), na.value = "black")
# ggsave(filename = "Veg-Proportion_Exotic-SurveyYear_RemRev.pdf", path = outdir, width = 14, height = 7)

ggplot(prop_exotic_df4, aes(x = RevRem, y = prop_exotic, fill=RevRem))+
  geom_violin()+
  geom_boxplot(outlier.shape = NA, width = 0.1, fill= "white")+
  geom_jitter(height = 0, width = 0.1, aes(colour = HILG))+
  facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
  scale_fill_manual(values = c("Revegetation"="#cc4c02","Remnant"="#41ab5d"))+
  theme_test()+
  scale_color_manual(values = ("HILG" = "red"), na.value = "black")
# ggsave(filename = "Veg-Proportion_Exotic-SurveyYear_RemRev-condensed.pdf", path = outdir, width = 7, height = 6)

# statsy things
prop_lmem <- LMEM_permute_anova(data = prop_exotic_df4,
                                response = "prop_exotic", fixed_effects = c("RevRem", "Survey_Year", "RevRem:Survey_Year"),
                                random_effects = c("WptID", "Ecosystem"),
                                nreps = 999)
prop_lmem$formula
# prop_exotic ~ RevRem + Survey_Year + RevRem:Survey_Year + (1 | WptID) + (1 | Ecosystem)
prop_lmem$Anova
#                      Chisq Df Pr(>Chisq)    
# RevRem              83.078  1  < 2.2e-16 ***
# Survey_Year        200.797  1  < 2.2e-16 ***
# RevRem:Survey_Year  60.566  1  7.113e-15 ***
prop_lmem$permutation_p_values
# RevRem        Survey_Year RevRem:Survey_Year 
#      0                  0                  0 

library(emmeans)
library(multcomp)
library(multcompView)

emm <- emmeans(model_exotic_visits, ~ RevRem * Survey_Year)
emm
pairs(emm, adjust = "tukey")
# Contrast                                                    estimate     SE     df t.ratio p.value
# Revegetation Survey_Year2015 - Remnant Survey_Year2015        0.5817 0.0348   94.8  16.699  <.0001
# Revegetation Survey_Year2015 - Revegetation Survey_Year2025   0.2918 0.0113 1277.3  25.899  <.0001
# Revegetation Survey_Year2015 - Remnant Survey_Year2025        0.5694 0.0349   95.0  16.335  <.0001
# Remnant Survey_Year2015 - Revegetation Survey_Year2025       -0.2899 0.0349   95.5  -8.304  <.0001
# Remnant Survey_Year2015 - Remnant Survey_Year2025            -0.0123 0.0213 1264.1  -0.579  0.9383
# Revegetation Survey_Year2025 - Remnant Survey_Year2025        0.2776 0.0349   95.7   7.946  <.0001
# 
# Degrees-of-freedom method: kenward-roger 
# P value adjustment: tukey method for comparing a family of 4 estimates 

cld(emm, adjust = "tukey")

# Indicator species analysis ----------------------------------------------
# https://www.stat.math.ethz.ch/CRAN/web/packages/indicspecies/vignettes/IndicatorSpeciesAnalysis.html  
# https://doi.org/10.1038/s41598-025-16501-8
# Remnant vs Reveg species

library(indicspecies)

# Cover/abundance table
# View(veg_species_wide_beta)

## Split data by Survey years  -------------------------------------------------
veg_species_wide_beta_2015 <- veg_species_wide_beta[metadata_s$Survey_Year == "2015",]
dim(veg_species_wide_beta_2015) # 76 686
metadata_s_2015 <- subset(metadata_s, Survey_Year == "2015")

veg_species_wide_beta_2025 <- veg_species_wide_beta[metadata_s$Survey_Year == "2025",]
dim(veg_species_wide_beta_2025) # 73 686
metadata_s_2025 <- subset(metadata_s, Survey_Year == "2025")

# Compute 2015
multipatt.obj_2015 <- multipatt(x = veg_species_wide_beta_2015, 
                                cluster = metadata_s_2015$RemRev,
                                control = how(nperm=999))
summary(multipatt.obj_2015)

# strassoc comparisoins (strength of associations)
strassoc.obj_2015 <- strassoc(veg_species_wide_beta_2015, 
                              cluster = metadata_s_2015$RemRev)

strassoc.obj_2015 <- as.data.frame(strassoc.obj_2015)
strassoc.obj_2015$Species <- rownames(strassoc.obj_2015)
summary(multipatt.obj_2015)

multipatt.obj_2015_sign <- multipatt.obj_2015$sign
multipatt.obj_2015_sign$Species <- rownames(multipatt.obj_2015_sign)

# Compute 2025
multipatt.obj_2025 <- multipatt(x = veg_species_wide_beta_2025, 
                                cluster = metadata_s_2025$RemRev,
                                control = how(nperm=999))
summary(multipatt.obj_2025)
summary(multipatt.obj_2025, indvalcomp=TRUE)

# strassoc comparisoins (strength of associations)
strassoc.obj_2025 <- strassoc(veg_species_wide_beta_2025, 
                              cluster = metadata_s_2025$RemRev)

strassoc.obj_2025 <- as.data.frame(strassoc.obj_2025)
strassoc.obj_2025$Species <- rownames(strassoc.obj_2025)
summary(multipatt.obj_2025)

multipatt.obj_2025_sign <- multipatt.obj_2025$sign
multipatt.obj_2025_sign$Species <- rownames(multipatt.obj_2025_sign)

## Split data by Reveg/rem treatments ------------------------------------------
veg_species_wide_beta_Reveg <- veg_species_wide_beta[metadata_s$RemRev == "Revegetation",]
dim(veg_species_wide_beta_Reveg) # 76 686
metadata_s_Reveg <- subset(metadata_s, RemRev == "Revegetation")

veg_species_wide_beta_Remnant <- veg_species_wide_beta[metadata_s$RemRev == "Remnant",]
dim(veg_species_wide_beta_Remnant) # 73 686
metadata_s_Remnant <- subset(metadata_s, RemRev == "Remnant")

# Compute Reveg
multipatt.obj_Reveg <- multipatt(x = veg_species_wide_beta_Reveg, 
                                 cluster = metadata_s_Reveg$Survey_Year,
                                 control = how(nperm=999))
summary(multipatt.obj_Reveg)

# strassoc comparisoins (strength of associations)
strassoc.obj_Reveg <- strassoc(veg_species_wide_beta_Reveg, 
                               cluster = metadata_s_Reveg$Survey_Year)

strassoc.obj_Reveg <- as.data.frame(strassoc.obj_Reveg)
strassoc.obj_Reveg$Species <- rownames(strassoc.obj_Reveg)
summary(multipatt.obj_Reveg)

multipatt.obj_Reveg_sign <- multipatt.obj_Reveg$sign
multipatt.obj_Reveg_sign$Species <- rownames(multipatt.obj_Reveg_sign)

# Compute Remnant
multipatt.obj_Remnant <- multipatt(x = veg_species_wide_beta_Remnant, 
                                   cluster = metadata_s_Remnant$Survey_Year,
                                   control = how(nperm=999))
summary(multipatt.obj_Remnant)
summary(multipatt.obj_Remnant, indvalcomp=TRUE)

# strassoc comparisoins (strength of associations)
strassoc.obj_Remnant <- strassoc(veg_species_wide_beta_Remnant, 
                                 cluster = metadata_s_Remnant$Survey_Year)

strassoc.obj_Remnant <- as.data.frame(strassoc.obj_Remnant)
strassoc.obj_Remnant$Species <- rownames(strassoc.obj_Remnant)
summary(multipatt.obj_Remnant)

multipatt.obj_Remnant_sign <- multipatt.obj_Remnant$sign
multipatt.obj_Remnant_sign$Species <- rownames(multipatt.obj_Remnant_sign)

### Diverging bar plot ---------------------------------------------------------
# 2015
combined_ind_2015 <- left_join(multipatt.obj_2015_sign, strassoc.obj_2015, by = "Species")

# get introduced species
veg_data_intr <- veg_data_long %>%
  select(`SCIENTIFIC NAME`, INTRODUCED) %>%
  distinct()
colnames(veg_data_intr) <- c("Species", "Exotic")

combined_ind_2015 <- left_join(combined_ind_2015, veg_data_intr, by = "Species")

class(multipatt.obj_2015_sign)
class(strassoc.obj_2015)

combined_ind_2015 <- combined_ind_2015 %>% 
  filter(p.value < 0.05) %>% 
  mutate(Species = factor(Species, levels = Species[order(Remnant, decreasing = FALSE)]),
         Native_Status = ifelse(Exotic == 1, "Non-native", "Native"))

combined_ind_2015$Exotic <- as.factor(combined_ind_2015$Exotic)

# colour labels
species_colors_2015 <- setNames(
  ifelse(combined_ind_2015$Native_Status == "Non-native", "orange", "darkgreen"),
  combined_ind_2015$Species
)

# Indicator species analysis plot
ggplot(combined_ind_2015, aes(x = Species, y = Remnant, fill = Native_Status)) +
  geom_col(colour= "black") +
  coord_flip() +
  scale_fill_manual(values = c("Non-native" = "orange", "Native" = "darkgreen"))+
  theme_test() +
  theme(axis.text.y = element_text(color = species_colors_2015[levels(factor(combined_ind_2015$Species))]))+
  labs(y = "Species-group association (Reference: Remnant)", x = "Species", title = "2015 Survey")
# ggsave(filename = "Veg-IndicatorSpecies-RemRev-2015.pdf", path = outdir, width = 9, height = 14)

# 2025
combined_ind_2025 <- left_join(multipatt.obj_2025_sign, strassoc.obj_2025, by = "Species")

# get introduced species
veg_data_intr <- veg_data_long %>%
  select(`SCIENTIFIC NAME`, INTRODUCED) %>%
  distinct()
colnames(veg_data_intr) <- c("Species", "Exotic")

combined_ind_2025 <- left_join(combined_ind_2025, veg_data_intr, by = "Species")

class(multipatt.obj_2025_sign)
class(strassoc.obj_2025)

combined_ind_2025 <- combined_ind_2025 %>% 
  filter(p.value < 0.05) %>% 
  mutate(Species = factor(Species, levels = Species[order(Remnant, decreasing = FALSE)]),
         Native_Status = ifelse(Exotic == 1, "Non-native", "Native"))

combined_ind_2025$Exotic <- as.factor(combined_ind_2025$Exotic)

# colour labels
species_colors_2025 <- setNames(
  ifelse(combined_ind_2025$Native_Status == "Non-native", "orange", "darkgreen"),
  combined_ind_2025$Species
)

# Indicator species analysis plot
ggplot(combined_ind_2025, aes(x = Species, y = Remnant, fill = Native_Status)) +
  geom_col(colour= "black") +
  coord_flip() +
  scale_fill_manual(values = c("Non-native" = "orange", "Native" = "darkgreen"))+
  theme_test() +
  theme(axis.text.y = element_text(color = species_colors_2025[levels(factor(combined_ind_2025$Species))]))+
  labs(y = "Species-group association", x = "Species", title = "2025 Survey")
# ggsave(filename = "Veg-IndicatorSpecies-RemRev-2025.pdf", path = outdir, width = 9, height = 14)

# Plot both together:
ggpubr::ggarrange(
  ggplot(combined_ind_2015, aes(x = Species, y = Remnant, fill = Native_Status)) +
    geom_col(colour= "black") +
    coord_flip() +
    scale_fill_manual(values = c("Non-native" = "orange", "Native" = "darkgreen"))+
    theme_test() +
    theme(axis.text.y = element_text(color = species_colors_2015[levels(factor(combined_ind_2015$Species))]))+
    labs(y = "Species-group association", x = "Species", title = "2015 Survey"),
  ggplot(combined_ind_2025, aes(x = Species, y = Remnant, fill = Native_Status)) +
    geom_col(colour= "black") +
    coord_flip() +
    scale_fill_manual(values = c("Non-native" = "orange", "Native" = "darkgreen"))+
    theme_test() +
    theme(axis.text.y = element_text(color = species_colors_2025[levels(factor(combined_ind_2025$Species))]))+
    labs(y = "Species-group association", x = "Species", title = "2025 Survey"),
  align = "hv", common.legend = TRUE
)
# ggsave(filename = "Veg-IndicatorSpecies-RemRev-2015-2025.pdf", path = outdir, width = 18, height = 14)

# Top 15 positive and negative difference values
# Revegetation: top 15 + bottom 15
combined_ind_2015_topbot <- combined_ind_2015 %>%
  arrange(Remnant) %>%
  slice(c(1:15, (n() - 14):n()))

# Remnant: top 15 + bottom 15
combined_ind_2025_topbot <- combined_ind_2025 %>%
  arrange(Remnant) %>%
  slice(c(1:15, (n() - 14):n()))

ggpubr::ggarrange(
  ggplot(combined_ind_2015, aes(x = Species, y = Remnant, fill = Native_Status)) +
    geom_col(colour= "black") +
    coord_flip()+
    scale_fill_manual(values = c("Non-native" = "orange", "Native" = "darkgreen"))+
    theme_test() +
    theme(axis.text.y = element_text(color = species_colors_2015[levels(factor(combined_ind_2015$Species))]))+
    labs(y = "Species-group association\n (Reveg indicators << | >> Remnant indicators)", x = "", title = "2025 surveys\n15 most strongly ± associated species", fill = "Native status"),
  ggplot(combined_ind_2025_topbot, aes(x = Species, y = Remnant, fill = Native_Status)) +
    geom_col(colour= "black") +
    coord_flip()+
    scale_fill_manual(values = c("Non-native" = "orange", "Native" = "darkgreen"))+
    theme_test() +
    theme(axis.text.y = element_text(color = species_colors_2025[levels(factor(combined_ind_2025_topbot$Species))]))+
    labs(y = "Species-group association\n (Reveg indicators << | >> Remnant indicators)", x = "", title = "2025 surveys\n15 most strongly ± associated species", fill = "Native status")
  ,
  align ="hv", common.legend = TRUE
)
# ggsave(filename = "Veg-IndicatorSpecies-SurveyYears-Remnant_Revegetation-TOP15.pdf", path = outdir, width = 12, height = 6)


## Indicators over time
veg_species_wide_beta_Rev <- veg_species_wide_beta[metadata_s$RemRev == "Revegetation",]
dim(veg_species_wide_beta_Rev) # 117 686
metadata_s_Rev <- subset(metadata_s, RemRev == "Revegetation")

veg_species_wide_beta_Rem <- veg_species_wide_beta[metadata_s$RemRev == "Remnant",]
dim(veg_species_wide_beta_Rem) # 32 686
metadata_s_Rem <- subset(metadata_s, RemRev == "Remnant")

# Compute Revegetation comparisons
multipatt.obj_Rev <- multipatt(x = veg_species_wide_beta_Rev, 
                               cluster = metadata_s_Rev$Survey_Year,
                               control = how(nperm=999))
summary(multipatt.obj_Rev)

# strassoc comparisoins (strength of associations)
strassoc.obj_Rev <- strassoc(veg_species_wide_beta_Rev, 
                             cluster = metadata_s_Rev$Survey_Year)

strassoc.obj_Rev <- as.data.frame(strassoc.obj_Rev)
strassoc.obj_Rev$Species <- rownames(strassoc.obj_Rev)

multipatt.obj_Rev_sign <- multipatt.obj_Rev$sign
multipatt.obj_Rev_sign$Species <- rownames(multipatt.obj_Rev_sign)

# Revegetation
combined_ind_Rev <- left_join(multipatt.obj_Rev_sign, strassoc.obj_Rev, by = "Species")

# get introduced species
veg_data_intr <- veg_data_long %>%
  select(`SCIENTIFIC NAME`, INTRODUCED) %>%
  distinct()
colnames(veg_data_intr) <- c("Species", "Exotic")

combined_ind_Rev2 <- left_join(combined_ind_Rev, veg_data_intr, by = "Species")

class(multipatt.obj_Rev_sign)
class(strassoc.obj_Rev)

combined_ind_Rev3 <- combined_ind_Rev2 %>% 
  filter(p.value < 0.05) %>% 
  mutate(Species = factor(Species, levels = Species[order(`2025`, decreasing = FALSE)]),
         Native_Status = ifelse(Exotic == 1, "Non-native", "Native"))

combined_ind_Rev3$Exotic <- as.factor(combined_ind_Rev3$Exotic)

# colour labels
species_colours_Rev3 <- setNames(
  ifelse(combined_ind_Rev3$Native_Status == "Non-native", "orange", "darkgreen"),
  combined_ind_Rev3$Species
)

# Indicator species analysis plot
ggplot(combined_ind_Rev3, aes(x = Species, y = `2025`, fill = Native_Status)) +
  geom_col(colour= "black") +
  coord_flip()+
  scale_fill_manual(values = c("Non-native" = "orange", "Native" = "darkgreen"))+
  theme_test() +
  # theme(axis.text.y = element_text(color = species_colors[levels(factor(combined_ind_2025$Species))]))+
  theme(axis.text.y = element_text(color = species_colours_Rev3[levels(factor(combined_ind_Rev3$Species))]))+
  labs(y = "Species-group association", x = "Species", title = "Revegetation sites")
# ggsave(filename = "Veg-IndicatorSpecies-SurveyYears-revegetation_only-short.pdf", path = outdir, width = 9, height = 7)

# Remnant
# Compute Remnant comparisons
multipatt.obj_Rem <- multipatt(x = veg_species_wide_beta_Rem, 
                               cluster = metadata_s_Rem$Survey_Year,
                               control = how(nperm=999))
summary(multipatt.obj_Rem)

# strassoc comparisoins (strength of associations)
strassoc.obj_Rem <- strassoc(veg_species_wide_beta_Rem, 
                             cluster = metadata_s_Rem$Survey_Year)

strassoc.obj_Rem <- as.data.frame(strassoc.obj_Rem)
strassoc.obj_Rem$Species <- rownames(strassoc.obj_Rem)

multipatt.obj_Rem_sign <- multipatt.obj_Rem$sign
multipatt.obj_Rem_sign$Species <- rownames(multipatt.obj_Rem_sign)

# Remnant
combined_ind_Rem <- left_join(multipatt.obj_Rem_sign, strassoc.obj_Rem, by = "Species")

# get introduced species
veg_data_intr <- veg_data_long %>%
  select(`SCIENTIFIC NAME`, INTRODUCED) %>%
  distinct()
colnames(veg_data_intr) <- c("Species", "Exotic")

combined_ind_Rem2 <- left_join(combined_ind_Rem, veg_data_intr, by = "Species")

class(multipatt.obj_Rem_sign)
class(strassoc.obj_Rem)

combined_ind_Rem3 <- combined_ind_Rem2 %>% 
  filter(p.value < 0.05) %>% 
  mutate(Species = factor(Species, levels = Species[order(`2025`, decreasing = FALSE)]),
         Native_Status = ifelse(Exotic == 1, "Non-native", "Native"))

combined_ind_Rem3$Exotic <- as.factor(combined_ind_Rem3$Exotic)

# colour labels
species_colours_Rem3 <- setNames(
  ifelse(combined_ind_Rem3$Native_Status == "Non-native", "orange", "darkgreen"),
  combined_ind_Rem3$Species
)

# Indicator species analysis plot
ggplot(combined_ind_Rem3, aes(x = Species, y = `2025`, fill = Native_Status)) +
  geom_col(colour= "black") +
  coord_flip()+
  scale_fill_manual(values = c("Non-native" = "orange", "Native" = "darkgreen"))+
  theme_test() +
  # theme(axis.text.y = element_text(color = species_colors[levels(factor(combined_ind_2025$Species))]))+
  theme(axis.text.y = element_text(color = species_colours_Rem3[levels(factor(combined_ind_Rem3$Species))]))+
  labs(y = "Species-group association (Reference: 2025 surveys)", x = "Species", title = "Remnant sites")
# ggsave(filename = "Veg-IndicatorSpecies-SurveyYears-remnant_only-short.pdf", path = outdir, width = 9, height = 7)

# Plot together
ggpubr::ggarrange(
  ggplot(combined_ind_Rev3, aes(x = Species, y = `2025`, fill = Native_Status)) +
    geom_col(colour= "black") +
    coord_flip()+
    scale_fill_manual(values = c("Non-native" = "orange", "Native" = "darkgreen"))+
    theme_test() +
    # theme(axis.text.y = element_text(color = species_colors[levels(factor(combined_ind_2025$Species))]))+
    theme(axis.text.y = element_text(color = species_colours_Rev3[levels(factor(combined_ind_Rev3$Species))]))+
    labs(y = "Species-group association\n (2015 indicators <- | -> 2025 indicators)", x = "Species", title = "Revegetation sites"),
  ggplot(combined_ind_Rem3, aes(x = Species, y = `2025`, fill = Native_Status)) +
    geom_col(colour= "black") +
    coord_flip()+
    scale_fill_manual(values = c("Non-native" = "orange", "Native" = "darkgreen"))+
    theme_test() +
    # theme(axis.text.y = element_text(color = species_colors[levels(factor(combined_ind_2025$Species))]))+
    theme(axis.text.y = element_text(color = species_colours_Rem3[levels(factor(combined_ind_Rem3$Species))]))+
    labs(y = "Species-group association\n (2015 indicators <- | -> 2025 indicators)", x = "Species", title = "Remnant sites")
  ,
  align ="hv", common.legend = TRUE
)
# ggsave(filename = "Veg-IndicatorSpecies-SurveyYears-Remnant_Revegetation.pdf", path = outdir, width = 13, height = 5)

# Top 15 positive and negative difference values
# Revegetation: top 15 + bottom 15
combined_ind_Rev_topbot <- combined_ind_Rev3 %>%
  arrange(`2025`) %>%
  slice(c(1:15, (n() - 14):n()))

# Remnant: top 15 + bottom 15
combined_ind_Rem_topbot <- combined_ind_Rem3 %>%
  arrange(`2025`) %>%
  slice(c(1:15, (n() - 14):n()))

ggpubr::ggarrange(
  ggplot(combined_ind_Rev_topbot, aes(x = Species, y = `2025`, fill = Native_Status)) +
    geom_col(colour= "black") +
    coord_flip()+
    scale_fill_manual(values = c("Non-native" = "orange", "Native" = "darkgreen"))+
    theme_test() +
    theme(axis.text.y = element_text(color = species_colours_Rev3[levels(factor(combined_ind_Rev_topbot$Species))]))+
    labs(y = "Species-group association\n 2015 indicators <<     |     >> 2025 indicators", x = "", title = "Revegetation sites\n15 most strongly ± associated species", fill = "Native status"),
  ggplot(combined_ind_Rem_topbot, aes(x = Species, y = `2025`, fill = Native_Status)) +
    geom_col(colour= "black") +
    coord_flip()+
    scale_fill_manual(values = c("Non-native" = "orange", "Native" = "darkgreen"))+
    theme_test() +
    theme(axis.text.y = element_text(color = species_colours_Rem3[levels(factor(combined_ind_Rem_topbot$Species))]))+
    labs(y = "Species-group association\n 2015 indicators <<     |     >> 2025 indicators", x = "", title = "Remnant sites\n15 most strongly ± associated species", fill = "Native status")
  ,
  align ="hv", common.legend = TRUE
)
# ggsave(filename = "Veg-IndicatorSpecies-SurveyYears-Remnant_Revegetation-TOP15.pdf", path = outdir, width = 12, height = 6)


# 

ggpubr::ggarrange(
  ggplot(combined_ind_Rev3, aes(x = Species, y = `2025`, fill = Native_Status)) +
    geom_col(colour= "black") +
    coord_flip()+
    scale_fill_manual(values = c("Exotic" = "orange", "Native" = "darkgreen"))+
    theme_test() +
    # theme(axis.text.y = element_text(color = species_colors[levels(factor(combined_ind_2025$Species))]))+
    theme(axis.text.y = element_text(color = species_colours_Rev3[levels(factor(combined_ind_Rev3$Species))]))+
    labs(y = "Species-group association (Reference: 2025 surveys)", x = "Species", title = "Revegetation sites"),
  ggplot(combined_ind_Rem3, aes(x = Species, y = `2025`, fill = Native_Status)) +
    geom_col(colour= "black") +
    coord_flip()+
    scale_fill_manual(values = c("Exotic" = "orange", "Native" = "darkgreen"))+
    theme_test() +
    # theme(axis.text.y = element_text(color = species_colors[levels(factor(combined_ind_2025$Species))]))+
    theme(axis.text.y = element_text(color = species_colours_Rem3[levels(factor(combined_ind_Rem3$Species))]))+
    labs(y = "Species-group association (Reference: 2025 surveys)", x = "Species", title = "Remnant sites")
  ,
  align ="hv", common.legend = TRUE, ncol=1
)
# ggsave(filename = "Veg-IndicatorSpecies-SurveyYears-Remnant_Revegetation-tall.pdf", path = outdir, width = 6.5, height = 10)





#### within veg comparisons ----------------------------------------------------
metadata_s$iLumpEco <- as.factor(metadata_s$iLumpEco)
eco_levels <- unique(metadata_s$tEcoName)
eco_levels
# functions
run_indicator_analysis <- function(eco, year = "2025") {
  
  idx <- metadata_s$tEcoName == eco & metadata_s$Survey_Year == year
  
  veg_sub  <- veg_species_wide_beta[idx, ]
  meta_sub <- metadata_s[idx, ]
  
  # Indicator species analysis
  multipatt.obj <- multipatt(x = veg_sub, cluster = meta_sub$RemRev, control = how(nperm = 999))
  
  # Strength of associations
  strassoc.obj <- strassoc(veg_sub, cluster = meta_sub$RemRev) |> as.data.frame()
  
  strassoc.obj$Species <- rownames(strassoc.obj)
  
  multipatt_sign <- multipatt.obj$sign
  multipatt_sign$Species <- rownames(multipatt_sign)
  
  combined_ind <- left_join(multipatt_sign, strassoc.obj, by = "Species")
  
  # Introduced / native status
  veg_data_intr <- veg_data_long %>%
    select(`SCIENTIFIC NAME`, INTRODUCED) %>%
    distinct()
  
  colnames(veg_data_intr) <- c("Species", "Exotic")
  
  combined_ind <- left_join(combined_ind, veg_data_intr, by = "Species")
  
  combined_ind <- combined_ind %>% 
    filter(p.value < 0.05) %>% 
    mutate(Species = factor(Species, levels = Species[order(Remnant)]), Native_Status = ifelse(Exotic == 1, "Exotic", "Native"), Exotic = as.factor(Exotic))
  
  # Colours for axis labels
  species_colors <- setNames(ifelse(combined_ind$Native_Status == "Exotic", "orange", "darkgreen"), combined_ind$Species)
  
  # Plot
  plot <- ggplot(combined_ind, aes(x = Species, y = Remnant, fill = Native_Status)) +
    geom_col(colour = "black") +
    coord_flip() +
    scale_fill_manual(values = c("Exotic" = "orange", "Native" = "darkgreen")) +
    theme_test() +
    theme(axis.text.y = element_text(color = species_colors[levels(combined_ind$Species)])) +
    labs(y = "Species-group association (Reference: Remnant)",x = "Species",title = paste("Eco", eco, "-", year, "Survey"))
  
  list(eco = eco,
       multipatt = multipatt.obj,
       combined = combined_ind,
       plot = plot)
}

indicator_results_2025 <- lapply(eco_levels, run_indicator_analysis, year = "2025")
names(indicator_results_2025) <- paste0("Eco", eco_levels)
indicator_results_2015 <- lapply(eco_levels, run_indicator_analysis, year = "2015")
names(indicator_results_2015) <- paste0("Eco", eco_levels)

eco_levels
ggpubr::ggarrange(
  indicator_results_2015$`EcoE. fasciculosa`$plot,
  indicator_results_2015$`EcoC. gracilis`$plot,
  indicator_results_2015$`EcoE. diversifolia`$plot,
  # indicator_results_2015$`EcoA. verticillata`$plot,
  indicator_results_2015$`EcoMixed mallee woodland`$plot,
  common.legend = TRUE, align = "hv"
)
# ggsave("Indicator-species-2015-surveys by ecosystemTYPE.pdf", path = outdir, height = 10, width = 15)

ggpubr::ggarrange(
  indicator_results_2025$`EcoE. fasciculosa`$plot,
  indicator_results_2025$`EcoC. gracilis`$plot,
  indicator_results_2025$`EcoE. diversifolia`$plot,
  # indicator_results_2025$`EcoA. verticillata`$plot,
  indicator_results_2025$`EcoMixed mallee woodland`$plot,
  common.legend = TRUE, align = "hv"
)
# ggsave("Indicator-species-2025-surveys by ecosystemTYPE.pdf", path = outdir, height = 10, width = 15)

# NEXT add detail about functional groups to see if different proportions of grasse and herbs are observed


# Comparisons to planting lists - WORK IN PROGRESS ------------------------------------------------

# Plot species as rows, treatments and sampling years as columns, and heatmap significant presence-absence outputs
# Prep planting data
planting_data <- read_excel("Survey_Data_2015.xlsx", sheet = "Planting lists")
# View(planting_data)
colnames(planting_data)

planting_data <- planting_data %>% 
  select(iWptID, Old_WptID, iPlanID, iYear_Planted, tNSXCode, SPECIES, iPlantsActual)

colnames(all_veg_data_v)
#  [1] "ProjectID"         "WptID"             "VisitID"           "NSXCODE"           "SCIENTIFIC NAME"   "INTRODUCED"       
#  [7] "tLF"               "tLFDesc"           "tAD"               "tADDesc"           "tLS"               "tLSDesc"          
# [13] "COVCODE"           "COVDESC"           "DeadCA"            "iPlantYear"        "iEcoID"            "tEcoName...18"    
# [19] "iLumpEco"          "tEcoName...20"     "tEcoName"          "tEcosystem_Desc"   "HILG"              "SITE_ID"          
# [25] "ogr_WptID"         "ogr_VisitI"        "ogr_PlanID"        "ogr_PolyID"        "Treatment"         "Control_site...9" 
# [31] "ogr_YearFi"        "ogr_Ecosys"        "ogr_Ecos_1"        "ogr_Eastin"        "ogr_Northi"        "Not_point"        
# [37] "Longitude"         "Latitude"          "Birdlife"          "Management"        "Control_site...20" "Survey_Year"      
# [43] "Site_ID" 

# Clean up (combine provenances etc)
planting_data2 <- planting_data %>%
  group_by(iWptID, iYear_Planted, iPlanID, SPECIES) %>%
  summarise(no_planted = sum(iPlantsActual))

planting_data2_joiner <- planting_data2
planting_data2_joiner$Survey_Year <- "Planting list"
planting_data2_joiner$presence <- 1

length(unique(planting_data2$iWptID))  # 60
length(unique(planting_data2$iPlanID)) # 55

# Prep survey data
colnames(all_veg_data_v)
colnames(planting_data2)

# get veg at site level
site_veg_data <- all_veg_data_v %>%
  distinct() %>%
  select(VisitID, WptID, VisitID, Survey_Year, `SCIENTIFIC NAME`, INTRODUCED, tEcoName, Birdlife, cover_midpoint) %>%
  group_by(Survey_Year, WptID, `SCIENTIFIC NAME`, INTRODUCED) %>%
  summarise(average_cover = mean(cover_midpoint)) %>%
  ungroup()

site_veg_met<- all_veg_data_v %>%
  distinct() %>%
  select(Survey_Year, WptID, tEcoName, iPlantYear) %>%
  distinct()

colnames(site_veg_data)
colnames(planting_data2)
intersect(colnames(site_veg_data), colnames(planting_data2))


unique(site_veg_data$`SCIENTIFIC NAME`)
unique(site_veg_data$Survey_Year)
planting_data3 <- left_join(site_veg_data, planting_data2, 
                            by = c("WptID" = "iWptID", 
                                   "SCIENTIFIC NAME" = "SPECIES"))

intersect(colnames(planting_data3), colnames(site_veg_met))

planting_data4 <- left_join(planting_data3, site_veg_met,
                            by = c("Survey_Year", "WptID"))
# View(planting_data4)

planting_data4$is_planted_sp <- ifelse(is.na(planting_data4$iPlanID) == "TRUE", 0, 1)
planting_data4$presence <- ifelse(planting_data4$average_cover > 0, 1, 0)
planting_data4$RemRev <- ifelse(planting_data4$iPlantYear == "Remnant", "Remnant", "Revegetation")

planting_data_combined <- planting_data4

planting_data_combined$Survey_treat <- paste0(planting_data_combined$Survey_Year, "_", planting_data_combined$RemRev)
planting_data_combined_REV <- subset(planting_data_combined, RemRev == "Revegetation")

# planting_data_combined_REVplanted <- subset(planting_data_combined_REV, is_planted_sp == 1)

# glimpse(planting_data_combined_REVplanted)
# Rows: 1,413
# Columns: 14
# $ Survey_Year       <fct> 2015, 2015, 2015, 2015, 2015, 2015, 2015, 2015, 2015, 2015, 2015, 2015, 2015, 2015, 2015, 2015, 2…
# $ WptID             <dbl> 576, 576, 576, 576, 576, 576, 576, 576, 576, 576, 576, 576, 576, 576, 576, 576, 577, 577, 577, 57…
# $ `SCIENTIFIC NAME` <chr> "Acacia dodonaeifolia", "Acacia myrtifolia", "Acacia paradoxa", "Acacia pycnantha", "Allocasuarin…
# $ INTRODUCED        <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, N…
# $ average_cover     <dbl> 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 0.500000, 9.416667, 0.50000…
# $ iYear_Planted     <dbl> 2014, 2014, 2014, 2014, 2014, 2014, 2014, 2014, 2014, 2014, 2014, 2014, 2014, 2014, 2014, 2014, 2…
# $ iPlanID           <dbl> 402, 402, 402, 402, 402, 402, 402, 402, 402, 402, 402, 402, 402, 402, 402, 402, 158, 158, 158, 15…
# $ no_planted        <dbl> 1861, 2554, 1025, 1852, 3499, 4020, 4683, 1197, 4538, 944, 3320, 300, 8443, 35, 700, 139, 25, 594…
# $ tEcoName          <chr> "E. fasciculosa", "E. fasciculosa", "E. fasciculosa", "E. fasciculosa", "E. fasciculosa", "E. fas…
# $ iPlantYear        <chr> "2014", "2014", "2014", "2014", "2014", "2014", "2014", "2014", "2014", "2014", "2014", "2014", "…
# $ is_planted_sp     <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
# $ presence          <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
# $ RemRev            <chr> "Revegetation", "Revegetation", "Revegetation", "Revegetation", "Revegetation", "Revegetation", "…
# $ Survey_treat      <chr> "2015_Revegetation", "2015_Revegetation", "2015_Revegetation", "2015_Revegetation", "2015_Reveget…


# 1) Should WptID have the species? ("WptID" should have a "PlanID" value with the "species" name present)
# 2) If so, does WptID actually have the species? ("WptID" should have a "presence" value of 1)
# 3) Therefore:
#       - If it should have the species, and it is found, assign correctly present
#       - If it should have the species, and it is not found, assign incorrectly absent
#       - If it should not have the species, and it is found, assign incorrectly present
#       - If it should not have the species, and it is not found, assign correctly absent

planting_data_combined_REV2 <- planting_data_combined_REV %>%
  mutate(presence_class = case_when(
    is_planted_sp == 1 & presence == 1 ~ "Correctly present",
    is_planted_sp == 1 & presence == 0 ~ "Incorrectly absent",
    is_planted_sp == 0 & presence == 1 ~ "Incorrectly present",
    is_planted_sp == 0 & presence == 0 ~ "Correctly absent",
    TRUE ~ NA_character_
  ))
table(planting_data_combined_REV2$presence_class)
# Correctly absent   Correctly present  Incorrectly absent Incorrectly present 
#               39                1405                   8                3869 

planting_data_combined_REV2 %>%
  count(Survey_treat, presence_class)
# # A tibble: 6 × 3
#   Survey_treat      presence_class          n
#   <chr>             <chr>               <int>
# 1 2015_Revegetation Correctly present     771
# 2 2015_Revegetation Incorrectly present  2042
# 3 2025_Revegetation Correctly absent       39
# 4 2025_Revegetation Correctly present     634
# 5 2025_Revegetation Incorrectly absent      8
# 6 2025_Revegetation Incorrectly present  1827

# Plot results
# species by Survey_treat
dropout_by_species <- planting_data_combined_REV2 %>%
  filter(is_planted_sp == 1) %>%
  count(`SCIENTIFIC NAME`, Survey_treat, presence_class) %>%
  group_by(`SCIENTIFIC NAME`, Survey_treat) %>%
  mutate(prop = n / sum(n))
# View(dropout_by_species)

ggplot(dropout_by_species %>% filter(presence_class == "Incorrectly absent"), 
       aes(x = reorder(`SCIENTIFIC NAME`, prop),  y = prop,  fill = Survey_treat)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(x = "Species", y = "Proportion incorrectly absent", fill = "Survey", title = "Planted species dropping out (incorrectly absent)") +
  theme_test()

species_presence_summary <- planting_data_combined_REV2 %>%
  count(`SCIENTIFIC NAME`, Survey_treat, presence_class)
species_presence_summary

ggplot(species_presence_summary, aes(x = `SCIENTIFIC NAME`, y = n, fill = presence_class)) +
  geom_col() +
  facet_grid(presence_class ~ Survey_treat, scales = "free_y") +
  coord_flip() +
  labs(x = "Species", y = "Number of records", fill = "Presence class", title = "Species presence outcomes by survey year") +
  theme_test()

species_presence_prop <- species_presence_summary %>%
  group_by(`SCIENTIFIC NAME`, Survey_treat) %>%
  mutate(prop = n / sum(n))
species_presence_prop

ggplot(species_presence_prop %>% filter(!presence_class == "Correctly absent")
       , aes(x = `SCIENTIFIC NAME`, y = prop, fill = presence_class )) +
  geom_col() +
  facet_grid(presence_class ~ Survey_treat, scales = "free_y") +
  coord_flip() +
  labs(x = "Species", y = "Proportion of sites where speces was found", fill = "Presence class", title = "Proportional presence outcomes by species") +
  theme_test()

species_dropout_index <- planting_data_combined_REV2 %>%
  filter(is_planted_sp == 1) %>%
  group_by(`SCIENTIFIC NAME`, Survey_treat) %>%
  summarise(
    n_sites = n(),
    dropout_rate = mean(presence == 0),
    .groups = "drop"
  )

# Presence based drop out rates for species
ggplot(
  species_dropout_index %>% filter(dropout_rate > 0),
  aes(x = reorder(`SCIENTIFIC NAME`, dropout_rate), y = dropout_rate)) +
  geom_col(position = "dodge",colour = "black") +
  coord_flip() +
  labs(x = "Species", y = "Dropout rate", fill = "Survey") +
  theme_test()+
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(aes(label = paste0("n = ", n_sites)),  hjust = -0.25) 

# Abundance based declines in planted species ----------------------------------

## Parametric bootstrap binomial glmem ------------------------------------------
library(lme4)
library(pbkrtest)

### Remnant veg comparisons across years  --------------------------------------

species_delta_df <- all_veg_data_v %>%
  distinct() %>%
  select(WptID, Survey_Year, VisitID, `SCIENTIFIC NAME`, INTRODUCED, tEcoName, Birdlife, cover_midpoint)

planting_data_large <- left_join(species_delta_df, planting_data2, 
                            by = c("WptID" = "iWptID", 
                                   "SCIENTIFIC NAME" = "SPECIES"))
site_veg_met<- all_veg_data_v %>%
  distinct() %>%
  select(Survey_Year, WptID, tEcoName, iPlantYear) %>%
  distinct()

planting_data_large <- left_join(planting_data_large, site_veg_met, 
                                 by = c("WptID", "Survey_Year", "tEcoName"))

planting_data_large$RemRev <- ifelse(planting_data_large$iPlantYear == "Remnant", "Remnant", "Revegetation")
planting_data_large$Survey_Year <- as.factor(planting_data_large$Survey_Year)

planting_data_combined_revonly <- subset(planting_data_large, RemRev == "Revegetation")
length(unique(planting_data_combined_revonly$`SCIENTIFIC NAME`)) # 446 unique species...

spxsp_split_reveg <- split(planting_data_combined_revonly, planting_data_combined_revonly$`SCIENTIFIC NAME`)
length(spxsp_split_reveg)

# analysis of species
spcmodel_listRemnanta <- lapply(spxsp_split_reveg, function(data) {
  
  tryCatch({
    
    full <- lmer(cover_midpoint ~ Survey_Year + (1 | WptID), data = data, 
                 control = lmerControl(optimizer = "bobyqa",  optCtrl = list(maxfun = 1e5)))
    
    null <- lmer(cover_midpoint ~ 1 + (1 | WptID), data = data, 
                  control = lmerControl(optimizer = "bobyqa",  optCtrl = list(maxfun = 1e5)))
    
    pb <- PBmodcomp(full, null, nsim = 999)

    print(pb)
    
    return(list(
      species = unique(data$`SCIENTIFIC NAME`),
      p_value = pb$test$p.value,
      full_model = full)
    )
    
  }, error = function(e) {
    message("Failed for species ", unique(data$tSpp), ": ", e$message)
    return(NULL)
  })
})

# Drop NULLs (failed models (ie. rare species))
spcmodel_listRemnantb <- Filter(Negate(is.null), spcmodel_listRemnanta)

View(spcmodel_listRemnantb[1])

# Turn list into a dataframe
df_tot_Remnantb <- data.frame()

for (species in names(spcmodel_listRemnantb)) {
  mod <- spcmodel_listRemnantb[[species]]
  
  p_val <- mod$p_value[2]
  estimate <- fixef(mod$full_model)["Survey_Year2025"]
  
  tmp <- data.frame(
    Species = species,
    Estimate = estimate,
    P_value = p_val
  )
  
  df_tot_Remnantb <- rbind(df_tot_Remnantb, tmp)
}

# Plotting and interpretation
df_tot_Remnantb$Significant <- ifelse(df_tot_Remnantb$P_value < 0.05, "P <0.05", "NS")
df_tot_Remnantb$Species <- factor(df_tot_Remnantb$Species, levels = df_tot_Remnantb$Species[order(df_tot_Remnantb$Estimate, decreasing = FALSE)])

df_tot_Remnantb$P_adj <- p.adjust(df_tot_Remnantb$P_value, method = "BH")
df_tot_Remnantb$Significance <- ifelse(df_tot_Remnantb$P_value < 0.05, "P < 0.05", "NS")
df_tot_Remnantb$Significance <- ifelse(df_tot_Remnantb$P_adj < 0.05, "Padj < 0.05", df_tot_Remnantb$Significance)
df_tot_Remnantb$Significance <- factor(df_tot_Remnantb$Significance, levels = c("Padj < 0.05", "P < 0.05", "NS"))


plot_Remnant <- ggplot(df_tot_Remnantb[df_tot_Remnantb$Significance == "Padj < 0.05" & df_tot_Remnantb$Estimate < 0,], aes(x = Species, y = Estimate, fill = Estimate))+
  geom_col(colour = "black")+coord_flip()+
  ggtitle("Remnant sites only")+
  theme_test()+
  labs(y = "Net change in percentage cover \nrelative to 2015 surveys", x = "Species", title="Remnant sites only")+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
plot_Remnant

# calculate absolute and relative change of species abundance from 2015
abundance_change <- planting_data_combined_REV2 %>%
  filter(is_planted_sp == 1) %>%
  select(WptID, `SCIENTIFIC NAME`, Survey_Year, average_cover) %>%
  pivot_wider(
    names_from = Survey_Year,
    values_from = average_cover,
    names_prefix = "abundance_"
  ) %>%
  mutate(
    abundance_2015 = ifelse(is.na(abundance_2015), 0, abundance_2015),
    abundance_2025 = ifelse(is.na(abundance_2025), 0, abundance_2025),
    absolute_change = abundance_2025 - abundance_2015,
    
    normalised_change = absolute_change / (abundance_2025+abundance_2015),
    
    relative_change = ifelse(abundance_2015 > 0, absolute_change / abundance_2015, NA_real_)
    # dropped_out = abundance_2025 < 0.1 & abundance_2015 >= 0.1,  # example threshold for dropout,
  )

abundance_dropout_summary <- abundance_change %>%
  group_by(`SCIENTIFIC NAME`) %>%
  summarise(
    n_sites = n(), # how many times the species is counted (across sites)
    # dropout_sites = sum(dropped_out, na.rm = TRUE),
    # dropout_rate = dropout_sites / n_sites,
    mean_relative_change = mean(relative_change, na.rm = TRUE),
    mean_normalised_change = mean(normalised_change, na.rm = TRUE)
  )
ggplot(abundance_dropout_summary, aes(x = reorder(`SCIENTIFIC NAME`, mean_normalised_change), y = mean_normalised_change))+
  geom_col()+
  coord_flip() +
  theme_bw()

# Management analysis ----------------------------------------------------------
grazing_info_only

unique(management_2025$`Domestic grazing`)

library(stringr)

management_2025 <- management_2025 %>%
  mutate(level = case_when(
    str_detect(`Domestic grazing`, "Low")    ~ "Low",
    str_detect(`Domestic grazing`, "Medium") ~ "Medium",
    str_detect(`Domestic grazing`, "High")   ~ "High",
    TRUE                            ~ NA_character_
  ))

