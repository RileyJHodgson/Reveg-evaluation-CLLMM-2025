# Analysis of Goyder - CLLMM revegetation project (bird data)
# RJH 2026

# Set working directory --------------------------------------------------------
setwd("/PATH/TO/DATA")

# libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(vegan)
library(stringr)

# create output directory
outdir <- "/PATH/TO/DATA/R_output_birds"
# dir.create(file.path(outdir))
outdir <- file.path(outdir)

# Read data
all_bird_data <- read_excel("Joint_Bird_Data.xlsx")
colnames(all_bird_data)
#SITE SPECIFIC METADATA
metadata_v <- read_excel("Goyder-soil_analysis-metadata.xlsx", sheet = "site-info")
colnames(metadata_v)

# Assuming bInOut is whether Birds are in or outside of the survey boundaries
all_bird_data <- subset(all_bird_data, bInOut == 1)
all_bird_data$RemRev <- ifelse(all_bird_data$iPlantYear == "Remnant", "Remnant", "Revegetation")
all_bird_data$daDate <- as.Date(all_bird_data$daDate, "%Y-%m-%d")
class(all_bird_data$daDate)
colnames(all_bird_data)
# [1] "Survey_Year" "daDate"      "Old_WptID"   "WptID"       "iVisitID"    "Treat_type"  "tNSXCode"    "tSpp"        "tName"       "bInOut"      "iPlantYear"  "tVisitNotes" "RemRev"     

# Get bird richness from each Site - Visit combination
all_bird_data2 <- all_bird_data %>%
  select(Survey_Year, WptID, iVisitID, Treat_type, tNSXCode, tSpp, tName, iPlantYear) %>%
  distinct()

nrow(all_bird_data)  # 2078
nrow(all_bird_data2) # 2070

head(all_bird_data2)
# # A tibble: 6 × 8
# Survey_Year WptID iVisitID Treat_type tNSXCode tSpp                         tName                  iPlantYear
#         <dbl> <dbl>    <dbl> <chr>      <chr>    <chr>                        <chr>                  <chr>     
# 1        2015   689      525 Remnant    C00273   Eolophus roseicapillus       Galah                  Remnant   
# 2        2015   689      525 Remnant    K04149   Pachycephala rufiventris     Rufous Whistler        Remnant   
# 3        2015   689      525 Remnant    U04126   Phylidonyris novaehollandiae New Holland Honeyeater Remnant   
# 4        2015   689      525 Remnant    C00361   Rhipidura albiscapa          Grey Fantail           Remnant   
# 5        2015   689      525 Remnant    S00465   Smicrornis brevirostris      Weebill                Remnant   
# 6        2015   689      525 Remnant    U00254   Trichoglossus haematodus     Rainbow Lorikeet       Remnant   

# Richness - alpha diversity ---------------------------------------------------
richness_data <- all_bird_data %>%
  group_by(Survey_Year, WptID, iVisitID, Treat_type, iPlantYear) %>%
  reframe(Richness = length(unique(tSpp)))
head(richness_data)

richness_data$iPlantYear <- factor(richness_data$iPlantYear, levels = c("2011", "2012", "2013", "2014", "2015", "Remnant"))
richness_data$Treat_type <- factor(richness_data$Treat_type, levels = c("Revegetated", "Remnant"))

richness_data_average <- richness_data %>%
  group_by(Survey_Year, Treat_type, WptID, iPlantYear) %>%
  summarise(mean_Richness = mean(Richness))
head(richness_data_average)

# plot
ggplot(richness_data, aes(x=iPlantYear, y=Richness, fill = iPlantYear))+
  geom_violin()+
  geom_boxplot(outlier.shape = NA, width = 0.15, fill = "white")+
  geom_jitter(width = 0.1, height = 0, alpha = 0.125)+
  theme_test()

ggplot(richness_data, aes(x=Treat_type, y=Richness, fill = Treat_type))+
  geom_violin()+
  geom_boxplot(outlier.shape = NA, width = 0.15, fill = "white")+
  geom_jitter(width = 0.1, height = 0, alpha = 0.125)+
  theme_test()+
  labs(y= "Bird richness", x = "Treatment")

library(ggbeeswarm)
ggplot(richness_data, aes(x=Treat_type, y=Richness, fill = Treat_type))+
  geom_violin()+
  geom_boxplot(outlier.shape = NA, width = 0.15, fill = "white")+
  geom_beeswarm(alpha = 0.25, cex = 1.75)+
  theme_test()+
  labs(y= "Bird richness", x = "Treatment")

### stats  ---------------------------------------------------------------------
library(lme4)
poisson_model <- glmer(Richness ~ Treat_type + (1|WptID), family = poisson(), data = richness_data)
car::Anova(poisson_model)
# Analysis of Deviance Table (Type II Wald chisquare tests)
# Response: Richness
#             Chisq Df Pr(>Chisq)   
# Treat_type 8.1269  1   0.004361 **

# Useful plot
ggplot(richness_data_average, aes(x=Treat_type, y=mean_Richness, fill = Treat_type))+
  geom_violin()+
  geom_boxplot(outlier.shape = NA, width = 0.15, fill = "white")+
  geom_jitter(width = 0.1, height = 0, alpha = 0.125)+
  theme_test()+
  labs(y= "Mean bird richness", x = "Treatment")+
  facet_grid(~Survey_Year)
# ggsave(filename = "Bird-Richness-SurveyYear_RemRev.pdf", path = outdir, width = 7, height = 6)

# Beta - Jaccard ordination differences across sites-visits ---------------------------
all_bird_data$One_ID_Var <- paste0(all_bird_data$Survey_Year, "_G", all_bird_data$WptID, "_", all_bird_data$iVisitID, "_py", all_bird_data$iPlantYear)

wide_bird_data <- all_bird_data %>%
  select(One_ID_Var, tSpp) %>%
  mutate(present = 1) %>%
  distinct() %>%
  pivot_wider(
    id_cols = One_ID_Var,             
    names_from = tSpp,
    values_from = present,
    values_fill = 0) %>% 
  as.data.frame()

rownames(wide_bird_data) <- wide_bird_data$One_ID_Var

wide_bird_data <- wide_bird_data %>%
  select(-One_ID_Var)
wide_bird_data[1:5,1:5]

# Site level only
all_bird_data$Site_ID_Var <- paste0(all_bird_data$Survey_Year, "_G", all_bird_data$WptID, "_py", all_bird_data$iPlantYear)

wide_bird_data_site <- all_bird_data %>%
  mutate(present = 1) %>%
  select(Site_ID_Var, iPlantYear, tSpp, present) %>%
  group_by(Site_ID_Var, tSpp) %>%
  distinct() %>%
  pivot_wider(
    id_cols = Site_ID_Var,             
    names_from = tSpp,
    values_from = present,
    values_fill = 0) %>% 
  as.data.frame()

rownames(wide_bird_data_site) <- wide_bird_data_site$Site_ID_Var

wide_bird_data_site <- wide_bird_data_site %>%
  select(-Site_ID_Var)
wide_bird_data_site[1:5,1:5]

# Plot
set.seed(123)
NMDS_Jacc <- metaMDS(wide_bird_data, distance = "jaccard", binary = TRUE)
NMDS_Jacc$stress # 0.1290578

plot(NMDS_Jacc)

# Plot site level
# Plot
set.seed(123)
NMDS_Jacc <- metaMDS(wide_bird_data_site, distance = "jaccard", binary = TRUE)
NMDS_Jacc$stress # 0.2199447

plot(NMDS_Jacc)

points_Jacc <- NMDS_Jacc$points

points_Jacc <- as.data.frame(points_Jacc)
points_Jacc$Site_ID_Var <- rownames(points_Jacc)

# Attach metadata

all_bird_site_beta <- all_bird_data %>% 
  select(Site_ID_Var, Survey_Year, WptID, iPlantYear, Treat_type) %>%
  distinct()
head(all_bird_site_beta)
head(points_Jacc)

bird_metadata_s_Jacc <- left_join(all_bird_site_beta, points_Jacc, by = "Site_ID_Var")
bird_metadata_s_Jacc$Survey_Year <- as.factor(bird_metadata_s_Jacc$Survey_Year)

colnames(bird_metadata_s_Jacc)
# "Site_ID"      "Survey_Year"  "iWptID"       "iPlantYear"   "iTreatID"     "iEcosystemID"
# "tEcoName"     "MDS1"         "MDS2"
head(bird_metadata_s_Jacc)

ggplot(bird_metadata_s_Jacc, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =Treat_type))+
  geom_point(size= 2)+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()
# ggsave(filename = "Bird-BETA-NMDS-Jaccard-SurveyYear_RemRev.pdf", path = outdir, width = 7, height = 6)

ggplot(bird_metadata_s_Jacc, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =iPlantYear))+
  geom_point(size= 2)+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()
# ggsave(filename = "Bird-BETA-NMDS-Jaccard-SurveyYear_PlantYear.pdf", path = outdir, width = 7, height = 6)

# stats
dist_mat_jaccard_bird <- vegdist(wide_bird_data, method = "jaccard")
# labels(dist_mat_bray_all)

# I need to control for two variables in my data: ecosystem type (tEcoName), and repeated sites (iWptID)
# metadata_s$strata_ecosys_wpt <- interaction(
#   metadata_s$tEcoName,
#   metadata_s$iWptID,
#   drop = TRUE
#   )
all_bird_visit_beta <- all_bird_data %>% 
  select(One_ID_Var, Survey_Year, WptID, iPlantYear, Treat_type) %>%
  distinct()

set.seed(123)
adonis2(dist_mat_jaccard_bird ~ Survey_Year * Treat_type, data = all_bird_visit_beta,
        strata = all_bird_visit_beta$WptID
        )
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Blocks:  strata 
# Permutation: free
# Number of permutations: 999
# adonis2(formula = dist_mat_jaccard_bird ~ Survey_Year * Treat_type, data = all_bird_visit_beta, strata = all_bird_visit_beta$WptID)
#                         Df SumOfSqs      R2       F Pr(>F)    
# Survey_Year              1    2.559 0.02058  6.3949  0.001 ***
# Treat_type               1    4.184 0.03365 10.4554  0.001 ***
# Survey_Year:Treat_type   1    1.543 0.01241  3.8568  0.001 ***
# Residual               290  116.047 0.93336                   
# Total                  293  124.334 1.00000  
ncol(all_bird_visit_beta)

# Species level comparisons ----------------------------------------------------
length(unique(all_bird_data$tSpp)) # 124

spxsp_df <- all_bird_data
spxsp_df$Presence <- 1

# Species names and common names
spxsp_df_spp <- all_bird_data %>% 
  select(tSpp, tName) %>%
  distinct()
  
metadata_v$Old_WptID <- as.numeric(metadata_v$Old_WptID)
spxsp_df <- left_join(spxsp_df, metadata_v, by = c("WptID" = "iWptID", "Old_WptID"))
colnames(spxsp_df)

# add all 0s
sampling_frame <- spxsp_df %>%
  distinct(Survey_Year, WptID, iVisitID)

species_pool <- spxsp_df %>%
  distinct(tSpp)

site_meta <- spxsp_df %>%
  distinct(
    WptID,
    Survey_Year, RemRev, tEcoName,
    Treat_type, Site_ID, Latitude, Longitude,
    Management)

spxsp_df_pa <- sampling_frame %>%
  tidyr::crossing(species_pool) %>%
  left_join(
    spxsp_df %>%
      select(Survey_Year, WptID, iVisitID, tSpp, Presence) %>%
      mutate(Presence = 1),
    by = c("Survey_Year", "WptID", "iVisitID", "tSpp")
  ) %>%
  mutate(Presence = if_else(is.na(Presence), 0L, 1L)) %>%
  left_join(site_meta, by = c("Survey_Year", "WptID"))

colnames(spxsp_df_pa)

# run functions on all species in loop.
library(lme4)

GLMEM_binomial_permute_anova <- function(data, response, fixed_effects, random_effects, nreps = 9999, seed = 123) {
  
  set.seed(seed)
  
  # Build fixed effects
  fixed_formula <- if (length(fixed_effects) == 1) {
    fixed_effects
  } else {
    paste(fixed_effects, collapse = " + ")
  }
  # Build random effects
  # random_formula <- paste("(1 | ", random_effects, ")", collapse = " + ")
  random_formula <- paste("(1 |", random_effects, ")", collapse = " + ")
  
  # Full model formula
  full_formula <- as.formula(paste(response, "~", fixed_formula, "+", random_formula))
  
  # Fit observed model
  withCallingHandlers(
    model <- glmer(full_formula, family = binomial(), data = data),
    warning = function(w) {
      if (grepl("singular", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
  # model <- glmer(full_formula, family = binomial(), data = data)
  anova_obs <- car::Anova(model)
  model_summary <- summary(model)
  fixed_names <- rownames(anova_obs)
  observed_chi <- anova_obs$Chisq
  
  # Storage for permutation stats
  permute_teststats <- matrix(
    NA, nrow = nreps, ncol = length(fixed_names),
    dimnames = list(NULL, fixed_names)
  )
  
  N <- nrow(data)
  
  for (i in seq_len(nreps)) {
    data_perm <- data
    data_perm[[response]] <- sample(data[[response]], N, replace = FALSE)
    
    mod_perm <- tryCatch(
      suppressWarnings(
        glmer(full_formula, family = binomial(), data = data_perm)
      ),
      error = function(e) NULL
    )
    # mod_perm <- tryCatch(
    #   glmer(full_formula, family = binomial(), data = data_perm),
    #   error = function(e) {
    #     message("Permutation ", i, " failed: ", e$message)
    #     NULL
    #   }
    # )
    
    if (!is.null(mod_perm)) {
      anova_perm <- tryCatch(
        car::Anova(mod_perm),
        error = function(e) {
          message("Anova failed on permutation ", i, ": ", e$message)
          NULL
        }
      )
      
      if (!is.null(anova_perm) && all(fixed_names %in% rownames(anova_perm))) {
        permute_teststats[i, ] <- anova_perm[fixed_names, "Chisq"]
      }
    }
  }
  # Drop failed permutations, keeping matrix structure
  permute_teststats <- permute_teststats[complete.cases(permute_teststats), , drop = FALSE]
  
  # Then calculate permutation p-values safely
  p_values <- sapply(seq_along(fixed_names), function(j) {
    nbout <- which(permute_teststats[, j] >= observed_chi[j])
    p_perm <- length(nbout) / nrow(permute_teststats)
    return(p_perm)
  })
  
  names(p_values) <- fixed_names
  
  return(list(
    formula = full_formula,
    model = model,
    Anova = anova_obs,
    Summary = model_summary,
    observed_chisq = observed_chi,
    permuted_chisq = permute_teststats,
    permutation_p_values = p_values
  ))
}

# 2025 data
spxsp_df_2025 <- subset(spxsp_df_pa, Survey_Year == "2025")

spcmodel_list2025 <- list()

for (species in unique(spxsp_df_2025$tSpp)) {
  data <- spxsp_df_2025 %>%
    filter(tSpp == species)
  
  # skip rare species
  if (sum(data$Presence == 1, na.rm = TRUE) < 10) next
  # skip if no variation in RemRev
  if (length(unique(data$RemRev[data$Presence == 1])) < 2) next
  
  gmodel <- GLMEM_binomial_permute_anova(data = data, 
                                response = "Presence", 
                                fixed_effects = "RemRev", 
                                random_effects = c("WptID"), 
                                nreps = 999, seed = 123
                                )

  # if (inherits(gmodel, "try-error")) next
  spcmodel_list2025[[species]] <- gmodel
}

# saveRDS(object = spcmodel_list2025, file = "spcmodel_list2025.RDS")
# spcmodel_list2025 <- readRDS("spcmodel_list2025.RDS")

df_tot_2025 <- data.frame()

for (species in seq_along(names(spcmodel_list2025))) {
  df <- as.data.frame(spcmodel_list2025[[species]]$Summary$coefficients)
  p_val <- spcmodel_list2025[[species]]$permutation_p_values
  
  tmp <- data.frame(
    Species = names(spcmodel_list2025)[[species]],
    Estimate = df$Estimate[2],  # assuming the 2nd row is the coefficient of interest
    P_value = p_val
  )
  df_tot_2025 <- rbind(df_tot_2025, tmp)
}

df_tot_2025 <- left_join(df_tot_2025, spxsp_df_spp, by = c("Species" = "tSpp"))
df_tot_2025$Species <- factor(df_tot_2025$Species, levels = df_tot_2025$Species[order(df_tot_2025$Estimate, decreasing = FALSE)])
df_tot_2025$tName <- factor(df_tot_2025$tName, levels = df_tot_2025$tName[order(df_tot_2025$Estimate, decreasing = FALSE)])
df_tot_2025$P_adj <- p.adjust(df_tot_2025$P_value, method = "BH")
df_tot_2025$Significance <- ifelse(df_tot_2025$P_value < 0.05, "P < 0.05", "NS")
df_tot_2025$Significance <- ifelse(df_tot_2025$P_adj < 0.05, "Padj < 0.05", df_tot_2025$Significance)
df_tot_2025$Significance <- factor(df_tot_2025$Significance, levels = c("Padj < 0.05", "P < 0.05", "NS"))

colnames(df_tot_2025)

# ggplot
df_tot_2025$Common_species <- paste0(df_tot_2025$tName, " | ", df_tot_2025$Species)
df_tot_2025$Common_species <- factor(df_tot_2025$Common_species, levels = df_tot_2025$Common_species[order(df_tot_2025$Estimate, decreasing = FALSE)])

# ggplot
birdchange_2025 <- ggplot(data = df_tot_2025, aes(x = tName, y = Estimate, fill = Significance))+
  geom_col()+coord_flip()+ggtitle("2025 surveys only")+
  scale_fill_manual(values = c("Padj < 0.05" = "darkblue", "P < 0.05" = "blue", "NS" = "lightgrey"))+
  theme_test()+
  labs(y = "Log odds ratio of species presence:\n Remnant (reference group) vs Revegetation", x = "Bird species")
birdchange_2025
birdchange_2025_sig <- ggplot(data = df_tot_2025[df_tot_2025$Significance == "P < 0.05" | df_tot_2025$Significance == "Padj < 0.05",], aes(x = tName, y = Estimate, fill = Significance))+
  geom_col()+coord_flip()+ggtitle("2025 surveys only")+
  scale_fill_manual(values = c("Padj < 0.05" = "darkblue", "P < 0.05" = "blue", "NS" = "lightgrey"))+
  theme_test()+
  labs(y = "Log odds ratio of species presence:\n Remnant (reference group) vs Revegetation", x = "Bird species")
birdchange_2025_sig

# 2015
spxsp_df_2015 <- subset(spxsp_df_pa, Survey_Year == "2015")

spcmodel_list2015 <- list()

for (species in unique(spxsp_df_2015$tSpp)) {
  data <- spxsp_df_2015 %>%
    filter(tSpp == species)
  
  # skip rare species
  if (sum(data$Presence == 1, na.rm = TRUE) < 10) next
  # skip if no variation in RemRev
  if (length(unique(data$RemRev[data$Presence == 1])) < 2) next
  
  gmodel <- GLMEM_binomial_permute_anova(data = data, 
                                         response = "Presence", 
                                         fixed_effects = "RemRev", 
                                         random_effects = c("WptID"), 
                                         nreps = 999, seed = 123
  )
  
  # if (inherits(gmodel, "try-error")) next
  spcmodel_list2015[[species]] <- gmodel
}

# saveRDS(object = spcmodel_list2015, file = "spcmodel_list2015.RDS")
# spcmodel_list2015 <- readRDS("spcmodel_list2015.RDS")

df_tot_2015 <- data.frame()

for (species in seq_along(names(spcmodel_list2015))) {
  df <- as.data.frame(spcmodel_list2015[[species]]$Summary$coefficients)
  p_val <- spcmodel_list2015[[species]]$permutation_p_values
  
  tmp <- data.frame(
    Species = names(spcmodel_list2015)[[species]],
    Estimate = df$Estimate[2],  # assuming the 2nd row is the coefficient of interest
    P_value = p_val
  )
  df_tot_2015 <- rbind(df_tot_2015, tmp)
}

df_tot_2015 <- left_join(df_tot_2015, spxsp_df_spp, by = c("Species" = "tSpp"))
df_tot_2015$Species <- factor(df_tot_2015$Species, levels = df_tot_2015$Species[order(df_tot_2015$Estimate, decreasing = FALSE)])
df_tot_2015$tName <- factor(df_tot_2015$tName, levels = df_tot_2015$tName[order(df_tot_2015$Estimate, decreasing = FALSE)])
df_tot_2015$P_adj <- p.adjust(df_tot_2015$P_value, method = "BH")
df_tot_2015$Significance <- ifelse(df_tot_2015$P_value < 0.05, "P < 0.05", "NS")
df_tot_2015$Significance <- ifelse(df_tot_2015$P_adj < 0.05, "Padj < 0.05", df_tot_2015$Significance)
df_tot_2015$Significance <- factor(df_tot_2015$Significance, levels = c("Padj < 0.05", "P < 0.05", "NS"))


# ggplot
df_tot_2015$Common_species <- paste0(df_tot_2015$tName, " | ", df_tot_2015$Species)
df_tot_2015$Common_species <- factor(df_tot_2015$Common_species, levels = df_tot_2015$Common_species[order(df_tot_2015$Estimate, decreasing = FALSE)])

# ggplot
birdchange_2015 <- ggplot(data = df_tot_2015, aes(x = tName, y = Estimate, fill = Significance))+
  geom_col()+coord_flip()+ggtitle("2015 surveys only")+
  scale_fill_manual(values = c("Padj < 0.05" = "darkblue", "P < 0.05" = "blue", "NS" = "lightgrey"))+
  theme_test()+
  labs(y = "Log odds ratio of species presence:\n Remnant (reference group) vs Revegetation", x = "Bird species")

birdchange_2015_sig <- ggplot(data = df_tot_2015[df_tot_2015$Significance == "P < 0.05" | df_tot_2015$Significance == "Padj < 0.05",], aes(x = tName, y = Estimate, fill = Significance))+
  geom_col()+coord_flip()+ggtitle("2015 surveys only")+
  scale_fill_manual(values = c("Padj < 0.05" = "darkblue", "P < 0.05" = "blue", "NS" = "lightgrey"))+
  theme_test()+
  labs(y = "Log odds ratio of species presence:\n Remnant (reference group) vs Revegetation", x = "Bird species")


# Reveg sites only
spxsp_df_Reveg <- subset(spxsp_df_pa, RemRev == "Revegetation")
spxsp_df_Reveg$Survey_Year <- as.factor(spxsp_df_Reveg$Survey_Year)
levels(spxsp_df_Reveg$Survey_Year)

spcmodel_listReveg <- list()

spxsp_split <- split(spxsp_df_Reveg, spxsp_df_Reveg$tSpp)
spcmodel_listReveg <- lapply(spxsp_split, function(data) {
  
  if (sum(data$Presence == 1, na.rm = TRUE) < 10) return(NULL)
  if (length(unique(data$Survey_Year[data$Presence == 1])) < 2) return(NULL)
  
  tryCatch(
    suppressMessages(
      GLMEM_binomial_permute_anova(
        data = data,
        response = "Presence",
        fixed_effects = "Survey_Year",
        random_effects = c("WptID"),
        nreps = 99,
        seed = 123
      )
    ),
    error = function(e) NULL
  )
})

# Drop NULLs
spcmodel_listReveg <- Filter(Negate(is.null), spcmodel_listReveg)

# saveRDS(object = spcmodel_listReveg, file = "spcmodel_listReveg.RDS")
# spcmodel_listReveg <- readRDS("spcmodel_listReveg.RDS")

df_tot_Reveg <- data.frame()

for (species in seq_along(names(spcmodel_listReveg))) {
  df <- as.data.frame(spcmodel_listReveg[[species]]$Summary$coefficients)
  p_val <- spcmodel_listReveg[[species]]$permutation_p_values
  
  tmp <- data.frame(
    Species = names(spcmodel_listReveg)[[species]],
    Estimate = df$Estimate[2],  # assuming the 2nd row is the coefficient of interest
    P_value = p_val,
    Significance = ifelse(p_val < 0.05, "P < 0.05", "NS")
  )
  df_tot_Reveg <- rbind(df_tot_Reveg, tmp)
}

df_tot_Reveg <- left_join(df_tot_Reveg, spxsp_df_spp, by = c("Species" = "tSpp"))
df_tot_Reveg$Species <- factor(df_tot_Reveg$Species, levels = df_tot_Reveg$Species[order(df_tot_Reveg$Estimate, decreasing = FALSE)])
df_tot_Reveg$tName <- factor(df_tot_Reveg$tName, levels = df_tot_Reveg$tName[order(df_tot_Reveg$Estimate, decreasing = FALSE)])
df_tot_Reveg$P_adj <- p.adjust(df_tot_Reveg$P_value, method = "BH")
df_tot_Reveg$Significance <- ifelse(df_tot_Reveg$P_value < 0.05, "P < 0.05", "NS")
df_tot_Reveg$Significance <- ifelse(df_tot_Reveg$P_adj < 0.05, "Padj < 0.05", df_tot_Reveg$Significance)
df_tot_Reveg$Significance <- factor(df_tot_Reveg$Significance, levels = c("Padj < 0.05", "P < 0.05", "NS"))

# ggplot
df_tot_Reveg$Common_species <- paste0(df_tot_Reveg$tName, " | ", df_tot_Reveg$Species)
df_tot_Reveg$Common_species <- factor(df_tot_Reveg$Common_species, levels = df_tot_Reveg$Common_species[order(df_tot_Reveg$Estimate, decreasing = FALSE)])

birdchange_rev <- ggplot(data = df_tot_Reveg, aes(x = tName, y = Estimate, fill = Significance))+
  geom_col()+coord_flip()+
  ggtitle("Revegetation only")+
  scale_fill_manual(values = c("Padj < 0.05" = "darkblue", "P < 0.05" = "blue", "NS" = "lightgrey"))+
  theme_test()+
  labs(y = "Log odds ratio of species presence:\n 2015 (reference group) vs 2025", x = "Bird species")

birdchange_rev_sig <- ggplot(data = df_tot_Reveg[df_tot_Reveg$Significance == "P < 0.05" | df_tot_Reveg$Significance == "Padj < 0.05",], aes(x = tName, y = Estimate, fill = Significance))+
  geom_col()+coord_flip()+
  scale_fill_manual(values = c("Padj < 0.05" = "darkblue", "P < 0.05" = "blue", "NS" = "lightgrey"))+
  ggtitle("Revegetation only")+
  theme_test()+
  labs(y = "Log odds ratio of species presence:\n 2015 (reference group) vs 2025", x = "Bird species")

table(df_tot_Reveg$Significance)

# Remnant sites only
spxsp_df_Remnant <- subset(spxsp_df_pa, RemRev == "Remnant")
spxsp_df_Remnant$Survey_Year <- as.factor(spxsp_df_Remnant$Survey_Year)
levels(spxsp_df_Remnant$Survey_Year)

spcmodel_listRemnant <- list()

spxsp_split <- split(spxsp_df_Remnant, spxsp_df_Remnant$tSpp)
spcmodel_listRemnant <- lapply(spxsp_split, function(data) {
  
  if (sum(data$Presence == 1, na.rm = TRUE) < 10) return(NULL)
  if (length(unique(data$Survey_Year[data$Presence == 1])) < 2) return(NULL)
  
  tryCatch(
    suppressMessages(
      GLMEM_binomial_permute_anova(
        data = data,
        response = "Presence",
        fixed_effects = "Survey_Year",
        random_effects = c("WptID"),
        nreps = 99,
        seed = 123
      )
    ),
    error = function(e) NULL
  )
})

# Drop NULLs
spcmodel_listRemnant <- Filter(Negate(is.null), spcmodel_listRemnant)

# saveRDS(object = spcmodel_listRemnant, file = "spcmodel_listRemnant.RDS")
# spcmodel_listRemnant <- readRDS("spcmodel_listRemnant.RDS")

df_tot_Remnant <- data.frame()

for (species in seq_along(names(spcmodel_listRemnant))) {
  df <- as.data.frame(spcmodel_listRemnant[[species]]$Summary$coefficients)
  p_val <- spcmodel_listRemnant[[species]]$permutation_p_values
  
  tmp <- data.frame(
    Species = names(spcmodel_listRemnant)[[species]],
    Estimate = df$Estimate[2],  # assuming the 2nd row is the coefficient of interest
    P_value = p_val,
    Significance = ifelse(p_val < 0.05, "P < 0.05", "NS")
  )
  df_tot_Remnant <- rbind(df_tot_Remnant, tmp)
}

# formate df for plotting
df_tot_Remnant <- left_join(df_tot_Remnant, spxsp_df_spp, by = c("Species" = "tSpp"))
df_tot_Remnant$Species <- factor(df_tot_Remnant$Species, levels = df_tot_Remnant$Species[order(df_tot_Remnant$Estimate, decreasing = FALSE)])
df_tot_Remnant$tName <- factor(df_tot_Remnant$tName, levels = df_tot_Remnant$tName[order(df_tot_Remnant$Estimate, decreasing = FALSE)])
df_tot_Remnant$P_adj <- p.adjust(df_tot_Remnant$P_value, method = "BH")
df_tot_Remnant$Significance <- ifelse(df_tot_Remnant$P_value < 0.05, "P < 0.05", "NS")
df_tot_Remnant$Significance <- ifelse(df_tot_Remnant$P_adj < 0.05, "Padj < 0.05", df_tot_Remnant$Significance)
df_tot_Remnant$Significance <- factor(df_tot_Remnant$Significance, levels = c("Padj < 0.05", "P < 0.05", "NS"))

# ggplot
df_tot_Remnant$Common_species <- paste0(df_tot_Remnant$tName, " | ", df_tot_Remnant$Species)
df_tot_Remnant$Common_species <- factor(df_tot_Remnant$Common_species, levels = df_tot_Remnant$Common_species[order(df_tot_Remnant$Estimate, decreasing = FALSE)])

birdchange_rem <- ggplot(data = df_tot_Remnant, aes(x = tName, y = Estimate, fill = Significance))+
  geom_col()+coord_flip()+
  ggtitle("Remnant only")+
  scale_fill_manual(values = c("Padj < 0.05" = "darkblue", "P < 0.05" = "blue", "NS" = "lightgrey"))+
  theme_test()+
  labs(y = "Log odds ratio of species presence:\n 2015 (reference group) vs 2025", x = "Bird species")

birdchange_rem_sig <- ggplot(data = df_tot_Remnant[df_tot_Remnant$Significance == "P < 0.05" | df_tot_Remnant$Significance == "Padj < 0.05",], aes(x = tName, y = Estimate, fill = Significance))+
  geom_col()+coord_flip()+
  scale_fill_manual(values = c("Padj < 0.05" = "darkblue", "P < 0.05" = "blue", "NS" = "lightgrey"))+
  ggtitle("Remnant only")+
  theme_test()+
  labs(y = "Log odds ratio of species presence:\n 2015 (reference group) vs 2025", x = "Bird species")

table(df_tot_Remnant$Significance)
# Padj < 0.05    P < 0.05          NS 
#           2           7          20 

# arrange all plots together 
ggpubr::ggarrange(birdchange_2015, birdchange_2025,
                  align="hv", labels = c("a)", "b)"))
# ggsave(filename = "Log-odds-bird-RemRev.pdf", path = outdir, width = 12, height = 10)
ggpubr::ggarrange(birdchange_2015_sig, birdchange_2025_sig,
                  align="hv", labels = c("a)", "b)"))
# ggsave(filename = "Log-odds-bird-RemRev-significant_only.pdf", path = outdir, width = 12, height =10)

ggpubr::ggarrange(birdchange_rev, birdchange_rem,
                  align="hv", labels = c("a)", "b)"))
# ggsave(filename = "Log-odds-bird-survey_years.pdf", path = outdir, width = 12, height = 10)

ggpubr::ggarrange(birdchange_rev_sig, birdchange_rem_sig,
                  align="hv", labels = c("a)", "b)"))
# ggsave(filename = "Log-odds-bird-survey_years-significant_only.pdf", path = outdir, width = 12, height =10)

# EVONET functional information (global database of species) -------------------
all_bird_data2$Species1 <- all_bird_data2$tSpp

# Key terms and definitions ----------------------------------------------------
# TROPIC LEVEL:                                                                 
#    Herbivore = species obtaining at least 70% of food resources from plants;   
#    Carnivore = species obtaining at least 70% of food resources by consuming live invertebrate or vertebrate animals; 
#    Scavenger = species obtaining at least 70% of food resources from carrion or refuse; 
#    Omnivore = species obtaining resources from multiple trophic level in roughly equal proportion
# TROPIC NICHE:
#    Frugivore = species obtaining at least 60% of food resources from fruit;   
#    Granivore = species obtaining at least 60% of food resources from seeds or nuts; 
#    Nectarivore =  species obtaining at least 60% of food resources from nectar; 
#    Herbivore = species obtaining at least 60% of food resources from other plant materials in non-aquatic systems, including leaves, buds, whole flowers etc.; 
#    Herbivore aquatic = species obtaining at least 60% of food resources from plant materials in aquatic systems, including algae and aquatic plant leaves; 
#    Invertivore = species obtaining at least 60% of food resources from invertebrates in terrestrial systems, including insects, worms, arachnids, etc.; 
#    Vertivore = species obtaining at least 60% of food resources from vertebrate animals in terrestrial systems, including mammals, birds, reptiles etc.; 
#    Aquatic Predator = species obtaining at least 60% of food resources from vertebrate and invertebrate animals in aquatic systems, including fish, crustacea, molluscs, etc; 
#    Scavenger = species obtaining at least 60% of food resources from carrion, offal or refuse; 
#    Omnivore = Species using multiple niches, within or across trophic levels, in relatively equal proportions 
# PRIMARY LIFESTYLE                                                             
#    Aerial = species spends much of the time in flight, and hunts or forages predominantly on the wing; 
#    Terrestrial = species spends majority of its time on the ground, where it obtains food while either walking or hopping (note this includes species that also wade in water with their body raised above the water); 
#    Insessorial = species spends much of the time perching above the ground, either in branches of trees and other vegetation (i.e. arboreal), or on other raised substrates including rocks, buildings, posts, and wires; 
#    Aquatic = species spends much of the time sitting on water, and obtains food while afloat or when diving under the water's surface; 
#    Generalist = species has no primary lifestyle because it spends time in different lifestyle classes

## Attatch Functions to data ---------------------------------------------------
bird_functional_full 
bird_functional <- read.csv(file = "//PATH/TO/DATA/Bird database stuff/AVONET/ELEData/TraitData/AVONET1_BirdLife.csv")
bird_functional_full <- bird_functional

bird_functional <- bird_functional %>% 
  select(Species1, Family1, Order1, Mass, Habitat, Habitat.Density, Trophic.Level, Trophic.Niche, Primary.Lifestyle)
head(bird_functional$Species1)

bird_functional_CLLMM <- all_bird_data2 %>%
  left_join(bird_functional, by = "Species1")

# sanity check table
bird_functional_CLLMM_SC <- bird_functional_CLLMM %>%
  select(tSpp, tName, Species1, Family1, Order1, Mass, Habitat, Habitat.Density, Trophic.Level, Trophic.Niche, Primary.Lifestyle) %>%
  distinct()
table(bird_functional_CLLMM_SC$Trophic.Niche)
table(bird_functional_CLLMM_SC$Primary.Lifestyle)

length(unique(bird_functional_CLLMM_SC$Trophic.Niche))
length(unique(bird_functional_CLLMM_SC$Primary.Lifestyle))

# Couple of NA species
# Extract these and assign them without subspecies data
bird_functional_CLLMM_SC_NA <- subset(bird_functional_CLLMM_SC, is.na(Primary.Lifestyle) == TRUE)
colnames(bird_functional_CLLMM_SC_NA)

bird_functional_CLLMM_SC_NA$Species1 <- str_extract(bird_functional_CLLMM_SC_NA$tSpp, "^[^ ]+ [^ ]+")
bird_functional_CLLMM_SC_NA <- bird_functional_CLLMM_SC_NA %>% 
  select(tSpp, tName, Species1)
head(bird_functional_CLLMM_SC_NA)
head(bird_functional)
bird_functional_CLLMM_SC_NA1 <- left_join(bird_functional_CLLMM_SC_NA, bird_functional, by ="Species1")
# head(bird_functional_CLLMM_SC_NA1)

bird_functional_CLLMM_SC_NA1
# Atlas of Living Australia names and synonyms

name_map <- c(
  # Change to accepted names
  "Eolophus roseicapillus"       = "Eolophus roseicapilla",
  "Cracticus tibicen"            = "Gymnorhina tibicen",
  "Lichenostomus virescens"      = "Gavicalis virescens",
  "Threskiornis molucca"         = "Threskiornis moluccus",
  "Lichenostomus penicillatus"   = "Ptilotula penicillata",
  "Chroicocephalus novaehollandiae" = "Larus novaehollandiae",
  "Calyptorhynchus funereus"     = "Zanda funerea",
  "Stigmatopelia chinensis"      = "Spilopelia chinensis",
  "Lichenostomus chrysops"       = "Caligavis chrysops",
  "Lichenostomus ornatus"        = "Ptilotula ornata",
  "Sugomel niger"                = "Sugomel nigrum",
  # and move subspecies up
  "Platycercus elegans 'adelaidae' (NC)" = "Platycercus elegans",
  "Northiella haematogaster haematogaster (NC)" = "Northiella haematogaster"
)

all_bird_data_updated <- all_bird_data2

all_bird_data_updated$Species1 <-
  ifelse(
    all_bird_data_updated$Species1 %in% names(name_map),
    name_map[all_bird_data_updated$Species1],
    all_bird_data_updated$Species1
  )

# Scale up to larger dataframe
bird_functional_CLLMM_Fixed <- all_bird_data_updated %>%
  left_join(bird_functional, by = "Species1")

table(bird_functional_CLLMM_Fixed$Trophic.Niche)
table(bird_functional_CLLMM_Fixed$Habitat)
table(bird_functional_CLLMM_Fixed$Primary.Lifestyle)

bird_functional_CLLMM_Fixed$Tropic_Niche_Lifestyle <- paste0(bird_functional_CLLMM_Fixed$Primary.Lifestyle, ": ", bird_functional_CLLMM_Fixed$Trophic.Level, "-", bird_functional_CLLMM_Fixed$Trophic.Niche)
table(bird_functional_CLLMM_Fixed$Tropic_Niche_Lifestyle)
length(table(bird_functional_CLLMM_Fixed$Tropic_Niche_Lifestyle))

# Bird functional diversity
# wide_bird_functions

bird_functional_CLLMM_Fixed$One_ID_Var <- paste0(bird_functional_CLLMM_Fixed$Survey_Year, "_G", bird_functional_CLLMM_Fixed$WptID, "_", bird_functional_CLLMM_Fixed$iVisitID, "_py", bird_functional_CLLMM_Fixed$iPlantYear)
bird_functional_CLLMM_Fixed$Site_ID_Var <- paste0(bird_functional_CLLMM_Fixed$Survey_Year, "_G", bird_functional_CLLMM_Fixed$WptID, "_py", bird_functional_CLLMM_Fixed$iPlantYear)

wide_bird_functions <- bird_functional_CLLMM_Fixed %>%
  select(One_ID_Var, Tropic_Niche_Lifestyle) %>%
  mutate(present = 1) %>%
  
  group_by(One_ID_Var, Tropic_Niche_Lifestyle) %>%
  reframe(
    FunG_count = sum(present)
  ) %>%
  # distinct() %>%
  pivot_wider(
    id_cols = One_ID_Var,             
    names_from = Tropic_Niche_Lifestyle,
    values_from = FunG_count,
    values_fill = 0) %>% 
  as.data.frame()

rownames(wide_bird_functions) <- wide_bird_functions$One_ID_Var
wide_bird_Funs <- wide_bird_functions %>%
  select(-One_ID_Var)
wide_bird_Funs[1:5,1:5]

library(vegan)
bird_fun_diversity_Richness <- specnumber(wide_bird_Funs)
bird_fun_diversity_Shannon <- diversity(wide_bird_Funs, index = "shannon")
bird_fun_diversity_InvSimpson <- diversity(wide_bird_Funs, index = "invsimpson")

One_ID_Var_frame <- bird_functional_CLLMM_Fixed %>%
  select(One_ID_Var, Survey_Year, WptID, iVisitID, Treat_type, iPlantYear) %>%
  distinct()

bird_fun_alphadiversity <- data.frame(
  One_ID_Var = names(bird_fun_diversity_Richness),
  Richness = bird_fun_diversity_Richness,
  Shannon = bird_fun_diversity_Shannon,
  InvSimpson = bird_fun_diversity_InvSimpson
)
bird_fun_alphadiversity$Effective_Shannons <- exp(bird_fun_alphadiversity$Shannon)

bird_fun_alphadiversity <- left_join(One_ID_Var_frame, bird_fun_alphadiversity, by = "One_ID_Var")
colnames(bird_fun_alphadiversity)
bird_fun_alphadiversity$Site_ID_Var <- paste0(bird_fun_alphadiversity$Survey_Year, "_G", bird_fun_alphadiversity$WptID, "_py", bird_fun_alphadiversity$iPlantYear)


mean_bird_fun_alphadiversity <- bird_fun_alphadiversity %>%
  group_by(Site_ID_Var, Survey_Year, Treat_type) %>%
  reframe( mean_Richness = mean(Richness), 
           mean_Shannon = mean(Shannon),
           mean_InvSimpson = mean(InvSimpson),
           mean_Effective_Shannons = mean(Effective_Shannons))

mean_bird_fun_alphadiversity$Treat_type <- factor(mean_bird_fun_alphadiversity$Treat_type, levels = c("Revegetated", "Remnant"))

library(ggpubr)

plot_fun_richness <- ggplot(mean_bird_fun_alphadiversity, aes(x = Treat_type, y = mean_Richness, fill = Treat_type))+
  geom_violin()+
  geom_boxplot(outlier.shape = NA, width = 0.1, fill = "white")+
  geom_jitter(height = 0, width = 0.15)+
  theme_test()+
  facet_grid(~Survey_Year)+
  xlab("")
plot_fun_shannons <- ggplot(mean_bird_fun_alphadiversity, aes(x = Treat_type, y = mean_Shannon, fill = Treat_type))+
  geom_violin()+
  geom_boxplot(outlier.shape = NA, width = 0.1, fill = "white")+
  geom_jitter(height = 0, width = 0.15)+
  theme_test()+
  facet_grid(~Survey_Year)+
  xlab("")
plot_fun_InvSimpson <- ggplot(mean_bird_fun_alphadiversity, aes(x = Treat_type, y = mean_InvSimpson, fill = Treat_type))+
  geom_violin()+
  geom_boxplot(outlier.shape = NA, width = 0.1, fill = "white")+
  geom_jitter(height = 0, width = 0.15)+
  theme_test()+
  facet_grid(~Survey_Year)+
  xlab("")
plot_fun_ENF <- ggplot(mean_bird_fun_alphadiversity, aes(x = Treat_type, y = mean_Effective_Shannons, fill = Treat_type))+
  geom_violin()+
  geom_boxplot(outlier.shape = NA, width = 0.1, fill = "white")+
  geom_jitter(height = 0, width = 0.15)+
  theme_test()+
  facet_grid(~Survey_Year)+
  xlab("")

ggarrange(
  plot_fun_richness,
  plot_fun_shannons,
  plot_fun_InvSimpson,
  plot_fun_ENF, 
  common.legend = TRUE, align = "hv"
)

fun_mod_rich_alpha <- lmer(Richness ~ Treat_type * Survey_Year + (1|Site_ID_Var), data =bird_fun_alphadiversity)
car::Anova(fun_mod_rich_alpha)

fun_mod_ENF_alpha <- lmer(Effective_Shannons ~ Treat_type + Survey_Year + (1|Site_ID_Var), data =bird_fun_alphadiversity)
car::Anova(fun_mod_ENF_alpha)

### Function Beta ords -----------------------------------------------------
colnames(bird_functional_CLLMM_Fixed)
# "Survey_Year"            "WptID"                  "iVisitID"               "Treat_type"             "tNSXCode"              
# "tSpp"                   "tName"                  "iPlantYear"             "Species1"               "Family1"               
# "Order1"                 "Mass"                   "Habitat"                "Habitat.Density"        "Trophic.Level"         
# "Trophic.Niche"          "Primary.Lifestyle"      "Tropic_Niche_Lifestyle"

# bird_functional_CLLMM_Fixed

bird_functional_CLLMM_Fixed$One_ID_Var <- paste0(bird_functional_CLLMM_Fixed$Survey_Year, "_G", bird_functional_CLLMM_Fixed$WptID, "_", bird_functional_CLLMM_Fixed$iVisitID, "_py", bird_functional_CLLMM_Fixed$iPlantYear)

wide_bird_functions <- bird_functional_CLLMM_Fixed %>%
  select(One_ID_Var, Tropic_Niche_Lifestyle) %>%
  mutate(present = 1) %>%
  
  group_by(One_ID_Var, Tropic_Niche_Lifestyle) %>%
  reframe(
    FunG_count = sum(present)
  ) %>%
  # distinct() %>%
  pivot_wider(
    id_cols = One_ID_Var,             
    names_from = Tropic_Niche_Lifestyle,
    values_from = FunG_count,
    values_fill = 0) %>% 
  as.data.frame()

rownames(wide_bird_functions) <- wide_bird_functions$One_ID_Var

wide_bird_Funs <- wide_bird_functions %>%
  select(-One_ID_Var)
wide_bird_Funs[1:5,1:5]

# Site level only
all_bird_data_updated2 <- all_bird_data_updated
bird_functional_CLLMM_Fixed$Site_ID_Var <- paste0(all_bird_data_updated2$Survey_Year, "_G", all_bird_data_updated2$WptID, "_py", all_bird_data_updated2$iPlantYear)

wide_bird_functions_SITE <- bird_functional_CLLMM_Fixed %>%
  select(One_ID_Var, Tropic_Niche_Lifestyle) %>%
  mutate(present = 1) %>%
  
  group_by(One_ID_Var, Tropic_Niche_Lifestyle) %>%
  reframe(
    FunG_count = sum(present)
  ) %>%
  # distinct() %>%
  pivot_wider(
    id_cols = One_ID_Var,             
    names_from = Tropic_Niche_Lifestyle,
    values_from = FunG_count,
    values_fill = 0) %>% 
  as.data.frame()

wide_bird_functions_SITE$Site_ID_Var <- sub("_[0-9]+_", "_", wide_bird_functions_SITE$One_ID_Var)

wide_bird_functions_SITE2 <- wide_bird_functions_SITE %>%
  group_by(Site_ID_Var) %>%
  summarise(
    across(
      where(is.numeric),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  as.data.frame()

rownames(wide_bird_functions_SITE2) <- wide_bird_functions_SITE2$Site_ID_Var

wide_bird_functions_SITE2 <- wide_bird_functions_SITE2 %>%
  select(-Site_ID_Var)
wide_bird_functions_SITE2[1:5,1:5]

# Plot 
set.seed(123)
NMDS_Jacc_fun <- metaMDS(wide_bird_functions_SITE2, distance = "jaccard", binary = TRUE)
NMDS_Jacc_fun$stress # 0.2155163

plot(NMDS_Jacc_fun)

points_Jacc_funs <- NMDS_Jacc_fun$points

points_Jacc_funs <- as.data.frame(points_Jacc_funs)
points_Jacc_funs$Site_ID_Var <- rownames(points_Jacc_funs)

# attatch metadata back
# all_bird_site_beta <- all_bird_data %>% 
#   select(Site_ID_Var, Survey_Year, WptID, iPlantYear, Treat_type) %>%
#   distinct()
all_bird_site_beta

head(all_bird_site_beta)
head(points_Jacc_funs)
colnames(points_Jacc_funs) <- c("MDS1_FUN", "MDS2_FUN", "Site_ID_Var")

bird_metadata_s_Jacc2 <- left_join(bird_metadata_s_Jacc, points_Jacc_funs, by = "Site_ID_Var")
bird_metadata_s_Jacc2$Survey_Year <- as.factor(bird_metadata_s_Jacc2$Survey_Year)

colnames(bird_metadata_s_Jacc2)
# "Site_ID"      "Survey_Year"  "iWptID"       "iPlantYear"   "iTreatID"     "iEcosystemID"
# "tEcoName"     "MDS1"         "MDS2"
head(bird_metadata_s_Jacc2)

ggplot(bird_metadata_s_Jacc2, aes(x = MDS1_FUN, y = MDS2_FUN, shape = Survey_Year, colour =Treat_type))+
  geom_point(size= 2)+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()
# ggsave(filename = "Bird-BETA-NMDS-Jaccard-FUNCTIONAL-SurveyYear_RemRev.pdf", path = outdir, width = 7, height = 6)

ggplot(bird_metadata_s_Jacc, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =iPlantYear))+
  geom_point(size= 2)+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()
# ggsave(filename = "Bird-BETA-NMDS-Jaccard-FUNCTIONAL-SurveyYear_PlantYear.pdf", path = outdir, width = 7, height = 6)


# Proportion of functional groups for identified species 

## Tropic Niche ----------------------------------------------------------------
bird_functional_CLLMM_Fixed_prop <- bird_functional_CLLMM_Fixed %>%
  select(One_ID_Var, Trophic.Niche, Survey_Year, Treat_type) %>%
  mutate(present = 1) %>%
  group_by(One_ID_Var, Trophic.Niche, Survey_Year, Treat_type) %>%
  summarise(
    FunG_count = sum(present),
    .groups = "drop_last"
  ) %>% 
  ungroup()

# dim(bird_functional_CLLMM_Fixed_prop)
# # bird_functional_CLLMM_Fixed_prop <- na.omit(bird_functional_CLLMM_Fixed_prop)
# dim(bird_functional_CLLMM_Fixed_prop)

bird_functional_CLLMM_Fixed_prop$Site_ID_Var <- sub("_[0-9]+_", "_", bird_functional_CLLMM_Fixed_prop$One_ID_Var)
head(bird_functional_CLLMM_Fixed_prop)

bird_functional_CLLMM_Fixed_prop2 <- bird_functional_CLLMM_Fixed_prop %>%
  group_by(Site_ID_Var, Trophic.Niche, Survey_Year, Treat_type) %>%
  reframe(
    mean_FunG_count = mean(FunG_count)
  ) %>% 
  ungroup() %>%
  group_by(Survey_Year, Treat_type, Trophic.Niche) %>%
  reframe(
    mean_count = mean(mean_FunG_count)
  ) %>% 
  ungroup() %>%
  group_by(Survey_Year, Treat_type) %>%   
  mutate(
    FunG_Proportion = mean_count / sum(mean_count)
  ) %>%
  ungroup()
head(bird_functional_CLLMM_Fixed_prop2)

# Sanity check ;) 
bird_functional_CLLMM_Fixed_prop2 %>%
  group_by(Survey_Year, Treat_type) %>%
  summarise(sum_prop = sum(FunG_Proportion))

# Recombine with metadata
meta_b_F <- bird_functional_CLLMM_Fixed %>%
  select(Survey_Year, Treat_type) %>%
  distinct()

bird_fun_Groups_stack <- left_join(meta_b_F, bird_functional_CLLMM_Fixed_prop2, by = c("Survey_Year", "Treat_type"))

head(bird_fun_Groups_stack)

Niche_colours <- c(
  "Aquatic predator"  = "#6A1B9A", 
  "Granivore"         = "#1B5E20", 
  "Herbivore aquatic" = "#B9F6CA", 
  "Invertivore"       = "#EF9A9A", 
  "Nectarivore"       = "#FFB300", 
  "Omnivore"          = "#64B5F6",
  "Vertivore"         = "#BF360C", 
  "NA: NA"            = "darkgrey"
)

bird_fun_Groups_stack$Treat_type <- factor(bird_fun_Groups_stack$Treat_type, levels = c("Revegetated", "Remnant"))

ggpubr::ggarrange(
  ggplot(bird_fun_Groups_stack, aes(x = Treat_type, y = FunG_Proportion, fill = Trophic.Niche)) +
    geom_col(colour = "black") +
    labs(x = "Treatment", y = "Proportion", fill = "Trophic niche") +
    theme_bw() +
    scale_fill_manual(values = Niche_colours)+
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ,
  ggplot(na.omit(bird_fun_Groups_stack), aes(x = Treat_type, y = Trophic.Niche, fill = FunG_Proportion)) +
    # geom_col(colour = "black") +
    geom_tile(colour = "black")+
    labs(x = "Treatment", y = "Trophic niche", fill = "Proportion") +
    theme_bw() +
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  align = "hv"
  )
# ggsave(filename = "Bird-Stack-FUNCTIONS-TrophicNiche-SurveyYear_RemRev.pdf", path = outdir, width = 13, height = 6)

ggpubr::ggarrange(
  ggplot(bird_fun_Groups_stack, aes(x = Treat_type, y = mean_count, fill = Trophic.Niche)) +
    geom_col(colour = "black") +
    labs(x = "Treatment", y = "Mean count", fill = "Trophic niche") +
    theme_bw() +
    scale_fill_manual(values = Niche_colours)+
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ,
  ggplot(na.omit(bird_fun_Groups_stack), aes(x = Treat_type, y = Trophic.Niche, fill = mean_count)) +
    # geom_col(colour = "black") +
    geom_tile(colour = "black")+
    labs(x = "Treatment", y = "Trophic niche", fill = "Mean count") +
    theme_bw() +
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  align = "hv"
)
# ggsave(filename = "Bird-Stack-FUNCTIONS-COUNT-TrophicNiche-SurveyYear_RemRev.pdf", path = outdir, width = 13, height = 6)

## Lifestyle -------------------------------------------------------------------
bird_functional_CLLMM_Fixed_prop_LS <- bird_functional_CLLMM_Fixed %>%
  select(One_ID_Var, Primary.Lifestyle, Survey_Year, Treat_type) %>%
  mutate(present = 1) %>%
  group_by(One_ID_Var, Primary.Lifestyle, Survey_Year, Treat_type) %>%
  summarise(
    FunG_count = sum(present),
    .groups = "drop_last"
  ) %>% 
  ungroup()

# dim(bird_functional_CLLMM_Fixed_prop)
# # bird_functional_CLLMM_Fixed_prop <- na.omit(bird_functional_CLLMM_Fixed_prop)
# dim(bird_functional_CLLMM_Fixed_prop)

bird_functional_CLLMM_Fixed_prop_LS$Site_ID_Var <- sub("_[0-9]+_", "_", bird_functional_CLLMM_Fixed_prop_LS$One_ID_Var)
head(bird_functional_CLLMM_Fixed_prop_LS)

bird_functional_CLLMM_Fixed_prop_LS2 <- bird_functional_CLLMM_Fixed_prop_LS %>%
  group_by(Site_ID_Var, Primary.Lifestyle, Survey_Year, Treat_type) %>%
  reframe(
    mean_FunG_count = mean(FunG_count)
  ) %>% 
  ungroup() %>%
  group_by(Survey_Year, Treat_type, Primary.Lifestyle) %>%
  reframe(
    mean_count = mean(mean_FunG_count)
  ) %>% 
  ungroup() %>%
  group_by(Survey_Year, Treat_type) %>%   
  mutate(
    FunG_Proportion = mean_count / sum(mean_count)
  ) %>%
  ungroup()
head(bird_functional_CLLMM_Fixed_prop_LS2)

# Sanity check ;) 
bird_functional_CLLMM_Fixed_prop_LS2 %>%
  group_by(Survey_Year, Treat_type) %>%
  summarise(sum_prop = sum(FunG_Proportion))

# Recombine with metadata
meta_b_F <- bird_functional_CLLMM_Fixed %>%
  select(Survey_Year, Treat_type) %>%
  distinct()

bird_fun_Groups_stack_Lifestyle <- left_join(meta_b_F, bird_functional_CLLMM_Fixed_prop_LS2, by = c("Survey_Year", "Treat_type"))

head(bird_fun_Groups_stack_Lifestyle)

Lifestyle_colours <- c(
  "Aerial"          = "#64B5F6",
  "Aquatic"         = "#B9F6CA", 
  "Generalist"      = "#6A1B9A", 
  "Insessorial"  = "#2E7D32", 
  "Terrestrial"     = "#FFB300", 
  "NA: NA"          = "darkgrey"
)

bird_fun_Groups_stack_Lifestyle$Treat_type <- factor(bird_fun_Groups_stack_Lifestyle$Treat_type, levels = c("Revegetated", "Remnant"))

ggpubr::ggarrange(
  ggplot(bird_fun_Groups_stack_Lifestyle, aes(x = Treat_type, y = FunG_Proportion, fill = Primary.Lifestyle)) +
    geom_col(colour = "black") +
    labs(x = "Treatment", y = "Proportion", fill = "Primary lifestyle") +
    theme_bw() +
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values = Lifestyle_colours)
  ,
  ggplot(na.omit(bird_fun_Groups_stack_Lifestyle), aes(x = Treat_type, y = Primary.Lifestyle, fill = FunG_Proportion)) +
    # geom_col(colour = "black") +
    geom_tile(colour = "black")+
    labs(x = "Treatment", y = "Primary lifestyle", fill = "Proportion") +
    theme_bw() +
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  align = "hv"
)
# ggsave(filename = "Bird-Stack-FUNCTIONS-PrimaryLifestyle-SurveyYear_RemRev.pdf", path = outdir, width = 13, height = 6)

ggpubr::ggarrange(
  ggplot(bird_fun_Groups_stack_Lifestyle, aes(x = Treat_type, y = mean_count, fill = Primary.Lifestyle)) +
    geom_col(colour = "black") +
    labs(x = "Treatment", y = "Mean count", fill = "Primary lifestyle") +
    theme_bw() +
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_fill_manual(values = Lifestyle_colours)
  ,
  ggplot(na.omit(bird_fun_Groups_stack_Lifestyle), aes(x = Treat_type, y = Primary.Lifestyle, fill = mean_count)) +
    # geom_col(colour = "black") +
    geom_tile(colour = "black")+
    labs(x = "Treatment", y = "Primary lifestyle", fill = "Mean count") +
    theme_bw() +
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  align = "hv"
)
# ggsave(filename = "Bird-Stack-FUNCTIONS-COUNT-PrimaryLifestyle-SurveyYear_RemRev.pdf", path = outdir, width = 13, height = 6)

## Trophic Level: Niche --------------------------------------------------------
bird_functional_CLLMM_Fixed_prop_TLN <- bird_functional_CLLMM_Fixed
bird_functional_CLLMM_Fixed_prop_TLN$Trophic.Level_Niche <- paste0(bird_functional_CLLMM_Fixed_prop_TLN$Trophic.Level, ": ", bird_functional_CLLMM_Fixed_prop_TLN$Trophic.Niche)

bird_functional_CLLMM_Fixed_prop_TLN <- bird_functional_CLLMM_Fixed_prop_TLN %>%
  select(One_ID_Var, Trophic.Level_Niche, Survey_Year, Treat_type) %>%
  mutate(present = 1) %>%
  group_by(One_ID_Var, Trophic.Level_Niche, Survey_Year, Treat_type) %>%
  summarise(
    FunG_count = sum(present),
    .groups = "drop_last"
  ) %>% 
  ungroup()

# dim(bird_functional_CLLMM_Fixed_prop)
# # bird_functional_CLLMM_Fixed_prop <- na.omit(bird_functional_CLLMM_Fixed_prop)
# dim(bird_functional_CLLMM_Fixed_prop)

bird_functional_CLLMM_Fixed_prop_TLN$Site_ID_Var <- sub("_[0-9]+_", "_", bird_functional_CLLMM_Fixed_prop_TLN$One_ID_Var)
head(bird_functional_CLLMM_Fixed_prop_TLN)

bird_functional_CLLMM_Fixed_prop_TLN2 <- bird_functional_CLLMM_Fixed_prop_TLN %>%
  group_by(Site_ID_Var, Trophic.Level_Niche, Survey_Year, Treat_type) %>%
  reframe(
    mean_FunG_count = mean(FunG_count)
  ) %>% 
  ungroup() %>%
  group_by(Survey_Year, Treat_type, Trophic.Level_Niche) %>%
  reframe(
    mean_count = mean(mean_FunG_count)
  ) %>% 
  ungroup() %>%
  group_by(Survey_Year, Treat_type) %>%   
  mutate(
    FunG_Proportion = mean_count / sum(mean_count)
  ) %>%
  ungroup()
head(bird_functional_CLLMM_Fixed_prop_TLN2)

# Sanity check ;) 
bird_functional_CLLMM_Fixed_prop_TLN2 %>%
  group_by(Survey_Year, Treat_type) %>%
  summarise(sum_prop = sum(FunG_Proportion))

# Recombine with metadata
meta_b_F <- bird_functional_CLLMM_Fixed %>%
  select(Survey_Year, Treat_type) %>%
  distinct()

bird_fun_Groups_stack_TLN <- left_join(meta_b_F, bird_functional_CLLMM_Fixed_prop_TLN2, by = c("Survey_Year", "Treat_type"))

head(bird_fun_Groups_stack_TLN)

unique(bird_fun_Groups_stack_TLN$Trophic.Level_Niche)
TLN_colours <- c(
  "Carnivore: Aquatic predator" = "#BF360C", "Herbivore: Granivore"         = "#1B5E20", "Omnivore: Invertivore" = "#0D47A1",
  "Carnivore: Invertivore"      = "#EF6C00", "Herbivore: Herbivore aquatic" = "#2E7D32", "Omnivore: Nectarivore" = "#1976D2",
  "Carnivore: Omnivore"         = "#FFB300", "Herbivore: Nectarivore"       = "#66BB6A", "Omnivore: Omnivore"    = "#64B5F6",
  "Carnivore: Vertivore"        = "#FFE0B2",  "Herbivore: Omnivore"         = "#B9F6CA", "NA: NA"                = "darkgrey"
  )

bird_fun_Groups_stack_TLN$Treat_type <- factor(bird_fun_Groups_stack_TLN$Treat_type, levels = c("Revegetated", "Remnant"))

ggpubr::ggarrange(
  ggplot(bird_fun_Groups_stack_TLN, aes(x = Treat_type, y = FunG_Proportion, fill = Trophic.Level_Niche)) +
    geom_col(colour = "black") +
    labs(x = "Treatment", y = "Proportion", fill = "Tropic Levels (Niche)") +
    theme_bw() +
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    scale_fill_manual(values = TLN_colours)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ,
  ggplot(bird_fun_Groups_stack_TLN, aes(x = Treat_type, y = Trophic.Level_Niche, fill = FunG_Proportion)) +
    # geom_col(colour = "black") +
    geom_tile(colour = "black")+
    labs(x = "Treatment", y = "Tropic Levels (Niche)", fill = "Proportion") +
    theme_bw() +
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  align = "hv"
)
# ggsave(filename = "Bird-Stack-FUNCTIONS-CombTrophicLevelNiche-SurveyYear_RemRev.pdf", path = outdir, width = 13, height = 6)

ggpubr::ggarrange(
  ggplot(bird_fun_Groups_stack_TLN, aes(x = Treat_type, y = mean_count, fill = Trophic.Level_Niche)) +
    geom_col(colour = "black") +
    labs(x = "Treatment", y = "Mean Count", fill = "Tropic Levels (Niche)") +
    theme_bw() +
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    scale_fill_manual(values = TLN_colours)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ,
  ggplot(bird_fun_Groups_stack_TLN, aes(x = Treat_type, y = Trophic.Level_Niche, fill = mean_count)) +
    # geom_col(colour = "black") +
    geom_tile(colour = "black")+
    labs(x = "Treatment", y = "Tropic Levels (Niche)", fill = "Mean Count") +
    theme_bw() +
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  align = "hv"
)
# ggsave(filename = "Bird-Stack-FUNCTIONS-COUNT-CombTrophicLevelNiche-SurveyYear_RemRev.pdf", path = outdir, width = 13, height = 6)

## Lifestyle: Niche --------------------------------------------------------
bird_functional_CLLMM_Fixed_prop_LifeNiche <- bird_functional_CLLMM_Fixed
bird_functional_CLLMM_Fixed_prop_LifeNiche$Trophic.Lifestyle_Niche <- paste0(bird_functional_CLLMM_Fixed_prop_LifeNiche$Primary.Lifestyle, ": ", bird_functional_CLLMM_Fixed_prop_LifeNiche$Trophic.Niche)

bird_functional_CLLMM_Fixed_prop_LifeNiche <- bird_functional_CLLMM_Fixed_prop_LifeNiche %>%
  select(One_ID_Var, Trophic.Lifestyle_Niche, Survey_Year, Treat_type) %>%
  mutate(present = 1) %>%
  group_by(One_ID_Var, Trophic.Lifestyle_Niche, Survey_Year, Treat_type) %>%
  summarise(
    FunG_count = sum(present),
    .groups = "drop_last"
  ) %>% 
  ungroup()

# dim(bird_functional_CLLMM_Fixed_prop)
# # bird_functional_CLLMM_Fixed_prop <- na.omit(bird_functional_CLLMM_Fixed_prop)
# dim(bird_functional_CLLMM_Fixed_prop)

bird_functional_CLLMM_Fixed_prop_LifeNiche$Site_ID_Var <- sub("_[0-9]+_", "_", bird_functional_CLLMM_Fixed_prop_LifeNiche$One_ID_Var)
head(bird_functional_CLLMM_Fixed_prop_LifeNiche)

bird_functional_CLLMM_Fixed_prop_LifeNiche2 <- bird_functional_CLLMM_Fixed_prop_LifeNiche %>%
  group_by(Site_ID_Var, Trophic.Lifestyle_Niche, Survey_Year, Treat_type) %>%
  reframe(
    mean_FunG_count = mean(FunG_count)
  ) %>% 
  ungroup() %>%
  group_by(Survey_Year, Treat_type, Trophic.Lifestyle_Niche) %>%
  reframe(
    mean_count = mean(mean_FunG_count)
  ) %>% 
  ungroup() %>%
  group_by(Survey_Year, Treat_type) %>%   
  mutate(
    FunG_Proportion = mean_count / sum(mean_count)
  ) %>%
  ungroup()
head(bird_functional_CLLMM_Fixed_prop_LifeNiche2)

# Sanity check ;) 
bird_functional_CLLMM_Fixed_prop_LifeNiche2 %>%
  group_by(Survey_Year, Treat_type) %>%
  summarise(sum_prop = sum(FunG_Proportion))

# Recombine with metadata
meta_b_F <- bird_functional_CLLMM_Fixed %>%
  select(Survey_Year, Treat_type) %>%
  distinct()

bird_fun_Groups_stack_LifeNiche <- left_join(meta_b_F, bird_functional_CLLMM_Fixed_prop_LifeNiche2, by = c("Survey_Year", "Treat_type"))

head(bird_fun_Groups_stack_LifeNiche)

unique(bird_fun_Groups_stack_LifeNiche$Trophic.Lifestyle_Niche)
Life_niche_colours <- c(
  "Insessorial: Granivore"   = "#4A148C", "Aerial: Aquatic predator" = "#E3F2FD", "Terrestrial: Aquatic predator" = "#BF360C", "Generalist: Invertivore" = "#8E0000", "Aquatic: Aquatic predator"  = "#00695C",  
  "Insessorial: Invertivore" = "#6A1B9A", "Aerial: Vertivore"        = "#0D47A1", "Terrestrial: Granivore"        = "#EF6C00", "Generalist: Omnivore"    = "#C62828", "Aquatic: Herbivore aquatic" = "#80CBC4",  
  "Insessorial: Nectarivore" = "#8E24AA", "Aerial: Invertivore"      = "#1976D2", "Terrestrial: Invertivore"      = "#FFB300", "Generalist: Vertivore"   = "#EF9A9A",
  "Insessorial: Omnivore"    = "#CE93D8", "Aerial: Omnivore"         = "#64B5F6", "Terrestrial: Omnivore"         = "#FFE0B2",  
  "Insessorial: Vertivore"   = "#F3E5F5", "NA: NA" = "darkgrey"
  )

bird_fun_Groups_stack_LifeNiche$Treat_type <- factor(bird_fun_Groups_stack_LifeNiche$Treat_type, levels = c("Revegetated", "Remnant"))

ggpubr::ggarrange(
  ggplot(bird_fun_Groups_stack_LifeNiche, aes(x = Treat_type, y = FunG_Proportion, fill = Trophic.Lifestyle_Niche)) +
    geom_col(colour = "black") +
    labs(x = "Treatment", y = "Proportion", fill = "Lifestyle: Tropic Niche") +
    theme_test() +
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    scale_fill_manual(values = Life_niche_colours)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ,
  ggplot(bird_fun_Groups_stack_LifeNiche, aes(x = Treat_type, y = Trophic.Lifestyle_Niche, fill = FunG_Proportion)) +
    # geom_col(colour = "black") +
    geom_tile(colour = "black")+
    labs(x = "Treatment", y = "Lifestyle: Tropic Niche", fill = "Proportion") +
    theme_test() +
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  align = "hv"
)
# ggsave(filename = "Bird-Stack-FUNCTIONS-CombLifestyleNiche-SurveyYear_RemRev.pdf", path = outdir, width = 13, height = 6)

ggpubr::ggarrange(
  ggplot(bird_fun_Groups_stack_LifeNiche, aes(x = Treat_type, y = mean_count, fill = Trophic.Lifestyle_Niche)) +
    geom_col(colour = "black") +
    labs(x = "Treatment", y = "Mean count", fill = "Lifestyle: Tropic Niche") +
    theme_test() +
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    scale_fill_manual(values = Life_niche_colours)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ,
  ggplot(bird_fun_Groups_stack_LifeNiche, aes(x = Treat_type, y = Trophic.Lifestyle_Niche, fill = mean_count)) +
    # geom_col(colour = "black") +
    geom_tile(colour = "black")+
    labs(x = "Treatment", y = "Lifestyle: Tropic Niche", fill = "Mean count") +
    theme_test() +
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  align = "hv"
)
# ggsave(filename = "Bird-Stack-FUNCTIONS-COUNT-CombLifestyleNiche-SurveyYear_RemRev.pdf", path = outdir, width = 13, height = 6)

## Niche: Lifestyle --------------------------------------------------------
bird_functional_CLLMM_Fixed_prop_LifeNicheREV <- bird_functional_CLLMM_Fixed
bird_functional_CLLMM_Fixed_prop_LifeNicheREV$Trophic.Lifestyle_NicheREV <- paste0(bird_functional_CLLMM_Fixed_prop_LifeNicheREV$Trophic.Niche, ": ",  bird_functional_CLLMM_Fixed_prop_LifeNicheREV$Primary.Lifestyle)

bird_functional_CLLMM_Fixed_prop_LifeNicheREV <- bird_functional_CLLMM_Fixed_prop_LifeNicheREV %>%
  select(One_ID_Var, Trophic.Lifestyle_NicheREV, Survey_Year, Treat_type) %>%
  mutate(present = 1) %>%
  group_by(One_ID_Var, Trophic.Lifestyle_NicheREV, Survey_Year, Treat_type) %>%
  summarise(
    FunG_count = sum(present),
    .groups = "drop_last"
  ) %>% 
  ungroup()

# dim(bird_functional_CLLMM_Fixed_prop)
# # bird_functional_CLLMM_Fixed_prop <- na.omit(bird_functional_CLLMM_Fixed_prop)
# dim(bird_functional_CLLMM_Fixed_prop)

bird_functional_CLLMM_Fixed_prop_LifeNicheREV$Site_ID_Var <- sub("_[0-9]+_", "_", bird_functional_CLLMM_Fixed_prop_LifeNicheREV$One_ID_Var)
head(bird_functional_CLLMM_Fixed_prop_LifeNicheREV)

bird_functional_CLLMM_Fixed_prop_LifeNicheREV2 <- bird_functional_CLLMM_Fixed_prop_LifeNicheREV %>%
  group_by(Site_ID_Var, Trophic.Lifestyle_NicheREV, Survey_Year, Treat_type) %>%
  reframe(
    mean_FunG_count = mean(FunG_count)
  ) %>% 
  ungroup() %>%
  group_by(Survey_Year, Treat_type, Trophic.Lifestyle_NicheREV) %>%
  reframe(
    mean_count = mean(mean_FunG_count)
  ) %>% 
  ungroup() %>%
  group_by(Survey_Year, Treat_type) %>%   
  mutate(
    FunG_Proportion = mean_count / sum(mean_count)
  ) %>%
  ungroup()
head(bird_functional_CLLMM_Fixed_prop_LifeNicheREV2)

# Sanity check ;) 
bird_functional_CLLMM_Fixed_prop_LifeNicheREV2 %>%
  group_by(Survey_Year, Treat_type) %>%
  summarise(sum_prop = sum(FunG_Proportion))

# Recombine with metadata
meta_b_F <- bird_functional_CLLMM_Fixed %>%
  select(Survey_Year, Treat_type) %>%
  distinct()

bird_fun_Groups_stack_LifeNicheREV <- left_join(meta_b_F, bird_functional_CLLMM_Fixed_prop_LifeNicheREV2, by = c("Survey_Year", "Treat_type"))

head(bird_fun_Groups_stack_LifeNicheREV)

unique(bird_fun_Groups_stack_LifeNicheREV$Trophic.Lifestyle_NicheREV)

NICHE_LIFE_colours <- c(
  "Aquatic predator: Aerial"      = "#DAA520",
  "Aquatic predator: Aquatic"     = "#FFDB58",
   "Aquatic predator: Terrestrial"= "#FFB300",
   
   "Granivore: Insessorial"       = "#1B5E20",
   "Granivore: Terrestrial"       = "#66BB6A",
   
   "Herbivore aquatic: Aquatic"   = "#FFE0B2",
   
   "Invertivore: Aerial"          = "#3E2723",
   "Invertivore: Generalist"      = "#6D4C41",
   "Invertivore: Insessorial"     = "#A1887F",
   "Invertivore: Terrestrial"     = "#D7CCC8",
   
   "Nectarivore: Insessorial"     = "#FF6F61",
   
   "Omnivore: Aerial"             = "#E3F2FD",
   "Omnivore: Generalist"         = "#0D47A1",
   "Omnivore: Insessorial"        = "#1976D2",
   "Omnivore: Terrestrial"        = "#64B5F6",
   
   "Vertivore: Generalist"        = "#BF360C",
   "Vertivore: Insessorial"       = "#EF6C00",
   
   "Vertivore: Aerial"            = "#4B0082",
   
   "NA: NA"                       = "darkgrey"
   )

bird_fun_Groups_stack_LifeNicheREV$Treat_type <- factor(bird_fun_Groups_stack_LifeNicheREV$Treat_type, levels = c("Revegetated", "Remnant"))

ggpubr::ggarrange(
  ggplot(bird_fun_Groups_stack_LifeNicheREV, aes(x = Treat_type, y = FunG_Proportion, fill = Trophic.Lifestyle_NicheREV)) +
    geom_col(colour = "black") +
    labs(x = "Treatment", y = "Proportion", fill = "Tropic Niche: Lifestyle") +
    theme_test() +
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    scale_fill_manual(values = NICHE_LIFE_colours)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ,
  ggplot(bird_fun_Groups_stack_LifeNicheREV, aes(x = Treat_type, y = Trophic.Lifestyle_NicheREV, fill = FunG_Proportion)) +
    # geom_col(colour = "black") +
    geom_tile(colour = "black")+
    labs(x = "Treatment", y = "Tropic Niche: Lifestyle", fill = "Proportion") +
    theme_test() +
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  align = "hv"
)
# ggsave(filename = "Bird-Stack-FUNCTIONS-Comb-Niche-To-Lifestyle-SurveyYear_RemRev.pdf", path = outdir, width = 13, height = 6)

ggpubr::ggarrange(
  ggplot(bird_fun_Groups_stack_LifeNicheREV, aes(x = Treat_type, y = mean_count, fill = Trophic.Lifestyle_NicheREV)) +
    geom_col(colour = "black") +
    labs(x = "Treatment", y = "Mean count", fill = "Lifestyle: Tropic Niche") +
    theme_bw() +
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    scale_fill_manual(values = NICHE_LIFE_colours)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ,
  ggplot(bird_fun_Groups_stack_LifeNicheREV, aes(x = Treat_type, y = Trophic.Lifestyle_NicheREV, fill = mean_count)) +
    # geom_col(colour = "black") +
    geom_tile(colour = "black")+
    labs(x = "Treatment", y = "Lifestyle: Tropic Niche", fill = "Mean count") +
    theme_bw() +
    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
  align = "hv"
)
# ggsave(filename = "Bird-Stack-FUNCTIONS-COUNT-Comb-Niche-To-Lifestyle-SurveyYear_RemRev.pdf", path = outdir, width = 13, height = 6)

# reorganise some plots 
ggpubr::ggarrange(ggplot(bird_fun_Groups_stack_Lifestyle, aes(x = Treat_type, y = FunG_Proportion, fill = Primary.Lifestyle)) +
                    geom_col(colour = "black") +
                    labs(x = "Treatment", y = "Proportion", fill = "Primary lifestyle") +
                    theme_bw() +
                    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
                    scale_fill_manual(values = Lifestyle_colours),
                  ggplot(bird_fun_Groups_stack, aes(x = Treat_type, y = FunG_Proportion, fill = Trophic.Niche)) +
                    geom_col(colour = "black") +
                    labs(x = "Treatment", y = "Proportion", fill = "Trophic niche") +
                    theme_bw() +
                    scale_fill_manual(values = Niche_colours)+
                    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
                    theme(axis.text.x = element_text(angle = 45, hjust = 1)),
                  ggplot(bird_fun_Groups_stack_Lifestyle, aes(x = Treat_type, y = mean_count, fill = Primary.Lifestyle)) +
                    geom_col(colour = "black") +
                    labs(x = "Treatment", y = "Mean count", fill = "Primary lifestyle") +
                    theme_bw() +
                    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
                    scale_fill_manual(values = Lifestyle_colours),
                  ggplot(bird_fun_Groups_stack, aes(x = Treat_type, y = mean_count, fill = Trophic.Niche)) +
                    geom_col(colour = "black") +
                    labs(x = "Treatment", y = "Mean count", fill = "Trophic niche") +
                    theme_bw() +
                    scale_fill_manual(values = Niche_colours)+
                    facet_grid(~Survey_Year, scales = "free_x", space = "free_x")+
                    theme(axis.text.x = element_text(angle = 45, hjust = 1))
                  , align = "hv")
# ggsave(filename = "Bird-Stack-FUNCTIONS-all-two_stackploats-SurveyYear_RemRev.pdf", path = outdir, width = 13, height = 6)

 ### stats ----------------------------------------------------------------------
bird_functional_CLLMM_Fixed_MODEL <- bird_functional_CLLMM_Fixed_prop
bird_functional_CLLMM_Fixed_MODEL$Survey_Year <- as.factor(bird_functional_CLLMM_Fixed_MODEL$Survey_Year)
# bird_functional_CLLMM_Fixed_MODEL$count <- 1

model_aquatic <- glmer(FunG_count ~ Treat_type + Survey_Year + (1|Site_ID_Var), family = poisson(), 
                       data = bird_functional_CLLMM_Fixed_MODEL[bird_functional_CLLMM_Fixed_MODEL$Trophic.Niche == "Aquatic predator",])
summary(model_aquatic)
#                       Estimate Std. Error z value Pr(>|z|)
# (Intercept)           -0.04742    0.45722  -0.104    0.917
# Treat_typeRevegetated  0.38619    0.46638   0.828    0.408
# Survey_Year2025        0.16063    0.22209   0.723    0.470

model_aquatic_sy <- glmer(FunG_count ~ Treat_type + (1|Survey_Year) + (1|Site_ID_Var), family = poisson(), 
                       data = bird_functional_CLLMM_Fixed_MODEL[bird_functional_CLLMM_Fixed_MODEL$Trophic.Niche == "Aquatic predator",])
summary(model_aquatic_sy) # not sig

model_graniv <- glmer(FunG_count ~ Treat_type + Survey_Year + (1|Site_ID_Var), family = poisson(), 
                      data = bird_functional_CLLMM_Fixed_MODEL[bird_functional_CLLMM_Fixed_MODEL$Trophic.Niche == "Granivore",])
summary(model_graniv)
#                       Estimate Std. Error z value Pr(>|z|)    
# (Intercept)             0.6561     0.1019   6.439  1.2e-10 ***
# Treat_typeRevegetated  -0.1785     0.1191  -1.499    0.134    
# Survey_Year2025         0.0156     0.1187   0.131    0.895    

model_graniv_sy <- glmer(FunG_count ~ Treat_type + (1|Survey_Year) + (1|Site_ID_Var), family = poisson(), 
                      data = bird_functional_CLLMM_Fixed_MODEL[bird_functional_CLLMM_Fixed_MODEL$Trophic.Niche == "Granivore",])
summary(model_graniv_sy) # not sig

model_invert <- glmer(FunG_count ~ Treat_type * Survey_Year + (1|Site_ID_Var), family = poisson(), 
                      data = bird_functional_CLLMM_Fixed_MODEL[bird_functional_CLLMM_Fixed_MODEL$Trophic.Niche == "Invertivore",])
summary(model_invert)
#                                       Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                            1.60995    0.09829  16.380  < 2e-16 ***
# Treat_typeRevegetated                 -0.81722    0.13201  -6.191 5.99e-10 ***
# Survey_Year2025                       -0.28372    0.14342  -1.978  0.04790 *  
# Treat_typeRevegetated:Survey_Year2025  0.58108    0.18585   3.127  0.00177 ** 

model_nect <- glmer(FunG_count ~ Treat_type + Survey_Year + (1|Site_ID_Var), family = poisson(), 
                    data = bird_functional_CLLMM_Fixed_MODEL[bird_functional_CLLMM_Fixed_MODEL$Trophic.Niche == "Nectarivore",])
summary(model_nect)
#                       Estimate Std. Error z value Pr(>|z|)  
# (Intercept)            0.31712    0.18319   1.731   0.0834 .
# Treat_typeRevegetated -0.04905    0.17513  -0.280   0.7794  
# Survey_Year2025        0.07714    0.18040   0.428   0.6690  

model_nect_sy <- glmer(FunG_count ~ Treat_type + (1|Survey_Year) + (1|Site_ID_Var), family = poisson(), 
                    data = bird_functional_CLLMM_Fixed_MODEL[bird_functional_CLLMM_Fixed_MODEL$Trophic.Niche == "Nectarivore",])
summary(model_nect_sy) # not sig

model_omni <- glmer(FunG_count ~ Treat_type + Survey_Year + (1|Site_ID_Var), family = poisson(), 
                    data = bird_functional_CLLMM_Fixed_MODEL[bird_functional_CLLMM_Fixed_MODEL$Trophic.Niche == "Omnivore",])
summary(model_omni)
#                         Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)            0.93739    0.08212  11.415   <2e-16 ***
#   Treat_typeRevegetated -0.19675    0.08758  -2.247   0.0247 *  
#   Survey_Year2025        0.06751    0.08486   0.796   0.4263    

model_vert <- glmer(FunG_count ~ Treat_type + Survey_Year + (1|Site_ID_Var), family = poisson(), 
                    data = bird_functional_CLLMM_Fixed_MODEL[bird_functional_CLLMM_Fixed_MODEL$Trophic.Niche == "Vertivore",])
summary(model_vert)
#                       Estimate Std. Error z value Pr(>|z|)
# (Intercept)            0.00335    0.29466   0.011    0.991
# Treat_typeRevegetated  0.07385    0.39342   0.188    0.851 
# Survey_Year2025       -0.04096    0.73603  -0.056    0.956

model_vert_sy <- glmer(FunG_count ~ Treat_type + (1|Survey_Year) + (1|Site_ID_Var), family = poisson(), 
                    data = bird_functional_CLLMM_Fixed_MODEL[bird_functional_CLLMM_Fixed_MODEL$Trophic.Niche == "Vertivore",])
summary(model_vert_sy) # not sig


# Network analysis with Vegetation:
# Bipartite networks to characterise associations between bird and plant communities based on species co-occurrence across survey sites

# birds
wide_bird_data_site[1:5,1:3]

# save bird wide site data_frame
# saveRDS(object = wide_bird_data_site, file = "wide_bird_data_site_nw.RDS")
library(stringr)
wide_bird_data_site_save <- wide_bird_data_site
wide_bird_data_site_save$sample_id_new <- rownames(wide_bird_data_site_save)

wide_bird_data_site_save$sample_id_new <- str_replace(wide_bird_data_site_save$sample_id_new,
  "^(\\d{4})_([^_]+)_.*$",
  "\\2_\\1"
)
rownames(wide_bird_data_site_save) <- wide_bird_data_site_save$sample_id_new
wide_bird_data_site_save <- wide_bird_data_site_save %>%
  select(-sample_id_new)
# saveRDS(object = wide_bird_data_site_save, file = "wide_bird_data_site_save_nw.RDS")

# ABD data -----
# nrow(bird_functional)
# colnames(bird_functional)
# read_excel("AVONET1_Birdlife.xlsx")

bird_functional_ABD <- read.csv(file = "/PATH/TO/DATA/Bird database stuff/Garnet_2015_Sdata/Australian_Bird_Data_Version_1.csv")

# bird_functional_ABD_FOOD <- bird_functional_ABD %>% 
#   select(matches("X4_"), matches("X5_"), matches("X6_"), matches("X3_"), matches("_Food_"))

bird_functional_ABD_NEST <- bird_functional_ABD %>% 
  select(matches("X4_"), matches("X5_"), matches("X6_"), matches("X3_"), matches("X29_"), matches("X71_"), matches("_Nest_location_"))

# bird_functional_ABD_POPDESC <- bird_functional_ABD %>% 
  # select(matches("X4_"), matches("X5_"), matches("X6_"), matches("X3_"), matches("X29_Population_description_4"))

### ABD Nesting --------------------------------------------------------------------
bird_nest_long <- bird_functional_ABD_NEST %>%
  pivot_longer(
    cols = contains("Nest"),
    names_to = "Nest_location",
    values_to = "Nest_location_1"
  ) %>%
  mutate(
    Nest_location = Nest_location %>%
      str_remove("^X[0-9]+_Nest_location_") %>%   
      str_remove("_12$") %>%      
      str_replace_all("_", " ") %>%
      str_to_sentence()           
  )
unique(bird_nest_long$Nest_location)
colnames(bird_nest_long)


bird_nesting_summary <- bird_nest_long %>%
  filter(Nest_location_1 == 1 ) %>%
  group_by(X4_Genus_name_2, X5_Species_name_2, X6_Subspecies_name_2, X3_Taxon_common_name_2, X29_Population_description_4, X71_South_Australia_7) %>%
  summarise(
    Nesting_levels = list(unique(Nest_location)),
    n_levels = n_distinct(Nest_location),
    Nest_Level = ifelse(n_levels > 1, "Generalist nester", Nest_location[1])
  )
table(bird_nesting_summary$Nest_Level)

bird_nesting_summary$Nesting_niche_Weighting <- 1/bird_nesting_summary$n_levels
bird_nesting_summary$Species <- paste0(bird_nesting_summary$X4_Genus_name_2, " ", bird_nesting_summary$X5_Species_name_2)
unique(bird_nesting_summary$Species)

colnames(bird_nesting_summary) <- c("Genusonly", "Speciesonly", "Subspecies", "Common_name", "Population_description", "Status_South_Australia", "All_nests", "Nesting_Level_n", "Nest_Level", "Nesting_Niche_Weighting", "Species")
bird_nesting_summary2 <- bird_nesting_summary %>%
  ungroup() %>%
  select(-c(Genusonly,Speciesonly))
bird_nesting_summary3 <- bird_nesting_summary2%>%
  select(Species, Subspecies, Common_name, Population_description, Status_South_Australia, Nest_Level, All_nests, Nesting_Niche_Weighting)

# Combine nesting wuth EVONET data
# Alluvial plots----------------------
# View(bird_functional_CLLMM_Fixed)

library(ggalluvial)
ggplot(data = bird_functional_CLLMM_Fixed,
       aes(axis1 = Primary.Lifestyle, axis2 = Trophic.Niche, axis3 = Habitat)) +
  geom_alluvium(aes()) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_test()+
  facet_grid(Survey_Year ~ Treat_type, scales = "free_y")


# Suppose bird_functional_CLLMM_Fixed is your dataframe

# Gather the axes into long format for easy scaling
bird_long_alluvial <- bird_functional_CLLMM_Fixed %>%
  pivot_longer(
    cols = c(Primary.Lifestyle, Trophic.Niche, Habitat),
    names_to = "axis",
    values_to = "stratum"
  )

# Calculate total weights per axis for relative abundance scaling
bird_long_alluvial <- bird_long_alluvial %>%
  group_by(Survey_Year, Treat_type, axis) %>%
  mutate(
    weight = 1  # or if you have an existing weight column, use that here
  ) %>%
  ungroup()

# Calculate sum of weights per axis (for scaling)
axis_sums_alluvial <- bird_long_alluvial %>%
  group_by(Survey_Year, Treat_type, axis) %>%
  summarise(total_weight = sum(weight))

# Join back and compute relative weight per axis
bird_long_alluvial <- bird_long_alluvial %>%
  left_join(axis_sums_alluvial, by = c("axis", "Survey_Year", "Treat_type")) %>%
  mutate(
    rel_weight = weight / total_weight
  )

# Now, pivot wider to recreate the format ggalluvial wants (if needed)
bird_wide_rel_alluvial <- bird_long_alluvial %>%
  select(-weight, -total_weight) %>%
  pivot_wider(
    names_from = axis,
    values_from = stratum
  )

bird_wide_rel_alluvial$Treat_type <- factor(bird_wide_rel_alluvial$Treat_type, levels = c("Revegetated", "Remnant"))

# Use manual fill scale
bird_wide_rel_alluvial$Primary.Lifestyle[bird_wide_rel_alluvial$Primary.Lifestyle == "Insessorial"] <- "Insessorial\n(Perching)"

# species still now found
# Rhipidura albiscapa --> Rhipidura fuliginosa
# Lichenostomus penicillatus --> Ptilotula penicillata  --> Ptilotula penicillatus
# Calyptorhynchus funereus --> Zanda funerea -->  Zanda funereus
# Stigmatopelia chinensis  --> Spilopelia chinensis --> Streptopelia chinensis
# Lichenostomus ornatus --> Ptilotula ornata --> Ptilotula ornatus
# "Cinclosoma castanotum" ---> "Cinclosoma castanotus"
# "Lalage sueurii" --> "Lalage tricolor"

bird_synonyms <- tibble::tribble(
  ~synonym,                       ~Species,
  "Rhipidura albiscapa",           "Rhipidura fuliginosa",
  "Lichenostomus penicillatus",    "Ptilotula penicillatus",
  "Ptilotula penicillata",         "Ptilotula penicillatus",
  "Calyptorhynchus funereus",      "Zanda funereus",
  "Zanda funerea",                 "Zanda funereus",
  "Stigmatopelia chinensis",       "Streptopelia chinensis",
  "Spilopelia chinensis",          "Streptopelia chinensis",
  "Lichenostomus ornatus",         "Ptilotula ornatus",
  "Ptilotula ornata",              "Ptilotula ornatus",
  "Cinclosoma castanotum",         "Cinclosoma castanotus",
  "Lalage sueurii",                "Lalage tricolor",
  "Eolophus roseicapilla",         "Eolophus roseicapillus",
  "Petroica boodang",              "Petroica multicolor"
)

bird_wide_rel_alluvial_2 <- bird_wide_rel_alluvial %>%
  left_join(bird_synonyms, by = c("tSpp" = "synonym")) %>%
  rename(tSpp_std = Species) %>%
  left_join(bird_synonyms, by = c("Species1" = "synonym")) %>%
  rename(Species1_std = Species) %>%
  mutate(
    ABD_preferred = coalesce(Species1_std, tSpp_std, Species1, tSpp)
  )

bird_wide_rel_alluvial3 <- bird_wide_rel_alluvial_2 %>%
  mutate(
    join_species = if_else(
      ABD_preferred %in% bird_nesting_summary3$Species,
      ABD_preferred,
      tSpp
    )
  )

bird_wide_rel_alluvial4 <- bird_wide_rel_alluvial3 %>%
  left_join(
    bird_nesting_summary3,
    by = c("join_species" = "Species")
  )

colnames(bird_wide_rel_alluvial4)
ggplot(data = bird_wide_rel_alluvial4,
       # aes(axis1 = Order1, axis2 = Habitat, axis3 = Primary.Lifestyle, axis4 = Trophic.Niche,
       aes(axis1 = Habitat, axis2 = Primary.Lifestyle, axis3 = Trophic.Niche, axis4 = Nest_Level,
           axis5 = Population_description,
           y = rel_weight)) +
  geom_alluvium(aes(fill = Order1)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  facet_grid(Survey_Year ~ Treat_type, scales = "free")+
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5),  # positions of axis2, axis3, axis4
                     labels = c("Habitat", "Primary Lifestyle", "Trophic Niche", "Nest Location", "Status")
  ) +
  theme_test()+
  ylab("Proportion")
# ggsave(filename = "Alluvial-bird-survey-treat2.pdf", width = 15, height =12)

bird_wide_rel_alluvial4 <- bird_wide_rel_alluvial4 %>%
  mutate(Population_description = recode(
    Population_description,
    "Australian"              = "Native\n(Australia)",
    "Endemic (breeding only)" = "Endemic\n(breeding only)",
    "Endemic (entirely Australian)" = "Endemic\n(entirely Australian)",
    "Non-breeding migrant"    = "Non-breeding migrant",
    "Introduced"              = "Introduced"
  ))

colours_for_all_stratums <- c(
  # axis 1 Habitat
  "Woodland"               = "#228B22",   
  "Forest"                 = "#2E8B57", 
  "Shrubland"              = "#00FA9A",    
  "Human Modified"         = "#32CD32",     
  "Grassland"              = "#6B8E23",   
  "Wetland"                = "#556B2F",   
  "Coastal"                = "#98FB98",   
  " "              = NA,  # spacer
  # Alt axis 1 order
  "Passeriformes"          = "#32CD32", 
  "Psittaciformes"         = "#00FA9A",  
  "Columbiformes"          = "#228B22",     
  "Pelecaniformes"         = "#556B2F",  
  "Charadriiformes"        = "#98FB98",    
  "Accipitriformes"        = "#2E8B57",   
  "Other orders"           = "#6B8E23",
  # "Suliformes"             = "#6B8E23",
  # "Cuculiformes"           = "#6B8E23",   
  # "Falconiformes"          = "#6B8E23",   
  # "Anseriformes"           = "#6B8E23", 
  # "Galliformes"            = "#6B8E23",
  # "Coraciiformes"          = "#6B8E23",   
  " "              = NA,  # spacer
  # axis 2 Primary lifestyle
  "Terrestrial"            = "#00003f",
  "Insessorial\n(Perching)"= "#1F4ED8",
  "Generalist"             = "#6EC6FF",
  "Aerial"                 = "#1CA7A6",
  "Aquatic"                = "#5A6FAE",
  " "              = NA,  # spacer
  # axis 3 Tropic niche
  "Granivore"              = "#DC143C",
  "Invertivore"            = "#B22222",
  "Nectarivore"            = "#FF6347",
  "Omnivore"               = "#FF7F50",
  "Aquatic predator"       = "#FF4500",
  "Vertivore"              = "#FF8C00",
  "Herbivore aquatic"      = "#F4A460",
  " "              = NA,  # spacer
  # axis 4 Nest location
  "Hollow"                 = "#FFD700",
  "Supported"              = "#FFA500",
  "Hanging"                = "#F4A460",
  "Burrow"                 = "#D2B48C",
  "Floating"               = "#B8860B",
  "Ground level"           = "#A0522D",
  "Generalist nester"      = "#8B4513",
  " "              = NA,  # spacer
  # axis 5 Status
  "Endemic\n(entirely Australian)" = "#1B9E77",  
  "Native\n(Australia)"                    = "#66C2A5", 
  "Endemic\n(breeding only)"       = "#A6D854",
  "Non-breeding migrant"          = "#7570B3",
  "Introduced"                    = "#D95F02",
  " "              = NA,  # spacer
  "NA" = "darkgrey"
)

new_order_table <- data.frame(
  Order1 = c("Passeriformes"  , "Psittaciformes" , "Columbiformes", "Pelecaniformes" ,
             "Charadriiformes", "Accipitriformes", "Suliformes",    "Cuculiformes",
             "Falconiformes",   "Anseriformes",    "Galliformes",   "Coraciiformes"),
  Order2 = c("Passeriformes",   "Psittaciformes" , "Columbiformes", "Pelecaniformes",
             "Charadriiformes",  "Accipitriformes","Other orders",  "Other orders", 
             "Other orders",     "Other orders",   "Other orders",  "Other orders")
  )

bird_wide_rel_alluvial5 <- left_join(bird_wide_rel_alluvial4, new_order_table, 
                                            by = "Order1")
# need to do "minor orders" columns
bird_wide_rel_alluvial5$count <- 1

bird_wide_rel_alluvial5$Population_description <- factor(bird_wide_rel_alluvial5$Population_description, 
                                                         levels = c("Non-breeding migrant"          ,
                                                                    "Endemic\n(entirely Australian)",
                                                                    "Endemic\n(breeding only)"      ,
                                                                    "Native\n(Australia)"           ,
                                                                    "Introduced"))

bird_wide_rel_alluvial5$Nest_Level <- factor(bird_wide_rel_alluvial5$Nest_Level, 
                                             levels = c("Generalist nester",
                                                        "Hanging"         ,
                                                        "Hollow"          ,
                                                        "Burrow"          ,
                                                        "Floating"        ,
                                                        "Supported"       ,
                                                        "Ground level"    ))

bird_wide_rel_alluvial5$Order2 <- factor(bird_wide_rel_alluvial5$Order2, 
                                         levels = c("Accipitriformes", "Charadriiformes", 
                                                    "Passeriformes"  , "Columbiformes"  ,
                                                    "Pelecaniformes" , "Psittaciformes" , 
                                                    "Other orders"))

# Plot again
ggplot(data = bird_wide_rel_alluvial5,
       # aes(axis1 = Order1, axis2 = Habitat, axis3 = Primary.Lifestyle, axis4 = Trophic.Niche,
       aes(axis1 = Order2, axis2 = Population_description, axis3 = Primary.Lifestyle, axis4 = Trophic.Niche, axis5 = Nest_Level,
           y = rel_weight)) +
  geom_alluvium(fill = "grey", alpha=0.5) +
  geom_stratum(aes(fill = after_stat(stratum)), alpha = 0.6) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  facet_grid(Survey_Year~Treat_type, scales = "free")+
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5 ),  # positions of axis2, axis3, axis4
                     labels = c("Order", "Status", "Primary Lifestyle", "Trophic Niche", "Nest Location")
  ) +
  scale_fill_manual(values = colours_for_all_stratums)+
  theme_test()+
  ylab("Proportion")
# ggsave(filename = "Alluvial-bird-survey-treat_5col-v2.pdf", width = 20, height =12, path = outdir)

bird_wide_rel_alluvial5_countsav <- bird_wide_rel_alluvial5 %>%
  group_by(Order2, Primary.Lifestyle, Trophic.Niche, Nest_Level, 
           Survey_Year, Treat_type, WptID, iVisitID) %>%
  pivot_longer(
    cols = c(Habitat, Primary.Lifestyle, Trophic.Niche, Nest_Level),
    names_to = "Axis",
    values_to = "Stratum"
  )

bird_alluvial_counts_visit <- bird_wide_rel_alluvial5 %>%
  select(Survey_Year, Treat_type, Site_ID_Var, tSpp, Order2,
         Primary.Lifestyle, Trophic.Niche, Nest_Level, Population_description) %>%
  group_by(Survey_Year, Treat_type, Site_ID_Var) %>%
  distinct() %>% mutate(count=1)

bird_alluvial_counts_visit2 <- bird_alluvial_counts_visit %>%
  select(Survey_Year, Treat_type, Site_ID_Var, tSpp, Order2,
         Primary.Lifestyle, Trophic.Niche, Nest_Level, Population_description, count) %>%
  group_by(Survey_Year, Treat_type,
           Order2,
           Primary.Lifestyle, Trophic.Niche, Nest_Level, Population_description,)%>%
  reframe(meancounts= mean(count),
          prop = mean(count)/sum(count))
  
bird_alluvial_counts_visit2 %>%
  group_by(Survey_Year, Treat_type, Order2, Primary.Lifestyle, 
           Trophic.Niche, Nest_Level, Population_description,)%>%
  reframe(total = sum(prop)) %>% View()
  


group_by(Survey_Year, Treat_type,
         Order2,
         Primary.Lifestyle, Trophic.Niche, Nest_Level, Population_description,)%>%
  bird_alluvial_props <- bird_alluvial_counts_visit %>%
  group_by(Treat_type, Survey_Year) %>%
  mutate(
    prop = n_samples / sum(n_samples)
  ) %>%
  ungroup()

ggplot(data = bird_wide_rel_alluvial5,
       # aes(axis1 = Order1, axis2 = Habitat, axis3 = Primary.Lifestyle, axis4 = Trophic.Niche,
       aes(axis1 = Order2, axis2 = Primary.Lifestyle, axis3 = Trophic.Niche, axis4 = Nest_Level, axis5 = Population_description,
           y = count)) +
  geom_alluvium(fill = "grey") +
  geom_stratum(aes(fill = after_stat(stratum)), alpha = 0.6) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  facet_grid(Survey_Year ~ Treat_type, scales = "free")+
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5),  # positions of axis2, axis3, axis4
                     labels = c("Order", "Primary Lifestyle", "Trophic Niche", "Nest Location", "Status")
  ) +
  scale_fill_manual(values = colours_for_all_stratums)+
  theme_test()+
  ylab("Average count")

## overlay functional group info onto NMDSs ------------------------------------
set.seed(123)
NMDS_Jacc <- metaMDS(wide_bird_data_site, distance = "jaccard", binary = TRUE)
NMDS_Jacc$stress # 0.2199447

plot(NMDS_Jacc)

# determine Species vectors...
sp_fit <- envfit(NMDS_Jacc, wide_bird_data_site, permutations = 999)
plot(sp_fit, p.max = 0.05, col = "red")
# sp_fit

# Niche groups
guild_membership_matrix <- bird_wide_rel_alluvial5 %>%
  distinct(tSpp, Trophic.Niche) %>%   # safety
  mutate(value = 1) %>%
  pivot_wider(
    names_from = Trophic.Niche,
    values_from = value,
    values_fill = 0
  ) %>%
  column_to_rownames("tSpp") %>%
  as.matrix()

dim(guild_membership_matrix)
all(rownames(guild_membership_matrix) %in% colnames(wide_bird_data_site))

guild_membership_matrix <- guild_membership_matrix[
  colnames(wide_bird_data_site),
  ,
  drop = FALSE
]
class(wide_bird_data_site)
class(guild_membership_matrix)
wide_bird_mat <- wide_bird_data_site %>%
  as.data.frame() %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

is.numeric(guild_membership_matrix)

guild_counts <- wide_bird_mat %*% guild_membership_matrix
guild_prop   <- guild_counts / rowSums(guild_counts)

# Lifestyle groups
guild_membership_matrix_LS <- bird_wide_rel_alluvial5 %>%
  distinct(tSpp, Primary.Lifestyle) %>%   # safety
  mutate(value = 1) %>%
  pivot_wider(
    names_from = Primary.Lifestyle,
    values_from = value,
    values_fill = 0
  ) %>%
  column_to_rownames("tSpp") %>%
  as.matrix()

dim(guild_membership_matrix_LS)
all(rownames(guild_membership_matrix_LS) %in% colnames(wide_bird_data_site))

guild_membership_matrix_LS <- guild_membership_matrix_LS[
  colnames(wide_bird_data_site),
  ,
  drop = FALSE
]
class(wide_bird_data_site)
class(guild_membership_matrix_LS)
wide_bird_mat_LS <- wide_bird_data_site %>%
  as.data.frame() %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

guild_counts_LS <- wide_bird_mat_LS %*% guild_membership_matrix_LS
guild_prop_LS   <- guild_counts_LS / rowSums(guild_counts_LS)

# Nesting groups
guild_membership_matrix_NL <- bird_wide_rel_alluvial5 %>%
  distinct(tSpp, Nest_Level) %>%   # safety
  mutate(value = 1) %>%
  pivot_wider(
    names_from = Nest_Level,
    values_from = value,
    values_fill = 0
  ) %>%
  column_to_rownames("tSpp") %>%
  as.matrix()

dim(guild_membership_matrix_NL)
all(rownames(guild_membership_matrix_NL) %in% colnames(wide_bird_data_site))

guild_membership_matrix_NL <- guild_membership_matrix_NL[
  colnames(wide_bird_data_site),
  ,
  drop = FALSE
]

wide_bird_mat_NL <- wide_bird_data_site %>%
  as.data.frame() %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

guild_counts_NL <- wide_bird_mat_NL %*% guild_membership_matrix_NL
guild_prop_NL <- guild_counts_NL / rowSums(guild_counts_NL)

# Population Status
guild_membership_matrix_PD <- bird_wide_rel_alluvial5 %>%
  distinct(tSpp, Population_description) %>%   # safety
  mutate(value = 1) %>%
  pivot_wider(
    names_from = Population_description,
    values_from = value,
    values_fill = 0
  ) %>%
  column_to_rownames("tSpp") %>%
  as.matrix()

dim(guild_membership_matrix_PD)
all(rownames(guild_membership_matrix_PD) %in% colnames(wide_bird_data_site))

guild_membership_matrix_PD <- guild_membership_matrix_PD[
  colnames(wide_bird_data_site),
  ,
  drop = FALSE
]

wide_bird_mat_PD <- wide_bird_data_site %>%
  as.data.frame() %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

guild_counts_PD <- wide_bird_mat_PD %*% guild_membership_matrix_PD
guild_prop_PD <- guild_counts_PD / rowSums(guild_counts_PD)

### examine dfs for vector arrows
# Determine Niche vectors
TN_fit <- envfit(
  NMDS_Jacc,
  guild_prop,
  permutations = 999
)
#                      NMDS1    NMDS2     r2 Pr(>r)    
# Granivore         -0.73202 -0.68128 0.2359  0.001 ***
# Invertivore       -0.63116  0.77566 0.2404  0.001 ***
# Nectarivore       -0.47036  0.88248 0.2330  0.001 ***
# Omnivore           0.66317 -0.74847 0.1532  0.001 ***
# Aquatic predator   0.96713 -0.25430 0.3696  0.001 ***
# Vertivore          0.46652 -0.88451 0.1339  0.002 ** 
# Herbivore aquatic  0.06136 -0.99812 0.0228  0.324    
# NA                -0.45762  0.88915 0.0210  0.353    

PS_fit <- envfit(
  NMDS_Jacc,
  guild_prop_LS,
  permutations = 999
)
PS_fit
#                            NMDS1    NMDS2     r2 Pr(>r)    
# Terrestrial              0.94945 -0.31391 0.3176  0.001 ***
# Insessorial\n(Perching) -0.87651  0.48138 0.6626  0.001 ***
# Generalist               0.94018 -0.34068 0.0436  0.115    
# Aerial                   0.75508 -0.65564 0.3958  0.001 ***
# Aquatic                  0.68175 -0.73158 0.1062  0.005 ** 
# NA                      -0.45762  0.88915 0.0210  0.380    

NL_fit <- envfit(
  NMDS_Jacc,
  guild_prop_NL,
  permutations = 999
)
#                      NMDS1    NMDS2     r2 Pr(>r)    
# Hollow            -0.46240 -0.88667 0.3521  0.001 ***
# Supported         -0.69463  0.71936 0.6774  0.001 ***
# Hanging           -0.68153 -0.73179 0.0980  0.010 ** 
# Burrow            -0.77727 -0.62917 0.0768  0.026 *  
# Floating           0.98200 -0.18886 0.2705  0.001 ***
# Ground level       0.87707 -0.48036 0.5080  0.001 ***
# Generalist nester  0.94743 -0.31998 0.5750  0.001 ***

PD_fit <- envfit(
  NMDS_Jacc,
  guild_prop_PD,
  permutations = 999
)
#                                   NMDS1    NMDS2     r2 Pr(>r)    
# Endemic\n(entirely Australian) -0.51506  0.85715 0.3523  0.001 ***
# Native\n(Australia)             0.04801 -0.99885 0.0376  0.159    
# Non-breeding migrant           -0.01406 -0.99990 0.0694  0.030 *  
# Introduced                      0.60214 -0.79839 0.2424  0.001 ***
# Endemic\n(breeding only)        0.79110 -0.61169 0.3694  0.001 ***

# fit to dfs
PS_fit_values <- as.data.frame(PS_fit$vectors$arrows)
PS_fit_values$Primary.Lifestyle <- rownames(PS_fit_values)
PS_fit_values$Sign <- ifelse(PS_fit$vectors$pvals <0.05, 1, 0)

TN_fit_values <- as.data.frame(TN_fit$vectors$arrows)
TN_fit_values$Tropic.Niche <- rownames(TN_fit_values)
TN_fit_values$Sign <-  ifelse(TN_fit$vectors$pvals <0.05, 1, 0)

NL_fit_values <- as.data.frame(NL_fit$vectors$arrows)
NL_fit_values$Nest_Level <- rownames(NL_fit_values)
NL_fit_values$Sign <-  ifelse(NL_fit$vectors$pvals <0.05, 1, 0)

PD_fit_values <- as.data.frame(PD_fit$vectors$arrows)
PD_fit_values$Pop_Status <- rownames(PD_fit_values)
PD_fit_values$Sign <-  ifelse(PD_fit$vectors$pvals <0.05, 1, 0)

# Replot earlier NMDSs
ggplot()+
  # Tropic niche vectors
  geom_segment(data = TN_fit_values[TN_fit_values$Sign==1,],colour = "red", aes(xend = NMDS1, yend=NMDS2, x=0, y=0))+
  geom_text(data = TN_fit_values[TN_fit_values$Sign==1,], colour = "red", aes(x = NMDS1, y=NMDS2, label=Tropic.Niche))+
  # Primary Lifestyle vectors
  geom_segment(data = PS_fit_values[PS_fit_values$Sign==1,], colour = "blue", aes(xend = NMDS1, yend=NMDS2, x=0, y=0))+
  geom_text(data = PS_fit_values[PS_fit_values$Sign==1,], colour = "blue", aes(x = NMDS1, y=NMDS2, label=Primary.Lifestyle))+
  # Nest Location vectors
  geom_segment(data = NL_fit_values[NL_fit_values$Sign==1,], colour = "green", aes(xend = NMDS1, yend=NMDS2, x=0, y=0))+
  geom_text(data = NL_fit_values[NL_fit_values$Sign==1,], colour = "darkgreen", aes(x = NMDS1, y=NMDS2, label=Nest_Level))+
  # Population Description
  geom_segment(data = PD_fit_values[PD_fit_values$Sign==1,], colour = "orange", aes(xend = NMDS1, yend=NMDS2, x=0, y=0))+
  geom_text(data = PD_fit_values[PD_fit_values$Sign==1,], colour = "darkorange", aes(x = NMDS1, y=NMDS2, label=Pop_Status))+
  # points
  geom_point(data = bird_metadata_s_Jacc, size= 2, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =Treat_type, hjust=1, vjust=1))+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()+
  labs(x= "NMDS1", y="NMDS2")

TN_NMDS <- ggplot()+
  # Tropic niche vectors
  geom_segment(data = TN_fit_values[TN_fit_values$Sign==1,],colour = "red", aes(xend = NMDS1, yend=NMDS2, x=0, y=0))+
  geom_text(data = TN_fit_values[TN_fit_values$Sign==1,], colour = "red", aes(x = NMDS1, y=NMDS2, label=Tropic.Niche))+
  # Points
  geom_point(data = bird_metadata_s_Jacc, size= 3, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =Treat_type))+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()+
  ggtitle(paste0("Stress: ", round(NMDS_Jacc$stress, 4)))+
  labs(x= "NMDS1", y="NMDS2")

PS_NMDS <- ggplot()+
  # Primary Lifestyle vectors
  geom_segment(data = PS_fit_values[PS_fit_values$Sign==1,], colour = "blue", aes(xend = NMDS1, yend=NMDS2, x=0, y=0))+
  geom_text(data = PS_fit_values[PS_fit_values$Sign==1,], colour = "blue", aes(x = NMDS1, y=NMDS2, label=Primary.Lifestyle))+
  # Points
  geom_point(data = bird_metadata_s_Jacc, size= 3, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =Treat_type))+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()+
  ggtitle(paste0("Stress: ", round(NMDS_Jacc$stress, 4)))+
  labs(x= "NMDS1", y="NMDS2")

NL_NMDS <- ggplot()+
  # Nest Location vectors
  geom_segment(data = NL_fit_values[NL_fit_values$Sign==1,], colour = "green", aes(xend = NMDS1, yend=NMDS2, x=0, y=0))+
  geom_text(data = NL_fit_values[NL_fit_values$Sign==1,], colour = "darkgreen", aes(x = NMDS1, y=NMDS2, label=Nest_Level))+
  # Points
  geom_point(data = bird_metadata_s_Jacc, size= 3, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =Treat_type))+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()+
  ggtitle(paste0("Stress: ", round(NMDS_Jacc$stress, 4)))+
  labs(x= "NMDS1", y="NMDS2")

PD_NMDS <- ggplot()+
  # Nest Location vectors
  geom_segment(data = PD_fit_values[PD_fit_values$Sign==1,], colour = "orange", aes(xend = NMDS1, yend=NMDS2, x=0, y=0))+
  geom_text(data = PD_fit_values[PD_fit_values$Sign==1,], colour = "darkorange", aes(x = NMDS1, y=NMDS2, label=Pop_Status))+
  # Points
  geom_point(data = bird_metadata_s_Jacc, size= 3, aes(x = MDS1, y = MDS2, shape = Survey_Year, colour =Treat_type))+
  scale_shape_manual(values =  c("2015" = 21, "2025" = 15))+
  theme_test()+
  ggtitle(paste0("Stress: ", round(NMDS_Jacc$stress, 4)))+
  labs(x= "NMDS1", y="NMDS2")

ggpubr::ggarrange(TN_NMDS, PS_NMDS, NL_NMDS, PD_NMDS, common.legend = TRUE)
# ggsave(filename = "Bird-BETA-NMDS-Jaccard-SurveyYear_PlantYear-ALLVECTORS.pdf", path = outdir, width = 14, height = 12)

# Additional ABD data stuff
### ABD Food and Diet --------------------------------------------------------------
bird_food_long <- bird_functional_ABD_FOOD %>%
  pivot_longer(
    cols = contains("Food"),
    names_to = "Food_type",
    values_to = "Food_present"
  ) %>%
  mutate(
    Food_type = Food_type %>%
      str_remove("^X[0-9]+_Food_") %>%   # drop X163_Food_
      str_remove("_10$") %>%             # drop _10
      str_replace_all("_", " ") %>%      # underscores to spaces
      str_to_sentence()                  # nice capitalisation
  )
unique(bird_food_long$Food_type)

# Step 1: Assign each Food_type a basic trophic category (e.g., Herbivore, Carnivore, etc.) for each row
# Step 2: Summarize per bird to check if multiple trophic categories exist
# ####### Replace Bird_ID with the actual bird ID column name
# Step 3: Join back the Trophic_Level to the original data if needed

bird_food_long2 <- bird_food_long %>%
  filter(Food_present == 1) %>%
  mutate(
    Basic_Trophic = case_when(
      Food_type %in% c("Foliage or herbs", "Fruit", "Seeds", "Corms or tubers") ~ "Herbivore",
      Food_type == "Nectar or pollen" ~ "Nectarivore",
      Food_type %in% c("Terrestrial vertebrates", "Terrestrial invertebrates", 
                       "Intertidal invertebrates", "Fish or invertebrates marine", 
                       "Fish or invertebrates inland waters") ~ "Carnivore",
      Food_type == "Carrion" ~ "Scavenger",
      TRUE ~ NA_character_
    )
  )

bird_trophic_summary <- bird_food_long2 %>%
  group_by(X4_Genus_name_2, X5_Species_name_2, X6_Subspecies_name_2, X3_Taxon_common_name_2) %>%
  summarise(
    trophic_levels = list(unique(Basic_Trophic)),
    n_levels = n_distinct(Basic_Trophic),
    Trophic_Level = ifelse(n_levels > 1, "Omnivore", Basic_Trophic[1])
  )

bird_food_long3 <- bird_food_long2 %>%
  left_join(select(bird_trophic_summary, X4_Genus_name_2, X5_Species_name_2, X6_Subspecies_name_2, X3_Taxon_common_name_2, Trophic_Level), by = c("X4_Genus_name_2", "X5_Species_name_2", "X6_Subspecies_name_2", "X3_Taxon_common_name_2")) %>%
  # Then apply your Trophic_Niche as before (if needed)
  mutate(
    Trophic_Niche = case_when(
      Food_type == "Fruit" ~ "Frugivore",
      Food_type == "Seeds" ~ "Granivore",
      Food_type %in% c("Foliage or herbs", "Corms or tubers") ~ "Herbivore",
      Food_type == "Nectar or pollen" ~ "Nectarivore",
      Food_type %in% c("Terrestrial invertebrates", "Intertidal invertebrates") ~ "Invertivore",
      Food_type %in% c("Fish or invertebrates marine", "Fish or invertebrates inland waters") ~ "Piscivore",
      Food_type == "Terrestrial vertebrates" ~ "Vertivore",
      Food_type == "Carrion" ~ "Scavenger",
      TRUE ~ NA_character_
    )
  )

bird_trophic_summary2 <- bird_food_long3 %>%
  group_by(X4_Genus_name_2, X5_Species_name_2, X6_Subspecies_name_2, X3_Taxon_common_name_2) %>%
  summarise(
    All_Trophic_Levels = list(unique(Basic_Trophic)),
    n_Trophic_levels = n_distinct(Basic_Trophic),
    Trophic_Level = ifelse(n_Trophic_levels > 1, "Omnivore", Basic_Trophic[1]),
    
    All_Trophic_Niches = list(unique(Trophic_Niche)),
    n_Trophic_niches = n_distinct(Trophic_Niche)
    )

bird_trophic_summary2$Trophic_niches_Weighting <- 1/bird_trophic_summary2$n_Trophic_niches
bird_trophic_summary2$Species <- paste0(bird_trophic_summary2$X4_Genus_name_2, " ", bird_trophic_summary2$X5_Species_name_2)

colnames(bird_trophic_summary2)
bird_trophic_summary3 <- bird_trophic_summary2 %>%
  select(Species, X6_Subspecies_name_2, X3_Taxon_common_name_2, Trophic_Level, All_Trophic_Niches, Trophic_niches_Weighting, X4_Genus_name_2, X5_Species_name_2) %>%
  ungroup() %>%
  as.data.frame()

nrow(bird_trophic_summary3)
length(unique(bird_trophic_summary3$Species))

colnames(bird_trophic_summary3) <- c("Species", "Subspecies", "Common_name", "Trophic_Level", "Trophic_Niche", "Trophic_Niche_Weighting", "Genusonly", "Speciesonly")
bird_trophic_summary3 <- bird_trophic_summary3 %>%
  select(-c(Genusonly,Speciesonly))

head(bird_trophic_summary3)
