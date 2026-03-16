# Permuted Linear Mixed Effects Models (LMEM) Toolkit
# RJH 2026

# This toolkit fits linear mixed effects models and evaluates fixed effects
# using permutation (randomisation) tests.
#
# Users can specify any number of fixed and random effects. Statistical
# significance is assessed by repeatedly permuting the response variable
# (default = 9,999 permutations) and refitting the same model to generate
# a null distribution of test statistics.
#
# Fixed effects are tested using permutation-based inference rather than
# asymptotic assumptions, making this approach robust to small sample sizes
# and non-normal residuals.
#
# IMPORTANT:
# - Interaction terms must be specified explicitly in the fixed effects.
# - Interactions should be written using ':' rather than '*'.
#
# Example:
# fixed_effects = c("Factor1", "Factor2", "Factor1:Factor2")
#
# Random effects are specified as grouping variables and are included as
# random intercepts.
#
# Example:
# random_effects = c("SiteID", "Block")
#
# The output includes:
# - The fitted LMEM
# - Observed test statistics for each fixed effect
# - Permutation-based null distributions
# - Permutation p-values for each fixed effect


# Function using Summary(LMEM_model)
LMEM_permute_summary <- function(data, response, fixed_effects, random_effects, nreps = 9999, seed = 123) {
  set.seed(seed)
  
  fixed_formula <- if(length(fixed_effects) == 1) {
    fixed_effects
  } else {
    paste(fixed_effects, collapse = " + ")
  }
  
  random_formula <- paste0("(1 | ", random_effects, ")", collapse = " + ")
  full_formula_str <- paste(response, "~", fixed_formula, "+", random_formula)
  full_formula <- as.formula(full_formula_str)
  
  model <- lmer(full_formula, data = data)
  coef_table <- summary(model)$coefficients
  
  fixed_coef_names <- rownames(coef_table)[-1]  # exclude intercept
  
  bootstrap_tvalues <- matrix(NA, nrow = nreps, ncol = length(fixed_coef_names))
  colnames(bootstrap_tvalues) <- fixed_coef_names
  N <- nrow(data)
  
  for(i in seq_len(nreps)) {
    data_boot <- data
    data_boot[[response]] <- sample(data[[response]], N, replace = FALSE)
    
    mod_boot <- tryCatch(
      lmer(full_formula, data = data_boot),
      error = function(e) NULL
    )
    if (!is.null(mod_boot)) {
      boot_coefs_i <- summary(mod_boot)$coefficients
      if(all(fixed_coef_names %in% rownames(boot_coefs_i))) {
        bootstrap_tvalues[i, ] <- boot_coefs_i[fixed_coef_names, "t value"]
      } else {
        bootstrap_tvalues[i, ] <- NA
      }
    } else {
      bootstrap_tvalues[i, ] <- NA
    }
  }
  
  bootstrap_tvalues <- bootstrap_tvalues[complete.cases(bootstrap_tvalues), ]
  observed_tvalues <- coef_table[fixed_coef_names, "t value"]
  
  p_values <- sapply(seq_along(fixed_coef_names), function(j) {
    boot_dist <- bootstrap_tvalues[, j]
    obs_val <- observed_tvalues[j]
    mean(abs(boot_dist) >= abs(obs_val))
  })
  names(p_values) <- fixed_coef_names
  
  return(list(
    formula = full_formula,
    model = model,
    summary = summary(model),
    observed_tvalues = observed_tvalues,
    bootstrap_tvalues = bootstrap_tvalues,
    permutation_p_values = p_values
  ))
}

# Function using car::Anova(LMEM_model)
LMEM_permute_anova <- function(data, response, fixed_effects, random_effects, nreps = 9999, seed = 123) {
  
  set.seed(seed)
  
  # Build fixed effects
  fixed_formula <- if (length(fixed_effects) == 1) {
    fixed_effects
  } else {
    paste(fixed_effects, collapse = " + ")
  }
  
  # Build random effects
  random_formula <- paste0("(1 | ", random_effects, ")", collapse = " + ")
  
  # Full model formula
  full_formula <- as.formula(
    paste(response, "~", fixed_formula, "+", random_formula)
  )
  
  # Fit observed model
  model <- lmer(full_formula, data = data)
  anova_obs <- car::Anova(model)
  
  fixed_names <- rownames(anova_obs)
  observed_chi <- anova_obs$Chisq
  
  # Storage for permutation stats
  permute_teststats <- matrix(
    NA, nrow = nreps, ncol = length(fixed_names),
    dimnames = list(NULL, fixed_names)
  )
  
  N <- nrow(data)
  
  # Permutation loop
  for (i in seq_len(nreps)) {
    
    data_perm <- data
    data_perm[[response]] <- sample(data[[response]], N, replace = FALSE)
    
    mod_perm <- tryCatch(
      lmer(full_formula, data = data_perm),
      error = function(e) NULL
    )
    
    if (!is.null(mod_perm)) {
      anova_perm <- tryCatch(
        car::Anova(mod_perm),
        error = function(e) NULL
      )
      
      if (!is.null(anova_perm) &&
          all(fixed_names %in% rownames(anova_perm))) {
        permute_teststats[i, ] <- anova_perm[fixed_names, "Chisq"]
      }
    }
  }
  
  # Drop failed permutations
  permute_teststats <- permute_teststats[complete.cases(permute_teststats), ]
  
  # Permutation p-values (one-sided, χ²)
  p_values <- sapply(seq_along(fixed_names), function(j) {
    nbout <- which(permute_teststats[, j] >= observed_chi[j])
    p_perm <- length(nbout) / length(permute_teststats[, j])
    return(p_perm)
  })
  
  names(p_values) <- fixed_names
  
  return(list(
    formula = full_formula,
    model = model,
    Anova = anova_obs,
    observed_chisq = observed_chi,
    permuted_chisq = permute_teststats,
    permutation_p_values = p_values
  ))
}

# Genearlised LMEM (binomial) using Summary(GLMEM_model) for effect size and ANOVA(GLMEM_model) to permute Chisq
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


# This function will give the pairwise t test results with multiple corrections
# it is intedned for use after permuted ANOVAs (aov()) with a single explanatory variable that requires multiple corrections
library(multcompView)

# Make sure levels of the group/treatment factor are as desired
pairwise_perm_cld <- function(df, response, treatment, nreps = 9999, p_adjust = "BH") {
  # Ensure response is numeric and treatment is factor
  df[[response]] <- as.numeric(as.character(df[[response]]))
  df[[treatment]] <- as.factor(df[[treatment]])
  
  trts <- levels(df[[treatment]])
  
  # All pairwise combinations
  pairs <- combn(trts, 2, simplify = FALSE)
  
  # Permutation t-test function
  perm_t_test <- function(y, group, nreps_inner) {
    obs <- t.test(y ~ group)$statistic
    r <- replicate(nreps_inner, {
      new_y <- sample(y)
      t.test(new_y ~ group)$statistic
    })
    pv <- sum(abs(r) >= abs(obs))/nreps_inner
    data.frame(t_stat = obs, pvalue = pv)
  }
  
  # Run all pairwise permutation tests
  pvals <- sapply(pairs, function(pair) {
    sub <- df[df[[treatment]] %in% pair, ]
    perm_t_test(sub[[response]], sub[[treatment]], nreps_inner = nreps)
  })
  
  # Format results
  pvals <- t(pvals)
  pvals <- as.data.frame(pvals)
  pvals$comparison <- sapply(pairs, paste, collapse = " vs ")
  pvals$pvalue_adj <- p.adjust(pvals$pvalue, method = p_adjust)
  pvals <- pvals[, c("comparison", "t_stat", "pvalue", "pvalue_adj")]
  pvals$significance <- ifelse(pvals$pvalue_adj < 0.05, "*", 
                               ifelse(pvals$pvalue_adj < 0.01, "**", 
                                      ifelse(pvals$pvalue_adj < 0.001, "***", "ns")))
  
  # Build significance matrix
  sig_mat <- matrix(TRUE, nrow = length(trts), ncol = length(trts), dimnames = list(trts, trts))
  for(i in seq_along(pairs)){
    t1 <- pairs[[i]][1]
    t2 <- pairs[[i]][2]
    sig_mat[t1, t2] <- sig_mat[t2, t1] <- pvals$pvalue_adj[i] < 0.05
  }
  
  # Generate compact letter display
  cld <- multcompLetters(sig_mat, threshold = 0.05)
  
  # Return results
  list(
    pairwise_comparisons = pvals,
    group_differences = cld$Letters
  )
}
