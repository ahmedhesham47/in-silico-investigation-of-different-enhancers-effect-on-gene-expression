library(mpra)
library(dplyr)

treatn <- function(fit, lfc = log2(1.2), trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05, 0.1)) {
  # Check fit
  if (!is(fit, "MArrayLM")) stop("fit must be an MArrayLM object")
  if (is.null(fit$coefficients)) stop("coefficients not found in fit object")
  if (is.null(fit$stdev.unscaled)) stop("stdev.unscaled not found in fit object")
  
  fit$lods <- NULL
  
  # Extract the coefficients and other information
  coefficients <- as.matrix(fit$coefficients)
  stdev.unscaled <- as.matrix(fit$stdev.unscaled)
  sigma <- fit$sigma
  df.residual <- fit$df.residual
  
  # Check if we have all necessary components
  if (is.null(coefficients) || is.null(stdev.unscaled) || is.null(sigma) || is.null(df.residual)) 
    stop("No data, or argument is not a valid lmFit object")
  if (max(df.residual) == 0) stop("No residual degrees of freedom in linear model fits")
  if (!any(is.finite(sigma))) stop("No finite residual standard deviations")
  
  # For trend calculation, if needed
  if (trend) {
    covariate <- fit$Amean
    if (is.null(covariate)) stop("Need Amean component in fit to estimate trend")
  } else {
    covariate <- NULL
  }
  
  sv <- squeezeVar(sigma^2, df.residual, covariate = covariate, robust = robust, winsor.tail.p = winsor.tail.p)
  fit$df.prior <- sv$df.prior
  fit$s2.prior <- sv$var.prior
  fit$s2.post <- sv$var.post
  
  df.total <- df.residual + sv$df.prior
  df.pooled <- sum(df.residual, na.rm = TRUE)
  df.total <- pmin(df.total, df.pooled)
  fit$df.total <- df.total
  
  lfc <- abs(lfc)
  acoef <- abs(coefficients)
  
  # Calculate the standard error
  se <- stdev.unscaled * sqrt(fit$s2.post)
  
  # Calculate the t-statistics
  tstat.right <- (acoef - lfc) / se
  tstat.left <- (acoef + lfc) / se
  fit$t <- array(0, dim(coefficients), dimnames = dimnames(coefficients))
  
  fit$p.value <- pt(tstat.right, df = df.total, lower.tail = FALSE) + pt(tstat.left, df = df.total, lower.tail = FALSE)
  
  # Add standard errors (SE) to the results
  fit$SE = se 
  
  tstat.right <- pmax(tstat.right, 0)
  fc.up <- (coefficients > lfc)
  fc.down <- (coefficients < -lfc)
  fit$t[fc.up] <- tstat.right[fc.up]
  fit$t[fc.down] <- -tstat.right[fc.down]
  fit$treat.lfc <- lfc
  
  return(fit)
}

counts_path_upstream = '/project2/kribelba_1515/saadawy/df_for_mpra_upstream.tsv'
counts_path_UPSTREAM = '/project2/kribelba_1515/saadawy/df_for_mpra.tsv'
counts_path_downstream = '/project2/kribelba_1515/saadawy/df_for_mpra_downstream.tsv'


df_upstream <- read.csv(counts_path_upstream, sep = '\t')
df_UPSTREAM <- read.csv(counts_path_UPSTREAM, sep = '\t')
df_downstream <- read.csv(counts_path_downstream, sep = '\t')

df_upstream <- df_upstream %>%
  group_by(Enhancer) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  rename_with(~ paste0(., "_upstream"), -Enhancer)

df_UPSTREAM <- df_UPSTREAM %>%
  group_by(Enhancer) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  rename_with(~ paste0(., "_UPSTREAM"), -Enhancer)

df_downstream <- df_downstream %>%
  group_by(Enhancer) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  rename_with(~ paste0(., "_downstream"), -Enhancer)

merged_df <- df_upstream %>%
  inner_join(df_UPSTREAM, by='Enhancer') %>%
  inner_join(df_downstream, by='Enhancer')

gDNA_counts <- merged_df %>% dplyr::select(starts_with("gDNA_rep")) %>% as.matrix()
mRNA_counts <- merged_df %>% dplyr::select(starts_with("mRNA_rep")) %>% as.matrix()

eid_vector <- as.character(merged_df$Enhancer)
rownames(gDNA_counts) <- rownames(mRNA_counts) <- eid_vector

colnames(gDNA_counts) <- gsub("gDNA_", "", colnames(gDNA_counts))
colnames(mRNA_counts) <- gsub("mRNA_", "", colnames(mRNA_counts))

mpra_obj <- MPRASet(
  RNA = mRNA_counts,
  DNA = gDNA_counts,
  eid = eid_vector,
)

n_upstream <- 4
n_UPSTREAM <- 5
n_downstream <- 3

total <- n_upstream + n_UPSTREAM + n_downstream

design <- data.frame(
  intcpt = rep(1, total),
  upstream = c(rep(TRUE, n_upstream), rep(FALSE, n_UPSTREAM + n_downstream)),
  UPSTREAM = c(rep(FALSE, n_upstream), rep(TRUE, n_UPSTREAM), rep(FALSE, n_downstream)),
  downstream = c(rep(FALSE, n_upstream + n_UPSTREAM), rep(TRUE, n_downstream))
)

group <- factor(c(rep(1, 4), rep(2, 5), rep(3, 3)))

design <- model.matrix(~ 0 + group)

output_dir <- "/project2/kribelba_1515/saadawy/mpra/integrative_design"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

fit <- mpralm(object = mpra_obj, design = design,
              aggregate = 'none', normalize = FALSE,
              model_type = 'indep_groups', block = NULL, plot = TRUE)

logr = compute_logratio(mpra_obj, aggregate='none')

tr <- treatn(fit)
result <- topTreat(tr, coef = 1, number = Inf)

result$uid <- rownames(result)
result$stdev.unscaled <- fit$stdev.unscaled[match(result$uid, rownames(fit$stdev.unscaled))]
result$SE <- tr$SE[match(result$uid, rownames(tr$SE))]
colnames(result)[length(result[1,])] <- 'SE'

out_path <- file.path(output_dir, paste0("integrative_conditions_mpra", ".csv"))
write.csv(result, out_path, row.names = FALSE)

# _____________________________________________________________

mpralm_fit_results <- list()

normalize_options <- c(TRUE, FALSE)
aggregation_options <- c('sum', 'mean')

for (norm in normalize_options) {
  for (agg in aggregation_options) {
    
    # Run mpralm with the current combination of parameters
    mpralm_fit <- mpralm(object = mpra_obj, design = design,
                         aggregate = agg, normalize = norm,
                         model_type = 'indep_groups', block = NULL, plot = TRUE)
    key <- paste0("agg_", agg, "_norm_", norm, "_model_indep_groups")
    mpralm_fit_results[[key]] <- mpralm_fit
    
    plot_obj <- recordPlot()  # Captures the last plot generated
    
    plot_filename <- paste0("mpralm_", 'agg_', agg, "_norm_", norm, "_model_type_", 'indep_groups', ".png")
    plot_path <- file.path(output_dir, plot_filename)
    
    png(plot_path, width = 800, height = 600)
    replayPlot(plot_obj)
    dev.off()
  }
}

for (norm in normalize_options) {
  for (agg in aggregation_options) {
    key <- paste0("agg_", agg, "_norm_", norm, "_model_indep_groups")
    
    fit <- mpralm(object = mpra_obj, design = design,
                  aggregate = agg, normalize = norm,
                  model_type = 'indep_groups', block = NULL, plot = TRUE)
    
    mpralm_fit_results[[key]] <- fit
    
    tr <- treatn(fit)
    result <- topTreat(tr, coef = 1, number = Inf)
    
    result$uid <- rownames(result)
    result$stdev.unscaled <- fit$stdev.unscaled[match(result$uid, rownames(fit$stdev.unscaled))]
    result$SE <- tr$SE[match(result$uid, rownames(tr$SE))]
    colnames(result)[length(result[1,])] <- 'SE'
    
    out_path <- file.path(output_dir, paste0(key, ".csv"))
    write.csv(result, out_path, row.names = FALSE)
    message("Saved result: ", out_path)
  }
}