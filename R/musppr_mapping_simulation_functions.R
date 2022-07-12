#' Simulate MPP data with effects from a single QTL, strain, and kinship
#'
#' This function ...
#' 
#' @export
#' @examples sim_mpp_data()
sim_mpp_data <- function (genoprobs, map, K = NULL, samples = NULL,
                          qtl_effect_size = 0.1, strain_effect_size = 0, kinship_effect_size = 0.1,
                          locus = NULL, vary_locus = TRUE, 
                          num_replicates = 1, num_sim = 1, 
                          M_ID = NULL, sample_method = c("uniform", "crp"), num_alleles = 8, 
                          scale_type = c("sample", "balanced"),
                          num_founders = 8, beta = NULL, sim_label = "sim_y") {
  
  sample_method <- sample_method[1]
  scale_type <- scale_type[1]
  if (is.null(samples)) {
    samples <- dimnames(genoprobs[[1]])[[1]]
  }
  num_genomes <- length(samples)
  
  ## Adjust effect sizes for replicates
  original_effects <- list(qtl_effect_size = qtl_effect_size, 
                           strain_effect_size = strain_effect_size,
                           kinship_effect_size = kinship_effect_size)
  noise_effect_size <- 1 - qtl_effect_size - strain_effect_size - kinship_effect_size
  noise_effect_size <- noise_effect_size/num_replicates
  scale_factor <- (1 - noise_effect_size)/sum(c(qtl_effect_size, strain_effect_size, kinship_effect_size))
  qtl_effect_size <- scale_factor * qtl_effect_size
  strain_effect_size <- scale_factor * strain_effect_size
  kinship_effect_size <- scale_factor * kinship_effect_size

  map_df <- map_list_to_df(map)
  
  if (!is.null(beta)) {
    num_alleles <- length(beta)
  }
  loci <- unlist(lapply(genoprobs, function(x) dimnames(x)[[3]]))
  if (is.null(locus)) {
    locus <- sample(loci, size = ifelse(vary_locus, num_sim, 1), replace = TRUE)
  }
  if (vary_locus) {
    locus_index <- 1:num_sim
  }
  else {
    locus_index <- rep(1, num_sim)
  }
  
  M <- NULL
  if (!is.null(M_ID)) {
    M <- model_matrix_from_ID(M_ID)
    num_alleles <- length(unique(unlist(strsplit(x = M_ID, 
                                                 split = ","))))
  }
  sim_matrix <- matrix(NA, nrow = num_genomes, ncol = num_sim)
  for (i in 1:num_sim) {
    ## LOCO K
    if (is.list(K)) {
      this_K <- K[[map_df[locus[i], "chr"]]]
    } else {
      this_K <- K
    }
    locus_matrix <- qtl2::pull_genoprobpos(genoprobs = genoprobs, marker = locus[locus_index[i]])
    this_sim <- sim_mpp_qtl(M = M, beta = beta, sample_method = sample_method, 
                            num_alleles = num_alleles, num_founders = num_founders, K = this_K[samples, samples],
                            qtl_effect_size = qtl_effect_size, strain_effect_size = strain_effect_size, 
                            kinship_effect_size = kinship_effect_size, noise_effect_size = noise_effect_size, 
                            scale_type = scale_type, locus_matrix = locus_matrix[samples,], 
                            num_sim = 1, sim_label = sim_label)
    sim_matrix[, i] <- this_sim$data
    if (i == 1) {
      rownames_holder <- rownames(this_sim$data)
    }
  }
  colnames(sim_matrix) <- paste(sim_label, 1:num_sim, sep = "_")
  rownames(sim_matrix) <- rownames_holder
  
  return(list(data = sim_matrix, 
              locus = as.character(locus), 
              locus_pos = as.numeric(map_df[locus, "pos"]), locus_chr = as.character(map_df[locus, "chr"]), 
              properties = list(num_alleles = num_alleles, 
                                sample_method = sample_method, 
                                num_replicates = num_replicates, num_founders = num_founders, 
                                qtl_effect_size = original_effects$qtl_effect_size, 
                                strain_effect_size = original_effects$strain_effect_size, 
                                kinship_effect_size = original_effects$kinship_effect_size,
                                M_ID = M_ID, vary_locus = vary_locus,
                                scale_type = scale_type)))
}


sim_mpp_qtl <- function (locus_matrix, 
                         M = NULL, sample_method = c("uniform", "crp"), 
                         scale_type = c("sample", "balanced"),
                         qtl_effect_size, beta = NULL, strain_effect_size, 
                         kinship_effect_size, K = NULL,
                         noise_effect_size, num_alleles = 8, num_founders = 8, num_sim, 
                         sim_label = "sim_y", 
                         ...) {
  
  sample_method <- sample_method[1]
  scale_type <- scale_type[1]
  strains <- rownames(locus_matrix)

  D <- locus_matrix
  
  ## QTL effect
  if (qtl_effect_size != 0) {
    qtl_effect <- sim_M_and_beta(num_alleles = num_alleles, 
                                 num_founders = num_founders, M = M, sample_method = sample_method, 
                                 beta = beta, ...)
  }
  else {
    qtl_effect <- list(M = diag(8), beta = rep(0, 8))
  }
  M <- qtl_effect$M
  raw_beta <- qtl_effect$beta
  if (qtl_effect_size != 0) {
    beta <- (raw_beta - mean(raw_beta))/sqrt(non_sample_var(raw_beta))
  }
  else {
    beta <- raw_beta
  }
  if (qtl_effect_size != 0) {
    if (scale_type == "sample") {
      DMB <- D %*% M %*% beta
      qtl_predictor <- DMB * sqrt(qtl_effect_size/non_sample_var(DMB)[1])
    } else {
      MB <- M %*% beta
      qtl_predictor <- D %*% MB * sqrt(qtl_effect_size/non_sample_var(MB)[1])
    }
  }
  else {
    qtl_predictor <- D %*% M %*% beta
  }

  scale_factor <- 1/sum(c(non_sample_var(qtl_predictor), strain_effect_size, kinship_effect_size, noise_effect_size))
  ## Strain effect
  if (strain_effect_size != 0) {
    strain_effect <- rnorm(n = nrow(locus_matrix))
    strain_effect <- (strain_effect - mean(strain_effect))/sqrt(non_sample_var(strain_effect))
    strain_predictor <- strain_effect * sqrt(strain_effect_size * scale_factor)
  }
  else {
    strain_predictor <- rep(0, nrow(D))
  }
  
  ## Kinship effect
  if (kinship_effect_size != 0 & !is.null(K)) {
    kinship_effect <- MASS::mvrnorm(n = 1, mu = rep(0, nrow(locus_matrix)), Sigma = K)
    kinship_effect <- (kinship_effect - mean(kinship_effect))/sqrt(non_sample_var(kinship_effect))
    kinship_predictor <- kinship_effect * sqrt(kinship_effect_size * scale_factor)
  } else {
    kinship_predictor <- rep(0, nrow(locus_matrix))
  }
  
  sim_data <- matrix(NA, nrow = nrow(D), ncol = num_sim)
  for (i in 1:num_sim) {
    scaled_resid <- calc_scaled_residual(noise_effect_size = noise_effect_size, n = nrow(D))
    sim_data[, i] <- qtl_predictor + strain_predictor + kinship_predictor + scaled_resid
  }
  colnames(sim_data) <- paste(sim_label, 1:ncol(sim_data), sep = "_")
  rownames(sim_data) <- strains

  return(list(data = sim_data, 
              properties = list(qtl_effect_size = qtl_effect_size, 
                                strain_effect_size = strain_effect_size, 
                                kinship_effect_size = kinship_effect_size,
                                num_alleles = num_alleles, 
                                sample_method = sample_method)))
}

sim_M_and_beta <- function (num_alleles = 8, num_founders = 8, M = NULL, 
                            sample_method = c("uniform", "crp"), beta = NULL, ...) {
  
  sample_method <- sample_method[1]
  if (is.null(M)) {
    M <- matrix(0, num_founders, num_alleles)
    M[cbind(sample(1:num_founders, num_alleles), 1:num_alleles)] <- 1
    if (sample_method == "uniform") {
      M[which(rowSums(M) == 0), ] <- t(rmultinom(num_founders - num_alleles, 1, rep(1/num_alleles, num_alleles)))
    }
    else if (sample_method == "crp") {
      for (i in which(rowSums(M) == 0)) {
        M[i, ] <- rmultinom(1, 1, colSums(M)/sum(M))
      }
    }
  }
  if (is.null(beta)) {
    beta <- rnorm(num_alleles)
  }
  effect <- list(M = M, beta = beta)
  effect
}

calc_scaled_residual <- function (noise_effect_size, n) {
  
  residual <- rnorm(n = n)
  residual <- (residual - mean(residual))/sqrt(non_sample_var(residual))
  residual <- residual * sqrt(noise_effect_size)
  residual
}

model_matrix_from_ID <- function (M_ID) {
  
  m <- as.numeric(unlist(strsplit(M_ID, ","))) + 1
  J <- length(m)
  M <- matrix(0, J, max(m))
  M[cbind(1:J, m)] <- 1
  M
}

#' Run scans of parametric bootstrap samples to estimate a positional confidence interval for a QTL
#'
#' This function ...
#' 
#' @export
#' @examples run_parboot()
run_parboot <- function(qtl_table, sim_data, genoprobs, 
                        K, map, num_samples = 100,
                        prob = seq(0.8, 0.99, by = 0.01)) {
  
  results_dat <- NULL
  for (i in 1:nrow(qtl_table)) {
    this_locus <- genoprobs[[qtl_table$chr[i]]][rownames(sim_data),,qtl_table$locus[i]]
    this_K <- K[[qtl_table$chr[i]]][rownames(sim_data), rownames(sim_data)]
    
    locus_fit <- fit1(pheno = sim_data[,qtl_table$lodcolumn[i]], 
                      kinship = this_K, 
                      genoprobs = this_locus)
      
    sigma2 <- mean(locus_fit$resid^2) * (nrow(sim_data)/(nrow(sim_data) - 7))
    Xb <- matrix(locus_fit$fitted, ncol = 1)
    Xb <- Xb[,rep(1, num_samples)]
    e <- sapply(1:num_samples, function(i) rnorm(n = nrow(sim_data), mean = 0, sd = sqrt(sigma2)))
    
    y_bs <- Xb + e
    rownames(y_bs) <- rownames(sim_data)
    colnames(y_bs) <- paste("bs_y", 1:num_samples, sep = "_")

    reduced_genoprobs <- list(genoprobs[[qtl_table$chr[i]]])
    names(reduced_genoprobs) <- qtl_table$chr[i]
    attr(reduced_genoprobs, "is_x_chr") <- attr(genoprobs, "is_x_chr")[qtl_table$chr[i]]
    attr(reduced_genoprobs, "crosstype") <- attr(genoprobs, "crosstype")
    attr(reduced_genoprobs, "alleles") <- attr(genoprobs, "alleles")
    attr(reduced_genoprobs, "alleleprobs") <- attr(genoprobs, "alleleprobs")
    attr(reduced_genoprobs, "class") <- attr(genoprobs, "class")
    
    bs_scans <- qtl2::scan1(genoprobs = reduced_genoprobs,
                            pheno = y_bs, 
                            kinship = NULL)
    peak_markers <- apply(bs_scans, 2, function(x) rownames(bs_scans)[which.max(x)])
    peak_pos <- map[[qtl_table$chr[i]]][peak_markers]

    ci_lo <- quantile(peak_pos, prob = (1 - prob)/2)
    ci_hi <- quantile(peak_pos, prob = prob + (1 - prob)/2)
    
    results_dat <- dplyr::bind_rows(results_dat,
                                    data.frame(lodcolumn = qtl_table$lodcolumn[i],
                                               chr = qtl_table$chr[i],
                                               pos = qtl_table$pos[i],
                                               lod = qtl_table$lod[i],
                                               locus = qtl_table$locus[i],
                                               true_locus = qtl_table$true_locus[i],
                                               true_chr = qtl_table$true_chr[i],
                                               true_pos = qtl_table$true_pos[i],
                                               prob = prob,
                                               ci_lo = ci_lo,
                                               ci_hi = ci_hi))
    
    print(paste("Locus", i, "out of", nrow(qtl_table), "done!", collapse = " "))
  }
  
  results_dat <- results_dat %>%
    dplyr::rowwise() %>%
    dplyr::mutate(covered = true_chr == chr & true_pos >= ci_lo & true_pos <= ci_hi,
                  ci_width = ci_hi - ci_lo) %>%
    dplyr::ungroup()
  
  results_dat
}

#' Run scans of Bayesian bootstrap samples to estimate a positional confidence interval for a QTL
#'
#' This function ...
#' 
#' @export
#' @examples run_bayesboot()
run_bayesboot <- function(qtl_table, sim_data, genoprobs, 
                          K, map, num_samples = 100,
                          prob = seq(0.8, 0.99, by = 0.01)) {
  
  ## Sampling weights
  n <- nrow(sim_data)
  weight_matrix <- matrix(NA, nrow = n, ncol = num_samples)
  for(i in 1:num_samples){
    random_num <- sort(runif(n - 1, min = 0, max = 1))
    random_num[n] <- 1
    weights <- c(random_num[1], diff(random_num))
    weight_matrix[,i] <- weights * n
  }
  rownames(weight_matrix) <- rownames(sim_data)
  
  results_dat <- NULL
  for (i in 1:nrow(qtl_table)) {
    this_locus <- genoprobs[[qtl_table$chr[i]]][rownames(sim_data),,qtl_table$locus[i]]
    this_K <- K[[qtl_table$chr[i]]][rownames(sim_data), rownames(sim_data)]
    
    reduced_genoprobs <- list(genoprobs[[qtl_table$chr[i]]])
    names(reduced_genoprobs) <- qtl_table$chr[i]
    attr(reduced_genoprobs, "is_x_chr") <- attr(genoprobs, "is_x_chr")[qtl_table$chr[i]]
    attr(reduced_genoprobs, "crosstype") <- attr(genoprobs, "crosstype")
    attr(reduced_genoprobs, "alleles") <- attr(genoprobs, "alleles")
    attr(reduced_genoprobs, "alleleprobs") <- attr(genoprobs, "alleleprobs")
    attr(reduced_genoprobs, "class") <- attr(genoprobs, "class")
    
    scans <- matrix(NA, nrow = dim(reduced_genoprobs[[1]])[3], ncol = num_samples)
    for (j in 1:num_samples) {
      scans[,j] <- qtl2::scan1(genoprobs = reduced_genoprobs,
                               pheno = sim_data[, i, drop = FALSE], 
                               kinship = NULL, ## weights argument does not seem to work when kinship is included
                               weights = weight_matrix[,j])
    }
    rownames(scans) <- dimnames(reduced_genoprobs[[1]])[[3]]
    
    peak_markers <- apply(scans, 2, function(x) rownames(scans)[which.max(x)])
    peak_pos <- map[[qtl_table$chr[i]]][peak_markers]
    
    ci_lo <- quantile(peak_pos, prob = (1 - prob)/2)
    ci_hi <- quantile(peak_pos, prob = prob + (1 - prob)/2)
    
    results_dat <- dplyr::bind_rows(results_dat,
                                    data.frame(lodcolumn = qtl_table$lodcolumn[i],
                                               chr = qtl_table$chr[i],
                                               pos = qtl_table$pos[i],
                                               lod = qtl_table$lod[i],
                                               locus = qtl_table$locus[i],
                                               true_locus = qtl_table$true_locus[i],
                                               true_chr = qtl_table$true_chr[i],
                                               true_pos = qtl_table$true_pos[i],
                                               prob = prob,
                                               ci_lo = ci_lo,
                                               ci_hi = ci_hi))
    
    print(paste("Locus", i, "out of", nrow(qtl_table), "done!", collapse = " "))
  }
  
  results_dat <- results_dat %>%
    dplyr::rowwise() %>%
    dplyr::mutate(covered = true_chr == chr & true_pos >= ci_lo & true_pos <= ci_hi,
                  ci_width = ci_hi - ci_lo) %>%
    dplyr::ungroup()
  
  results_dat
}

#' Run scans of parametric permutation samples to estimate a positional confidence interval for a QTL
#'
#' This function ...
#' 
#' @export
#' @examples run_parperm()
run_parperm <- function(qtl_table, sim_data, genoprobs, 
                        K, map, num_samples = 100,
                        prob = seq(0.8, 0.99, by = 0.01)) {
  
  results_dat <- NULL
  for (i in 1:nrow(qtl_table)) {
    this_locus <- genoprobs[[qtl_table$chr[i]]][rownames(sim_data),,qtl_table$locus[i]]
    this_K <- K[[qtl_table$chr[i]]][rownames(sim_data), rownames(sim_data)]
    
    locus_fit <- fit1(pheno = sim_data[,qtl_table$lodcolumn[i]], 
                      kinship = this_K, 
                      genoprobs = this_locus)
    
    Xb <- matrix(locus_fit$fitted, ncol = 1)
    Xb <- Xb[,rep(1, num_samples)]
    e <- sapply(1:num_samples, function(i) sample(locus_fit$resid))
    
    y_perm <- Xb + e
    rownames(y_perm) <- rownames(sim_data)
    colnames(y_perm) <- paste("perm_y", 1:num_samples, sep = "_")
    
    reduced_genoprobs <- list(genoprobs[[qtl_table$chr[i]]])
    names(reduced_genoprobs) <- qtl_table$chr[i]
    attr(reduced_genoprobs, "is_x_chr") <- attr(genoprobs, "is_x_chr")[qtl_table$chr[i]]
    attr(reduced_genoprobs, "crosstype") <- attr(genoprobs, "crosstype")
    attr(reduced_genoprobs, "alleles") <- attr(genoprobs, "alleles")
    attr(reduced_genoprobs, "alleleprobs") <- attr(genoprobs, "alleleprobs")
    attr(reduced_genoprobs, "class") <- attr(genoprobs, "class")
    
    perm_scans <- qtl2::scan1(genoprobs = reduced_genoprobs,
                              pheno = y_perm, 
                              kinship = NULL)
    peak_markers <- apply(perm_scans, 2, function(x) rownames(perm_scans)[which.max(x)])
    peak_pos <- map[[qtl_table$chr[i]]][peak_markers]
    
    ci_lo <- quantile(peak_pos, prob = (1 - prob)/2)
    ci_hi <- quantile(peak_pos, prob = prob + (1 - prob)/2)
    
    results_dat <- dplyr::bind_rows(results_dat,
                                    data.frame(lodcolumn = qtl_table$lodcolumn[i],
                                               chr = qtl_table$chr[i],
                                               pos = qtl_table$pos[i],
                                               lod = qtl_table$lod[i],
                                               locus = qtl_table$locus[i],
                                               true_locus = qtl_table$true_locus[i],
                                               true_chr = qtl_table$true_chr[i],
                                               true_pos = qtl_table$true_pos[i],
                                               prob = prob,
                                               ci_lo = ci_lo,
                                               ci_hi = ci_hi))
    
    print(paste("Locus", i, "out of", nrow(qtl_table), "done!", collapse = " "))
  }
  
  results_dat <- results_dat %>%
    dplyr::rowwise() %>%
    dplyr::mutate(covered = true_chr == chr & true_pos >= ci_lo & true_pos <= ci_hi,
                  ci_width = ci_hi - ci_lo) %>%
    dplyr::ungroup()
  
  results_dat
}

#' Run scans of parametric permutation samples (with the kinship matrix) to estimate a positional confidence interval for a QTL
#'
#' This function ...
#' 
#' @export
#' @examples run_parperm_kinship()
run_parperm_kinship <- function(qtl_table, sim_data, genoprobs, 
                                K, map, num_samples = 100,
                                prob = seq(0.8, 0.99, by = 0.01)) {
  
  n <- nrow(sim_data)
  e_index <- sapply(1:num_samples, function(i) sample(1:n))
  results_dat <- NULL
  for (i in 1:nrow(qtl_table)) {
    this_locus <- genoprobs[[qtl_table$chr[i]]][rownames(sim_data),,qtl_table$locus[i]]
    this_K <- K[[qtl_table$chr[i]]][rownames(sim_data), rownames(sim_data)]
    
    locus_fit <- fit1(pheno = sim_data[,qtl_table$lodcolumn[i]], 
                      kinship = this_K, 
                      genoprobs = this_locus)
    
    Xb <- matrix(locus_fit$fitted, ncol = 1)
    Xb <- Xb[,rep(1, num_samples)]
    e <- sapply(1:num_samples, function(i) locus_fit$resid[e_index[,i]])
    
    y_perm <- Xb + e
    rownames(y_perm) <- rownames(sim_data)
    colnames(y_perm) <- paste("perm_y", 1:num_samples, sep = "_")
    
    reduced_genoprobs <- list(genoprobs[[qtl_table$chr[i]]])
    names(reduced_genoprobs) <- qtl_table$chr[i]
    attr(reduced_genoprobs, "is_x_chr") <- attr(genoprobs, "is_x_chr")[qtl_table$chr[i]]
    attr(reduced_genoprobs, "crosstype") <- attr(genoprobs, "crosstype")
    attr(reduced_genoprobs, "alleles") <- attr(genoprobs, "alleles")
    attr(reduced_genoprobs, "alleleprobs") <- attr(genoprobs, "alleleprobs")
    attr(reduced_genoprobs, "class") <- attr(genoprobs, "class")
    
    perm_scans <- matrix(NA, nrow = dim(reduced_genoprobs[[1]])[3], ncol = num_samples)
    for (j in 1:num_samples) {
      perm_K <- this_K[rownames(this_K)[e_index[,j]], rownames(this_K)[e_index[,j]]]
      rownames(perm_K) <- colnames(perm_K) <- rownames(this_K)
      perm_scans[,j] <- qtl2::scan1(genoprobs = reduced_genoprobs,
                                    pheno = y_perm[, j, drop = FALSE], 
                                    kinship = perm_K)
    }
    rownames(perm_scans) <- dimnames(reduced_genoprobs[[1]])[[3]]
    
    peak_markers <- apply(perm_scans, 2, function(x) rownames(perm_scans)[which.max(x)])
    peak_pos <- map[[qtl_table$chr[i]]][peak_markers]
    
    ci_lo <- quantile(peak_pos, prob = (1 - prob)/2)
    ci_hi <- quantile(peak_pos, prob = prob + (1 - prob)/2)
    
    results_dat <- dplyr::bind_rows(results_dat,
                                    data.frame(lodcolumn = qtl_table$lodcolumn[i],
                                               chr = qtl_table$chr[i],
                                               pos = qtl_table$pos[i],
                                               lod = qtl_table$lod[i],
                                               locus = qtl_table$locus[i],
                                               true_locus = qtl_table$true_locus[i],
                                               true_chr = qtl_table$true_chr[i],
                                               true_pos = qtl_table$true_pos[i],
                                               prob = prob,
                                               ci_lo = ci_lo,
                                               ci_hi = ci_hi))
    
    print(paste("Locus", i, "out of", nrow(qtl_table), "done!", collapse = " "))
  }
  
  results_dat <- results_dat %>%
    dplyr::rowwise() %>%
    dplyr::mutate(covered = true_chr == chr & true_pos >= ci_lo & true_pos <= ci_hi,
                  ci_width = ci_hi - ci_lo) %>%
    dplyr::ungroup()
  
  results_dat
}

#' Run scans of subsamples of data to estimate a positional confidence interval for a QTL
#'
#' This function ...
#' 
#' @export
#' @examples run_subsample()
run_subsample <- function(qtl_table, sim_data, genoprobs, 
                          K, map, num_draws = 100,
                          sample_perc = 2/3, 
                          prob = seq(0.8, 0.99, by = 0.01)) {
  
  results_dat <- NULL
  for (i in 1:nrow(qtl_table)) {
    this_locus <- genoprobs[[qtl_table$chr[i]]][rownames(sim_data),,qtl_table$locus[i]]
    this_K <- K[[qtl_table$chr[i]]][rownames(sim_data), rownames(sim_data)]
  
    reduced_genoprobs <- list(genoprobs[[qtl_table$chr[i]]])
    names(reduced_genoprobs) <- qtl_table$chr[i]
    attr(reduced_genoprobs, "is_x_chr") <- attr(genoprobs, "is_x_chr")[qtl_table$chr[i]]
    attr(reduced_genoprobs, "crosstype") <- attr(genoprobs, "crosstype")
    attr(reduced_genoprobs, "alleles") <- attr(genoprobs, "alleles")
    attr(reduced_genoprobs, "alleleprobs") <- attr(genoprobs, "alleleprobs")
    attr(reduced_genoprobs, "class") <- attr(genoprobs, "class")
    
    scans <- matrix(NA, nrow = dim(reduced_genoprobs[[1]])[3], ncol = num_draws)
    for (j in 1:num_draws) {
      this_sample <- sample(rownames(this_K), size = ceiling(sample_perc*nrow(sim_data)))
      scans[,j] <- qtl2::scan1(genoprobs = reduced_genoprobs,
                               pheno = sim_data[this_sample, i, drop = FALSE], 
                               kinship = this_K)
    }
    rownames(scans) <- dimnames(reduced_genoprobs[[1]])[[3]]
    
    peak_markers <- apply(scans, 2, function(x) rownames(scans)[which.max(x)])
    peak_pos <- map[[qtl_table$chr[i]]][peak_markers]
    
    ci_lo <- quantile(peak_pos, prob = (1 - prob)/2)
    ci_hi <- quantile(peak_pos, prob = prob + (1 - prob)/2)
    
    results_dat <- dplyr::bind_rows(results_dat,
                                    data.frame(lodcolumn = qtl_table$lodcolumn[i],
                                               chr = qtl_table$chr[i],
                                               pos = qtl_table$pos[i],
                                               lod = qtl_table$lod[i],
                                               locus = qtl_table$locus[i],
                                               true_locus = qtl_table$true_locus[i],
                                               true_chr = qtl_table$true_chr[i],
                                               true_pos = qtl_table$true_pos[i],
                                               prob = prob,
                                               ci_lo = ci_lo,
                                               ci_hi = ci_hi))
    
    print(paste("Locus", i, "out of", nrow(qtl_table), "done!", collapse = " "))
  }
  
  results_dat <- results_dat %>%
    dplyr::rowwise() %>%
    dplyr::mutate(covered = true_chr == chr & true_pos >= ci_lo & true_pos <= ci_hi,
                  ci_width = ci_hi - ci_lo) %>%
    dplyr::ungroup()
  
  results_dat
}

#' Evaluate mapping performance when using LOD drop or Bayesian credible support intervals
#'
#' This function ...
#' 
#' @export
#' @examples eval_mapping_results()
eval_mapping_results <- function(thresh = seq(6, 9, by = 0.5), 
                                 lod_drop = seq(1, 3, by = 0.1),
                                 bci = seq(0.8, 0.95, by = 0.01),
                                 use_lod_drop = TRUE,
                                 scans, map, true_qtl) {
  
  rate_dat <- NULL
  
  if (use_lod_drop) {
    interval <- lod_drop
  } else {
    interval <- bci
  }
  for (i in 1:length(interval)) {
    if (use_lod_drop) {
      lodpeaks <- qtl2::find_peaks(scans, map = map, threshold = min(thresh), drop = interval[i])
    } else {
      lodpeaks <- qtl2::find_peaks(scans, map = map, threshold = min(thresh), prob = interval[i])
    }
    
    for (j in 1:length(thresh)) {
      detected_qtl <- lodpeaks %>%
        dplyr::filter(lod > thresh[j]) %>%
        dplyr::left_join(true_qtl %>%
                           dplyr::rename(true_chr = chr,
                                         true_pos = pos)) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(ci_width = ci_hi - ci_lo,
                      covered = true_chr == chr & true_pos >= ci_lo & true_pos <= ci_hi) %>%
        dplyr::ungroup()
      
      detection_rate <- sum(detected_qtl$covered)/nrow(true_qtl)
      true_positive_rate <- sum(detected_qtl$covered)/nrow(detected_qtl)
      false_positive_rate <- sum(!detected_qtl$covered)/nrow(detected_qtl)
      
      reduced_qtl <- detected_qtl %>%
        dplyr::filter(true_chr == chr & abs(true_pos - pos) < 20)
      
      coverage_rate <- sum(reduced_qtl$covered)/nrow(reduced_qtl)
      mean_ci_width <- mean(reduced_qtl$ci_width)
      
      if (use_lod_drop) {
        rate_dat <- bind_rows(rate_dat,
                              data.frame(thresh = thresh[j], lod_drop = lod_drop[i],
                                         detection_rate = detection_rate, true_positive_rate = true_positive_rate,
                                         false_positive_rate = false_positive_rate, coverage_rate = coverage_rate,
                                         mean_ci_width = mean_ci_width))
      } else {
        rate_dat <- bind_rows(rate_dat,
                              data.frame(thresh = thresh[j], bci = bci[i],
                                         detection_rate = detection_rate, true_positive_rate = true_positive_rate,
                                         false_positive_rate = false_positive_rate, coverage_rate = coverage_rate,
                                         mean_ci_width = mean_ci_width))
      }
    }
  }
  
  rate_dat
}

#' Evaluate mapping performance with null simulations
#'
#' This function ...
#' 
#' @export
#' @examples eval_null_mapping_results()
eval_null_mapping_results <- function(thresh = seq(6, 10, by = 0.5), scans, map) {
  
  count_dat <- NULL
  
  lodpeaks <- qtl2::find_peaks(scans, map = map, threshold = min(thresh))
  
  for (i in 1:length(thresh)) {
    count_dat <- dplyr::bind_rows(count_dat,
                                  data.frame(thresh = thresh[i], 
                                             count = lodpeaks %>%
                                               dplyr::filter(lod > thresh[i]) %>%
                                               nrow))
  }
  
  count_dat
}

#' Evaluate mapping performance when using parametric bootstrap performance
#'
#' This function ...
#' 
#' @export
#' @examples eval_mapping_parboot_results()
eval_mapping_parboot_results <- function(thresh = seq(6, 9, by = 0.5), 
                                         prob = seq(0.8, 0.95, by = 0.01),
                                         parboot_results, 
                                         num_true_qtl) {
  
  rate_dat <- NULL
  
  for (i in 1:length(thresh)) {
    for (j in 1:length(prob)) {
      detected_qtl <- parboot_results %>%
        dplyr::filter(lod > thresh[i],
                      prob == prob[j])
      detection_rate <- sum(detected_qtl$covered)/num_true_qtl
      true_positive_rate <- sum(detected_qtl$covered)/nrow(detected_qtl)
      false_positive_rate <- sum(!detected_qtl$covered)/nrow(detected_qtl)
      mean_ci_width = mean(detected_qtl$ci_width)
      
      reduced_qtl <- detected_qtl %>%
        dplyr::filter(true_chr == chr & abs(true_pos - pos) < 20)
      
      coverage_rate <- sum(reduced_qtl$covered)/nrow(reduced_qtl)

      rate_dat <- dplyr::bind_rows(rate_dat,
                                   data.frame(thresh = thresh[i], prob = prob[j],
                                              detection_rate = detection_rate, 
                                              true_positive_rate = true_positive_rate,
                                              false_positive_rate = false_positive_rate, 
                                              coverage_rate = coverage_rate,
                                              mean_ci_width = mean_ci_width))
    }
  }
  
  rate_dat
}

#' Pull FWER thresholds from maximum LOD scores
#'
#' This function ...
#' 
#' @export
#' @examples pull_fwer_thresh()
pull_fwer_thresh <- function(maxlod, 
                             right_side = TRUE, 
                             fwer = 0.05,
                             use_gev = TRUE) {
  
  if (right_side) {
    fwer <- 1 - fwer
  }
  if (use_gev) {
    ## Fit GEV from maximum LODs of permutations
    evd_pars <- evd::fgev(maxlod)$estimate
    ## Pull threshold as quantile from GEV
    thresh <- evd::qgev(p = fwer, 
                        loc = evd_pars["loc"], 
                        scale = evd_pars["scale"], 
                        shape = evd_pars["shape"])
  } else {
    ## Pull empirical p-values from permutation maximum LOD scores
    thresh <- quantile(maxlod, p = fwer)
  }
  
  as.numeric(thresh)
}

