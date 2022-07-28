#' Simulate data and evaluate heritability estimation
#'
#' This function ...
#' 
#' @export
#' @examples eval_sim_h2()
eval_sim_h2 <- function(sim_h2, n, n_per_genome = 1,
                        K, K_fit = NULL,
                        intercept = 0,
                        method = c("qtl2", "sommer", "miqtl"),
                        use_rint = FALSE) {
  
  if (is.null(K_fit)) { K_fit <- K }
  method <- method[1]
  
  ## Simulation
  u <- t(MASS::mvrnorm(n = n, mu = rep(0, nrow(K)), Sigma = K))
  e <- sapply(1:n, function(i) rnorm(n = nrow(K)))
  
  ## Adjust effect sizes for replicates
  if (sim_h2 > 0) {
    noise_effect_size <- 1 - sim_h2
    noise_effect_size <- noise_effect_size/n_per_genome
    scale_factor <- (1 - noise_effect_size)/sim_h2
    sim_h2 <- scale_factor * sim_h2
  }
  
  y <- sapply(1:n, function(i) intercept + u[,i] * sqrt(sim_h2/non_sample_var(u[,i])) + e[,i] * sqrt((1 - sim_h2)/non_sample_var(e[,i])))
  if (use_rint) {
    y <- apply(y, 2, function(x) rint(x))
  }
  rownames(y) <- rownames(K)
  
  if (method == "qtl2") {
    h2 <- qtl2::est_herit(pheno = y, kinship = K_fit)
  } else if (method == "sommer") {
    fit_sommer <- function(sim_y, i, K) {
      
      fit <- sommer::mmer(y ~ 1,
                          random = ~ sommer::vs(id, Gu = K),
                          rcov = ~ units, 
                          data = data.frame(y = sim_y[,i], id = rownames(sim_y)),
                          verbose = FALSE,
                          date.warning = FALSE)
      h2 <- sommer::vpredict(fit, h2 ~ (V1)/(V1 + V2))
      h2$Estimate
    }
    
    h2 <- sapply(1:n, function(i) fit_sommer(sim_y = y, i = i, K = K_fit))
  } else if (method == "miqtl") {
    eigen.K <- eigen(K_fit)
    fit_miqtl <- function(sim_y, i, K, eigen.K) {
      
      fit <- miqtl::lmmbygls(y ~ 1, 
                             data = data.frame(y = sim_y[,i], SUBJECT.NAME = rownames(sim_y)),
                             K = K, eigen.K = eigen.K, use.par = "h2.REML")
      fit$h2
    }
    h2 <- sapply(1:n, function(i) fit_miqtl(sim_y = y, i = i, K = K_fit, eigen.K = eigen.K))
  }
  
  h2
}

#' Simulate data with strain/genome replicates and evaluate heritability estimation
#'
#' This function ...
#' 
#' @export
#' @examples eval_sim_h2_with_reps()
eval_sim_h2_with_reps <- function(sim_h2, n_sims, n_per_strain, K_strains, K_strains_fit = NULL,
                                  intercept = 0, method = c("sommer", "miqtl"),
                                  use_rint = FALSE, do_Kasv_fit = FALSE) {
  
  if (is.null(K_strains_fit)) { K_strains_fit <- K_strains }
  if (do_Kasv_fit & method == "sommer") {
    cat("Cannot currently use the full ASV form of the kinship matrix (mouse-level matrix with replicates) with sommer")
    cat("do_Kasv_fit essentially set to FALSE")
  }
  method <- method[1]
  
  u <- t(MASS::mvrnorm(n = n_sims, mu = rep(0, nrow(K_strains)), Sigma = K_strains))
  u <- u[rep(1:nrow(K_strains), each = n_per_strain),]
  e <- sapply(1:n_sims, function(i) rnorm(n = nrow(K_strains) * n_per_strain))
  
  y <- sapply(1:n_sims, function(i) intercept + u[,i] * sqrt(sim_h2/non_sample_var(u[,i])) + e[,i] * sqrt((1 - sim_h2)/non_sample_var(e[,i])))
  if (use_rint) {
    y <- apply(y, 2, function(x) rint(x))
  }
  rownames(y) <- paste(rep(rownames(K_strains), each = n_per_strain), rep(1:n_per_strain, times = nrow(K_strains)), sep = "_")
  
  if (method == "sommer") {
    fit_sommer <- function(sim_y, i, K) {
      
      fit <- sommer::mmer(y ~ 1,
                          random = ~ sommer::vs(strain, Gu = K),
                          rcov = ~ units, 
                          data = data.frame(y = sim_y[,i], 
                                            id = rownames(sim_y), 
                                            strain = gsub(x = rownames(sim_y), 
                                                          pattern = "_[0-9]+", 
                                                          replacement = "")),
                          verbose = FALSE,
                          date.warning = FALSE)
      
      h2 <- sommer::vpredict(fit, h2 ~ (V1)/(V1 + V2))
      h2$Estimate
    }
    
    h2 <- sapply(1:n_sims, function(i) fit_sommer(sim_y = y, i = i, K = K_strains_fit))
  } else if (method == "miqtl") {
    ind_to_strain_data <- data.frame(SUBJECT.NAME = rownames(y), 
                                     strain = gsub(x = rownames(y), 
                                                   pattern = "_[0-9]+", 
                                                   replacement = "")) %>%
      mutate(strain = factor(strain, levels = colnames(K_strains_fit)))
    Z <- model.matrix(miqtl:::process.random.formula(geno.id = "strain"), 
                      data = ind_to_strain_data)
    K_ind <- Z %*% tcrossprod(K_strains_fit, Z)
    rownames(K_ind) <- colnames(K_ind) <- as.character(ind_to_strain_data[,"SUBJECT.NAME"])
    if (do_Kasv_fit) {
      K_ind <- make_Kasv(K_ind)
    }
    
    eigen.K <- eigen(K_ind)
    
    fit_miqtl <- function(sim_y, i, K, eigen.K) {
      
      fit <- miqtl::lmmbygls(y ~ 1, 
                             data = data.frame(y = sim_y[,i], SUBJECT.NAME = rownames(sim_y)),
                             K = K, eigen.K = eigen.K, use.par = "h2.REML")
      fit$h2
    }
    h2 <- sapply(1:n_sims, function(i) fit_miqtl(sim_y = y, i = i, K = K_ind, eigen.K = eigen.K))
  }

  h2
}

#' Simulate data with additive and strain/genome heritability components and evaluate heritability estimation
#'
#' This function ...
#' 
#' @export
#' @examples eval_sim_h2_sommer_strainvar()
eval_sim_h2_sommer_strainvar <- function(sim_h2_add_prop, h2_total,
                                         n_sims, n_per_strain, K_strains,  K_strains_fit = NULL,
                                         intercept = 0,
                                         use_rint = FALSE) {
  
  if (is.null(K_strains_fit)) { K_strains_fit <- K_strains }
  
  u_add <- t(MASS::mvrnorm(n = n_sims, mu = rep(0, nrow(K_strains)), Sigma = K_strains))
  u_add <- u_add[rep(1:nrow(K_strains), each = n_per_strain),]
  u_strain <- sapply(1:n_sims, function(i) rnorm(n = nrow(K_strains)))
  u_strain <- u_strain[rep(1:nrow(K_strains), each = n_per_strain),]
  e <- sapply(1:n_sims, function(i) rnorm(n = nrow(K_strains) * n_per_strain))
  
  y <- sapply(1:n_sims, function(i) intercept + u_add[,i] * sqrt((sim_h2_add_prop * h2_total)/non_sample_var(u_add[,i])) + u_strain[,i] * sqrt(((1 - sim_h2_add_prop) * h2_total)/non_sample_var(u_strain[,i])) + e[,i] * sqrt((1 - h2_total)/non_sample_var(e[,i])))
  if (use_rint) {
    y <- apply(y, 2, function(x) rint(x))
  }
  rownames(y) <- paste(rep(rownames(K_strains), each = n_per_strain), rep(1:n_per_strain, times = nrow(K_strains)), sep = "_")
  
  fit_sommer <- function(sim_y, i, K) {
    
    fit <- sommer::mmer(y ~ 1,
                        random = ~ sommer::vs(strain, Gu = K) + strain,
                        rcov = ~ units, 
                        data = data.frame(y = sim_y[,i], 
                                          id = rownames(sim_y), 
                                          strain = gsub(x = rownames(sim_y), 
                                                        pattern = "_[0-9]+", 
                                                        replacement = "")),
                        verbose = FALSE,
                        date.warning = FALSE)
    h2_add <- sommer::vpredict(fit, h2_add ~ (V1)/(V1 + V2 + V3))
    h2_strain <- sommer::vpredict(fit, h2_strain ~ (V2)/(V1 + V2 + V3))
    c(h2_add$Estimate, h2_strain$Estimate)
  }
  
  h2 <- t(sapply(1:n_sims, function(i) fit_sommer(sim_y = y, i = i, K = K_strains_fit), simplify = "array"))
  h2
}

#' Scale kinship matrix to average semivariance 
#'
#' This function ...
#' 
#' @export
#' @examples make_Kasv()
make_Kasv <- function(Z) {
  
  Kbar <- scale(Z, scale = FALSE) %*% t(scale(Z, scale = FALSE))
  Kasv <- Kbar/(psych::tr(Kbar)/(nrow(Z) - 1))
  rownames(Kasv) <- colnames(Kasv) <- rownames(Z)
  Kasv
}

#' Simulate data with additive and strain/genome heritability components and evaluate heritability estimation in an additive-only model
#'
#' This function ...
#' 
#' @export
#' @examples eval_sim_h2_miqtl_strainvar_add_only_fit()
eval_sim_h2_miqtl_strainvar_add_only_fit <- function(sim_h2_add_prop, h2_total,
                                                     n_sims, n_per_strain, K_strains,  K_strains_fit = NULL,
                                                     intercept = 0,
                                                     use_rint = FALSE,
                                                     do_Kasv_fit = FALSE) {
  
  if (is.null(K_strains_fit)) { K_strains_fit <- K_strains }
  
  u_add <- t(MASS::mvrnorm(n = n_sims, mu = rep(0, nrow(K_strains)), Sigma = K_strains))
  u_add <- u_add[rep(1:nrow(K_strains), each = n_per_strain),]
  u_strain <- sapply(1:n_sims, function(i) rnorm(n = nrow(K_strains)))
  u_strain <- u_strain[rep(1:nrow(K_strains), each = n_per_strain),]
  e <- sapply(1:n_sims, function(i) rnorm(n = nrow(K_strains) * n_per_strain))
  
  y <- sapply(1:n_sims, function(i) intercept + u_add[,i] * sqrt((sim_h2_add_prop * h2_total)/non_sample_var(u_add[,i])) + u_strain[,i] * sqrt(((1 - sim_h2_add_prop) * h2_total)/non_sample_var(u_strain[,i])) + e[,i] * sqrt((1 - h2_total)/non_sample_var(e[,i])))
  if (use_rint) {
    y <- apply(y, 2, function(x) rint(x))
  }
  rownames(y) <- paste(rep(rownames(K_strains), each = n_per_strain), rep(1:n_per_strain, times = nrow(K_strains)), sep = "_")
  
  ind_to_strain_data <- data.frame(SUBJECT.NAME = rownames(y), 
                                   strain = gsub(x = rownames(y), 
                                                 pattern = "_[0-9]+", 
                                                 replacement = "")) %>%
    mutate(strain = factor(strain, levels = colnames(K_strains_fit)))
  Z <- model.matrix(miqtl:::process.random.formula(geno.id = "strain"), 
                    data = ind_to_strain_data)
  K_ind <- Z %*% tcrossprod(K_strains_fit, Z)
  rownames(K_ind) <- colnames(K_ind) <- as.character(ind_to_strain_data[,"SUBJECT.NAME"])
  
  if (do_Kasv_fit) {
    K_ind <- make_Kasv(K_ind)
  }
  
  eigen.K <- eigen(K_ind)
  
  fit_miqtl <- function(sim_y, i, K, eigen.K) {
    
    fit <- miqtl::lmmbygls(y ~ 1, 
                           data = data.frame(y = sim_y[,i], SUBJECT.NAME = rownames(sim_y)),
                           K = K, eigen.K = eigen.K, use.par = "h2.REML")
    fit$h2
  }
  h2 <- sapply(1:n_sims, function(i) fit_miqtl(sim_y = y, i = i, K = K_ind, eigen.K = eigen.K))
  h2
}


