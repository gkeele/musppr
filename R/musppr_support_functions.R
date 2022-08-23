#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

# Conversion function taken from qtl2convert
map_list_to_df <- function (map_list, chr_column = "chr", pos_column = "pos", marker_column = "marker") {
  
  nmar <- vapply(map_list, length, 1)
  markers <- unlist(lapply(map_list, names))
  result <- data.frame(chr = rep(names(map_list), nmar), pos = unlist(map_list), 
                       marker = markers, stringsAsFactors = FALSE)
  rownames(result) <- markers
  names(result)[1] <- chr_column
  names(result)[2] <- pos_column
  if (is.null(marker_column)) 
    result <- result[, -3, drop = FALSE]
  else names(result)[3] <- marker_column
  result
}

#' Scale kinship matrix to average semivariance 
#'
#' This function takes a kinship matrix and calculates the average semivariance form.
#' 
#' @param K DEFAULT: NULL. Haplotype-based kinship matrix, as calculated with \code{qtl2::calc_kinship()}.
#' @param Z DEFAULT: NULL. SNP allele count design matrix. Provide Z if using a SNP-based kinship matrix.
#' @export
#' @examples make_Kasv()
make_Kasv <- function(K = NULL,
                      Z = NULL) {
  
  if (!is.null(Z)) {
    Kbar <- scale(Z, scale = FALSE) %*% t(scale(Z, scale = FALSE))
  } else if (!is.null(K)) {
    P <- diag(nrow(K)) - (1/nrow(K))*matrix(1, nrow = nrow(K), ncol = 1) %*% matrix(1, nrow = 1, ncol = nrow(K))
    Kbar <- P %*% K %*% t(P)
  }
  Kasv <- Kbar/(psych::tr(Kbar)/(nrow(Kbar) - 1))
  
  if (!is.null(Z)) {
    rownames(Kasv) <- colnames(Kasv) <- rownames(Z)
  } else if (!is.null(K)) {
    rownames(Kasv) <- colnames(Kasv) <- rownames(K)
  }
  
  Kasv
}

non_sample_var <- function (x)  {
  
  var_x <- var(x) * ((length(x) - 1)/length(x))
  var_x
}

incidence_matrix <- function (fact) {
  
  m <- diag(nlevels(fact))[fact, ]
  colnames(m) <- levels(fact)
  return(m)
}

rint <- function (x) {
  x = rank(x, na.last = "keep", ties.method = "average")/(sum(!is.na(x)) + 1)
  qnorm(x)
}
