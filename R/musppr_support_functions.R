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

non_sample_var <- function (x)  {
  
  var_x <- var(x) * ((length(x) - 1)/length(x))
  var_x
}

incidence_matrix <- function (fact) {
  
  m <- diag(nlevels(fact))[fact, ]
  colnames(m) <- levels(fact)
  return(m)
}
