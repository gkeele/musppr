#' CC-RIX grid plot
#'
#' This function ...
#' 
#' @export
#' @examples plot_ccrix_grid()
plot_ccrix_grid <- function(cc_strains,
                            ccrix_set,
                            low = "white",
                            high = "black", 
                            include_cc = TRUE,
                            grid_line_width = 0.5,
                            text_size = 6) {
  
  ccrix_grid <- matrix(0, nrow = length(cc_strains), ncol = length(cc_strains))
  rownames(ccrix_grid) <- colnames(ccrix_grid) <- cc_strains
  
  for (i in 1:length(ccrix_set)) {
    these_strains <- unlist(strsplit(ccrix_set[i], split = "x"))
    ccrix_grid[these_strains[1], these_strains[2]] <- ccrix_grid[these_strains[2], these_strains[1]] <- 1
  }
  if (include_cc) {
    for (i in 1:length(cc_strains)) {
      ccrix_grid[i, i] <- 1
    }
  }
  
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat, diag = TRUE)] <- NA
    return(cormat)
  }
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  ccrix_grid_long_dat <- reshape2::melt(get_upper_tri(ccrix_grid), na.rm = TRUE)
  blank_dat <- reshape2::melt(get_lower_tri(matrix(1, nrow = length(cc_strains), ncol = length(cc_strains))), na.rm = TRUE)
  
  g <- ggplot2::ggplot(data = ccrix_grid_long_dat, 
                       ggplot2::aes(x = Var1, y = Var2, fill = value)) + 
    ggplot2::geom_tile(color = "gray", size = grid_line_width) + 
    ggplot2::geom_tile(data = blank_dat, fill = "white", size = 0.1) + 
    ggplot2::xlab("") + 
    ggplot2::ylab("") +
    ggplot2::scale_fill_gradient(low = low, high = high, limit = c(0, 1)) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::guides(fill = "none") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x.top = ggplot2::element_text(angle = 90, size = text_size, vjust = 0.5, hjust = 1), 
                   axis.text.y = ggplot2::element_text(size = text_size),
                   plot.title = ggplot2::element_text(hjust = 0.5, size = 16),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
    coord_fixed()
  
  g
}

#' Genome scan plot
#'
#' This function ...
#' 
#' @export
#' @examples plot_genome_scan()
plot_genome_scan <- function(qtl_scan,
                             lodcolumn = 1,
                             med_scan = NULL,
                             med_highlight_dat = NULL,
                             map,
                             bgcolor = "white",
                             altbgcolor = "white",
                             color_alpha = 0.6,
                             altcolor_alpha = 0.6,
                             color = "cadetblue",
                             altcolor = "cadetblue",
                             medcol = "gray90",
                             outlinecol = NA,
                             hlines_col = "gray",
                             hlines_lwd = 0.5,
                             added_hline = NULL,
                             annot_dat,
                             target_symbol = NULL,
                             chr = c(1:19, "X"),
                             main = "",
                             my_y_max = NULL,
                             rug_col = "black") {
  
  include_chr <- chr
  if (is.null(med_scan)) {
    y_max <- max(qtl_scan)
  } else {
    med_scan <- med_scan %>% 
      filter(chr %in% chr)
    y_max <- max(c(qtl_scan, med_scan$lod))
  }
  if (!is.null(my_y_max)) {
    y_max <- max(y_max, my_y_max)
  }
  plot(qtl_scan, map, lodcolumn,
       ylim = c(0, y_max), 
       bgcolor = bgcolor, 
       altbgcolor = altbgcolor, 
       col = NA, 
       gap = 0,
       chr = chr,
       main = main,
       hlines_col = hlines_col,
       hlines_lwd = hlines_lwd)
  map_df <- map_list_to_df(map)
  j <- 1
  for (i in chr) {
    x <- xpos_scan1(map = map, thepos = map_df %>% filter(chr == i) %>% pull(pos), thechr = map_df %>% filter(chr == i) %>% pull(chr), gap = 0, chr = chr)
    y <- qtl_scan[map_df %>% filter(chr == i) %>% pull(marker), lodcolumn]
    polygon(x = c(x[1], x, x[length(x)]),
            y = c(0, y, 0), 
            border = NA,
            col = ifelse(j %% 2 != 0, 
                         scales::alpha(color, color_alpha), 
                         scales::alpha(altcolor, altcolor_alpha)))
    j <- j + 1
  }
  if (!is.null(med_scan)) {
    points(x = qtl2::xpos_scan1(map, 
                                thechr = med_scan$chr, 
                                thepos = med_scan$pos, 
                                gap = 0,
                                chr = chr), 
           y = med_scan$lod,
           col = medcol, 
           pch = 20)
  }
  
  plot(qtl_scan, map, lodcolumn,
       ylim = c(0, y_max), 
       bgcolor = NA, 
       altbgcolor = NA,
       hlines = NULL,
       hlines_col = NA,
       hlines_lwd = NA,
       col = outlinecol, 
       lwd = 0.5,
       gap = 0,
       add = TRUE,
       chr = chr)
  if (!is.null(med_highlight_dat)) {
    for (i in 1:nrow(med_highlight_dat)) {
      points(x = qtl2::xpos_scan1(map, 
                                  thechr = med_scan %>% filter(symbol == med_highlight_dat$symbol[i]) %>% pull(chr), 
                                  thepos = med_scan %>% filter(symbol == med_highlight_dat$symbol[i]) %>% pull(pos), 
                                  gap = 0,
                                  chr = chr), 
             y = med_scan %>% filter(symbol == med_highlight_dat$symbol[i]) %>% pull(lod),
             col = med_highlight_dat$col[i], 
             bg = med_highlight_dat$bgcol[i],
             cex = med_highlight_dat$cex[i],
             pch = med_highlight_dat$pch[i])
    }
  }
  if (!is.null(target_symbol)) {
    for (i in 1:length(target_symbol)) {
      rug(x = qtl2::xpos_scan1(map, 
                               thechr = annot_dat %>% filter(symbol == target_symbol[i]) %>% pull(chr),
                               thepos = annot_dat %>% filter(symbol == target_symbol[i]) %>% pull(middle),
                               gap = 0,
                               chr = chr),
          col = rug_col,
          lwd = 3)
    }
  }
  if (!is.null(added_hline)) {
    abline(h = added_hline, lty = 2)
  }
}

pull_half_circle_coord <- function(width, 
                                   x_center = NULL,
                                   y_center = 0,
                                   interval_length = 500) {
  
  rad <- width/2
  if (is.null(x_center)) { x_center <- rad }
  
  x <- seq(x_center - rad, x_center + rad, length.out = interval_length)
  
  y = sqrt(rad^2 - (x - x_center)^2) + y_center
  circle_coord <- list(y = y, x = x)
  circle_coord
}

pull_quarter_circle_coord <- function(width, 
                                      x_center = NULL,
                                      y_center = 0,
                                      interval_length = 500) {
  
  rad <- width/8
  if (is.null(x_center)) { x_center <- rad }
  
  x <- seq(x_center - rad, x_center, length.out = interval_length)
  
  y = sqrt(rad^2 - (x - x_center)^2) + y_center
  circle_coord <- list(y = y, x = x)
  circle_coord
}

plot_chromatid <- function(intervals, 
                           cols,
                           width = 0.1,
                           y_shift = NULL,
                           x_shift = 0,
                           add_border = TRUE,
                           border_lwd = 1.5,
                           end_type = c("rounded", "circle", "square")) {
  
  end_type <- end_type[1]
  if (is.null(y_shift)) { y_shift <- width/2 }
  
  intervals <- intervals/(1/(1 - y_shift*2)) + y_shift
  for (i in 1:(length(intervals) - 1)) {
    polygon(x = c(x_shift, x_shift, x_shift + width, x_shift + width), 
            y = c(intervals[i], intervals[i + 1], intervals[i + 1], intervals[i]), 
            col = cols[i], 
            border = NA)
  }
  if (add_border) {
    segments(x0 = x_shift, x1 = x_shift, y0 = y_shift, y1 = 1 - y_shift, lwd = border_lwd)
    segments(x0 = x_shift + width, x1 = x_shift + width, y0 = y_shift, y1 = 1 - y_shift, lwd = border_lwd)
  }
  if (end_type == "rounded") {
    segments(x0 = x_shift, x1 = x_shift + width, y0 = y_shift, y1 = y_shift, lwd = 1, col = cols[1])
    segments(x0 = x_shift, x1 = x_shift + width, y0 = 1 - y_shift, y1 = 1 - y_shift, lwd = 1, col = cols[length(cols)])
    
    rounded_cap <- pull_quarter_circle_coord(width = width, x_center = width/8, y_center = 0)
    polygon(x = c(x_shift + rounded_cap$x, x_shift + width - rev(rounded_cap$x)), 
            y = c(y_shift - rounded_cap$y, y_shift - rev(rounded_cap$y)),
            col = cols[1], 
            border = NA)
    polygon(x = c(x_shift + rounded_cap$x, x_shift + width - rev(rounded_cap$x)), 
            y = c(1 - y_shift + rounded_cap$y, 1 - y_shift + rev(rounded_cap$y)),
            col = cols[1], 
            border = NA)
    
    if (add_border) {
      points(x = c(x_shift + rounded_cap$x, x_shift + width - rev(rounded_cap$x)), 
             y = c(y_shift - rounded_cap$y, y_shift - rev(rounded_cap$y)), 
             type = "l", lwd = border_lwd)
      points(x = c(x_shift + rounded_cap$x, x_shift + width - rev(rounded_cap$x)), 
             y = c(1 - y_shift + rounded_cap$y, 1 - y_shift + rev(rounded_cap$y)), 
             type = "l", lwd = border_lwd)
    }
  } else if (end_type == "circle") {
    segments(x0 = x_shift, x1 = x_shift + width, y0 = y_shift, y1 = y_shift, lwd = 1, col = cols[1])
    segments(x0 = x_shift, x1 = x_shift + width, y0 = 1 - y_shift, y1 = 1 - y_shift, lwd = 1, col = cols[length(cols)])
    
    circle_cap <- pull_half_circle_coord(width = width, x_center = width/2, y_center = 0)
    polygon(x = x_shift + circle_cap$x, 
            y = y_shift - circle_cap$y, 
            col = cols[1], 
            border = NA)
    polygon(x = x_shift + circle_cap$x, 
            y = 1 - y_shift + circle_cap$y, 
            col = cols[length(cols)], 
            border = NA)
    
    if (add_border) {
      points(x = circle_cap$x + x_shift, 
             y = y_shift - circle_cap$y,
             type = "l", lwd = border_lwd)
      points(x = circle_cap$x + x_shift, 
             y = 1 - y_shift + circle_cap$y, 
             type = "l", lwd = border_lwd)
    }
  } else {
    if (add_border) {
      segments(x0 = x_shift, x1 = x_shift + width, y0 = y_shift, y1 = y_shift, lwd = border_lwd)
      segments(x0 = x_shift, x1 = x_shift + width, y0 = 1 - y_shift, y1 = 1 - y_shift, lwd = border_lwd)
    }
  }
}

#' Draw cartoon chromosome pair for the founder strains of the CC, CC-RIX, and DO
#'
#' This function ...
#' 
#' @export
#' @examples plot_founder_chromosomes(col = qtl2::CCcolors[1], width = 0.15, gap = 0.3)
plot_founder_chromosomes <- function(col,
                                     width = 0.15,
                                     gap = 0.3,
                                     y_shift = NULL,
                                     x_shift = 0,
                                     add_border = TRUE,
                                     border_lwd = 1.5,
                                     end_type = c("rounded", "circle", "square")) {
  
  end_type <- end_type[1]
  plot.new() 
  ## First chromatid
  plot_chromatid(intervals = c(0, 1),
                 cols = col,
                 width = width,
                 y_shift = y_shift,
                 x_shift = x_shift,
                 add_border = add_border,
                 border_lwd = border_lwd,
                 end_type = end_type)
  ## Second chromatid
  plot_chromatid(intervals = c(0, 1),
                 cols = col,
                 width = width,
                 y_shift = y_shift,
                 x_shift = x_shift + gap,
                 add_border = add_border,
                 border_lwd = border_lwd,
                 end_type = end_type)
}

#' Draw cartoon chromosome pair for CC strains
#'
#' This function ...
#' 
#' @export
#' @examples 
#' set.seed(10)
#' plot_cc_chromosomes(intervals = c(0, sort(runif(n = 15, min = 0, max = 1)), 1), 
#'                     col = sample(qtl2::CCcolors, 16, replace = TRUE), width = 0.15, gap = 0.3)
plot_cc_chromosomes <- function(intervals, 
                                cols,
                                width = 0.15,
                                gap = 0.3,
                                y_shift = NULL,
                                x_shift = 0,
                                add_border = TRUE,
                                border_lwd = 1.5,
                                end_type = c("rounded", "circle", "square")) {
  
  end_type <- end_type[1]
  plot.new() 
  ## First chromatid
  plot_chromatid(intervals = intervals,
                 cols = cols,
                 width = width,
                 y_shift = y_shift,
                 x_shift = x_shift,
                 add_border = add_border,
                 border_lwd = border_lwd,
                 end_type = end_type)
  ## Second chromatid
  plot_chromatid(intervals = intervals,
                 cols = cols,
                 width = width,
                 y_shift = y_shift,
                 x_shift = x_shift + gap,
                 add_border = add_border,
                 border_lwd = border_lwd,
                 end_type = end_type)
}

#' Draw cartoon chromosome pair for DO or CC-RIX
#'
#' This function ...
#' 
#' @export
#' @examples 
#' set.seed(10)
#' ## CC-RIX
#' plot_do_chromosomes(c(0, sort(runif(15, 0, 1)), 1), c(0, sort(runif(15, 0, 1)), 1), 
#'                     cols1 = sample(qtl2::CCcolors, 16, replace = TRUE), 
#'                     cols2 = sample(qtl2::CCcolors, 16, replace = TRUE), width = 0.15, gap = 0.3)
#' ## DO
#' plot_do_chromosomes(c(0, sort(runif(50, 0, 1)), 1), c(0, sort(runif(50, 0, 1)), 1), 
#'                     cols1 = sample(qtl2::CCcolors, 51, replace = TRUE), 
#'                     cols2 = sample(qtl2::CCcolors, 51, replace = TRUE), width = 0.15, gap = 0.3)
plot_do_chromosomes <- function(intervals1,
                                intervals2,
                                cols1,
                                cols2,
                                width = 0.15,
                                gap = 0.3,
                                y_shift = NULL,
                                x_shift = 0,
                                add_border = TRUE,
                                border_lwd = 1.5,
                                end_type = c("rounded", "circle", "square")) {
  end_type = end_type[1]
  plot.new() 
  ## First chromatid
  plot_chromatid(intervals = intervals1,
                 cols = cols1,
                 width = width,
                 y_shift = y_shift,
                 x_shift = x_shift,
                 add_border = add_border,
                 border_lwd = border_lwd,
                 end_type = end_type)
  ## Second chromatid
  plot_chromatid(intervals = intervals2,
                 cols = cols2,
                 width = width,
                 y_shift = y_shift,
                 x_shift = x_shift + gap,
                 add_border = add_border,
                 border_lwd = border_lwd,
                 end_type = end_type)
}



