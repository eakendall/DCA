ggplot_tornado <- function(dat,
                           baseline_output = NA,
                           annotate_scale = 0,
                           title="",
                           multiple = F){
  
  if (all(class(dat) != "tornado")) stop("Input data must be tornado class data frame.")
  if (length(baseline_output) != 1) stop("Input baseline_output must be length one.")
  
  output_name <- attr(dat, "output_name")
  
  if (is.na(baseline_output)) {
    
    baseline_output <- attr(dat, "baseline")[output_name]
  }
  
  dat$baseline <- unlist(baseline_output, use.names = FALSE)
  
  # don't strictly need this
  # order output columns as decending and ascending
  datplot <-
    dat[ ,c(output_name, "baseline")] %>%
    dplyr::mutate("min" = apply(., 1, min),
                  "max" = apply(., 1, max)) %>%
    dplyr::select(min, max)
  
  datplot <- cbind(dat, datplot)
  NAMES <- as.character(datplot$names)
  
  
  # order by length of bars
    datplot$names = factor(as.character(datplot$names),
                           levels = rev(unique(datplot$names[order(abs(datplot$min - datplot$max), decreasing = T)])))
    plotorder <- order(abs(datplot$min - datplot$max))
    
    # if("All of the above" %in% levels(datplot$names)) 
    #   datplot$names = factor(as.character(datplot$names),
    #                          levels = c(levels(datplot$names)[startsWith(levels(datplot$names),"All of the above")],
    #                                     (levels(datplot$names))[!startsWith(levels(datplot$names),"All of the above")]))
  
  # check if parameter values are provided
  if (all(NAMES %in% names(datplot))) {
    barLabels <- datplot[, NAMES] %>%
      OpenMx::diag2vec()
  }else{barLabels <- ""}
  
  # shift annotation left or right
  nudge <- (with(datplot, eval(parse(text = output_name)) > baseline) - 0.5) * annotate_scale
  
  ggplot2::ggplot(datplot,
                  aes(names, ymin = min, ymax = max, colour = val)) +
    geom_linerange(size = 5) +
    coord_flip() +
    xlab("") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = dat$baseline, linetype = "dashed") +
    theme_bw() +
    theme(axis.text = element_text(size = 10),
          # legend.position = "none",
          legend.title = element_blank()) +
    annotate("text", x = datplot$names[plotorder], y = datplot[plotorder, output_name] + nudge, label = barLabels[plotorder]) +
    ggtitle(title)
}
