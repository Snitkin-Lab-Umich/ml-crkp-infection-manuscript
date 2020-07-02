# boxplot similar to base R

library(ggplot2)

geom_boxplot2 <- function(mapping = NULL, data = NULL, stat = "boxplot", position = "dodge2", 
                          ..., outlier.colour = NULL, outlier.color = NULL, outlier.fill = NULL, 
                          outlier.shape = 1, outlier.size = 1.5, outlier.stroke = 0.5, 
                          outlier.alpha = NULL, notch = FALSE, notchwidth = 0.5, varwidth = FALSE, 
                          na.rm = FALSE, show.legend = NA, inherit.aes = TRUE,
                          linetype = "dashed"){
  list(
    geom_boxplot(mapping = mapping, data = data, stat = stat, position = position,
                 outlier.colour = outlier.colour, outlier.color = outlier.color, 
                 outlier.fill = outlier.fill, outlier.shape = outlier.shape, 
                 outlier.size = outlier.size, outlier.stroke = outlier.stroke, 
                 outlier.alpha = outlier.alpha, notch = notch, 
                 notchwidth = notchwidth, varwidth = varwidth, na.rm = na.rm, 
                 show.legend = show.legend, inherit.aes = inherit.aes, 
                 linetype = linetype, ...),
    stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.shape = 1) ,
    stat_boxplot(geom = "errorbar", width=0.5, aes(ymin = ..ymax..)) ,
    stat_boxplot(geom = "errorbar", width=0.5, aes(ymax = ..ymin..)) ,
    theme_bw(), 
    theme(plot.title = element_text(hjust = 0.5,  
                                    size = 14,
                                    face = "bold"),
          panel.border = element_rect(linetype = "solid",
                                      colour = "black", fill = "NA", size = 0.5))
  )
}
