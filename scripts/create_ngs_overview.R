library(googlesheets4)
library(dplyr)
library(ggplot2)
library(ggrepel)

seq_method_table <- read_sheet("https://docs.google.com/spreadsheets/d/1fS6rmYjRAU26Y2en21hKoWyDn8jqiHMlsDs7bxZyA6Y/edit?usp=sharing"
)

vendors <- unique(seq_method_table$vendor)
vendor_colors <- RColorBrewer::brewer.pal(length(vendors), "Paired")
names(vendor_colors) <- vendors

plot_seq_method <- function(df, x, y, shape = NULL, xlab = x, ylab = y,
                            shape_lab = shape,
                            vendor_colors) {
  p <- ggplot(df,
         aes(x = log10(.data[[x]]),
             y = log10(.data[[y]]),
             shape = if (is.null(shape)) NULL else .data[[shape]],
             color = vendor,
             label = type)) + 
    xlab(xlab) +
    ylab(ylab) +
    geom_point(size = 2.5) +
    scale_color_manual(values = vendor_colors) +
    theme_classic() + 
    geom_text_repel(max.overlaps = 5,
                    color = 'darkgrey', box.padding = 0.3) +
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16)) 
  
  if (!is.null(shape)) {
    p <- p + labs(shape = shape_lab)
  }
  
  return(p)
}

p <- seq_method_table |>
  filter(vendor != "Sanger") |>
  plot_seq_method(x = "yield_high_gb",
                  y = "price_per_gb_min",
                  xlab = "log10(throughput in Gb)",
                  ylab = "log10(cost per Gb in $)",
                  vendor_colors = vendor_colors)

png("docs/assets/images/price_vs_yield.png",
    res = 300, width = 7.5, height = 6, units = "in")
print(p)
dev.off()

p <- seq_method_table |>
  plot_seq_method(x = "yield_high_gb",
                  y = "price_per_gb_min",
                  xlab = "log10(throughput in Gb)",
                  ylab = "log10(cost per Gb in $)",
                  vendor_colors = vendor_colors)

png("docs/assets/images/price_vs_yield_incl_sanger.png",
    res = 300, width = 7.5, height = 6, units = "in")
print(p)
dev.off()

p <- seq_method_table |>
  plot_seq_method(x = "read_length_max",
                  y = "yield_high_gb",
                  xlab = "log10(max read length)",
                  ylab = "log10(throughput (Gb))",
                  shape = "read_type",
                  shape_lab = "read type",
                  vendor_colors = vendor_colors)

png("docs/assets/images/yield_vs_max_length.png",
    res = 300, width = 7.5, height = 6, units = "in")
print(p)
dev.off()

p <- seq_method_table |>
  plot_seq_method(x = "read_length_max",
                  y = "price_per_gb_min",
                  xlab = "log10(max read length)",
                  ylab = "log10(cost per Gb in $)",
                  shape = "read_type",
                  shape_lab = "read type",
                  vendor_colors = vendor_colors)

png("docs/assets/images/price_vs_max_length.png",
    res = 300, width = 7.5, height = 6, units = "in")
print(p)
dev.off()
  