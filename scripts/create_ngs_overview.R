library(googlesheets4)
library(dplyr)
library(ggplot2)
library(ggrepel)

contr_list <- read_sheet("https://docs.google.com/spreadsheets/d/1fS6rmYjRAU26Y2en21hKoWyDn8jqiHMlsDs7bxZyA6Y/edit?usp=sharing"
)

p <- contr_list |> ggplot(
  aes(x = log10(yield_high_gb),
      y = log10(price_per_gb_min),
      color = vendor,
      label = type)) + 
  xlab("log10(throughput (Gb))") +
  ylab("log10(cost per Gb ($))") +
  geom_point(size = 2.5) +
  scale_color_brewer(palette = 'Paired') +
  theme_classic() + 
  geom_text_repel(max.overlaps = 5,
                  color = 'darkgrey', box.padding = 0.3) +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16)) 

png("docs/assets/images/price_vs_yield.png",
    res = 300, width = 7.5, height = 6, units = "in")
print(p)
dev.off()
