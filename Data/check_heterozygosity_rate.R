## check_heterozygosity_rate.R
## Produces a histogram of missing for individuals and loci 

library(ggplot2)

# Read data into R 
het <- read.table("R_check.het", head=TRUE)

het$HET_RATE = (het$"N.NM." - het$"O.HOM.")/het$"N.NM."

# plink.imiss histogram
hist.het <- ggplot(het, aes(x=HET_RATE)) +
  geom_histogram(binwidth = 0.001, col = "black", fill="tomato") + 
  labs(title = "Heterozygosity Rates") + 
  xlab("Heterozygosity Rate") + 
  ylab("Frequency") + 
  theme_bw()


show(hist.het)

# Save plots 
ggsave(plot=hist.het, filename="heterozygosity.png")
