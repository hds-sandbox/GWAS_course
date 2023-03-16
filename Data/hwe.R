## MAF_check.R
## Produces a histogram of MAF globally and for extreme values

library(ggplot2)

# Read data into R 
hwe<-read.table (file="plink.hwe", header=TRUE)
hwe_zoom<-read.table (file="plinkzoomhwe.hwe", header=TRUE)

# maf_freq histogram
hist.hwe <- ggplot(hwe, aes(x=hwe[,9])) +
  geom_histogram(col = "black", fill="tomato") + 
  labs(title = "HWE distribution") + 
  xlab("HWE p-value") + 
  ylab("Frequency") +
  theme_bw()

# maf_freq histogram
hist.hwe_below_threshold <- ggplot(hwe_zoom, aes(x=hwe_zoom[,9])) +
  geom_histogram(binwidth = 0.0000015, col = "black", fill="tomato") + 
  labs(title = "HWE distribution for strongly deviating SNPs only") + 
  xlab("HWE p-value") + 
  ylab("Frequency") +
  theme_bw()

# Save plots 
ggsave(plot=hist.hwe, filename="histhwe.png")
ggsave(plot=hist.hwe_below_threshold, filename="histhwe_below_threshold.png")
