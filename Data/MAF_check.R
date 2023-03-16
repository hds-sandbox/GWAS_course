## MAF_check.R
## Produces a histogram of MAF globally and for extreme values

library(ggplot2)

# Read data into R 
maf_freq <- read.table("MAF_check.frq", header =TRUE, as.is=T)

# maf_freq histogram
hist.maf <- ggplot(maf_freq, aes(x=maf_freq[,5])) +
  geom_histogram(col = "black", fill="tomato") + 
  labs(title = "MAF distribution") + 
  xlab("MAF") + 
  ylab("Frequency") +
  theme_bw()

# Save plots 
ggsave(plot=hist.maf, filename="histmaf.png")


