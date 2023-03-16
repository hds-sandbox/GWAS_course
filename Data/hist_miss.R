## hist_miss.R
## Produces a histogram of missing for individuals and loci 

library(ggplot2)

# Read data into R 
indmiss<-read.table(file="plink.imiss", header=TRUE)
snpmiss<-read.table(file="plink.lmiss", header=TRUE)

# plink.imiss histogram
hist.imiss <- ggplot(indmiss, aes(x=indmiss[,6])) +
  geom_histogram(binwidth = 0.001, col = "black", fill="tomato") + 
  labs(title = "Frequency of missingness rates in individuals") + 
  xlab("Frequency") + 
  ylab("Missingness in Individuals") + 
  theme_bw()

# plink.lmiss histogram
hist.lmiss <- ggplot(snpmiss, aes(x=snpmiss[,5])) +
  geom_histogram(binwidth = 0.005, col = "black", fill="tomato") + 
  labs(title = "Frequency of missingness rates in SNPs") + 
  xlab("Frequency") + 
  ylab("Missingness in SNPs") + 
  theme_bw()

show(hist.lmiss)

# Save plots 
ggsave(plot=hist.imiss, filename="histimiss.png")
ggsave(plot=hist.lmiss, filename="histlmiss.png")
