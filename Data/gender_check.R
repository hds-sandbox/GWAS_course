## hist_miss.R
## Produces a histogram of missing for individuals and loci 

library(ggplot2)

# Read data into R
gender <- read.table("plink.sexcheck", header=T,as.is=T)
male=subset(gender, gender$PEDSEX==1)
female=subset(gender, gender$PEDSEX==2)

# plink.imiss histogram
hist.gender <- ggplot(gender, aes(x=gender[,6])) +
  geom_histogram(col = "black", fill="tomato") + 
  labs(title = "All sexes") + 
  xlab("F") + 
  theme_bw()

# plink.imiss histogram
hist.male <- ggplot(male, aes(x=male[,6])) +
  geom_histogram(col = "black", fill="tomato") + 
  labs(title = "Male") + 
  xlab("F") + 
  theme_bw()


# plink.imiss histogram
hist.female <- ggplot(female, aes(x=female[,6])) +
  geom_histogram(col = "black", fill="tomato") + 
  labs(title = "Female") + 
  xlab("F") + 
  theme_bw()


# Save plots 
ggsave(plot=hist.gender, filename="histgender.png")
ggsave(plot=hist.male, filename="histmale.png")
ggsave(plot=hist.female, filename="histfemale.png")
