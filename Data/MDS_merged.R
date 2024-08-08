## MDS_merged.R
## MDS on our samples + anchoring from 1000 Genomes

library(ggplot2)

# Read data into R 
data<- read.table(file="MDS_merge2.mds",header=TRUE)
race<- read.table(file="racefile.txt",header=TRUE)
datafile<- merge(data,race,by=c("IID","FID"))

scatter.mds <- ggplot(datafile, aes(x=datafile[,4], y=datafile[,5], color=datafile[,14])) +
  geom_point() +
  xlab("MD Component 1") + 
  ylab("MD Component 2") +
  labs(color="Population") +
  theme_bw()

# Save plots 
ggsave(plot=scatter.mds, filename="MDS.png")

