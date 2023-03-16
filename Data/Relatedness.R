## Relatedness.R
## Produces a histogram of missing for individuals and loci 

library(ggplot2)

# Read data into R 
relatedness <- read.table("pihat_min0.2.genome", header=T)
relatedness_zoom <- read.table("zoom_pihat.genome", header=T)

# Relatedness plot
plot.relatedness <- ggplot(relatedness) +
  geom_point(aes(x=Z0, y=Z1, col=RT)) + 
  xlim(0,1) + 
  ylim(0,1) +
  labs(title = "Relatedness", col = "Relationship") + 
  xlab("Z0") + 
  ylab("Z1") + 
  theme_bw()

# Zoomed relatedness plot
plot.relatedness_zoom <- ggplot(relatedness_zoom) +
  geom_point(aes(x=Z0, y=Z1, col=RT)) + 
  xlim(0,0.005) + 
  ylim(0,1) +
  labs(title = "Relatedness", col = "Relationship") + 
  xlab("Z0") + 
  ylab("Z1") + 
  theme_bw()

# Relatedness histogram
hist.relatedness <- ggplot(relatedness, aes(x=relatedness[,10])) +
  geom_histogram(binwidth = 0.02, col = "black", fill="tomato") + 
  labs(title = "Histogram of relatedness (Pi Hat)") + 
  xlab("Pi Hat") + 
  ylab("Frequency") + 
  theme_bw()

show(hist.relatedness)

# Save plots 
ggsave(plot=plot.relatedness, filename="relatedness.png")
ggsave(plot=plot.relatedness_zoom, filename="relatedness_zoom.png")
ggsave(plot=hist.relatedness, filename="histrelatedness.png")


