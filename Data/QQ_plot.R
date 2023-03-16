options(warn=-1)

# Install and load GWAS package qqman
install.packages("qqman", repos = "http://cran.us.r-project.org")
library("qqman") 

# QQ plot for --logistic
results_log <- read.table("logistic_results.assoc_2.logistic", head=TRUE)
jpeg("QQ-Plot_logistic.jpeg")
qq(results_log$P, main = "Q-Q plot of GWAS p-values (log) using --logistic")
dev.off()

# QQ plot for --assoc
results_as <- read.table("assoc_results.assoc", head=TRUE)
jpeg("QQ-Plot_assoc.jpeg")
qq(results_as$P, main = "Q-Q plot of GWAS p-values (log) using --assoc")
dev.off()