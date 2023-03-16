options(warn=-1)

# Install and load GWAS package qqman
install.packages("qqman", repos = "http://cran.us.r-project.org")
library("qqman") 

# Manhattan plot using --logistic
results_log <- read.table("logistic_results.assoc_2.logistic", head=TRUE)
jpeg("Logistic_manhattan.jpeg")
manhattan(results_log,chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot: logistic")
dev.off()

# Manhattan plot using --assoc
results_as <- read.table("assoc_results.assoc", head=TRUE)
jpeg("assoc_manhattan.jpeg")
manhattan(results_as,chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot: assoc")
dev.off()  