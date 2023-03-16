---
title: Genome-Wide Association Studies
summary: A brief description of my document.
date: 2000-01-01
hide:
  - navigation
  - toc
---

<!--
# Put above to hide navigation (left), toc (right) or footer (bottom)

hide:
  - navigation 
  - toc
  - footer 

# You should hide the navigation if there are no subsections
# You should hide the Table of Contents if there are no important titles
-->

<center>
# Genome-Wide Association Studies
</center>

**Updated:** 15 September 2022

This course is an introduction to the method of Genome-Wide Association Studies (GWAS), which quantifies the statistical 
association between a genetic variant and a phenotype (often on disease traits). This course will not focus on using any 
particular software, instead explaining *why* the given analyses are done from a statistical and biological perspective. 

<br>

<!-- AUTHORS -->
<!-- Format: [author name](link to author page){ .md-button } -->
<h4>Authors</h4>

[Samuele Soraggi](https://www.linkedin.com/in/samuelesoraggi/){ .md-button }
[Conor O'Hare](https://github.com/conor-ohare){ .md-button }

<hr>

<!-- OVERVIEW OF COURSE -->
!!! abstract "Overview"
    &#128172; **Syllabus:**    
        1. Understand what is a GWAS is and why we use it  
        2. Statistics of GWAS (regression coefficients, P-values, statistical power, Bayes factors)  
        3. Genetic relatedness and population structure  
        4. Confounding and covariates in GWAS  
        5. Haplotypes, linkage disequilibrum, imputation, fine-mapping  
        6. Linear mixed models and heritability  
        7. Summary statistics and meta-analysis  
        8. Mendelian randomization  
        9. Success and critisism of GWAS  
        10. Human genetics research at FIMM  

    &#128368; **Total Time Estimation:** X hours  

    &#128193; **Supporting Materials:**  [Original course from the University of Helsinki](https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course) | 
                                         [PLINK documentation](https://zzz.bwh.harvard.edu/plink/) | 
                                         [An Introduction to Statistical Learning (for further statistics explanations)](https://www.statlearning.com/)

    &#128203; **License:** [Tutorial Content is licensed under Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/)


!!! warning "Course Requirements"
    - Knowledge of R. It is recommended that you have at least followed our workshop [From Excel to R](https://github.com/Center-for-Health-Data-Science/FromExceltoR)
    - Basic statistics and mathematics skills

<hr>

We do GWAS because a statistical association between a particular physical region of the genome and the phenotype  

- [x] can point to biological mechanisms affecting the phenotype,  
- [x] can allow prediction of the phenotype from genomic information.  

These results may further benefit  

- [x] medicine by leading to molecular or environmental interventions against harmful phenotypes,  
- [x] biotechnology by improving the ways we utilize microbes, plants or animals,  
- [x] forensics by more accurate identification of an individual from a DNA sample,  
- [x] biogeographic ancestry inference of individuals, populations and species,  
- [x] our understanding of the role of natural selection and other evolutionary forces in the living world.  

The genome of an individual remain (nearly) constant throughout the individual’s lifetime. This is a truly remarkable property compared to, e.g., other molecular sources of information (such as metabolomics, metagenomics, transcriptomics, proteomics or epigenomics) or environmental factors that may vary widely across time. Therefore, the genome seems an ideal starting point for scientific research: it needs to be measured only once for an individual and there is no reverse causation from the phenotype to genome (with cancer as an important exception).

### Acknowledgements

- [Center for Health Data Science](https://heads.ku.dk/), University of Copenhagen
- [Matti Pirinen, PhD](https://www.mv.helsinki.fi/home/mjxpirin/), University of Helsinki
