# IMAS
##Scope
A balance of alternatively spliced (AS) transcript isoforms is critical for a normal phenotype, and dysregulated expression might result in phenotypic and clinical alterations in certain individuals and groups. For example, the expression ratio (denoted as PSI: Percent Spliced In) of AS transcript isoforms or AS exons is modulated in genes that are associated with many complex diseases, such as cancer and diabetes. Moreover, molecular factors that aﬀect the expression ratio, such as genetic variations (sQTLs) and methylation (deﬁned as me-sQTLs), give additional insight into the molecular mechanisms that explain the regulatory elements in splicing machinery. 
The Bioconductor R package IMAS oﬀers two components.

First, IMAS provides an integrative analysis tool to detect diﬀerential AS events (e.g., exon skipping, intron retention, and 3- and 5-prime alternative splicing sites) using paired and junction reads that are created from high-throughput RNA-seq data.

Second, IMAS links AS exons that are diﬀerentially expressed across groups, deﬁned as PSI, to molecular factors (eg, SNPs and DNA methylation). IMAS can identify sQTLs and me-sQTLs to interpret their molecular functions. In addition, IMAS can incorporate the clinical information of a patient and determine whether exon skipping or inclusion aﬀects the clinical outcomes of disease groups.

Statistical methods (linear regression model or generalized linear mixed model in sQTLs) are used to identify signiﬁcant diﬀerences in PSI values across groups, sQTLs, me-sQTLs, and clinical outcome associations. IMAS can accept mapped bam ﬁles, which may be created from several mapping tools. IMAS analyzes all paired-end and junction reads. Moreover, all functions of this package accept a simple FPKM matrix dataset of transcripts (output ﬁle format of assemble tools, such as cuﬄinks) through IVAS, which is an R/bioconductor package for identifying genetic variants that aﬀect alternative splicing patterns with the FPKM matrix dataset. 

Finally, IMAS provides a function for visualizing all identiﬁed AS transcript isoforms, sQTLs, me-sQTLs, and clinical outcomes.

##Installation
To install IMAS package, start R and enter:
``` r
source("https://bioconductor.org/biocLite.R")
biocLite("IMAS")
```
To view vignette of IMAS, start R and enter:
``` r
browseVignettes("IMAS")
```
The document provides an implementation detail with example data sets.
