\name{timeSeq.sort}
\alias{timeSeq.sort}
\title{
Sort NDPE Genes by Kullback Leibler Distance Ratios or P-values
}
\description{
Sort all genes in the dataset by their Kullback Leibler distance ratios or p-values.
}
\usage{
timeSeq.sort(genenames,npde,pde,table,count,pvalue)
}
\arguments{
  \item{genenames}{a vector of gene names.}
  \item{npde}{a vector of Kullback Leibler distance ratios or p-values for NPDE genes.}
  \item{pde}{a vector of Kullback Leibler distance ratios or p-values for PDE genes.}
  \item{table}{gene expression values.}
  \item{count}{the number of exons or replicates for each gene.}
  \item{pvalue}{logical:indicating if p-values are returned.}
}
\value{
A list with components
  \item{npde.list}{dataframe of NPDE genes sorted by KLRs or p-values.}
  \item{pde.list}{dataframe of PDE genes sorted by KLRs or p-values.}
  \item{table1}{gene expression values for each gene, corresponding to npde.list.}
  \item{table2}{gene expression values for each gene, corresponding to pde.list.}
}
\references{
Sun, Xiaoxiao, David Dalpiaz, Di Wu, Jun S. Liu, Wenxuan Zhong, and Ping Ma. "Statistical inference for time course RNA-Seq data using a negative binomial mixed-effect model." BMC Bioinformatics, 17(1):324, 2016.

Chong Gu. Model diagnostics for smoothing spline ANOVA models. Canadian Journal of Statistics, 32(4):347-358, 2004.

Chong Gu. Smoothing spline ANOVA models. Springer, second edition, 2013.

Chong Gu and Ping Ma. Optimal smoothing in nonparametric mixed-effect models. Annals of Statistics, pages 1357-1379, 2005a.

Wood (2001) mgcv:GAMs and Generalized Ridge Regression for R. R News 1(2):20-25
}
\author{
Fan Gao and Xiaoxiao Sun
}
