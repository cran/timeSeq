\name{timeSeq.sort}
\alias{timeSeq.sort}
\title{
Sort NDPE Genes by Kullback Leibler Distance Ratios
}
\description{
Sort all genes in the expression dataset by their Kullback Leibler distance ratios.
}
\usage{
timeSeq.sort(genenames, NPDE.ratio, PDE.ratio, table, count)
}
\arguments{
  \item{genenames}{Vector of gene names}
  \item{NPDE.ratio}{Vector of Kullback Leibler distance ratios for NPDE genes}
  \item{PDE.ratio}{Vector of Kullback Leibler distance ratios for PDE genes}
  \item{table}{Gene expression values for each gene.}
  \item{count}{Number of exons for each gene.}
}
\value{
A list with components
  \item{NPDE_list}{Dataframe for NPDE genes sorted by KLRs}
  \item{PDE_list}{Dataframe for PDE genes sorted by KLRs}
  \item{table1}{Gene expression values for each gene, corresponding to NPDE_list}
  \item{table2}{Gene expression values for each gene, corresponding to PDE_list}
}
\references{
Chong Gu. Model diagnostics for smoothing spline ANOVA models. Canadian Journal of Statistics, 32(4):347-358, 2004.

Chong Gu. Smoothing spline ANOVA models. Springer, second edition, 2013.

Chong Gu and Ping Ma. Optimal smoothing in nonparametric mixed-effect models. Annals of Statistics, pages 1357-1379, 2005a.
}
\author{
Fan Gao and Xiaoxiao Sun
}
