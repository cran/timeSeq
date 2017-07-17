\name{timeSeq.heatmap}
\alias{timeSeq.heatmap}
\title{
Heatmap of the Most Significant NDPE Genes
}
\description{
Plots the heatmap for the significant NDPE genes.
}
\usage{
timeSeq.heatmap(timeSeq.obj, n)
}
\arguments{
\item{timeSeq.obj}{an object returned by timeSeq function}
\item{n}{the number of the most significant NPDE genes. It must be a positive integer.}
}
\references{
Sun, Xiaoxiao, David Dalpiaz, Di Wu, Jun S. Liu, Wenxuan Zhong, and Ping Ma. "Statistical inference for time course RNA-Seq data using a negative binomial mixed-effect model." BMC Bioinformatics, 17(1):324, 2016.

Chong Gu. Model diagnostics for smoothing spline ANOVA models. Canadian Journal of Statistics, 32(4):347-358, 2004.

Chong Gu. Smoothing spline ANOVA models. Springer, second edition, 2013.

Chong Gu and Ping Ma. Optimal smoothing in nonparametric mixed-effect models. Annals of Statistics, pages 1357-1379, 2005a.
}
\author{
Fan Gao and Xiaoxiao Sun
}

\examples{

data(simulate.dt)
attach(simulate.dt)
model.fit <- timeSeq(data.count, group.label, gene.names, exon.level = FALSE)
timeSeq.heatmap(model.fit, n = 10)

}

