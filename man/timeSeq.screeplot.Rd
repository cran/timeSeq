\name{timeSeq.screeplot}
\alias{timeSeq.screeplot}
\title{
Scree Plot of Kullback Leibler Distance Ratios}
\description{
Plot the scree plot for all genes in the dataset.
}
\usage{
timeSeq.screeplot(timeSeq.obj, type = c("barplot", "lines"))
}
\arguments{
  \item{timeSeq.obj}{an object returned by timeSeq function
}  
  \item{type}{type of plot: "barplot" for bar plot and "lines" for line graph.
}
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
model.fit <- timeSeq(data.count, group.label, gene.names, reads, exon.level = FALSE,  n.cores = 1)
timeSeq.screeplot(model.fit, "lines")

}
