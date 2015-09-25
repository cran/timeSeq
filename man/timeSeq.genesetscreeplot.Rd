\name{timeSeq.genesetscreeplot}
\alias{timeSeq.genesetscreeplot}
\title{
Scree Plot of Kullback-Leibler Distance Ratios for Gene Sets
}
\description{
Plot the scree plot for all gene sets in the dataset.
}
\usage{
timeSeq.genesetscreeplot(timeSeq.geneset)
}
\arguments{
\item{timeSeq.geneset}{an object returned by timeSeq.geneset.}
}
\value{
  \item{NPDE.list}{a dataframe for expression values.}
}
\references{
Chong Gu. Model diagnostics for smoothing spline ANOVA models. Canadian Journal of Statistics, 32(4):347-358, 2004.

Chong Gu. Smoothing spline ANOVA models. Springer, second edition, 2013.

Chong Gu and Ping Ma. Optimal smoothing in nonparametric mixed-effect models. Annals of Statistics, pages 1357-1379, 2005a.
}
\author{
Fan Gao and Xiaoxiao Sun
}
\examples{
#### Simulation
set.seed(1)
geneset.data <- matrix(nrow=90, ncol=16)
mu.produce.exp <- function(x,exon.sd)
{exp(10*(1e3*(x^11*(1-x)^6)+100*(x^3*(1-x)^10)) + 5 + exon.sd)}
ti1 <- c(1:8)/10
ti2 <- c(8:1)/10
nu <- 20000
pp <- 1
exon_length <- c()
for(i in 1:10){
  for(j in 1:3){
    exon_length[pp] <- c(0.1,0.25,0.4)[j]
    exon_sd <- rnorm(1, sd=exon_length[pp])
    mu1 <- mu.produce.exp(ti1, exon_sd)
    mu2 <- mu.produce.exp(ti2, exon_sd)
    geneset.data[(i-1)*3 + j, 1:8] <- rnbinom(8,nu,mu=mu1)
    geneset.data[(i-1)*3 + j, 9:16] <- rnbinom(8,nu,mu=mu2)
    pp <- pp + 1
  }
}

for(i in 11:30){
  for(j in 1:3){
    exon_length[pp] <- c(0.1,0.25,0.4)[j]
    exon_sd <- rnorm(1, sd=exon_length[pp])
    mu1 <- mu.produce.exp(ti1, exon_sd)
    geneset.data[(i-1)*3 + j, 1:8] <- rnbinom(8,nu,mu=mu1)
    geneset.data[(i-1)*3 + j, 9:16] <- rnbinom(8,nu,mu=mu1)
    pp <- pp + 1
  }
}

### Gene names should match with the data. 
gene.names=paste("g",rep(1:30, each=3),sep="")

### create some random gene sets
gene.sets=vector("list",3)
for(i in 1:3){gene.sets[[i]]=paste("g",(1:30)[((i-1)*10+1):(i*10)],sep="")}
geneset.names=paste("set",as.character(1:3),sep="")
group.label <- c(rep(1,8), rep(2,8))

##### Model fitting (If you want to run 'timeSeq.geneset' by yourself, just remove the 
##### "## " from the head of two lines below.)
## geneset.fit <-  timeSeq.geneset(geneset.data, group.label, gene.sets, gene.names, 
##   geneset.names, n_cores = 2, offset=FALSE)
  

###Load the object returned by timeSeq.geneset
data(object_by_timeSeq.geneset)

### Scree plot to check cutoff values.
scree.plot <- timeSeq.genesetscreeplot(geneset.fit)

### Treat top 10 gene sets as NPDE gene sets. 
scree.plot[1 : 10,]
}
