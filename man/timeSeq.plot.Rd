\name{timeSeq.plot}
\alias{timeSeq.plot}
\title{
Display the Significant NDPE Genes
}
\description{
Plot the expression values for the significant genes.
}
\usage{
timeSeq.plot(timeSeq.obj, N)
}
\arguments{
\item{timeSeq.obj}{an object returned by timeSeq function}
\item{N}{the number of significant NPDE genes you would like to plot. It must be a positive integer and should be no greater than 10.}
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
###Data Simulation
mu.produce.exp <- function(x,exon.sd)
{exp(10*(1e3*(x^11*(1-x)^6)+100*(x^3*(1-x)^10)) + 5 + exon.sd)}
ti1 <- c(1:8)/10
ti2 <- c(8:1)/10
nu <- 1000

set.seed(1)
settingC.exon <- matrix(nrow=30, ncol=16)
pp <- 1
exon_length <- c()
for(i in 1:5){
  for(j in 1:3){
    exon_length[pp] <- sample(seq(0.1,0.4,length=3), 1)
    exon_sd <- rnorm(1, sd=exon_length[pp])
    mu1 <- mu.produce.exp(ti1, exon_sd)
    mu2 <- mu.produce.exp(ti2, exon_sd)
    p1 <- nu/(mu1+nu)
    p2 <- nu/(mu2+nu)
    settingC.exon[(i-1)*3 + j, 1:8] <- rnbinom(8,nu,p1)
    settingC.exon[(i-1)*3 + j, 9:16] <- rnbinom(8,nu,p2)
    pp <- pp + 1
  }
}

for(i in 6:10){
  for(j in 1:3){
    exon_length[pp] <- sample(seq(0.1,0.4,length=3), 1)
    exon_sd <- rnorm(1, sd=exon_length[pp])
    mu1 <- mu.produce.exp(ti1, exon_sd)
    p1 <- nu/(mu1+nu)
    settingC.exon[(i-1)*3 + j, 1:8] <- rnbinom(8,nu,p1)
    settingC.exon[(i-1)*3 + j, 9:16] <- rnbinom(8,nu,p1)
    pp <- pp + 1
  }
}

gene.names <- paste("gene", rep(1:10, each=3), sep="")
group.label <- c(rep(1,8), rep(2,8))

##### Model fitting (If you want to run 'timeSeq' by yourself, just remove the 
##### "## " from the head of two lines below.)
## model.fit <- timeSeq(settingC.exon, group.label, gene.names, exon.length=exon_length, 
##   reads=NULL, n_cores = 2, offset=FALSE)

###Load the object returned by timeSeq
data(object_by_timeSeq)

###Plot 
timeSeq.plot(model.fit, 6)
}