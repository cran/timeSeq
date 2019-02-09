# Sort NDPE and PDE Genes by Kullback Leibler Distance Ratios
# Author: Fan Gao
# Created 25 Seq 2015. Last modified 10 Oct 2016.

timeSeq.sort <- function(genenames,npde,pde,table,count,pvalue) {
    list.length <- length(genenames)
    if(pvalue){
        npde.order <- order(npde)
        npde.list <- data.frame(genenames=genenames[npde.order],pvalues=npde[npde.order],count=count[npde.order])
        table1 <- table[npde.order,,]
        pde.order <- order(pde)
        pde.list <- data.frame(genenames=genenames[pde.order],pvalues=pde[pde.order],count=count[pde.order])
        table2 <- table[pde.order,,]
    } else {
        npde.order <- order(-npde)
        npde.list <- data.frame(genenames=genenames[npde.order],ratios=npde[npde.order],count=count[npde.order])
        table1 <- table[npde.order,,]
        pde.order <- order(-pde)
        pde.list <- data.frame(genenames=genenames[pde.order],ratios=pde[pde.order],count=count[pde.order])
        table2 <- table[pde.order,,]
    }
    out <- list(npde.list=npde.list, 
                pde.list=pde.list,
                table1=table1,
                table2=table2)
    out
}

