#   Heatmap of the Most Significant NDPE Genes
#   Author: Fan Gao
#   Created 25 Seq 2015. Last modified 10 Oct 2016.

timeSeq.heatmap <- function(timeSeq.obj,n) {
    if (n<=1) {
		cat("n must be a positive integer.\n")
		return(cat("ERROR!"))
	}		
	obj <- timeSeq.obj
	N <- length(obj$count)
	if (n>N) n <- N
	graphics.off()
    ans <- obj$sorted
    m <- ans$table1[1,(1:ans$npde.list$count[1]),]
    now <- ans$npde.list$count[1]
    for (i in c(2:n)) {
        m <- rbind(m,ans$table1[i,(1:ans$npde.list$count[i]),])
        now <- now + ans$npde.list$count[i]
    }
    m <- data.frame(m)
    rownames(m)[1:ans$npde.list$count[1]]=paste(as.character(ans$npde.list$genenames[1])," _",1:ans$npde.list$count[1],sep="")
    now <- ans$npde.list$count[1]
    for (i in c(2:n)) {
        rownames(m)[(now+1):(now+ans$npde.list$count[i])]=paste(as.character(ans$npde.list$genenames[i])," _",1:ans$npde.list$count[i],sep="")     
        now <- now + ans$npde.list$count[i]
    }	 
    colnames(m) <- paste("Time",1:(obj$group1.length+obj$group2.length),sep="")
    annotation <- data.frame(Group=factor(c(rep(0,obj$group1.length),rep(1,obj$group2.length)),labels=c("Group1", "Group2")))
    rownames(annotation) <- colnames(m)
    var1 <- c("navy", "darkgreen")
    names(var1) <- c("Group1","Group2")   
    ann.colors <- list(var1=var1)
    pheatmap(m,scale="row",show_rownames=TRUE,
             show_colnames=TRUE,cluster_rows=FALSE,
             cluster_cols=FALSE,
             main=paste("Heatmap of Top ",n," significant NPDE genes",sep=""),
             display_numbers=FALSE,annotation=annotation)
    return(m)
}
