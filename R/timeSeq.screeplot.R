#   Scree Plot of Kullback Leibler Distance Ratios
#   Author: Fan Gao
#   Created 25 Seq 2015. Last modified 10 Oct 2016.

timeSeq.screeplot <- function(timeSeq.obj,type=c("barplot","lines")) {
	graphics.off()
    obj.sorted <- timeSeq.obj$sorted
    obj.sorted <- obj.sorted$npde.list
    gene.length <- sum(!is.na(obj.sorted[,2]))
    gene.index <- !(is.na(obj.sorted[,2]))
    if(timeSeq.obj$pvalue){
        if (type=="barplot") {
            obj.sorted$genenames <- factor(obj.sorted$genenames, 
                                           levels=unique(as.character(obj.sorted$genenames)))    
            barchart(pvalues~genenames,data=obj.sorted, 
                     main="P-values of Genes",xlab="Genes",ylab="P-values",
                     col = "aliceblue")    
        } else if (type=="lines") {
            plot(x=1:gene.length,y=obj.sorted$pvalues[gene.index],
                 type="b",pch=21,col="red",xaxt="n",lty = 2, 
                 main="P-values of Genes",xlab="Genes",ylab = "P-values")    
            axis(1,1:gene.length,obj.sorted$genenames[gene.index],col.axis="blue")
	    } else {
		    cat("The type of plot must be bar plot or line graph.\n")
		    return(cat("ERROR!"))
	}
    } else {
        if (type=="barplot") {
            obj.sorted$genenames <- factor(obj.sorted$genenames, 
                                           levels=unique(as.character(obj.sorted$genenames)))    
            barchart(ratios~genenames,data=obj.sorted, 
                     main="Ratios of Genes",xlab="Genes",ylab="Ratios",
                     col = "aliceblue")    
        } else if (type=="lines") {
            plot(x=1:gene.length,y=obj.sorted$ratios[gene.index],
                 type="b",pch=21,col="red",xaxt="n",lty = 2, 
                 main="Ratios of Genes",xlab="Genes",ylab = "Ratios")    
            axis(1,1:gene.length,obj.sorted$genenames[gene.index],col.axis="blue")
	    } else {
		    cat("The type of plot must be bar plot or line graph.\n")
		    return(cat("ERROR!"))
	    }

    }
}
