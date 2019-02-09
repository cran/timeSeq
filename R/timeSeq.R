#	Fit negative binomial mixed effect model for each gene
#	Authors: Xiaoxiao Sun & Fan Gao
#	Created 25 Seq 2015. Last modified 08 February 2019.

timeSeq <- function(data.count,group.label,gene.names,exon.length=NULL,exon.level=FALSE,pvalue=TRUE) {
    lev.label <- unique(group.label)
    if(length(lev.label)>2){
        stop("Currently, only two experimental conditions are allowed.")
    }
    if(exon.level & pvalue){
        cat("No p-values for exon level data!","\n")
    }
    if(is.factor(gene.names)){
        gene.names <- as.character(gene.names)
    }
    group.length <- length(group.label)
    group.1 <- length(group.label[group.label==lev.label[1]])
    group.2 <- length(group.label[group.label==lev.label[2]])
    gene.list <- unique(gene.names)
    npde <- NULL
    pde <- NULL
    for(genenames in gene.list){
        ind <- which(gene.list%in%genenames)
        if (exon.level){
            gene <- matrix(data.count[gene.names==genenames,],ncol=group.length)
            n <- dim(gene)[1]
            if (n<1) return(sprintf("Gene: %s has no data!",genenames))
            colnames(gene) <- group.label
            rownames(gene) <- NULL
            m.dt <- melt(gene)
            if (!is.null(exon.length)){
                el <-  rep(exon.length[gene.names==genenames],group.length)
                el <- log(el)
                if(sum(el<=0)){
                    el[which(el<=0)] <- 0.1
                }
                ndt <- data.frame(id=as.factor(m.dt[,1]),g=as.factor(m.dt[,2]),t=c(rep(1:group.1,each=n),rep(1:group.2,each=n)),leng=el,count=m.dt[,3])
                ran <- mkran(~1|id,ndt)
                randm <- list(z=el*ran$z,sigma=ran$sigma,init=ran$init)
                gss.npde <- try(gssanova(count~t+g+t:g,data=ndt,family="nbinomial",random=randm),silent=TRUE)
                gss.pde <- try(gssanova(count~t+g,data=ndt,family="nbinomial",random=randm),silent=TRUE)
            } else {
                ndt <- data.frame(id=as.factor(m.dt[,1]),g=as.factor(m.dt[,2]),t=c(rep(1:group.1,each=n),rep(1:group.2,each=n)),count=m.dt[,3])
                gss.npde <- try(gssanova(count~t+g+t:g,data=ndt,family="nbinomial",random=~1|id),silent=TRUE)
                gss.pde <- try(gssanova(count~t+g,data=ndt,family="nbinomial",random=~1|id),silent=TRUE)
            }
            if (("try-error" %in% class(gss.npde))) npde[ind] <- NA else {
                pj.int <- try(project(gss.npde,c("t","g")),silent=TRUE)
                if (("try-error" %in% class(pj.int))) npde[ind] <- NA else npde[ind] <- pj.int$ratio
            }
            if (("try-error" %in% class(gss.pde))) {
                pde[ind] <- NA
            } else {
                pj.add <- try(project(gss.pde,c("t")),silent=TRUE)
                if (("try-error" %in% class(pj.add))) pde[ind] <- NA else pde[ind] <- pj.add$ratio
            }
        } else {
            gene <- matrix(data.count[gene.names==genenames,],ncol=group.length)
            n <- dim(gene)[1]
            if (n < 1) return(sprintf("Gene: %s has no data!",genenames))
            colnames(gene) <- group.label
            rownames(gene) <- NULL
            m.dt <- melt(gene)
            ndt <- data.frame(id=as.factor(m.dt[, 1]),g=as.factor(m.dt[,2]),t=c(rep(1:group.1,each=n),rep(1:group.2,each=n)),count=m.dt[, 3])
            if(!pvalue){
                gss.npde <- try(gssanova(count~t+g+t:g,data=ndt,family="nbinomial"),silent=TRUE)
                gss.pde <- try(gssanova(count~t+g,data=ndt,family="nbinomial"),silent=TRUE)
                if (("try-error" %in% class(gss.npde))) npde[ind] <- NA else {
                    pj.int <- try(project(gss.npde,c("t","g")),silent=TRUE)
                    if (("try-error" %in% class(pj.int))) npde[ind] <- NA else npde[ind] <- pj.int$ratio
                }
                if (("try-error" %in% class(gss.pde))) pde[ind] <- NA else {
                    pj.add <- try(project(gss.pde, c("t")),silent=TRUE)
                    if (("try-error" %in% class(pj.add))) pde[ind] <- NA else pde[ind] <- pj.add$ratio
                }
            } else {
                gam.npde <- try(gam(count~ti(t)+g+ti(t,g,bs=c("cr","fs")),family=nb(),data=ndt),silent=TRUE)
                gam.pde <- try(gam(count~ti(t)+g,family=nb(),data=ndt),silent=TRUE)
                if (("try-error" %in% class(gam.npde))) npde[ind] <- NA else {
                    npde[ind] <- summary(gam.npde)$s.table[2,4]
                }
                if (("try-error" %in% class(gam.pde))) pde[ind] <- NA else {
                    pde[ind] <- summary(gam.pde)$p.table[2,4]
                }
            }
        }
    }
	count <- numeric(length(gene.list))
	max.length <- 0
	for (i in c(1:length(gene.list))) {
	    count[i] <- dim(matrix(data.count[gene.names==gene.list[i],],ncol=group.length))[1]
		if (count[i]>max.length) max.length <- count[i]
	}
	table <- array(0,c(length(gene.list),max.length,group.length))
	for (i in c(1:length(gene.list))) {
	    m <- unlist(data.count[gene.names==gene.list[i],1:group.length])
		m <- t(matrix(m,c(count[i],group.length)))
		for (j in c(1:count[i])) table[i,j,] <- m[,j]
	}
    sorted <- timeSeq.sort(gene.list,npde,pde,table,count,pvalue)
   	out <- list(sorted=sorted,count=count,NPDE=npde,PDE=pde,
               genenames=gene.list,table=table,data=data.count,
               gene.names=gene.names,exon.length=exon.length,
		       group.label=group.label,group.length=group.length,
               group1.length=group.1,group2.length=group.2,
               exon.level=exon.level,pvalue=pvalue)
	class(out) <- "timeSeq.obj"
	return(out)
}
