
#	Fit negative binomial mixed effect model for each gene
#	Authors: Xiaoxiao Sun & Fan Gao
#	Created 25 Seq 2015. Last modified 10 Oct 2016.

timeSeq = function(data.count, group.label, gene.names, reads = NULL, exon.length = NULL, exon.level = TRUE, 
	n.cores = NULL, offset = TRUE, iterations=10, p.values=FALSE) {

    lev.label = unique(group.label)
    if(length(lev.label)>2){
        stop("Currently, only two experimental conditions are allowed.")
    }

    if(is.factor(gene.names)){
    	gene.names = as.character(gene.names)
    }

    group.length = length(group.label)
    group.1 = length(group.label[group.label == lev.label[1]])
    group.2 = length(group.label[group.label == lev.label[2]])
    List = unique(gene.names)

    if(is.null(reads)){
        lib.size = rep(1, group.length)
    } else {
    	lib.size = (reads/mean(reads))
    }


	desired.cores = n.cores
	if (length(desired.cores) == 0) {
	    max.cores = detectCores()
	    desired.cores = max.cores - 1		
	} 
    cl = makeCluster(desired.cores)
    registerDoParallel(cl)

	fits = foreach(genenames = List, .packages = c("reshape", "gss")) %dopar% {

		if (exon.level) {
			gene = matrix(data.count[gene.names == genenames, ], ncol = group.length)
			n = dim(gene)[1]
			if (n < 1) return(sprintf("Gene: %s has no data!", genenames))
			colnames(gene) = group.label
			rownames(gene) = NULL

			m.dt = melt(gene)
			offset.lib = rep(lib.size, each = n)
			if (!is.null(exon.length)) {
				el =  rep(exon.length[gene.names == genenames], group.length)
				el = log(el)
                if(sum(el<=0)){
                	el[which(el<=0)] = 0.1
                }

				ndt = data.frame(id = as.factor(m.dt[, 1]), g = as.factor(m.dt[, 2]), t = c(rep(1:group.1, each = n), rep(1:group.2, each = n)), leng = el, count = m.dt[, 3])
				ran <- mkran(~1|id, ndt)
                randm = list(z=el*ran$z, sigma=ran$sigma, init=ran$init)

		        if(offset){
        			fit1 = try(gssanova(count ~ t + g + t:g, offset = offset.lib, data = ndt, family = "nbinomial", random = randm, skip.iter=TRUE), silent = TRUE)
		            fit2 = try(gssanova(count ~ t + g, data = ndt, offset = offset.lib, family = "nbinomial", random = randm, skip.iter=TRUE), silent = TRUE)
		        } else {
          			fit1 = try(gssanova(count ~ t + g + t:g, data = ndt, family = "nbinomial", random = randm, skip.iter=TRUE), silent = TRUE)
			        fit2 = try(gssanova(count ~ t + g, data = ndt, family = "nbinomial", random = randm, skip.iter=TRUE), silent = TRUE)
        		}
			
			} else {
				ndt <- data.frame(id = as.factor(m.dt[, 1]), g = as.factor(m.dt[, 2]), t = c(rep(1:group.1, each = n), rep(1:group.2, each = n)), 
				count = m.dt[, 3])
		        if(offset){
        			fit1 = try(gssanova(count ~ t + g + t:g, offset = offset.lib, data = ndt, family = "nbinomial", random = ~1 | id, skip.iter=TRUE), silent = TRUE)
		            fit2 = try(gssanova(count ~ t + g, offset = offset.lib, data = ndt, family = "nbinomial", random = ~1 | id, skip.iter=TRUE), silent = TRUE)
        		} else {
		            fit1 = try(gssanova(count ~ t + g + t:g, data = ndt, family = "nbinomial", random = ~1 | id, skip.iter=TRUE), silent = TRUE)
		            fit2 = try(gssanova(count ~ t + g, data = ndt, family = "nbinomial", random = ~1 | id, skip.iter=TRUE), silent = TRUE)
        		}				
			}

			if (("try-error" %in% class(fit1))) NPDE.ratio = NA
			else {
				pj.tc.temp = try(project(fit1, c("t", "g")), silent = TRUE)
				if (("try-error" %in% class(pj.tc.temp))) NPDE.ratio = NA else NPDE.ratio = pj.tc.temp$ratio
			}
			if (("try-error" %in% class(fit2))) {
				PDE.ratio = NA
				perm.nu = NA
				perm.eta = NA
			} else {
				perm.nu = fit2$nu 
				perm.eta = fit2$eta
				pj.tg.temp = try(project(fit2, c("t")), silent = TRUE)
				if (("try-error" %in% class(pj.tg.temp))) PDE.ratio = NA else PDE.ratio = pj.tg.temp$ratio
			}
		} else {
			gene = matrix(data.count[gene.names == genenames, ], ncol=group.length)
			n = dim(gene)[1]
			if (n < 1) return(sprintf("Gene: %s has no data!", genenames))
			colnames(gene) = group.label
			rownames(gene) = NULL
			m.dt = melt(gene)
			offset.lib = rep(lib.size, each = n)
	
			ndt = data.frame(id = as.factor(m.dt[, 1]), g = as.factor(m.dt[, 2]), t = c(rep(1:group.1, each = n), rep(1:group.2, each = n)), count = m.dt[, 3])
	        if(offset){
	  	        fit1 = try(gssanova(count ~ t + g + t:g, offset = offset.lib, data = ndt, family = "nbinomial",  skip.iter=TRUE), silent = TRUE)
		        fit2 = try(gssanova(count ~ t + g, offset = offset.lib, data = ndt, family = "nbinomial", skip.iter=TRUE), silent = TRUE)
        	} else {
		        fit1 = try(gssanova(count ~ t + g + t:g, data = ndt, family = "nbinomial", skip.iter=TRUE), silent = TRUE)
          	    fit2 = try(gssanova(count ~ t + g, data = ndt, family = "nbinomial", skip.iter=TRUE), silent = TRUE)
        	}

			if (("try-error" %in% class(fit1))) NPDE.ratio = NA
			else {
				pj.tc.temp = try(project(fit1, c("t", "g")), silent = TRUE)
				if (("try-error" %in% class(pj.tc.temp))) NPDE.ratio = NA else NPDE.ratio = pj.tc.temp$ratio
			}
			if (("try-error" %in% class(fit2))){
				PDE.ratio = NA
				perm.nu = NA
				perm.eta = NA
			} else {
				perm.nu = fit2$nu 
				perm.eta = fit2$eta
				pj.tg.temp = try(project(fit2, c("t")), silent = TRUE)
				if (("try-error" %in% class(pj.tg.temp))) PDE.ratio = NA else PDE.ratio = pj.tg.temp$ratio
			}
		}
		list(gene.name=genenames, NPDE.ratio=NPDE.ratio, PDE.ratio=PDE.ratio, nu=perm.nu, fitted=perm.eta) 
	}

	stopCluster(cl)
	
	if (p.values) {
		if (iterations <= 1) {
			cat("Iterations must be a positive integer greater than 1.\n")
			return(cat("ERROR!"))
		}	
		
		pvalue = array(0, c(length(List), 1))
		perm.values = list()
	
		  for (i in 1 : length(List)){
			genenames = List[i]
			gene = matrix(data.count[gene.names == genenames, ], ncol=group.length)
			out = timeSeq.permutate(gene, genenames, gene.names, group.length, group.1, group.2, group.label, fits[[i]], 
									n.cores, iterations, lib.size, exon.length, exon.level, offset)					
			pvalue[i, ] = c(out[[1]])
			perm.values[[i]] = out[[2]]
		  }
		

		pvalue = data.frame(List, pvalue)
		colnames(pvalue) = c("genename", "NPDE")
	} else pvalue = NULL


	count = numeric(length(List))
	max_length = 0
	for (i in c(1:length(List))) {
		count[i] = dim(matrix(data.count[gene.names == List[i], ], ncol=group.length))[1]
		if (count[i] > max_length) max_length = count[i]
	}

	table = array(0, c(length(List), max_length, group.length))
	for (i in c(1 : length(List))) {
		m = unlist(data.count[gene.names == List[i], 1:group.length])
		m = t(matrix(m, c(count[i], group.length)))
		for (j in c(1:count[i])) table[i, j, ] = m[, j]
	}

    NPDE.ratio = unlist(lapply(fits, "[[", "NPDE.ratio"))
    PDE.ratio = unlist(lapply(fits, "[[", "PDE.ratio"))
	###Sort by ratios
	sorted = timeSeq.sort(List, NPDE.ratio, PDE.ratio, table, count)

	out = list(sorted = sorted, 
	           count = count,
	           NPDE.ratio = NPDE.ratio,
               PDE.ratio = PDE.ratio, 
               genenames = List, 
               table = table, 
               data = data.count, 
               gene.names = gene.names, 
		       lib.size = lib.size, 
		       exon.length = exon.length,
		       group.label = group.label,
               group.length = group.length, 
               group1.length = group.1, 
               group2.length = group.2, 
               exon.level = exon.level,
               p.values = p.values,
               pvalue = pvalue
		       )
               
	class(out) = "timeSeq.obj"
	return(out)
}
