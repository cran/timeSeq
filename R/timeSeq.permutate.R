#	Permutation for Significance Testing (only for NPDE genes)
#	Authors: Fan Gao & Xiaoxiao Sun
#	Created 25 Seq 2015. Last modified 10 Oct 2016.

timeSeq.permutate = function(gene.initial, gene.name, gene.names, group.length, group1.length, group2.length, group.label, 
	ratio.initial, n.cores = NULL, iterations = 10, lib.size, exon.length, exon.level, offset) {
  
	if (iterations <= 1) {
		cat("Iterations must be a positive integer greater than 1.\n")
		return(cat("ERROR!"))
	}	
  
	desired.cores = n.cores
	if (length(desired.cores) == 0) {
	    max.cores = detectCores()
	    desired.cores = max.cores - 1		
	} 	
	
	cl.perm = makeCluster(desired.cores)
	registerDoParallel(cl.perm)	

    fitted = ratio.initial$fitted
    permnu = ratio.initial$nu

	perm.fits = foreach(k = 1 : iterations, .packages = c("reshape", "gss")) %dopar% {  
    	if (exon.level) {
            rnb.dt = rnbinom(length(fitted), prob=(exp(fitted))/(1+exp(fitted)), size=permnu)
			gene = gene.initial
			n = dim(gene)[1]
			if (n < 1) return(sprintf("Gene: %s has no data!", gene.name))
			colnames(gene) = group.label
			rownames(gene) = NULL
			m.dt = melt(gene)
			offset.lib = rep(lib.size, each = n)
			if (!is.null(exon.length)) {
				el =  rep(exon.length[gene.names == gene.name], group.length)
				el = log(el)
                if(sum(el<=0)){
                	el[which(el<=0)] = 0.1
                }

				ndt = data.frame(id = as.factor(m.dt[, 1]), g = as.factor(m.dt[, 2]), t = c(rep(1:group1.length, each = n), rep(1:group2.length, each = n)), leng = el, count = rnb.dt)
				ran <- mkran(~1|id, ndt)
                randm = list(z=el*ran$z, sigma=ran$sigma, init=ran$init)

		        if(offset){
        			fit1 = try(gssanova(count ~ t + g + t:g, offset = offset.lib, data = ndt, family = "nbinomial", random = randm, skip.iter=TRUE), silent = TRUE)
		        } else {
          			fit1 = try(gssanova(count ~ t + g + t:g, data = ndt, family = "nbinomial", random = randm, skip.iter=TRUE), silent = TRUE)
        		}
			
			} else {
				ndt <- data.frame(id = as.factor(m.dt[, 1]), g = as.factor(m.dt[, 2]), t = c(rep(1:group1.length, each = n), rep(1:group2.length, each = n)), 
				count =rnb.dt)
		        if(offset){
        			fit1 = try(gssanova(count ~ t + g + t:g, offset = offset.lib, data = ndt, family = "nbinomial", skip.iter=TRUE), silent = TRUE)
        		} else {
		            fit1 = try(gssanova(count ~ t + g + t:g, data = ndt, family = "nbinomial", skip.iter=TRUE), silent = TRUE)
        		}				
			}

			if (("try-error" %in% class(fit1))) NPDE.ratio = NA
			else {
				pj.tc.temp = try(project(fit1, c("t", "g")), silent = TRUE)
				if (("try-error" %in% class(pj.tc.temp))) NPDE.ratio = NA else NPDE.ratio = pj.tc.temp$ratio
			}
			
		} else {

            rnb.dt = rnbinom(length(fitted), prob=(exp(fitted))/(1+exp(fitted)), size=permnu) 
			gene = gene.initial
			n = dim(gene)[1]
			if (n < 1) return(sprintf("Gene: %s has no data!", gene.name))
			colnames(gene) = group.label
			rownames(gene) = NULL
			m.dt = melt(gene)
			offset.lib = rep(lib.size, each = n)
			ndt = data.frame(id = as.factor(m.dt[, 1]), g = as.factor(m.dt[, 2]), t = c(rep(1:group1.length, each = n), rep(1:group2.length, each = n)), count = rnb.dt)

	        if(offset){
	  	        fit1 = try(gssanova1(count ~ t + g + t:g, offset = offset.lib, data = ndt, family = "nbinomial", skip.iter=TRUE, nu=permnu, alpha=1.4), silent = TRUE)
        	} else {
		        fit1 = try(gssanova1(count ~ t + g + t:g, data = ndt, family = "nbinomial", skip.iter=TRUE, nu=permnu), silent = TRUE)
        	}

			if (("try-error" %in% class(fit1))) NPDE.ratio = NA
			else {
				pj.tc.temp = try(project(fit1, c("t", "g")), silent = TRUE)
				if (("try-error" %in% class(pj.tc.temp))) NPDE.ratio = NA else NPDE.ratio = pj.tc.temp$ratio
			}
		}	
		list(NPDE.perm = NPDE.ratio)
    }
  
	stopCluster(cl.perm)

    NPDE.state = ratio.initial$NPDE.ratio

    perm.fits = unlist(lapply(perm.fits, "[[", "NPDE.perm"))
    
    if ((!is.na(NPDE.state)) && (sum(!is.na(as.numeric(perm.fits))) >= 1)) 
		pvalue.NPDE = sum(as.numeric(perm.fits[!is.na(as.numeric(perm.fits))]) > NPDE.state) / sum(!is.na(as.numeric(perm.fits)))
    else pvalue.NPDE = NA
   	
    out = list(pvalue.NPDE, perm.fits)
    return(out)
}

