#	Fit negative binomial mixed effect model for each gene
#	Authors: Xiaoxiao Sun & Fan Gao
#	Created 25 Seq 2015. Last modified 10 July 2017 (more stable).

timeSeq = function(data.count, group.label, gene.names, exon.length = NULL, exon.level = TRUE) {
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
    NPDE.ratio = NULL
    PDE.ratio = NULL
    for(genenames in List){
          ind = which(List %in% genenames)
          if (exon.level){
             gene = matrix(data.count[gene.names == genenames, ], ncol = group.length)
             n = dim(gene)[1]
             if (n < 1) return(sprintf("Gene: %s has no data!", genenames))
             colnames(gene) = group.label
             rownames(gene) = NULL
             m.dt = melt(gene)
             if (!is.null(exon.length)){
                 el =  rep(exon.length[gene.names == genenames], group.length)
                 el = log(el)
                 if(sum(el<=0)){
                     el[which(el<=0)] = 0.1
                 }
                 ndt = data.frame(id = as.factor(m.dt[, 1]), g = as.factor(m.dt[, 2]), t = c(rep(1:group.1, each = n), rep(1:group.2, each = n)), leng = el, count = m.dt[, 3])
                 ran <- mkran(~1|id, ndt)
                 randm = list(z=el*ran$z, sigma=ran$sigma, init=ran$init)
                 fit1 = try(gssanova(count ~ t + g + t:g, data = ndt, family = "nbinomial", random = randm), silent = TRUE)
                 fit2 = try(gssanova(count ~ t + g, data = ndt, family = "nbinomial", random = randm), silent = TRUE)
            } else {
                ndt <- data.frame(id = as.factor(m.dt[, 1]), g = as.factor(m.dt[, 2]), t = c(rep(1:group.1, each = n), rep(1:group.2, each = n)),
                count = m.dt[, 3])
                fit1 = try(gssanova(count ~ t + g + t:g, data = ndt, family = "nbinomial", random = ~1 | id), silent = TRUE)
                fit2 = try(gssanova(count ~ t + g, data = ndt, family = "nbinomial", random = ~1 | id), silent = TRUE)
            }
            if (("try-error" %in% class(fit1))) NPDE.ratio[ind] = NA
            else {
                pj.tc.temp = try(project(fit1, c("t", "g")), silent = TRUE)
                if (("try-error" %in% class(pj.tc.temp))) NPDE.ratio[ind] = NA else NPDE.ratio[ind] = pj.tc.temp$ratio
            }
            if (("try-error" %in% class(fit2))) {
                PDE.ratio[ind] = NA
            } else {
                pj.tg.temp = try(project(fit2, c("t")), silent = TRUE)
                if (("try-error" %in% class(pj.tg.temp))) PDE.ratio[ind] = NA else PDE.ratio[ind] = pj.tg.temp$ratio
            }
          } else {
            gene = matrix(data.count[gene.names == genenames, ], ncol=group.length)
            n = dim(gene)[1]
            if (n < 1) return(sprintf("Gene: %s has no data!", genenames))
            colnames(gene) = group.label
            rownames(gene) = NULL
            m.dt = melt(gene)
            ndt = data.frame(id = as.factor(m.dt[, 1]), g = as.factor(m.dt[, 2]), t = c(rep(1:group.1, each = n), rep(1:group.2, each = n)), count = m.dt[, 3])
            fit1 = try(gssanova(count ~ t + g + t:g, data = ndt,family = "nbinomial"), silent = TRUE)
            fit2 = try(gssanova(count ~ t + g, data = ndt,family = "nbinomial"), silent = TRUE)
            if (("try-error" %in% class(fit1))) NPDE.ratio[ind] = NA
            else {
                pj.tc.temp = try(project(fit1, c("t", "g")), silent = TRUE)
            if (("try-error" %in% class(pj.tc.temp))) NPDE.ratio[ind] = NA else NPDE.ratio[ind] = pj.tc.temp$ratio
            }
            if (("try-error" %in% class(fit2))){
                PDE.ratio[ind] = NA
            } else {
                pj.tg.temp = try(project(fit2, c("t")), silent = TRUE)
            if (("try-error" %in% class(pj.tg.temp))) PDE.ratio[ind] = NA else PDE.ratio[ind] = pj.tg.temp$ratio
            }
         }
    }

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

	  sorted = timeSeq.sort(List, NPDE.ratio, PDE.ratio, table, count)
	  out = list(sorted = sorted,
	             count = count,
	             NPDE.ratio = NPDE.ratio,
               PDE.ratio = PDE.ratio,
               genenames = List,
               table = table,
               data = data.count,
               gene.names = gene.names,
		           exon.length = exon.length,
		           group.label = group.label,
               group.length = group.length,
               group1.length = group.1,
               group2.length = group.2,
               exon.level = exon.level)
	   class(out) = "timeSeq.obj"
	   return(out)
}
