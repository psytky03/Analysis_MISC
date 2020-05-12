# testGeneOverlap
# a function for simple hypergeometric enrichement test
# inputs are vector of observed IDs, reference pathway IDs, and universe IDs.
# return df with p.value and odd.ratio
testGeneOverlap <- function(observedIDs, pathwayIDs, universeIDs, side = "greater"){
  obs_size <- length(observedIDs)
  pat_size <- length(pathwayIDs)
  side <- "less"
  obs_pat_size <- length(base::union(observedIDs, pathwayIDs))
  univ_size <- length(universeIDs)
  shared_size <- length(base::intersect(observedIDs, pathwayIDs))
  
  cont.tbl <- matrix(c(univ_size - obs_pat_size, 
                       pat_size - shared_size, 
                       obs_size - shared_size,
                       shared_size), 
                     ncol=2)
  colnames(cont.tbl) <- c('notPathway', 'inPathway')
  rownames(cont.tbl) <- c('notObs', 'inObs')
  
  print(cont.tbl)
  res.fisher <- try(fisher.test(cont.tbl, alternative= side), silent=TRUE)
  res.fisher
  
  res <- data.frame(pval =double(), odds.ratio=double())
  
  if(is.list(res.fisher)) {
    res[1,]$odds.ratio <- res.fisher$estimate
    res[1,]$pval <- res.fisher$p.value
  } else {
    res[1,]$odds.ratio <- .0
    res[1,]$pval <- 1.
  }
  return(res)
}

# From https://rdrr.io/bioc/GeneOverlap/man/GeneOverlap.html
# The problem of gene overlap testing can be described by a hypergeometric distribution 
# where one gene list A defines the number of white balls in the urn and the 
# other gene list B defines the number of white balls in the draw. 
# Assume the total number of genes is n, the number of genes in A is a and the number of genes in B is b. 
# If the intersection between A and B is t, the probability density of seeing t can be calculated as:
#   
#   dhyper(t, a, n - a, b)
# 
# without loss of generality, we can assume b <= a. So the largest possible value for t is b. 
# Therefore, the p-value of seeing intersection t is:
#   
#   sum(dhyper(t:b, a, n - a, b))
# 
# The Fisher's exact test forms this problem slightly different but its calculation is also
# based on the hypergeometric distribution. It starts by constructing a contingency table:
# 
# matrix(c(n - union(A,B), setdiff(A,B), setdiff(B,A), intersect(A,B)), nrow=2)
# 
# It therefore tests the independence between A and B and is conceptually more straightforward. 
# It is better to illustrate a concept using some example.  
# Let's assume we have a genome of size 200 and two gene lists with 70 and 30 genes each. 
# If the intersection between the two is 10, the hypergeometric way to calculate the p-value is:
# 
# fisher.test(matrix(c(110, 60, 20, 10), nrow=2), alternative="greater")