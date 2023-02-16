
library(ape)

as.dist(Disimil_Matrix)   # Cambiar Disimil_Matrix por la matriz de distancia 

NJ <- ape::nj(as.matrix(Disimil_Matrix))


f <- function(njs) njs(as.dist(Disimil_Matrix))
myBoots <- boot.phylo(NJ, Disimil_Matrix, f, B= 100)


plot.phylo(NJ, type = "p", use.edge.length = TRUE, edge.color = "turquoise4", font = 3, root.edge = TRUE)
nodelabels(myBoots, cex = .5)
