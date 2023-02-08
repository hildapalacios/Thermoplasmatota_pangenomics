PosArgs <- as.character(commandArgs(trailingOnly = TRUE)) #Reading positional argument 



salida_gethoms <- PosArgs[1]
archivo_especies <- PosArgs[2]

clusterverse <- list.files(path = salida_gethoms, pattern=".faa")
genomeverse <- readLines(archivo_especies)

pan_matrix <- matrix(nrow = length(clusterverse), ncol = length(genomeverse))
rownames(pan_matrix) <- clusterverse
colnames(pan_matrix) <- genomeverse

for (i in 1:nrow(pan_matrix)) {
  cur_clust <- paste(salida_gethoms,clusterverse[i],sep = "/")
  cur_headers <- system(paste("grep '>'",cur_clust),intern = TRUE)
  for (j in 1:ncol(pan_matrix)) {
    cur_genome <- genomeverse[j]
    if(any(base::grepl(cur_genome,cur_headers))){
      pan_matrix[i,j] <- 1
    }else{
      pan_matrix[i,j] <- 0
    }
  }
  print(nrow(pan_matrix)-i)
}

write.csv(pan_matrix, file = "binary_pangenome_matrix.csv")


# Calcular matriz de distancia 
Entrada_matriz <- t(pan_matrix)

#Distance_matrix

library(ape)
library(ade4)
   
