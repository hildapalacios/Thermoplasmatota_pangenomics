#leer matriz

library(tidyverse)
library(ape)


setwd("/home/arturo/FiloArqueaBBYs_KEGG/")
unfiltered_binary_pangenome_matrix <- read.csv("binary_pangenome_matrix.csv")


#arturo@Microbio:~/FiloArqueaBBYs_KEGG$ scp "Filogenia345.txt" hilda@beagle.local:/home/hilda/Descargas/ 



# 1) añadir taxonomia
# 2) criterios de inclusión para reducir redundancia
# 3) aqui va la chamba de get_homologues
# 4) script para matriz de presencia-ausencia
# 5) Matriz de distancia
# 6) Filogenia 


#### FILTRAR MATRIZ (NO SINGLETONES)

binary_pangenome_matrix_plus_sum <- unfiltered_binary_pangenome_matrix[,-1] %>%
  replace(is.na(.), 0) %>%
  mutate(sum= rowSums(.))
rownames(binary_pangenome_matrix_plus_sum)<- unfiltered_binary_pangenome_matrix[,1]

plot(hist(binary_pangenome_matrix_plus_sum$sum))

#Filtrar singletones

porcentaje=5


ngenomas= ncol(binary_pangenome_matrix_plus_sum)-1
num_gen_perc = (ngenomas*(porcentaje))/100

binary_pm_No_singletones<- binary_pangenome_matrix_plus_sum %>% 
  filter(sum > num_gen_perc)

binary_pm_No_singletones_No_Hypothetical <- binary_pm_No_singletones %>% 
  filter(!grepl("hypothetical", rownames(binary_pm_No_singletones)))

plot(hist(binary_pm_No_singletones$sum))


table(as.character(binary_pangenome_matrix_plus_sum$sum)) %>% sort()

binary_pangenome_matrix <- binary_pm_No_singletones_No_Hypothetical[,!(colnames(binary_pangenome_matrix)=="sum")]


#Estimación de matriz de distancias

Distance_matrix <- stats::dist(t(binary_pangenome_matrix))

Distance_matrix.tb <- as.data.frame(as.matrix(Distance_matrix)) 
colnames(Distance_matrix.tb) <- colnames(binary_pangenome_matrix)
rownames(Distance_matrix.tb) <- colnames(binary_pangenome_matrix)

Distance_matrix.tb <- Distance_matrix.tb[-1, -1]

Distance_matrix.tb <- Distance_matrix.tb[-grep("tmp", rownames(Distance_matrix.tb)), -grep("tmp", rownames(Distance_matrix.tb))]


NJ <- ape::njs(as.dist(Distance_matrix.tb))

ape::write.tree(NJ, "Filogenia_nosingle_nohypo.txt")



