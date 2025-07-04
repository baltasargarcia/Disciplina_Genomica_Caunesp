# Remove todos os objetos da memória
rm(list = ls())

# Leitura do arquivo com as soluções "verdadeiras" (do conjunto completo)
setwd("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/PBLUP/All")
sol <- read.table("solutions", header = TRUE, row.names = NULL)
head(sol)
table(sol$trait.effect)

# Seleciona apenas os EBVs (effect == 3)
ebvs <- sol[sol$trait.effect == 3, ]
head(ebvs)
colnames(ebvs) <- c("Trait", "Effect", "Id_new", "EBV", "EBV_SE")

# Lê o arquivo renadd para obter os IDs originais
code_ids <- read.table("renadd03.ped", col.names = c(
  "Id_new", "Parent1", "Parent2", "3-NofParents", "YOB",
  "Par_genotype", "N_records", "N_progeny1", "N_progeny2", "Original_ID"
))
head(code_ids)

# Faz o merge com os IDs originais
ebvs_ids <- merge(ebvs, code_ids, by = "Id_new")
head(ebvs_ids)

# Seleciona colunas finais: Id_new, Original_ID, EBV e SE
ebvs_final <- ebvs_ids[, c("Id_new", "Original_ID", "EBV", "EBV_SE")]
head(ebvs_final)

# Exporta para um arquivo
write.table(ebvs_final, "EBVs_All.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

# --------------------------
# LOOP para os 5 grupos de validação cruzada (CV1 a CV5)
# --------------------------
rm(list = ls())
for(i in 1:5) {
  # Define diretório do CV
  setwd(paste0("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/PBLUP/CV", i))
  
  # Lê o arquivo de soluções do CV
  sol <- read.table("solutions", header = TRUE, row.names = NULL)
  table(sol$trait.effect)
  
  # Seleciona os EBVs (effect == 3)
  ebvs <- sol[sol$trait.effect == 3, ]
  colnames(ebvs) <- c("Trait", "Effect", "Id_new", "EBV", "EBV_SE")
  
  # Lê os IDs originais
  code_ids <- read.table("renadd03.ped", col.names = c(
    "Id_new", "Parent1", "Parent2", "3-NofParents", "YOB",
    "Par_genotype", "N_records", "N_progeny1", "N_progeny2", "Original_ID"
  ))
  
  # Junta EBVs com IDs originais
  ebvs_ids <- merge(ebvs, code_ids, by = "Id_new")
  
  # Seleciona colunas finais
  ebvs_final <- ebvs_ids[, c("Id_new", "Original_ID", "EBV", "EBV_SE")]
  
  # Exporta resultado para cada grupo
  write.table(ebvs_final,
              file = paste0("EBVs_cv", i, ".txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE)
}
