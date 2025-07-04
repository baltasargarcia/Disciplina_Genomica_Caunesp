# Limpa o ambiente de trabalho
rm(list = ls())

# Carrega o pacote necessário
require(dplyr)

# Define o diretório de trabalho
setwd("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica")

# Lê o arquivo com os fenótipos dos animais genotipados
feno <- read.table("Feno_geno.txt")
head(feno)

# Adiciona uma nova coluna com grupos de validação cruzada (CV)
# 5 grupos com aproximadamente 200 animais cada
# Proporção: 80% referência / 20% validação
set.seed(123456)
grupo <- sample(c(rep(1:5, 200)))  # Vetor com 1000 indivíduos
feno$V6 <- grupo                   # Adiciona como nova coluna no dataset
head(feno)

# Visualiza a distribuição dos grupos
plot(table(feno$V6))

# Salva o dataset com os grupos de CV adicionados
write.table(feno, "CV_groups.txt", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Cria cópias do conjunto para cada grupo de validação
val1 <- feno
val2 <- feno
val3 <- feno
val4 <- feno
val5 <- feno

# Mascara os fenótipos (coluna V5) do grupo de validação com o valor -999
# Esse valor é interpretado como "missing" nos programas de predição
val1$V5[val1$V6 == 1] <- -999
val2$V5[val2$V6 == 2] <- -999
val3$V5[val3$V6 == 3] <- -999
val4$V5[val4$V6 == 4] <- -999
val5$V5[val5$V6 == 5] <- -999

# Salva os conjuntos de dados com validação (genotipados apenas)
write.table(val1, "Data_g1.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ")
write.table(val2, "Data_g2.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ")
write.table(val3, "Data_g3.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ")
write.table(val4, "Data_g4.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ")
write.table(val5, "Data_g5.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ")

# Prepara os conjuntos de dados para ssGBLUP, incluindo indivíduos sem genótipos
feno_no_geno <- read.table("Feno_no_geno.txt")
feno_no_geno$V6 <- 0  # Grupo 0 para não genotipados

# Junta os dados genotipados + não genotipados para cada grupo de validação
val1_ssg <- rbind(val1, feno_no_geno)
val2_ssg <- rbind(val2, feno_no_geno)
val3_ssg <- rbind(val3, feno_no_geno)
val4_ssg <- rbind(val4, feno_no_geno)
val5_ssg <- rbind(val5, feno_no_geno)

# Salva os arquivos para ssGBLUP
write.table(val1_ssg, "Data_g1_ssg.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ")
write.table(val2_ssg, "Data_g2_ssg.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ")
write.table(val3_ssg, "Data_g3_ssg.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ")
write.table(val4_ssg, "Data_g4_ssg.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ")
write.table(val5_ssg, "Data_g5_ssg.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ")
