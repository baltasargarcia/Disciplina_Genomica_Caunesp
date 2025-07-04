# Limpa o ambiente de trabalho
rm(list = ls())

# Carrega os pacotes necessários
require(data.table)
require(gdata)

# Define o diretório de trabalho
setwd("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica")

# Lê o arquivo com os fenótipos e efeitos
data <- read.table("Desafio_SRS.txt", header = TRUE)
head(data)

# Prepara um novo dataframe com colunas selecionadas
data2 <- data.frame(data$ID, data$Age, data$Tag_weight, data$Mortal_bin, data$Mortal_DDM)
head(data2)

# Verifica se todos os animais do fenótipo estão presentes no pedigree
ped <- read.table("pedigree.txt")
head(ped)
presentes <- ped[ped$V1 %in% data2$data.ID, ]

# Gera dois conjuntos de fenótipos: com e sem genótipos disponíveis
fam <- read.table("SRS.fam")
feno_geno <- data2[data2$data.ID %in% fam$V1, ]
feno_no_geno <- data2[!data2$data.ID %in% fam$V1, ]

# Salva os arquivos com os fenótipos
write.fwf(feno_geno, "Feno_geno.txt", colnames = FALSE, rownames = FALSE, quote = FALSE, sep = " ")
write.fwf(feno_no_geno, "Feno_no_geno.txt", colnames = FALSE, rownames = FALSE, quote = FALSE, sep = " ")
write.fwf(data2, "Feno_All.txt", colnames = FALSE, rownames = FALSE, quote = FALSE, sep = " ")
write.fwf(ped, "Pedigree_final.txt", colnames = FALSE, rownames = FALSE, quote = FALSE, sep = " ")

# Prepara o arquivo de genótipos com controle de qualidade usando PLINK
# 1. Filtro por taxa de chamada dos SNPs (call rate >= 0.9)
system("plink --bfile SRS --geno 0.1 --chr-set 29 --make-bed --out SRS_gen90")

# 2. Filtro por frequência alélica mínima (MAF >= 0.01)
system("plink --bfile SRS_gen90 --maf 0.01 --chr-set 29 --make-bed --out SRS_gen90_maf01")

# 3. Filtro por equilíbrio de Hardy-Weinberg (HWE, p > 1.49e-06)
# (Valor obtido por 0.05 / número de SNPs: 0.05 / 33542)
system("plink --bfile SRS_gen90_maf01 --hwe 1.490668e-06 --chr-set 29 --make-bed --out SRS_gen90_maf01_hwe")

# 4. Filtro por taxa de chamada dos indivíduos (call rate >= 0.9)
system("plink --bfile SRS_gen90_maf01_hwe --mind 0.1 --chr-set 29 --make-bed --out SRS_final")

# Converte genótipos de ACTG para notação 012 (PLINK recodeA)
system("plink --bfile SRS_final --chr-set 29 --recodeA --out SRS_geno")

# Lê o arquivo de genótipos recodificado
geno <- read.table("SRS_geno.raw", skip = 1, sep = " ")

# Seleciona apenas o ID e as colunas de genótipos (removendo colunas auxiliares)
geno2 <- geno[, c(2, 7:ncol(geno))]

# Substitui valores ausentes (NA) por 5 (BLUPF90 interpreta 5 como dado ausente)
geno2[is.na(geno2)] <- 5

# Separa IDs e genótipos
ids <- data.frame(geno2[, 1])                # Apenas IDs
geno3 <- geno2[, 2:ncol(geno2)]              # Apenas genótipos

# Junta os genótipos em uma única string por indivíduo
geno4 <- do.call(paste, c(geno3, sep = ""))

# Reagrupa os IDs com os genótipos concatenados
geno_final <- cbind(ids, geno4)

# Salva o arquivo final de genótipos em colunas de largura fixa
write.fwf(geno_final, "Genotipos_final.txt", quote = FALSE, sep = " ", colnames = FALSE)
