rm(list = ls())

# Instalar o pacote caso necessário
# install.packages("qqman")
library(qqman)

# Definindo o diretório de trabalho
setwd("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia3/Pratica/GCTA")

########################
# GWAS usando GCTA
########################

# 1. Criando a matriz de relação genômica (GRM)
system("gcta64 --bfile IPN_final --autosome-num 29 --make-grm --out GRM_IPN")

# 2. Gerando GRM final com limiar de relacionamento (bK)
system("gcta64 --grm GRM_IPN --autosome-num 29 --make-bK 0.05 --out GRM_final")

# 3. Estimando variância genética com REML
system("gcta64 --grm GRM_final --autosome-num 29 --pheno Phen_ddm.txt --qcovar Cov.txt --reml --out Reml_challenge")

# 4. Executando GWAS com MLMA-LOCO (leave-one-chromosome-out)
system("gcta64 --bfile IPN_final --autosome-num 29 --pheno Phen_ddm.txt --qcovar Cov.txt --grm GRM_final --mlma-loco --out GWAS_ddm")

# 5. Lendo os resultados do GWAS
MLMA_IPN <- read.table("GWAS_ddm.loco.mlma", header = TRUE)
head(MLMA_IPN)

# 6. Preparando o dataframe para os gráficos
m_plot <- data.frame(
  SNP_ID = MLMA_IPN$SNP,
  CHR = MLMA_IPN$Chr,
  BPP = MLMA_IPN$bp,
  P_VALUE = MLMA_IPN$p,
  FREQ = MLMA_IPN$Freq
)

########################
# Gráficos: Manhattan e QQ plot
########################

# Gráfico de Manhattan
manhattan(m_plot,
          chr = "CHR",
          bp = "BPP",
          snp = "SNP_ID",
          p = "P_VALUE",
          ylim = c(-0.25, 8),
          suggestiveline = -log10((0.05)/(31812/44)),  # linha sugestiva
          genomewideline = -log10(0.05/31812),         # linha genômica
          col = c("black", "red", "gray"),
          main = "GWAS for IPN resistance (Day of death)"
)

# Gráfico QQ
qq(MLMA_IPN$p,
   main = "QQ plot for GWAS (Day of death)",
   col = "darkblue",
   pch = 20,
   cex = 1.1)
abline(0, 1, col = "red", lty = 2)  # linha de referência

##############################
# Cálculo da Variância Explicada por SNP
##############################

# Calculando a variância explicada (%)
MLMA_IPN$var <- 100 * (MLMA_IPN[, 7]^2) * (2 * MLMA_IPN[, 6] * (1 - MLMA_IPN[, 6])) /
  sum((MLMA_IPN[, 7]^2) * (2 * MLMA_IPN[, 6] * (1 - MLMA_IPN[, 6])), na.rm = TRUE)

# Soma total da variância explicada
sum(MLMA_IPN$var, na.rm = TRUE)

# Reorganizando os dados por ordem de cromossomo e posição
MLMA_IPN2 <- MLMA_IPN[order(MLMA_IPN$Chr), ]
MLMA_IPN2$row <- 1:nrow(MLMA_IPN2)

# Ordenando pela variância explicada (maior para menor)
MLMA_IPN3 <- MLMA_IPN2[order(MLMA_IPN2$var, decreasing = TRUE), ]
head(MLMA_IPN3, 10)

# Filtrando SNPs significativos com p-valor < 10^-4.16
sig_snps_pvalor <- MLMA_IPN3[-log10(MLMA_IPN3$p) >= 4.16, ]

# Exportando os SNPs significativos
write.table(sig_snps_pvalor, "Sig_SNPs.txt", quote = FALSE, col.names = TRUE, row.names = FALSE)

# Anotação manual: investigar gene senp5 em Chr5 (posição ~86.195.983 ± 200.000) no genoma Omyk_1.0
# Referência: https://doi.org/10.1534/g3.119.400463
