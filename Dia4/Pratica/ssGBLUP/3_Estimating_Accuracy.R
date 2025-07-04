# Limpa o ambiente
rm(list = ls())

# Leitura dos EBVs do conjunto completo e dos cinco grupos de validação cruzada
ebvs_all <- read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/ssGBLUP/All/EBVs_All.txt", header = TRUE)
ebvs_cv1 <- read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/ssGBLUP/CV1/EBVs_cv1.txt", header = TRUE)
ebvs_cv2 <- read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/ssGBLUP/CV2/EBVs_cv2.txt", header = TRUE)
ebvs_cv3 <- read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/ssGBLUP/CV3/EBVs_cv3.txt", header = TRUE)
ebvs_cv4 <- read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/ssGBLUP/CV4/EBVs_cv4.txt", header = TRUE)
ebvs_cv5 <- read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/ssGBLUP/CV5/EBVs_cv5.txt", header = TRUE)

# Leitura do arquivo de fenótipos
all_pheno <- read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/Feno_geno.txt")

# Leitura e filtragem dos indivíduos de validação de cada grupo
cv1_val <- all_pheno[all_pheno$V1 %in% read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/ssGBLUP/CV1/Data_g1_ssg.txt")[, "V1"], ]
cv1_val <- cv1_val[read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/ssGBLUP/CV1/Data_g1_ssg.txt")[, "V6"] == 1, ]

cv2_val <- all_pheno[all_pheno$V1 %in% read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/ssGBLUP/CV2/Data_g2_ssg.txt")[, "V1"], ]
cv2_val <- cv2_val[read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/ssGBLUP/CV2/Data_g2_ssg.txt")[, "V6"] == 2, ]

cv3_val <- all_pheno[all_pheno$V1 %in% read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/ssGBLUP/CV3/Data_g3_ssg.txt")[, "V1"], ]
cv3_val <- cv3_val[read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/ssGBLUP/CV3/Data_g3_ssg.txt")[, "V6"] == 3, ]

cv4_val <- all_pheno[all_pheno$V1 %in% read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/ssGBLUP/CV4/Data_g4_ssg.txt")[, "V1"], ]
cv4_val <- cv4_val[read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/ssGBLUP/CV4/Data_g4_ssg.txt")[, "V6"] == 4, ]

cv5_val <- all_pheno[all_pheno$V1 %in% read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/ssGBLUP/CV5/Data_g5_ssg.txt")[, "V1"], ]
cv5_val <- cv5_val[read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/ssGBLUP/CV5/Data_g5_ssg.txt")[, "V6"] == 5, ]

# ------------------------------------------------------------------------------
# Função para cálculo de:
# - Acurácia 1: correlação entre EBV e fenótipo, ajustada pela herdabilidade
# - Acurácia 2: baseada na covariância entre EBVs completos e EBVs CV
# - Viés: coeficiente da regressão linear fenótipo ~ EBV (método LR)
#
# Referências:
# - Yoshida et al. (2018): https://doi.org/10.1534/g3.117.300499
# - Legarra & Reverter (2018): https://doi.org/10.1186/s12711-018-0426-6
# ------------------------------------------------------------------------------

# Parâmetros estimados(PBLUP)
h2 <- 0.38080
varg <- 23.875

# Função para calcular métricas por grupo
calcula_metricas <- function(cv_ebv, cv_val, grupo) {
  ebvs_merged <- merge(ebvs_all, cv_ebv, by = "Original_ID")
  colnames(ebvs_merged) <- c("Original_ID", "Id_all", "EBV_all", "SE_all",
                             "Id_cv", paste0("EBV_cv", grupo), paste0("SE_cv", grupo))
  
  colnames(cv_val) <- c("Original_ID", "Idade", "Peso", "SB", "DDM")
  dados_final <- merge(ebvs_merged, cv_val, by = "Original_ID")
  
  EBV_cv <- dados_final[[paste0("EBV_cv", grupo)]]
  
  cor_val <- cor(EBV_cv, dados_final$DDM)                # Correlação fenótipo ~ EBV
  cov_val <- cov(dados_final$EBV_all, EBV_cv)            # Covariância EBV_all vs EBV_cv
  acc1 <- cor_val / sqrt(h2)                             # Acurácia 1: Yoshida et al. (2018)
  acc2 <- sqrt(cov_val / varg)                           # Acurácia 2: Legarra & Reverter (2018)
  bias <- summary(lm(DDM ~ EBV_cv, data = dados_final))$coefficients[2, 1]  # Viés (slope)
  
  return(c(acc1, acc2, bias))
}

# Cálculo das métricas para cada grupo
metrics_cv1 <- calcula_metricas(ebvs_cv1, cv1_val, 1)
metrics_cv2 <- calcula_metricas(ebvs_cv2, cv2_val, 2)
metrics_cv3 <- calcula_metricas(ebvs_cv3, cv3_val, 3)
metrics_cv4 <- calcula_metricas(ebvs_cv4, cv4_val, 4)
metrics_cv5 <- calcula_metricas(ebvs_cv5, cv5_val, 5)

# Consolidando as métricas
all <- rbind(metrics_cv1, metrics_cv2, metrics_cv3, metrics_cv4, metrics_cv5)
mean_acc1 <- mean(all[,1]); sd_acc1 <- sd(all[,1])
mean_acc2 <- mean(all[,2]); sd_acc2 <- sd(all[,2])
mean_bias <- mean(all[,3]); sd_bias <- sd(all[,3])

# Criando a tabela final
final <- data.frame(
  CV_group = c("CV1", "CV2", "CV3", "CV4", "CV5", "Média", "Desvio"),
  Acurácia1 = c(all[,1], mean_acc1, sd_acc1),
  Acurácia2 = c(all[,2], mean_acc2, sd_acc2),
  Viés      = c(all[,3], mean_bias, sd_bias)
)

# Exibe a tabela
print(final)

# Salva em arquivo txt
write.table(final,
            file = "C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/ssGBLUP/Acc_ssGBLUP.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
