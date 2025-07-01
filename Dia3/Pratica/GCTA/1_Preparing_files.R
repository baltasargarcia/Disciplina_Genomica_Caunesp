rm(list=ls())

# Lendo os fenótipos
setwd("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia3/Pratica/GCTA")
phen<-read.table("Data.txt")
head(phen)

# Mantendo apenas os animais desafiados que possuem genótipos (749)
ind_gen<-read.table("IPN.fam")
phen_gen<-phen[phen$V1 %in% ind_gen$V2,]
head(phen_gen)

# Preparando os arquivos de covariáveis e fenótipos
cov<-data.frame(phen_gen$V1, phen_gen$V1, phen_gen$V2, phen_gen$V3)
phen_bin<-data.frame(phen_gen$V1, phen_gen$V1, phen_gen$V4)
phen_ddm<-data.frame(phen_gen$V1, phen_gen$V1, phen_gen$V5)
write.table(cov, "Cov.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(phen_bin, "Phen_bin.txt", row.names = F, col.names = F, quote = F, sep = "\t")
write.table(phen_ddm, "Phen_ddm.txt", row.names = F, col.names = F, quote = F, sep = "\t")

# Preparando os arquivos de genótipos
# Mantendo apenas os genótipos dos animais com fenótipos
write.table(data.frame(phen_ddm[,1:2]), "Only_phen.txt", row.names = F, col.names = F, quote = F, sep = "\t")
system("plink --bfile IPN --chr-set 29 --keep Only_phen.txt --make-bed --out IPN_gen")

# Realizando controle de qualidade
# Call-rate dos SNPs ≥ 90%
system("plink --bfile IPN_gen --geno 0.1 --chr-set 29 --make-bed --out IPN_gen90")

# Frequência alélica mínima (MAF) ≥ 0.05
system("plink --bfile IPN_gen90 --maf 0.05 --chr-set 29 --make-bed --out IPN_gen90_maf05")

0.05/33181  # HWE corrigido por Bonferroni

# Equilíbrio de Hardy-Weinberg (HWE)
system("plink --bfile IPN_gen90_maf05 --hwe 1.506886e-06 --chr-set 29 --make-bed --out IPN_gen90_maf05_hwe")

# Call-rate dos indivíduos ≥ 90%
system("plink --bfile IPN_gen90_maf05_hwe --mind 0.1 --chr-set 29 --make-bed --out IPN_final")
