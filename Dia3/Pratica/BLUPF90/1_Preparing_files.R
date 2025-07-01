rm(list=ls())

#install.packages("gdata")
require(gdata)

# Lendo os fenótipos
setwd("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia3/Pratica/BLUPF90")
phen<-read.table("Data.txt")
head(phen)

# Mantendo apenas os animais desafiados que possuem genótipos (749)
ind_gen<-read.table("IPN.fam")
phen_gen<-phen[phen$V1 %in% ind_gen$V2,]
head(phen_gen)

# Preparando o arquivo de fenótipos (espaco entre colunas)
write.fwf(phen_gen, "Fenotipos.txt", rownames = F, colnames = F, quote = F, sep = " ")

# Preparando os arquivos de genótipos
# Mantendo apenas os genótipos dos animais com fenótipos
write.table(data.frame(phen_gen[,c(1,1)]), "Only_phen.txt", row.names = F, col.names = F, quote = F, sep = "\t")
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
system("plink --bfile IPN_gen90_maf05_hwe --mind 0.1 --chr-set 29 --make-bed --recodeA --out IPN_final")

#Preparando o arquivo de genotipos para o formato aceito pelo BLUPF90
geno<-read.table("IPN_final.raw",skip=1,sep=" ")

#Selectiong only ID and genotypes (removing unnecessary columns)
geno2<-geno[,c(2,7:ncol(geno))]
str(geno2)

#Replacing NA by 5 (BLUPF90 reads 5 as missing genotype)
geno2[is.na(geno2)] <- 5

#Separating IDs from genotypes
#Only Ids
ids<-data.frame(geno2[,1])

#Only genotypes
geno3<-geno2[,c(2:ncol(geno2))]

#Removing space between SNPs 
my_cols<-colnames(geno3)
geno4<-do.call(paste, c(geno3[my_cols], sep = ""))

#Joining IDs and Genotypes again
geno_final<-cbind(ids,geno4)

write.fwf(geno_final,"Genotipos_final.txt",sep=" ",quote = F,colnames = F)

#Preparando mapa
mapa<-read.table("IPN_final.bim")
mapa2<-data.frame(1:nrow(mapa),mapa$V1,mapa$V4)
mapa2
colnames(mapa2)<-c("SNP_ID","CHR","POS")
write.fwf(mapa2,"Mapa.txt",sep=" ",quote = F,colnames = T)
