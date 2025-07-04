# Clase practica:Control de calidad de genotipos y análisis de poblaciones

#Limpieza del ambiente de trabajo
rm(list=ls())

# Camino al directorio de trabajo

setwd("C:/Users/vitoo/OneDrive/Documentos/Disciplina_Genômica/Dia2B")
getwd()


### Archivos de genotipo de tambaqui compuestos por 600 individuos y 192 SNPs (9 familias)

### Importante: el archivo ejecutable plink.exe y los archivos de genotipos deben estar dentro del directorio de trabajo.



##Parte 1: Control de calidad de genotipos


#Filtro_1 - elimina individuos con más del 10% de datos faltantes.


system("plink --bfile genotipos_poblaciones --aec --mind 0.1 --out genotipos_poblaciones_filtered_mind --make-bed")



#Filtro_2 - elimina SNPs que fallan en más del 10% de los individuos.


system("plink --bfile genotipos_poblaciones_filtered_mind --aec --geno 0.1 --out genotipos_poblaciones_filtered_geno --make-bed")



#Filtro_3 - elimina SNPs cuya frecuencia del alelo menor (Minor Allele Frequency) < 5%


system("plink --bfile genotipos_poblaciones_filtered_geno --aec --maf 0.05 --out genotipos_poblaciones_filtered_maf --make-bed")





#Filtrar HWE por familia 

#Capturar el ID de la familia (columna 1) y el ID del individuo (columna 2) de cada archivo .fam y crear un archivo (ids_*.txt).

fam <- read.table("genotipos_poblaciones_filtered_maf.fam", header=FALSE)
pops <- unique(fam$V1)
for (pop in pops) {
  subset <- fam[fam$V1 == pop, c(1,2)]
  write.table(subset, file=paste0("ids_", pop, ".txt"), 
              row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")
}


#Utilizo estos archivos para realizar el filtro de HWE por familia, mediante la opción --keep.

# Captura archivos de IDs generados previamente.
id_files <- list.files(pattern = "^ids_.*\\.txt$")
#Para cada archivo de IDs, generar y ejecutar el comando de PLINK.
for (ids_file in id_files) {
  pop <- sub("ids_(.*)\\.txt", "\\1", ids_file)
  pop <- trimws(gsub("[[:space:]]", "_", pop))
  out_prefix <- paste0("genotipos_", pop, "_hwe_filtered")
  
  cmd <- paste("plink --bfile genotipos_poblaciones_filtered_maf --keep", ids_file, "--hwe 1e-6 --make-bed --out", out_prefix)
  cat(cmd, "\n")
  system(cmd)
}



##Usar exactamente el mismo panel de SNPs para los análisis de diversidad y estructura.

# 1️ ler todos archivos bim   
bim_files <- list.files(pattern = "^genotipos_.*_hwe_filtered\\.bim$")
snp_lists <- lapply(bim_files, function(file) {
  bim <- read.table(file, header=FALSE)
  bim$V2  # coluna dos IDs de SNP
})

# 2️ Encontrar SNPs comuns a todos
snps_comuns <- Reduce(intersect, snp_lists)

# 3️ Guardar lista de SNPs comunes
delete.na <- function(x) x[!is.na(x)]
snps_comuns <- delete.na(snps_comuns)
write.table(snps_comuns, file="snps_comunes.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)



# Archivo de SNPs comunes generado previamente: snps_comunes.txt

# Carrega archivos de genótipos filtrados por HWE
bfiles <- list.files(pattern = "^genotipos_.*_hwe_filtered\\.bed$")

# Filtra los SNPs en común de cada familia mediante la opción --extract.
for (bed in bfiles) {
  prefix <- sub("\\.bed", "", bed)
  cmd <- paste("plink --bfile", prefix, "--extract snps_comunes.txt --make-bed --out", paste0(prefix, "_final"))
  cat("Executando:", cmd, "\n")
  system(cmd)
}


###Generar para cada família individualmente

##Número de loci polimórficos → Número de loci polimórficos

##Média de MAF (frequência alélica menor) → Media del MAF (frecuencia alélica menor)

##Heterozigosidade esperada (He) → Heterocigosidad esperada (He)

##Heterozigosidade observada (Ho) → Heterocigosidad observada (Ho)


##Etapas:
  
#1. Generar estadísticas básicas en PLINK para cada población (--freq y --hardy)


# Loop de 1 a 9 para cada família

for (i in 1:9) {
  base <- paste0("genotipos_", i, "_hwe_filtered_final")
  out <- paste0("genotipos_", i, "_stats")
  
  # Comando PLINK para generar freq e hardy
  cmd <- paste(
    "plink",
    "--bfile", base,
    "--freq",
    "--hardy",
    "--out", out
  )
  
   system(cmd)
  
  cat("Arquivos gerados para família", i, ":", out, ".frq e .hwe\n")
}


#Archivo .frq: Frecuencia alélica

  
#Columna Significado
#CHR	Número del cromosoma donde se encuentra el SNP. 
#SNP	Identificador del SNP (nombre o código único).
#A1	Alelo menor (minor allele), es decir, el alelo con menor frecuencia en la muestra.
#A2	Alelo mayor (major allele), es decir, el otro alelo observado en ese locus.
#MAF	Minor Allele Frequency (frecuencia del alelo menor)valor entre 0 y 0.5.	
#NCHROBS	Número de cromosomas observados (es decir, 2 veces el número de individuos genotipados).



#Archivo .hwe: Prueba de equilibrio Hardy-Weinberg

  
#Columna	Significado
#CHR	
#SNP
#TEST	"ALL" (se refiere a toda la muestra o subgrupo analizado).
#A1
#A2	
#GENO	AA = homocigotos para A1; AB = heterocigotos; BB = homocigotos para A2	
#O(HET)	Proporción observada de individuos heterocigotos (AB) para ese SNP.
#E(HET)	Proporción esperada de heterocigotos 
#P	Valor p de la prueba de Hardy-Weinberg.


##2. Script en R para procesar los archivos de salida .frq y .hwe 


pops <- 1:9  

# Criar data.frame final
resultados <- data.frame(Populacao = character(), 
                         LociPolimorficos = integer(),
                         Media_MAF = numeric(), 
                         Media_He = numeric(), 
                         Media_Ho = numeric(),
                         stringsAsFactors = FALSE)

for (i in pops) {
  # Leitura de arquivos PLINK
  frq <- read.table(paste0("genotipos_", i, "_stats.frq"), header = TRUE)
  hwe <- read.table(paste0("genotipos_", i, "_stats.hwe"), header = TRUE, check.names = FALSE)
  
  # Remover SNPs fixos
  frq_filt <- frq[frq$MAF > 0 & frq$MAF < 1, ]
  loci_polimorficos <- nrow(frq_filt)
  media_maf <- mean(frq_filt$MAF, na.rm = TRUE)
  
  # Filtrar apenas TEST == "ALL"
  hwe_filt <- hwe[grepl("ALL", hwe$TEST), ]
  
  # Acessar colunas com nomes especiais
  col_ho <- grep("O\\(HET\\)", colnames(hwe_filt), value = TRUE)
  col_he <- grep("E\\(HET\\)", colnames(hwe_filt), value = TRUE)
  
  media_ho <- mean(hwe_filt[[col_ho]], na.rm = TRUE)
  media_he <- mean(hwe_filt[[col_he]], na.rm = TRUE)
  
  resultados <- rbind(resultados,
                      data.frame(Populacao = as.character(i),
                                 LociPolimorficos = loci_polimorficos,
                                 Media_MAF = media_maf,
                                 Media_He = media_he,
                                 Media_Ho = media_ho))
}

# Visualizar resultados
print(resultados)

write.csv(resultados, file = "tabela_diversidade_por_familia.csv", row.names = FALSE)


#Gráficos comparativos

#install.packages("ggplot2")
#install.packages('reshape2')

library(ggplot2)
library(reshape2)


# Gráficas comparativas

# Convertir a formato largo para ggplot (praticidad)
resultados_long <- melt(resultados, id.vars = "Populacao")

# Gráfico de barras

jpeg(filename="graficas_comparativas.jpg",
     width=11,
     height=9,units = "in",res = 600)

ggplot(resultados_long, aes(x = Populacao, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = "Diversidade genética por família",
       x = "Família", y = "Valor") +
  theme_minimal()

dev.off()




###Corregir los valores de p de Hardy-Weinberg por FDR (Benjamini-Hochberg)

###Comparar el número de loci dentro y fuera del equilibrio en cada familia


pops <- 1:9

# Tabla resumen por población
resumen_hwe <- data.frame(Poblacion = character(),
                          SNPs_totales = integer(),
                          SNPs_en_equilibrio = integer(),
                          SNPs_fuera_equilibrio = integer(),
                          stringsAsFactors = FALSE)

for (i in pops) {
  hwe <- read.table(paste0("genotipos_", i, "_stats.hwe"), header = TRUE, check.names = FALSE)
  
  # Filtrar sólo las líneas con TEST = \"ALL\" o \"ALL(NP)\"
  hwe_filtrado <- hwe[grepl("ALL", hwe$TEST), ]
  
  # Corregir los p-valores con el método de FDR (Benjamini-Hochberg)
  hwe_filtrado$FDR <- p.adjust(hwe_filtrado$P, method = "BH")
  
  # Contar SNPs
  total <- nrow(hwe_filtrado)
  en_eq <- sum(hwe_filtrado$FDR > 0.05, na.rm = TRUE)
  fuera_eq <- sum(hwe_filtrado$FDR <= 0.05, na.rm = TRUE)
  
  resumen_hwe <- rbind(resumen_hwe,
                       data.frame(Poblacion = as.character(i),
                                  SNPs_totales = total,
                                  SNPs_en_equilibrio = en_eq,
                                  SNPs_fuera_equilibrio = fuera_eq))
}

# Ver resultados
print(resumen_hwe)


write.csv(resumen_hwe, "resumen_HWE_FDR.csv", row.names = FALSE)


resumen_largo <- melt(resumen_hwe[, c("Poblacion", "SNPs_en_equilibrio", "SNPs_fuera_equilibrio")],
                      id.vars = "Poblacion")


jpeg(filename="graficos_comparativos_HWE.jpg",
     width=11,
     height=9,units = "in",res = 600)

ggplot(resumen_largo, aes(x = Poblacion, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "SNPs en o fuera del equilibrio de Hardy-Weinberg (FDR ≤ 0.05)",
       x = "Familia", y = "Número de SNPs", fill = "Categoría") +
  theme_minimal()

dev.off()




##Calcular el FIS con el paquete hierfstat.

#Realizar bootstrap sobre los loci para estimar intervalos de confianza y significancia.

#Construir una tabla con el FIS observado, el intervalo de confianza y si el valor es significativo o no (basado en el IC del 95%).

#Limpieza del ambiente de trabajo
rm(list=ls())

#Instalar paquetes (si es necesario):
  
 
#install.packages("hierfstat")

#install.packages("adegenet")

library(hierfstat)

library(adegenet)


#Passo 1: unir archivos de las familias y obtener archivos ped/map con Plink

# Crear el archivo de lista con los nombres base de los archivos .bed/.bim/.fam
familias <- 1:9
nombres <- paste0("genotipos_", familias, "_hwe_filtered_final")
writeLines(nombres, "lista_familias.txt")

# Combinar todos los archivos usando PLINK
comando_merge <- paste(
  "plink",
  "--merge-list lista_familias.txt",
  "--make-bed",
  "--out genotipos_combinados"
)
cat("Ejecutando:", comando_merge, "\n")
system(comando_merge)

# Generar archivos PED y MAP desde el archivo combinado
comando_ped <- paste(
  "plink",
  "--bfile genotipos_combinados",
  "--recode",
  "--out genotipos_combinados"
)
cat("Ejecutando:", comando_ped, "\n")
system(comando_ped)

cat("✅ Proceso completo. Archivos combinados y recodificados generados con éxito.\n")



#ped to dat 


# Cargar las bibliotecas necesarias
library(data.table)
library(stringr)
library(adegenet)

# Archivo base
prefijo <- "genotipos_combinados"

# Leer archivos de PLINK
fam <- fread(paste0(prefijo, ".fam"), header = FALSE)
bim <- fread(paste0(prefijo, ".bim"), header = FALSE)
ped <- fread(paste0(prefijo, ".ped"), header = FALSE)

# Nombres de los loci
loci <- bim$V2
n_loci <- length(loci)
n_ind <- nrow(ped)

# Genotipos: columnas a partir de la 7
geno_raw <- ped[, 7:ncol(ped)]

# Tabla de conversión de alelos
alelos_unicos <- unique(unlist(geno_raw))
alelos_unicos <- alelos_unicos[alelos_unicos != "0"]
tabla_conversion <- data.frame(
  letra = sort(na.omit(unique(alelos_unicos))),
  codigo = sprintf("%02d", 1:length(unique(alelos_unicos)))
)
tabla_conversion <- rbind(data.frame(letra = "0", codigo = "00"), tabla_conversion)

# Función para recodificar
recodificar_alelo <- function(alelo) {
  cod <- tabla_conversion$codigo[match(alelo, tabla_conversion$letra)]
  ifelse(is.na(cod), "00", cod)
}

# Construir matriz de genótipos codificados
genos_codificados <- matrix(nrow = n_ind, ncol = n_loci)

for (i in 1:n_loci) {
  a1 <- as.character(geno_raw[[2*i - 1]])
  a2 <- as.character(geno_raw[[2*i]])
  cod1 <- recodificar_alelo(a1)
  cod2 <- recodificar_alelo(a2)
  genos_codificados[, i] <- paste0(cod1, cod2)
}

# Identificador da población
poblacion <- fam$V1
matriz_final <- cbind(poblacion, genos_codificados)

# Calcular el número de poblaciones
n_pops <- length(unique(poblacion))

# Archivo de saída
archivo_salida <- "genotipos_fstat_convertido.dat"
con <- file(archivo_salida, open = "wt")

# Linea 1: n_pops, n_ind, n_loci, 3
writeLines(paste(n_pops, n_ind, n_loci, 3), con)

# Lineas 2 até n_loci+1: nomes dos loci
writeLines(loci, con)

# Genotipos de los indivíduos
apply(matriz_final, 1, function(x) writeLines(paste(x, collapse = " "), con))

# Cerrar
close(con)

# Mostrar tabla de conversión
cat("✅ Arquivo .dat gerado corretamente:", archivo_salida, "\n")
cat("📋 Alelos convertidos:\n")
print(tabla_conversion)


##check archivo dat


# Leer el archivo .dat
archivo <- "genotipos_fstat_convertido.dat"
lineas <- readLines(archivo)

# 🔸 Metainformación (linha 1 tem 4 valores: npop, nind, nloci, 3)
meta <- str_split(lineas[1], "\\s+")[[1]]
n_pops <- as.numeric(meta[1])
n_ind <- as.numeric(meta[2])
n_loci <- as.numeric(meta[3])

# 🔸 Loci (linhas 2 até 1+n_loci)
loci <- lineas[2:(1 + n_loci)]

# 🔸 Genótipos (demais linhas)
genotipos <- lineas[(2 + n_loci):length(lineas)]

# 🔸 Separar por espaço
datos_split <- str_split(genotipos, "\\s+")

# ⚠️ Verificação de número de colunas
comprimentos <- sapply(datos_split, length)
if (any(comprimentos != (1 + n_loci))) {
  stop("❌ Número de colunas inconsistente com loci + coluna de população")
}

# 🔸 Construir matriz
matriz <- do.call(rbind, datos_split)
colnames(matriz) <- c("pop", loci)

datos <- as.data.frame(matriz)

# Convertir a factores para contar
datos$pop <- as.factor(datos$pop)
n_pops <- length(levels(datos$pop))

# Contar alelos únicos por locus
# Ex: "0103" → alelos: "01" y "03"
conteo_alelos <- sapply(datos[, -1], function(col) {
  col <- as.character(col)
  col <- col[col != "0000"]  # remover datos faltantes
  alelos1 <- substr(col, 1, 2)
  alelos2 <- substr(col, 3, 4)
  length(unique(c(alelos1, alelos2)))
})

# Mostrar resumen
cat("📊 RESUMEN DEL ARCHIVO .DAT\n")
cat("🔹 Número de poblaciones:", n_pops, "\n")
cat("🔹 Número de individuos  :", n_ind, "\n")
cat("🔹 Número de loci        :", n_loci, "\n")
cat("🔹 Promedio de alelos por locus:", round(mean(conteo_alelos), 2), "\n")
cat("   (rango:", min(conteo_alelos), "-", max(conteo_alelos), ")\n")



###Executar hierfstat


# 🔸 Leer archivo .dat
archivo <- "genotipos_fstat_convertido.dat"
lineas <- readLines(archivo)

# 🔸 Extraer metadatos (linha 1: n_pops, n_ind, n_loci, 3)
meta <- str_split(lineas[1], "\\s+")[[1]]
n_pops <- as.numeric(meta[1])
n_ind <- as.numeric(meta[2])
n_loci <- as.numeric(meta[3])

# 🔸 Nombres de loci (linhas 2 até 1 + n_loci)
loci <- lineas[2:(1 + n_loci)]

# 🔸 Genótipos (demais linhas)
genotipos <- lineas[(2 + n_loci):length(lineas)]

# 🔸 Separar por espaço
datos_split <- str_split(genotipos, "\\s+")

# 🔸 Validar número de colunas (1 coluna pop + n_loci loci)
comprimentos <- sapply(datos_split, length)
if (any(comprimentos != (1 + n_loci))) {
  stop("❌ Número de colunas inconsistente com loci + coluna de população")
}

# 🔸 Construir matriz
matriz <- do.call(rbind, datos_split)
colnames(matriz) <- c("pop", loci)


# Convertir a data.frame con valores numéricos
datos <- as.data.frame(matriz)
datos[, -1] <- lapply(datos[, -1], as.numeric)
datos$pop <- as.factor(datos$pop)

# Estadísticas básicas
res <- basic.stats(datos)
res

# Bootstrap para IC de FIS por población
set.seed(123)
boot_fis <- boot.ppfis(datos, nboot = 1000, quant = c(0.05, 0.95))

# Combinar resultados
tabla_resultados <- data.frame(
  Poblacion = as.character(unique(datos$pop)),
  FIS_observado = round(colMeans(res$Fis, na.rm = TRUE), 4),
  IC_inferior = round(boot_fis$fis.ci$ll, 4),
  IC_superior = round(boot_fis$fis.ci$hl, 4),
  Significativo = NA
)

# Determinar significancia
tabla_resultados$Significativo <- with(tabla_resultados, ifelse(
  IC_inferior > 0, "Positivo (*)",
  ifelse(IC_superior < 0, "Negativo (*)", "No significativo")
))

# Mostrar resultados
print(tabla_resultados)

write.csv(tabla_resultados, file = "fis_por_familia.csv", row.names = FALSE)

cat("✅ Tabla de FIS por familia guardada como 'fis_por_familia.csv'\n")



##PCA colorido por población

#Paso 1: Cargar los paquetes

# Instalar paquetes si es necesario
#if (!require("SNPRelate")) install.packages("BiocManager"); BiocManager::install("SNPRelate")

library(SNPRelate)

bed.fn <- "genotipos_combinados.bed"
bim.fn <- "genotipos_combinados.bim"
fam.fn <- "genotipos_combinados.fam"
gds.fn <- "genotipos_combinados.gds"



snpgdsBED2GDS(
  bed.fn = "genotipos_combinados.bed",
  bim.fn = "genotipos_combinados.bim",
  fam.fn = "genotipos_combinados.fam",
  out.gdsfn = file.path(tempdir(), "genotipos_combinados.gds"),
  verbose = TRUE
)

# Cargar la librería

library(ggplot2)

# Ruta del archivo GDS generado (en carpeta temporal)
gds_path <- file.path(tempdir(), "genotipos_combinados.gds")

# Abrir archivo GDS
genofile <- snpgdsOpen(gds_path)

# Ejecutar PCA
pca <- snpgdsPCA(genofile, autosome.only = FALSE)

# Leer archivo .fam para obtener IDs y poblaciones
fam <- read.table("genotipos_combinados.fam", header = FALSE)
colnames(fam) <- c("Poblacion", "Individuo", "Padre", "Madre", "Sexo", "Fenotipo")

# Confirmar correspondencia de IDs
ids_gds <- pca$sample.id
fam_ordenado <- fam[match(ids_gds, fam$Individuo), ]

# Verificación de seguridad
stopifnot(identical(as.character(fam_ordenado$Individuo), as.character(ids_gds)))

# Crear data frame con los resultados de PCA
pca_df <- data.frame(
  Individuo = fam_ordenado$Individuo,
  Poblacion = as.factor(fam_ordenado$Poblacion),
  PC1 = pca$eigenvect[, 1],
  PC2 = pca$eigenvect[, 2]
)

# Porcentaje de varianza explicada por cada componente
var_exp <- round(pca$varprop[1:2] * 100, 2)



jpeg(filename="PCA_poblaciones.jpg",
     width=11,
     height=9,units = "in",res = 600)

# Gráfico de dispersión
library(RColorBrewer)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Poblacion)) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)"),
    title = "Análisis de Componentes Principales (PCA)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
dev.off()



#Estimar Fst

# 🧼 Paso 0: Limpiar el entorno y cargar la librería
rm(list = ls())
gc()
library(hierfstat)

# 📥 Paso 1: Leer el archivo .dat exportado para FSTAT
archivo <- "genotipos_fstat_convertido.dat"
lineas <- readLines(archivo)

# 📊 Paso 2: Extraer metadatos correctamente (4 valores en la línea 1)
cabecera <- strsplit(lineas[1], "\\s+")[[1]]
num_pops <- as.numeric(cabecera[1])
num_ind  <- as.numeric(cabecera[2])
num_loci <- as.numeric(cabecera[3])

# 🧬 Paso 3: Leer nombres de loci y genotipos
nombres_loci <- lineas[2:(1 + num_loci)]
datos_gen <- lineas[(2 + num_loci):length(lineas)]

# 🔹 Separar por espacios
split_gen <- strsplit(datos_gen, "\\s+")
matriz_gen <- do.call(rbind, split_gen)

# ✔️ Validar número de columnas esperadas
if (ncol(matriz_gen) != (1 + num_loci)) {
  stop(paste0("❌ Esperado ", 1 + num_loci, " columnas, pero se encontraron ", ncol(matriz_gen)))
}

# ✅ Asignar nombres de columnas
colnames(matriz_gen) <- c("pop", nombres_loci)

# 🔄 Transformar a data.frame y tipos adecuados
datos_fstat <- as.data.frame(matriz_gen)
datos_fstat[] <- lapply(datos_fstat, as.character)
datos_fstat[, -1] <- lapply(datos_fstat[, -1], as.numeric)
datos_fstat$pop <- as.factor(datos_fstat$pop)

# Pronto para usar com hierfstat


resumen_stats <- basic.stats(datos_fstat)
fst_global <- resumen_stats$overall["Fst"]
cat("✅ Valor Fst global:", round(fst_global, 4), "\n")

# 🔸 Paso 6: Calcular Fst por pares de poblaciones
fst_pares <- pairwise.WCfst(datos_fstat,diploid=TRUE)
cat("📊 Fst por pares:\n")
print(round(fst_pares, 4))



###A análise DAPC (Discriminant Analysis of Principal Components, ou Análisis Discriminante de Componentes Principales)

rm(list = ls())

# ===============================
# ETAPA 1: Cargar paquetes necesários
# ===============================
#install.packages("adegenet")
#install.packages("tidyverse")
library(adegenet)
library(tidyverse)

# ===============================
# ETAPA 2: Leer archivo .raw generado con PLINK
system("plink --bfile genotipos_combinados --recodeA --out genotipos_raw")
# ===============================
# ============================================
# ETAPA 1: Leer archivo .raw do PLINK
# ============================================
genlight_obj <- read.PLINK("genotipos_raw.raw")

# ============================================
# ETAPA 2: Atribuir las 9 famílias como grupos
# ============================================

fam <- read.table("genotipos_combinados.fam", header = FALSE)

#"Los nombres de los individuos en el objeto genlight provienen de la columna V2 (por ejemplo: 'CAU_001')"
#"Las familias están en la columna V1 (valores del 1 al 9)"

pop(genlight_obj) <- fam$V1[match(indNames(genlight_obj), fam$V2)]

# Verificar

cat("Distribuição de indivíduos por família:\n")
print(table(pop(genlight_obj)))

# ============================================
# ETAPA 3: Análisis de PCA para estimar los PC ideales
# ============================================

pca_res <- glPca(genlight_obj, nf = 100)

# Variância explicada

explained <- pca_res$eig / sum(pca_res$eig) * 100
cum_var <- cumsum(explained)

# Graficar la varianza acumulada (entre 80% y 90%)

plot(cum_var, type = "b", pch = 19,
     xlab = "Componentes principais (PC)", ylab = "% Variância acumulada",
     main = "PCA – Variância explicada acumulada")

# Mostrar los primeros 70 componentes principales (PCs)

cat("Variância acumulada por PC (primeiros 70):\n")
print(cum_var[1:70])

# Definir el número de componentes principales (PCs) en base al gráfico

n_pcs <- 70  # Ajuste se necessário

# ============================================
# ETAPA 4: Executar DAPC com las 9 famílias
# ============================================

n_da <- length(unique(pop(genlight_obj))) - 1  # máx. = 8

dapc_familias <- dapc(genlight_obj, pop(genlight_obj), n.pca = n_pcs, n.da = n_da)

# ============================================
# ETAPA 5: Visualizar resultado
# ============================================

png(filename="DAPC_graph.png",
     width=800, height=600 )


scatter(dapc_familias, scree.da = TRUE, posi.da = "bottomright",
        main = "DAPC – Famílias 1 a 9")


dev.off()




#Estimar la matriz de parentesco genético (GRM) entre los 600 individuos de las 9 poblaciones para el control de cruzamientos.

# 📦 Instalar paquetes necesarios
#install.packages("AGHmatrix")
#install.packages("pheatmap")

# 📚 Cargar librerías
library(AGHmatrix)
library(pheatmap)


# Leer el archivo .raw exportado desde PLINK
datos <- read.table("genotipos_raw.raw", header = TRUE)

# Extraer solo la matriz de genotipos (eliminando FID, IID, PAT, MAT, SEX, PHENOTYPE)
genotipos <- datos[, -(1:6)]

# Convertir a matriz numérica (requerido por AGHmatrix)
geno_matrix <- as.matrix(genotipos)
mode(geno_matrix) <- "numeric"  # asegura que todo sea numérico

# Calcular la matriz de parentesco aditivo (VanRaden 2008)
G <- Gmatrix(SNPmatrix = geno_matrix,
             missingValue = NA,       # Cambiar a -9 si tus datos faltantes están codificados así
             maf = 0.05,
             method = "VanRaden")

# Asignar nombres de individuos
rownames(G) <- datos$IID
colnames(G) <- datos$IID

# Visualizar una parte de la matriz
round(G[1:5, 1:5], 3)

# Guardar la matriz en archivo CSV
write.csv(G, "matriz_parentesco_genomico.csv")

# 🔹 Crear un mapa de calor (heatmap)
pheatmap(G,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Mapa de Calor - Parentesco Genómico",
         fontsize = 6)


#Parentesco esperado	Valor aproximado (VanRaden)
#Idênticos (clones)	~1.00
#Irmãos completos	~0.50
#Meio-irmãos	~0.25
#Não aparentados	~0.00
#Relações negativas	< 0


##Clasificación y visualización de parentesco

# Transformar la matriz G en formato de pares
G_long <- as.data.frame(as.table(G))
colnames(G_long) <- c("Ind1", "Ind2", "Parentesco")

# Eliminar diagonales (comparación consigo mismo) y duplicados
G_long <- G_long[G_long$Ind1 != G_long$Ind2, ]
G_long <- G_long[!duplicated(t(apply(G_long[, 1:2], 1, sort))), ]

# Clasificar los pares según el grado de parentesco
G_long$Categoria <- cut(G_long$Parentesco,
                        breaks = c(-Inf, 0.1, 0.35, 0.65, Inf),
                        labels = c("No emparentados", "Medios hermanos", "Hermanos completos", "Posibles clones"))

# Contar el número de pares en cada categoría
tabla_categorias <- table(G_long$Categoria)
print(tabla_categorias)


##Gráfico de barras – Conteo por categoría

jpeg(filename="bar_GRM.jpg",
     width=11,
     height=9,units = "in",res = 600)

# Definir uma margem extra para o topo do gráfico
margen_superior <- max(tabla_categorias) * 0.15

# Aumentar a margem esquerda e afastar o título do eixo Y
par(mar = c(5, 7, 4, 2))     # aumenta a margem esquerda (segundo valor)
par(mgp = c(4, 1, 0))        # afasta o texto do eixo (primeiro valor)

# Definir margem superior para o topo do gráfico
margen_superior <- max(tabla_categorias) * 0.15

# Criar o gráfico e capturar as posições das barras

bp <- barplot(tabla_categorias,
              col = "red",
              main = "Número de pares por categoría de parentesco",
              ylab = "Número de pares",
              xlab = "Categoría de parentesco",
              las = 1,
              ylim = c(0, max(tabla_categorias) + margen_superior))

# Adicionar os valores sobre as barras
text(x = bp,
     y = tabla_categorias + max(tabla_categorias) * 0.02,
     labels = tabla_categorias,
     cex = 0.9,
     font = 2)

dev.off()


jpeg(filename="histo_GRM.jpg",
     width=11,
     height=9,units = "in",res = 600)

hist(G_long$Parentesco,
     breaks = 50,
     col = "lightgreen",
     main = "Distribución de valores de parentesco genómico",
     xlab = "Parentesco",
     ylab = "Frecuencia")

dev.off()



# graficar por población


# 📦 Instalar paquetes necesarios
#install.packages("ggplot2")
#install.packages("pheatmap")
#install.packages("AGHmatrix")

# 📚 Cargar librerías
library(ggplot2)
library(pheatmap)
library(AGHmatrix)

# ──────────────────────────────────────────────────────────────
# 1. Leer el archivo .fam y crear info_populaciones.csv
# ──────────────────────────────────────────────────────────────

# Leer archivo .fam (columna 1: población; columna 2: ID del individuo)
fam <- read.table("genotipos_combinados.fam", header = FALSE)

# Crear archivo info_populaciones.csv
info_pop <- data.frame(ID = fam$V2, Poblacion = fam$V1)
write.csv(info_pop, "info_populaciones.csv", row.names = FALSE)
cat("✅ Archivo 'info_populaciones.csv' creado automáticamente a partir de genotipos_combinados.fam\n")

# ──────────────────────────────────────────────────────────────
# 2. Calcular la matriz de parentesco genómico (VanRaden)
# ──────────────────────────────────────────────────────────────

# Leer archivo .raw exportado desde PLINK
datos <- read.table("genotipos_raw.raw", header = TRUE)

# Extraer genotipos y convertir a matriz numérica
genotipos <- datos[, -(1:6)]
geno_matrix <- as.matrix(genotipos)
mode(geno_matrix) <- "numeric"

# Calcular matriz G con método de VanRaden
G <- Gmatrix(SNPmatrix = geno_matrix,
             missingValue = NA,
             maf = 0.05,
             method = "VanRaden")

# Asignar nombres de individuos
rownames(G) <- datos$IID
colnames(G) <- datos$IID

# ──────────────────────────────────────────────────────────────
# 3. Construir G_long (pares) y asociar poblaciones
# ──────────────────────────────────────────────────────────────

# Convertir matriz G a tabla larga (pares)
G_long <- as.data.frame(as.table(G))
colnames(G_long) <- c("Ind1", "Ind2", "Parentesco")

# Eliminar diagonales y duplicados
G_long <- G_long[G_long$Ind1 != G_long$Ind2, ]
G_long <- G_long[!duplicated(t(apply(G_long[, 1:2], 1, sort))), ]

# Leer archivo info_populaciones
info_pop <- read.csv("info_populaciones.csv", stringsAsFactors = FALSE)

# Asociar población al individuo 1
G_long <- merge(G_long, info_pop, by.x = "Ind1", by.y = "ID", all.x = TRUE)
colnames(G_long)[ncol(G_long)] <- "Pob1"

# Asociar población al individuo 2
G_long <- merge(G_long, info_pop, by.x = "Ind2", by.y = "ID", all.x = TRUE)
colnames(G_long)[ncol(G_long)] <- "Pob2"

# ──────────────────────────────────────────────────────────────
# 4. Clasificar por parentesco y filtrar intra-población
# ──────────────────────────────────────────────────────────────

# Filtrar pares dentro de la misma población
G_intra <- subset(G_long, Pob1 == Pob2)
G_intra$Poblacion <- G_intra$Pob1

# Clasificar cada par según el valor de parentesco
G_intra$Categoria <- cut(G_intra$Parentesco,
                         breaks = c(-Inf, 0.1, 0.35, 0.65, Inf),
                         labels = c("No emparentados", "Medios hermanos", "Hermanos completos", "Posibles clones"))

# ──────────────────────────────────────────────────────────────
# 5. Contar pares por categoría y población
# ──────────────────────────────────────────────────────────────

# Tabla cruzada: población × categoría de parentesco
tabla_por_pob <- table(G_intra$Poblacion, G_intra$Categoria)
df_plot <- as.data.frame(tabla_por_pob)
colnames(df_plot) <- c("Poblacion", "Categoria", "Frecuencia")

# ──────────────────────────────────────────────────────────────
# 6. Gráfico de barras agrupadas por población
# ──────────────────────────────────────────────────────────────

jpeg(filename="barras_poblacion.jpg",
     width=11,
     height=9,units = "in",res = 600)

ggplot(df_plot, aes(x = Poblacion, y = Frecuencia, fill = Categoria)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Pares por categoría de parentesco en cada población",
       x = "Población",
       y = "Número de pares") +
  scale_fill_manual(values = c("gray80", "skyblue", "orange", "tomato")) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

























