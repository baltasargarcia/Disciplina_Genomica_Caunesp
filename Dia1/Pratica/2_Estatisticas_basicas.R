rm(list=ls())
require(ggplot2)
require(ggpubr)
#install.packages("GGally")
require(GGally)


#Carpeta de trabajo
setwd("C:/Users/Usuario/Desktop")

#Transformar los datos de .xls para .txt

#install.packages("readxl")
library(readxl)

dados <- read_excel("Desafio_IPN.xlsx")

write.table(dados, file = "Desafio_IPN.txt", sep = "\t", row.names = FALSE)

#Leyendo el archivo de datos
data<-read.table("Desafio_IPN.txt",h=T)
head(data)

#estadisticas descriptivas basicas
mean(data$Mortal_weight)
sd(data$Mortal_weight)
summary(data$Mortal_weight)

mean(data$Mortal_DDM)
sd(data$Mortal_DDM)
summary(data$Mortal_DDM)


#plots basicos
plot(data$Mortal_weight)
plot(data$Mortal_DDM)

hist(data$Mortal_weight)
hist(data$Mortal_DDM)

#plot normalidad 
#compara la distribución de tus datos con una distribución normal — si los puntos caen cerca de una línea recta, indica que los datos son aproximadamente normales.
#El eje x muestra los cuantiles teóricos de la distribución normal (valores esperados bajo normalidad).
#El eje y muestra los cuantiles muestrales de tus datos.

qqnorm(data$Mortal_weight)
qqline(data$Mortal_weight)
qqnorm(data$Mortal_DDM)
qqline(data$Mortal_DDM)

#boxplots
boxplot(data$Mortal_weight)
boxplot(data$Mortal_DDM)

#Correlacion entre dos variables
cor(data$Tag_weight,data$Tag_length)
cor.test(data$Tag_weight,data$Tag_length)

#Plot de correlaci?n
ggplot(data) +
  aes(x = Tag_weight, y = Tag_length) +
  geom_point(colour = "blue") +
  theme_minimal()+
geom_smooth(mapping = aes(x = Tag_weight, y = Tag_length))

#Edición de datos (removiendo datos)
data[data$Tag_weight > 4.9,]<-NA

#Plot de correlaci?n
ggplot(data) +
  aes(x = Tag_weight, y = Tag_length) +
  geom_point(colour = "blue") +
  theme_minimal()+
  geom_smooth(mapping = aes(x = Tag_weight, y = Tag_length))+
  stat_cor(method = "pearson", label.x = 3, label.y = 8)

#Plot de correlacion con v?rias variables
data2<-data[,c(10,11,14)]
ggpairs(data2)

#Diagonal: densidades de cada variable.
#Debajo de la diagonal: gráficos de dispersión (scatterplots) que muestran la relación entre pares de variables.
#Encima de la diagonal: coeficientes de correlación
  
