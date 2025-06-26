rm(list=ls())
require(ggplot2)
require(ggpubr)
#install.packages("GGally")
require(GGally)
#Transformar los datos de .xls para .txt

#Carpeta de trabajo
setwd("C:/Users/Usuario/Desktop")

#Leyendo el archivo de datos
data<-read.table("Desafio_IPN.txt",h=T)
head(data)

#estad?sticas descriptivas b?sicas
mean(data$Mortal_weight)
sd(data$Mortal_weight)
summary(data$Mortal_weight)

mean(data$Mortal_DDM)
sd(data$Mortal_DDM)
summary(data$Mortal_DDM)


#plots b?sicos
plot(data$Mortal_weight)
plot(data$Mortal_DDM)

hist(data$Mortal_weight)
hist(data$Mortal_DDM)

#plot normalidad
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

#Edici?n de datos (removiendo datos)
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
  
