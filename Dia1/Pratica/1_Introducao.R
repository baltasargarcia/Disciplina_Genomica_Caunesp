#############################
# Introdução ao R
#############################

##### Diretório de trabalho #####

# Qual é o diretório atual?
getwd() 

# Criando um novo diretório
dir.create("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia1/Pratica/Test") 

# Definindo o diretório de trabalho
setwd("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia1/Pratica/Test")                      
getwd() 

##### Instalação e carregamento de pacotes #####

# Instalar o pacote (descomente se necessário)
# install.packages("ggplot2")

# Carregar o pacote
library(ggplot2) 

##### Tipos de objetos em R #####

# atribuir um valor númerico a um objeto
a <- 3
a
b = 5
b

# criação de vetores (concatenar elementos em um vetor)
x <- c(1, 2, 3, 4, 5, 6) 
x
x2 <- 1:6
x2
# No R, colocamos texto (strings) entre aspas (' ' ou " ") para indicar que são valores de texto e não nomes de variáveis ou números.
w <- c('an1', 'an2', 'an3')
w

# Verificando o tipo dos objetos
class(x)
class(w)

# Simulando valores com distribuição normal
#cria um vetor y com 6 números gerados aleatoriamente de uma distribuição normal com média 10 e desvio padrão 5.
y <- rnorm(6, mean = 10, sd = 5)
y

# Ajuda sobre funções
help(rnorm)
?rnorm

##### Matrizes #####

# Criando matrizes com colunas e linhas
#cbind = combina vetores como colunas
#rbind = combina vetores como linhas
X <- cbind(x, y) # Atenção às letras maiúsculas e minúsculas
X2 <- rbind(x, y)
X
X2

# Outra maneira de criar uma matriz
X2 <- matrix(nrow = 6, ncol = 2, data = c(x, y))
X2

##### Fatores e caracteres #####

x3 <- c("low", "medium", "low", "high")
x3

x4 <- as.factor(x3)
x4
# as.factor = texto em variável categórica
x5 <- factor(c("AA", "AT", "TT", "TT", "AT", "AT", "AA", "TT"))
x5

class(x3)
class(x4)

##### Operadores matemáticas básicas #####

# vetores
x <- 1:10
y <- 11:20
# O ; permite executar múltiplos comandos em uma linha
x; y

#operações com vetores
z <- x + y
z
x * x
x^2
x / 2

# Condições lógicas em vetores
z2 <- x > 5
z2
z3 <- x == 3
z3
z4 <- x != 3
z4

##### Estrutura dos objetos #####

str(x)
str(z)
str(a)
str(x3)
str(x4)

##### Listar e remover objetos #####

# Listar objetos no ambiente
ls()

# Remover objetos específicos
rm(x2, z2)
ls()

# Remover todos os objetos
rm(list = ls())
ls()

##### Funções úteis #####

# Definir semente para reprodutibilidade
set.seed(12345)

# Simular 1000 valores com média 100 e desvio padrão 20
x <- rnorm(1000, mean = 100, sd = 20)

# Estatísticas básicas

#média
mean(x)

#desvio padrão
sd(x)

# variância amostral de x, mede o quanto os valores se dispersam em torno da média, usando a fórmula: soma dos quadrados dos desvios dividida por (n - 1).
var(x)

#retorna um vetor com o valor mínimo e o valor máximo de x.
range(x)

min(x)
max(x)

#calcula os quantis de x para as probabilidades especificadas: 10%, 30%, 50% e 90%.
quantile(x, prob = c(0.1, 0.3, 0.5, 0.9))

summary(x)

# Sumário de fator
x5 <- factor(c("AA", "AT", "TT", "TT", "AT", "AT", "AA", "TT"))
summary(x5)

# Plot da densidade
plot(density(x))

##### Simulação de dados #####

#tamanho da amostra
N <- 30

# gera 30 números aleatórios de uma distribuição normal padrão, (média 0, desvio padrão 1)
x1 <- rnorm(N)

#30 valores de uma distribuição binomial com size = 1 (sucesso ou fracasso, ex: 1 = infectado; 0 = resistente), com probabilidade de sucesso 0.2
x2 <- rbinom(n = N, size = 1, prob = 0.2)

# Tabela de frequências
table(x2)

# Criar matriz com os dados simulados
X <- cbind(x1, x2)
head(X)
tail(X)

##### Data frames #####

# Criar data frame

myData <- data.frame(Peso = x1, sex = x2)

#dimensão do meu data frame
dim(myData)

head(myData)

#estrutura do meu data frame
str(myData)

# Acessar colunas
myData$Peso
myData$sex

##### Escrita e leitura de arquivos #####

getwd()

# Escrever arquivo de texto
write.table(myData, file = "myData.txt", row.names = FALSE, col.names = TRUE, sep = "\t")

# Limpar ambiente
rm(list = ls())
ls()
getwd()

# Ler arquivo
myData <- read.table("myData.txt", sep = "\t", header = TRUE)
head(myData)
str(myData)
