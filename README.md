# 🧬 Tópicos Especiais: Seleção Genômica na Aquicultura

Este repositório reúne os **materiais didáticos da disciplina _"Tópicos Especiais: Genômica Aplicada à Conservação e Melhoramento Genético de Peixes"_**, voltada a estudantes e profissionais interessados em integrar ferramentas genômicas a programas de melhoramento genético de organismos aquáticos.

---

## 💻 Softwares Necessários

Para otimizar as atividades práticas do curso, solicitamos que todos os participantes instalem previamente os seguintes softwares:

### 1. R (base) e RStudio
- **R base**: [Download](https://brieger.esalq.usp.br/CRAN/)  
  (Escolha o sistema operacional correspondente)
- **RStudio**: [Download](https://posit.co/download/rstudio-desktop/)  
  (Interface mais amigável e funcional)

## Pacotes a serem utilizados dentro do ambiente R (solicitamos que a instalação seja realizada no RStudio previamente à Aula 2)
install.packages("ggplot2")
library(ggplot2)

install.packages('reshape2')
library(reshape2)

install.packages("devtools")
library(devtools)
install_github("jgx65/hierfstat")
library("hierfstat")

options(
  repos = c(
    zkamvar = "https://zkamvar.r-universe.dev",
    CRAN = "https://cloud.r-project.org"
  )
)
install.packages("adegenet")
library(adegenet)

install.packages("data.table")
library(data.table)

install.packages("stringr")
library(stringr)

if (!require("SNPRelate")) install.packages("BiocManager"); BiocManager::install("SNPRelate")
library(SNPRelate)

> ⚠️ O RStudio depende do R base. Instale o R base antes do RStudio.


---

### 2. PLINK
- Site oficial: [https://www.cog-genomics.org/plink/](https://www.cog-genomics.org/plink/)
- Baixe a versão estável mais recente ("Stable") de acordo com seu sistema operacional.
- Após o download (`.zip`), extraia o arquivo `plink.exe`.  
  **Não é necessário instalar** — basta mantê-lo em uma pasta acessível.

---

### 3. Pacotes BLUPF90
- Disponíveis em: [https://nce.ads.uga.edu/html/projects/programs/](https://nce.ads.uga.edu/html/projects/programs/)
- Baixe os arquivos diretamente na pasta **"Temp"**.
- Módulos necessários:
  - `renumf90`
  - `airemlf90`
  - `blupf90`
  - `preGSf90`

> 📁 Os arquivos são executáveis e não requerem instalação.

---

### 4. Biblioteca `libiomp5md.dll`
- Essencial para o funcionamento dos programas BLUPF90.
- Copie esta biblioteca para **todas as pastas** onde serão executados os pacotes BLUPF90.
- [Download aqui](https://pt.dll-files.com/download/3a7902626cddec83a3da541a96118b46/libiomp5md.dll.html?c=ZFZVcSs5RE9jOTAwMjdpeGlWS3dkUT09)

---

### 5. TextPad *(opcional – somente para Windows)*
- [Download aqui](https://www.textpad.com/download)
- Recomendado para abrir arquivos de texto grandes, que podem ser lentos com o Bloco de Notas padrão do Windows.

---

## 📂 Material das Aulas Práticas

Os dados utilizados nas atividades estão organizados em **pastas por dia de aula**.  
Ao fazer o download do material, **mantenha a estrutura original das pastas**, pois ela é essencial para a execução correta dos scripts e análises.

---

## Dúvidas

baltasar.garcia@unesp.br
mastrochirico.filho@unesp.br
