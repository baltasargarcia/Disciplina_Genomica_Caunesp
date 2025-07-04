# Limpa o ambiente
rm(list = ls())

# Carregando pacotes necessários
library(ggplot2)
library(dplyr)
library(magrittr)
library(patchwork)
library(plotly)

# ============================================================
# 1. Leitura dos dados
# ============================================================

# Lê as tabelas dos três métodos
pblup   <- read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/PBLUP/Acc_PBLUP.txt", header = TRUE)
gblup   <- read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/GBLUP/Acc_GBLUP.txt", header = TRUE)
ssgblup <- read.table("C:/Post_Doc_UNESP/Cursos_congressos_discip/Disciplina_Caunesp_2025/Dia4/Pratica/ssGBLUP/Acc_ssGBLUP.txt", header = TRUE)

# Mantém apenas os cinco grupos de validação (remove média e desvio)
pblup   <- pblup[1:5, ]
gblup   <- gblup[1:5, ]
ssgblup <- ssgblup[1:5, ]

# Junta os três métodos em um único dataframe
final <- rbind(pblup, gblup, ssgblup)
final$Metodo <- factor(c(rep("PBLUP", 5), rep("GBLUP", 5), rep("ssGBLUP", 5)),
                       levels = c("PBLUP", "GBLUP", "ssGBLUP"))

# Garante a ordem dos grupos de validação
final$CV_group <- factor(final$CV_group, levels = c("CV1", "CV2", "CV3", "CV4", "CV5"))

# ============================================================
# 2. Gráficos de barras para Acurácia1 e Acurácia2 (por grupo)
# ============================================================

grafico_acc1 <- ggplot(final, aes(x = CV_group, y = Acurácia1, fill = Metodo)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(x = "Grupo de Validação", y = "Acurácia 1", fill = "Método") +
  ggtitle("Acurácia 1 por grupo de validação") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1))

grafico_acc2 <- ggplot(final, aes(x = CV_group, y = Acurácia2, fill = Metodo)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(x = "Grupo de Validação", y = "Acurácia 2", fill = "Método") +
  ggtitle("Acurácia 2 por grupo de validação") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1))

# Exibe os dois gráficos lado a lado
grafico_acc1 + grafico_acc2 + plot_layout(ncol = 2)

# ============================================================
# 3. Resumo: médias, desvios e aumentos percentuais
# ============================================================

# Função para calcular médias, desvios e aumento percentual em relação ao PBLUP
calc_resumo <- function(coluna_acuracia) {
  resumo <- final %>%
    group_by(Metodo) %>%
    summarise(media = mean(.data[[coluna_acuracia]]),
              desvio = sd(.data[[coluna_acuracia]])) %>%
    mutate(aumento = ifelse(Metodo == "PBLUP", NA, (media - media[Metodo == "PBLUP"]) / media[Metodo == "PBLUP"] * 100))
  return(resumo)
}

resumo_ac1 <- calc_resumo("Acurácia1")
resumo_ac2 <- calc_resumo("Acurácia2")

# ============================================================
# 4. Gráficos de barras para médias + % aumento
# ============================================================

grafico1 <- ggplot(resumo_ac1, aes(x = Metodo, y = media, fill = Metodo)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_errorbar(aes(ymin = media - desvio, ymax = media + desvio), width = 0.2) +
  geom_text(aes(label = ifelse(!is.na(aumento), paste0(sprintf("%.1f", aumento), "%"), "")),
            vjust = -3, size = 7, fontface = "bold") +
  theme_minimal() +
  labs(x = "Método", y = "Acurácia 1 média", fill = "Método") +
  ggtitle("Acurácia 1 média por método") +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1))

grafico2 <- ggplot(resumo_ac2, aes(x = Metodo, y = media, fill = Metodo)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_errorbar(aes(ymin = media - desvio, ymax = media + desvio), width = 0.2) +
  geom_text(aes(label = ifelse(!is.na(aumento), paste0(sprintf("%.1f", aumento), "%"), "")),
            vjust = -3, size = 7, fontface = "bold") +
  theme_minimal() +
  labs(x = "Método", y = "Acurácia 2 média", fill = "Método") +
  ggtitle("Acurácia 2 média por método") +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16, face = "bold"),
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1))

# Exibe os dois gráficos lado a lado
grafico1 + grafico2 + plot_layout(ncol = 2)

# ============================================================
# 5. Gráfico de linhas do viés por grupo de validação
# ============================================================

# Calcula a distância do viés ao valor ideal (1.0)
final$distancia_do_ideal <- round(final$Viés - 1, 3)

grafico_vies <- ggplot(final, aes(x = CV_group, y = Viés, color = Metodo, group = Metodo)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_text(aes(label = distancia_do_ideal), vjust = -1, size = 4, fontface = "bold") +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = 5.3, y = 1.03, label = "Viés ideal", color = "red", fontface = "bold", size = 5) +
  scale_color_manual(values = c("PBLUP" = "#E69F00", "GBLUP" = "#56B4E9", "ssGBLUP" = "#009E73")) +
  theme_minimal() +
  labs(title = "Viés estimado por grupo de validação e método",
       x = "Grupo de Validação (CV)",
       y = "Viés",
       color = "Método") +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
  coord_cartesian(ylim = c(0.5, 1.3))

print(grafico_vies)

# ============================================================
# 6. Gráfico interativo com plotly
# ============================================================

grafico_interativo <- plot_ly(final,
                              x = ~CV_group,
                              y = ~Viés,
                              color = ~Metodo,
                              type = 'scatter',
                              mode = 'lines+markers+text',
                              text = ~paste("Viés:", round(Viés, 3),
                                            "<br>Método:", Metodo,
                                            "<br>Grupo:", CV_group,
                                            "<br>Dist. do ideal:", round(Viés - 1, 3)),
                              hoverinfo = 'text') %>%
  layout(title = "Viés por Grupo de Validação",
         yaxis = list(title = "Viés", range = c(0.5, 1.3)),
         xaxis = list(title = "Grupo de Validação"),
         shapes = list(list(type = "line",
                            x0 = 0, x1 = 1,
                            xref = "paper",
                            y0 = 1, y1 = 1,
                            line = list(color = "red", dash = "dash"))))

grafico_interativo
