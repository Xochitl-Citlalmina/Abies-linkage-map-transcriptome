#####Poisson

# Cargar bibliotecas necesarias
library(dplyr)
library(ggplot2)
library(MASS)  # Para la función de Binomial Negativa

# Leer los datos
datos <- read.csv("allLGFig.txt", sep = "\t")

# Crear una columna para los intervalos de 10 cM
datos <- datos %>%
  mutate(intervalo = floor(POS_avg / 20) * 20)  # Agrupa por intervalos de 10 cM

# Agrupar por CHR y intervalo, y contar la cantidad de marcadores
frecuencias <- datos %>%
  group_by(CHR, intervalo) %>%
  summarise(frecuencia = n()) %>%
  ungroup()

# Función para ajustar la distribución de Poisson y Binomial Negativa
ajustar_distribuciones <- function(frecuencias) {
  # Poisson: estimar lambda como la media de las frecuencias observadas
  lambda_poisson <- mean(frecuencias$frecuencia)
  frecuencias$esperada_poisson <- dpois(frecuencias$frecuencia, lambda = lambda_poisson) * sum(frecuencias$frecuencia)
  
  # Binomial Negativa: estimar parámetros size y mu
  fit_nb <- fitdistr(frecuencias$frecuencia, "negative binomial")
  size <- fit_nb$estimate["size"]
  mu <- fit_nb$estimate["mu"]
  frecuencias$esperada_binomial_neg <- dnbinom(frecuencias$frecuencia, size = size, mu = mu) * sum(frecuencias$frecuencia)
  
  return(frecuencias)
}

# Aplicar el ajuste de distribuciones para cada CHR
frecuencias_ajustadas <- frecuencias %>%
  group_by(CHR) %>%
  do(ajustar_distribuciones(.))

# Graficar las frecuencias observadas y las curvas ajustadas de Poisson y Binomial Negativa
ggplot(frecuencias_ajustadas, aes(x = intervalo, y = frecuencia)) +
  geom_bar(stat = "identity", aes(fill = as.factor(CHR)), position = "dodge", alpha = 0.6) +  # Frecuencia observada
  geom_line(aes(y = esperada_poisson, color = "Poisson"), size = 1) +  # Curva Poisson
  geom_line(aes(y = esperada_binomial_neg, color = "Negative Binomial"), size = 1) +  # Curva Binomial Negativa
  facet_wrap(~ CHR, scales = "free_y") +  # Un gráfico por cada grupo de ligamiento
  labs(x = "Interval of 20 cM", y = "Frecuency", fill = "Linkage Group", color = "Distribution") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, NA))  # Asegura que los ejes Y sean consistentes

# Especificar la ruta y nombre del archivo
output_file <- "bin20cM.jpg"

# Guardar la gráfica en formato JPEG con alta calidad
ggsave(
  filename = output_file,          # Nombre del archivo
  plot = last_plot(),              # La última gráfica que se generó
  device = "jpeg",                 # Formato del archivo (JPEG)
  path = "C:/Users/Usuario/Mi unidad/Posdoc/ResultadosLepMap/Marzo_24/Figuras_finales/",  # Ruta a la carpeta donde se guardará
  dpi = 300,                       # Resolución en puntos por pulgada (DPI), ajustable para alta calidad
  width = 12, height = 8,          # Ajuste de las dimensiones en pulgadas
  units = "in"                     # Unidades de las dimensiones ("in" para pulgadas)
)

