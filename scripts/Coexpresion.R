# Instalar paquetes necesarios si no los tienes
if (!requireNamespace("WGCNA", quietly = TRUE)) install.packages("WGCNA")

# Cargar librería
library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# ===== PASO 1: Cargar expresión cruda =====
expresion_raw <- read.csv("conteo_data_sin_duplicados.csv", header = TRUE, row.names = 1, check.names = FALSE)

# ===== PASO 2: Filtrar genes poco expresados =====
# Conservamos genes expresados en al menos la mitad de las muestras con al menos 10 cuentas
filtro <- rowSums(expresion_raw >= 10) >= (ncol(expresion_raw) / 2)
expresion_filtrada <- expresion_raw[filtro, ]
cat("Genes después del filtrado:", nrow(expresion_filtrada), "\n")

# ===== PASO 3: Transponer para WGCNA (muestras en filas, genes en columnas) =====
datExpr <- t(expresion_filtrada)

# ===== PASO 4: Verificar calidad de los datos =====
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  cat("Se eliminaron muestras o genes con muchos valores faltantes o constantes.\n")
}

# ===== PASO 5: Clustering de muestras (opcional, para detección de outliers) =====
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Clustering jerárquico de muestras", xlab="", sub="", cex.lab=1.2)

# Remover muestra outlier
datExpr <- datExpr[!rownames(datExpr) %in% "SS01_15", ]
sampleTree2 <- hclust(dist(datExpr), method = "average")
plot(sampleTree2, main = "Clustering sin muestra outlier (SS01_15)", xlab="", sub="")

# === PASO 4: Elegir soft-thresholding power ===
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Visualización
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit (R²)",
     type="n", main = "Seleccionar power")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="red")
abline(h=0.8, col="blue")  # target: R² > 0.8

# === PASO 5: Elegir un power apropiado (ajusta según tu gráfica) ===
softPower <- 16  # <- cambia si ves otro mejor en tu gráfica

# Guardamos temporalmente los nombres originales de muestra
muestra_nombres <- rownames(datExpr)

# Convertimos a numeric
datExpr <- as.data.frame(apply(datExpr, 2, as.numeric))

# Restauramos los nombres de las muestras correctamente
rownames(datExpr) <- muestra_nombres

# === PASO 6: Construcción de red y detección de módulos ===
net <- blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "unsigned",   
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  verbose = 3
)


# === PASO 7: Resultados: colores y asignación de módulos ===
moduleColors <- labels2colors(net$colors)
table(moduleColors)  # número de genes por módulo

# === PASO 8: Visualización de dendrograma con módulos ===
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Módulos detectados", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)

# === PASO 9: Guardar genes y sus módulos ===
genes_modulos <- data.frame(GENE = colnames(datExpr), modulo = moduleColors)
write.table(genes_modulos, file = "genes_por_modulo_WGCNA.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# === PASO 10 (opcional): Guardar la red completa ===
save(net, file = "red_WGCNA_modulos.RData")

cat("✅ Análisis WGCNA terminado. Módulos detectados y guardados.\n")

####Unir con mapa

# Cargar tabla de genes y su módulo (de WGCNA)
modulos <- read.table("genes_por_modulo_WGCNA.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Cargar mapa de ligamiento
mapa <- read.table("combined_map_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Unir por nombre del gen
modulos_mapa <- merge(modulos, mapa, by = "GENE")

# Ver primeras líneas
head(modulos_mapa)

# Guardar la tabla unificada
write.table(modulos_mapa, file = "modulos_con_posicion.txt", sep = "\t", quote = FALSE, row.names = FALSE)

library(ggplot2)

# Histograma: número de genes por módulo y grupo de ligamiento (CHR)
ggplot(modulos_mapa, aes(x = CHR, fill = modulo)) +
  geom_bar(position = "stack") +
  labs(title = "Distribución de genes por módulo y grupo de ligamiento",
       x = "Grupo de ligamiento (CHR)", y = "Número de genes") +
  theme_minimal()

# Calcular distancia promedio entre genes del mismo módulo en cada CHR
distancias_modulo <- do.call(rbind, lapply(split(modulos_mapa, modulos_mapa$modulo), function(df_mod) {
  do.call(rbind, lapply(split(df_mod, df_mod$CHR), function(df_chr) {
    if (nrow(df_chr) >= 2) {
      dists <- dist(sort(df_chr$POS_avg))
      data.frame(
        MODULO = unique(df_chr$modulo),
        CHR = unique(df_chr$CHR),
        N_genes = nrow(df_chr),
        dist_promedio = mean(dists)
      )
    }
  }))
}))

# Ver distancias promedio entre genes coexpresados
print(distancias_modulo)

# Guardar resultados
write.table(distancias_modulo, file = "distancia_promedio_por_modulo_y_CHR.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#############CORRELACIONES

# Requiere:
# - datExpr: expresión con muestras en filas y genes en columnas
# - modulos_mapa: tabla con GENE, CHR, POS_avg

library(dplyr)

# Asegurar solo genes comunes
genes_comunes <- intersect(modulos_mapa$GENE, colnames(datExpr))
mod_df <- modulos_mapa %>% filter(GENE %in% genes_comunes)
expr_df <- datExpr[, mod_df$GENE]
mod_df <- mod_df[match(colnames(expr_df), mod_df$GENE), ]

# Función para obtener pares y correlaciones dentro de cada CHR
correlacion_por_chr <- function(chr_name, data_expr, data_map) {
  # Subconjunto por CHR
  sub_map <- data_map %>% filter(CHR == chr_name)
  if (nrow(sub_map) < 2) return(NULL)
  
  sub_expr <- data_expr[, sub_map$GENE]
  if (ncol(sub_expr) < 2) return(NULL)
  
  # Correlación de expresión
  cor_mat <- cor(sub_expr, method = "pearson")
  
  # Matriz de distancias físicas
  dist_mat <- as.matrix(dist(sub_map$POS_avg))
  colnames(dist_mat) <- rownames(dist_mat) <- sub_map$GENE
  
  # Obtener pares únicos
  results <- data.frame()
  for (i in 1:(nrow(dist_mat)-1)) {
    for (j in (i+1):ncol(dist_mat)) {
      g1 <- rownames(dist_mat)[i]
      g2 <- colnames(dist_mat)[j]
      dist <- dist_mat[i, j]
      cor_val <- cor_mat[g1, g2]
      
      # Clasificación por distancia
      bin <- ifelse(dist <= 1, "≤1cM",
                    ifelse(dist <= 5, "≤5cM", NA))
      if (!is.na(bin)) {
        results <- rbind(results, data.frame(
          CHR = chr_name,
          gene1 = g1,
          gene2 = g2,
          distancia_cM = dist,
          correlacion = cor_val,
          bin = bin
        ))
      }
    }
  }
  return(results)
}

# Aplicar a cada CHR
CHRs <- unique(mod_df$CHR)
res_chr <- do.call(rbind, lapply(CHRs, correlacion_por_chr, data_expr = expr_df, data_map = mod_df))

# Promedio por CHR y bin
resumen <- res_chr %>%
  group_by(CHR, bin) %>%
  summarise(
    n_pares = n(),
    promedio_correlacion = mean(correlacion, na.rm = TRUE),
    .groups = "drop"
  )

# Mostrar resultados
print(resumen)

# Guardar
write.table(resumen, file = "correlacion_por_CHR_1cM_5cM.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(res_chr, file = "correlaciones_pares_por_CHR.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Si no está cargado, volver a cargar resumen:
resumen <- read.table("correlacion_por_CHR_1cM_5cM.txt", header = TRUE, sep = "\t")

library(ggplot2)

# Barplot
ggplot(resumen, aes(x = CHR, y = promedio_correlacion, fill = bin)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Correlación promedio entre genes por CHR",
       subtitle = "Para pares de genes a ≤1 cM y ≤5 cM",
       x = "Grupo de ligamiento (CHR)",
       y = "Correlación de expresión (Pearson)",
       fill = "Distancia genética") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Asegurarte de tener estas tablas listas:
# - expr_df: expresión de los genes (solo los presentes en el mapa), columnas = genes
# - mod_df: data frame con columnas GENE, modulo, CHR, POS_avg

# Convertir modulo a factor para facilitar gráficos
mod_df$modulo <- as.factor(mod_df$modulo)

# Función para un módulo en un CHR
correlacion_mod_chr <- function(chr, modulo, data_expr, data_map) {
  sub_map <- data_map %>% filter(CHR == chr, modulo == modulo)
  if (nrow(sub_map) < 2) return(NULL)
  
  sub_expr <- data_expr[, sub_map$GENE]
  if (ncol(sub_expr) < 2) return(NULL)
  
  cor_mat <- cor(sub_expr, method = "pearson")
  dist_mat <- as.matrix(dist(sub_map$POS_avg))
  colnames(dist_mat) <- rownames(dist_mat) <- sub_map$GENE
  
  results <- data.frame()
  for (i in 1:(nrow(dist_mat)-1)) {
    for (j in (i+1):ncol(dist_mat)) {
      g1 <- rownames(dist_mat)[i]
      g2 <- colnames(dist_mat)[j]
      dist <- dist_mat[i, j]
      cor_val <- cor_mat[g1, g2]
      bin <- ifelse(dist <= 1, "≤1cM",
                    ifelse(dist <= 5, "≤5cM", NA))
      if (!is.na(bin)) {
        results <- rbind(results, data.frame(
          CHR = chr,
          modulo = modulo,
          gene1 = g1,
          gene2 = g2,
          distancia_cM = dist,
          correlacion = cor_val,
          bin = bin
        ))
      }
    }
  }
  return(results)
}

# Aplicarlo para todos los módulos y CHRs
mod_chr_list <- expand.grid(
  CHR = unique(mod_df$CHR),
  modulo = unique(mod_df$modulo),
  stringsAsFactors = FALSE
)

res_mod_chr <- do.call(rbind, mapply(
  correlacion_mod_chr,
  chr = mod_chr_list$CHR,
  modulo = mod_chr_list$modulo,
  MoreArgs = list(data_expr = expr_df, data_map = mod_df),
  SIMPLIFY = FALSE
))

# Resumen
resumen_mod_chr <- res_mod_chr %>%
  group_by(modulo, CHR, bin) %>%
  summarise(
    n_pares = n(),
    promedio_correlacion = mean(correlacion, na.rm = TRUE),
    .groups = "drop"
  )

# Guardar
write.table(resumen_mod_chr, file = "correlacion_por_modulo_CHR.txt", sep = "\t", quote = FALSE, row.names = FALSE)

library(ggplot2)

# Convertir a factor para orden correcto
resumen_mod_chr$modulo <- factor(resumen_mod_chr$modulo)
resumen_mod_chr$CHR <- factor(resumen_mod_chr$CHR)

ggplot(resumen_mod_chr, aes(x = CHR, y = promedio_correlacion, fill = bin)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ modulo, scales = "free_y") +
  labs(title = "Correlación de expresión por módulo y grupo de ligamiento",
       subtitle = "Promedio entre pares de genes a ≤1 cM y ≤5 cM",
       x = "CHR", y = "Correlación (Pearson)", fill = "Distancia genética") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

####################ANALISIS DE MODULOS POR CHR##########

library(dplyr)
library(ggplot2)

# Cargar archivo de texto separado por tabulaciones
datos <- read.table("modulos_con_posicion.txt", header = TRUE, sep = "\t")

# Contar genes por módulo y grupo cromosómico
modulo_chr <- datos %>%
  group_by(modulo, CHR) %>%
  summarise(cantidad_genes = n()) %>%
  ungroup()

# Visualizar resultados
ggplot(modulo_chr, aes(x = factor(CHR), y = cantidad_genes, fill = modulo)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Grupo de ligamiento (CHR)", y = "Cantidad de genes", fill = "Módulo") +
  theme_minimal()

# Crear una tabla de contingencia módulo vs CHR
tabla_contingencia <- table(datos$modulo, datos$CHR)

# Realizar prueba Chi-cuadrado
chi_resultado <- chisq.test(tabla_contingencia)

# Ver resultados
chi_resultado







