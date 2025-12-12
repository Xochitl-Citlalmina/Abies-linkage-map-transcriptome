
library(dplyr)
# Cargar datos
cor_data <- read.table("correlaciones_pares_por_CHR.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Calcular p-valores para las correlaciones
# Vamos a estimar los p-valores a partir de los coeficientes de correlación

n_muestras <- 10  

# Transformar r a p-valor (t-distribución)
cor_data$pvalue <- 2 * pt(-abs(cor_data$correlacion * sqrt((n_muestras - 2) / (1 - cor_data$correlacion^2))), df = n_muestras - 2)

# Ajustar por FDR
cor_data$fdr <- p.adjust(cor_data$pvalue, method = "fdr")

# Filtrar correlaciones positivas y significativas
signif_cor <- cor_data %>%
  filter(correlacion > 0, fdr < 0.05)

# Ver cuántos pares quedaron
cat("Número de pares significativos positivos:\n")
nrow(signif_cor)

# Obtener genes únicos involucrados
genes_involucrados <- unique(c(signif_cor$gene1, signif_cor$gene2))
cat("Número de genes únicos correlacionados:\n")
length(genes_involucrados)

# leer mapa
mapa <- read.table("combined_map_data.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Obtener CHR de los genes correlacionados
genes_chr <- mapa %>%
  filter(GENE %in% genes_involucrados) %>%
  select(GENE, CHR)

# Resumen por CHR
resumen_chr <- genes_chr %>%
  group_by(CHR) %>%
  summarise(n_genes = n_distinct(GENE))

# Ver resultado
print(resumen_chr)

# Guardar archivos
write.table(signif_cor, file = "correlaciones_significativas_positivas.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(genes_chr, file = "genes_correlacionados_por_CHR.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Instalar igraph si no lo tienes
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
library(igraph)

# Crear grafo
red <- graph_from_data_frame(
  signif_cor[, c("gene1", "gene2", "correlacion")],
  directed = FALSE
)

# Opcional: asignar CHR a los nodos
genes_chr_unicos <- genes_chr[!duplicated(genes_chr$GENE), ]
V(red)$CHR <- genes_chr_unicos$CHR[match(V(red)$name, genes_chr_unicos$GENE)]

# Colorear por CHR
colores_chr <- rainbow(length(unique(V(red)$CHR)))
names(colores_chr) <- unique(V(red)$CHR)
V(red)$color <- colores_chr[V(red)$CHR]

# Tamaño de nodo según número de conexiones (grado)
V(red)$size <- degree(red) * 3 + 5

# Dibujar red
plot(
  red,
  vertex.label = V(red)$name,
  vertex.label.cex = 0.7,
  edge.width = signif_cor$correlacion * 5,  # grosor según fuerza
  layout = layout_with_fr,
  main = "Co-expression network (positive significant correlations)"
)

######## Saber que muestras son sanas y dañadas

# 1. Cargar expresión y metadata
expresion <- read.csv("conteo_data_sin_duplicados.csv", header = TRUE, row.names = 1, check.names = FALSE)  # o .txt con sep = "\t"

# Transponer si tiene genes en filas
expresion <- t(expresion)

# 2. Cargar metadata
meta <- read.csv("muestra_info.csv", header = TRUE, stringsAsFactors = FALSE)

# 3. Cargar genes correlacionados
correlados <- read.table("correlaciones_significativas_positivas.txt", header = TRUE, sep = "\t")
genes_cor <- unique(c(correlados$gene1, correlados$gene2))

# 4. Subset de la matriz de expresión con esos genes
expr_sub <- expresion[, colnames(expresion) %in% genes_cor]

# 5. Revisar si los genes tienen expresión en muestras dañadas o sanas

# Agregar etiqueta de estado (sano/dañado)
meta_estado <- meta[, c("muestra", "estado")]
rownames(meta_estado) <- meta_estado$muestra

# Crear resumen por gen
resumen_genes <- data.frame(GENE = colnames(expr_sub))

resumen_genes$estado <- sapply(colnames(expr_sub), function(gene) {
  expr_g <- expr_sub[, gene]
  muestras_expresadas <- names(expr_g[expr_g > 0])
  estados <- unique(meta_estado[muestras_expresadas, "estado"])
  if (length(estados) == 1) {
    return(estados)
  } else if (length(estados) > 1) {
    return("ambos")
  } else {
    return("no_expresado")
  }
})

# Resultado final
print(resumen_genes)

# Guardar
write.table(resumen_genes, file = "genes_estado_presencia.txt", sep = "\t", row.names = FALSE, quote = FALSE)











