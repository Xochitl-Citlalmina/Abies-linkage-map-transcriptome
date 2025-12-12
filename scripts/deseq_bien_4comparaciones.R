library(DESeq2)
library(ggplot2)
library("pheatmap")
library(qqman)

# Leer y preparar countData
countData <- read.csv("conteo_data_sin_duplicados.csv", header = TRUE)
rownames(countData) <- countData$GENE
countData <- countData[ , -1]  # quitar columna GENE
### Leer y preparar metadata
metaData <- read.csv("muestra_info.csv", header = TRUE)

# Revisa los niveles de condición y estado
unique(metaData$condicion)
unique(metaData$estado)

# Asegura formato uniforme
metaData$condicion <- trimws(tolower(metaData$condicion))
metaData$estado <- trimws(tolower(metaData$estado))

# Filtrar muestras con condición "alta"
samples_alta <- metaData$muestra[metaData$condicion == "alta"]

# Verifica que hay suficientes muestras
print(samples_alta)

# Subconjunto de conteos y metadata
countData_alta <- countData[ , colnames(countData) %in% samples_alta]
metaData_alta <- metaData[match(colnames(countData_alta), metaData$muestra), ]

# Verificar coincidencia de orden
stopifnot(all(colnames(countData_alta) == metaData_alta$muestra))

###### Dañado vs sano en alta concentración de ozono

samples_alta <- metaData$muestra[metaData$condicion == "alta"]
countData_alta <- countData[ , samples_alta]
metaData_alta <- metaData[metaData$muestra %in% samples_alta, ]
metaData_alta <- metaData_alta[match(colnames(countData_alta), metaData_alta$muestra), ]

dds_alta <- DESeqDataSetFromMatrix(countData = countData_alta,
                                   colData = metaData_alta,
                                   design = ~ estado)
dds_alta <- dds_alta[rowSums(counts(dds_alta)) >= 10, ]
dds_alta <- DESeq(dds_alta)
res_alta <- results(dds_alta, contrast = c("estado", "dañado", "sano"))

####### 2. Dañado vs sano en concentración moderada de ozono

samples_moderada <- metaData$muestra[metaData$condicion == "moderada"]
countData_moderada <- countData[ , samples_moderada]
metaData_moderada <- metaData[metaData$muestra %in% samples_moderada, ]
metaData_moderada <- metaData_moderada[match(colnames(countData_moderada), metaData_moderada$muestra), ]

dds_moderada <- DESeqDataSetFromMatrix(countData = countData_moderada,
                                       colData = metaData_moderada,
                                       design = ~ estado)
dds_moderada <- dds_moderada[rowSums(counts(dds_moderada)) >= 10, ]
dds_moderada <- DESeq(dds_moderada)
res_moderada <- results(dds_moderada, contrast = c("estado", "dañado", "sano"))

######## 3. Dañado vs sano sin importar concentración

dds_estado <- DESeqDataSetFromMatrix(countData = countData,
                                     colData = metaData,
                                     design = ~ estado)
dds_estado <- dds_estado[rowSums(counts(dds_estado)) >= 10, ]
dds_estado <- DESeq(dds_estado)
res_estado <- results(dds_estado, contrast = c("estado", "dañado", "sano"))

###### 4. Sanos y dañados bajo diferentes concentraciones (alta vs moderada)

dds_full <- DESeqDataSetFromMatrix(countData = countData,
                                   colData = metaData,
                                   design = ~ estado + condicion)
dds_full <- dds_full[rowSums(counts(dds_full)) >= 10, ]
dds_full <- DESeq(dds_full)

# Comparación entre condiciones
res_condicion <- results(dds_full, contrast = c("condicion", "alta", "moderada"))

# Comparación entre estado ajustada por condición
res_estado_ajustado <- results(dds_full, contrast = c("estado", "dañado", "sano"))

#### guardar en csv
write.csv(as.data.frame(res_alta), "res_alta_dañado_vs_sano.csv")
write.csv(as.data.frame(res_moderada), "res_moderada_dañado_vs_sano.csv")
write.csv(as.data.frame(res_estado), "res_todo_dañado_vs_sano.csv")
write.csv(as.data.frame(res_condicion), "res_condicion_alta_vs_moderada.csv")


##### Heatmap

library(pheatmap)
library(RColorBrewer)
library(DESeq2)

# Obtener los genes más significativos
topgenes <- head(order(res_alta$padj, na.last = NA), 20)
vsd <- vst(dds_alta, blind = FALSE)

# Crear anotación de columnas
anno_col <- data.frame(Group = colData(dds_alta)$estado)
rownames(anno_col) <- colnames(dds_alta)

# Traducir etiquetas a inglés
anno_col$Group <- factor(anno_col$Group,
                         levels = c("sano", "dañado"),
                         labels = c("healthy", "damaged"))

# Paleta de colores para la anotación
ann_colors <- list(Group = c(healthy = "#00D9D9", damaged = "#E60000"))

# Generar el heatmap
pheatmap(assay(vsd)[topgenes, ],
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = anno_col,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize = 10,
         fontsize_row = 8,
         fontsize_col = 9,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Damaged vs. healthy in high ozone")

library(grid)

# Guardar el heatmap en un objeto
p <- pheatmap(assay(vsd)[topgenes, ],
              scale = "row",
              cluster_rows = TRUE,
              cluster_cols = TRUE,
              annotation_col = anno_col,
              annotation_colors = ann_colors,
              show_rownames = TRUE,
              show_colnames = TRUE,
              fontsize = 10,
              fontsize_row = 8,
              fontsize_col = 9,
              color = colorRampPalette(c("blue", "white", "red"))(50),
              main = "Damaged vs. healthy in high ozone")

# Exportar como TIFF
tiff("heatmap_high_O3.tiff", width = 2000, height = 2000, res = 300)
grid::grid.newpage()
grid::grid.draw(p$gtable)
dev.off()

#### heatmap otros 3 
# Función general para generar y exportar heatmaps
crear_heatmap <- function(dds_obj, res_obj, nombre_archivo, titulo_plot) {
  # Selección de los 20 genes más significativos
  topgenes <- head(order(res_obj$padj, na.last = NA), 80)
  vsd <- vst(dds_obj, blind = FALSE)
  
  # Anotaciones en inglés
  anno_col <- data.frame(Group = colData(dds_obj)$estado)
  rownames(anno_col) <- colnames(dds_obj)
  anno_col$Group <- factor(anno_col$Group,
                           levels = c("sano", "dañado"),
                           labels = c("healthy", "damaged"))
  ann_colors <- list(Group = c(healthy = "#00D9D9", damaged = "#E60000"))
  
  # Crear objeto de pheatmap sin imprimir aún
  p <- pheatmap(assay(vsd)[topgenes, ],
                scale = "row",
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                annotation_col = anno_col,
                annotation_colors = ann_colors,
                show_rownames = TRUE,
                show_colnames = TRUE,
                fontsize = 10,
                fontsize_row = 8,
                fontsize_col = 9,
                color = colorRampPalette(c("blue", "white", "red"))(50),
                main = titulo_plot)
  
  # Exportar a TIFF
  tiff(nombre_archivo, width = 2000, height = 2000, res = 300)
  grid.newpage()
  grid.draw(p$gtable)
  dev.off()
}
### Heatmap 1: Damaged vs. healthy in high ozone
crear_heatmap(dds_alta, res_alta, "heatmap_high_O3.tiff", "Damaged vs. healthy in high ozone")
#### Heatmap 2: Damaged vs. healthy in moderate ozone
crear_heatmap(dds_moderada, res_moderada, "heatmap_moderate_O3.tiff", "Damaged vs. healthy in moderate ozone")

### Heatmap 3: Damaged vs. healthy (all samples)
crear_heatmap(dds_estado, res_estado, "heatmap_all_samples.tiff", "Damaged vs. healthy (all samples)")

### Heatmap 4: High vs. moderate ozone (ajustado por estado)
crear_heatmap_condicion <- function(dds_obj, res_obj, nombre_archivo, titulo_plot) {
  topgenes <- head(order(res_obj$padj, na.last = NA), 80)
  vsd <- vst(dds_obj, blind = FALSE)
  
  anno_col <- data.frame(Condition = colData(dds_obj)$condicion)
  rownames(anno_col) <- colnames(dds_obj)
  anno_col$Condition <- factor(anno_col$Condition,
                               levels = c("moderada", "alta"),
                               labels = c("moderate O3", "high O3"))
  ann_colors <- list(Condition = c(`moderate O3` = "#91D1C2", `high O3` = "#FF914D"))
  
  p <- pheatmap(assay(vsd)[topgenes, ],
                scale = "row",
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                annotation_col = anno_col,
                annotation_colors = ann_colors,
                show_rownames = TRUE,
                show_colnames = TRUE,
                fontsize = 10,
                fontsize_row = 8,
                fontsize_col = 9,
                color = colorRampPalette(c("blue", "white", "red"))(50),
                main = titulo_plot)
  
  tiff(nombre_archivo, width = 2000, height = 2000, res = 300)
  grid.newpage()
  grid.draw(p$gtable)
  dev.off()
}

crear_heatmap_condicion(dds_full, res_condicion, "heatmap_condicion_O3.tiff",
                        "High vs. moderate ozone (all trees)")




