#!/usr/bin/env bash
# config.sh (fallback if you don't want YAML)

SCRIPTS_DIR="./scripts"
OUT_ROOT="./runs"
RUN_NAME="latest"

COUNT_MATRIX_CSV="./data/conteo_data_sin_duplicados.csv"
METADATA_CSV="./data/muestra_info.csv"
LINKAGE_MAP_TSV="./data/combined_map_data.txt"
EXPRESSION_LONG_TSV="./data/allgenesorder.txt"

DO_DESEQ2=1
DO_WGCNA_COEXP=1
DO_COEXP_FDR=1
DO_CORRELACION_CORTAS=1

# Optional patch parameters
N_MUESTRAS=10
EXPRESSION_THRESHOLD=500
CM_BINS_JSON='[1,5]'

# Optional concat
DO_CONCAT_EXPRESSION=0
EXPRESSION_FILES_DIR="./data/expression_files"

# Optional BAM counting
DO_COUNT_FROM_BAM=0
BAM_DIR="./data/bams"
SAMTOOLS_BIN="samtools"
SAMPLES_JSON='["DC01_15_sw10L50","SC01_15_sw10L50"]'

# Optional mapping steps
DO_BUSQUEDA_LG=0
LOCUS_LIST_TXT="./data/locusLG12.txt"
BLAST_TSV="./data/e-30.blastn.txt"

DO_COORDENADAS_GEN=0
ORDER_MAPPED_TSV="./data/order12_10PS.mapped"
RESULTADOS_LG_TSV="./data/resultadoslg12.txt"

DO_EXTRAER_POS_GEN=0
COORDS_DIR="./data/coordenadas_genesTransc"

DO_OBTENER_LOCUS_GEN_POS=0
DO_POISSON=0
