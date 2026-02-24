# Integrating high-density linkage maps and transcriptomics reveals ozone‐response expression modules in a peri-urban forest tree

This repository documents the workflow used to construct a saturated linkage map in *Abies religiosa* and integrate it with transcriptomic data to identify genes involved in ozone stress response.

## Quick start (recommended)

This repo includes a small runner that executes the main downstream analyses (DESeq2 → WGCNA/co-expression → FDR/network → Pearson-by-distance) **without editing your original scripts**. The runner stages your inputs into a run directory using the exact filenames expected by your scripts.

### 1) Install the environment

```bash
conda env create -f environment.yml
conda activate abies-linkage-transcriptome
```

### 2) Configure inputs

Edit `config.yaml` and set the paths to your files:

- `count_matrix_csv`: `conteo_data_sin_duplicados.csv`
- `metadata_csv`: `muestra_info.csv`
- `linkage_map_tsv`: `combined_map_data.txt`
- `expression_long_tsv`: `allgenesorder.txt` (long format: `expresion`, `GENE`, `muestra`)

### 3) Run

```bash
./run_pipeline.sh -c config.yaml
```

Outputs and logs are written to:

```text
runs/<run_name>/
  logs/
  (script outputs)
```

> Tip: you can turn steps on/off in `config.yaml` using `do_*` flags.

---

## Repository layout

```text
.
├── scripts/                 # analysis scripts (.sh, .R, .py)
├── run_pipeline.sh          # pipeline runner (stages inputs, runs steps)
├── config.yaml              # config (paths + step toggles)
├── environment.yml          # conda environment
├── CITATION.cff             # citation metadata (GitHub “Cite this repository”)
└── LICENSE
```

---

## Inputs expected by downstream analyses (runner)

The runner creates a working directory and stages these filenames (symlink/copy):

| Staged filename in run dir | Meaning |
|---|---|
| `conteo_data_sin_duplicados.csv` | count matrix (rows=genes, cols=samples) |
| `muestra_info.csv` | sample metadata (must include `muestra`, `estado`, `condicion`) |
| `combined_map_data.txt` | linkage map table (must include `GENE`, `CHR`, `POS_avg`) |
| `allgenesorder.txt` | long expression table (columns: `expresion`, `GENE`, `muestra`) |

---

## Main scripts (downstream)

### Differential expression (DESeq2)
- **Script:** `scripts/deseq_bien_4comparaciones.R`
- **Inputs (expected filenames):**  
  - `conteo_data_sin_duplicados.csv`  
  - `muestra_info.csv`
- **Outputs:** DESeq2 results tables and diagnostic plots (see script for filenames).

### Co-expression (WGCNA) + correlation by CHR and distance bins
- **Script:** `scripts/Coexpresion.R`
- **Inputs (expected filenames):**  
  - `conteo_data_sin_duplicados.csv`  
  - `combined_map_data.txt`
- **Key outputs (examples):**  
  - `genes_por_modulo_WGCNA.txt`  
  - `modulos_con_posicion.txt`  
  - `correlaciones_pares_por_CHR.txt`  
  - `correlacion_por_CHR_1cM_5cM.txt`

### FDR correction + co-expression network
- **Script:** `scripts/Coexp_FDR.R`
- **Inputs (expected filenames):**  
  - `correlaciones_pares_por_CHR.txt`  
  - `combined_map_data.txt`  
  - `muestra_info.csv` (used for presence-by-condition summaries)
- **Key outputs (examples):**  
  - `correlaciones_significativas_positivas.txt`  
  - `genes_correlacionados_por_CHR.txt`  
  - `genes_estado_presencia.txt`

### Pearson correlations by distance bins (Python)
- **Script:** `scripts/correlacion_cortas.py`
- **Inputs (expected filenames):**  
  - `allgenesorder.txt`  
  - `combined_map_data.txt`
- **Outputs:** binned correlation summaries + plots (see script for filenames).  
  The runner can optionally patch `expression_threshold` and the list of `cm` bins via `config.yaml`.

---

## Manual workflow (full map + annotation pipeline)

The sections below preserve the original “manual” commands used in the study.  
**Important:** paths like `~/lep-map3/bin` are system-specific; update them for your installation.

### 🧬 Linkage Map Construction using Lep-MAP3

#### 1. ParentCall2
```bash
java -cp ~/lep-map3/bin ParentCall2 data=pedigreeAll_PR.txt vcfFile=populations.snps.vcf removeNonInformative=1 halfSibs=1 > data.call
```

#### 2. Filtering
```bash
java -cp ~/lep-map3/bin Filtering2 data=data.call removeNonInformative=1 dataTolerance=0.001 heterozygoteRate=0.05 > data_fil.call
```

#### 3. Identification of Linkage Groups
```bash
java -cp ~/lep-map3/bin SeparateChromosomes2 data=data_fil.call lodLimit=30 sizeLimit=20 distortionLod=1 numThreads=12 > map30.txt
```

#### 4. Add Unique Markers to LGs
```bash
java -cp ~/lep-map3/bin JoinSingles2All map=map30.txt data=data_fil.call lodLimit=10 lodDifference=5 > map30_js.txt
```

#### 5. Construct the Linkage Map
```bash
java -cp ~/lep-map3/bin OrderMarkers2 map=map30_js.txt data=data_fil.call selfingPhase=1 chromosome=1 numThreads=16 numMergeIterations=15 numPolishIterations=5 proximityScale=100 > order1PS.txt
```

#### 6. Generate Dot Plot
```bash
java -cp ~/lep-map3/bin LMPlot order1PS.txt > order1PS.dot
```

---

### 🧮 LOD Score Matrix with Lep-MAP2

#### 1. Convert and Transpose
```bash
awk -f simpleConvert.awk dataAll_fil.call | cut -f 3- | java -cp ~/lep_map/bin Transpose > data_filter.linkage
```

#### 2. Compute LOD Scores
```bash
java -cp ~/lep_map/bin OrderMarkers data=data_filter.linkage evaluateOrder=order1_2.txt improveOrder=0 computeLODScores=1 > lod1_2.txt
```

#### 3. Visualize LOD Scores with Gnuplot
```gnuplot
set terminal png size 2000,2000
set output "lod1PS.png"
set pm3d map
set zrange [0:35]
splot "< cut -f 2- lod1PS.txt | awk '(NR>2)'" matrix
```

---

### 🧬 Marker Matching and BLAST Analysis

#### 1. Extract Loci from VCF and Sequences
```r
library(vcfR); library(dplyr)
lg1 <- read.vcfR("populations.snps1_10.recode.vcf")
locusLG1 <- as.data.frame(lg1@fix)[, "CHROM", drop=FALSE]
write.table(locusLG1, "locusLG1.txt", sep="\t", row.names=FALSE, col.names=FALSE)
```

```bash
seqtk subseq populations.loci.fa locusLG1.txt > locusLG1.fa ... locusLG2.fa ...
cat *.fa > allLG.fa
```

#### 2. BLAST Against Transcriptome
```bash
makeblastdb -in GCAT_AB-RNA-1.0.16.fa -dbtype nucl -parse_seqids
blastn -task blastn -query allLG.fa -db GCAT_AB-RNA-1.0.16.fa -evalue 1e-30 -out e-30.blastn.txt -outfmt 6
```

---

## 🔬 Functional Annotation and Expression

This section describes how loci from the linkage map were matched to transcriptomic sequences and annotated functionally.

### 1. Identify matched loci from BLAST results
```bash
awk '{print $1}' e-30.blastn.txt | uniq > loci_match.txt
```

### 2. Extract sequences of matched loci from the genome
```bash
xargs samtools faidx allLG.fa < loci_match.txt >> matches_loci.fa
```

### 3. Extract sequences of target transcriptome genes
```bash
xargs samtools faidx GCAT_AB-RNA-1.0.16.fa < genes_exp_diferencial2.txt >> matches_genes_exp_dif.fa
```

### 4. Concatenate both locus and gene sequences
```bash
cat matches_loci.fa matches_genes_exp_dif.fa > matches_duplicated_loci.fa
```

### 5. Remove duplicated FASTA entries
```bash
awk '/^>/'"{if(seq) print seq; print; seq=\"\"; next} {seq=seq $0} END{if(seq) print seq}' matches_duplicated_loci.fa | awk '!seen[$0]++' > matches_duplicated_loci_nodup.fasta
```

### 6. Annotate matched loci using InterProScan
```bash
./interproscan.sh -cpu 12 -t n -i matches_loci.fa -b lociLG_blast -goterms
```

---

## 🧬 Integration of Expression Data (manual)

> These helper scripts may require editing if they contain fixed paths. The runner provides portable alternatives for several steps.

#### 7. Extract coordinates of loci and matched genes
```bash
bash scripts/coordenadas_gen.sh
```

#### 8. Collect expression data per sample (add sample ID column)
```bash
bash scripts/agregar_col.sh <sample_file.tsv>
```

#### 9. Concatenate all expression files
```bash
bash scripts/concatenar.sh
```

#### 10. Extract gene name and map position
```bash
bash scripts/extraer_pos_gen.sh
```

#### 11. Combine gene–locus–position data
```bash
python scripts/obtener_locus_gen_pos.py
```

---

## 📊 Genetic map visualization
```bash
./genetic_mapper.pl --bar --map=allLGFig.txt > all.svg
```

---

## 📈 Marker distribution (Poisson / NB)
```bash
Rscript scripts/poisson_ultimo.R
```

---

## Data availability

Raw sequencing data are available at NCBI under BioProject **PRJNA1314304**
(*Genotyping-by-sequencing for linkage map construction in Abies religiosa*).

This repository hosts the analysis pipeline, a small `data/example/` dataset for testing, and sequence-defined resources that enable reuse of the linkage map across studies by matching loci via sequence similarity (rather than dataset-specific marker names).

### Linkage-map loci and annotations
- `data/linkage_map_loci/matches_loci.fa` — marker-associated locus sequences used in the linkage map (FASTA; STACKS locus IDs retained)
- `data/annotation/lociLG_blast.gff3` — functional annotations for mapped loci (GFF3; InterPro-based)

**Typical uses:** (i) within-*A. religiosa* map-to-map comparisons by sequence matching; (ii) map-to-assembly anchoring/validation by aligning scaffolds/contigs to `matches_loci.fa`.

## 🧾 Citation

If you use this pipeline, please cite our manuscript:

> Granados-Aguilar X. et al. (202X). *Integrating high-density linkage maps and transcriptomics reveals ozone‐response expression modules in a peri-urban forest tree*. [In review].

And/or cite the software release via **GitHub → “Cite this repository”** (powered by `CITATION.cff`).
