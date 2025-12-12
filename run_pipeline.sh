#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./run_pipeline.sh              # uses config.yaml if available, else config.sh
#   ./run_pipeline.sh -c config.yaml
#
# Notes:
# - This runner *stages* the inputs into a working directory so your existing scripts
#   (which expect fixed filenames like conteo_data_sin_duplicados.csv) can run unmodified.

CONFIG=""
WORKDIR=""
while getopts ":c:w:h" opt; do
  case "$opt" in
    c) CONFIG="$OPTARG" ;;
    w) WORKDIR="$OPTARG" ;;
    h)
      echo "Usage: $0 [-c config.yaml|config.sh] [-w workdir]"
      exit 0
      ;;
    \?)
      echo "Unknown option: -$OPTARG" >&2
      exit 2
      ;;
  esac
done

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Default config choice
if [[ -z "${CONFIG}" ]]; then
  if [[ -f "${ROOT_DIR}/config.yaml" ]]; then
    CONFIG="${ROOT_DIR}/config.yaml"
  elif [[ -f "${ROOT_DIR}/config.sh" ]]; then
    CONFIG="${ROOT_DIR}/config.sh"
  else
    echo "ERROR: No config found. Create config.yaml (preferred) or config.sh in repo root." >&2
    exit 1
  fi
fi

# ---------- helpers ----------
log() { echo "[$(date '+%F %T')] $*"; }

die() { echo "ERROR: $*" >&2; exit 1; }

have() { command -v "$1" >/dev/null 2>&1; }

# Try symlink; fallback to copy (useful on systems without symlink perms)
stage_file() {
  local src="$1"
  local dst="$2"
  [[ -n "$src" ]] || die "stage_file: empty src for $dst"
  [[ -f "$src" ]] || die "Missing input file: $src"
  if ln -sf "$src" "$dst" 2>/dev/null; then
    :
  else
    cp -f "$src" "$dst"
  fi
}

# ---------- load config ----------
load_from_yaml() {
  local yaml_path="$1"
  python3 - "$yaml_path" <<'PY'
import sys, os
cfg_path = sys.argv[1]
try:
    import yaml  # type: ignore
except Exception:
    sys.stderr.write("NO_PYYAML\n")
    sys.exit(13)

with open(cfg_path, "r", encoding="utf-8") as f:
    cfg = yaml.safe_load(f) or {}

# Minimal flattening: only top-level keys (strings, ints, bools).
# Nested dicts are exported as JSON strings.
def bash_escape(s: str) -> str:
    return s.replace("\\", "\\\\").replace('"', '\\"')

for k, v in cfg.items():
    key = str(k).upper()
    if isinstance(v, (str, int, float)):
        print(f'export {key}="{bash_escape(str(v))}"')
    elif isinstance(v, bool):
        print(f'export {key}={"1" if v else "0"}')
    elif v is None:
        print(f'export {key}=""')
    else:
        import json
        print(f'export {key}="{bash_escape(json.dumps(v))}"')
PY
}

if [[ "${CONFIG}" == *.yml || "${CONFIG}" == *.yaml ]]; then
  if ! have python3; then die "python3 is required to read YAML config"; fi
  set +e
  EXPORTS="$(load_from_yaml "${CONFIG}")"
  status=$?
  set -e
  if [[ $status -eq 13 ]]; then
    die "PyYAML is not installed. Install with: python3 -m pip install pyyaml  (or use config.sh instead)."
  elif [[ $status -ne 0 ]]; then
    die "Failed to parse YAML config: ${CONFIG}"
  fi
  eval "${EXPORTS}"
else
  # shellcheck source=/dev/null
  source "${CONFIG}"
fi

# ---------- defaults ----------
: "${SCRIPTS_DIR:=${ROOT_DIR}/scripts}"
: "${OUT_ROOT:=${ROOT_DIR}/runs}"
: "${RUN_NAME:=latest}"

if [[ -z "${WORKDIR}" ]]; then
  WORKDIR="${OUT_ROOT}/${RUN_NAME}"
fi

mkdir -p "${WORKDIR}"
mkdir -p "${WORKDIR}/logs"

log "Repo root: ${ROOT_DIR}"
log "Config: ${CONFIG}"
log "Workdir: ${WORKDIR}"

# ---------- stage required inputs ----------
# These scripts expect *exact* filenames; we stage them accordingly.
if [[ -n "${COUNT_MATRIX_CSV:-}" ]]; then
  stage_file "${COUNT_MATRIX_CSV}" "${WORKDIR}/conteo_data_sin_duplicados.csv"
fi
if [[ -n "${METADATA_CSV:-}" ]]; then
  stage_file "${METADATA_CSV}" "${WORKDIR}/muestra_info.csv"
fi
if [[ -n "${LINKAGE_MAP_TSV:-}" ]]; then
  stage_file "${LINKAGE_MAP_TSV}" "${WORKDIR}/combined_map_data.txt"
fi
if [[ -n "${EXPRESSION_LONG_TSV:-}" ]]; then
  stage_file "${EXPRESSION_LONG_TSV}" "${WORKDIR}/allgenesorder.txt"
fi

# ---------- optional: build allgenesorder.txt from a folder ----------
if [[ "${DO_CONCAT_EXPRESSION:-0}" == "1" ]]; then
  [[ -n "${EXPRESSION_FILES_DIR:-}" ]] || die "DO_CONCAT_EXPRESSION=1 but EXPRESSION_FILES_DIR is empty"
  [[ -d "${EXPRESSION_FILES_DIR}" ]] || die "EXPRESSION_FILES_DIR does not exist: ${EXPRESSION_FILES_DIR}"

  log "Concatenating expression files from: ${EXPRESSION_FILES_DIR}"
  out="${WORKDIR}/allgenesorder.txt"
  rm -f "$out"

  first=1
  shopt -s nullglob
  files=( "${EXPRESSION_FILES_DIR}"/*.txt "${EXPRESSION_FILES_DIR}"/*.tsv )
  shopt -u nullglob
  [[ ${#files[@]} -gt 0 ]] || die "No .txt/.tsv files found in ${EXPRESSION_FILES_DIR}"

  for f in "${files[@]}"; do
    if [[ $first -eq 1 ]]; then
      cat "$f" >> "$out"
      first=0
    else
      tail -n +2 "$f" >> "$out"
    fi
  done
  log "Wrote: ${out}"
fi

# ---------- optional: count per-contig from BAMs (portable replacement for count_genes_bamfile.sh) ----------
if [[ "${DO_COUNT_FROM_BAM:-0}" == "1" ]]; then
  [[ -n "${BAM_DIR:-}" ]] || die "DO_COUNT_FROM_BAM=1 but BAM_DIR is empty"
  [[ -d "${BAM_DIR}" ]] || die "BAM_DIR does not exist: ${BAM_DIR}"
  [[ -n "${SAMPLES_JSON:-}" ]] || die "DO_COUNT_FROM_BAM=1 but SAMPLES_JSON is empty (list of sample IDs)"
  have "${SAMTOOLS_BIN:-samtools}" || die "samtools not found (set SAMTOOLS_BIN in config)"

  log "Counting contigs from BAMs (portable step)"
  python3 - <<'PY' "${BAM_DIR}" "${WORKDIR}" "${SAMTOOLS_BIN:-samtools}" "${SAMPLES_JSON}"
import sys, os, json, subprocess, pathlib, collections
bam_dir, workdir, samtools, samples_json = sys.argv[1:]
samples = json.loads(samples_json)

outdir = pathlib.Path(workdir) / "bam_counts"
outdir.mkdir(parents=True, exist_ok=True)

for s in samples:
    bam = pathlib.Path(bam_dir) / f"{s}.bam"
    if not bam.exists():
        raise SystemExit(f"Missing BAM: {bam}")
    # Stream alignments -> count RNAME (3rd field in SAM)
    cmd = [samtools, "view", str(bam)]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    counts = collections.Counter()
    assert p.stdout is not None
    for line in p.stdout:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3:
            counts[parts[2]] += 1
    rc = p.wait()
    if rc != 0:
        raise SystemExit(f"samtools view failed for {bam} (exit {rc})")

    out = outdir / f"{s}.countgenes.tsv"
    with out.open("w", encoding="utf-8") as fh:
        for contig, n in counts.most_common():
            fh.write(f"{n}\t{contig}\n")

print(f"Wrote per-sample count lists to: {outdir}")
PY
fi

# ---------- run DESeq2 ----------
if [[ "${DO_DESEQ2:-1}" == "1" ]]; then
  have Rscript || die "Rscript not found"
  [[ -f "${SCRIPTS_DIR}/deseq_bien_4comparaciones.R" ]] || die "Missing script: ${SCRIPTS_DIR}/deseq_bien_4comparaciones.R"
  log "Running DESeq2 script..."
  ( cd "${WORKDIR}" && Rscript "${SCRIPTS_DIR}/deseq_bien_4comparaciones.R" ) |& tee "${WORKDIR}/logs/deseq2.log"
fi

# ---------- run WGCNA + correlations ----------
if [[ "${DO_WGCNA_COEXP:-1}" == "1" ]]; then
  have Rscript || die "Rscript not found"
  [[ -f "${SCRIPTS_DIR}/Coexpresion.R" ]] || die "Missing script: ${SCRIPTS_DIR}/Coexpresion.R"
  log "Running Coexpresion.R (WGCNA + correlations)..."
  ( cd "${WORKDIR}" && Rscript "${SCRIPTS_DIR}/Coexpresion.R" ) |& tee "${WORKDIR}/logs/coexpresion.log"
fi

# ---------- run FDR over correlations (patch n_muestras if requested) ----------
if [[ "${DO_COEXP_FDR:-1}" == "1" ]]; then
  have Rscript || die "Rscript not found"
  [[ -f "${SCRIPTS_DIR}/Coexp_FDR.R" ]] || die "Missing script: ${SCRIPTS_DIR}/Coexp_FDR.R"

  tmp="${WORKDIR}/Coexp_FDR.patched.R"
  cp -f "${SCRIPTS_DIR}/Coexp_FDR.R" "$tmp"

  if [[ -n "${N_MUESTRAS:-}" ]]; then
    log "Patching Coexp_FDR.R: n_muestras <- ${N_MUESTRAS}"
    # Replace the first occurrence of "n_muestras <- <number>"
    perl -0777 -i -pe 's/n_muestras\s*<-\s*\d+\s*/n_muestras <- '"${N_MUESTRAS}"'\n/s' "$tmp"
  fi

  log "Running Coexp_FDR (FDR + network)..."
  ( cd "${WORKDIR}" && Rscript "$tmp" ) |& tee "${WORKDIR}/logs/coexp_fdr.log"
fi

# ---------- run Pearson correlations by distance bins (Python) ----------
if [[ "${DO_CORRELACION_CORTAS:-1}" == "1" ]]; then
  have python3 || die "python3 not found"
  [[ -f "${SCRIPTS_DIR}/correlacion_cortas.py" ]] || die "Missing script: ${SCRIPTS_DIR}/correlacion_cortas.py"

  tmp="${WORKDIR}/correlacion_cortas.patched.py"
  cp -f "${SCRIPTS_DIR}/correlacion_cortas.py" "$tmp"

  if [[ -n "${EXPRESSION_THRESHOLD:-}" ]]; then
    log "Patching correlacion_cortas.py: expression_threshold=${EXPRESSION_THRESHOLD}"
    perl -0777 -i -pe 's/expression_threshold\s*=\s*\d+/expression_threshold='"${EXPRESSION_THRESHOLD}"'/s' "$tmp"
  fi
  if [[ -n "${CM_BINS_JSON:-}" ]]; then
    log "Patching correlacion_cortas.py: cm bins=${CM_BINS_JSON}"
    python3 - "$tmp" "$CM_BINS_JSON" <<'PY'
import sys, json, re, pathlib
path, bins_json = sys.argv[1], sys.argv[2]
bins = json.loads(bins_json)
text = pathlib.Path(path).read_text(encoding="utf-8")
# Replace: for cm in [1, 5]:
text = re.sub(r"for\s+cm\s+in\s+\[[^\]]+\]\s*:", f"for cm in {bins}:", text)
pathlib.Path(path).write_text(text, encoding="utf-8")
PY
  fi

  log "Running correlacion_cortas.py..."
  ( cd "${WORKDIR}" && python3 "$tmp" ) |& tee "${WORKDIR}/logs/correlacion_cortas.log"
fi

# ---------- optional: LG BLAST lookup (portable replacement for BusquedaLGotro.sh) ----------
if [[ "${DO_BUSQUEDA_LG:-0}" == "1" ]]; then
  [[ -n "${LOCUS_LIST_TXT:-}" ]] || die "DO_BUSQUEDA_LG=1 but LOCUS_LIST_TXT is empty"
  [[ -n "${BLAST_TSV:-}" ]] || die "DO_BUSQUEDA_LG=1 but BLAST_TSV is empty"
  stage_file "${LOCUS_LIST_TXT}" "${WORKDIR}/locusLG.txt"
  stage_file "${BLAST_TSV}" "${WORKDIR}/blast.txt"
  out="${WORKDIR}/resultados_lg.txt"
  log "Running BLAST locus lookup -> ${out}"
  python3 - <<'PY' "${WORKDIR}/locusLG.txt" "${WORKDIR}/blast.txt" "${out}"
import sys
locus, blast, out = sys.argv[1:]
wanted = set(x.strip() for x in open(locus, encoding="utf-8") if x.strip())
with open(blast, encoding="utf-8") as f, open(out, "w", encoding="utf-8") as o:
    for line in f:
        if not line.strip(): 
            continue
        parts = line.rstrip("\n").split()
        if parts and parts[0] in wanted:
            # mimic your awk '{print $1, $2}'
            o.write(f"{parts[0]}\t{parts[1]}\n")
PY
fi

# ---------- optional: map<->gene coordinate join (portable replacement for coordenadas_gen.sh) ----------
if [[ "${DO_COORDENADAS_GEN:-0}" == "1" ]]; then
  [[ -n "${ORDER_MAPPED_TSV:-}" ]] || die "DO_COORDENADAS_GEN=1 but ORDER_MAPPED_TSV is empty"
  [[ -n "${RESULTADOS_LG_TSV:-}" ]] || die "DO_COORDENADAS_GEN=1 but RESULTADOS_LG_TSV is empty"
  stage_file "${ORDER_MAPPED_TSV}" "${WORKDIR}/order.mapped"
  stage_file "${RESULTADOS_LG_TSV}" "${WORKDIR}/resultadoslg.txt"
  out="${WORKDIR}/coordenadas_genes.txt"
  log "Joining order.mapped + resultadoslg.txt -> ${out}"
  python3 - <<'PY' "${WORKDIR}/order.mapped" "${WORKDIR}/resultadoslg.txt" "${out}"
import sys
order_path, res_path, out_path = sys.argv[1:]
# Build dict from resultadoslg (expects at least 4 cols, tab-delimited)
d = {}
with open(res_path, encoding="utf-8") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 1:
            d[parts[0]] = parts  # keep full row
with open(order_path, encoding="utf-8") as f, open(out_path, "w", encoding="utf-8") as o:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 4 and parts[0] in d:
            # echo -e "$col1\t$col3\t$col4\t${diccionario[$col1]}"
            o.write("\t".join([parts[0], parts[2], parts[3], *d[parts[0]]]) + "\n")
PY
fi

# ---------- optional: extract gene + position columns from coordinate files (portable replacement for extraer_pos_gen.sh) ----------
if [[ "${DO_EXTRAER_POS_GEN:-0}" == "1" ]]; then
  [[ -n "${COORDS_DIR:-}" ]] || die "DO_EXTRAER_POS_GEN=1 but COORDS_DIR is empty"
  [[ -d "${COORDS_DIR}" ]] || die "COORDS_DIR does not exist: ${COORDS_DIR}"
  outdir="${WORKDIR}/posGen"
  mkdir -p "${outdir}"
  log "Extracting columns 2,3,5 from *.txt in ${COORDS_DIR} -> ${outdir}"
  shopt -s nullglob
  for f in "${COORDS_DIR}"/*.txt; do
    bn="$(basename "$f" .txt)"
    awk -F'\t' '{print $2"\t"$3"\t"$5}' "$f" > "${outdir}/${bn}_posGen.txt"
  done
  shopt -u nullglob
fi

# ---------- optional scripts (run only if present) ----------
if [[ "${DO_OBTENER_LOCUS_GEN_POS:-0}" == "1" ]]; then
  if [[ -f "${SCRIPTS_DIR}/obtener_locus_gen_pos.py" ]]; then
    log "Running obtener_locus_gen_pos.py"
    ( cd "${WORKDIR}" && python3 "${SCRIPTS_DIR}/obtener_locus_gen_pos.py" ) |& tee "${WORKDIR}/logs/obtener_locus_gen_pos.log"
  else
    log "WARNING: obtener_locus_gen_pos.py not found in ${SCRIPTS_DIR} (skipping)"
  fi
fi

if [[ "${DO_POISSON:-0}" == "1" ]]; then
  if [[ -f "${SCRIPTS_DIR}/poisson_ultimo.R" ]]; then
    log "Running poisson_ultimo.R"
    ( cd "${WORKDIR}" && Rscript "${SCRIPTS_DIR}/poisson_ultimo.R" ) |& tee "${WORKDIR}/logs/poisson.log"
  else
    log "WARNING: poisson_ultimo.R not found in ${SCRIPTS_DIR} (skipping)"
  fi
fi

log "Done. Outputs are in: ${WORKDIR}"
