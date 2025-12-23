#!/usr/bin/env bash
set -euo pipefail

FASTQ_DIR="${FASTQ_DIR:-/path/to/fastqs}"
SAMPLE_PREFIX="${SAMPLE_PREFIX:-sample_}"
FIRST_SAMPLE="${FIRST_SAMPLE:-1}"
LAST_SAMPLE="${LAST_SAMPLE:-22}"

# FASTQ filename templates for integer sample IDs.
# Default matches Illumina bcl2fastq style: 1_S1_L001_R1_001.fastq.gz
R1_TEMPLATE="${R1_TEMPLATE:-{i}_S{i}_L001_R1_001.fastq.gz}"
R2_TEMPLATE="${R2_TEMPLATE:-{i}_S{i}_L001_R2_001.fastq.gz}"

# Optional: provide a sample list file (one integer per line) instead of FIRST_SAMPLE..LAST_SAMPLE
SAMPLE_LIST_FILE="${SAMPLE_LIST_FILE:-}"

# Taxon list files (one taxon per line, comments allowed with #)
HOST_TAXA_FILE="${HOST_TAXA_FILE:-configs/hosts/host_taxa.txt}"
PATHOGEN_TAXA_FILE="${PATHOGEN_TAXA_FILE:-configs/panels/pathogen_taxa.txt}"
VIRUS_TAXA_FILE="${VIRUS_TAXA_FILE:-configs/panels/virus_taxa.txt}"

THREADS="${THREADS:-14}"
OUTDIR="${OUTDIR:-confirm_out_competitive_v12}"

MIN_MAPQ="${MIN_MAPQ:-20}"

STRONG_MIN_MAPQ20="${STRONG_MIN_MAPQ20:-500}"
STRONG_MIN_BREADTH_PCT="${STRONG_MIN_BREADTH_PCT:-0.10}"
STRONG_MIN_MEAN_DEPTH="${STRONG_MIN_MEAN_DEPTH:-0.005}"

WEAK_MIN_MAPQ20="${WEAK_MIN_MAPQ20:-50}"
WEAK_MIN_BREADTH_PCT="${WEAK_MIN_BREADTH_PCT:-0.01}"

MASK_REF="${MASK_REF:-1}"
FILTER_LOW_COMPLEX_READS="${FILTER_LOW_COMPLEX_READS:-1}"
BBDUK_ENTROPY="${BBDUK_ENTROPY:-0.70}"
BBDUK_WINDOW="${BBDUK_WINDOW:-50}"
BBDUK_K="${BBDUK_K:-5}"

MM2_SPLIT_DIR="${MM2_SPLIT_DIR:-${OUTDIR}/tmp}"

MM2_HOST_OPTS="${MM2_HOST_OPTS:- -x sr -I 1G}"
MM2_PATH_OPTS="${MM2_PATH_OPTS:- -x sr -I 1G}"


# Load taxa lists from files (skip blank lines and comments)
load_taxa_file() {
  local f="$1"
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: taxa file not found: ${f}" >&2
    exit 10
  fi
  awk '{gsub(/\r/,""); if($0 ~ /^[[:space:]]*#/){next} if($0 ~ /^[[:space:]]*$/){next} print $0}' "${f}"
}

mapfile -t HOST_TAXA < <(load_taxa_file "${HOST_TAXA_FILE}")
mapfile -t PATHOGEN_TAXA < <(load_taxa_file "${PATHOGEN_TAXA_FILE}")

# Virus taxa file is optional; if absent or empty, skip virus downloading
VIRUS_TAXA=()
if [[ -f "${VIRUS_TAXA_FILE}" ]]; then
  mapfile -t VIRUS_TAXA < <(load_taxa_file "${VIRUS_TAXA_FILE}" || true)
fi

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing tool: $1" >&2; exit 1; }; }
need minimap2; need samtools; need pigz; need datasets; need dataformat; need unzip; need awk; need sed; need python3
if [[ "${MASK_REF}" == "1" ]]; then need dustmasker; fi
if [[ "${FILTER_LOW_COMPLEX_READS}" == "1" ]]; then need bbduk.sh; fi

mkdir -p "${OUTDIR}"/{logs,tmp,refs,host_refs,pathogen_refs,host_combined,pathogen_combined,nonhost_fastq,bams,metrics}
mkdir -p "${MM2_SPLIT_DIR}"

log() { echo "[$(date +'%F %T')] $*" | tee -a "${OUTDIR}/logs/run.log" >&2; }

get_assembly_acc() {
  local taxon="$1" level="$2" source="$3"
  datasets summary genome taxon "${taxon}" --assembly-level "${level}" --assembly-source "${source}" \
    --limit 1 --as-json-lines \
  | dataformat tsv genome --fields accession \
  | awk 'NR==2{print $1}'
}

pick_best_assembly_host() {
  local taxon="$1" acc=""
  for combo in "chromosome RefSeq" "chromosome all" "complete RefSeq" "scaffold RefSeq" "contig RefSeq"; do
    local level="${combo% *}" source="${combo#* }"
    acc="$(get_assembly_acc "${taxon}" "${level}" "${source}" || true)"
    [[ -n "${acc}" && "${acc}" != "null" ]] && { echo "${acc}"; return 0; }
  done
  echo ""
}

pick_best_assembly_pathogen() {
  local taxon="$1" acc=""
  for combo in "complete RefSeq" "complete GenBank" "complete all" "chromosome RefSeq" "chromosome all"; do
    local level="${combo% *}" source="${combo#* }"
    acc="$(get_assembly_acc "${taxon}" "${level}" "${source}" || true)"
    [[ -n "${acc}" && "${acc}" != "null" ]] && { echo "${acc}"; return 0; }
  done
  echo ""
}

download_genome_fna_by_assembly() {
  local taxon="$1" acc="$2" out_fna="$3" manifest="$4"
  if [[ -s "${out_fna}" ]]; then
    log "REF exists, skip: ${out_fna}"
    echo -e "${taxon}\t${acc}\t${out_fna}\tOK(SKIP)" >> "${manifest}"
    return 0
  fi
  local tmpzip="${OUTDIR}/refs/${acc}.zip"
  local tmpdir="${OUTDIR}/refs/${acc}"
  rm -rf "${tmpdir}"; mkdir -p "${tmpdir}"
  log "Downloading genome for ${taxon} (assembly ${acc})"
  datasets download genome accession "${acc}" --include genome --filename "${tmpzip}" >/dev/null
  unzip -q "${tmpzip}" -d "${tmpdir}"
  local fna
  fna="$(find "${tmpdir}/ncbi_dataset/data" -name "*_genomic.fna" -o -name "genomic.fna" | head -n 1)"
  if [[ -z "${fna}" ]]; then
    log "WARN: cannot find genomic.fna for ${taxon} (${acc})"
    echo -e "${taxon}\t${acc}\tNA\tFAILED" >> "${manifest}"
    return 1
  fi
  cp "${fna}" "${out_fna}"
  echo -e "${taxon}\t${acc}\t${out_fna}\tOK" >> "${manifest}"
  log "Saved: ${out_fna}"
}

download_host_refs() {
  local outdir="${OUTDIR}/host_refs/individual"; mkdir -p "${outdir}"
  local manifest="${OUTDIR}/host_refs/host_manifest.tsv"
  echo -e "taxon\tassembly_accession\tfasta_path\tstatus" > "${manifest}"
  log "Downloading host references"
  for taxon in "${HOST_TAXA[@]}"; do
    local out="${outdir}/${taxon// /_}.fna"
    if [[ -s "${out}" ]]; then
      log "REF exists, skip: ${out}"
      echo -e "${taxon}\tNA\t${out}\tOK(SKIP)" >> "${manifest}"
      continue
    fi
    local acc; acc="$(pick_best_assembly_host "${taxon}")"
    if [[ -z "${acc}" ]]; then
      log "WARN: no host assembly found for ${taxon}"
      echo -e "${taxon}\tNA\tNA\tFAILED" >> "${manifest}"
      continue
    fi
    download_genome_fna_by_assembly "${taxon}" "${acc}" "${out}" "${manifest}" || true
  done
}

download_pathogen_refs() {
  local outdir="${OUTDIR}/pathogen_refs/individual"; mkdir -p "${outdir}"
  local manifest="${OUTDIR}/pathogen_refs/pathogen_manifest.tsv"
  echo -e "taxon\tassembly_accession\tfasta_path\tstatus" > "${manifest}"
  log "Downloading curated pathogen references"
  for taxon in "${PATHOGEN_TAXA[@]}"; do
    local out="${outdir}/${taxon// /_}.fna"
    if [[ -s "${out}" ]]; then
      log "REF exists, skip: ${out}"
      echo -e "${taxon}\tNA\t${out}\tOK(SKIP)" >> "${manifest}"
      continue
    fi
    local acc; acc="$(pick_best_assembly_pathogen "${taxon}")"
    if [[ -z "${acc}" ]]; then
      log "WARN: no assembly found for ${taxon}, skipping"
      echo -e "${taxon}\tNA\tNA\tFAILED" >> "${manifest}"
      continue
    fi
    download_genome_fna_by_assembly "${taxon}" "${acc}" "${out}" "${manifest}" || true
  done

  log "Downloading curated virus references"
  for taxon in "${VIRUS_TAXA[@]}"; do
    local out="${outdir}/${taxon// /_}.fna"
    if [[ -s "${out}" ]]; then
      log "REF exists, skip: ${out}"
      echo -e "${taxon}\tNA\t${out}\tOK(SKIP)" >> "${manifest}"
      continue
    fi
    local tmpzip="${OUTDIR}/refs/virus_${taxon// /_}.zip"
    local tmpdir="${OUTDIR}/refs/virus_${taxon// /_}"
    rm -rf "${tmpdir}"; mkdir -p "${tmpdir}"
    log "Downloading virus genome for ${taxon}"
    datasets download virus genome taxon "${taxon}" --complete-only --include genome --filename "${tmpzip}" >/dev/null || true
    unzip -q "${tmpzip}" -d "${tmpdir}" || true
    local fna
    fna="$(find "${tmpdir}/ncbi_dataset/data" -name "genomic.fna" -o -name "*genomic*.fna" | head -n 1 || true)"
    if [[ -z "${fna}" ]]; then
      log "WARN: cannot find virus genomic.fna for ${taxon}, skipping"
      echo -e "${taxon}\tNA\tNA\tFAILED" >> "${manifest}"
      continue
    fi
    cp "${fna}" "${out}"
    echo -e "${taxon}\tNA\t${out}\tOK" >> "${manifest}"
    log "Saved: ${out}"
  done
}

combine_and_prefix() {
  local input_glob="$1" out_fna="$2"
  : > "${out_fna}"
  mapfile -t files < <(ls -1 ${input_glob} 2>/dev/null | sort || true)
  [[ ${#files[@]} -eq 0 ]] && { log "ERROR: no FASTAs matched: ${input_glob}"; exit 1; }
  for ref in "${files[@]}"; do
    local ref_name; ref_name="$(basename "${ref}" .fna)"
    awk -v P="${ref_name}" '
      BEGIN{n=0}
      /^>/{
        n++;
        printf(">%s|seq%06d\n", P, n);
        next
      }
      {print}
    ' "${ref}" >> "${out_fna}"
  done
}

prepare_refs() {
  local host_comb="${OUTDIR}/host_combined/combined_host.fna"
  local path_comb="${OUTDIR}/pathogen_combined/combined_pathogens.fna"
  local path_use="${OUTDIR}/pathogen_combined/combined_pathogens.use.fna"

  [[ -s "${host_comb}" ]] || { log "Combining host FASTAs"; combine_and_prefix "${OUTDIR}/host_refs/individual/*.fna" "${host_comb}"; }
  [[ -s "${path_comb}" ]] || { log "Combining pathogen FASTAs"; combine_and_prefix "${OUTDIR}/pathogen_refs/individual/*.fna" "${path_comb}"; }

  if [[ "${MASK_REF}" == "1" ]]; then
    if [[ ! -s "${path_use}" ]]; then
      log "Masking combined pathogen reference with dustmasker"
      dustmasker -in "${path_comb}" -out "${path_use}" -outfmt fasta > "${OUTDIR}/logs/dustmasker.log" 2>&1 || true
      if [[ ! -s "${path_use}" ]]; then
        log "WARN dustmasker failed, using unmasked reference"
        cp "${path_comb}" "${path_use}"
      fi
    fi
  else
    [[ -s "${path_use}" ]] || cp "${path_comb}" "${path_use}"
  fi
}

get_fastqs_for_i() {
  local i="$1"
  local r1="${FASTQ_DIR}/$(python3 - <<PY
i=int("${i}")
print("${R1_TEMPLATE}".format(i=i))
PY
)"
  local r2="${FASTQ_DIR}/$(python3 - <<PY
i=int("${i}")
print("${R2_TEMPLATE}".format(i=i))
PY
)"
  [[ -f "${r1}" && -f "${r2}" ]] && echo -e "${r1}	${r2}" && return 0
  return 1
}

host_subtract() {
  local sample="$1" r1="$2" r2="$3"
  local host_ref="${OUTDIR}/host_combined/combined_host.fna"
  local sdir="${OUTDIR}/bams/${sample}"; mkdir -p "${sdir}"
  local hostbam="${sdir}/host.unsorted.bam"
  local nh_r1="${OUTDIR}/nonhost_fastq/${sample}.nonhost_R1.fastq.gz"
  local nh_r2="${OUTDIR}/nonhost_fastq/${sample}.nonhost_R2.fastq.gz"
  local mm2log="${OUTDIR}/logs/${sample}.host_mm2.log"
  local split="${MM2_SPLIT_DIR}/mm2split.${sample}.host"

  if [[ -s "${nh_r1}" && -s "${nh_r2}" && -s "${hostbam}" ]]; then
    echo -e "${nh_r1}\t${nh_r2}"
    return 0
  fi

  log "Host mapping ${sample}"
  minimap2 -a ${MM2_HOST_OPTS} -t "${THREADS}" --secondary=no --split-prefix="${split}" \
    "${host_ref}" "${r1}" "${r2}" 2> "${mm2log}" \
    | samtools view -b -o "${hostbam}" -

  if ! samtools view -H "${hostbam}" >/dev/null 2>&1; then
    log "ERROR host mapping failed for ${sample}"
    tail -n 80 "${mm2log}" >&2 || true
    exit 2
  fi

  samtools flagstat "${hostbam}" > "${sdir}/host.flagstat.txt" || true

  log "Extract non-host read pairs (both mates unmapped) ${sample}"
  samtools fastq -f 12 -F 256 -n \
    -1 >(pigz -p "${THREADS}" > "${nh_r1}") \
    -2 >(pigz -p "${THREADS}" > "${nh_r2}") \
    -0 /dev/null -s /dev/null \
    "${hostbam}" >/dev/null

  echo -e "${nh_r1}\t${nh_r2}"
}

filter_low_complex_reads() {
  local sample="$1" in1="$2" in2="$3"
  if [[ "${FILTER_LOW_COMPLEX_READS}" != "1" ]]; then
    echo -e "${in1}\t${in2}"
    return 0
  fi
  local out1="${OUTDIR}/nonhost_fastq/${sample}.nonhost.clean_R1.fastq.gz"
  local out2="${OUTDIR}/nonhost_fastq/${sample}.nonhost.clean_R2.fastq.gz"
  if [[ -s "${out1}" && -s "${out2}" ]]; then
    echo -e "${out1}\t${out2}"
    return 0
  fi
  log "Low-complexity filter (bbduk) ${sample}"
  bbduk.sh in1="${in1}" in2="${in2}" out1="${out1}" out2="${out2}" \
    entropy="${BBDUK_ENTROPY}" entropywindow="${BBDUK_WINDOW}" entropyk="${BBDUK_K}" \
    stats="${OUTDIR}/metrics/${sample}.bbduk_stats.txt" overwrite=t \
    2> "${OUTDIR}/logs/${sample}.bbduk.log" || true
  if [[ ! -s "${out1}" || ! -s "${out2}" ]]; then
    log "WARN bbduk outputs missing for ${sample}, using unfiltered reads"
    echo -e "${in1}\t${in2}"
    return 0
  fi
  echo -e "${out1}\t${out2}"
}

map_to_pathogens() {
  local sample="$1" r1="$2" r2="$3"
  local ref="${OUTDIR}/pathogen_combined/combined_pathogens.use.fna"
  local sdir="${OUTDIR}/bams/${sample}"; mkdir -p "${sdir}"
  local bam="${sdir}/combined.bam"
  local mm2log="${OUTDIR}/logs/${sample}.path_mm2.log"
  local split="${MM2_SPLIT_DIR}/mm2split.${sample}.path"

  if [[ -s "${bam}" ]]; then
    echo "${bam}"
    return 0
  fi

  log "Competitive mapping ${sample}"
  minimap2 -a ${MM2_PATH_OPTS} -t "${THREADS}" --secondary=no --split-prefix="${split}" \
    "${ref}" "${r1}" "${r2}" 2> "${mm2log}" \
    | samtools sort -@ "${THREADS}" -o "${bam}" -

  if ! samtools view -H "${bam}" >/dev/null 2>&1; then
    log "ERROR pathogen mapping failed for ${sample}"
    tail -n 120 "${mm2log}" >&2 || true
    exit 3
  fi
  samtools index "${bam}"
  samtools flagstat "${bam}" > "${sdir}/combined.flagstat.txt" || true
  echo "${bam}"
}

metrics_for_sample() {
  local sample="$1" bam="$2"
  local mdir="${OUTDIR}/metrics/${sample}"; mkdir -p "${mdir}"

  log "Metrics ${sample}"

  samtools view -F 0x904 "${bam}" | cut -f 3 | sed 's/|seq[0-9]\+$//' | sort | uniq -c \
    | awk '{print $2"\t"$1}' \
    > "${mdir}/counts_primary.tsv"

  samtools view -F 0x904 -q "${MIN_MAPQ}" "${bam}" | cut -f 3 | sed 's/|seq[0-9]\+$//' | sort | uniq -c \
    | awk '{print $2"\t"$1}' \
    > "${mdir}/counts_primary_mapq${MIN_MAPQ}.tsv"

  samtools coverage -q "${MIN_MAPQ}" --ff "UNMAP,SECONDARY,QCFAIL,SUPPLEMENTARY" "${bam}" \
    > "${mdir}/coverage_by_contig.tsv" || true

  awk -v S="${sample}" '
    BEGIN{OFS="\t"}
    NR==1{next}
    {
      split($1,a,"|"); ref=a[1];
      len=($3-$2+1);
      covbases=$5;
      meandepth=$7;
      tot_len[ref]+=len;
      tot_covbases[ref]+=covbases;
      tot_depth_x_len[ref]+=meandepth*len;
    }
    END{
      for(ref in tot_len){
        breadth = (tot_len[ref]>0 ? (tot_covbases[ref]/tot_len[ref])*100.0 : 0.0);
        mean_depth = (tot_len[ref]>0 ? (tot_depth_x_len[ref]/tot_len[ref]) : 0.0);
        printf "%s\t%s\t%.6f\t%.6f\n", S, ref, mean_depth, breadth
      }
    }
  ' "${mdir}/coverage_by_contig.tsv" | sort -k2,2 > "${mdir}/depth_breadth_by_ref.tsv"

  {
    cut -f1 "${mdir}/counts_primary.tsv" 2>/dev/null || true
    cut -f1 "${mdir}/counts_primary_mapq${MIN_MAPQ}.tsv" 2>/dev/null || true
    cut -f2 "${mdir}/depth_breadth_by_ref.tsv" 2>/dev/null || true
  } | sort -u > "${mdir}/ref_names.all.txt"

  sort -k1,1 "${mdir}/counts_primary.tsv" > "${mdir}/_c1.tsv" 2>/dev/null || true
  sort -k1,1 "${mdir}/counts_primary_mapq${MIN_MAPQ}.tsv" > "${mdir}/_c2.tsv" 2>/dev/null || true

  local outrows="${mdir}/summary_rows.tsv"
  : > "${outrows}"
  while read -r ref; do
    mp=$(awk -v r="${ref}" '$1==r{print $2}' "${mdir}/_c1.tsv" 2>/dev/null || true); mp=${mp:-0}; mp=$((mp+0))
    mq=$(awk -v r="${ref}" '$1==r{print $2}' "${mdir}/_c2.tsv" 2>/dev/null || true); mq=${mq:-0}; mq=$((mq+0))
    db=$(awk -v r="${ref}" '$2==r{print $3"\t"$4}' "${mdir}/depth_breadth_by_ref.tsv" 2>/dev/null || true)
    md=$(echo "${db}" | awk '{print $1}'); md=${md:-0}
    br=$(echo "${db}" | awk '{print $2}'); br=${br:-0}

    label=$(awk -v mq="${mq}" -v br="${br}" -v md="${md}" \
      -v smq="${STRONG_MIN_MAPQ20}" -v sbr="${STRONG_MIN_BREADTH_PCT}" -v smd="${STRONG_MIN_MEAN_DEPTH}" \
      -v wmq="${WEAK_MIN_MAPQ20}" -v wbr="${WEAK_MIN_BREADTH_PCT}" '
      BEGIN{
        if(mq+0 >= smq+0 && br+0 >= sbr+0 && md+0 >= smd+0) {print "STRONG_SIGNAL"; exit}
        if(mq+0 >= wmq+0 && br+0 >= wbr+0) {print "WEAK_SIGNAL"; exit}
        print "NO_SIGNAL"
      }')

    printf "%s\t%s\t%d\t%d\t%.6f\t%.6f\t%s\n" "${sample}" "${ref}" "${mp}" "${mq}" "${md}" "${br}" "${label}" >> "${outrows}"
  done < "${mdir}/ref_names.all.txt"

  sort -k4,4nr "${outrows}" -o "${outrows}"
}

main() {
  log "START host-subtraction + competitive confirmatory mapping v12"
  log "FASTQ_DIR=${FASTQ_DIR} SAMPLES=${FIRST_SAMPLE}..${LAST_SAMPLE}"
  log "OUTDIR=${OUTDIR} THREADS=${THREADS} MIN_MAPQ=${MIN_MAPQ}"
  log "MASK_REF=${MASK_REF} FILTER_LOW_COMPLEX_READS=${FILTER_LOW_COMPLEX_READS}"
  log "STRONG thresholds: mapq20>=${STRONG_MIN_MAPQ20} breadth>=${STRONG_MIN_BREADTH_PCT}% mean_depth>=${STRONG_MIN_MEAN_DEPTH}"
  log "WEAK thresholds: mapq20>=${WEAK_MIN_MAPQ20} breadth>=${WEAK_MIN_BREADTH_PCT}%"

  download_host_refs
  download_pathogen_refs
  prepare_refs

  local summary="${OUTDIR}/confirm_summary_competitive.tsv"
  echo -e "sample\tref_name\tmapped_primary\tmapped_primary_mapq${MIN_MAPQ}\tmean_depth\tbreadth_covered_pct\tsignal_label" > "${summary}"

  if [[ -n "${SAMPLE_LIST_FILE}" ]]; then
    mapfile -t _SAMPLES < "${SAMPLE_LIST_FILE}"
  else
    mapfile -t _SAMPLES < <(seq "${FIRST_SAMPLE}" "${LAST_SAMPLE}")
  fi
  for i in "${_SAMPLES[@]}"; do
    fq="$(get_fastqs_for_i "${i}")" || { log "WARN missing FASTQs for sample ${i}, skip"; continue; }
    r1="$(echo "${fq}" | awk '{print $1}')"
    r2="$(echo "${fq}" | awk '{print $2}')"
    sample="${SAMPLE_PREFIX}${i}"

    nh="$(host_subtract "${sample}" "${r1}" "${r2}")"
    nh1="$(echo "${nh}" | awk '{print $1}')"
    nh2="$(echo "${nh}" | awk '{print $2}')"
    if [[ ! -s "${nh1}" || ! -s "${nh2}" ]]; then
      log "WARN nonhost reads empty for ${sample}, skip"
      continue
    fi

    clean="$(filter_low_complex_reads "${sample}" "${nh1}" "${nh2}")"
    c1="$(echo "${clean}" | awk '{print $1}')"
    c2="$(echo "${clean}" | awk '{print $2}')"

    bam="$(map_to_pathogens "${sample}" "${c1}" "${c2}")"
    metrics_for_sample "${sample}" "${bam}"

    cat "${OUTDIR}/metrics/${sample}/summary_rows.tsv" >> "${summary}"
  done

  log "DONE: ${summary}"
}

main "$@"
