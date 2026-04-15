#!/bin/bash
set -Eeuo pipefail

########################################
# CONFIGURATION
########################################

USER_NAME="BCG2026_Cristianelli_E"
EXAM_DIR="/home/BCG2026_exam"
USER_DIR="${EXAM_DIR}/${USER_NAME}"
WORK_DIR="$HOME/project_Ennio"
THREADS=12

BED="${WORK_DIR}/chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed"
REF="${WORK_DIR}/chr20.fa"
SUMMARY_TSV="${WORK_DIR}/trio_summary.tsv"

########################################
# PREPARE WORKING DIRECTORY
########################################

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

ln -sf "${EXAM_DIR}/trios.txt" .
ln -sf "${EXAM_DIR}/list_disorders.txt" .
ln -sf "${EXAM_DIR}/chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed" .
ln -sf "${EXAM_DIR}/chr20.fa" .
ln -sf "${EXAM_DIR}/chr20.fa.fai" .
ln -sf "${EXAM_DIR}"/chr20*.bt2 .
ln -sf "${USER_DIR}/mode_inherithance.tsv" .

########################################
# INITIALIZE SUMMARY TABLE
########################################

printf "trio\tmode\ttotal_variants\tvariants_in_bed\tvariants_after_inheritance\tfinal_vep_filtered\n" > "$SUMMARY_TSV"

########################################
# PROCESS ALL TRIOS
########################################

for TRIO in trio_1 trio_2 trio_3 trio_4 trio_5
do
    echo "=========================================="
    echo "Processing ${TRIO}"
    echo "=========================================="

    mkdir -p "$TRIO"
    cd "$TRIO"

    MODE="NA"
    total_variants="NA"
    variants_in_bed="NA"
    variants_after_inheritance="NA"
    final_vep_filtered="NA"

    ########################################
    # LINK FASTQ FILES
    ########################################

    ln -sf "${USER_DIR}/${TRIO}"/*.fq.gz .

    ########################################
    # FIND SAMPLE IDS
    ########################################

    mapfile -t IDS < <(ls *.targets_R1.fq.gz | sed 's/\.targets_R1\.fq\.gz//' | sort)

    if [ "${#IDS[@]}" -ne 3 ]; then
        echo "ERROR: expected 3 samples in ${TRIO}, found ${#IDS[@]}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$TRIO" "$MODE" "$total_variants" "$variants_in_bed" "$variants_after_inheritance" "$final_vep_filtered" >> "$SUMMARY_TSV"
        exit 1
    fi

    trio_line=$(awk -v a="${IDS[0]}" -v b="${IDS[1]}" -v c="${IDS[2]}" '
        NR > 1 && (($1==a || $1==b || $1==c) && ($2==a || $2==b || $2==c) && ($3==a || $3==b || $3==c)) {
            print
            exit
        }' ../trios.txt)

    if [ -z "$trio_line" ]; then
        echo "ERROR: could not find matching trio in trios.txt for ${TRIO}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$TRIO" "$MODE" "$total_variants" "$variants_in_bed" "$variants_after_inheritance" "$final_vep_filtered" >> "$SUMMARY_TSV"
        exit 1
    fi

    # trios.txt columns: child father mother
    child_id=$(echo "$trio_line" | awk '{print $1}')
    father_id=$(echo "$trio_line" | awk '{print $2}')
    mother_id=$(echo "$trio_line" | awk '{print $3}')

    echo "Child : $child_id"
    echo "Father: $father_id"
    echo "Mother: $mother_id"

    ########################################
    # MAPPING
    ########################################

    bowtie2 -x ../chr20 \
      -1 "${child_id}.targets_R1.fq.gz" \
      -2 "${child_id}.targets_R2.fq.gz" \
      --rg-id child --rg SM:child -p "$THREADS" \
      | samtools view -Sb - \
      | samtools sort -o child.bam

    bowtie2 -x ../chr20 \
      -1 "${father_id}.targets_R1.fq.gz" \
      -2 "${father_id}.targets_R2.fq.gz" \
      --rg-id father --rg SM:father -p "$THREADS" \
      | samtools view -Sb - \
      | samtools sort -o father.bam

    bowtie2 -x ../chr20 \
      -1 "${mother_id}.targets_R1.fq.gz" \
      -2 "${mother_id}.targets_R2.fq.gz" \
      --rg-id mother --rg SM:mother -p "$THREADS" \
      | samtools view -Sb - \
      | samtools sort -o mother.bam

    samtools index child.bam
    samtools index father.bam
    samtools index mother.bam

    ########################################
    # QC
    ########################################

    fastqc *.bam

    qualimap bamqc -bam child.bam  --feature-file "$BED" --outdir child
    qualimap bamqc -bam father.bam --feature-file "$BED" --outdir father
    qualimap bamqc -bam mother.bam --feature-file "$BED" --outdir mother

    multiqc . -f -o . -n "multiqc_${TRIO}.html"

    ########################################
    # VARIANT CALLING
    ########################################

    freebayes -f "$REF" -t "$BED" child.bam father.bam mother.bam \
        -m 20 -C 5 -Q 10 --min-coverage 10 \
        | bcftools sort -Ov -o "${TRIO}.vcf"

    bgzip -f "${TRIO}.vcf"
    bcftools index -f "${TRIO}.vcf.gz"

    if [ ! -s "${TRIO}.vcf.gz" ]; then
        echo "ERROR: ${TRIO}.vcf.gz was not created"
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$TRIO" "$MODE" "$total_variants" "$variants_in_bed" "$variants_after_inheritance" "$final_vep_filtered" >> "$SUMMARY_TSV"
        exit 1
    fi

    ########################################
    # SAMPLE ORDER FOR FILTERING
    ########################################

    printf "child\nfather\nmother\n" > samples.txt

    echo "Sample order after bcftools -S:"
    reordered_samples=$(mktemp)
    bcftools view -S samples.txt "${TRIO}.vcf.gz" -Ov -o "$reordered_samples"
    bcftools query -l "$reordered_samples"
    rm -f "$reordered_samples"

    ########################################
    # COUNTS BEFORE FILTERS
    ########################################

    total_variants=$(bcftools view -H "${TRIO}.vcf.gz" | wc -l)
    variants_in_bed=$(bcftools view -H -R "$BED" "${TRIO}.vcf.gz" | wc -l)

    ########################################
    # INHERITANCE FILTERING
    ########################################

    if [[ "$TRIO" == "trio_1" || "$TRIO" == "trio_4" || "$TRIO" == "trio_5" ]]; then
        MODE="autosomal recessive"
        GT_FILTER='GT[0]="1/1" && GT[1]!="0/0" && GT[2]!="0/0"'

    elif [[ "$TRIO" == "trio_2" ]]; then
        MODE="autosomal dominant inherited"
        GT_FILTER='GT[0]!="0/0" && GT[1]!="0/0" && GT[2]="0/0"'

    elif [[ "$TRIO" == "trio_3" ]]; then
        MODE="autosomal dominant de novo"
        GT_FILTER='GT[0]!="0/0" && GT[1]="0/0" && GT[2]="0/0"'

    else
        echo "ERROR: unknown trio ${TRIO}"
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$TRIO" "$MODE" "$total_variants" "$variants_in_bed" "$variants_after_inheritance" "$final_vep_filtered" >> "$SUMMARY_TSV"
        exit 1
    fi

    echo "Mode: ${MODE}"
    echo "GT filter: ${GT_FILTER}"

    variants_after_inheritance=$(
        bcftools view -R "$BED" "${TRIO}.vcf.gz" \
          | bcftools view -S samples.txt \
          | bcftools view -i "$GT_FILTER" -H \
          | wc -l
    )

    bcftools view -R "$BED" "${TRIO}.vcf.gz" \
        | bcftools view -S samples.txt \
        | bcftools view -i "$GT_FILTER" -Ov \
        > "${TRIO}.cand.vcf"

    echo "Candidate count for ${TRIO}: $(grep -vc '^#' "${TRIO}.cand.vcf" || true)"

    if [ ! -s "${TRIO}.cand.vcf" ]; then
        echo "WARNING: ${TRIO}.cand.vcf is empty after inheritance filtering"
        final_vep_filtered=0
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$TRIO" "$MODE" "$total_variants" "$variants_in_bed" "$variants_after_inheritance" "$final_vep_filtered" >> "$SUMMARY_TSV"
        cd ..
        continue
    fi

    if [ "$(grep -vc '^#' "${TRIO}.cand.vcf" || true)" -eq 0 ]; then
        echo "WARNING: ${TRIO}.cand.vcf contains no candidate variants"
        final_vep_filtered=0
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$TRIO" "$MODE" "$total_variants" "$variants_in_bed" "$variants_after_inheritance" "$final_vep_filtered" >> "$SUMMARY_TSV"
        cd ..
        continue
    fi

    ########################################
    # ANNOTATION
    ########################################

    vep -i "${TRIO}.cand.vcf" \
      -o "${TRIO}.vep_annotated.vcf" \
      --format vcf \
      --vcf \
      --cache \
      --offline \
      --dir_cache /data/vep_cache \
      --species homo_sapiens \
      --assembly GRCh38 \
      --refseq \
      --cache_version 115 \
      --use_given_ref \
      --mane \
      --pick_allele \
      --af \
      --af_1kg \
      --max_af \
      --sift b \
      --polyphen b

    if [ ! -s "${TRIO}.vep_annotated.vcf" ]; then
        echo "WARNING: VEP did not produce ${TRIO}.vep_annotated.vcf"
        final_vep_filtered="NA"
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$TRIO" "$MODE" "$total_variants" "$variants_in_bed" "$variants_after_inheritance" "$final_vep_filtered" >> "$SUMMARY_TSV"
        cd ..
        continue
    fi

    filter_vep -i "${TRIO}.vep_annotated.vcf" \
      -o "${TRIO}.vep_filtered.vcf" \
      --format vcf \
      --filter "(IMPACT is HIGH or IMPACT is MODERATE) and (not MAX_AF or MAX_AF < 0.0001)"

    if [ ! -s "${TRIO}.vep_filtered.vcf" ]; then
        echo "WARNING: filter_vep did not produce ${TRIO}.vep_filtered.vcf"
        final_vep_filtered="NA"
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$TRIO" "$MODE" "$total_variants" "$variants_in_bed" "$variants_after_inheritance" "$final_vep_filtered" >> "$SUMMARY_TSV"
        cd ..
        continue
    fi

    final_vep_filtered=$(grep -vc '^#' "${TRIO}.vep_filtered.vcf" || true)

    ########################################
    # COVERAGE TRACKS
    ########################################

    bedtools genomecov -ibam child.bam  -bg -trackline -trackopts 'name="child"'  -max 100 > childCov.bg
    bedtools genomecov -ibam father.bam -bg -trackline -trackopts 'name="father"' -max 100 > fatherCov.bg
    bedtools genomecov -ibam mother.bam -bg -trackline -trackopts 'name="mother"' -max 100 > motherCov.bg

    ########################################
    # WRITE SUMMARY ROW
    ########################################

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$TRIO" "$MODE" "$total_variants" "$variants_in_bed" "$variants_after_inheritance" "$final_vep_filtered" >> "$SUMMARY_TSV"

    cd ..
done

echo "All trios processed."
echo "Summary written to: $SUMMARY_TSV"