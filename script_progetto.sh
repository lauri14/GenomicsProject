#####################################
# Part I: getting the data
BCG2026_Cristianelli_E@159.149.160.7

# I will make a new folder to store the data
mkdir case594_summary

# I will then move to that folder
cd case594_summary
# to get the data I will use the ln -s command


ln -s /home/BCG2022_genomics_exam/uni.* .	
ln -s /home/BCG2022_genomics_exam/case594_* .	
ln -s /home/BCG2022_genomics_exam/targetsPad100.bed . 
ln -s /home/BCG2022_genomics_exam/worked/594/samples.txt



bowtie2 -U case594_child.fq.gz -x uni --rg-id "child" --rg "SM:child"  | samtools view -Sb | samtools sort -o child.bam 
bowtie2 -U case594_father.fq.gz -x uni --rg-id "father" --rg "SM:father"  | samtools view -Sb | samtools sort -o father.bam 
bowtie2 -U case594_mother.fq.gz -x uni --rg-id "mother" --rg "SM:mother"  | samtools view -Sb | samtools sort -o mother.bam 

mv case594_child.fq.gz child.fq.gz
mv case594_father.fq.gz father.fq.gz
mv case594_mother.fq.gz mother.fq.gz

samtools index child.bam
samtools index father.bam
samtools index mother.bam

fastqc *.bam 	
qualimap bamqc -bam child.bam --feature-file targetsPad100.bed -outdir child	
qualimap bamqc -bam father.bam --feature-file targetsPad100.bed -outdir father
qualimap bamqc -bam mother.bam --feature-file targetsPad100.bed -outdir mother


multiqc .

freebayes -f /home/BCG2022_genomics_exam/universe.fasta -m 20 -C 5 -Q 10 --min-coverage 10 mother.bam child.bam father.bam  > case594.vcf

bgzip case594.vcf
bcftools index case594.vcf.gz


bcftools view -R targetsPad100.bed case594.vcf.gz \ | bcftools view -S samples.txt \ | bcftools view -i 'GT[2]="0/1" && ((GT[0]="0/0" && GT[1]="0/0") || (GT[0]="0/1" || GT[1]="0/1"))' \ | bcftools filter -i 'QUAL>20' -Ov -o case594.cand.vcf

vep  -i case594.cand.vcf  -o case594.vep_annotated.vcf  --vcf  --cache   --offline   --assembly GRCh37  --dir_cache /data/vep_cache --use_given_ref --mane  --pick_allele --af --af_1kg  --af_gnomade --max_af --sift b --polyphen b

# Final filter: high impact + rare
filter_vep -i case594.vep_annotated.vcf   -o case594.vep_filtered.vcf   --filter "IMPACT is HIGH and (not MAX_AF or MAX_AF < 0.0001)"

#####################################
# Part VI coverage tracks
ad de novo :
bcftools view -R targetsPad100.bed case594.vcf.gz | bcftools view -S samples.txt | bcftools view -i 'GT[0]!="RR" && GT[1]="RR" && GT[2]="RR"' | bcftools filter -i 'QUAL>20' -Ov -o case594finale.cand.vcf

ref ref 
bcftools view -R targetsPad100.bed case594.vcf.gz | bcftools view -S samples.txt | bcftools view -i 'GT[0]="RR"  && GT[1]="RA" && GT[2]="RA"' | bcftools filter -i 'QUAL>20' -Ov -o case594.cand.vcf

# here I use -trackline and -trackopts to set the name of the tracks to be displayed at ucsc. It will not work otherwise
bedtools genomecov -ibam father.bam -bg -trackline -trackopts 'name="father"' -max 100 > fatherCov.bg
bedtools genomecov -ibam mother.bam -bg -trackline -trackopts 'name="mother"' -max 100 > motherCov.bg
bedtools genomecov -ibam child.bam -bg -trackline -trackopts 'name="child"' -max 100 > childCov.bg





#####
#!/bin/bash

#Make folder HW3
#mkdir hw3

#cd hw3

#IMPORT FILES
#Link all files that start with case... to my folder
ln -s /home/BCG2022_genomics_exam/uni.* .

#Link all files that start with uni. to my folder
ln -s /home/BCG2022_genomics_exam/universe* .

#"' .' tells you directory current - where this command gets run"
ln -s /home/BCG2022_genomics_exam/targetsPad100.bed .

list=(452 594)

for i in "${list[@]}"
do

ln -s /home/BCG2022_genomics_exam/case${i}* .

#CREATE OR EDIT TEXT FILE IN NANO
nano samplescase${i}.txt

# child
bowtie2 -U "case${i}_child.fq.gz"  -x uni --rg-id "childcase${i}" --rg "SM:childcase${i}" -p 12 | samtools view -Sb |samtools sort -o childcase${i}.bam

#father
bowtie2 -U "case${i}_father.fq.gz"  -x uni --rg-id "fathercase${i}" --rg "SM:fathercase${i}" -p 12 | samtools view -Sb |samtools sort -o fathercase${i}.bam

#mother
bowtie2 -U "case${i}_mother.fq.gz"  -x uni --rg-id "mothercase${i}" --rg "SM:mothercase${i}" -p 12 | samtools view -Sb |samtools sort -o mothercase${i}.bam

fastqc *case${i}.bam

#child
qualimap bamqc -bam "childcase${i}.bam" --feature-file targetsPad100.bed --outdir "childcase${i}"

#father
qualimap bamqc -bam "fathercase${i}.bam" --feature-file targetsPad100.bed --outdir "fathercase${i}"

#mother
qualimap bamqc -bam "mothercase${i}.bam" --feature-file targetsPad100.bed --outdir "mothercase${i}"

multiqc *case${i} .

mv multiqc_report.html "multiqc_report_case${i}.html"

mv multiqc_report_1.html "multiqc_report_case${i}.html"

freebayes -f universe.fasta -m 20 -C 5 -Q 10 -q 10 --min-coverage 10  "childcase${i}.bam" "fathercase${i}.bam" "mothercase${i}.bam" > "case${i}.vcf"

bedtools genomecov -ibam fathercase${i}.bam -bg -trackline -trackopts 'name="fathercase${i}"' -max 100 > fatherCovcase${i}.bg
bedtools genomecov -ibam mothercase${i}.bam -bg -trackline -trackopts 'name="mothercase${i}"' -max 100 > motherCovcase${i}.bg
bedtools genomecov -ibam childcase${i}.bam -bg -trackline -trackopts 'name="childcase${i}"' -max 100 > childCovcase${i}.bg

bgzip case${i}.vcf
    
bcftools index case${i}.vcf.gz

echo "before the if done"

if [ "${i}" -eq 452 ]; then
# Autosomal recessive

bcftools view -R targetsPad100.bed case${i}.vcf.gz | bcftools view -S samplescase452.txt | bcftools view -i 'GT[0]="AA" && GT[1]="RA" && GT[2]="RA"' | bcftools filter -i 'QUAL>20' -Ov -o case${i}.cand.vcf

echo " pipe OK"

elif [ "${i}" -eq 594 ]; then
# De Novo Autosomal Dominant
bcftools view -R targetsPad100.bed case${i}.vcf.gz | bcftools view -S samplescase594.txt | bcftools view -i 'GT[0]!="RR" && GT[1]="RR" && GT[2]="RR"' | bcftools filter -i 'QUAL>20' -Ov -o case${i}.cand.vcf
fi

done
