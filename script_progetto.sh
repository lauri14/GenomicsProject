BCG2026_Cristianelli_E@159.149.160.7

#i rename the files basing on their role
#files from Cristianelli Ennio
mv HG00451.targets_R1.fq.gz child1.fq.gz
mv HG00451.targets_R2.fq.gz child2.fq.gz
mv HG00452.targets_R1.fq.gz father1.fq.gz
mv HG00452.targets_R2.fq.gz father2.fq.gz
mv HG00453.targets_R1.fq.gz mother1.fq.gz
mv HG00453.targets_R2.fq.gz mother2.fq.gz
#files from Andaloro Laura
mv HG00448.targets_R1.fq.gz child1.fq.gz
mv HG00448.targets_R2.fq.gz child2.fq.gz
mv HG00449.targets_R1.fq.gz father1.fq.gz
mv HG00449.targets_R2.fq.gz father2.fq.gz
mv HG00450.targets_R1.fq.gz mother1.fq.gz
mv HG00450.targets_R2.fq.gz mother2.fq.gz


bowtie2 -1 child1.fq.gz -2 child2.fq.gz  -x ../chr20 --rg-id "child" --rg "SM:child"  | samtools view -Sb | samtools sort -o child.bam 
bowtie2 -1 father1.fq.gz -2 father2.fq.gz -x ../chr20 --rg-id "father" --rg "SM:father"  | samtools view -Sb | samtools sort -o father.bam 
bowtie2 -1 mother1.fq.gz -2 mother2.fq.gz -x ../chr20 --rg-id "mother" --rg "SM:mother"  | samtools view -Sb | samtools sort -o mother.bam 

samtools index child.bam
samtools index father.bam
samtools index mother.bam

fastqc *.bam 	
qualimap bamqc -bam child.bam --feature-file ../chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed -outdir child	
qualimap bamqc -bam father.bam --feature-file ../chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed -outdir father
qualimap bamqc -bam mother.bam --feature-file ../chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed -outdir mother


multiqc .


num=$(basename "$PWD" | grep -o '[0-9]\+')
freebayes -f ../chr20.fa -m 20 -C 5 -Q 10 --min-coverage 5 mother.bam child.bam father.bam > trio${num}.vcf

bgzip trio${num}.vcf
bcftools index trio${num}.vcf.gz


# single line version of the if block
if [ "$num" -eq 1 ] || [ "$num" -eq 4 ] || [ "$num" -eq 5 ]; then bcftools view -R ../chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed trio${num}.vcf.gz | bcftools view -S ../samples.txt | bcftools view -i 'GT[0]="AA" && GT[1]="RA" && GT[2]="RA"' | bcftools filter -i 'QUAL>20' -Ov -o trio${num}.cand.vcf; elif [ "$num" -eq 3 ]; then bcftools view -R ../chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed trio${num}.vcf.gz | bcftools view -S ../samples.txt | bcftools view -i 'GT[0]!="RR" && GT[1]="RR" && GT[2]="RR"' | bcftools filter -i 'QUAL>20' -Ov -o trio${num}.cand.vcf; elif [ "$num" -eq 2 ]; then bcftools view -R ../chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed trio${num}.vcf.gz | bcftools view -S ../samples.txt | bcftools view -i 'GT[0] != "RR" && (GT[1] != "RR" || GT[2] != "RR")' | bcftools filter -i 'QUAL>20' -Ov -o trio${num}.cand.vcf; fi




##########
#if block for inherithance conditions


if [ "${num}" -eq 1 ]||[ "${num}" -eq 4 ]|ls|[ "${num}" -eq 5 ]; then 
# Autosomal recessive case
bcftools view -R ../chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed  trio${num}.vcf.gz | bcftools view -S ../samples.txt | bcftools view -i 'GT[0]="AA" && GT[1]="RA" && GT[2]="RA"' | bcftools filter -i 'QUAL>20' -Ov -o trio${num}.cand.vcf

elif [ "${num}" -eq  3 ]; then
# De Novo Autosomal Dominant case
bcftools view -R ../chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed  trio${num}.vcf.gz | bcftools view -S ../samples.txt | bcftools view -i 'GT[0]!="RR" && GT[1]="RR" && GT[2]="RR"' | bcftools filter -i 'QUAL>20' -Ov -o trio${num}.cand.vcf
elif [ "${num}" -eq  2 ]; then
#inherited Autosomal Domminant case
bcftools view -R ../chr20_ILMN_Exome_2.0_Plus_Panel.hg38_padded.bed  trio${num}.vcf.gz | bcftools view -S ../samples.txt | bcftools view -i 'GT[0] != "RR" && (GT[1] != "RR" || GT[2] != "RR")' | bcftools filter -i 'QUAL>20' -Ov -o trio${num}.cand.vcf

fi
