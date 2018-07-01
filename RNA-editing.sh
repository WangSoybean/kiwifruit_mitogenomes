#create bam files
#!/bin/bash
set -eu -o pipefail

genome=RZmt.fasta

hisat2-build $genome mt
for n in `cat name`
do
reads1=${n}_1.fastq.gz
reads2=${n}_2.fastq.gz
hisat2 -x mt -1 $reads1 -2 $reads2 -S $n.sam -p 16 2>mt.log
samtools view -bS -F 12 -q 20 -@ 10 $n.sam | samtools sort -@ 10 - | samtools rmdup - --reference $genome $n.F12.sort.rmdup.bam
done
#create vcf files
freebayes -f RZmt.fasta Aa-1-1_clean.F12.sort.rmdup.bam > Aa-1-1_clean.F12.sort.rmdup.bam.vcf
freebayes -f DHmt.fasta AH-1-1_clean.F12.sort.rmdup.bam > AH-1-1_clean.F12.sort.rmdup.bam.vcf
freebayes -f MHmt.fasta Ae-1-1_clean.F12.sort.rmdup.bam > Ae-1-1_clean.F12.sort.rmdup.bam.vcf
#create RNA editing information
#!/bin/bash
set -eu -o pipefail

for n in *.F12.sort.rmdup.bam.vcf ; do cat $n | vcffilter -f "QUAL > 20" -f "DP > 10" > $n.filer ; done

for n in Aa*_clean.F12.sort.rmdup.bam.vcf.filer ; do vcftools --vcf $n --bed RZmt.cds.gff.bed --recode --out ${n/.vcf/}.RZ.cds ; done
for n in AH*_clean.F12.sort.rmdup.bam.vcf.filer ; do vcftools --vcf $n --bed DHmt1.cds.gff.bed --recode --out ${n/.vcf/}.DH1.cds ; done
for n in AH*_clean.F12.sort.rmdup.bam.vcf.filer ; do vcftools --vcf $n --bed DHmt2.cds.gff.bed --recode --out ${n/.vcf/}.DH2.cds ; done
for n in Ae*_clean.F12.sort.rmdup.bam.vcf.filer ; do vcftools --vcf $n --bed MHmt.cds.gff.bed --recode --out ${n/.vcf/}.MH.cds ; done
rm *.log

for m in Aa*_clean.F12.sort.rmdup.bam.filer.RZ.cds.recode.vcf ; do for n in `cat RZmt.cds.gff.bed | cut -f 2-4,6 | sed 's/\t/#/g'` ; do b=${n%%#*} ; b1=${n#*#}; e=${b1%%#*} ; name=${b1#*#}; cat $m | awk -v v0=$n -v v1=$b -v v2=$e -v v3=$name -F '\t|:' 'BEGIN{OFS="\t"}{split($19,a,",");if($2>=v1 && $2<=v2){print v0,$1,$2,$4,$5,$6,$17,$18,$19,a[1]/$18}}' ; done | awk -F '#|\t' 'BEGIN{OFS="\t"}{if($4~/+/){print $3,$6-$1,$15,$7,$8,$9,$10,$11,$12,$13,$14}else{print $3,$2-$6+1,$15,$7,$8,$9,$10,$11,$12,$13,$14}}' | sed 's/\t\t\t/\t/' | sed 's/\t\t/\t/' | sed 's/\t/_/' | sed 's/:/_/' | sort -t _ -k 1,1 -k 3n,3 > $m.editing ; done  #sed '1i\#Gene\tLoca\tRef\tAlt\tQual\tGenotype\tDP\tAD\tRatio'

for m in AH*_clean.F12.sort.rmdup.bam.filer.DH1.cds.recode.vcf ; do for n in `cat DHmt1.cds.gff.bed  | cut -f 2-4,6 | sed 's/\t/#/g'` ; do b=${n%%#*} ; b1=${n#*#}; e=${b1%%#*} ; name=${b1#*#}; cat $m | awk -v v0=$n -v v1=$b -v v2=$e -v v3=$name -F '\t|:' 'BEGIN{OFS="\t"}{split($19,a,",");if($2>=v1 && $2<=v2){print v0,$1,$2,$4,$5,$6,$17,$18,$19,a[1]/$18}}' ; done | awk -F '#|\t' 'BEGIN{OFS="\t"}{if($4~/+/){print $3,$6-$1,$15,$7,$8,$9,$10,$11,$12,$13,$14}else{print $3,$2-$6+1,$15,$7,$8,$9,$10,$11,$12,$13,$14}}' | sed 's/\t\t\t/\t/' | sed 's/\t\t/\t/' | sed 's/\t/_/' | sed 's/:/_/' | sort -t _ -k 1,1 -k 3n,3 > $m.editing ; done

for m in AH*_clean.F12.sort.rmdup.bam.filer.DH2.cds.recode.vcf ; do for n in `cat DHmt2.cds.gff.bed  | cut -f 2-4,6 | sed 's/\t/#/g'` ; do b=${n%%#*} ; b1=${n#*#}; e=${b1%%#*} ; name=${b1#*#}; cat $m | awk -v v0=$n -v v1=$b -v v2=$e -v v3=$name -F '\t|:' 'BEGIN{OFS="\t"}{split($19,a,",");if($2>=v1 && $2<=v2){print v0,$1,$2,$4,$5,$6,$17,$18,$19,a[1]/$18}}' ; done | awk -F '#|\t' 'BEGIN{OFS="\t"}{if($4~/+/){print $3,$6-$1,$15,$7,$8,$9,$10,$11,$12,$13,$14}else{print $3,$2-$6+1,$15,$7,$8,$9,$10,$11,$12,$13,$14}}' | sed 's/\t\t\t/\t/' | sed 's/\t\t/\t/' | sed 's/\t/_/' | sed 's/:/_/' | sort -t _ -k 1,1 -k 3n,3 > $m.editing ; done

for m in Ae*_clean.F12.sort.rmdup.bam.filer.MH.cds.recode.vcf ; do for n in `cat MHmt.cds.gff.bed  | cut -f 2-4,6 | sed 's/\t/#/g'` ; do b=${n%%#*} ; b1=${n#*#}; e=${b1%%#*} ; name=${b1#*#}; cat $m | awk -v v0=$n -v v1=$b -v v2=$e -v v3=$name -F '\t|:' 'BEGIN{OFS="\t"}{split($19,a,",");if($2>=v1 && $2<=v2){print v0,$1,$2,$4,$5,$6,$17,$18,$19,a[1]/$18}}' ; done | awk -F '#|\t' 'BEGIN{OFS="\t"}{if($4~/+/){print $3,$6-$1,$15,$7,$8,$9,$10,$11,$12,$13,$14}else{print $3,$2-$6+1,$15,$7,$8,$9,$10,$11,$12,$13,$14}}' | sed 's/\t\t\t/\t/' | sed 's/\t\t/\t/' | sed 's/\t/_/' | sed 's/:/_/' | sort -t _ -k 1,1 -k 3n,3 > $m.editing ; done
# annotation of RNA editing
java -jar snpEff.jar 1.RZmt_RNA_editing 1.vcf.files/4.compare.RNA.editing/Aa-1-1_clean.F12.sort.rmdup.bam.filer.RZ.cds.recode.vcf > 1.vcf.files/4.compare.RNA.editing/Aa-1-1_clean.F12.sort.rmdup.bam.filer.RZ.cds.recode.vcf.ann
cp snpEff_genes.txt RZ.gene.editing
java -jar snpEff.jar 2.MHmt_RNA_editing 1.vcf.files/4.compare.RNA.editing/Ae-1-1_clean.F12.sort.rmdup.bam.filer.MH.cds.recode.vcf > 1.vcf.files/4.compare.RNA.editing/Ae-1-1_clean.F12.sort.rmdup.bam.filer.MH.cds.recode.vcf.ann
cp snpEff_genes.txt MH.gene.editing
java -jar snpEff.jar 3.DHmt1_RNA_editing 1.vcf.files/4.compare.RNA.editing/AH-1-1_clean.F12.sort.rmdup.bam.filer.DH1.cds.recode.vcf > 1.vcf.files/4.compare.RNA.editing/AH-1-1_clean.F12.sort.rmdup.bam.filer.DH1.cds.recode.vcf.ann
cp snpEff_genes.txt DH1.gene.editing
java -jar snpEff.jar 4.DHmt2_RNA_editing 1.vcf.files/4.compare.RNA.editing/AH-1-1_clean.F12.sort.rmdup.bam.filer.DH2.cds.recode.vcf > 1.vcf.files/4.compare.RNA.editing/AH-1-1_clean.F12.sort.rmdup.bam.filer.DH2.cds.recode.vcf.ann
cp snpEff_genes.txt DH2.gene.editing



