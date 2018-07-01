spades.py -o spades1_RZ_k127 --pe1-1 RZ_R1.fq.gz --pe1-2 RZ_R2.fq.gz --mp1-1 RZ_3000_R1.fq.gz --mp1-2 RZ_3000_R2.fq.gz -m 250 -t 20 -k 127
spades.py -o spades1_MH_k127 -1 KW_R1.fq.gz -2 KW_R2.fq.gz --pacbio m54169_170713_065211.subreads.wsb.fasta -m 250 -t 20 -k 127
canu -p DHcp -d DH -pacbio-raw ../raw_data/filtered_subreads.fastq genomeSize=0.9m gnuplotImageFormat=svg
quast -o quast -t 10 *fa

