samtools faidx E.coli_K12_MG1655.fa

bwa index E.coli_K12_MG1655.fa

bwa mem -t 24 -R '@RG\tID:K12\tSM:K12' \
    E.coli_K12_MG1655.fa SRR1770413_1.fastq SRR1770413_2.fastq \
    | samtools view -Shu - \
    | sambamba sort /dev/stdin -o SRR1770413.raw.bam
sambamba markdup SRR1770413.raw.bam SRR1770413.bam

freebayes -f E.coli_K12_MG1655.fa --ploidy 1 SRR1770413.bam >SRR1770413.vcf



awk 'BEGIN {RS=">";ORS=""} length($0) >= 1000 {print ">"$0}' rawdata/Montipora_assembly.fa > rawdata/Montipora_assembly_1000.fa

samtools faidx rawdata/Montipora_assembly_1000.fa
bwa index rawdata/Montipora_assembly_1000.fa


#### TRY ON MY DATAA

# take subset of 10000 reads
seqtk sample -s100 BC_R1_trimmed.fq 10000 > BCsub1.fq
seqtk sample -s100 BC_R2_trimmed.fq 10000 > BCsub2.fq

seqtk sample -s100 NC_R1_trimmed.fq 10000 > NCsub1.fq
seqtk sample -s100 NC_R2_trimmed.fq 10000 > NCsub2.fq

seqtk sample -s100 ND_R1_trimmed.fq 10000 > NDsub1.fq
seqtk sample -s100 ND_R2_trimmed.fq 10000 > NDsub2.fq


#alternative map reads to reference: RC 11.21.16
bwa mem -t 24 -R '@RG\tID:BC\tSM:BC' \
    rawdata/Montipora_assembly_1000.fa procdata/BCsub1.fq procdata/BCsub2.fq \
    | samtools view -Shu - \
    | sambamba sort /dev/stdin -o procdata/BCsub_raw.bam
bwa mem -t 24 -R '@RG\tID:NC\tSM:NC' \
    rawdata/Montipora_assembly_1000.fa procdata/NCsub1.fq procdata/NCsub2.fq \
    | samtools view -Shu - \
    | sambamba sort /dev/stdin -o procdata/NCsub_raw.bam
bwa mem -t 24 -R '@RG\tID:ND\tSM:ND' \
    rawdata/Montipora_assembly_1000.fa procdata/NDsub1.fq procdata/NDsub2.fq \
    | samtools view -Shu - \
    | sambamba sort /dev/stdin -o procdata/NDsub_raw.bam


ls procdata/*sub_raw.bam > procdata/sub_bam.list
bamtools merge -list procdata/sub_bam.list > procdata/cat_sub.bam

bamtools header -in procdata/cat_sub.bam | grep -P "@RG"

samtools index procdata/cat_sub.bam

freebayes -f rawdata/Montipora_assembly_1000.fa -b procdata/BCsub_raw.bam --pooled-continuous -p 28 --use-best-n-alleles 3 -v procdata/BCsub_raw.vcf

freebayes -f rawdata/Montipora_assembly_1000.fa -b procdata/BCsub_raw.bam --pooled-continuous -p 28 --use-best-n-alleles 3 -v procdata/BCsub_raw.vcf

