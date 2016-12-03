#map reads to reference: RC 11.21.16
bwa mem -t 24 -R '@RG\tID:BC\tSM:BC' \
    rawdata/Montipora_assembly_1000.fa procdata/BC_R1_trimmed.fq procdata/BC_R2_trimmed.fq \
    | samtools view -Shu - \
    | sambamba sort /dev/stdin -o procdata/BC__raw.bam
bwa mem -t 24 -R '@RG\tID:NC\tSM:NC' \
    rawdata/Montipora_assembly_1000.fa procdata/NC_R1_trimmed.fq procdata/NC_R2_trimmed.fq \
    | samtools view -Shu - \
    | sambamba sort /dev/stdin -o procdata/NC__raw.bam
bwa mem -t 24 -R '@RG\tID:ND\tSM:ND' \
    rawdata/Montipora_assembly_1000.fa procdata/ND_R1_trimmed.fq procdata/ND_R2_trimmed.fq \
    | samtools view -Shu - \
    | sambamba sort /dev/stdin -o procdata/ND__raw.bam


echo "procdata/BC__raw.bam
procdata/NC__raw.bam
procdata/ND__raw.bam" > procdata/__bam.list

bamtools merge -list procdata/__bam.list > procdata/cat__raw.bam

#bamtools header -in procdata/cat__raw.bam | grep -P "@RG"

#samtools sort procdata/cat__raw.bam -o procdata/catsort.bam
samtools index procdata/cat__raw.bam

freebayes -f rawdata/Montipora_assembly_1000.fa -b procdata/cat__raw.bam --pooled-continuous -p 28 --use-best-n-alleles 3 -v procdata/cat__raw.vcf

# can freebayes run parallel?



# this script did not work when rawdata/Montipora_assembly.fa was used. retry with Montipora_assembly_1000.fa
# script WORKS when using Montipora_assembly_1000.fa!


samtools faidx data/ref/Montipora_assembly_1000.fa
bwa index data/ref/Montipora_assembly_1000.fa

# Create input files for popoolation analysis
samtools mpileup -B data/bam/BC__raw.bam data/bam/NC__raw.bam > data/popool/BCNC.mpileup
java -ea -Xmx7g -jar ~/popoolation2_1201/mpileup2sync.jar --input data/popool/BCNC.mpileup --output data/popool/BCNC.sync --fastq-type sanger --min-qual 20 --threads 16

perl ~/popoolation2_1201/snp-frequency-diff.pl --input data/popool/BCNC.sync --output-prefix BCNC --min-count 6 --min-coverage 50 --max-coverage 200
