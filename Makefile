data/vcf/cat.vcf: data/bam/cat.bam data/ref/Montipora_assembly_1000.fa.fai
	samtools index data/bam/cat.bam
	# Variant calling using freebayes on merged bam file
	freebayes -f data/ref/Montipora_assembly_1000.fa -b data/bam/cat.bam --pooled-continuous -p 28 --use-best-n-alleles 3 -v $@

data/bam/cat.bam: data/bam/BC__raw.bam data/bam/NC__raw.bam data/bam/ND__raw.bam
	ls data/bam/*raw.bam > data/bam/__bam.list
	# Merge bam files from each pooled library together
	bamtools merge -list data/bam/__bam.list > $@

data/bam/BC__raw.bam: data/trimmed_fq/BC_R1_trimmed.fq data/trimmed_fq/BC_R2_trimmed.fq data/ref/Montipora_assembly_1000.fa.bwt
	# Map reads to reference
	bwa mem -t 24 -R '@RG\tID:BC\tSM:BC' \
	data/ref/Montipora_assembly_1000.fa $(word 1,$^) $(word 2,$^) \
	| samtools view -Shu - \
	| sambamba sort /dev/stdin -o $@

data/bam/NC__raw.bam: data/trimmed_fq/NC_R1_trimmed.fq data/trimmed_fq/NC_R2_trimmed.fq data/ref/Montipora_assembly_1000.fa.bwt
	# Map reads to reference
	bwa mem -t 24 -R '@RG\tID:NC\tSM:NC' \
	data/ref/Montipora_assembly_1000.fa $(word 1,$^) $(word 2,$^) \
	| samtools view -Shu - \
	| sambamba sort /dev/stdin -o $@

data/bam/ND__raw.bam: data/trimmed_fq/ND_R1_trimmed.fq data/trimmed_fq/ND_R2_trimmed.fq data/ref/Montipora_assembly_1000.fa.bwt
	# Map reads to reference
	bwa mem -t 24 -R '@RG\tID:ND\tSM:ND' \
	data/ref/Montipora_assembly_1000.fa $(word 1,$^) $(word 2,$^) \
	| samtools view -Shu - \
	| sambamba sort /dev/stdin -o $@

data/trimmed_fq/%_R1_trimmed.fq: data/raw_fq/%_R1.fq.gz
	# Trim adapters from read 1 (may contain indexed adapter sequence at 3' end)
	cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG -e 0.10 -O 3 $< -o $@

data/trimmed_fq/%_R2_trimmed.fq: data/raw_fq/%_R2.fq.gz
	# Trim adapters from read 2 (may contain universal adapter sequence (reverse complemented) at 3' end)
	cutadapt -a GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -e 0.10 -O 3 $< -o $@

data/ref/Montipora_assembly_1000.fa.bwt: data/ref/Montipora_assembly_1000.fa
	# Index reference assembly
	bwa index $<

data/ref/Montipora_assembly_1000.fa.fai: data/ref/Montipora_assembly_1000.fa
	samtools faidx $<
