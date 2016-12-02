

## Move files from permanent archive into project directory and rename
mkdir ~/McRAD/rawdata
cp ~/RCunning-ezRAD-32557555/BC-39972135/BC_S1_L001_R1_001.fastq.gz ~/McRAD/rawdata/BC_R1.fq.gz
cp ~/RCunning-ezRAD-32557555/BC-39972135/BC_S1_L001_R2_001.fastq.gz ~/McRAD/rawdata/BC_R2.fq.gz
cp ~/RCunning-ezRAD-32557555/NC-39977071/NC_S2_L001_R1_001.fastq.gz ~/McRAD/rawdata/NC_R1.fq.gz
cp ~/RCunning-ezRAD-32557555/NC-39977071/NC_S2_L001_R2_001.fastq.gz ~/McRAD/rawdata/NC_R2.fq.gz
cp ~/RCunning-ezRAD-32557555/ND-39981033/ND_S3_L001_R1_001.fastq.gz ~/McRAD/rawdata/ND_R1.fq.gz
cp ~/RCunning-ezRAD-32557555/ND-39981033/ND_S3_L001_R2_001.fastq.gz ~/McRAD/rawdata/ND_R2.fq.gz

## Run fastqc for quality summary of raw reads
#cd ~/McRAD
#fastqc rawdata/*.fq.gz

## Trim adapters from read 1 (may contain indexed adapter sequence at 3' end)
cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
-e 0.10 -O 3 rawdata/BC_R1.fq.gz -o procdata/BC_R1_trimmed.fq
cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
-e 0.10 -O 3 rawdata/NC_R1.fq.gz -o procdata/NC_R1_trimmed.fq
cutadapt -a GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
-e 0.10 -O 3 rawdata/ND_R1.fq.gz -o procdata/ND_R1_trimmed.fq

## Trim adapters from read 2 (may contain universal adapter sequence (reverse complemented) at 3' end)
cutadapt -a GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
-e 0.10 -O 3 rawdata/BC_R2.fq.gz -o procdata/BC_R2_trimmed.fq
cutadapt -a GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
-e 0.10 -O 3 rawdata/NC_R2.fq.gz -o procdata/NC_R2_trimmed.fq
cutadapt -a GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
-e 0.10 -O 3 rawdata/ND_R2.fq.gz -o procdata/ND_R2_trimmed.fq

## Index reference assembly
samtools faidx rawdata/Montipora_assembly.fa
bwa index rawdata/Montipora_assembly.fa
picard-tools CreateSequenceDictionary R=rawdata/Montipora_assembly.fa O=rawdata/Montipora_assembly.dict

## Map reads to reference
bwa mem -t 24 rawdata/Montipora_assembly.fa procdata/BC_R1_trimmed.fq procdata/BC_R2_trimmed.fq -t 16 -a -M -T 10 -R "@RG\tID:BC\tSM:BC\tPL:Illumina" 2> procdata/BC.bwa.log | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@16 -q 1 -SbT rawdata/Montipora_assembly.fa - > procdata/BC.bam 2>procdata/BC.bam.log
bwa mem rawdata/Montipora_assembly.fa procdata/NC_R1_trimmed.fq procdata/NC_R2_trimmed.fq -t 16 -a -M -T 10 -R "@RG\tID:NC\tSM:NC\tPL:Illumina" 2> procdata/NC.bwa.log | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@16 -q 1 -SbT rawdata/Montipora_assembly.fa - > procdata/NC.bam 2>procdata/NC.bam.log
bwa mem rawdata/Montipora_assembly.fa procdata/ND_R1_trimmed.fq procdata/ND_R2_trimmed.fq -t 16 -a -M -T 10 -R "@RG\tID:ND\tSM:ND\tPL:Illumina" 2> procdata/ND.bwa.log | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@16 -q 1 -SbT rawdata/Montipora_assembly.fa - > procdata/ND.bam 2>procdata/ND.bam.log




# is it not working because too many contigs in reference? try filtering out small contigs
awk 'BEGIN {RS=">";ORS=""} length($0) >= 500 {print ">"$0}' rawdata/Montipora_assembly.fa > rawdata/Montipora_assembly_500.fa











#Sort each bam file
samtools sort procdata/BC.bam -o procdata/BC_sort.bam
samtools sort procdata/NC.bam -o procdata/NC_sort.bam
samtools sort procdata/ND.bam -o procdata/ND_sort.bam

#run diagnostics on alignments
#you are looking to see that a good chunk of alignments were properly paired (R1 and R2 mapped to same region in reference)
samtools flagstat procdata/BC_sort.bam > procdata/BC_sort.flagstat
samtools flagstat procdata/NC_sort.bam > procdata/NC_sort.flagstat
samtools flagstat procdata/ND_sort.bam > procdata/ND_sort.flagstat

#run filters on .bam mappings
#these settings will remove alignments that were not mapped, not properly paired, and that had mapping qualities < 30
#Extract properly paired (-f 0x02) AND exclude un-mapped reads (-F 0x04) (PM = paired mapped)
samtools view -f 0x02 -F 0x04 -q30 -b procdata/BC_sort.bam > procdata/BC_sortPM.bam
samtools view -f 0x02 -F 0x04 -q30 -b procdata/NC_sort.bam > procdata/NC_sortPM.bam
samtools view -f 0x02 -F 0x04 -q30 -b procdata/ND_sort.bam > procdata/ND_sortPM.bam

#run flagstat diagnostics again to see impact of filters and adjust if needed
samtools flagstat procdata/BC_sortPM.bam > procdata/BC_sortPM.flagstat
samtools flagstat procdata/NC_sortPM.bam > procdata/NC_sortPM.flagstat
samtools flagstat procdata/ND_sortPM.bam > procdata/ND_sortPM.flagstat

##### MERGE bam files & Review Mapping Stats #####

#Create list of BAM files
ls procdata/*_sort.bam > procdata/bamlist_sort.list
ls procdata/*_sortPM.bam > procdata/bamlist_sortPM.list
#merge bamfiles together (this is the appropriate time to remove duplicates)
bamtools merge -list procdata/bamlist_sort.list > procdata/cat_sort.bam
bamtools merge -list procdata/bamlist_sortPM.list > procdata/cat_sortPM.bam
#validate that both read groups (RGs) are in merged file
bamtools header -in procdata/cat_sortPM.bam  | grep -P "@RG" > procdata/RG_validate.log 
#MAPPING STATS - All individuals 
more procdata/*_sort.flagstat > procdata/flagstat_pre-filter.txt
more procdata/*_sortPM.flagstat > procdata/flagstat_post-filter.txt

#index merged .bam
samtools index procdata/cat_sort.bam
samtools index procdata/cat_sortPM.bam

#Use GATK to realign mappings around INDELS (this makes mapping of INDELs consistent among multiple alignments - so they stack properly
#For example, say ind1 and ind2 both have a deletion of 1 of the A's in this sequence,
#CCTAATG - reference
#CCT*ATG - ind1 .bam mapped it like this
#CCTA*TG - ind2 .bam mapped it like this
#The realigner will go back through and compare alignments to resolve these differences so they line up correctly:
#CCTAATG - reference
#CCTA*TG - ind1 (realigned)
#CCTA*TG - ind2 (realigned)
#this is really important when we are doing variant calling, so this is read as 1 variant site and not 2.


#First, you'll need to create your realigner intervals
java -jar $GATK -R rawdata/Montipora_assembly.fa -T RealignerTargetCreator -I procdata/cat_sortPM.bam -o procdata/cat_realigner.intervals -L
### Can't get this to work....

# Try alternate approach to fixing indels
samtools calmd -Abr procdata/cat_sortPM.bam rawdata/Montipora_assembly.fa > procdata/cat_sortPM_baq.bam  
samtools index procdata/cat_sortPM_baq.bam
# This worked.....

########### Variant Calling with FreeBayes ##############################################
#take the time to look at freebayes --help to familiarize yourself with all the settings.
#for example, the -O option will filter input alignments just like we did using samtools view above. 
#so as you can see there are many ways (and many places) for filtering poor mappings.  

#we can use 3 variant calling methods for pooled samples: J (discrete), K (continuous), and Naive.  In short, -K (pooled-continous) is the best method that generates frequency-based calls for all variants. J-mode calculates genotypes from pooled data and looks for variants between pools (this is over-simplified and produces very few SNPS - highly underestimated), and Naive takes way too long and calls a SNP for each site and counts #s of each variant (it will likely overestimate SNPs)
#K-mode Continuous-pooling (no-ploidy set and ignores genotype outcome). Generates frequency-based calls for all variants passing input thresholds
#This is a command for using pooled samples. To run on individuals remove the "--pooled-continuous" option
freebayes -f rawdata/Montipora_assembly.fa -b procdata/cat_sort.bam -v procdata/cat_sort.vcf -p 28 --pooled-continuous
freebayes --fasta-reference rawdata/Montipora_assembly.fa -b procdata/cat_sortPM_baq.bam > test.vcf

#-X no multi-nucleotide polymorphisms
#-O Use stringent input base and mapping quality filters (Equivalent to -m 30 -q 20 -R 0 -S 0)
# -m 5 = mapping quality threshold
# -q = quality score of base (PHRED) to be called
# -E = complex gap 
# -C = min number observations of ALT allele (per individual)
# -G = at least N observations of ALT allele (for total population)
# -z = mismatch threshold (proportion) 
# -Q = mismatch base quality 


############### VcfLib Tools to work with VCF files ##################
#Now that you have the vcf file you can use vcffilter to reduce the vcf file to include whatever you want.
#Personally, I like to bring my raw.vcf into R, where I can test various filters, and then output a vcf that includes only variant sites passing all criteria

##### remove complex and MNPs
vcffilter -f "TYPE = snp | TYPE = ins | TYPE = del & NUMALT < 4" cat_fbraw.vcf > cat_snps-indels.vcf 

#filter by genotype qualities
vcffilter -f "QUAL > 20" cat_fbraw.vcf > cat_q20.vcf 

### look at just biallelic snps
vcffilter -f "TYPE = snp & NUMALT < 2" cat_fbraw.vcf > cat_bisnps.vcf 

#Filter to include only SNPs with total Depth > 20
vcffilter -f "DP > 20" cat_fbraw.vcf  > cat_20x.vcf 


################ What to do once you have properly filtered VCF file? ############################
#if you have individually-barcoded data you can use PGDspider to convert your vcf file into a number 
#of formats for use with other programs (e.g., Structure, Arlequin, GenePop, Bayescan, etc.)
#if you have pooled samples you have less options, but the best one is PoPoolation2, which you can use 
#to calculate locus-by-locus pairwise Fst comparisons (and p-values) for each population pair. 
#I have written an R script that will create the necessary .sync input file for PoPoolation and help
#you do quality control filtering on your vcf file output. It is by no means a "program", its a mess of 
#commands, but I'd be happy to show you how to modify it for your data. I know that you can do the same
#things in command line or python, but I know R better so thats what I used.

#If you have pooled samples, come talk to me and I'll give you the R script, and additional PoPoolation script.

