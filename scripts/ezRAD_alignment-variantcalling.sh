#!/bin/bash

#Alignment, Variant Calling, VCF output
#you will need the following programs installed:
#samtools, bwa, picard-tools, bamtools, freebayes, GATK, vcftools


#first you will need your reference to be in a fasta file
#index reference
samtools faidx reference.fasta
bwa index reference.fasta
java -jar picard-tools-1.108/CreateSequenceDictionary.jar R=reference.fasta O=reference.dict


######################################
##Use BWA to map reads to reference

#first declare an array of all your file names. For example, if you have 3 individuals, you might have files names like this:
#ind1.R1.fq ind1.R2.fq ind2.R1.fq ind2.R2.fq #ind3.R1.fq ind3.R2.fq
#so your array will be the prefix, like this:
declare -a NAMES=("ind1" "ind2" "ind3")

#then you can run a loop for all of your individuals at once. 
#the bwa mem command runs the paired-read alignment; the "@RG..." part adds a read group to each alignment (identifier for when you merge alignments)
#the samtools command then converts the .sam to .bam and does the first filter - which removes alignments that 
#had mapping quality scores of 0 (which means it didn't map to anything)
for i in "${NAMES[@]}"
do
bwa mem reference.fasta $i.R1.fq $i.R2.fq -t 16 -a -M -T 10 -R "@RG\tID:$i\tSM:$i\tPL:Illumina" 2> bwa.$i.log | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@16 -q 1 -SbT reference.fasta - > $i.bam 2>$i.bam.log
done
#use bwa -h to see what each option means

#OR if the loop isn't working for you just list each command separately
bwa mem chaetodon_reference.fasta ind1.R1.fq ind1.R2.fq -t 16 -a -M -T 10 -R "@RG\tID:ind1\tSM:ind1\tPL:Illumina" 2> bwa.ind1.log | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@16 -q 1 -SbT reference.fasta - > ind1.bam 2>ind1.bam.log
bwa mem chaetodon_reference.fasta ind2.R1.fq ind2.R2.fq -t 16 -a -M -T 10 -R "@RG\tID:ind2\tSM:ind2\tPL:Illumina" 2> bwa.ind2.log | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@16 -q 1 -SbT reference.fasta - > ind2.bam 2>ind2.bam.log
bwa mem chaetodon_reference.fasta ind3.R1.fq ind3.R2.fq -t 16 -a -M -T 10 -R "@RG\tID:ind3\tSM:ind3\tPL:Illumina" 2> bwa.ind3.log | mawk '!/\t[2-9].[SH].*/' | mawk '!/[2-9].[SH]\t/' | samtools view -@16 -q 1 -SbT reference.fasta - > ind3.bam 2>ind3.bam.log


#Now, sort each bam file
samtools sort ind1.bam ind1_sort
samtools sort ind2.bam ind2_sort
samtools sort ind3.bam ind3_sort

#run diagnostics on alignments
#you are looking to see that a good chunk of alignments were properly paired (R1 and R2 mapped to same region in reference)
samtools flagstat ind1_sort.bam > ind1_sort.flagstat
samtools flagstat ind2_sort.bam > ind2_sort.flagstat
samtools flagstat ind3_sort.bam > ind3_sort.flagstat


#run filters on .bam mappings
#these settings will remove alignments that were not mapped, not properly paired, and that had mapping qualities < 30
#Extract properly paired (-f 0x02) AND exclude un-mapped reads (-F 0x04) (PM = paired mapped)
samtools view -f 0x02 -F 0x04 -q30 -b ind1_sort.bam > ind1_sortPM.bam
samtools view -f 0x02 -F 0x04 -q30 -b ind2_sort.bam > ind2_sortPM.bam
samtools view -f 0x02 -F 0x04 -q30 -b ind3_sort.bam > ind3_sortPM.bam

#run flagstat diagnostics again to see impact of filters and adjust if needed
samtools flagstat ind1_sortPM.bam > ind1_sortPM.flagstat
samtools flagstat ind2_sortPM.bam > ind2_sortPM.flagstat
samtools flagstat ind3_sortPM.bam > ind3_sortPM.flagstat


##### MERGE bam files & Review Mapping Stats #####

#Create list of BAM files
ls *_sort.bam > bamlist_sort.list
ls *_sortPM.bam > bamlist_sortPM.list
#merge bamfiles together (this is the appropriate time to remove duplicates)
bamtools merge -list bamlist_sort.list > cat_sort.bam
bamtools merge -list bamlist_sortPM.list > cat_sortPM.bam
#validate that both read groups (RGs) are in merged file
bamtools header -in cat_sortPM.bam  | grep -P "@RG" > RG_validate.log 
#MAPPING STATS - All individuals 
more *_sort.flagstat > flagstat_pre-filter.txt
more *_sortPM.flagstat > flagstat_post-filter.txt

#index merged .bam
samtools index cat_sort.bam
samtools index cat_PM.bam

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
java -jar GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar -R reference.fasta \
 -T RealignerTargetCreator  \
 -I indcat_PM.bam -o cat_realigner.intervals -nt 16
#you will notice that GATK runs internal filters that are pretty equivalent to the ones we did above.
#Once you're familar with the filters, and tested how increasing or decreasing thresholds impacts your output,
#you can just as easily run this and skip the above filtering.


#Indel Realigner
java -XX:MaxPermSize=64g -jar /home/jw2/programs/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar -T IndelRealigner  -R reference.fasta \
     -I indcat_sort.bam  \
     -targetIntervals ind_sort_realigner.intervals \
     -o indcat_sort_realigned.bam
#Realigner NOTES:
#you might get some warnings about intervals that have too many reads, e.g., 	
	#IndelRealigner - Not attempting realignment in interval outFile_10032_654:9-458 because there are too many reads.
#in my experience these are usually PhiX contamination or messy areas with really high coverage.
#so its good that it ignores them.

#index realigned .bam
samtools index cat_sort_realigned.bam

#ALTERNATIVE to fixing INDELS
#if GATK isn't working for you, samtools has an alternative called BAQ (Base Alignment Quality), 
#which helps avoid the INDEL artifact in low-coverage multi-sample SNP calling
#The BAQ strategy is invoked by default in mpileup. To make other SNP callers take advantage of BAQ, you could run:
samtools calmd -Abr cat_sortPM.bam reference.fasta > cat_sortPM_baq.bam  
samtools index cat_sortPM_baq.bam

#########################################################################################
#Variant calling can be done using various tools, but freebayes is a good choice


########### Variant Calling with FreeBayes ##############################################
#take the time to look at freebayes --help to familiarize yourself with all the settings.
#for example, the -O option will filter input alignments just like we did using samtools view above. 
#so as you can see there are many ways (and many places) for filtering poor mappings.  

#we can use 3 variant calling methods for pooled samples: J (discrete), K (continuous), and Naive.  In short, -K (pooled-continous) is the best method that generates frequency-based calls for all variants. J-mode calculates genotypes from pooled data and looks for variants between pools (this is over-simplified and produces very few SNPS - highly underestimated), and Naive takes way too long and calls a SNP for each site and counts #s of each variant (it will likely overestimate SNPs)
#K-mode Continuous-pooling (no-ploidy set and ignores genotype outcome). Generates frequency-based calls for all variants passing input thresholds
#This is a command for using pooled samples. To run on individuals remove the "--pooled-continuous" option
freebayes -f reference.fasta -0 -E 3 -C 4 -G 10 -z .1 -X -u -n 4 --pooled-continuous --min-coverage 10 --min-repeat-entropy 1 -V --populations poplist -b cat_sort_realigned.bam -v cat_fbraw.vcf

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
