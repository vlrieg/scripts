#!/bin/bash

# This script is supposed to be based on https://software.broadinstitute.org/gatk/documentation/article.php?id=6483
# Run like `sbatch --mem=8G -c10 ./variant-calling.sh ERR2679009`

module load bwa
module load samtools
module load jdk/1.8.0_45-fasrc01

# Exit immediately if any command returns a failing exit status
set -e

#PvP01 vivax reference genome
reference=/data/wraycompute/vdp5/reference_data/PVP01.fasta
ref=${reference}

##### 1. set up variables #####
base=$1
shift

echo "analyzing data set $base"
echo "creating work dir $base"
mkdir -p $base
cd $base

Read="../${base}.fastq.gz"
echo "read for ${base} is ${Read}"
#if [ ! -f ${Read} ] ; then
#echo 'read file ${Read} not found! quitting"
#exit 1
#fi

#R1="../${base}_1.fastq.gz"
#echo "forward read for $base is $R1"

#if [ ! -f ${R1} ] ; then
#echo "forward read file ${R1} not found! quitting"
#exit 1
#fi
    
#R2="../${base}_2.fastq.gz"
#echo "reverse read for $base is $R2"

#if [ ! -f ${R2} ] ; then
#echo "reverse read file ${R2} not found! quitting"
#exit 1
#fi


RUN_PICARD="java -jar -Xmx7g /data/wraycompute/malaria/Applications/picard/picard.jar"
RUN_GATK3="java -jar /home/vdp5/src/GenomeAnalysisTK.jar"

################### mapping ###################

function Map_Fastq ()
{
echo "mapping fastq.gz files: start"
bwa mem -M -t 4 -v 2 -A 2 -L 15 -U 9 -T 75 -k 19 -w 100 -d 100 -r 1.5 -c 10000 -B 4 -O 6 -E 1 ${ref} ${Read} | samtools view  -Sb - | samtools sort - -o ${base}.bam
echo "mapping fastq.gz files: done"
}
Map_Fastq

################### variant processing ###################

#put the output in my directory for now
# valdirectory=/data/wraycompute/val/variants
# valdir=${valdirectory}

# Add read groups 
# necessary step for downstream processing. The validation stringency is necessary for malaria since its genome is full of contigs
function ReadGroups ()
{
echo "Adding read group info: start"
${RUN_PICARD} AddOrReplaceReadGroups INPUT=${base}.bam OUTPUT=${base}-readgroups.bam RGID=1 RGLB=${base} RGPL=illumina RGPU=unit1 RGSM=${base}  VALIDATION_STRINGENCY=LENIENT
echo "Adding read group info: done"
}
ReadGroups #call function

# this step is necessary for GATK, but note that it is NOT for tools like freebayes since they do their own duplicates processing.
function MarkDuplicates ()
{
echo "Mark Duplicates & Indexing: start"
${RUN_PICARD} MarkDuplicates INPUT=${base}-readgroups.bam OUTPUT=${base}-dedup_reads.bam METRICS_FILE=metrics.txt  VALIDATION_STRINGENCY=LENIENT
rm -rf readgroups.bam
${RUN_PICARD} BuildBamIndex INPUT=${base}-dedup_reads.bam  VALIDATION_STRINGENCY=LENIENT

echo "Mark Duplicates & Indexing: Done!"
}
MarkDuplicates #call the function!


# initial variant processing, using 4 cores. Note that we call full variants here, not the GVCF (which is what we want as the final output)
function RawVariants ()
{
echo "GATK3 raw variants: start..."
${RUN_GATK3} -T HaplotypeCaller -R ${ref} -I ${base}-dedup_reads.bam -o ${base}-raw_variants.vcf -nct 4
echo "GATK3 raw variants: done!"
}
RawVariants #call the function

# we split up the variant types since SNPs and indels need different processing
function RawSnps ()
{
echo "GATK3 raw SNPs: start..."
${RUN_GATK3} -T SelectVariants -R ${ref} -V ${base}-raw_variants.vcf -selectType SNP -o ${base}-raw_snps.vcf
echo "GATK3 raw SNPs: done!"
}
RawSnps #call the function


function RawIndels ()
{
echo "GATK3 raw Indels: start..."
${RUN_GATK3} -T SelectVariants -R ${ref} -V ${base}-raw_variants.vcf -selectType INDEL -o ${base}-raw_indels.vcf
rm -rf ${base}-raw_variants.vcf
echo "GATK3 raw Indels: done!"
}
RawIndels #call the function

#basic filters
function FilterSnp ()
{
echo "GATK3 filter SNPs: starting..."
${RUN_GATK3} -T VariantFiltration -R ${ref} -V ${base}-raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o ${base}-filtered_snps.vcf
echo "GATK3 filter SNPs: done!"
}
FilterSnp #call the function

function FilterIndel ()
{
echo "GATK3 filter indels: starting..."
${RUN_GATK3} -T VariantFiltration -R ${ref} -V ${base}-raw_indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o ${base}-filtered_indels.vcf
echo "GATK3 filter indels: done!!!!!"
}
FilterIndel # call the function

# Here is our first recalibration step. We use the hard filtered SNPs set as the training set. It's not an issue that this set isn't totally accurate because we're going to do more rounds of processing later.
function BaseRecalibrator1 ()
{
echo "GATK3 base recalibrator round 1: start"
${RUN_GATK3} -T BaseRecalibrator -R ${ref} -I ${base}-dedup_reads.bam -knownSites ${base}-filtered_snps.vcf -knownSites ${base}-filtered_indels.vcf -o ${base}-recal_data.table
${RUN_GATK3} -T PrintReads -R ${ref} -I ${base}-dedup_reads.bam -BQSR ${base}-recal_data.table -o ${base}-recal_reads.bam
echo "GATK3 base recalibrator round 1: done"
}
BaseRecalibrator1 #call the function

#Round 2 
# we call the second set here using the recalculated reads (variants2)
function SecondSet()
{
echo 'GATK3 second set: start'
${RUN_GATK3} -T HaplotypeCaller -nct 4 -R ${ref} -I ${base}-recal_reads.bam -o ${base}-raw_variants2.vcf
echo "GATK3 second set: done"
}
SecondSet #call the function

#some cleanup
function CleanUp()
{
rm -rf *.idx ${base}-*.idx ${base}-raw_snps.vf ${base}-raw_indels.vcf ${base}-recal_reads.bam ${base}-recal_data.table ${base}-post_recal_data.table ${base}-filtered_snps.vcf ${base}-filtered_indels.vcf
echo "cleaned up some intermediate files"
}
CleanUp #call the function

function RawSnps2()
{
echo "GATK3 call raw snps round 2: start"
${RUN_GATK3} -T SelectVariants -R ${ref} -V ${base}-raw_variants2.vcf -selectType SNP -o ${base}-raw_snps.vcf
echo "GATK3 call raw snps round 2: end"
}
RawSnps2 #call the function

function RawIndels2()
{
echo "GATK3 call raw Indels round 2: start"
${RUN_GATK3} -T SelectVariants -R ${ref} -V ${base}-raw_variants2.vcf -selectType INDEL -o ${base}-raw_indels.vcf
rm -rf ${base}-raw_variants2.vcf
echo "GATK3 call raw Indels round 2: end"
}
RawIndels2

function FilterSnps2()
{
echo "GATK3 call filter snps round 2: start"
${RUN_GATK3} -T VariantFiltration -R ${ref} -V ${base}-raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o ${base}-filtered_snps.vcf
echo "GATK3 call filter snps round 2: done"
}
FilterSnps2 #call the function

function FilterIndels2()
{
echo "GATK3 call filter indels round 2: start"
${RUN_GATK3} -T VariantFiltration -R ${ref} -V ${base}-raw_indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o ${base}-filtered_indels.vcf
echo "GATK3 call filter indels round 2: done"
}
FilterIndels2 # call the function

# really important to note that we use the dedup_reads file here as the input, not the printed reads from the last step. I have not tested what would happen if we used the recalculated reads here instead, but none of the guides I read seemed to do so.
function Recalibrate2()
{
echo "GATK3 recalibrate round 2: start"
${RUN_GATK3} -T BaseRecalibrator -R ${ref} -I ${base}-dedup_reads.bam  -knownSites ${base}-filtered_snps.vcf -knownSites ${base}-filtered_indels.vcf -o ${base}-recal_data.table
${RUN_GATK3} -T PrintReads -R ${ref} -I ${base}-dedup_reads.bam -BQSR ${base}-recal_data.table -o ${base}-recal_reads.bam
echo "GATK3 recalibrate round 2: done"
}
Recalibrate2 # call the function

#Round 3
function ThirdSet()
{
echo "GATK3 third set: start"
${RUN_GATK3} -T HaplotypeCaller -nct 4 -R ${ref} -I ${base}-recal_reads.bam -o ${base}-raw_variants3.vcf
echo "GATK3 third set: done"
}
ThirdSet #call function

function RawSnps3()
{
echo "calling raw SNPs round 3: start"
${RUN_GATK3} -T SelectVariants -R ${ref} -V ${base}-raw_variants3.vcf -selectType SNP -o ${base}-raw_snps.vcf
echo "calling raw SNPs round 3: done"
}
RawSnps3 #call the function


function RawIndels3()
{
echo "calling raw indels round 3: start"
${RUN_GATK3} -T SelectVariants -R ${ref} -V ${base}-raw_variants3.vcf -selectType INDEL -o ${base}-raw_indels.vcf
rm -rf ${base}-raw_variants3.vcf
echo "calling raw indels round 3: done"
}
RawIndels3 #call the function


function FilterSnps3()
{
echo "filter snps round 3: start"
${RUN_GATK3} -T VariantFiltration -R ${ref} -V ${base}-raw_snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 4.0' --filterName "basic_snp_filter" -o ${base}-filtered_snps.vcf
echo "filter snps round 3: done"
}
FilterSnps3 #call function

function FilterIndels3()
{
echo "filter indels round 3: start"
${RUN_GATK3} -T VariantFiltration -R ${ref} -V ${base}-raw_indels.vcf --filterExpression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filterName "basic_indel_filter" -o ${base}-filtered_indels.vcf
echo "filter indels round 3: done"
}
FilterIndels3

# really important to note that we use the dedup_reads file here as the input, not the printed reads from the last step. I have not tested what would happen if we used the recalculated reads here instead, but none of the guides I read seemed to do so.
function Recalibrate3()
{
echo "recalibration round 3: start"
${RUN_GATK3} -T BaseRecalibrator -R ${ref} -I ${base}-dedup_reads.bam  -knownSites ${base}-filtered_snps.vcf -knownSites ${base}-filtered_indels.vcf -o ${base}-recal_data.table
${RUN_GATK3} -T PrintReads -R ${ref} -I ${base}-dedup_reads.bam -BQSR ${base}-recal_data.table -o ${base}-recal_reads.bam
echo "recalibration round 3: done"
}
Recalibrate3 #call function


#Last Step
function Output()
{
echo "outputting GVCF file: start"
${RUN_GATK3} -T HaplotypeCaller -nct 4 -R ${ref} -I ${base}-recal_reads.bam --emitRefConfidence GVCF -o ${base}.g.vcf
echo "outputting GVCF file: done"
}
Output #call function


#some cleanup
function FinalCleanUp()
{
rm -rf *.bam *.bai ${base}-raw_snps.vcf ${base}-raw_indels.vcf ${base}-recal_reads.bam ${base}-recal_data.table ${base}-post_recal_data.table ${base}-filtered_snps.vcf ${base}-filtered_indels.vcf
rm -rf ${base}-filtered_indels.vcf.idx ${base}-filtered_snps.vcf.idx ${base}-raw_indels.vcf.idx ${base}-raw_snps.vcf.idx ${base}-raw_variants3.vcf.idx ${base}-raw_variants2.vcf.idx metrics.txt ../${base}_1.fastq.gz ../${base}_2.fastq.gz
echo "clean up large intermediate files: DONE!"
}
FinalCleanUp #call the function


# This discusses handling primary and secondary alignments:
# https://www.biostars.org/p/206396/

# Base Quality Score Recalibration (BQSR)
# https://software.broadinstitute.org/gatk/documentation/article?id=2801

# GATK4 Tutorial on combining gvcfs for joint calling
# https://software.broadinstitute.org/gatk/documentation/article?id=11813
