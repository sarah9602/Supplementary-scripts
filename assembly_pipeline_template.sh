#!/bin/bash
#
#SBATCH -J assembly_SAMPLE
#SBATCH -o %x-%j.o ## x denotes job name and j denotes job ID
#SBATCH -e %x-%j.e
#SBATCH -p PARTITION-NAME ## Unique to HPC system
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH -t 1-00:00:00
#SBATCH --mem=0
#SBATCH --mail-type=END
#SBATCH --mail-user=sarahjjohnson@ou.edu

module load MEGAHIT/1.2.9-foss-2019a-Python-3.8.0

## Create a variable to your working directory
working=/scratch/sjj9602/Assembly_SAMPLE
samp=SAMPLE
mkdir -p $working
cd $working

#############---------- Megahit Assembly
## Add continue flag in case a directory already exists
megahit --k-min 21 --k-max 91 --k-step 10 -r INFILE -o Megahit --out-prefix $samp -t 8 --memory 0.9 --continue
mv Megahit/$samp.contigs.fa .

#############---------- Filter reads < 1000 nt in length
cat $samp.contigs.fa|seqkit seq -m 1000 -w 0 > $samp.contigs.filt.fa

#############---------- Read mapping
module load Bowtie2/2.3.5.1-GCC-8.3.0

echo "building bowtie2 index"
bowtie2-build $samp.contigs.filt.fa $samp.ref 
## Input file should be fastq - align filtered fastq data to contig reference
echo "contig alignment"
## For ancient samples add --very-sensitive flag and -N 1 to allow for 
## additional mismatch for presence of ancient DNA damage
bowtie2 -p 5 --very-sensitive -N 1 -x $samp.ref -U INFILE | \
samtools view -@ 5 -Sb - | \
samtools sort -@ 5 -o $samp.sorted.bam -
samtools index -@ 5 $samp.sorted.bam

#############---------- Contig binning
module load MetaBAT/2.12.1

jgi_summarize_bam_contig_depths --outputDepth $samp.depth.txt $samp.sorted.bam
mkdir -p Bins
echo "running metabat2"
metabat2 -t 0 -m 1500 -a $samp.depth.txt -i $samp.contigs.filt.fa -o Bins/$samp.bin

