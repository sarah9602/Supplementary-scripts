# Assembly pipeline
Process reads prior to running pipeline
Remove adapters, remove PCR duplicates, filter out host reads

Establish variables for where your samples are located
```bash
sampdir=/PATH/TO/PROCESSED/READS
```

## ------- Assembly - submit each sample through assembly pipeline -------
See [assembly_pipeline_template.sh](https://github.com/sarah9602/Supplementary-scripts/blob/main/assembly_pipeline_template.sh) using HPC with SLURM queing system
Process is as follows:
**Megahit**
-k-min 21 
-k-max 91 
-k-step 10
-memory 0.9 
-t 8 
### --- Filter out contigs < 1000 nt long with seqkit
### --- Readmapping with Bowtie2
1. bowtie2-build contigs.fa file
2. Align processed short fastq reads to the bt2 contig database for the corresponding sample
3. Feed the ouput to samtools to conver to bam, sort, and index the alignment file.
### --- Metabat binning
jgi_summarize_bam_contig_depths to output depth file
bin with -m 1500 (This is the minimum allowed amount)

### Read list of files to feed through assembly_pipeline_template.sh
The output will be individual assembly pipeline files for each processed metagnome.
```bash
ls $sampdir/*PROCESSED.fastq.gz| while read infile; do samp=`echo $infile|awk -F"/" '{print $NF}'|awk -F"." '{print $1}'`; sed -e "s|INFILE|$infile|g" -e "s/SAMPLE/$samp/g" assembly_pipeline_template.sh > $samp.assembly.sh; done
```
Then run all the individual pipeline files so they run on the HPC simultaneously.
```bash
ls *assembly.sh|while read in; do sbatch $in; done
```

### Move to resulting output from working directory to output directory
```bash
mkdir -p $sampdir/Assembly/Binning $sampdir/Assembly/Readmapping $sampdir/Assembly/Megahit
mv /scratch/sjj9602/Assembly_*/Bins/* $sampdir/Assembly/Binning
mv /scratch/sjj9602/Assembly_*/*depth.txt $sampdir/Assembly/Binning
mv /scratch/sjj9602/Assembly_*/*sorted.bam* $sampdir/Assembly/Readmapping
mv /scratch/sjj9602/Assembly_*/*contigs.filt.fa $sampdir/Assembly/Megahit
```

## ------- Assembly Stats -------
```bash
bbmap/statswrapper.sh $sampdir/Assembly/Megahit/*contigs.filt.fa format=5 > $sampdir/Assembly/post_filt_assembly_stats.txt
```

############################################
## Assembly quality check
############################################

### ------- CheckM -------
This doesn't require high performance computing. Binning is output directory that holds contig bins with .fa extension.CheckM_output is created upon command.
**See [checkm](https://github.com/Ecogenomics/CheckM?tab=readme-ov-file) for installation instructions.**

**OR Load conda environment**
[checkm.yml](https://github.com/sarah9602/Supplementary-scripts/checkm.yml)
```bash
conda env create --name checkm checkm.yml
```

**Run checkm**
```bash
path=$sampdir/Assembly
conda activate checkm
## This step takes a long time
checkm lineage_wf -f checkm_output.txt -t 20 -x fa ./Binning ./CheckM_output 2>&1 | tee checkm.lineage.log
## Depending on how many bins, this takes ~ 5 mins or less
checkm qa -o 2 -t 20 -f checkm_stats --tab_table CheckM_output/lineage.ms CheckM_output 2>&1 | tee checkm.qa.log

## Pasoli parameters
## Completeness
## Medium quality - 50% <= x <= 90%
## High quality - > 90%
## All < 5% contamination
cat checkm_stats|awk -F"\t" '$6 >= "50"'|awk -F"\t" '$6 < "90"'|awk -F"\t" '$7 <= "5"'|awk -F"\t" '{print $1}' > med_qual.txt
cat checkm_stats|awk -F"\t" '$6 >= "90"'|awk -F"\t" '$7 <= "5"'|awk -F"\t" '{print $1}' > high_qual.txt 

mkdir HighMediumQuality/
cat med_qual.txt|while read in; do mv Binning/$in.fa HighMediumQuality/; done
cat high_qual.txt|while read in; do mv Binning/$in.fa HighMediumQuality/; done
```
############################################
## Assign Taxonomy
############################################
Can be with phylophlan or kraken2 - have not determined best route yet

**Install [PhyloPhlAn](https://github.com/biobakery/phylophlan)]**

## ------- Phylophlan -------
phylophlan_assign_sgbs -i MediumQuality/ -d /media/projects/References/phylophlan_databases/SGB.Jan21 -n 1 --nproc 20
phylophlan_assign_sgbs -i HighQuality/ -d /media/projects/References/phylophlan_databases/SGB.Jan21 -n 1 --nproc 20
## From Grunhall et al - limit MASH distance to < 5%
phylophlan_draw_metagenomic -i ethiopian_mags/ethiopian_mags.tsv --map tutorial_ethiopia__mag2meta.tsv -f png --verbose 2>&1 | tee phylophlan_draw_metagenomic.log
phylophlan_draw_metagenomic -i chomper_mags.1.tsv --map chomper_mag2meta.tsv -f pdf --verbose 2>&1|tee phylophlan_draw_metagenomic.log

## ------- Kraken2 -------
**See [kraken.sbatch](https://github.com/sarah9602/Anson-St-Ancestors-Oral-Microbiome/blob/main/Kraken-Analysis/run_kraken-bracken.sh) script**


## ------- Damage ------- Optional - have not tested
conda activate damage
path=/media/projects/johnson/Assembly/Binninb
ls 
mapDamage -i $path/ReadMapping/$sample.sorted.bam -r $path/Binning/$sample.bins/High/$sample.bin.13.fa -d $sample.bin.13.fa.MapDamage_results --rescale --rescale-out=$sample.mapped_sorted_rmdup.rescaled.bam]

**See [pyDamage](https://github.com/maxibor/pydamage) also**




