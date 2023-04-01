# Transcriptomic Pipeline
## Sam R. Hirst

***

# Quality check and trimming
We can use [TrimGalore](https://github.com/FelixKrueger/TrimGalore) to produce some QC files and trim the raw files
```
#Trim Galore has its own conda envrionment 

F1BLOOD=/path/to/raw/RNA/fileF1
R1BLOOD=/path/to/raw/RNA/fileR1
F2BLOOD=/path/to/raw/RNA/fileF2
R2BLOOD=/path/to/raw/RNA/fileR2

OUTPUT=/shares_bgfs/margres_lab/Snakes/Transcriptomes/Cruber/CLP2635-Blood/trimmed-Blood

start=`date +%s`

source ~/.bashrc

conda activate trim-galore


#Trimming and FastQC
trim_galore --paired $F1BLOOD $R1BLOOD $F2BLOOD $R2BLOOD -o $OUTPUT

end=`date +%s`
printf "Execution time: $((end-start))s"
```
# Mapping
Next, we need to map our trimmed files to our reference genome
For this, we will use [Hisat2](https://github.com/DaehwanKimLab/hisat2)
```

F1BLOOD=/path/to/trimmed/RNA/fileF1
R1BLOOD=/path/to/trimmed/RNA/fileR1
F2BLOOD=/path/to/trimmed/RNA/fileF2
R2BLOOD=/path/to/trimmed/RNA/fileR2

source ~/.bashrc
conda activate hisat2

#Start by building DB for Hisat2 - can use same DB for all transcriptomes of an individual
hisat2-build /shares_bgfs/margres_lab/Snakes/Genomes/Cruber/CLP2635/Assembly/RunPurge/Cruber_CLP2635_v2.purge_dups.asm.fasta genome

#Now, you can map your transcriptomes to that DB, producing sam files
hisat2 -x genome -1 $F1 $F2 -2 $R1 $R2 -S hisat2_blood.sam

#Bam files are better
module load apps/samtools/1.3.1

samtools view -Sb hisat2_blood.sam > hisat2_blood.bam
#Sort bam file
samtools sort -n hisat2_blood.bam -o hisat2_blood_sorted.bam
```
# Mapping to estimate expression
This also uses hisat2, but the process is a bit different
You will need the reference genome assembly and annotation IN GTF

```
#Build index first
#Building genome index with annotation - First generate ss and exon files from gtf
#############

hisat2_extract_splice_sites.py $ANNO > Cruber_genome.ss

hisat2_extract_exons.py $ANNO > Cruber_Genome.exon

#Now generate annotated index

hisat2-build -p 50 --ss Cruber_genome.ss --exon Cruber_Genome.exon -f $GENOME Cruber_genome_Anno

#Now we can run hisat2:
hisat2 -p 50 -k 10 --dta -x Cruber_genome_Anno -1 $CLP2635VGF -2 $CLP2635VGR -S /path_to_output/CLP2635_VenomGland.sam

#convert to Sam and Index

module load apps/samtools/1.3.1

samtools sort -@ 50 -o /Path_to_Output/Cruber_CLP2635_VG.bam /Path_to_Input/SAM/Cruber_CLP2635_VG.sam
samtools index /shares_bgfs/margres_lab/Snakes/Transcriptomes/Cruber/Analyses/Hisat2/BAM/Cruber_CLP2635_VG.bam
```
# Estimating expression
This is done using Stringtie

```
stringtie -p 50 -G $ANNO -e -o /Path_to_output/CLP2635_VenomGland.gtf -B -A /Path_to_output/CLP2635_VenomGland_abund.tab /Path_to/BAM/CLP2635_VenomGland.bam

```
The abund.tab file is most relevant. Final column shows transcript per million (TPM) for each annotated gene. 
