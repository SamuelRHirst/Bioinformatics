# HiFi Genome Assembly and Annotation
## Sam R. Hirst

***

This markdown document walks through the general pipeline followed for [PacBio HiFi](https://www.pacb.com/smrt-science/smrt-sequencing/hifi-reads-for-highly-accurate-long-read-sequencing/) genome assembly and annotation. I use my genome for *Crotalus ruber* (Red Diamond Rattlesnake CLP2635) as an example 

# Quality Check and Stats
We can use [NanoPlot](https://github.com/wdecoster/NanoPlot) to check out some general statisitics on our long-read data such as the average/median read lengths and N50 before assembly.
```
NanoPlot -t 23 -o NanoPlot_res --fasta demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_reads.fasta  demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_reads.fasta
```
# Assembly
After testing other methodologies, I chose to use [hifiasm](https://hifiasm.readthedocs.io/en/latest/index.html) as the primary assembler as it was both the fastest and produced the best assembly.
-t is number of CPUs, -o specifies prefix of output file. hifiasm outputs the primary assembly in `gfa` format, we can convert this to fasta format using `awk`
```
hifiasm -o Cruber_CLP2635_v2.asm -t23 ~/Run1/demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_reads.fasta ~/Run2/demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_reads.fasta
awk '$1 ~/S/ {print ">"$2"\n"$3}' Cruber_CLP2635_v2.asm.p_ctg.gfa > Cruber_CLP2635_v2.asm.fasta
```
## OPTIONAL: Additional Purge_dups
You can increase genome assembly quality by running Purge_dups post assembly.  HiFiasm already does Purge_dups, but it apparently is not as quality as running it independently

```
#set the target assembly and reads to variables so I don't have to type them
#repeatedly
ASM="~/Cruber/CLP2635/Assembly/Run3/Cruber_CLP2635_v2.asm.fasta"
PBDATA="~/Cruber/CLP2635/Raw/Run2/demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_combined_reads.fasta"


#generate paf, read depth histogram, base-level read depth
conda activate
minimap2 -xmap-pb $ASM $PBDATA | gzip -c - > ${ASM}.paf.gz
conda deactivate

/home/h/hirsts/purge_dups/bin/pbcstat ${ASM}.paf.gz #(produces PB.base.cov and PB.stat files)

/home/h/hirsts/purge_dups/bin/calcuts PB.stat > cutoffs 2>calcuts.log

#split assembly and self-align
/home/h/hirsts/purge_dups/bin/split_fa $ASM > ${ASM}.split

conda activate
minimap2 -xasm5 -DP ${ASM}.split ${ASM}.split | gzip -c - > ${ASM}.split.self.paf.gz
conda deactivate

#purge haplotigs and overlaps
/home/h/hirsts/purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov ${ASM}.split.self.paf.gz > dups.bed 2> purge_dups.log

#Get purged primary and haplotig sequences from draft assembly
/home/h/hirsts/purge_dups/bin/get_seqs -e dups.bed $ASM
```
Final file name will be purged.fa.  I eventually change this to Cruber_CLP2635_v2.purge_dups.asm.fasta
# General Assembly Stats
I used [PyPi](https://pypi.org/project/assembly-stats/)

Will show on both initial Hifiasm assembly and the assembly with additional Purge_dups run
```
#Run on hifiasm only assembly
source ~/.bashrc
python /home/h/hirsts/compute/assembly_stats/assembly_stats.py ~/Cruber/CLP2635/Assembly/Run3/Cruber_CLP2635_v2.asm.fasta
```
  "Contig Stats": 
    "L10": 6,
    "L20": 16,
    "L30": 30,
    "L40": 45,
    "L50": 67,
    "N10": 20420946,
    "N20": 13483913,
    "N30": 11415746,
    "N40": 8893022,
    "N50": 6076602,
    "gc_content": 39.816707174710174,
    "longest": 29977177,
    "mean": 1122078.5326388888,
    "median": 174182.5,
    "sequence_count": 1440,
    "shortest": 13316,
```
#Run on post-purge_dups assembly
source ~/.bashrc
python /home/h/hirsts/compute/assembly_stats/assembly_stats.py ~/Cruber/CLP2635/Assembly/RunPurge/purged.fa
```
  "Contig Stats": 
  "L10": 6,
  "L20": 16,
  "L30": 29,
  "L40": 44,
  "L50": 65,
  "N10": 20420946,
  "N20": 13483913,
  "N30": 11437693,
  "N40": 9038368,
  "N50": 6254503,
  "gc_content": 39.77977003198217,
  "longest": 29977177,
  "mean": 1408698.393428064,
  "median": 274293.5,
  "sequence_count": 1126,
  "shortest": 17104,
  "total_bps": 1586194391


From here on out, I only use the additional Purge_dups assembly
# BUSCO
To assess the quality of your genome assembly, you will want to run [BUSCO](https://busco.ezlab.org/). Ideally, you will have BUSCO completeness at least >90%
I use [BUSCO Lineages](https://busco.ezlab.org/list_of_lineages.html) vertebrata_odb10 and sauropsida_odb10
```
conda activate BUSCO

asm="/shares_bgfs/margres_lab/Snakes/Genomes/Cruber/CLP2635/Assembly/RunPurge/purged.fa"

busco -i $asm -l sauropsida_odb10 -c 16 -m genome -r --out busco_purge_sauro
busco -i $asm -l vertebrata_odb10 -c 16 -m genome -r --out busco_purge_vert
```
 
 ## Additional Stats
 This is based on the [Proposed standards and metrics for defining genome assembly quality](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8081667/) from Rhie et al. 2021. Nature
 
Start with Phred based Accuracy (QV) and K-mer Completeness using [MERQURY](https://github.com/marbl/merqury)
 ```
#Start by creating Meryl counts of the high quality PacBio HiFi data
conda activate merqury

asm=~/Cruber/CLP2635/Assembly/RunPurge/Cruber_CLP2635_v2.purge_dups.asm.fasta
hapasm=~/Cruber/CLP2635/Assembly/RunPurge/hap.fa

#create meryl database for Merqury#
meryl count k=21 /shares_bgfs/margres_lab/Snakes/Genomes/Cruber/CLP2635/Raw/Run2/demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_combined_reads.fasta output pacbio.meryl

#Run Merqury using the meryl database, your primary purge output and haplotig output
merqury.sh pacbio.meryl $asm $hapasm merqury_pacbio
 
##QV ~/Cruber/CLP2635/Assembly/RunPurge/Assembly_stats/kmer/merqury_pacbio.qv contains QV value ##
##Kmer Completness ~/Cruber/CLP2635/Assembly/RunPurge/Assembly_stats/kmer contains K-mer completeness ##
 ```
False Duplication %
```
conda activate merqury

$MERQURY/eval/false_duplications.sh merqury_pacbio.Cruber_CLP2635_v2.purge_dups.asm.spectra-cn.hist
```
        

Genome Tools for Estimating Genome Size and Calculating NG50
You will need [Genometools](http://genometools.org/) [Genomescope](http://qb.cshl.edu/genomescope/) and [Jellyfish](https://github.com/gmarcais/Jellyfish)
They also are in Conda
```
###Genome Tools for Estimating Genome Size and Calculating NG50###
#Jellyfish and genomescope and genometools are installed in their own conda environments. Make sure you get kmer-jellyfish and genometools-genometools as well within their respective conda environments

READS=~/Cruber/CLP2635/Raw/Run2/demultiplex.bc1012_BAK8A_OA--bc1012_BAK8A_OA.hifi_combined_reads.fasta
ASM=~/Cruber/CLP2635/Assembly/RunPurge/Cruber_CLP2635_v2.purge_dups.asm.fasta

#compute k-mer distribution with jellyfish
conda activate jellyfish
jellyfish count -C -m 21 -s 1000000000 -t 20 $READS -o reads.jf
jellyfish histo -t 20 reads.jf > reads.histo
conda deactivate

#estimate genome size with GenomeScope2
#https://github.com/tbenavi1/genomescope2.0
conda activate genomescope
genomescope2 -i reads.histo -o genomescope_reads -k 21 > genomescope_out
conda deactivate

SIZE=$(sed -n -e 's/^.*len://p' genomescope_out)

#Compute NG50 (and other stats)
conda activate genometools
gt seqstat -contigs -genome $SIZE $ASM > ${ASM}_assembly.stats
conda deactivate
```
# Annotation
If you have reference genomes of closely-related taxa, I highly reccommend using [GeMoMa](http://www.jstacs.de/index.php/GeMoMa). Otherwise, you may want to try other annotators such as MAKER, BRAKER2, or Funannotate
## Repeat Masking
This is needed for most annotators. Hard masking turns repeat elements into N. Soft masking turns repeat elements to lower case letters.
I used [Repeatmodeler](https://www.repeatmasker.org/)
```
#Repeatmodeler has its own conda environment - also part of repeatmodeler
#xsmall will enable softmasking 

#Start with RepeatModeler, which will create a repeat database of your assembly

###-UPDATE YOUR INPUTS
INPUT=/path/to/your/assembly/.fasta
DB=DB_NAME

###-SETUP DETECTION OF REPEATS
#echo "Build database ..."
BuildDatabase -name $DB $INPUT

echo "Run RepeatModeler ..."
RepeatModeler -pa 50 -database $DB > out.log

ln -s RM_*/consensi.fa.classified ./

#RepeatMasker will now actually modify your FASTA file for repeats. Annotators prefer different types of masking - hard or soft. GEMOMA prefers soft

echo "Run RepeatMasker ..."
RepeatMasker -pa 50 -xsmall -gff -lib consensi.fa.classified -dir MaskerOutput $INPUT

```
## Gene Annotation
I used [GeMoMa](http://www.jstacs.de/index.php/GeMoMa) which was reccommended to me by a geneticist with the Smithsonian Data Science Lab (Dr. Paul Frandsen)
It relies heavily on references, so having closely related taxa with good reference genomes is necesarry
It can also use Transcriptomes as reference (Other tutorial available for Transcriptome data processing)
```
#GEMOMA has its own conda envrionment#
#References can be obtained from GenBank

#References
acgff=~/Reference/Acarol/GCF_000090745.1_AnoCar2.0_genomic_edited.gff
acfna=~/Reference/Acarol/GCF_000090745.1_AnoCar2.0_genomic.fna
cagff=~/Reference/Cadam/MASTER_Cadam_genes_cds_mRNA_30Jun21_edited.gff3
cafna=~/Reference/Cadam/Cadamanteus_3dDNAHiC_1.2.fasta
ctgff=~/Reference/Ctigris/GCF_016545835.1_ASM1654583v1_genomic_edited.gff
ctfna=~/Reference/Ctigris/GCF_016545835.1_ASM1654583v1_genomic.fna
tegff=~/Reference/Telegans/GCF_009769535.1_rThaEle1.pri_genomic_edited.gff
tefna=~/Reference/Telegans/GCF_009769535.1_rThaEle1.pri_genomic.fna

genome=~/Cruber/CLP2635/Assembly/RunPurge/Cruber_CLP2635_v2.purge_dups.asm.fasta

# -Xmx is max memory and is needed for many java programs. 
# All references were obtained from NCBI. RNA seq data from BAM files is also incoorperated with r=mapped because they are BAM files

conda activate gemoma
GeMoMa GeMoMaPipeline threads=50 -Xmx250G outdir=final_annot GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=true t=$genome s=own i=acarol a=$acgff g=$acfna s=own i=cadam a=$cagff g=$cafna s=own i=ctigris a=$ctgff g=$ctfna s=own i=telegans a=$tegff g=$tefna r=MAPPED ERE.m=$blood ERE.m=$heart ERE.m=$gonad ERE.m=$kideny ERE.m=$liver ERE.m=$VG
```
Next, we can add more useful annotations using [Interproscan](https://www.ebi.ac.uk/interpro/search/sequence/) and [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

You first need to obtain a Uniprot protein database of reference organisms for Blast. I chose Acarol Cadam Ctigris Ohannah Telegans and combined their uniprot files into a single FASTA file combined_uniprot.fasta

Interproscan
```
module load apps/jdk/15.0.2
module load apps/python/3.8.5

#The edited predicted proteins fasta just had to do with removing the * in the file.  Interproscan doesn't read * #

./interproscan.sh -i ~/Cruber/CLP2635/Annotation/GeMoMa/misc_out/predicted_proteins_edited.fasta -b output_GEMOMA/ -cpu 48
```
BLAST
```
conda activate blast

#start by making blastdb and specify protein
makeblastdb -in combined_uniprot.fasta -out Uniprot_combined -dbtype prot

#default out format is 6, which we want it to be for our next step. Also, we are using the protein FASTA output from the annotation of GeMoMa.  Some will add * as stop codons, but blastp doesn't like that, so be sure to edit your protein FASTA file to elminate the *

blastp -query ~/Cruber/CLP2635/Annotation/GeMoMa/misc_out/predicted_proteins_edited.fasta -db Uniprot_combined -max_target_seqs 2 -max_hsps 1 -evalue 1e-6 -outfmt 6 -out Cruber.outfmt6.blastp -num_threads 40
conda deactivate
```
Add BLAST and Interproscan outputs to annotation using [AGAT](https://github.com/NBISweden/AGAT)

```
##Add BLAST and Interproscan output using agat

Annotation=~/Cruber/CLP2635/Annotation/GeMoMa/misc_out/Cruber_annotation_GeMoMa.gff
Inter=~/Cruber/CLP2635/Annotation/interproscan/interproscan-5.52-86.0/output_GEMOMA/predicted_proteins_edited.fasta.tsv
blast=~/Cruber/CLP2635/Annotation/blast/UNIPROT/Cruber.outfmt6.blastp
blastdb=~/Cruber/CLP2635/Annotation/blast/UNIPROT/combined_uniprot.fasta

conda activate agat

agat_sp_manage_functional_annotation.pl -f $Annotation -i $Inter -b $blast -d $blastdb -o Functional

# This combines the information from interproscan and blast into your annotation GFF to give functional information, including mRNA product info
```



