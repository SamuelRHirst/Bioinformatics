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
I used [GeMoMa](http://www.jstacs.de/index.php/GeMoMa)
It relies heavily on references, so having closely related taxa with good reference genomes is necesarry
It can also use Transcriptomes as reference (Other tutorial available for Transcriptome data processing)
```
#GEMOMA has its own conda envrionment#
#References can be obtained from GenBank

#References
#Anolis Carolinensis
acgff=~/Reference/Acarol/GCF_000090745.1_AnoCar2.0_genomic_edited.gff
acfna=~/Reference/Acarol/GCF_000090745.1_AnoCar2.0_genomic.fna
#Crotalus adamanteus
cagff=~/Reference/Cadam/MASTER_Cadam_genes_cds_mRNA_30Jun21_edited.gff3
cafna=~/Reference/Cadam/Cadamanteus_3dDNAHiC_1.2.fasta
#Crotalus tigris
ctgff=~/Reference/Ctigris/GCF_016545835.1_ASM1654583v1_genomic_edited.gff
ctfna=~/Reference/Ctigris/GCF_016545835.1_ASM1654583v1_genomic.fna
#Thamnophis elegans
tegff=~/Reference/Telegans/GCF_009769535.1_rThaEle1.pri_genomic_edited.gff
tefna=~/Reference/Telegans/GCF_009769535.1_rThaEle1.pri_genomic.fna

#Reference Genome
genome=~/Cruber/CLP2635/Assembly/RunPurge/Cruber_CLP2635_v2.purge_dups.asm.fasta

# -Xmx is max memory and is needed for many java programs. 
# All references were obtained from NCBI. RNA seq data from BAM files is also incoorperated with r=mapped because they are BAM files

#I reccommend the following flags. For o=true, you will have an reference output for each reference, o=false will be a single output
conda activate gemoma
GeMoMa GeMoMaPipeline threads=50 -Xmx250G outdir=final_annot GeMoMa.Score=ReAlign AnnotationFinalizer.r=NO o=false t=$genome s=own i=acarol a=$acgff g=$acfna s=own i=cadam a=$cagff g=$cafna s=own i=ctigris a=$ctgff g=$ctfna s=own i=telegans a=$tegff g=$tefna r=MAPPED ERE.m=$blood ERE.m=$heart ERE.m=$gonad ERE.m=$kideny ERE.m=$liver ERE.m=$VG
```
Next, we can add more useful annotations using [Interproscan](https://www.ebi.ac.uk/interpro/search/sequence/) and [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

Note on April 3 2023, I tried rerunning this part, but Interproscan was giving me issues so I just used BLAST. That is likely sufficient for our goal here. 

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
