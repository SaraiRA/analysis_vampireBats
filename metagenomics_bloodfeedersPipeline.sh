#!/bin/bash

# CALL MODULES
#module avail 
module load AdapterRemoval/v2.2.0 
module load perl/v5.28.1
module load prinseq/v0.20.4 
module load bwa/v0.7.15
module load htslib/v1.6
module load samtools/v1.6 
module load bedtools/v20170719
module load jellyfish/v1.1.11 
module load java/v11.0.1-jdk
module load fastqc/v0.11.8a
module load bcftools/v1.4
module load vcftools/v0.1.14

module load kraken/v2.0.7 
module load python/v2.7.12
module load bowtie2/v2.2.9
module load muscle/v3.8.31
module load blast+/v2.2.28
module load samtools/v0.1.19
module load RAxML/v8.2.11
module load diamond/v0.8.31
module load MetaPhlAn/v2.6.0


## Define the variables
PROJECT=/groups/hologenomics/sarai/data/batProject
RAWREADS=$PROJECT/rawReads
FASTQ=$PROJECT/fastqs
ADAPTER=$PROJECT/adapRem
NOHOST=$PROJECT/noHost
NOHOSTCOPY=$PROJECT/noHost_copy
MAP=$PROJECT/mapping
PREY=$PROJECT/preys
UNPREY=$PROJECT/unpreys
KRAKEN=$PROJECT/kraken2
TAPIRGENOME=$PROJECT/tapir_genome

# references
DIREF=/groups/hologenomics/data/genomes
BATGENOME=/groups/hologenomics/sarai/data/my_referenceGenomes/vampire/vampire_bat_10Aug2015_WXtkA.fasta
COWGENOME=${DIREF}/cow/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fasta
PIGGENOME=${DIREF}/pig/susScr11.fasta
SHEEPGENOME=${DIREF}/sheep/oviAri4.fasta
DONKEYGENOME=${DIREF}/donkey/GCF_001305755.1_ASM130575v1_genomic.fasta
HORSEGENOME=${DIREF}/horse/equCab2.masked.fasta
CHICKENGENOME=${DIREF}/chicken/galGal5.fasta
RHINOGENOME=${DIREF}/white_rhino/cerSim1.fasta
KRAKEN2DB="/groups/hologenomics/data/db/Kraken2-standard-DB/"


############### MERGING ###################
#Merging files of the same lane 
#Dont run with the rest of the pipeline

echo "Merge file of the same line"

## LINE
for f in $(ls *_L006_R1_001.fastq.gz)
do bn=$(basename $f _L006_R1_001.fastq.gz); cat "$bn"_L006_R2_0*.fastq.gz > "$bn"_L006_R2_ME.fastq.gz
done


############# QUALITY ####################
# fastqc for quality comparision after clean the data
mkdir fastqs_fastqc
fastqc *


##### ADAPTERS #####

## Remove the illumina adapters P5 and P7 (AdapterRemoval). 
# Yeah, here we are removing the adapters

#INPUT: fastq files

echo "Remove adapters"
# -p: no error if existing 
mkdir -p $ADAPTER && cd $ADAPTER

## LINE
if [ ! -e .adap.done ]; then
	# For each fastq pair, do cutadapt
	for f in $FASTQ/*_L006_R1_ME.fastq.gz
		# Figure out the name of the sample.
		do
		bn=$(basename $f _L006_R1_ME.fastq.gz)
    		
		## Run cutadapt to verify that all the sequences have the TE target. 
    		# -a : check for sequence on the 3' end of read1. It is N, so any base will do.
		# -G : check for sequence on the 5' end of read2. It is the primer for LINE, ^ means it should be in the star. 
  		# --no-trim: do not remove the primer sequences.
    		# --no-indels: do not account for short insertions and deletions.
    		# -e : mismatch rate, here it is 5%, 1 mismatch in the primer sequence.
    		# -o : output name for read 1
  		# -p : output name for read 2
    		# $f and ${f/R1/R2}: read 1 and read 2 input files respectively.
   
 		## Run Adapter removal upon the successful completion of cutadapt
		# collapse: Combined into a single consensus sequence pair aligments
		# Output: output_paired.collapsed containing merged reads,
		# output_paired.collapsed.truncated containing merged reads that have been trimmed due to the --trimns or --trimqualities options.
		# The sequencing primers are not specified, since AdapterRemoval knows the Illumina sequencing primers.
		echo "AdapterRemoval  --qualitybase 33 --qualitybase-output 33 --qualitymax 45 --gzip --mm 3 --minlength 25 --trimns --trimqualities --minquality 10 \
			--collapse --basename ${bn}_primer_noAdap --file1 $f --file2 ${f/R1/R2}"
  	done | xsbatch -c 1 --mem-per-cpu=2G -R -J adap --
 	touch .adap.done
fi

#Wait until the jobs above finish 
echo "Waiting"
read dummy



############### REMOVE DUPLICATES #################
# Remove exact duplicates 

echo "Remove duplicated reads" 
cd $ADAPTER

### COLLAPSED READS 
#INPUT: *_primer_noAdap.collapsed.gz and  *_primer_noAdap.collapsed.truncated.gz

if [ ! -e .RD.done ]; then
	## Run prinseq
	# -derep: Type of duplicates to filter. 1 (exact duplicate)
	# -derep_min: This option specifies the number of allowed duplicates. If you want to remove sequence duplicates that occur more than x times
	# -out_format 3, fastq and names
	# uncompress gz file for prinseq
  	for f in $ADAPTER/*_noAdap.collapsed.gz
  		do
    		bn=$(basename $f _noAdap.collapsed.gz)
    		echo "( gzip -dc <(cat $f $ADAPTER/${bn}_noAdap.collapsed.truncated.gz) |  prinseq-lite -verbose -fastq stdin -derep 1 -derep_min 2 -out_format 3 -out_good ${bn}_noAdapCollapsed_noDupPrinSeq -out_bad ${bn}_noAdapCollapsed_DupPrinSeq; gzip $ADAPTER/${bn}_noAdapCollapsed_noDupPrinSeq.fastq $ADAPTER/${bn}_noAdapCollapsed_DupPrinSeq.fastq)"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J prinseq -R --max-array-jobs=10 --
  touch .RD.done
fi

### UNCOLLAPSED READS 
#INPUT: *_primer_noAdap.collapsed.gz and  *_primer_noAdap.collapsed.truncated.gz

cd $ADAPTER
if [ ! -e .RDun.done ]; then
	## Run prinseq
	# -derep: Type of duplicates to filter. 1 (exact duplicate)
	# -derep_min: This option specifies the number of allowed duplicates. If you want to remove sequence duplicates that occur more than x times
	# -out_format 3, fastq and names
  	for f in $ADAPTER/*_noAdap.pair1.truncated.gz
  		do
    		bn=$(basename $f _noAdap.pair1.truncated.gz)
    		echo "(gunzip $f $ADAPTER/${bn}_noAdap.pair2.truncated.gz; prinseq-lite -fastq $ADAPTER/${bn}_noAdap.pair1.truncated -fastq2 $ADAPTER/${bn}_noAdap.pair2.truncated -derep 1 -derep_min 2 -out_format 3 -out_good ${bn}_noAdapUncollapsed_noDupPrinSeq -out_bad ${bn}_noAdapUncollapsed_DupPrinSeq)"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J prinseq_un -R --max-array-jobs=10 --
  touch .RDun.done
fi

#Wait until the jobs above finish 
echo "Waiting"
read dummy

gzip *noAdapUncollapsed*

#Wait until the jobs above finish 
echo "Waiting"
read dummy


################ REMOVE HOST ###############

# Map to the host reference genomes using BWA mem. 
# Use mem so that we can soft clip the primer sequences in the beginning of the read.

# Sort the mapping by coordinates using 
echo "Map to reference genome"
# -p: no error if existing 
mkdir -p $NOHOST && cd $NOHOST

### MAPPING ####
### COLLAPSED READS 

#INPUT: *_noAdapCollapsed_noDupPrinSeq

cd $NOHOST
if [ ! -e .map.done ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
	for f in $ADAPTER/*_noAdapCollapsed_noDupPrinSeq.fastq.gz  
		do     
		bn=$(basename $f _noAdapCollapsed_noDupPrinSeq.fastq.gz)
		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
		echo "(bwa mem -M $BATGENOME $f | samtools view -b -| samtools sort - -o ${bn}_MapBat_collapsed.markdup.bam ; samtools view -b -f 4 $NOHOST/${bn}_MapBat_collapsed.markdup.bam > $NOHOST/${bn}_UNMap_Bat_collapsed.markdup.bam)" 
	done | xsbatch -c 1 --mem-per-cpu=10G -J map -R --max-array-jobs=10 --   
	touch .map.done
fi

### UNCOLLAPSED READS 

if [ ! -e .UNcollapsed.map.done ]; then   
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
	for f in $ADAPTER/*_primer_noAdapUncollapsed_noDupPrinSeq_1.fastq 
		do     
		bn=$(basename $f _primer_noAdapUncollapsed_noDupPrinSeq_1.fastq)     
		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
		echo "(bwa mem -M $BATGENOME $f $ADAPTER/${bn}_primer_noAdapUncollapsed_noDupPrinSeq_2.fastq | samtools view -b -| samtools sort - -o ${bn}_MapBat_uncollapsed.markdup.bam ; samtools view -b -f 4 $NOHOST/${bn}_MapBat_uncollapsed.markdup.bam > $NOHOST/${bn}_UNMap_Bat_uncollapsed.markdup.bam)" 
	done | xsbatch -c 1 --mem-per-cpu=10G -J Unmap -R --max-array-jobs=10 --   
	touch .UNcollapsed.map.done
fi


#Wait until the jobs above finish 
echo "Waiting"
read dummy

### CONVERTING BAM TO FASTQ ###

### COLLAPSED READS 
cd $NOHOST
if [ ! -e .recoverfq.done ]; then for f in $NOHOST/*_UNMap_Bat_collapsed.markdup.bam; do bn=$(basename $f _UNMap_Bat_collapsed.markdup.bam); echo "(bedtools bamtofastq -i $f -fq $NOHOST/${bn}_UNMap_Bat_collapsed.markdup.fastq; gzip $NOHOST/${bn}_UNMap_Bat_collapsed.markdup.fastq )"; done | xsbatch -c 1 --mem-per-cpu=5G -J fq -R --max-array-jobs=10 --; touch .recoverfq.done; fi

### UNCOLLAPSED READS 
cd $NOHOST
if [ ! -e .ucollap.recoverfq.done ]; then for f in $NOHOST/*_UNMap_Bat_uncollapsed.markdup.bam; do bn=$(basename $f _UNMap_Bat_uncollapsed.markdup.bam); echo "(bamToFastq -i $f -fq $NOHOST/${bn}_UNMap_Bat_uncollapsed.markdup.fastq; gzip $NOHOST/${bn}_UNMap_Bat_uncollapsed.markdup.fastq )"; done | xsbatch -c 1 --mem-per-cpu=5G -J uncollpased_fq -R --max-array-jobs=10 --; touch .uncollap.recoverfq.done; fi

#Wait until the jobs above finish 
echo "Waiting"
read dummy

###MERGING#### 
cd $NOHOST
#Merging files of the same lane 
#Dont run with the rest of the pipeline

echo "Merge file of the same line"

for f in $(ls *_UNMap_Bat_collapsed.markdup.fastq.gz)
do bn=$(basename $f _UNMap_Bat_collapsed.markdup.fastq.gz); cat $f ${bn}*fastq.gz > ${bn}_UNMap_BatME.markdup.fastq.gz
done

#Wait until the jobs above finish 
echo "Waiting"
read dummy

################ MAPPING AND GENOME ASSEMBLIES OF DESCRIBED PREYS #################
#Mapping the rest of the reads from host, to every prey separatly

echo "Assemble the preys"
# -p: no error if existing 
mkdir -p $PREY && cd $PREY


### COW ###
PREYCOW=$PREY/cow
mkdir -p $PREYCOW
cd $PREYCOW
if [ ! -e .map.cow ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $NOHOST/*_primer_UNMap_BatME.markdup.fastq.gz 
  		do
    		bn=$(basename $f _primer_noAdap.collapsed.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M $COWGENOME $f | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}_MAPcow_UNMAPbatME.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J mapCow -R --max-array-jobs=10 --
  touch .map.cow
fi

### PIG ###
PREYPIG=$PREY/pig
mkdir -p $PREYPIG
cd $PREYPIG
if [ ! -e .map.pig ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $NOHOST/*_primer_UNMap_BatME.markdup.fastq.gz 
  		do
    		bn=$(basename $f _primer_noAdap.collapsed.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M $PIGGENOME $f | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}_MAPpig_UNMAPbatME.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J mapPig -R --max-array-jobs=10 --
  touch .map.pig
fi

### SHEEP ###
PREYSHEEP=$PREY/sheep
mkdir -p $PREYSHEEP
cd $PREYSHEEP
if [ ! -e .map.sheep ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $NOHOST/*_primer_UNMap_BatME.markdup.fastq.gz 
  		do
    		bn=$(basename $f _primer_noAdap.collapsed.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M $SHEEPGENOME $f | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}_MAPsheep_UNMAPbatME.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J mapSheep -R --max-array-jobs=10 --
  touch .map.sheep
fi

### DONKEY ###
PREYDONKEY=$PREY/donkey
mkdir -p $PREYDONKEY
cd $PREYDONKEY
if [ ! -e .map.donkey ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $NOHOSTCOPY/*_primer_UNMap_BatME.markdup.fastq.gz 
  		do
    		bn=$(basename $f _primer_noAdap.collapsed.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M $DONKEYGENOME $f | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}_MAPdonkey_UNMAPbatME.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J mapDonkey -R --max-array-jobs=10 --
  touch .map.donkey
fi

### HORSE ###
PREYHORSE=$PREY/horse
mkdir -p $PREYHORSE
cd $PREYHORSE
if [ ! -e .map.horse ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $NOHOST/*_primer_UNMap_BatME.markdup.fastq.gz 
  		do
    		bn=$(basename $f _primer_noAdap.collapsed.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M $HORSEGENOME $f | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}_MAPhorse_UNMAPbatME.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J mapHorse -R --max-array-jobs=10 --
  touch .map.horse
fi

### CHICKEN ###
PREYCHICKEN=$PREY/chicken
mkdir -p $PREYCHICKEN
cd $PREYCHICKEN
if [ ! -e .map.chicken ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $NOHOST/*_primer_UNMap_BatME.markdup.fastq.gz 
  		do
    		bn=$(basename $f _primer_noAdap.collapsed.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M $CHICKENGENOME $f | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}_MAPchicken_UNMAPbatME.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J mapChicken -R --max-array-jobs=10 --
  touch .map.chicken
fi

### RHINO ###
PREYRHINO=$PREY/rhino
mkdir -p $PREYRHINO
cd $PREYRHINO
if [ ! -e .map.rhino ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $NOHOSTCOPY/*_primer_UNMap_BatME.markdup.fastq.gz 
  		do
    		bn=$(basename $f _primer_noAdap.collapsed.gz)
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M $RHINOGENOME $f | samtools sort -n -O bam - | samtools fixmate -r -p -m - - | samtools sort - | samtools markdup -r - ${bn}_MAPrhino_UNMAPbatME.markdup.bam )"
  	done | xsbatch -c 1 --mem-per-cpu=10G -J mapRhino -R --max-array-jobs=10 --
  touch .map.rhino
fi

################ GENOME ASSEMBLY OF TAPIR  #################
#Mapping the rest of the reads from host, to every prey by pipe

echo "Assembly tapir genome"
# -p: no error if existing 
mkdir -p $TAPIRGENOME && cd $TAPIRGENOME

#merge with samtools 
samtools merge tapir_8samples.bam $PREY/rhino/TOG_KLEY_116_CAGCTA_primer_UNMap_BatME.markdup.fastq.gz_MAPrhino_UNMAPbatME.markdup.bam $PREY/rhino/TOG_KLEY_121_GACGAC_primer_UNMap_BatME.markdup.fastq.gz_MAPrhino_UNMAPbatME.markdup.bam $PREY/rhino/TOG_KLEY_25_TGTGAC_primer_UNMap_BatME.markdup.fastq.gz_MAPrhino_UNMAPbatME.markdup.bam   $PREY/rhino/TOG_KLEY_29_GACACT_primer_UNMap_BatME.markdup.fastq.gz_MAPrhino_UNMAPbatME.markdup.bam   $PREY/rhino/TOG_KLEY_90_ACGCAT_primer_UNMap_BatME.markdup.fastq.gz_MAPrhino_UNMAPbatME.markdup.bam   $PREY/rhino/TOG_KLEY_94_ACATAC_primer_UNMap_BatME.markdup.fastq.gz_MAPrhino_UNMAPbatME.markdup.bam   $PREY/rhino/TOG_KLEY_54_TAGATG_primer_UNMap_BatME.markdup.fastq.gz_MAPrhino_UNMAPbatME.markdup.bam   $PREY/rhino/TOG_KLEY_54_TAGATG_primer_UNMap_BatME.markdup.fastq.gz_MAPrhino_UNMAPbatME.markdup.bam   

# create samtools index
samtools index tapir_8samples.bam

#Consensus sequences 
samtools mpileup -uf cerSim1 tapir_8samples.bam | bcftools view -cg - | bcftools/vcfutils.pl vcf2fq > consensus_tapir_8samples.bam.fasta


samtools mpileup -ABuf cerSim1.fasta tapir_8samples.bam  | bcftools call -cOz --pval-threshold 0.99 > mapped.vcf.gz
tabix mapped.vcf.gz 
cat reference.nt | bcftools consensus mapped.vcf.gz > mapped.fastq

################ RECOVERING UNMAPPED READS  #################
#Mapping the rest of the reads from host, to every prey by pipe

echo "Recovering unmapped reads"
# -p: no error if existing 
mkdir -p $UNPREY && cd $UNPREY;
cd $UNPREY;
if [ ! -e .map.done ]; then
	for f in $NOHOSTCOPY/*_primer_UNMap_BatME.markdup.fastq.gz 
  		do
    		bn=$(basename $f _primer_noAdap.collapsed.gz);
    		# Run bwa mem and then sort then sort the bam by coordinates.
		# M: mark shorter split hits as secondary
    		echo "(bwa mem -M $COWGENOME $f | samtools view -b -| samtools sort - -o ${bn}_MapCow_UNMapBat.markdup.bam ; samtools view -b -f 4 $UNPREY/${bn}_MapCow_UNMapBat.markdup.bam > $UNPREY/${bn}_Map_UNMapBatCow.markdup.bam; bedtools bamtofastq -i $UNPREY/${bn}_Map_UNMapBatCow.markdup.bam -fq $UNPREY/${bn}_Map_UNMapBatCow.markdup.fastq; gzip $UNPREY/${bn}_Map_UNMapBatCow.markdup.fastq; bwa mem -M $PIGGENOME $UNPREY/${bn}_Map_UNMapBatCow.markdup.fastq.gz | samtools view -b -| samtools sort - -o ${bn}_MapPig_UNMapBatCow.markdup.bam; samtools view -b -f 4 $UNPREY/${bn}_MapPig_UNMapBatCow.markdup.bam > $UNPREY/${bn}_Map_UNMapBatCowPig.markdup.bam; bedtools bamtofastq -i $UNPREY/${bn}_Map_UNMapBatCowPig.markdup.bam -fq $UNPREY/${bn}_Map_UNMapBatCowPig.markdup.fastq; gzip $UNPREY/${bn}_Map_UNMapBatCowPig.markdup.fastq; bwa mem -M $SHEEPGENOME $UNPREY/${bn}_Map_UNMapBatCowPig.markdup.fastq.gz | samtools view -b -| samtools sort - -o ${bn}_MapSheep_UNMapBatCowPig.markdup.bam; samtools view -b -f 4 $UNPREY/${bn}_MapSheep_UNMapBatCowPig.markdup.bam > $UNPREY/${bn}_Map_UNMapBatCowPigSheep.markdup.bam; bedtools bamtofastq -i $UNPREY/${bn}_Map_UNMapBatCowPigSheep.markdup.bam -fq $UNPREY/${bn}_Map_UNMapBatCowPigSheep.markdup.fastq; gzip  $UNPREY/${bn}_Map_UNMapBatCowPigSheep.markdup.fastq; bwa mem -M $DONKEYGENOME $UNPREY/${bn}_Map_UNMapBatCowPigSheep.markdup.fastq.gz | samtools view -b -| samtools sort - -o ${bn}_MapDonkey_UNMapBatCowPigSheep.markdup.bam; samtools view -b -f 4 $UNPREY/${bn}_MapDonkey_UNMapBatCowPigSheep.markdup.bam > $UNPREY/${bn}_Map_UNMapBatCowPigSheepDonkey.markdup.bam; bedtools bamtofastq -i $UNPREY/${bn}_Map_UNMapBatCowPigSheepDonkey.markdup.bam -fq $UNPREY/${bn}_Map_UNMapBatCowPigSheepDonkey.markdup.fastq; gzip $UNPREY/${bn}_Map_UNMapBatCowPigSheepDonkey.markdup.fastq; bwa mem -M $HORSEGENOME $UNPREY/${bn}_Map_UNMapBatCowPigSheepDonkey.markdup.fastq.gz | samtools view -b -| samtools sort - -o ${bn}_MapHorse_UNMapBatCowPigSheepDonkey.markdup.bam; samtools view -b -f 4 $UNPREY/${bn}_MapHorse_UNMapBatCowPigSheepDonkey.markdup.bam > $UNPREY/${bn}_Map_UNMapBatCowPigSheepDonkeyHorse.markdup.bam; bedtools bamtofastq -i $UNPREY/${bn}_Map_UNMapBatCowPigSheepDonkeyHorse.markdup.bam -fq $UNPREY/${bn}_Map_UNMapBatCowPigSheepDonkeyHorse.markdup.fastq; gzip  $UNPREY/${bn}_Map_UNMapBatCowPigSheepDonkeyHorse.markdup.fastq; bwa mem -M $CHICKENGENOME $UNPREY/${bn}_Map_UNMapBatCowPigSheepDonkeyHorse.markdup.fastq.gz | samtools view -b -| samtools sort - -o ${bn}_MapChicken_UNMapBatCowPigSheepDonkeyHorse.markdup.bam; samtools view -b -f 4 $UNPREY/${bn}_MapChicken_UNMapBatCowPigSheepDonkeyHorse.markdup.bam > $UNPREY/${bn}_Map_UNMapBatCowPigSheepDonkeyHorseChicken.markdup.bam; bedtools bamtofastq -i $UNPREY/${bn}_Map_UNMapBatCowPigSheepDonkeyHorseChicken.markdup.bam -fq $UNPREY/${bn}_Map_UNMapBatCowPigSheepDonkeyHorseChicken.markdup.fastq; gzip  $UNPREY/${bn}_Map_UNMapBatCowPigSheepDonkeyHorseChicken.markdup.fastq; bwa mem -M $RHINOGENOME $UNPREY/${bn}_Map_UNMapBatCowPigSheepDonkeyHorseChicken.markdup.fastq.gz | samtools view -b -| samtools sort - -o ${bn}_MapRhino_UNMapBatCowPigSheepDonkeyHorseChicken.markdup.bam; samtools view -b -f 4 $UNPREY/${bn}_MapRhino_UNMapBatCowPigSheepDonkeyHorseChicken.markdup.bam > $UNPREY/${bn}_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.bam; bedtools bamtofastq -i $UNPREY/${bn}_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.bam -fq $UNPREY/${bn}_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq; gzip  $UNPREY/${bn}_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq)";
  	done | xsbatch -c 1 --mem-per-cpu=70G -J recover -R --max-array-jobs=10 --
  touch .map.done
fi

################ IDENTITY UNMAPPED READS  #################
#Mapping the rest of the reads from host, to every prey by pipe

##### KRAKEN #####
echo "Running kraken"
# -p: no error if existing 
mkdir -p $KRAKEN2 && cd $KRAKEN2


if [ ! -e .kraken ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $UNPREY/*_primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq.gz
  		do
    		bn=$(basename $f _primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq.gz)
    		# Run kraken, report, mpa report
    		echo "(gunzip $f; kraken2 --db $KRAKEN2DB --threads 8 --preload  --output $KRAKEN2/${bn}.kraken $UNPREY/${bn}_primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq; kraken2-report --db $KRAKEN2DB  $KRAKEN2/${bn}.kraken  >  $KRAKEN2/${bn}.kraken.report; gzip  $UNPREY/${bn}_primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq; kraken2-mpa-report --db $KRAKEN2DB  $KRAKEN2/${bn}.kraken > $KRAKEN2/${bn}.kraken.mpa.report)"
  	done | xsbatch -c 1 --mem-per-cpu=30G -J kraken -R --max-array-jobs=10 --
  touch .kraken
fi



##### KRAKEN #####
KRAKEN2V2=$PROJECT/kraken2_v2
echo "Running kraken"
# -p: no error if existing 
mkdir -p $KRAKEN2V2 && cd $KRAKEN2V2


if [ ! -e .kraken ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $UNPREY/*_primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq.gz
  		do
    		bn=$(basename $f _primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq.gz)
    		# Run kraken, report, mpa report
    		echo "(gunzip $f; kraken2 --db $KRAKEN2DB --threads 8 --preload --unclassified-out $KRAKEN2V2/${bn}.unclassified.kraken --output $KRAKEN2V2/${bn}.kraken $UNPREY/${bn}_primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq; kraken2-report --db $KRAKEN2DB  $KRAKEN2V2/${bn}.kraken  >  $KRAKEN2V2/${bn}.kraken.report; gzip  $UNPREY/${bn}_primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq; kraken2-mpa-report --db $KRAKEN2DB  $KRAKEN2V2/${bn}.kraken > $KRAKEN2V2/${bn}.kraken.mpa.report)"
  	done | xsbatch -c 1 --mem-per-cpu=40G -J kraken -R --max-array-jobs=10 --
  touch .kraken
fi


##### KRAKEN #####
# -minimum-base-quality 35

KRAKEN2V3=$PROJECT/kraken2_k35
echo "Running kraken"
# -p: no error if existing 
mkdir -p $KRAKEN2V3 && cd $KRAKEN2V3


if [ ! -e .kraken ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $UNPREY/*_primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq.gz
  		do
    		bn=$(basename $f _primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq.gz)
    		# Run kraken, report, mpa report
    		echo "(gunzip $f; kraken2 --db $KRAKEN2DB --threads 8 --preload --minimum-base-quality 35 --unclassified-out $KRAKEN2V3/${bn}.unclassified.kraken.fastq --output $KRAKEN2V3/${bn}.kraken $UNPREY/${bn}_primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq; kraken2-report --db $KRAKEN2DB  $KRAKEN2V3/${bn}.kraken  >  $KRAKEN2V3/${bn}.kraken.report; gzip  $UNPREY/${bn}_primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq; kraken2-mpa-report --db $KRAKEN2DB  $KRAKEN2V3/${bn}.kraken > $KRAKEN2V3/${bn}.kraken.mpa.report)"
  	done | xsbatch -c 1 --mem-per-cpu=40G -J kraken -R --max-array-jobs=10 --
  touch .kraken
fi


##### KRAKEN #####
# -minimum-base-quality 35
# generate report 


UNPREY2=/groups/hologenomics/sarai/data/batProject/unpreys_2

KRAKEN2V3R=$PROJECT/kraken2_k35Report
echo "Running kraken"
# -p: no error if existing 
mkdir -p $KRAKEN2V3R && cd $KRAKEN2V3R


if [ ! -e .kraken ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $UNPREY2/*_primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq.gz
  		do
    		bn=$(basename $f _primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq.gz)
    		# Run kraken, report, mpa report
    		echo "(gunzip $f; kraken2 --db $KRAKEN2DB --threads 8 --preload --minimum-base-quality 35 --unclassified-out $KRAKEN2V3R/${bn}.unclassified.kraken.fastq --report $KRAKEN2V3R/${bn}.report.kraken --output $KRAKEN2V3R/${bn}.kraken $UNPREY2/${bn}_primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq; gzip  $UNPREY2/${bn}_primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq)"
  	done | xsbatch -c 1 --mem-per-cpu=40G -J kraken -R --max-array-jobs=10 --
  touch .kraken
fi

#Check to have report kraken 
kraken2 --report test.report --db $KRAKEN2DB --threads 8 --preload --minimum-base-quality 35 --output test.kraken /groups/hologenomics/sarai/data/batProject/unpreys/TOG_KLEY_121_GACGAC_primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq

#################FASTQC#######################3
module load java/v11.0.1-jdk
module load fastqc/v0.11.8

FASTQC=$PROJECT/fastqc
mkdir -p $FASTQC && cd $FASTQC


if [ ! -e .fastqc ]; then
	## Run BWA mem to map
	# Discard any reads that are orphans of a pair, so keep only valid pairs
  	for f in $UNPREY/*_primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq.gz
  		do
    		bn=$(basename $f _primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq.gz)
    		# Run kraken, report, mpa report
    		echo "(gunzip $f; kraken2 --db $KRAKEN2DB --threads 8 --preload  --output $KRAKEN2/${bn}.kraken $UNPREY/${bn}_primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq; kraken2-report --db $KRAKEN2DB  $KRAKEN2/${bn}.kraken  >  $KRAKEN2/${bn}.kraken.report; gzip  $UNPREY/${bn}_primer_UNMap_BatME.markdup.fastq.gz_Map_UNMapBatCowPigSheepDonkeyHorseChickenRhino.markdup.fastq; kraken2-mpa-report --db $KRAKEN2DB  $KRAKEN2/${bn}.kraken > $KRAKEN2/${bn}.kraken.mpa.report)"
  	done | xsbatch -c 1 --mem-per-cpu=5G -J kraken -R --max-array-jobs=10 --
  touch .fastqc
fi

##### KRONA #####
#Visualization of kraken results
# Run in personal computer

# to run the krona plots
for file in *.kraken
do
base=`basename $file .kraken`
echo "ktImportTaxonomy -q 2 -t 3 ${file} -o ${base}_kraken_krona.html"
#echo ${file}
done



