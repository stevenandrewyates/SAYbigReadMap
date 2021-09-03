##################################################
######            Details                   ######
##################################################

# author:   Steven Yates
# contact:  steven.yates@usys.ethz.ch
# Year:     2021
# Citation: TBA

##################################################
######            Description               ######
##################################################

# This script will align many fastq files against 
# a genome.

# It assumes all of the fastq data split already
# using "01_SAYsplit.sh" in this repository. 

# the data should be stored in a folder, e.g "SPLIT".
# Which contains sub folders for each sample (original
# fastq file). The script willl then make an output 
# directory called SAMand then populate it folders 
# for each sample, that conta sorted BAM files 
# aligned using bowtie2.

##################################################
######              Usage                   ######
##################################################

# This requires the following inputs
 
# -FG is the directory with the GENOME bowtie library
# -S is the directory to put the split data

# bash 02_SAYmanyReadMap.sh -S split2 -G GENOME

##################################################
######              Script                  ######
##################################################

while getopts S:G: flag
do
    case "${flag}" in
        G) GENOME=${OPTARG};;
        S) SPLIT=${OPTARG};;
    esac
done

# load bowtie2
module load gcc/4.8.2 bowtie2/2.2.3
module load gdc
module load samtools/1.7

mkdir SAM

GEN=$(ls $GENOME | grep '1.bt2' | grep -v rev | sed 's/.1.bt2//g')

for x in $(ls $SPLIT | sed 's/\..*//g' | uniq);
	do echo mkdir SAM/$x;done | sh

for x in $(ls $SPLIT | sed 's/\..*//g' | uniq )
	do for y in $(ls $SPLIT/${x}.1);
		do echo "bsub \"bowtie2 --no-unal -x $GENOME/$GEN -1 $SPLIT/${x}.1/${y} -2 $SPLIT/${x}.2/${y} | samtools view -bSu - | samtools sort -o SAM/${x}/$y.bam - \"";
	done;
done | sh


