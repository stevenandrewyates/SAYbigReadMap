##################################################

# author:   Steven Yates
# contact:  steven.yates@usys.ethz.ch
# Year:     2021
# Citation: TBA

##################################################
######            Description               ######
##################################################

# This script combines bam files using samtools merge

# It assumes all of the fastq data has already been 
# split and mapped using :

# "01_SAYsplit.sh" &
# "02_SAYmanyReadMap.sh" &
# "03_SAYmergeStepOne.sh

# in this repository. 

# using the aforementioned scripts an original fastq 
# file will be divided into chunks of one million reads,
# then aligned to the genome and converted into bam
# files. This script then combines the first step of 
# merged bams. 


##################################################
######              Usage                   ######
##################################################

# This requires the following inputs
 

# -S is the directory with the bam file data

# bash 04_SAYmergeStepTwo.sh -S SAM

##################################################
######              Script                  ######
##################################################

module load gdc
module load samtools/1.7


while getopts S: flag
do
    case "${flag}" in

        S) SAM=${OPTARG};;
    esac
done

for f in $(ls $SAM );do echo "ls -1 -d $SAM/${f}/* | grep 1M > $SAM/$f/BAM.list"; done | sh
mkdir BAMsorted
for f in $(ls $SAM );do echo "bsub samtools merge BAMsorted/$f.sorted.bam -b $SAM/$f/BAM.list"; done | sh


