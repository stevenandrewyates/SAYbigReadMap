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
# " 02_SAYmanyReadMap.sh"

# in this repository. 

# using the aforementioned scripts an original fastq 
# file will be divided into chunks of one million reads,
# then aligned to the genome and converted into bam
# files. In this step it will start making sorted- 
# merged bam files. 

# the split files will be named aa, ab, ac ... cz etc.
# To optimise merging not all files will be merged at once.
# Because having hundreds of files open is tricky, even
# for a computer. So it will combine all files that 
# start with "a" together and all files that start 
# with "b" together. The new merged bams will be 
# post-fixed with 1M.sort.bam. This way they can be 
# merged later (in the next script). However do not 
# run the next part, until this one has finished. 
# Sorry you have to wait 

 
##################################################
######              Usage                   ######
##################################################

# This requires the following inputs
 

# -S is the directory with the bam file data

# bash 03_SAYmergeStepOne.sh -S SAM

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


for f in $(ls $SAM); 
	do for x in $(ls $SAM/$f | cut -c 1 | uniq);
		do echo bsub samtools merge $SAM/$f/${x}1M.sort.bam $SAM/$f/${x}\*;
		done;
	done | sh






