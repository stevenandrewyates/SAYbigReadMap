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

# A script to split fastq data files into 1 million 
# read chunks

# This script will make directory (for example split)
# then it will make a new folder for each fastq file 
# in the split directory.

# it will then populate the split folder with chunks 
# of 1 million reads. 

# These will be called aa, ab, ac etc... 
##################################################
######              Usage                   ######
##################################################

# This requires the following inputs
 
# -F is the directory with the FASTQ data
# -S is the directory to put the split data
# -Z zipped or not, use "TRUE" or "FALSE"

# bash 01_SAYsplit.sh -F FASTQ -S SPLIT -Z FALSE


##################################################
######              Script                  ######
##################################################

while getopts F:S:Z: flag
do
    case "${flag}" in
        F) FASTQ=${OPTARG};;
        S) SPLIT=${OPTARG};;
        Z) ZIP=${OPTARG};;
    esac
done

echo "$ZIP"

echo "making directory $SPLIT"
mkdir $SPLIT
# you need to make a folder for the files which will be split
# then you need to split the files

echo "populating it with:" 
for x in $(ls $FASTQ | sed 's/.fastq.*//g');
	do echo mkdir $SPLIT/$x; 
	echo "made directory $SPLIT/$x"
	if [ $ZIP == "FALSE" ]
		then 
		echo "cat $FASTQ/${x}.fastq | split -l 4000000 - $SPLIT/$x/";
		fi
	if [ $ZIP == "TRUE" ]
		then 
		echo "zcat $FASTQ/${x}.fastq.gz | split -l 4000000 - $SPLIT/$x/";
		fi
done | sh




