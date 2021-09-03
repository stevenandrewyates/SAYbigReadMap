##################################################

# author:   Steven Yates
# contact:  steven.yates@usys.ethz.ch
# Year:     2021
# Citation: TBA

##################################################
######            Description               ######
##################################################

# get summary data after read mapping
 
##################################################
######              Usage                   ######
##################################################

# This requires the following inputs
 

# -S is the directory with the split fastq file data

# bash 02.1_SAYhowManyMapped.sh -S split2/

# it is essential that you add the "/" to directory!

# This script just outputs the data to the terminal.
# If you want to save it redirect to an output file

# bash 02.1_SAYhowManyMapped.sh -S split2/ > Summary.txt

##################################################
######              Script                  ######
##################################################

while getopts S: flag
do
    case "${flag}" in

        S) SPLIT=${OPTARG};;
    esac
done

# Check the output

# fnd how many fatsq files you have
FOUNDFILES=$(find $SPLIT -type f | wc -l)
# make sure it's an integer
declare -i FOUNDFILES

# setup a variable to divide by two
two=2

# now divide the number of fastq files by two
let h=$FOUNDFILES/$two
declare -i h

# find how many files were successfully completed
SUCCESS=$(grep "Successfully completed" lsf.* | wc -l)
declare -i SUCCESS

echo "I found $FOUNDFILES fastq files in the directory $SPLIT"
echo "and I found $SUCCESS Successfully completed lsf files."

if [ $h == $SUCCESS ]
	then
	echo "The number of paired fastq file($h) matches the number of Successfully completed ($SUCCESS) lsf files"
	echo "Happy days."
	else
	echo "The number of paired fastq file($h) does not matches the number of Successfully completed ($SUCCESS) lsf files"
	echo "Bad times. Maybe check the that all the lsf files are from read mapping?"
	fi
	

grep "overall" lsf.* | sed 's/%.*//g' | sed 's/.*://g' > OVERALL.txt	
OVER=$(awk '{ total += $1; count++ } END { print total/count }' OVERALL.txt)
	
echo "The mean mapping rate over all samples was $OVER."	
echo "You can find the overall mapping rates in the file \"OVERALL.txt\"" 
echo "Don't forget to write this information to a file!"
