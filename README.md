# SAYbigReadMap
A pipeline for mapping millions of DNA reads to a genome using bsub on Euler. 

# Description

Imagine the sequencing gods have been kind and furnished you with a whole NOVAseq of data, great. But how best to map all of that data efficiently? For one you can use a local machine and just let it run, or you could distribute it across a bsub system. This walkthough guides the process of distributing jobs using bsub.

The idea is simple split the data into small chunks and submit many small jobs. The advantage is the cluster will distribute the work easily and quickly. Also you don't run the risk of asking for not enough time and losing all your efforts.

# Data preparation

To begin we will download some data as an example. In this case we will work with an old favourite of mine red clover (*Trifolium pratense*). To make life easy we will take adavntage of previous [GitHub repository](https://github.com/stevenandrewyates/SAYReadMappingDNA) to download and install the reference genome. Which will download and install the red clover genome from [Ensembl plants](https://plants.ensembl.org/Trifolium_pratense/Info/Index).
```
cd $SCRATCH
mkdir RED
cd RED
module load git
git clone https://github.com/stevenandrewyates/SAYReadMappingDNA
git clone https://github.com/stevenandrewyates/SAYbigReadMap
sh SAYReadMappingDNA/01_DownloadGenome.sh -f ftp://ftp.ensemblgenomes.org/pub/plants/release-51/fasta/trifolium_pratense/dna/Trifolium_pratense.Trpr.dna_rm.toplevel.fa.gz
```
Next we need some data. For this we will another use a previous [repository](https://github.com/stevenandrewyates/SAYEulerDataManagement) and download data from two whole genome sequencing runs. This will be done using the SRA toolkit and downloads around 100 million reads per sample.
```
git clone https://github.com/stevenandrewyates/SAYEulerDataManagement
sh SAYEulerDataManagement/01_DownloadSRAtoolkit.sh
mkdir FASTQ
$HOME/sratoolkit.2.11.0-centos_linux64/bin/fastq-dump --split-files -X 91000000 -Z SRR1532811 | paste - - - - - - - - | tee >(cut -f 1-4 | tr \"\\t\" \"\\n\" > FASTQ/SRR1532811.1.fastq) | cut -f 5-8 | tr \"\\t\" \"\\n\" > FASTQ/SRR1532811.2.fastq
$HOME/sratoolkit.2.11.0-centos_linux64/bin/fastq-dump --split-files -X 91000000 -Z SRR847533 | paste - - - - - - - - | tee >(cut -f 1-4 | tr \"\\t\" \"\\n\" > FASTQ/SRR847533.1.fastq) | cut -f 5-8 | tr \"\\t\" \"\\n\" > FASTQ/SRR847533.2.fastq
```


Before starting with the data analysis we will first chop the data up into chunks of one million reads. This simply uses the (`split`) command which splits the input fastq files every four million lines (`split -l 4000000`). Being as each *read* has four lines per fastq file, this will split into chunks of one million reads. Besides this we use a (`for`) loop to iterate over all files and make sub-folders ((`mkdir`)). In case you're wondering I use my prefered method of having the loops (`echo`) the commands, so I can check them before passing ((`|`)) them to shell ((`sh`))

```
mkdir SPLIT
for x in $(ls $FASTQ | sed 's/.fastq.*//g');
	do echo mkdir $SPLIT/$x; 
	echo "cat $FASTQ/${x}.fastq | split -l 4000000 - $SPLIT/$x/";
done | sh
```

Splitting the data can be done with one script:
```
bash SAYbigReadMap/01_SAYsplit.sh -F FASTQ -S SPLIT -Z FALSE
```
The only thing added is if the data is *gzip* or not. This is specified in the (`-Z`) option either TRUE or FALSE.

# Read mapping

Time for the fun stuff read mapping. Like in other repositories I'll be using my prefered aligner [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). First we need to load bowtie2 and SAMtools. Then make a directory (`mkdir`) to output the data. In this example we'll call it **SAM**. FYI, you need SAMtools v.1.7 from the *gdc* module because the default (SAMtools v.1.2) doesn't handle pipes (|), which are used later.

```
module load gcc/4.8.2 bowtie2/2.2.3
module load gdc
module load samtools/1.7
mkdir SAM
```

Before starting we will make sub-folders for each fastq file. This uses a (`for`) loop like before to find the files. It then uses (`sed`) to remove any post-fixes (anything after a full stop) and get a unique (`uniq`) list of them. The aim being to find the sample names, and make only one SAM folder. In reality you have two files per sample so removing the ".1.fastq" and ".2.fastq" for forwards and reverse reads then removing duplicates does this.

```
for x in $(ls SPLIT | sed 's/\..*//g' | uniq);
	do echo mkdir SAM/$x;done | sh
```

Finally read mapping! For this a double (`for`) loop will find all the fastq files. These are then passed to (`echo`) which writes the bowtie2 commands. In addition two pipes ((`|`)) are used to convert the *SAM* files into  *BAM* file then a sorted *BAM* file: using SAMtools.

```
for x in $(ls SPLIT | sed 's/\..*//g' | uniq )
	do for y in $(ls SPLIT/${x}.1);
		do echo "bsub \"bowtie2 --no-unal -x GENOME/$GEN -1 $SPLIT/${x}.1/${y} -2 $SPLIT/${x}.2/${y} | samtools view -bSu - | samtools sort -o SAM/${x}/$y.bam - \"";
	done;
done | sh
```

For ease of use the whole thing can be done with one script:
```
bash SAYbigReadMap/02_SAYmanyReadMap.sh -S SPLIT -G GENOME
```
You can probably guess what the variables **SPLIT** and **GENOME** are?

Maybe you're suspicious if things worked? So you can get information about whether everything worked. Forgive me for not going into detail but I prepared a script which handles this. By default it writes the output to the terminal. If you want to save the information (besides copy and paste) you can redirect it to a file

```
bash SAYbigReadMap/02.1_SAYhowManyMapped.sh -S SPLIT/
bash SAYbigReadMap/02.1_SAYhowManyMapped.sh -S SPLIT/ > Summary_data.txt
```

# Merging BAM files

So far so good, you have mapped all the data. Not so good it's all split up. The split files will be named aa, ab, ac ... cz etc. To optimise merging, not all files will be merged at once. Because having hundreds of files open is tricky, even for a computer. So it will combine all files that start with "a" together and all files that start with "b" together, etc.,. The new merged *BAM*s will be post-fixed with 1M.sort.bam. This way they can be merged after. However do not run the next part, until this one has finished. Sorry you have to wait. The double (`for`) loop should look familiar to before. Also we are using (`SAMtools merge`) for this. The trick is using a wildcard (*) to match all files starting with "a", "b", "c" or "d". 

```
for f in $(ls SAM); 
	do for x in $(ls SAM/$f | cut -c 1 | uniq);
		do echo bsub samtools merge SAM/$f/${x}1M.sort.bam SAM/$f/${x}\*;
		done;
	done | sh
```

Like before you can do this easily with:

```
bash  SAYbigReadMap/03_SAYmergeStepOne.sh -S SAM
```

The very last step is combining all of those merged bams from before. Easy. First we make a list of files, per sample. Then pass this to (`SAMtools merge`). Again (`for`) loops are used to iterate over the files. To put all the finished files in one place we make the directory *BAMsorted* and write the data to it.

```
for f in $(ls $SAM );do echo "ls -1 -d $SAM/${f}/* | grep 1M > $SAM/$f/BAM.list"; done | sh
mkdir BAMsorted
for f in $(ls $SAM );do echo "bsub samtools merge BAMsorted/$f.sorted.bam -b $SAM/$f/BAM.list"; done | sh
```

Again this can all be done with one script

```
bash  SAYbigReadMap/04_SAYmergeStepTwo.sh -S SAM
```

Voila.

# All together

For simplicity I put everything together at the end. 
```
bash SAYbigReadMap/01_SAYsplit.sh -F FASTQ -S SPLIT -Z FALSE
bash SAYbigReadMap/02_SAYmanyReadMap.sh -S SPLIT -G GENOME
bash SAYbigReadMap/02.1_SAYhowManyMapped.sh -S SPLIT/
bash SAYbigReadMap/02.1_SAYhowManyMapped.sh -S SPLIT/ > Summary_data.txt
bash  SAYbigReadMap/03_SAYmergeStepOne.sh -S SAM
bash  SAYbigReadMap/04_SAYmergeStepTwo.sh -S SAM
```
You probably want to call some polymorphisms? So why not use the [following](https://github.com/stevenandrewyates/SAYVariantDetection]):
```
git clone https://github.com/stevenandrewyates/SAYVariantDetection
sh SAYVariantDetection/01_simpleVCF.sh -b BAMsorted -g GENOME/genome.fasta -o output.vcf -d 20
```
