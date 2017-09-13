# TB-ML-RD
Machine learning approach for *in silico* RD analysis of WGS samples of *M.tuberculosis*
## Files 
**ML_RD.py**: the standart pipeline for prediction of RD deletions and lineage identification. <br>
**model_cat.pickle**: the pretrained model of CatBoost classifier (CatBoost technology developed by YANDEX LLC) used for RD classiication <br>
**RD_data.csv**: pre-defined list of RD's used for analysis and lineage prediction. <br>
**ref/**: folder with bowtie2 index of H37Rv reference genome <br>
## TB-ML-RD requires
Python 2.7 <br>
Pysam <br>
CatBoost <br>
SAMtools <br>
bowtie2 <br>
## Usage
```
python ML_RD.py -f filename -i input_folder -o output_folder -m mode
```
**Options**: <br>
```
-f FILENAME, --file=FILENAME  Specify filename
-i INPUT_FOLDER, --input=INPUT_FOLDER Specify input folder [Default: running directory]
-m MODE, --mode=MODE  interaction mode: if bam, process starts directly from
                        analyzing bam file, bam file should be sorted and
                        indexed; if fastq, files will be processed with
                        bowtie2 first, will look for files _1.fastq.gz and
                        _2.fastq.fz; for fastq --refi and --reff should be
                        specified
-o OUTPUT, --output=OUTPUT  Specify outuput folder [Default: running directory]
--refi=REFI Path to folder and basename  [Default: ref/NC_000962.3]
--reff=REFF Path to reference fasta [Default: ref/NC_000962.3.fasta]
```
**Examples**:<br>
Starting from .bam file:<br>
```
python ML_RD.py -f example -i example/ -m bam

```
Starting from .fastq.gz files:<br>
```
python ML_RD.py -f example -i example/ -m fastq
```

