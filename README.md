see also [Snakefile](http://gitlab.naktuinbouw.net/bioinformatics/gt-seq/blob/master/new_pipeline/Snakefile) and [config.yaml](http://gitlab.naktuinbouw.net/bioinformatics/gt-seq/blob/master/new_pipeline/config.yaml)


#### merging reads with Pear 


#### Demultiplexing

 - **script:**
   - [RG_tagger.py](http://gitlab.naktuinbouw.net/bioinformatics/gt-seq/blob/master/new_pipeline/RG_tagger.py)
 - **purpose:**
   - Splitting up the reads belonging to different species
   - tagging reads with the sample they belong to (read group). 
 - **input:**
   - fastq file coming from a GT_seq run
   - xlsx/xlsm/csv file with barcodes and corresponding sample ID's
 - **output:**
   - fastq file with the readgroup tags now attached.


#### Mapping the reads with bwa mem

#### Adding readgroups

 - **script:**
   - [make_header.py](http://gitlab.naktuinbouw.net/bioinformatics/gt-seq/blob/master/new_pipeline/make_header.py)
 - **purpose:**
   - including the readgroups to the header of a SAM/BAM file
 - **input:**
   - xlsx/xlsm/csv file with barcodes and corresponding sample ID's
   - sam/bam withought readgroups in header
 - **output:**
   - file with just the sam/bam header now including the readgroups. 


#### Sorting and reheadering the sam/bam
 
#### Combining the reads for each sample with mpileup

#### Calculating the likelyhood of a informative SNP with bcftools call

#### Filtering vcf

 - **script:**
    - vcf_filter.py
    - [vcf_locationselect.py](http://gitlab.naktuinbouw.net/bioinformatics/gt-seq/blob/master/new_pipeline/vcf_locationselect.py)
 - **purpose:**
    - indicate locations previously selected
    - indicate other locations with informative SNP
 - **input:**
    - gff3 indicating the SNP positions
    - called vcf file 
 - **output:**
    - vcf file with filters


#### Removing the uninformative locations from the vcf file

 - **script:**
    - vcf_filter.py
    - [vcf_locationselect.py](http://gitlab.naktuinbouw.net/bioinformatics/gt-seq/blob/master/new_pipeline/vcf_locationselect.py)
 - **purpose:**
    - removing all the positions we are not interested in.
 - **input:**
    - vcf file with filters
 - **output:**
    - vcf file only containing interesting calls

#### Creating pics for each SNP

 - **script:**
   - [GTseq_SummaryFigures_vcf.py](http://gitlab.naktuinbouw.net/bioinformatics/gt-seq/blob/master/new_pipeline/GTseq_SummaryFigures_vcf.py)
 - **purpose:**
   - create images showing quality of the SNP's among all samples in the gt-seq library
 - **input:**
   - called vcf file 
   - TODO nate's genos
 - **output:**
   - for each SNP a plot
   - text file with statistics for each SNP


![dag.svg](new_pipeline/dag.png)   