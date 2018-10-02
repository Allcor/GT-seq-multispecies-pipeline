import os

'''
 ---  variables:  ---
'''

configfile: "config.yaml"

MAIN_DIR = os.getcwd()+'/'
SCRIPTS_DIR = "/bin/gt-seq/new_pipeline/"

BASE_NAME = MAIN_DIR+config['project_name']


rule all:
  input:
    #expand(BASE_NAME+'.{species}.merged.bam', species=config['species']),
    expand(BASE_NAME+'.{species}.filter.vcf.gz', species=config['species']),
    #expand(BASE_NAME+'.{species}.merged_species.vcf.gz', species=config['species']),
    expand(MAIN_DIR+"{species}_graphs/"+config['project_name']+'.{species}_SummaryData.txt', species=config['species']),
    expand(BASE_NAME+'.{species}.gff3', species=config['species']),
    expand(MAIN_DIR+"vcf_counts_{species}.csv", species=config['species']),
    expand(BASE_NAME+'.unassembled.{side}.fastq.gz', side=['forward','reverse']),
    BASE_NAME+".assembled.fastq.gz",


rule pear:
  input:
    config['fastq_forward'],
    config['fastq_reverse'],
  params:
    BASE_NAME
  threads:
    8
  output:
    temp(BASE_NAME+".assembled.fastq"),
    temp(BASE_NAME+".discarded.fastq"),
    temp(expand(BASE_NAME+'.unassembled.{side}.fastq', side=['forward','reverse']))
  shell:
    "pear -j {threads} -f {input[0]} -r {input[1]} -o {params}"


rule zip_reads_assembled:
  input:
    BASE_NAME+".assembled.fastq",
  output:
    BASE_NAME+".assembled.fastq.gz",
  shell:
    "pigz {input}"


rule zip_reads_unassembled:
  input:
    BASE_NAME+'.unassembled.{side}.fastq'
  output:
    BASE_NAME+'.unassembled.{side}.fastq.gz'
  shell:
    "pigz {input}"


rule read_split_assembled:
  input:
    BASE_NAME+".assembled.fastq",
    config['barcode_file']
  params:
    SCRIPTS_DIR+"RG_tagger.py",
    BASE_NAME+'.assembled.fastq'
  output:
    temp(expand(BASE_NAME+'.assembled.{species}.fastq', species=config['species'])),
    temp(BASE_NAME+'.assembled.nogroup.fastq')
  shell:
    "python3.6 {params[0]} -f {input[0]} -i {input[1]} -o {params[1]} -s -u"


rule read_split:
  input:
    BASE_NAME+".unassembled.{side}.fastq",
    config['barcode_file']
  params:
    SCRIPTS_DIR+"RG_tagger.py",
    BASE_NAME+'.{side}.fastq'
  output:
    temp(expand(BASE_NAME+'.{{side}}.{species}.fastq', species=config['species'])),
    temp(BASE_NAME+'.{side}.nogroup.fastq')
  shell:
    "python3.6 {params[0]} -f {input[0]} -i {input[1]} -o {params[1]} -s -u"


rule bwa_map_assembled:
  input:
    lambda wildcards: config["reference_genome"][wildcards.species],
    BASE_NAME+'.assembled.{species}.fastq',
  threads:
    10
  output:
    temp(BASE_NAME+'.assembled.{species}.bam')
  shell:
    "bwa mem -C -t {threads} {input} | samtools view -Shb - > {output}"


rule bwa_map:
  input:
    lambda wildcards: config["reference_genome"][wildcards.species],
    expand(BASE_NAME+'.{side}.{{species}}.fastq', side=['forward','reverse'])
  threads:
    10
  output:
    temp(BASE_NAME+'.unassembled.{species}.bam')
  shell:
    "bwa mem -C -t {threads} {input} | samtools view -Shb - > {output}"


rule reheader:
  input:
    config['barcode_file'],
    BASE_NAME+'.{assembled}.{species}.bam'
  params:
    SCRIPTS_DIR+"make_header.py"
    " -s {species}"
  output:
    temp(BASE_NAME+'.{assembled}.{species}.sorted.bam')
  shell:
    "python3.6 {params} -i {input[0]} -b {input[1]} | samtools reheader - {input[1]} | samtools sort - > {output}"


rule merge:
  input:
    expand(BASE_NAME+'.{assembled}.{{species}}.sorted.bam', assembled=['unassembled', 'assembled'])
  threads:
    10
  output:
    temp(BASE_NAME+'.{species}.merged.bam')
  shell:
    "samtools merge -c -@ {threads} {output} {input}"


rule check_bam:
  #TODO using the merged bam instead of only the assembled
  input:
    #BASE_NAME+'.assembled.{species}.sorted.bam',
    BASE_NAME+'.{species}.merged.bam',
    config['barcode_file']
  params:
    SCRIPTS_DIR+"mapping_stats.py"
  output:
    BASE_NAME+'.{species}.checked.bam',
    BASE_NAME+'.{species}.mapping_qual_log.xlsx'
  shell:
    "python3.6 {params} -b {input[0]} -l {output[1]} -o {output[0]}"


rule bedfile:
  input:
    lambda wildcards: config["gff3_file"][wildcards.species]
  params:
    '"PCR_product"',
    '-v "filtered=True"'
  output:
    temp(BASE_NAME+'.{species}.bed')
  shell:
    "grep {params[0]} {input} | grep {params[1]} | cut -f 1,4,5 > {output}"


rule mpileup:
  input:
    lambda wildcards: config["reference_genome"][wildcards.species],
    BASE_NAME+'.{species}.checked.bam',
    BASE_NAME+'.{species}.bed'
  output:
    temp(BASE_NAME+'.{species}.vcf.gz')
  shell:
    "samtools mpileup -vuf {input[0]} -l {input[2]} -a -d 1000000 -L 1000000 --output-tags DP,AD,ADF,ADR,SP {input[1]} | bgzip -c > {output}"


rule tabix:
  input:
    BASE_NAME+'.{species}.vcf.gz'
  output:
    temp(BASE_NAME+'.{species}.vcf.gz.tbi')
  shell:
    "tabix -p vcf {input}"


rule calling:
  input:
    BASE_NAME+'.{species}.vcf.gz'
  threads:
    3
  output:
    temp(BASE_NAME+'.{species}.called.vcf.gz')
  shell:
    "bcftools call -m --threads {threads} {input[0]} | bgzip -c > {output}"


rule tabix_calls:
  input:
    BASE_NAME+'.{species}.called.vcf.gz'
  output:
    temp(BASE_NAME+'.{species}.called.vcf.gz.tbi')
  shell:
    "tabix -p vcf {input}"


rule filter_vcf:
  input:
    #BASE_NAME+'.{species}.concat.vcf.gz'
    BASE_NAME+'.{species}.called.vcf.gz',
    lambda wildcards: config["gff3_file"][wildcards.species]
  params:
    "--local-script vcf_gtseq_filters.py",
    "sq", #site quality < 30
    "avg-dps --avg-depth-per-sample 5",
    "spread", # minimal 2 samples for each allele (ignored with less then 10 samples, DP must be more then 100)
    "gff3_notarget --gff3_target",
    "gff3_ignore --gff3_ignore"
  output:
    temp(BASE_NAME+'.{species}.filtered.vcf.gz')
  shell:
    "vcf_filter.py --no-short-circuit {params[0]} {input[0]} {params[1]} {params[2]} {params[3]} {params[4]} {input[1]} {params[5]} {input[1]} | bgzip -c > {output}"


rule add_provenance:
  input:
    BASE_NAME+'.{species}.filtered.vcf.gz',
    lambda wildcards: config["gff3_file"][wildcards.species]
  params:
    SCRIPTS_DIR+"vcf_provenance.py",
    #TODO: have the command of filter_vcf in the meta info as well.
    "vcf_gtseq_filters.py"
  output:
    BASE_NAME+'.{species}.filter.vcf.gz'
  shell:
    "python3.6 {params[0]} -v {input[0]} -g {input[1]} -l {params[1]} | bgzip -c > {output}"


rule cleanup_vcf:
  input:
    BASE_NAME+'.{species}.filter.vcf.gz'
  params:
    "failed_filtering"
  output:
    BASE_NAME+'.{species}.clean.vcf.gz'
  shell:
    "vcf_filter.py --no-filtered --local-script vcf_gtseq_filters.py {input} {params} | bgzip -c > {output}"


rule add_to_gff3:
  input:
    BASE_NAME+'.{species}.clean.vcf.gz',
    lambda wildcards: config["gff3_file"][wildcards.species]
  params:
    SCRIPTS_DIR+"add_discovered.py"
  output:
    BASE_NAME+'.{species}.gff3'
    #TODO: put the gff3 in a better place.
  shell:
    "python3.6 {params} -g {input[1]} -v {input[0]} -o {output}"


rule merge_samples:
  input:
    BASE_NAME+'.{species}.clean.vcf.gz',
    config['barcode_file']
  params:
    SCRIPTS_DIR+"merge_samples.py"
  output:
    BASE_NAME+'.{species}.merged_species.vcf.gz',
    BASE_NAME+'.merge_{species}_samples.log'
  shell:
    "python3.6 {params} -v {input[0]} -i {input[1]} -l {output[1]}| bgzip -c > {output[0]}"


'''
 ---  genos:  ---
 #TODO: should be removed if made obsolete
'''

if 'genos_folder' in config:
  GENOFOLDER = config['genos_folder']
else:
  GENOFOLDER = MAIN_DIR+"genos/"

filename_parts = ['i7_name','i5_name','Plate_name']

def fetch_barcode_info(infile, folder, filename_parts):
    #names as created in GTseq_BarcodeSplit
    from openpyxl import load_workbook
    filenames = []
    wb = load_workbook(infile, data_only=True)
    info_sheet = wb['Library']
    info_headers = []
    for column in info_sheet.iter_cols(max_row=1):
        info_headers.append(column[0].value)
    for row in info_sheet.iter_rows(min_row=2):
        info_dict = {info_headers[i]:item.value for i,item in enumerate(row)}
        # rows with empty cells are not allowed
        if all(x for x in info_dict.values()):
            filenames.append(folder+info_dict['Species'].lower()+'/'+'_'.join(
                ['_'.join(info_dict[x].split()) for x in filename_parts]
            ))
    return filenames

FILES = fetch_barcode_info(config['barcode_file'], GENOFOLDER, filename_parts)

rule barcode_split:
    input:
        meta = config['barcode_file'],
        #TODO will miss products that are to long to be merged
        fastq = BASE_NAME+".assembled.fastq"
    params:
        GENOFOLDER
    output:
        temp(expand('{sample}.fastq', sample=FILES))
    shell:
        "python3.6 /bin/gt-seq/new_pipeline/GTseq_BarcodeSplit_naktuinbouw.py -i {input.meta} -f {input.fastq} -d {params}"

rule genotyper:
    input:
        barcode = lambda wildcards: config["probeseq_file"][wildcards.sample.split('/')[-2]],
        fastq = "{sample}.fastq"
    output:
        "{sample}.genos"
    shell:
        '/bin/gt-seq/GTseq_Genotyper_v3.pl "{input.barcode}" "{input.fastq}" > {output}'


'''
 ---  figures:  ---
 #TODO: figure out how to make off-target plots without the geno's the geno's
'''


rule figures:
  input:
    lambda wildcards: config["gff3_file"][wildcards.species],
    BASE_NAME+'.{species}.clean.vcf.gz',
    expand('{sample}.genos', sample=FILES)
  params:
    SCRIPTS_DIR+"GTseq_SummaryFigures_vcf.py",
    MAIN_DIR+"{species}_graphs/"+config['project_name']+".{species}", #the start of the path used for the figures
    GENOFOLDER+'{species}/'
  output:
    MAIN_DIR+"{species}_graphs/"+config['project_name']+".{species}_SummaryData.txt"
  shell:
    "python3.6 {params[0]} -v {input[1]} -n {params[1]} -a -e {params[2]} -g {input[0]}"


rule sample_figures:
    input:
        BASE_NAME+'.{species}.clean.vcf.gz',
        config['barcode_file'],
        BASE_NAME+'.{species}.mapping_qual_log.xlsx'
    params:
        SCRIPTS_DIR+"vcf_stats.py",
        MAIN_DIR+"{species}_graphs/"
    output:
        MAIN_DIR+"vcf_counts_{species}.csv"
    shell:
        "python3.6 {params[0]} -v {input[0]} -i {input[1]} -o {output[0]} -c {input[2]} -p {params[1]} -g"