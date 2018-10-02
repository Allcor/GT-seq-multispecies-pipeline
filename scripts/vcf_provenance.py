#!/usr/local/bin/python3.6
"""
When using tools like samtools the vcf file will receive header lines containing
for instance version and parameters used in producing the file.
For vcf_filter.py this is not the case.
"""

from subprocess import Popen,PIPE
from argparse import ArgumentParser
import sys
import gzip
import datetime


def argument_parser():
    parser = ArgumentParser(description='Options for fastq RG tagger')
    parser.add_argument(
        '-v', '--vcf',
        required=True,
        help="the filtered vcf file to add the vcf_filter.py provenance to"
    )
    parser.add_argument(
        '-g', '--gff3',
        help="the gff3 file used to add filters with the vcf_filter.py"
    )
    parser.add_argument(
        '-o', '--output',
        help="the location to write the changed vcf file to. (default:stdout)",
    )
    parser.add_argument(
        '-l', '--local',
        help="the local script containing the custum filters"
    )
    parser.add_argument(
        '-c', '--command',
        help="the command used to filter the vcf"
    )
    return parser.parse_args()


def main(args):
    #create the meta lines
    pipe = Popen("pip3 freeze | grep PyVCF", shell=True, stdout=PIPE).stdout
    date = datetime.datetime.now().strftime("%a %b %d %X %Y")
    version = pipe.read().decode("utf-8")
    meta = "##filter_vcfVersion=" + version.strip() + "\n"
    ###bcftools_callCommand=call -m --threads 3 /nfs/data/Big_Data/L15-217_framboos/results/arlo/gt_seq_22-02-18.framboos.vcf.gz; Date=Wed Apr 25 09:24:29 2018
    if args.local:
        meta += "##exampleCommand=vcf_filter.py --no-short-circuit --local-script {} [filters below]; Date={}\n".format(args.local,date)
    if args.gff3:
        meta += "##gff3_file:" + args.gff3 + "\n"
    #set the outfile
    if args.output:
        outfile = open(args.output,'w')
    else:
        outfile = sys.stdout
    #write the file with vcf_filter meta info to outfile
    meta_added = False
    with gzip.open(args.vcf, 'rb') as infile:
        for line in infile:
            if line.startswith(b"##INFO") and not meta_added:
                outfile.write(meta)
                meta_added = True
            outfile.write(line.decode("utf-8"))
    outfile.close()


if __name__ == '__main__':
    args = argument_parser()
    main(args)