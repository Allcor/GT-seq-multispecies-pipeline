#!/usr/bin/python3

"""
adds the regroups to the header of a sam/bam

Input is
 1: sam/bam withought readgroups in header
 2: xlsx file with barcodes and corresponding sample ID's
Output is
 1: header.sam, all the Read groups added to the header of input sam/bam
"""

from argparse import ArgumentParser
from openpyxl import load_workbook
import pysam
import re
import RG_tagger

def argument_parser():
    parser = ArgumentParser(description='Options for fastq RG tagger')
    parser.add_argument(
        '-b', '--bam',
        required=True,
        help="the bam file that will be used to establish the read groups"
    )
    parser.add_argument(
        '-i', '--info',
        required=True,
        help="the .csv or .xlsx file with information on the read groups"
    )
    parser.add_argument(
        '-o', '--output',
        help="""
            the location where the new header.sam file will be created. 
            If not specified sys.stdout will be used
        """
    )
    parser.add_argument(
        '-s', '--species',
        help="the species of the readgroups to add, default is all",
        nargs='+',
        default=None
    )
    return parser.parse_args()


def make_readgroups_list(library_info_lines,species_list=None):
    readgroups_list = []
    for library_info in library_info_lines:
        readgroup_dict = RG_tagger.make_readgroup(library_info)
        if not species_list or library_info['Species'].lower() in species_list:
            readgroups_list.append(readgroup_dict)
    return readgroups_list


def main(args):
    # establish the readgroups
    info_filetype = args.info.split('.')[-1]
    if info_filetype == "csv":
        library_info_list = RG_tagger.obtain_barcode_info_csv(args.info)
    elif info_filetype in ["xlsx","xlsm"]:
        library_info_list = RG_tagger.obtain_barcode_info_xlsx(args.info)
    else:
        raise Exception("readgroup file extention .{} is not supported, needs to be .csv or .xlsx".format(
            info_filetype
        ))

    # open the sam/bam file
    alignment_filetype = args.bam.split('.')[-1]
    if alignment_filetype == "bam":
        infile = pysam.AlignmentFile(args.bam, 'rb')
    elif alignment_filetype == "sam":
        infile = pysam.AlignmentFile(args.bam, 'r')
    else:
        raise Exception("alignment file extension is not supported, needs to be .sam or .bam")

    header = infile.header.to_dict()
    header['RG'] = make_readgroups_list(library_info_list, args.species)

    # open the outfile
    if args.output:
        outfile = pysam.AlignmentFile(args.output, "w", header=header)
    else:
        outfile = pysam.AlignmentFile("-", "w", header=header)


if __name__ == '__main__':
    args = argument_parser()
    main(args)