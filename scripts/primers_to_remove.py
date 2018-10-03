#!/usr/bin/python3
"""
Script reads gff3 and returns a list with all the primers that can be removed and the reason why.

input:
  gff3

output:
  csv with the primers that are disabled in the gff3
"""

from gff3_manager import gff3_file
from argparse import ArgumentParser
import sys

def argument_parser():
    parser = ArgumentParser(description='Options for fastq RG tagger')
    parser.add_argument(
        '-g', '--gff3',
        required=True,
        help="the gff3 file writen with the primers as gff3 features"
    )
    parser.add_argument(
        '-o', '--output',
        help="the location to write the disabled primers to. (default:stdout)"
    )
    return parser.parse_args()

def main(args):
    # read the gff3 and put it in memory
    gff3 = gff3_file(args.gff3)
    # get all disabled primers
    disabled_snp_list = []
    for group in [x for x in gff3.sets if not isinstance(x, str)]:
        if group.pcr_product and group.are_all_snp_disabled():
            disabled_snp_list.append(group)
    # write the outfile
    if args.output:
        outfile = open(args.output, 'w')
    else:
        outfile = sys.stdout
    for item in disabled_snp_list:
        left_primer_id = ','.join([x.attributes.id for x in item.forward_primer])
        left_primer_name = ','.join([x.attributes.name for x in item.forward_primer])
        right_primer_id = ','.join([x.attributes.id for x in item.reverse_primer])
        right_primer_name = ','.join([x.attributes.name for x in item.reverse_primer])
        for snp in item.snp:
            snp_id = snp.attributes.id
            snp_reasons = ','.join([' '.join(x) for x in snp.attributes.history])
        outfile.write('{}\n'.format(';'.join([
            left_primer_id,
            left_primer_name,
            right_primer_id,
            right_primer_name,
            snp_id,
            snp_reasons
        ])))
    outfile.close()

if __name__ == '__main__':
    args = argument_parser()
    main(args)