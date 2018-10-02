#!/usr/bin/python3
"""
Script to add undescovered SNP of a vcf file to the gff3 that have not previously been targeted.

input:
vcf file (called and filtered)
gff3 file

output:
new gff3 file with new snp added
"""

from gff3_manager import gff3_file
from argparse import ArgumentParser
import vcf

def argument_parser():
    parser = ArgumentParser(description='Options for fastq RG tagger')
    parser.add_argument(
        '-g', '--gff3',
        required=True,
        help="the gff3 file writen with the primers as gff3 features"
    )
    parser.add_argument(
        '-v', '--vcf',
        required=True,
        help="the .vcf file that has SNP not yet described in the gff3"
    )
    parser.add_argument(
        '-o', '--output',
        help="the location to write the changed gff3 file to. (default:stdout)"
    )
    return parser.parse_args()

def get_new_snp(vcf_file):
    """
    Gets the positions of the new snp in a vcf file
    :param vcf_file: py_vcf file
    :return: list of new snp
    """
    new_snp = []
    for loci in vcf_file:
        if "gff3_notarget" in loci.FILTER:
            new_snp.append(loci)
    return(new_snp)

def main(args):
    # read the gff3 and put it in memory
    gff3 = gff3_file(args.gff3)
    # discover the new SNP from the vcf file
    vcf_file = vcf.Reader(filename=args.vcf)
    snp_to_add = get_new_snp(vcf_file)
    # add the new SNP to the gff3
    for snp in snp_to_add:
        gff3.add_vcf_call(snp)
    # write the new gff3
    gff3.write_gff3(args.output)

if __name__ == '__main__':
    args = argument_parser()
    main(args)