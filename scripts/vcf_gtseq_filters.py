#!/usr/local/bin/python3.6

import vcf
import sys

# TODO: set pithonpath with install
sys.path.insert(0, '/bin/gt-seq/new_pipeline/')
import gff3_manager

class IsFeature(vcf.filters.Base):
    'Filter calls that are not SNP features in [--gff3_target](above gff3_file)'

    name = 'gff3_notarget'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument(
            '--gff3_target',
            help='vcf file containing the calls made for vcf being filtered'
        )

    def __init__(self, args):
        gff3 = gff3_manager.gff3_file(args.gff3_target)
        self.calls = gff3.give_active_snp()

    def __call__(self, record):
        if record.is_indel and not record.is_deletion:
            record_id = record.CHROM + '_' + str(record.POS) + '_' + record.REF + '_indel'
        else:
            record_id = record.CHROM+'_'+str(record.POS)+'_'+record.REF
        if record_id not in self.calls :
            return True

    def filter_name(self):
        return self.name


class InFile(vcf.filters.Base):
    'Filter calls that are noted in [--gff3_ignore](above gff3_file) as safe to ignore'

    name = 'gff3_ignore'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument(
            '--gff3_ignore',
            help='vcf file containing the calls made for vcf being filtered'
        )

    def __init__(self, args):
        gff3 = gff3_manager.gff3_file(args.gff3_ignore)
        self.calls = gff3.give_deactivated_snp()

    def __call__(self, record):
        if record.is_indel and not record.is_deletion:
            record_id = record.CHROM + '_' + str(record.POS) + '_' + record.REF + '_indel'
        else:
            record_id = record.CHROM+'_'+str(record.POS)+'_'+record.REF
        if record_id in self.calls :
            return True

    def filter_name(self):
        return self.name


class Spread(vcf.filters.Base):
    'Has a minimum of [--min_calls](see number) on each allele among samples'

    name = 'spread'

    @classmethod
    def customize_parser(self, parser):
        parser.add_argument(
            '--min_calls',
            type=int,
            default=2,
            help="the minimal amount of significant allele calls among all samples"
        )

    def __init__(self, args):
        self.threshold = args.min_calls

    def __call__(self, record):
        if len(record.samples) > (self.threshold*5):
            # if there are not enough classes, skip this step.
            nums = {'0/0':0,'0/1':0,'1/1':0}
            max_DP = max([i.data.DP for i in record.samples])
            for i in record.samples:
                if i.gt_nums in nums:
                    # if the depth is low it should not be taken into account
                    if i.data.DP >= max_DP*0.2:
                        nums[i.gt_nums] += 1
            if any(item < self.threshold for item in nums.values()):
                return True


class RemoveCalls(vcf.filters.Base):
    'Checks the existing filters, some filter combinations are allowed.'

    name = 'failed_filtering'

    def __call__(self, record):
        if 'gff3_notarget' in record.FILTER and len(record.FILTER) > 1:
            return True

    def filter_name(self):
        return self.name