#!/usr/bin/python3
"""
Script to make adjustments to the gff3 with SNP and Primers for GT-seq.

When a gff3 file is loaded it attempts to recognise custom values used for the gt_seq pipeline.
Loading should work for legacy gff3 files. If a gff3 is writen it will be writen in the newest layout.

While the main use is to load it as module, it can be used for making changes when run.

Input:
 - gff3 to adjust,
 - changes to be made (csv, tsv, json)

Output:
 - new gff3 with changes
"""

from argparse import ArgumentParser
#from openpyxl import load_workbook
#from pyfasta import Fasta

import sys
import datetime


def argument_parser():
    parser = ArgumentParser(description='Options for fastq RG tagger')
    parser.add_argument(
        '-g', '--gff3',
        required=True,
        help="the gff3 file writen with the primers as gff3 features"
    )
    parser.add_argument(
        '-i', '--info',
        help="the .tsv file with SNP to filter"
    )
    parser.add_argument(
        '-o', '--output',
        help="the location to write the changed gff3 file to. (default:stdout)"
    )
    return parser.parse_args()


class GFF3_attributes():
    """
    Attributes of a gff3_line (feature)

    Includes custom values: discovered, inspected, filter, filtered and allele
    """

    history_key = 'filter'
    active_key = 'filtered'
    allele_key = 'allele'
    discovered_key = 'discovered'
    validated_key = 'inspected'

    disabled = ['False','false']
    enabled = ['True','true','new_find']

    def __init__(self,owner):
        self._owner_line = owner
        self.id = ""
        self.name = ""
        self.alias = ""
        self.parent = ""
        self.target = ""
        self.gap = ""
        self.derives_from = ""
        self.note = ""
        self.dbxref = ""
        self.ontology_term = ""
        self.is_circular = ""
        self.discovered = None
        self.validated = None
        self.active = None
        self.history = []
        self.allele = []
        self.other = {}

    def parse(self,att_line):
        """
        function to read the attributes of a gff3 line.
        all the default entries are there,
        :param att_line: the attributes of a gff3 line
        :return: self
        """
        att_line = {item.split('=')[0]:item.split('=')[1] for item in att_line.split(';')}
        if 'ID' in att_line:
            self.id = att_line.pop('ID',None)
        if 'Name' in att_line:
            self.name = att_line.pop('Name',None)
        if 'Alias' in att_line:
            self.alias = att_line.pop('Alias',None)
        if 'Parent' in att_line:
            self.parent = att_line.pop('Parent',None)
        if 'Target' in att_line:
            self.target = att_line.pop('Target',None)
        if 'Gap' in att_line:
            self.gap = att_line.pop('Gap',None)
        if 'Derives_from' in att_line:
            self.derives_from = att_line.pop('Derives_from',None)
        if 'Note' in att_line:
            self.note = att_line.pop('Note',None)
        if 'Dbxref' in att_line:
            self.dbxref = att_line.pop('Dbxref',None)
        if 'Ontology_term' in att_line:
            self.ontology_term = att_line.pop('Ontology_term',None)
        if 'Is_circular' in att_line:
            self.is_circular = att_line.pop('Is_circular',None)
        if self.allele_key in att_line:
            self.allele = att_line.pop(self.allele_key,None).split(',')
        if self.discovered_key in att_line:
            discovered_val = att_line.pop(self.discovered_key,None)
            if discovered_val in self.enabled:
                self.discovered = True
            elif discovered_val in self.disabled:
                self.discovered = False
        if self.validated_key in att_line:
            validated_val = att_line.pop(self.validated_key,None)
            if validated_val in self.enabled:
                self.validated = True
            elif validated_val in self.disabled:
                self.validated = False
        if self.active_key in att_line:
            active_val = att_line.pop(self.active_key,None)
            if active_val in self.enabled:
                self.active = False
            elif active_val in self.disabled:
                self.active = True
            else:
                pass#print("{} is not a valid entry for {}".format(active_val,self.active_key))
        else:
            #TODO this will make every element withought existing value active.
            self.active = True
        rest_dict = {}
        for key,value in att_line.items():
            if self.history_key in key:
                parts = att_line[key].split(',')
                if len(parts) == 3:
                    self.add_history(*parts)
                else:
                    reason = '.'.join(parts[2:])
                    self.add_history(parts[0],parts[1],reason)
            else:
                rest_dict[key] = value
        if rest_dict != {}:
            #TODO: maybe add this to a logger instead
            #print("attributes put in other: {}".format(', '.join(att_line.keys())))
            self.other = att_line
        return self

    def compose(self):
        """
        function to constitute the attributes into a string to write in a gff3 file
        :return: str
        """
        return_lib = self.other
        return_lib.update({
            'ID' : self.id,
            'Name' : self.name,
            'Alias' : self.alias,
            'Parent' : self.parent,
            'Target' : self.target,
            'Gap' : self.gap,
            'Derives_from' : self.derives_from,
            'Note' : self.note,
            'Dbxref' : self.dbxref,
            'Ontology_term' : self.ontology_term,
            'Is_circular' : self.is_circular,
            'allele' : ','.join([a.strip() for a in self.allele])
        })
        if self.discovered != None:
            if self.discovered:
                return_lib[self.discovered_key] =  self.enabled[0]
            else:
                return_lib[self.discovered_key] = self.disabled[0]
        if self.validated != None:
            if self.validated:
                return_lib[self.validated_key] =  self.enabled[0]
            else:
                return_lib[self.validated_key] = self.disabled[0]
        if self.active != None:
            if self.active:
                return_lib[self.active_key] = self.disabled[0]
            else:
                return_lib[self.active_key] = self.enabled[0]
        for i,item in enumerate(self.history):
            return_lib["{}_{}".format(self.history_key,i)] = ','.join(item)
        return_str = ';'.join([key + '=' + value for key, value in return_lib.items() if value])
        return return_str

    def add_history(self,date,what,note):
        """
        add a change/filter to the feature
        :param date: the date the change was made
        :param what: who/what made the change
        :param note: reason change was made
        :return: changes self
        """
        note = '.'.join(note.split(','))
        self.history.append([date,what,note])

    def disable_feature(self,reason,source="gff3_maniger"):
        """
        disables this feature and tells the containing set
        :param source: who/what made the change
        :param reason: reason change was made
        :return: changes self
        """
        date = datetime.datetime.now().strftime("%Y-%m-%d")
        self.add_history(date,source,reason)
        self.active = False
        if self._owner_line.type == 'SNP':
            self._owner_line._owner_set.all_snp_disabled()

    def check_note_for_history(self):
        """
        because the SNP have been made inactive differently in the case of rasberry,
        this is here to add those to the history.
        :return:
        """
        testrun_notes = [
            "multiple loci suspected",
            "suspected multicopy, poor performance",
            "fixed allele 1",
            "very poor amplification",
            "very poor amplification, high off target percent",
            "poor amplification, maybe redesign",
            "mono-allele 1?",
            "redesign primer",
            "most of target",
            "poor performance",
            "poor performance, primers off target",
            "off target amp",
            "mono-allele 1",
            "mono-allele 2 and off target",
            "Nate said it is a mess",
            "off target amp",
            "mono-allele 1 and off target"
        ]
        if self.note == "No primers made by primer3":
            self.add_history("2018-2-12","Nate","primers were not made for this sequence variation")
            self.note = "sequence variant selected by GBS-SNP-selection"
        elif self.note == "Removed by nate, close to other SNP":
            self.add_history("2018-2-19","Nate","Primers designed for this SNP were taken out, were to close to other SNP")
            self.note = "sequence variant selected by GBS-SNP-selection"
        elif self.note == "Predicted to form hetrodymer":
            self.add_history("2018-2-19","Nate","Predicted to form hetrodymer")
            self.note = "sequence variant selected by GBS-SNP-selection"
        elif self.note == "no valid primer pair could be made for this position":
            self.note = "sequence variant selected by GBS-SNP-selection"
        elif self.note in testrun_notes:
            self.add_history("2018-2-23","Thomas",self.note)
            self.note = "sequence variant selected by GBS-SNP-selection"
        #check if any were missed.
        if self.active and self.note != "sequence variant selected by GBS-SNP-selection":
            pass #print(self.note)


class GFF3_line():
    """
    Each line of a gff3 file has 9 tab separated fields.
    the last one, attributes, can have special fields.
    In case of SNP, this contains the filter.
    """
    def __init__(self, set_owner, gff3_line=None):
        self._owner_set = set_owner
        self.seqid = ""
        self.source = ""
        self.type = ""
        self.start = 0
        self.end = 0
        self.score = ""
        self.strand = ""
        self.phase = ""
        self.attributes = GFF3_attributes(self)
        if gff3_line:
            self.parse(gff3_line)

    def parse(self,gff3_line):
        """
        parse a line of a gff3 file
        :param gff3_line: line of a gff3 file
        :return: self
        """
        split_line = gff3_line.strip().split('\t')
        self.seqid = split_line[0]
        self.source = split_line[1]
        self.type = split_line[2]
        self.start = int(split_line[3])
        self.end = int(split_line[4])
        self.score = split_line[5]
        self.strand = split_line[6]
        self.phase = split_line[7]
        self.attributes.parse(split_line[8])
        return self

    def compose(self):
        """
        constitute the string back together that can be writen to a gff3 file
        :return: string
        """
        return_str = "\t".join([
            self.seqid,
            self.source,
            self.type,
            str(self.start),
            str(self.end),
            self.score,
            self.strand,
            self.phase,
            self.attributes.compose()
        ])
        return return_str

    def is_new_snp(self,seqid,pos,allele):
        """
        fill in the fields required for a new SNP
        :param seqid: string, id of sequence containing SNP
        :param pos: intiger, position of SNP
        :param allele: list, alleles found on this position.
        :return:
        """
        self.seqid = seqid
        self.source = "gff3_manager"
        self.type = "SNP"
        self.start = pos
        self.end = pos
        self.score = "."
        self.strand = "+"
        self.phase = "."
        self.attributes.id = seqid+"_"+str(pos)
        self.attributes.note = "new sequence variant found after sequencing"
        self.attributes.allele = allele
        self.attributes.active = True
        self.attributes.discovered = True
        self.attributes.validated = False

    def primer_start_fix(self):
        #TODO this function will not be used anymore, remove?
        """
        By mistake the start of primers were not on the correct position.
        this function corrects this.
        :return: adjusting self.start
        """
        if self.type in ["forward_primer", "reverse_primer", "PCR_product"]:
            self.start += 1
        if self.type == "region" and self.source == "Primer3":
            # this is the region containing the primers
            self.start += 1


class GT_seq_location():
    """
    For each species we made a GT-seq primer set there should also be a gff3
    with detailed information on each primer set and the interesting SNP it contains.

    This class represents one of those primer and SNP combinations.
    """

    def __init__(self):
        self._flanking_region = None
        self.flanking_region = None
        self.gt_seq_region = []
        self.pcr_product = []
        self.forward_primer = []
        self.reverse_primer = []
        self.snp = []

    def __repr__(self):
        if self.flanking_region and self.flanking_region.attributes:
            site_id = self.flanking_region.attributes.id
        else:
            site_id = 'unnamed site'
        return "GT-seq {}".format(site_id)

    def compose(self):
        #TODO: don't have to check all are disabled all the time
        self.are_all_snp_disabled()
        #TODO: should have to remove redundant reagions only once
        if self._flanking_region:
            #redundant region object is present. The redundant regions need to be removed.
            self.remove_redundant_regions()
        return_list = []
        if self._flanking_region:
            return_list.append(self._flanking_region.compose())
        if self.flanking_region:
            return_list.append(self.flanking_region.compose())
        if self.gt_seq_region:
            for feature in self.gt_seq_region:
                return_list.append(feature.compose())
        if self.pcr_product:
            for feature in self.pcr_product:
                return_list.append(feature.compose())
        if self.forward_primer:
            for feature in self.forward_primer:
                return_list.append(feature.compose())
        if self.reverse_primer:
            for feature in self.reverse_primer:
                return_list.append(feature.compose())
        if self.snp:
            for feature in self.snp:
                return_list.append(feature.compose())
        return '\n'.join(return_list)+'\n'

    def remove_redundant_regions(self):
        """
        Se it applies to Sequence Ontology the flanking region and pcr product are double defined.
        This function removes the redundant regions.
        :return: self._flanking_region and self.gt_seq_region
        """
        self.flanking_region.attributes.id = self._flanking_region.attributes.id
        self.flanking_region.attributes.parent = ''
        for feature in self.pcr_product:
            feature.attributes.id = feature.attributes.parent
            feature.attributes.parent = ''
        self._flanking_region = None
        self.gt_seq_region = []
        if self.pcr_product:
            snp_parent = self.pcr_product[0].attributes.id
        else:
            snp_parent = self.flanking_region.attributes.id
        for snp in self.snp:
            snp.attributes.parent = snp_parent


    def remove_region_names(self):
        """
        regions had the same name as a SNP, made it hard to look up.
        :return: changes self
        """
        for region in self.gt_seq_region:
            region.attributes.name = ''

    def get_main_snp(self):
        """
        checks the SNP associated to this group and returns one that can be seen as most important.
        :return:
        """
        for snp in self.snp:
            if snp.attributes.active:
                if not snp.attributes.discovered:
                    return snp

    def are_all_snp_disabled(self):
        """
        checks if any interesting SNP are still in the amplified region.
        :return: true if primers don't amplify any interesting snp.
        """
        for snp in self.snp:
            if snp.attributes.active:
                return False
        self.primers_are_useless()
        return True

    def add_snp(self,snp_line):
        """
        add a SNP to the location
        :param: snp_line: GFF3_line object
        :return: changes self.snp
        """
        #TODO: make sure there are no cases anymore where there are multiple primer sets
        parent = self.pcr_product[0].attributes.id
        snp_line.attributes.parent = parent
        self.snp.append(snp_line)

    def primers_are_useless(self):
        """
        deal with the fact the primers are not amplifying anything interesting
        :return:
        """
        #TODO: send a message telling these primers can be taken out.
        for feature in self.gt_seq_region:
            if feature.attributes.active:
                feature.attributes.disable_feature("has no interesting sequence variation")
        for feature in self.pcr_product:
            if feature.attributes.active:
                feature.attributes.disable_feature("has no interesting sequence variation")
        for feature in self.forward_primer:
            if feature.attributes.active:
                feature.attributes.disable_feature("has no interesting sequence variation")
        for feature in self.reverse_primer:
            if feature.attributes.active:
                feature.attributes.disable_feature("has no interesting sequence variation")


class gff3_file():
    """
    class holding the gff3 document itself.

    should be able to loop on each gt-seq group
    but each feature should also be reachable with it's id directly as well
    """

    def __init__(self,gff3_file=None):
        self.sets = []
        self.features_id = {}
        self.features_name = {}
        self.sequence_regions = {}
        if gff3_file:
            self.read_gff3(gff3_file)

    def read_gff3(self,gff3_file):
        """
        read a gff3 file and read the GT-seq feature sets
        :param gff3_file: used for GT-seq SNP and Primer annotation.
        :return: populate self.sets with multiple GT_seq_location objects
        """
        with open(gff3_file) as infile:
            set = None
            for line in infile:
                if line[0] == '#':
                    if line[:3] == '###' and set:
                        self.sets.append(set)
                        set = None
                    if line.startswith("##sequence-region"):
                        splitline = line.split()
                        self.sequence_regions[splitline[1]] = line
                    #TODO: properly deal with comment lines.
                    self.sets.append(line)
                else:
                    line = GFF3_line(set,line)
                    #adding the feature individually
                    self.features_id[line.attributes.id] = line
                    if line.attributes.name:
                        if line.attributes.name in self.features_name:
                            #TODO: find a way to handle features that have the same name.
                            pass#print(line.attributes.id, line.attributes.name,  self.features_name[line.attributes.name].attributes.id)
                        else:
                            self.features_name[line.attributes.name] = line
                    #adding the set of features
                    if line.type == "region" and not line.attributes.parent:
                        #this feature has been deemed redundant and is not used in recent versions of the gff3,
                        if set:
                            #this is the first element of a set,
                            # old set needs to be added to the list and a new set created
                            self.sets.append(set)
                            set = GT_seq_location()
                        else:
                            set = GT_seq_location()
                            #if the set is none, it was also during init, and we need to set the owner_set again
                            line._owner_set = set
                        set._flanking_region = line
                    elif line.type == "flanking_region":
                        if set and set.flanking_region:
                            # this can also be the first element of a set,
                            # if the set already has a flanking region
                            # old set needs to be added to the list and a new set created
                            self.sets.append(set)
                            set = GT_seq_location()
                        else:
                            set = GT_seq_location()
                            #if the set is none, it was also during init, and we need to set the owner_set again
                            line._owner_set = set
                        set.flanking_region = line
                    elif line.type == "region" and line.attributes.parent:
                        set.gt_seq_region.append(line)
                    elif line.type == "PCR_product":
                        set.pcr_product.append(line)
                    elif line.type == "forward_primer":
                        set.forward_primer.append(line)
                    elif line.type == "reverse_primer":
                        set.reverse_primer.append(line)
                    elif line.type == "SNP":
                        set.snp.append(line)
                    else:
                        pass#print("line of type {} not added.".format(line.type))
            if set:
                # there was no '###' at the end of the file so the last set needs to be added.
                self.sets.append(set)

    def write_gff3(self,gff3_file=None):
        """
        write the gff3 document on stdout or given file
        :param gff3_file: path to write the gff3 file to
        :return: to doc or stdout
        """
        # write the new gff3
        if gff3_file:
            outfile = open(gff3_file, 'w')
        else:
            outfile = sys.stdout
        for set in self.sets:
            if isinstance(set, GT_seq_location):
                outfile.write(set.compose())
            else:
                outfile.write(set)
        outfile.close()

    def set_containing_position(self,seqid,position,position2=None):
        """
        find the primer set containing a certain position.
        :param seqid: string
        :param position: int
        :return: GT_seq_location
        """
        #TODO: have a way to find the position without looping all the sets
        for set in self.sets:
            if isinstance(set, GT_seq_location) and set.flanking_region.seqid == seqid:
                # -1 in case of counts starting at 0
                start = set.flanking_region.start-1
                end = set.flanking_region.end
                if start <= position and end >= position:
                    if position2: #if position 2 is present, it must also fall between these bounds
                        if start <= position2 and end >= position2:
                            return set
                    else:
                        return set
        #print("this position is not inside any of the targeted regions")
        return False

    def add_vcf_call(self,VCF_call):
        """
        find the group the snp belongs to
        establish the values needed for a new feature
        and add the feature to the SNP of the group.
        :return:
        """
        seqid = VCF_call.CHROM
        pos = VCF_call.POS
        allele = [x.__str__() for x in VCF_call.alleles]
        containing_set = self.set_containing_position(seqid,pos)
        if containing_set:
            new_snp = GFF3_line(containing_set)
            new_snp.is_new_snp(seqid,pos,allele)
            containing_set.add_snp(new_snp)
        else:
            print("this snp is not amplified in any target region of the primers")

    def give_active_snp(self):
        """
        loop on all the SNP and return those that are enabled
        :return: list of id with allele
        """
        genome_features = set()
        for key,value in self.features_id.items():
            if value.type == "SNP" and value.attributes.active:
                ref_seq = value.attributes.allele[0]
                if value.attributes.discovered and not value.attributes.validated:
                    # snp that have been added to the gff3 not used for designing the primers
                    # should have been validated before they can be called targets.
                    # this is a edge case that should never happen in practise.
                    pass
                elif all(len(ref_seq)==len(x) for x in value.attributes.allele):
                    # if all the alleles are the same length it is a SNP not a indel
                    genome_features.add(value.attributes.id + '_' + ref_seq)
                else:
                    genome_features.add(value.attributes.id +'_'+ ref_seq + '_indel')
        return genome_features

    def give_deactivated_snp(self):
        """
        loop on all the SNP and return those that are disabled
        :return: list of id with allele
        """
        genome_features = set()
        for key,value in self.features_id.items():
            if value.type == "SNP" and not value.attributes.active:
                ref_seq = value.attributes.allele[0]
                if all(len(ref_seq) == len(x) for x in value.attributes.allele):
                    genome_features.add(value.attributes.id + '_' + ref_seq)
                else:
                    genome_features.add(value.attributes.id + '_' + ref_seq + '_indel')
        return genome_features

    def disable_SNP(self,snp_id,reason,source="gff3_maniger"):
        """
        Put the SNP on inactive and check if this has consequences.
        :param snp_id: string with id or name of the SNP feature
        :param reason: motivation for the change
        :return:
        """
        snp_to_change = None
        date = datetime.datetime.now().strftime("%Y-%m-%d")
        # fetch the snp object
        if snp_id in self.features_id:
            snp_to_change = self.features_id[snp_id]
        # disable it if found
        if snp_to_change:
            snp_to_change.attributes.add_history(date,source,reason)
        else:
            pass#print("SNP with id {} was not found".format(snp_id))
        # see if now the primers should be removed.
        snp_set = snp_to_change._owner_set
        snp_set.are_all_snp_disabled()


def read_changes_tsv(tsv_file):
    """
    simple tab separated file with snp_name, reason of removal,
    :param tsv_file:
    :return: list with SNP to remove
    """
    changes = {}
    with open(tsv_file, 'r') as info_file:
        for info in info_file:
            split_info = info.strip().split('/t')
            changes[split_info[0]] = split_info[1]
    return changes


def main(args):
    # read the gff3 and put it in memory
    gff3 = gff3_file(args.gff3)
    # read the SNP to filter
    if args.info:
        changes = read_changes_tsv(args.info)
        for snp, reason in changes.items():
            gff3.disable_SNP(snp,reason)
    # write the new gff3
    gff3.write_gff3(args.output)

if __name__ == '__main__':
    args = argument_parser()
    main(args)