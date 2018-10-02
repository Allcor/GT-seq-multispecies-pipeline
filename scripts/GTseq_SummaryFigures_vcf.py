#!/usr/local/bin/python3.6
# GTseq_SummaryFigures_v3.py
# by Nate Campbell
# produce summary figures for GTseq libraries using GTseq_Genotyper_v3 output formatted files.
# Also outputs summary data in text format for further analysis.
"""
edited by Arlo Hoogeveen

Adds a plot from vcf file where penos also has data. Will make plot of undiscovered SNP as well
"""

import math
import argparse
import vcf
from os import listdir
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


def argument_parse():
    parser = argparse.ArgumentParser(description='Options for BarcodeSplit script')
    parser.add_argument(
        "-e", "--genos",
        dest="genos_dir",
        help="the path to directory containing .genos files for library"
    )
    parser.add_argument(
        "-g", "--gff3",
        dest="gff_file",
        help="the path to the gff3 file with information on SNP locations and primer targets."
    )
    parser.add_argument(
        "-f", "--filter",
        dest="filter_mode",
        default=2,
        help="filter modes, all features:0, unfiltered:1, unfiltered and new SNP:2, new SNP only:3"
    )
    parser.add_argument(
        "-v", "--vcf",
        dest="vcf_file",
        help="the path to the vcf file containing varients of the GT-seq run."
    )
    parser.add_argument(
        "-n", "--name",
        dest="library_name",
        help="type the library name"
    )
    parser.add_argument(
        "-a", "--all",
        dest="all_plots",
        action='store_true',
        help="Print separate plots for all the assays",
    )
    return parser.parse_args()


def read_genome_features(gff_file):
    genome_features = {}
    with open(gff_file, 'r') as gff:
        keys = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        for line in gff:
            split_line = line.strip().split('\t')
            if split_line[0][0] != '#' and len(split_line) == 9:
                #if this is not working, it's not a gff3 file.
                gff_line = {keys[i]:x for i,x in enumerate(split_line)}
                if gff_line['type'] == 'SNP':
                    loci_name = gff_line['seqid']+'_'+gff_line['start']
                    attr = {value.split('=')[0]:value.split('=')[1] for value in gff_line['attributes'].split(';')}
                    # TODO: if the filter changes this might need to be adjusted.
                    if attr['filtered'] == 'False':
                        attr['filtered'] = None
                    genome_features[loci_name] = attr
    return genome_features


def make_assay_list(flist,assaylist):
    #open top file and create assay list...initialize dictionary of loci with their corresponding percentage of on-target reads
    inds = float(len(flist))
    OT_Dict = {}
    StDEV_Dict = {}
    OTP_Dict = {}
    StDEV2_Dict = {}
    with open(flist[0]) as f:
        lineNo = 0
        for line in f:
            lineNo += 1
            if lineNo > 1:
                stuff = line.split(',')
                assaylist.append(stuff[0])
                OT_Dict[stuff[0]] = float(0)
                StDEV_Dict[stuff[0]] = float(0)
                OTP_Dict[stuff[0]] = float(0)
                StDEV2_Dict[stuff[0]] = float(0)

    for loci in assaylist:
        for genos in flist:
            with open(genos) as g:
                for line in g:
                    info = line.split(',')
                    if loci in info[0] and len(info[0]) == len(loci):
                        OT_Dict[info[0]] = OT_Dict[info[0]] + float(info[10])
                        OTP_Dict[info[0]] = OTP_Dict[info[0]] + float(info[9])

    #Get the mean of percenage of OT reads for each locus and mean of percentage
    #forward primer reads containing probe sequences for each locus...
    for loci in assaylist:
        OT_Dict[loci] = OT_Dict[loci] / inds
        OTP_Dict[loci] = OTP_Dict[loci] / inds

    #Calculate standard deviations at each locus...
    for loci in assaylist:
        for genos in flist:
            with open(genos) as g:
                for line in g:
                    info = line.split(',')
                    if loci in info[0] and len(info[0]) == len(loci):
                        variance = (float(info[10]) - OT_Dict[loci])**2
                        variance2 = (float(info[9]) - OTP_Dict[loci])**2
                        StDEV_Dict[loci] = StDEV_Dict[loci] + variance
                        StDEV2_Dict[loci] = StDEV2_Dict[loci] + variance2

    for loci in assaylist:
        StDEV_Dict[loci] = (math.sqrt(StDEV_Dict[loci] / inds))
        StDEV2_Dict[loci] = (math.sqrt(StDEV2_Dict[loci] / inds))

    return (OT_Dict, StDEV_Dict, OTP_Dict, StDEV2_Dict, assaylist)


def make_assay_list_vcf(vcf_file):
    """
    ##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
    ##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">
    ##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">
    ##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
    ##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
    ##INFO=<ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
    ##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
    ##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
    ##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">
    ##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
    ##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">
    ##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
    ##FORMAT=<ID=ADF,Number=R,Type=Integer,Description="Allelic depths on the forward strand">
    ##FORMAT=<ID=ADR,Number=R,Type=Integer,Description="Allelic depths on the reverse strand">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##INFO=<ID=ICB,Number=1,Type=Float,Description="Inbreeding Coefficient Binomial test (bigger is better)">
    ##INFO=<ID=HOB,Number=1,Type=Float,Description="Bias in the number of HOMs number (smaller is better)">
    ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
    ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
    ##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
    ##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
    :param vcf_file:
    :return:
    """
    assaylist = []
    OTP_sample_dict = {}
    OffTarget_sample_dict = {}
    OTP_totals = {}
    vcf_reader = vcf.Reader(filename=vcf_file)
    # getting the on target reads for each site
    for record in vcf_reader:
        site_id = record.CHROM+'_'+str(record.POS)
        if record.num_called != 0:
            if site_id not in assaylist:
                assaylist.append(site_id)
            OTP_sample_dict[site_id] = {call.sample:call.data.DP for call in record.samples}
            OffTarget_sample_dict[site_id] = {call.sample:sum(call.data.AD[2:]) for call in record.samples}
    # total on target reads on all sites (for each sample)
    for sample in vcf_reader.samples:
        total_OTP = 0
        for key,value in OTP_sample_dict.items():
            total_OTP += value[sample]
        OTP_totals[sample] = total_OTP

    OT_Dict = {}
    OTP_Dict = {}
    StDEV_Dict = {}
    StDEV2_Dict = {}
    inds = 0
    # calculating the average on target rates
    for site in assaylist:
        OT = 0.0
        OTP = 0.0
        items = 0
        for key,value in OTP_sample_dict[site].items():
            if value != 0:
                items += 1
                OT += value / OTP_totals[key] * 100
                OTP += (value-OffTarget_sample_dict[site][key]) / value
        OT_Dict[site] = OT/items
        OTP_Dict[site] = OTP/items

        StDEV = 0.0
        StDEV2 = 0.0
        for key, value in OTP_sample_dict[site].items():
            if value != 0:
                OT = value / OTP_totals[key] * 100
                OTP = (value-OffTarget_sample_dict[site][key]) / value
                variance = (OT - OT_Dict[site])**2
                variance2 = (OTP - OTP_Dict[site])**2
                StDEV += variance
                StDEV2 += variance2
        StDEV_Dict[site] = StDEV / items
        StDEV2_Dict[site] = StDEV2 / items
        if items > inds:
            inds = items

    return (OT_Dict, StDEV_Dict, OTP_Dict, StDEV2_Dict, assaylist, inds)

def make_sorted_assay_list(OT_Dict,OTP_Dict,StDEV_Dict,StDEV2_Dict):
    # Get sorted list of standard deviations for error bars...
    Sorted_OT = sorted(OT_Dict.values())
    Sorted_OTkeys = sorted(OT_Dict, key=OT_Dict.get)
    Sorted_stDEV = []
    Sorted_OTP = sorted(OTP_Dict.values())
    Sorted_OTPkeys = sorted(OTP_Dict, key=OTP_Dict.get)
    Sorted_stDEV2 = []

    for loci in Sorted_OTkeys:
        Sorted_stDEV.append(StDEV_Dict[loci])
    for loci in Sorted_OTPkeys:
        Sorted_stDEV2.append(StDEV2_Dict[loci])

    return(Sorted_OT, Sorted_OTkeys, Sorted_stDEV, Sorted_OTP, Sorted_OTPkeys, Sorted_stDEV2)


def make_locus_graph(f_out, figID, subpos, assaylist, OT_Dict, StDEV_Dict):
    #plot read distribution bar graph using means and standard deviation for each locus...
    plt.figure(figID)
    bar_width = 1
    left = 0
    AssayNum = float(len(assaylist))
    L_AvOTP = 100 / AssayNum

    plt.subplot(subpos)
    for loci in assaylist:
        plt.bar(left, OT_Dict[loci], width=bar_width, bottom=None, hold=None,
        edgecolor='b', yerr=StDEV_Dict[loci], color='b', error_kw=dict(ecolor='black', lw=0.5,
        alpha=0.5, capsize=0.5))
        print(loci, OT_Dict[loci], StDEV_Dict[loci])
        f_out.write(loci + '\t' + str(OT_Dict[loci]) + '\t' + str(StDEV_Dict[loci]) + '\n')
        left = left + bar_width

    plt.plot([0, left], [L_AvOTP, L_AvOTP], 'r-', linewidth=1.0)
    plt.xlabel('Loci', fontsize=8)
    plt.ylabel('% of On-Target Reads', fontsize=8)
    plt.title('Read distribution (unsorted)', fontsize=10)


def make_percentage_graph(f_out, figID, subpos, assaylist, OTP_Dict, StDEV2_Dict):
    #Create bar graph of percentage (reads with forward primer AND probe / reads with fwd primer)...
    plt.figure(figID)
    bar_width = 1
    left = 0

    plt.subplot(subpos)
    for loci in assaylist:
        plt.bar(left, OTP_Dict[loci], width=bar_width, bottom=None, hold=None,
        edgecolor='b', yerr=StDEV2_Dict[loci], color='b', error_kw=dict(ecolor='black', lw=0.5,
        alpha=0.5, capsize=0.5))
        print(loci, OTP_Dict[loci], StDEV2_Dict[loci])
        f_out.write(loci + '\t' + str(OTP_Dict[loci]) + '\t' + str(StDEV2_Dict[loci]) + '\n')
        left = left + bar_width

    plt.xlabel('Loci', fontsize=8)
    plt.ylabel('% On-Target Primers', fontsize=8)
    plt.title('Primers On-Target (unsorted)', fontsize=10)


def make_value_graph(f_out, figID, subpos, assaylist, Sorted_OT, Sorted_OTkeys, Sorted_stDEV):
    #populate read distribution bar graph using sorted data...
    plt.figure(figID)
    bar_width = 1
    left = 0
    AssayNum = float(len(assaylist))
    L_AvOTP = 100 / AssayNum
    Assays = int(AssayNum)

    plt.subplot(subpos)
    for x in range(0, Assays):
        plt.bar(left, Sorted_OT[x], width=bar_width, bottom=None, hold=None,
        edgecolor='green', yerr=Sorted_stDEV[x], color='green', error_kw=dict(ecolor='black', lw=0.5, alpha=0.5, capsize=0.5))
        print(Sorted_OTkeys[x], Sorted_OT[x], Sorted_stDEV[x])
        f_out.write(Sorted_OTkeys[x] + '\t' + str(Sorted_OT[x]) + '\t' + str(Sorted_stDEV[x]) + '\n')
        left = left + bar_width

    plt.plot([0, left], [L_AvOTP, L_AvOTP], 'r-', linewidth=1.0)
    plt.xlabel('Loci', fontsize=8)
    plt.title('Read distribution (Sorted)', fontsize=10)


def make_sorted_percentage_graph(f_out, figID, subpos, assaylist, Sorted_OTP, Sorted_OTPkeys, Sorted_stDEV2):
    #Create sorted bar graph of percentage (reads with forward primer AND probe / reads with fwd primer)...
    plt.figure(figID)
    left = 0
    bar_width = 1
    AssayNum = float(len(assaylist))
    Assays = int(AssayNum)

    plt.subplot(subpos)
    for x in range(0, Assays):
        plt.bar(left, Sorted_OTP[x], width=bar_width, bottom=None, hold=None,
        color='red', edgecolor='red', yerr=Sorted_stDEV2[x], error_kw=dict(ecolor='black', lw=0.5,
        alpha=0.5, capsize=0.5))
        print(Sorted_OTPkeys[x], Sorted_OTP[x], Sorted_stDEV2[x])
        f_out.write(Sorted_OTPkeys[x] + '\t' + str(Sorted_OTP[x]) + '\t' + str(Sorted_stDEV2[x]) + '\n')
        left = left + bar_width

    plt.xlabel('Loci', fontsize=8)
    plt.title('Primers On-Target (sorted)', fontsize=10)


def make_loci_plot(file):
    lineNo = 0
    g = open(file)
    for line in g:
        lineNo = lineNo + 1
        if lineNo > 1:
            info = line.split(',')
            xarr = info[1].split('=')
            yarr = info[2].split('=')
            x = int(round(float(xarr[1])))
            y = int(round(float(yarr[1])))
            ratio = float(info[3])
            sum_xy = x + y
            scale = 30.0
            if sum_xy < 10:
                color = 'yellow'
            elif sum_xy == 0:
                color = 'yellow'
            elif ratio > 10:
                color = 'red'
            elif ratio < 0.1:
                color = 'blue'
            elif ratio < 2 and ratio > 0.5:
                color = 'purple'
            else:
                color = 'yellow'
            plt.subplot(223)
            plt.scatter(x, y, c='none', s=scale, label=color,alpha=0.4, edgecolors=color)
            plt.subplot(224)
            plt.scatter(x, y, c='none', s=scale, label=color,alpha=0.4, edgecolors=color)


def genos_plot(flist,loci,f_out,plot221,plot222):
    """
    Makes plot for sequence variant data coming from nate's genos files.
    :param flist: list with genos files
    :param loci: the sequence variant to make a plot from
    :param f_out: the logfile
    :param plot221: the plot with counts
    :param plot222: the plot with percentage
    :return: changed plot221 and plot222
    """
    xmax = 0
    ymax = 0
    fmax = 0
    fmax2 = 0
    A1_corr = float(0)
    A2_corr = float(0)
    gt_per = float(0)
    gt_inds = float(0)
    inds = float(len(flist))
    for genos in flist:
        g = open(genos)
        for line in g:
            info = line.split(',')
            if loci in info[0] and len(info[0]) == len(loci):
                size = len(info)
                if size == 13:
                    A1_corr = float(info[6])
                    A2_corr = float(info[7])
                xarr = info[1].split('=')
                yarr = info[2].split('=')
                x = int(round(float(xarr[1])))
                y = int(round(float(yarr[1])))
                ratio = float(info[3])
                sum_xy = x + y
                AF_div = float(xarr[1]) + float(yarr[1])
                if AF_div == 0:
                    AF_div = 0.1
                p_A2 = float(yarr[1])/AF_div*100
                scale = 50.0
                if x > xmax:
                    xmax = x
                if y > ymax:
                    ymax = y
                if sum_xy > fmax2:
                    fmax2 = sum_xy
                if sum_xy < 10:
                    color = 'yellow'
                elif sum_xy == 0:
                    color = 'yellow'
                elif ratio > 10:
                    color = 'red'
                    gt_inds = gt_inds + 1
                elif ratio < 0.1:
                    color = 'blue'
                    gt_inds = gt_inds + 1
                elif ratio < 2 and ratio > 0.5:
                    color = 'purple'
                    gt_inds = gt_inds + 1
                else:
                    color = 'yellow'
                plot221.scatter(x, y, c=color, s=scale, label=color,alpha=0.4, edgecolors='none')
                plot222.scatter(sum_xy, p_A2, c=color, s=scale, label=color,alpha=0.4, edgecolors='none')
    if xmax > ymax:
        fmax = xmax
    else:
        fmax = ymax
    gt_per = gt_inds / inds * 100
    gt_per = str(gt_per)
    A1_corr = str(A1_corr)
    A2_corr = str(A2_corr)
    text = ' : Corrections [' + A1_corr[:3] + ', ' + A2_corr[:3] + '] : GT% = ' + gt_per[:5]
    print(loci + text)
    f_out.write(loci + text + '\n')
    #Create first subplot XY scatter A1 vs A2 counts...
    plot221.grid(True)
    plot221.axis([-5, fmax, -5, fmax])
    plot221.plot([0, 10], [10, 0], 'y-', linewidth=2.0)
    plot221.plot([9, 10000], [1, 1000], 'r-', linewidth=2.0)
    plot221.plot([1000, 1], [10000, 9], 'b-', linewidth=2.0)
    plot221.plot([5, 10000], [5, 10000], 'm-', linewidth=2.0)
    plot221.plot([6.6, 10000], [3.3, 5000], 'k-', linewidth=2.0)
    plot221.plot([3.3, 5000], [6.6, 10000], 'k-', linewidth=2.0)
    plot221.set_title('A1 vs A2 counts (nate)', fontsize=10)
    plot221.set_xlabel('A1 counts', fontsize=8)
    plot221.set_ylabel('A2 counts', fontsize=8)

    #Create second suplot % A2 counts...
    plot222.axis([-5, fmax2, -5, 105])
    plot222.plot([10, 10], [0, 100], 'y-', linewidth=2.0)
    plot222.plot([0, 10000], [10, 10], 'r-', linewidth=2.0)
    plot222.plot([0, 10000], [90, 90], 'b-', linewidth=2.0)
    plot222.plot([10, 10000], [66.6, 66.6], 'k-', linewidth=2.0)
    plot222.plot([10, 10000], [33.3, 33.3], 'k-', linewidth=2.0)
    plot222.set_title('% allele 2 (nate)', fontsize=10)
    plot222.set_xlabel('Counts', fontsize=8)
    plot222.set_ylabel('% A2 counts', fontsize=8)

def make_locus_plots(loci,fig_num,f_out,assaylist,Sorted_OTkeys,Sorted_OT,Sorted_stDEV,Sorted_OTPkeys,Sorted_OTP,Sorted_stDEV2,flist,Lname):
    xmax = 0
    ymax = 0
    fmax = 0
    fmax2 = 0
    A1_corr = float(0)
    A2_corr = float(0)
    gt_per = float(0)
    gt_inds = float(0)
    # this is incorrect, reflects the read groups, this now contains all the species from the library.
    # inds = len(loci.samples)
    # ugly fix, if there are no reads (for this locus) don't count it.
    inds = len([x for x in loci.samples if x.data.DP != 0])
    if inds == 0:
        inds = 1
    # TODO, should be solved by not adding those samples in the sam header
    site_id = loci.CHROM + '_' + str(loci.POS)
    if flist != [] and site_id in assaylist:
        figsize = [6.4,6.8]
        subadjust = {'top':0.84, 'wspace':0.3, 'hspace':0.4}
        plot221 = plt.subplot(321)
        plot222 = plt.subplot(322)
        genos_plot(flist, site_id, f_out, plot221, plot222)
        plot221 = plt.subplot(323)
        plot222 = plt.subplot(324)
        plot223 = plt.subplot(325)
        plot224 = plt.subplot(326)
    else:
        figsize = [6.4,4.8]
        subadjust = {'top':0.8, 'wspace':0.3, 'hspace':0.4}
        plot221 = plt.subplot(221)
        plot222 = plt.subplot(222)
        plot223 = plt.subplot(223)
        plot224 = plt.subplot(224)
    for sample in loci.samples:
        #size = len(info)
        #if size == 13:
        #    A1_corr = float(info[6])
        #    A2_corr = float(info[7])
        #xarr = info[1].split('=')
        if isinstance(sample.data.AD,list):
            xarr = [loci.alleles[0], sample.data.AD[0]]
            try:
                yarr = [loci.alleles[1].sequence,sample.data.AD[1]]
            except AttributeError:
                yarr = [loci.alleles[1].type, sample.data.AD[1]]
        else:
            xarr = [loci.alleles[0],sample.data.AD]
            yarr = ['.',0]
        #yarr = info[2].split('=')
        x = int(round(float(xarr[1])))
        y = int(round(float(yarr[1])))
        if yarr[1] != 0:
            ratio = float(xarr[1]/yarr[1])
        else:
            ratio = xarr[1]
        sum_xy = x + y
        AF_div = float(xarr[1]) + float(yarr[1])
        if AF_div == 0:
            AF_div = 0.1
        p_A2 = float(yarr[1])/AF_div*100
        scale = 50.0
        if x > xmax:
            xmax = x
        if y > ymax:
            ymax = y
        if sum_xy > fmax2:
            fmax2 = sum_xy
        if sum_xy < 10:
            color = 'yellow'
        elif sum_xy == 0:
            color = 'yellow'
        elif ratio > 10:
            color = 'red'
            gt_inds = gt_inds + 1
        elif ratio < 0.1:
            color = 'blue'
            gt_inds = gt_inds + 1
        elif ratio < 2 and ratio > 0.5:
            color = 'purple'
            gt_inds = gt_inds + 1
        else:
            color = 'yellow'
        plot221.scatter(x, y, c=color, s=scale, label=color,alpha=0.4, edgecolors='none')
        plot222.scatter(sum_xy, p_A2, c=color, s=scale, label=color,alpha=0.4, edgecolors='none')
    if xmax > ymax:
        fmax = xmax
    else:
        fmax = ymax
    gt_per = gt_inds / inds * 100
    gt_per = str(gt_per)
    A1_corr = str(A1_corr)
    A2_corr = str(A2_corr)
    text = ' : Corrections [' + A1_corr[:3] + ', ' + A2_corr[:3] + '] : GT% = ' + gt_per[:5]
    print(site_id + text)
    f_out.write(site_id + text + '\n')
    #Create first subplot XY scatter A1 vs A2 counts...
    plot221.grid(True)
    plot221.axis([-5, fmax, -5, fmax])
    plot221.plot([0, 10], [10, 0], 'y-', linewidth=2.0)
    plot221.plot([9, 10000], [1, 1000], 'r-', linewidth=2.0)
    plot221.plot([1000, 1], [10000, 9], 'b-', linewidth=2.0)
    plot221.plot([5, 10000], [5, 10000], 'm-', linewidth=2.0)
    plot221.plot([6.6, 10000], [3.3, 5000], 'k-', linewidth=2.0)
    plot221.plot([3.3, 5000], [6.6, 10000], 'k-', linewidth=2.0)
    plot221.set_title('A1 vs A2 counts', fontsize=10)
    plot221.set_xlabel('A1 counts', fontsize=8)
    plot221.set_ylabel('A2 counts', fontsize=8)

    #Create second suplot % A2 counts...
    plot222.axis([-5, fmax2, -5, 105])
    plot222.plot([10, 10], [0, 100], 'y-', linewidth=2.0)
    plot222.plot([0, 10000], [10, 10], 'r-', linewidth=2.0)
    plot222.plot([0, 10000], [90, 90], 'b-', linewidth=2.0)
    plot222.plot([10, 10000], [66.6, 66.6], 'k-', linewidth=2.0)
    plot222.plot([10, 10000], [33.3, 33.3], 'k-', linewidth=2.0)
    plot222.set_title('% allele 2', fontsize=10)
    plot222.set_xlabel('Counts', fontsize=8)
    plot222.set_ylabel('% A2 counts', fontsize=8)

    #continue with bar graph subplots...
    AssayNum = float(len(assaylist))
    Assays = int(AssayNum)
    L_AvOTP = 100 / AssayNum
    bar_width = 1
    left = 0
    for x in range(0, Assays):
        bar_color='grey'
        ebar_color='grey'
        opacity=0.5
        if site_id in Sorted_OTkeys[x] and len(site_id) == len(Sorted_OTkeys[x]):
            bar_color='red'
            ebar_color='black'
            opacity=1
            text = text + '\nPercentage of On-Target Reads = ' + str(Sorted_OT[x])[:4]
        plot223.bar(left, Sorted_OT[x], width=bar_width, bottom=None,
            alpha=opacity, yerr=Sorted_stDEV[x], color=bar_color, edgecolor=bar_color,
            error_kw=dict(ecolor=ebar_color, lw=0.5, alpha=opacity, capsize=0.5))
        left = left + bar_width

    plot223.plot([0, left], [L_AvOTP, L_AvOTP], 'r-', linewidth=1.0)
    plot223.set_xlabel('Loci', fontsize=8)
    plot223.set_ylabel('Percentage of On-Target Reads', fontsize=8)
    plot223.set_title('Read distribution among loci', fontsize=10)

    left = 0
    for x in range(0, Assays):
        bar_color='grey'
        ebar_color='grey'
        opacity=0.5
        if site_id in Sorted_OTPkeys[x] and len(site_id) == len(Sorted_OTPkeys[x]):
            bar_color='red'
            ebar_color='black'
            opacity=1
            text = text + '\nOn-target Primer = ' + str(Sorted_OTP[x])[:4]
        plot224.bar(left, Sorted_OTP[x], width=bar_width, bottom=None,
            alpha=opacity, yerr=Sorted_stDEV2[x], color=bar_color, edgecolor=bar_color,
            error_kw=dict(ecolor=ebar_color, lw=0.5, alpha=opacity, capsize=0.5))
        left = left + bar_width

    plot224.set_xlabel('Loci', fontsize=8)
    plot224.set_ylabel('Percentage On-Target Primers', fontsize=8)
    plot224.set_title('Primers On-Target', fontsize=10)

    if loci.FILTER:
        text += "\nlocus has filter: {}".format(loci.FILTER)
        #text += "\nLocation was not called"

    plt.figure(1)
    plt.suptitle(site_id + text, fontsize=12)
    fig = plt.gcf()
    fig.set_size_inches(figsize)
    fig.subplots_adjust(**subadjust)
    if loci.is_snp:
        filename = "{}_{}.png".format(Lname,site_id)
    else:
        filename = "{}_{}.{}.png".format(Lname,site_id,loci.var_type)
    fig.savefig(filename.format(fig_num), format='png', dpi=300)
    plt.clf()


def main(args):

    # ask for sequence file if not given.
    if args.library_name:
        Lname = args.library_name
    else:
        print('type the library name *use single quotes*')
        Lname = input()

    fout_name = Lname + '_SummaryData.txt'
    f_out = open(fout_name, 'w')
    f_out.write(Lname)
    f_out.write('GTseq Summary Data\n\n')

    assaylist = []
    flist = []

    if args.genos_dir:
        path = args.genos_dir
        # filter for .genos files in directory and add to "flist"...
        list1 = listdir(path)
        for i in list1:
            if '.genos' in i:
                j = path + '/' + i
                flist.append(j)
        # open top file and create assay list...initialize dictionary of loci with their corresponding percentage of on-target reads
        (OT_Dict, StDEV_Dict, OTP_Dict, StDEV2_Dict, assaylist) = make_assay_list(flist, assaylist)
    else:
        (OT_Dict, StDEV_Dict, OTP_Dict, StDEV2_Dict, assaylist, inds) = make_assay_list_vcf(args.vcf_file)

    (Sorted_OT, Sorted_OTkeys, Sorted_stDEV, Sorted_OTP, Sorted_OTPkeys, Sorted_stDEV2) = make_sorted_assay_list(OT_Dict, OTP_Dict, StDEV_Dict, StDEV2_Dict)

    #### making the summery plots ####

    plt.rc('font', **{'size': 5})

    text='Read Distribution data (sorted by locus name)'
    print(text+'\n')
    f_out.write(text+'\n')
    make_locus_graph(f_out, 1, 221, assaylist, OT_Dict, StDEV_Dict)

    text='Read Distribution data (sorted by value)'
    print(text+'\n')
    f_out.write('\n'+text+'\n')
    make_value_graph(f_out, 1, 222, assaylist, Sorted_OT, Sorted_OTkeys, Sorted_stDEV)

    text='Primer Reads On-Target (reads with forward primer AND probe / reads with fwd primer)*100'
    print(text+'\n')
    f_out.write('\n'+text+'\n')
    make_percentage_graph(f_out, 1, 223, assaylist, OTP_Dict, StDEV2_Dict)

    text='Primer Reads On-Target (sorted) '
    text+='(reads with forward primer AND probe / reads with fwd primer)*100'
    print(text+'\n')
    f_out.write('\n'+text+'\n')
    make_sorted_percentage_graph(f_out, 1, 224, assaylist, Sorted_OTP, Sorted_OTPkeys, Sorted_stDEV2)

    plt.figure(1)
    plt.suptitle(Lname.split('/')[-1], fontsize=12)
    plt.subplots_adjust(top=0.9, wspace=0.3, hspace=0.4)
    plt.savefig(Lname + '_summary.png', dpi=300, format='png')
    plt.clf()

    #### making the summary library ####
    """
    print('Library Summary\n')
    f_out.write('\nLibrary Summary\n')
    #open each .genos file and define x and y coordinates for plotting Raw-reads vs. GT% graph...

    # Initialize variables...
    AssayNum = float(len(assaylist))
    xmax = 0
    xmax2 = 0
    num90 = float(0)
    aveOTP = float(0)
    for genos in flist:
        g = open(genos)
        genper = float(0)
        otreads = 0
        rawreads = 0
        NA = 0
        lineNo2 = 0
        for line in g:
            lineNo2 = lineNo2 + 1
            if lineNo2 == 1:
                info = line.split(',')
                RRarr = info[1].split(':')
                OTarr = info[2].split(':')
                OTParr = info[3].split(':')
                otreads = int(OTarr[1]) / 1000
                rawreads = int(RRarr[1]) / 1000
                aveOTP = aveOTP + float(OTParr[1])
            elif lineNo2 > 1:
                info2 = line.split(',')
                if 'NA' in info2[5]:
                    NA = NA + 1
        if rawreads > xmax:
            xmax = rawreads
        if otreads > xmax2:
            xmax2 = otreads
        NA = float(NA)
        genper = (1 - (NA / AssayNum)) * 100
        scale = 50.0
        if genper >= 90:
            num90 = num90 + 1
        # print(genos,otreads,genper)
        plt.subplot(221)
        plt.scatter(rawreads, genper, c='orange', s=scale, label='black', alpha=0.4, edgecolors='none')
        plt.subplot(222)
        plt.scatter(otreads, genper, c='black', s=scale, label='black', alpha=0.4, edgecolors='none')
    
    per90 = num90 / inds * 100
    aveOTP = aveOTP / inds

    perSTR = str(per90)
    aveOTPSTR = str(aveOTP)

    text1 = 'Samples in library: ' + str(int(inds))
    text2 = 'Samples over 90% GT: ' + str(int(num90))
    text3 = 'Percentage over 90% GT: ' + perSTR[:4] + '%'
    text4 = 'Average OT-Percentage: ' + aveOTPSTR[:4] + '%'
    print(text1, text2, text3, text4)
    f_out.write(text1 + '\n')
    f_out.write(text2 + '\n')
    f_out.write(text3 + '\n')
    f_out.write(text4 + '\n')

    font = {
        'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 10,
    }

    plt.figure(1)
    plt.subplot(221)
    plt.plot([0, xmax], [90, 90], 'r-', linewidth=2.0)
    plt.grid(True)
    plt.axis([-5, xmax, -5, 105])
    plt.xlabel('Raw Reads (K)')
    plt.ylabel('Genotyping Percentage')
    plt.text(0, 30, text1, fontdict=font)
    plt.text(0, 20, text2, fontdict=font)
    plt.text(0, 10, text3, fontdict=font)
    plt.text(0, 0, text4, fontdict=font)

    plt.subplot(222)
    plt.plot([0, xmax2], [90, 90], 'r-', linewidth=2.0)
    plt.grid(True)
    plt.axis([-5, xmax2, -5, 105])
    plt.xlabel('On-Target Reads (K)')
    plt.ylabel('Genotyping Percentage')
    plt.text(0, 30, text1, fontdict=font)
    plt.text(0, 20, text2, fontdict=font)
    plt.text(0, 10, text3, fontdict=font)
    plt.text(0, 0, text4, fontdict=font)

    print('Per-locus Summary\n')
    f_out.write('\nPer-locus Summary\n')

    #Gather data from all or first 100 .genos files and plot data from all loci onto a single graph...
    end = 100
    if inds < 100:
         end = int(inds)
    print('Plotting all loci for %s samples...\nKeepin it 100...\nLots of data points\nThis will take a couple minutes\n' % end)
    
    for i in range(0, end):
        make_loci_plot(flist[i])
    
    plt.subplot(223)
    plt.grid(True)
    plt.axis([-5, 500, -5, 500])
    plt.plot([0, 10], [10, 0], 'y-', linewidth=2.0)
    plt.plot([9, 10000], [1, 1000], 'r-', linewidth=2.0)
    plt.plot([1000, 1], [10000, 9], 'b-', linewidth=2.0)
    plt.plot([5, 10000], [5, 10000], 'm-', linewidth=2.0)
    plt.plot([6.6, 10000], [3.3, 5000], 'k-', linewidth=2.0)
    plt.plot([3.3, 5000], [6.6, 10000], 'k-', linewidth=2.0)
    plt.title('All Loci Plot')
    plt.xlabel('A1 counts')
    plt.ylabel('A2 counts')

    plt.subplot(224)
    plt.grid(True)
    plt.axis([-5, 100, -5, 100])
    plt.plot([0, 10], [10, 0], 'y-', linewidth=2.0)
    plt.plot([9, 10000], [1, 1000], 'r-', linewidth=2.0)
    plt.plot([1000, 1], [10000, 9], 'b-', linewidth=2.0)
    plt.plot([5, 10000], [5, 10000], 'm-', linewidth=2.0)
    plt.plot([6.6, 10000], [3.3, 5000], 'k-', linewidth=2.0)
    plt.plot([3.3, 5000], [6.6, 10000], 'k-', linewidth=2.0)
    plt.title('All Loci Plot (Zoom)')
    plt.xlabel('A1 counts')
    plt.ylabel('A2 counts')
    plt.figure(1)
    plt.suptitle(Lname, fontsize=12)
    plt.savefig(Lname + '_002.png', dpi=300, format='png')
    plt.clf()
    """

    #Create summary data for each locus including allele counts scatter plot, allele frequency scatter plot, % of OT reads bar graph
    #with locus bar highlighted, and % of reads on target with locus bar highlighted...
    plt.figure(1)
    fig_num = 2

    if args.all_plots:
        print('Done...\n\nNow creating plots for each locus...\n')
        vcf_reader = vcf.Reader(filename=args.vcf_file)
        #read the positions from the gff3
        if args.gff_file:
            filtermode = args.filter_mode
            known_sites = read_genome_features(args.gff_file)
        else:
            filtermode = 3
            known_sites = {}
        # initialising the check
        for value in known_sites.values():
            value['in_vcf'] = False
        for loci in vcf_reader:
            loci_name = loci.CHROM+'_'+str(loci.POS)
            if loci_name in known_sites:
                # location present, this can ether be filtered or not
                snp_info = known_sites[loci_name]
                filtered = snp_info['filtered']
                if filtermode == 0:
                    # all locations in the feature file
                    fig_num = fig_num + 1
                    make_locus_plots(loci, fig_num, f_out, assaylist, Sorted_OTkeys, Sorted_OT, Sorted_stDEV,
                                     Sorted_OTPkeys, Sorted_OTP, Sorted_stDEV2, flist, Lname)
                elif filtermode <= 2 and not filtered:
                    # only unfiltered locations in the feature file
                    fig_num = fig_num + 1
                    make_locus_plots(loci, fig_num, f_out, assaylist, Sorted_OTkeys, Sorted_OT, Sorted_stDEV,
                                     Sorted_OTPkeys, Sorted_OTP, Sorted_stDEV2, flist, Lname)
                known_sites[loci_name]['in_vcf'] = True
            elif filtermode >= 2: #and not loci.FILTER
                # location not present in feature file but has valid SNP
                fig_num = fig_num + 1
                make_locus_plots(loci, fig_num, f_out, assaylist, Sorted_OTkeys, Sorted_OT, Sorted_stDEV,
                                 Sorted_OTPkeys, Sorted_OTP, Sorted_stDEV2, flist, Lname)
        missed_SNP = [key for key,value in known_sites.items() if value['in_vcf'] == False]
        if known_sites != {}:
            print("figures for {} were not made.".format(', '.join(missed_SNP)))

    f_out.close()
    summary = 'All done.\nSummary data: %s\nSummary figures: are saved in the current folder'
    print(summary % (fout_name))
    print('Convert figures to .pdf using the convert command **$ convert *png LXXXX.pdf**')

if __name__ == '__main__':
    args = argument_parse()
    main(args)