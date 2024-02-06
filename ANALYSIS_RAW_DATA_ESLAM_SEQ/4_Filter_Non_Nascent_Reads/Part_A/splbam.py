#! /usr/bin/env python3 -u


"""Split a BAM file containing alignments
from a SLAM/TUC/TL-seq experiment into
labelled/unlabelled BAM files, to use as input
for PulseR (counts). Additionally outputs info
on mismatch rates, etc.

NOTE: We use the notation labelled/unlabelled, but
      new/old (RNA) would be less ambiguous, cf.
      biochemical separation.

NOTE: Hard coded library type RF-FIRSTSTRAND!

NOTE: Currently using GRAND-SLAM snpdata - default
      format (P Value = 0.001)
"""

import os
import sys
import logging
import argparse
import shlex
import pandas as pd
import numpy as np
import pysam as ps

from tabulate import tabulate

from collections import defaultdict

import utils

logger = logging.getLogger(__name__)


default_num_cpus = 48
default_mem = '400G'


def get_base_vec(base, strand):
    # report mismatch following JACUSA output - w/o coverage
    base = str(base).upper()
    vec=['A', 'C', 'G', 'T', base]
    vec.sort()
    idx = np.array([0]*4)
    idx[vec.index(base)] = 1
    if strand == '-':
        idx = idx[::-1]
    return ','.join(idx.astype(str))


def get_mismatches(contig, bam_file):
    # Find all mismatches
    bam = ps.AlignmentFile(bam_file, "rb").fetch(contig=contig)
    bed_entries = []
    # also add info on all mismatches
    mismatch_details_dict = defaultdict(int)

    for r in bam:
        query_id = r.query_name
        seq = r.query_alignment_sequence
        flag = r.flag
        qstart = r.query_alignment_start
        rlen = r.query_length
        # sort strand, based on read pair for library ** RF-FIRSTSTRAND **
        strand = '+'
        if (r.is_read1 and not r.is_reverse) or (r.is_read2 and r.is_reverse):
            strand = '-'
        # pre-filtering based on some featureCounts selected params
        # this is only for estimation of actual conversion rates
        used = True
        if r.is_supplementary or r.is_unmapped or r.mate_is_unmapped or not r.is_proper_pair:
            used = False
        # get all mismatches (in lowercase, only mismatches) - ignore INDELs
        mismatches = [m for m in r.get_aligned_pairs(with_seq=True) if m[2] and m[2].islower()]
        for m in mismatches:
            base_qual = r.query_alignment_qualities[m[0]-qstart]
            ref = m[2]
            # if anti-sense, we report A->G as T->C, with strand=-
            # reference is reversed, and base is reported accordingly in get_base_vec
            if strand == '-':
                ref = ref.translate(str.maketrans("ACGTacgtNnXx", "TGCAtgcaNnXx"))
            entry = {
                'contig': contig,
                'start': m[1],
                'end': m[1]+1,
                'name': query_id,
                'score': used,
                'strand': strand,
                'bases11': get_base_vec(seq[m[0]-qstart], strand),
                'ref': str(ref.upper()),
                'base_qual': base_qual,
                'm_pos': m[0]-qstart,
                'qstart': qstart,
                'rlen': rlen,
                'read1': r.is_read1,
                'samflag': flag,
            }
            bed_entries.append(entry)

        # mismatchdetails
        # but we skip unused reads
        if used:
            # For inserts, deletions, skipping either query or reference position may be None
            mismatches = [m for m in r.get_aligned_pairs(with_seq=True) if m[0] and m[2]]
            for m in mismatches:
                mpos = m[0]-qstart
                genomic = m[2].upper()
                read = seq[m[0]-qstart]
                if r.is_read2:
                    mpos += 150
                if (r.is_read1 and strand == '+') or (r.is_read2 and strand == '-'):
                    genomic = genomic.translate(str.maketrans("ACGTacgtNnXx", "TGCAtgcaNnXx"))
                    read = read.translate(str.maketrans("ACGTacgtNnXx", "TGCAtgcaNnXx"))
                key = '{}:{}:{}'.format(genomic, read, mpos)
                mismatch_details_dict[key] += 1

    bed = pd.DataFrame(bed_entries)
#    print("End of get_mismatches()")
    return (mismatch_details_dict, bed)


def get_mismatch_details(details, filename1, filename2):

    # first sort mismatch details list of dict into one dict
    mismatch_details = defaultdict(int)
    for d in details:
        for key, val in d.items():
            mismatch_details[key] += val

    # get coverage
    coverage = defaultdict(int)
    for key, val in mismatch_details.items():
        genomic, read, pos = key.split(':')
        ckey = '{}:{}'.format(genomic, pos)
        if genomic == 'N' or read == 'N':
            continue
        coverage[ckey] += val

    # and convert to final format
    sl = []
    for key, val in mismatch_details.items():
        genomic, read, pos = key.split(':')
        if genomic == read or genomic == 'N' or read == 'N':
            continue
        cov = coverage['{}:{}'.format(genomic, pos)]
        e = {'Category': 'Any',
            'Genomic': genomic,
            'Read': read,
            'Position': int(pos),
            'Coverage': cov,
            'Mismatches': val}
        sl.append(e)
    mismatch_details_df = pd.DataFrame(sl)
    mismatch_details_df = mismatch_details_df[['Category', 'Genomic', 'Read', 'Position', 'Coverage', 'Mismatches']]
    mismatch_details_df.sort_values(by=['Genomic', 'Read', 'Position'], inplace=True)
    mismatch_details_df.to_csv(filename1,
                               sep='\t',
                               index=False,
                               compression='gzip')

    # now sum everything to get overall` mismatches, but split by read
    mismatch_details_df.loc[mismatch_details_df['Position']<150, 'Orientation'] = 'First'
    mismatch_details_df.loc[mismatch_details_df['Position']>149, 'Orientation'] = 'Second'

    cov = mismatch_details_df[['Genomic', 'Orientation', 'Coverage']].drop_duplicates().groupby(['Genomic', 'Orientation']).sum()
    cov = pd.DataFrame(cov.to_records())

    mis = mismatch_details_df[['Genomic', 'Read', 'Orientation', 'Mismatches']].groupby(['Genomic', 'Read', 'Orientation']).sum()
    mis = pd.DataFrame(mis.to_records())

    mismatches_counts = mis.merge(cov, how='left', on=['Genomic', 'Orientation'])
    mismatches_counts['Category'] = 'Any'
    mismatches_counts = mismatches_counts[['Category', 'Orientation', 'Genomic', 'Read', 'Coverage', 'Mismatches']]
    mismatches_counts.sort_values(by=['Orientation', 'Genomic', 'Read'], inplace=True)

    mismatches_counts.to_csv(filename2,
                             sep='\t',
                             index=False,
                             compression='gzip')

    return mismatches_counts


def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="""Split a BAM file containing alignments from
                                     a SLAM/TUC/TL-seq experiment into labelled/unlabelled BAM files.
                                     Requires the MD tag. All reads are kept incl. unmapped,
                                     secondary/supplementary, etc. these can be filtered out when
                                     counting (featureCounts). Optionally, SNPs can be subtracted.""")

    parser.add_argument('bam', help="""The input BAM file (full path).""")

    parser.add_argument('outdir_mm', help="""The output directory (mismatch information).""")

    parser.add_argument('name', help="""The output base name without extension.""")

    parser.add_argument('-s', '--subtract', help="""SNPs to be subtracted (GS default format)""", type=str)

    parser.add_argument('--vcf', help="""Use this flag if SNPs are in VCF format""", action='store_true')

    parser.add_argument('-ref', '--ref-base', help="""Conversion reference base.""",
                        choices=['A', 'C', 'G', 'T'], default='T')

    parser.add_argument('-bc', '--base-change', help="""Conversion base (substitution/mismatch).""",
                        choices=['A', 'C', 'G', 'T'], default='C')

    parser.add_argument('-q', '--base-qual', help="The minimum base quality for any given mismatch (default: 20).",
                        type=int, default=20)

    parser.add_argument('--trim5p', help="The number bases to trim at the 5' ends of reads (default: 0).",
                        type=int, default=0)

    parser.add_argument('--trim3p', help="The number bases to trim at the 3' ends of reads (default: 0).",
                        type=int, default=0)

    parser.add_argument('--overwrite', help='''If this flag is present, then existing files
        will be overwritten.''', action='store_true')

    parser.add_argument('-t', '--tmp', help="""Optional argument: where to write
        temporary files. If not specified, programs-specific tmp will be used.""", default=None)


    utils.add_sbatch_options(parser, num_cpus=default_num_cpus, mem=default_mem)
    utils.add_logging_options(parser)
    args = parser.parse_args()
    utils.update_logging(args)

    msg = "[splbam]: {}".format(' '.join(sys.argv))
    logger.info(msg)

    # if using slurm, submit the script
    if args.use_slurm:
        cmd = "{}".format(' '.join(shlex.quote(s) for s in sys.argv))
        utils.check_sbatch(cmd, args=args)
        return

    # check output path
    exist = utils.check_files_exist([args.outdir_mm],
                                    raise_on_error=True,
                                    logger=logger)

    # check that all files exist
    input_files = [args.bam]
    if args.subtract:
        input_files.append(args.subtract)
    exist = utils.check_files_exist(input_files,
                                    raise_on_error=True,
                                    logger=logger)

    print("bam input: {}".format(args.bam))
    print("outdir_mm input: {}".format(args.outdir_mm))
    print("name input: {}".format(args.name))

    # create the output files - mismatches
    mismatch_details_filename = '{}.mismatchDetails.tab.gz'.format(args.name)
    mismatch_details_filename = os.path.join(args.outdir_mm, mismatch_details_filename)

    mismatch_filename = '{}.mismatches.tab.gz'.format(args.name)
    mismatch_filename = os.path.join(args.outdir_mm, mismatch_filename)

    mismatch_filename_final = '{}.mismatches-used.tab.gz'.format(args.name)
    mismatch_filename_final = os.path.join(args.outdir_mm, mismatch_filename_final)

    mismatch_all_filename = '{}.mismatch_all_final_list.tab.gz'.format(args.name)
    mismatch_all_filename = os.path.join(args.outdir_mm, mismatch_all_filename)

    trueconversions_filename = '{}.trueconversions_list.tab'.format(args.name)
    trueconversions_filename = os.path.join(args.outdir_mm, trueconversions_filename)

    print("mismatche details output: {}".format(mismatch_details_filename))
    print("mismatches output: {}".format(mismatch_filename))
    print("all mismatches output: {}".format(mismatch_filename_final))
    print("all mismatches output: {}".format(mismatch_all_filename))
    print("True conversions output: {}".format(trueconversions_filename))

    # and check if exist
    out_files = [mismatch_details_filename, mismatch_filename, mismatch_filename_final, mismatch_all_filename, trueconversions_filename]

    all_out_exists = all([os.path.exists(of) for of in out_files])
    if args.overwrite or not all_out_exists:
        pass
    else:
        msg = "All output files {} already exist. Skipping call.".format(out_files)
        logger.warning(msg)
        return

    msg = "Getting all mismatches"
    logger.info(msg)

    # first get all mismatches
    print("Start collecting mismatches")
    bam = ps.AlignmentFile(args.bam, "rb")
    SN = (SQ['SN'] for SQ in bam.header['SQ'])
    bam.close()

    mismatch_and_details = utils.apply_parallel_iter(SN,
                                                     args.num_cpus,
                                                     get_mismatches,
                                                     args.bam,
                                                     progress_bar=False,
                                                     backend='multiprocessing')

    all_details = [a for a, b in mismatch_and_details]
    mismatch_count = get_mismatch_details(all_details,
                                          mismatch_details_filename,
                                          mismatch_filename)

    all_mismatches = [b for a, b in mismatch_and_details]
    all_mismatches = pd.concat(all_mismatches)

#    print("save all_mismatches_before_filtering.txt")
#    with open('all_mismatches_before_filtering.txt', 'w') as f:
#        f.write(tabulate(all_mismatches))

    # get the conversion of interest
    print("get the conversion of interest")
    bc = get_base_vec(args.base_change, '+')
    m_bc = all_mismatches['bases11'] == bc
    m_ref = all_mismatches['ref'] == args.ref_base
    all_mismatches = all_mismatches[m_ref & m_bc]

    # filter base quality - no offset of 33 needs to be subtracted
    print("filter base quality - no offset of 33 needs to be subtracted")
    m_qual = all_mismatches['base_qual']>=args.base_qual

    # below we keep track of discarded mismatches to adjust rates
    print("below we keep track of discarded mismatches to adjust rates")
    discarded = all_mismatches[~m_qual].copy()
    m = (discarded['read1']==True) & (discarded['score']==True)
    discarded_first = discarded[m].shape[0]
    m = (discarded['read1']==False) & (discarded['score']==True)
    discarded_second = discarded[m].shape[0]

    all_mismatches = all_mismatches[m_qual]

    # discard mismatches found at read ends
    print("discard mismatches found at read ends")
    m_trim5p = all_mismatches['m_pos'] < (args.trim5p - all_mismatches['qstart'])
    m_trim3p = all_mismatches['m_pos'] >= (all_mismatches['rlen'] - args.trim3p)
    discarded = all_mismatches[m_trim5p | m_trim3p].copy()
    m = (discarded['read1']==True) & (discarded['score']==True)
    discarded_first += discarded[m].shape[0]
    m = (discarded['read1']==False) & (discarded['score']==True)
    discarded_second += discarded[m].shape[0]

    all_mismatches = all_mismatches[~m_trim5p & ~m_trim3p]

    # remove all SNPs from what remains
    print("remove all SNPs from what remains")
    if args.subtract:
        # currently GRAND-SLAM snpdata default format
        if args.vcf:
            snps = utils.fmt_convert(args.subtract)
        else:
            snps = pd.read_csv(args.subtract, sep='\t')
        snps = snps.Location.unique()
        # add field
        all_mismatches['Location'] = all_mismatches[['contig', 'start']].apply(lambda x: ':'.join([str(s) for s in x]), axis=1)

        discarded = all_mismatches[all_mismatches.Location.isin(snps)].copy()
        m = (discarded['read1']==True) & (discarded['score']==True)
        discarded_first += discarded[m].shape[0]
        m = (discarded['read1']==False) & (discarded['score']==True)
        discarded_second += discarded[m].shape[0]

        all_mismatches = all_mismatches[~all_mismatches.Location.isin(snps)]

    # adjust final mismatch counts
    print("adjust final mismatch counts")
    m = (mismatch_count.Orientation=='First') & (mismatch_count.Genomic=='A') & (mismatch_count.Read=='G')
    mismatch_count.loc[m, 'Mismatches'] = mismatch_count.loc[m, 'Mismatches'] - discarded_first
    n = (mismatch_count.Orientation=='Second') & (mismatch_count.Genomic=='T') & (mismatch_count.Read=='C')
    mismatch_count.loc[n, 'Mismatches'] = mismatch_count.loc[n, 'Mismatches'] - discarded_second
    mismatch_count = mismatch_count[m | n]
    mismatch_count.to_csv(mismatch_filename_final,
                          sep='\t',
                          index=False,
                          compression='gzip')

    print("output all final mismatches")
    all_mismatches.to_csv(mismatch_all_filename,
                             sep='\t',
                             index=False,
                             compression='gzip')

    print("output all true conversions")
    true_conversions = all_mismatches.name.unique()
    with open(trueconversions_filename, 'w') as f:
        for item in true_conversions:
            f.write("%s\n" % item)


    print("Python script finished")



if __name__ == '__main__':
    main()
