#!/mnt/software/unstowable/anaconda/bin/python2
##!/usr/bin/env python2
"""Chimera finding based on BWA-MEM split/secondary alignment
"""

import sys
import os
#import gzip
import argparse
import regex

import pysam
import chimeras2
# https://pypi.python.org/pypi/cigar/0.1
from cigar import Cigar


# anything below treated as singletons
MIN_CHIM_DIST = 50


def query_aln_seq(queryseq, cigar):
    """Infer aligned bit of query sequence from cigar

    >>> query_aln_seq("ACGTACGT", "8M")
    'ACGTACGT'
    >>> query_aln_seq("ACGTACGT", "1M7S")
    'A'
    >>> query_aln_seq("ACGTACGT", "6S2M")
    'GT'
    >>> query_aln_seq("ACGTACGT", "6S2M")
    'GT'
    >>> query_aln_seq("ACGTACGT", "2S4M2S")
    'GTAC'
    >>> query_aln_seq("ACGTACGT", "2S1M2I1M2S")
    'GTAC'
    """

    cig = chimeras2.cigar_to_list(cigar)

    start = 0
    for (op_len, op_type) in cig:
        if op_type == 'S':
            start += op_len
        else:
            break

    end = len(queryseq)
    for (op_len, op_type) in cig[::-1]:
        if op_type == 'S':
            end -= op_len
        else:
            break

    assert end > 0 and start < end
    assert cigar2querylen(cig) == end-start, ("{} != {}".format(cigar2querylen(cig), end-start))

    return queryseq[start:end]


def cigar2querylen(cigar):
    """Determine length of query sequence from cigar string
    """

    cig = chimeras2.cigar_to_list(cigar)
    qlen = 0
    for (op_len, op_type) in cig:
        assert not isinstance(op_type, int), (
            "Do not understand cigar {}. Might be from pysam?".format(op_type))
        if op_type in ['M', 'I', '=', 'X']:
            # FIXME check SAM spec for more
            qlen += op_len
    return qlen


def cigar2reflen(cigar):
    """infer ref length from cigar struct or string

    >>> cigar2reflen("10M")
    10
    >>> cigar2reflen("10M1S")
    10
    >>> cigar2reflen("1S10M")
    10
    >>> cigar2reflen("1S10M1S")
    10
    >>> cigar2reflen("1S5M1D5M1S")
    11
    >>> cigar2reflen("1S5M1I5M1S")
    10
    """

    cig = chimeras2.cigar_to_list(cigar)

    rlen = 0
    #ref_consuming = 'MDN=X'
    for (op_len, op_type) in cig:
        if op_type in Cigar.ref_consuming_ops:
            rlen += op_len
    return rlen


def read_is_primary(read):
    """see sam spec: For each read/contig in a SAM file, it is required
    that one and only one line associated with the read satisfies
    'FLAG & 0x900 == 0'. This line is called the primary line of the read.
    """
    return bool(read.flag & 0x900 == 0)


def find_sa_chimeras(bam, fh_out,
                     sa_singletons_only=False,# anything with dist below min_chim_dist but >0
                     min_chim_dist=MIN_CHIM_DIST,# anythin below min_chim_dist is called sa_singleton (split alignment singleton)
                     max_nm_perc=100,
                     ign_string=None,
                     ign_max_err=0,
                     add_nm=False):
    """Find chimeras from BWA-MEM split mappings (defined as chimeras here)
    """

    debug = False
    if debug:
        num_chimeras = 0

    fh_in = pysam.Samfile(bam)
    for (num_reads, read) in enumerate(fh_in):

        if debug and (num_chimeras > 100 or num_reads) > 10000:
            sys.stderr.write("DEBUG break\n")
            break
        if read.is_unmapped or read.is_qcfail or read.is_duplicate:
            continue
        if not read_is_primary(read):
            continue

        qname = read.qname
        qseq = read.seq
        # pysam 0.7.7 uses seq and not query_sequence (clipped anyway?)
        # pysam 0.7.7 uses qname not query_name

        if ign_string:
            if regex.search("({}){{e<={}}}".format(ign_string, ign_max_err), qseq, regex.BESTMATCH):
                #sys.stderr.write("DEBUG: ignoring {}\n".format(qseq))
                continue

        # previously ignored MQ0. Should be decided downstream.
        # Still valuable in large genomes.
        #
        # if not read.mapping_quality > 0:
        # pysam 0.7.7 has mapq and not mapping_quality
        #if not read.mapq > 0:
        #    continue

        # finding split alignments via SA tag in primary alignment
        tags = dict(read.tags)
        if not tags.has_key('SA'):
            continue

        assert 'NM' in tags
        perc_nm = tags['NM']*100.0/float(len(
            query_aln_seq(read.seq, read.cigarstring)))
        if perc_nm > max_nm_perc:
            continue

        ori_cigar = list(Cigar(read.cigarstring).items())

        # no need to determine clip site here and later. overlap/proximity of
        # mapping determines chimera already
        #
        # skip if clip site can't be determined
        #if not clip_site(ori_cigar):
        #    sys.stderr.write(
        #        "WARN: can't determine clip site (or clip too small) for {} in {}\n".format(
        #            ori_cigar, qname))
        #    continue

        # SA == supplementary alignment
        # for definition of SA tag see:
        # https://sourceforge.net/p/samtools/mailman/message/30853577/
        # chr,strandPos,CIGAR,mapQ,NM;
        num_valid_sa = 0
        for (_sa_num, sa_tag) in enumerate(tags['SA'].rstrip(";").split(';')):
            sa_tag = dict(zip(['chrom', 'pos', 'strand', 'cigar', 'mq', 'nm'],
                              sa_tag.split(",")))
            assert sa_tag['strand'] in ['+', '-']
            for k in ['pos', 'mq', 'nm']:
                sa_tag[k] = int(sa_tag[k])

            # previously ignored MQ0. Should be decided downstream.
            # Still valuable in large genomes.
            #if int(sa_tag['mq']) == 0:
            #    continue
            
            sa_cigar = list(Cigar(sa_tag['cigar']).items())
            #if not clip_site(sa_cigar):
            #    sys.stderr.write(
            #        "WARN: can't determine clip site (or clip too small) for {} in {} (SA)\n".format(
            #            sa_cigar, qname))
            #    continue

            # clips have to be on opposite sites
            #if clip_site(ori_cigar) == clip_site(sa_cigar):
            #    sys.stderr.write(
            #        "WARN: clip on identical sites for {}: {} and {}\n".format(
            #            qname, ori_cigar, sa_cigar))
            #    continue

            perc_nm = sa_tag['nm']*100.0/float(len(
                query_aln_seq(read.seq, sa_tag['cigar'])))
            if perc_nm > max_nm_perc:
                continue

            # indirect testing of cigar2rlen()
            # pysam 0.7.7 has aend and not reference_end
            #assert read.pos + cigar2rlen(ori_cigar) == read.reference_end
            assert read.pos + cigar2reflen(ori_cigar) == read.aend

            chim = chimeras2.Chimera()
            chim.qname = qname

            # default is to assign primary to left and SA to right and then sort

            chim.left.flag = read.flag
            chim.left.rname = fh_in.getrname(read.tid)
            # pysam 0.7.7 uses tid not reference_id
            chim.left.pos = read.pos
            chim.left.aend = read.aend
            chim.left.mapq = read.mapq
            chim.left.cigar = read.cigarstring
            chim.left.seq = query_aln_seq(qseq, chim.left.cigar)
            if add_nm:
                chim.left.nm = tags['NM']
            
            chim.right.flag = 0 if sa_tag['strand'] == '+' else 16
            chim.right.rname = sa_tag['chrom']
            chim.right.pos = sa_tag['pos']
            chim.right.aend = sa_tag['pos'] + cigar2reflen(sa_cigar)
            chim.right.mapq = sa_tag['mq']
            chim.right.cigar = sa_tag['cigar']
            chim.right.seq = query_aln_seq(qseq, chim.right.cigar)
            if add_nm:
                chim.right.nm = sa_tag['nm']

            chim.sort_halves()
            if chim.order == "invalid":
                # likely revcomp vs non-revcomp
                continue

            # no overlap allowed ever
            if chim.dist2d() < 0:
                continue
            # sa_singleton and wanted?
            if chim.dist2d() < min_chim_dist:
                if not sa_singletons_only:
                    continue
            else:
                if sa_singletons_only:
                    continue

            #sys.stderr.write(
            #    "DEBUG chim before sanity_check = {} (dist2d={})\n".format(
            #        chim, chim.dist2d()))
            chim.sanity_check()

            num_valid_sa += 1
            if num_valid_sa > 1:
                sys.stderr.write(
                    "WARN: More than one valid SA found for {}\n".format(qname))

            fh_out.write("{}\n".format(chim))
            if debug:
                num_chimeras += 1

    fh_in.close()


def main():
    """main
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--bam",
                        required=True,
                        help="Input BAM file")
    parser.add_argument("-o", "--chim",
                        default="-",
                        help="Output chimera file (- for stdout; default)")
    parser.add_argument("--sa-singletons",
                        action="store_true",
                        help="Predict split alignment (fake) singletons only (invert logic)")
    parser.add_argument("--add-nm",
                        action="store_true",
                        help="Add number of mismatches to ouput (likely breaks downstream scripts)")
    parser.add_argument("--min-chim-dist",
                        default=MIN_CHIM_DIST, type=int,
                        help="Minimum required distance of chimeric ends (default {})".format(MIN_CHIM_DIST))
    DEF_MAX_NM_PERC = 100
    parser.add_argument("--max-nm-perc",
                        default=DEF_MAX_NM_PERC, type=int,
                        help="Maximum percent of mismatches allowed in either end (default {})".format(DEF_MAX_NM_PERC))
    parser.add_argument("--ign_string",
                        help="sequence string to ignore e.g. linker")
    default = 0
    parser.add_argument("--ign_max_err", default=default, type=int,
                        help="maximum number of mismatches allowed in ignore string (see above; default={})".format(default))
    args = parser.parse_args()


    if not os.path.exists(args.bam):
        parser.error("Non-existant file {}".format(args.bam))
        sys.exit(1)
    if os.path.exists(args.chim) and args.chim != "-":
        parser.error("Refusing to overwrite file {}".format(args.chim))
        sys.exit(1)

    if args.chim == "-":
        chimfh = sys.stdout
    else:
        chimfh = open(args.chim, 'w')

    find_sa_chimeras(args.bam, chimfh,
                     sa_singletons_only=args.sa_singletons,
                     min_chim_dist=args.min_chim_dist,
                     max_nm_perc=args.max_nm_perc,
                     ign_string=args.ign_string,
                     ign_max_err=args.ign_max_err,
                     add_nm=args.add_nm)

    if chimfh != sys.stdout:
        chimfh.close()


if __name__ == "__main__":
    main()
