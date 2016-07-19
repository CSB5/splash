"""
All positions are one based half open (like bed coordinates or python slices)
"""

import os
import sys

# https://pypi.python.org/pypi/cigar/0.1
from cigar import Cigar


__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2016 Genome Institute of Singapore"
__license__ = "The MIT License (MIT)"


def cigar_to_list(cigar):
    """Convencience function that converts list, string or cigarstruct cigar to
    list of cigarstruct items
    """
    if isinstance(cigar, list):
        cig = cigar
    elif not isinstance(cigar, Cigar):
        cig = list(Cigar(cigar).items())
    else:
        cig = list(cigar)
    return cig


def cigar_first_match_pos(cigar):
    """Return first aligned read position

    >>> cigar_first_match_pos("55S9M1D33M")
    55
    >>> cigar_first_match_pos("67S21M9S")
    67
    >>> cigar_first_match_pos("30M62S")
    0
    >>> cigar_first_match_pos("10H30M")
    0
    """
    cig = cigar_to_list(cigar)
    pos = 0
    for (op_len, op_type) in cig:
        if op_type in 'M=X':
            return pos
        elif op_type == 'H':
            continue
        else:
            pos += op_len


def clip_site(cigar, minclip=10):
    """which contains larger clip.

    previously did not allow for end clippings in matching site, e.g. 50M1S.


    >>> c = list(Cigar("55S9M1D33M").items())
    >>> clip_site(c)
    'L'
    >>> c = list(Cigar("67S21M9S").items())
    >>> clip_site(c)
    'L'
    >>> c = list(Cigar("30M62S").items())
    >>> clip_site(c)
    'R'
    >>> clip_site("32S21M18S")
    'L'
    """

    cig = cigar_to_list(cigar)

    lclip = rclip = 0
    if cig[0][1] in 'HS':
        lclip = cig[0][0]
    if cig[-1][1] in 'HS':
        rclip = cig[-1][0]

    if lclip > rclip and lclip >= minclip:
        return 'L'
    elif lclip < rclip and rclip >= minclip:
        return 'R'
    else:
        return False


class ChimericHalf(object):
    r"""Describes one chimera half

    http://stackoverflow.com/questions/8834916/how-can-i-include-special-characters-tab-newline-in-a-python-doctest-result-s

    >>> chimhalf = ChimericHalf()
    >>> s = ['0', 'abc', '9', '100', '60', '90M', 'None']
    >>> chimhalf.from_split_str(s)
    >>> str(chimhalf)
    '0\tabc\t9\t100\t60\t90M\tNone'
    >>> chimhalf.is_reverse()
    False
    """


    def __init__(self):
        """init"""
        self.flag = None# sample flag
        self.rname = None# ref name
        self.pos = None# 0 based
        self.aend = None# half open
        self.mapq = None
        self.cigar = None
        self.seq = None
        self.nm = None# temp hack for shenyang to get nm reported

        
    def __str__(self):
        """string representation
        """

        # temp hack for shenyang to get nm reported
        if self.nm is None:
            return "\t".join([str(x) for x in self.flag, self.rname, self.pos,
                              self.aend, self.mapq, self.cigar, self.seq])
        else:
            return "\t".join([str(x) for x in self.flag, self.rname, self.pos,
                              self.aend, self.mapq, self.cigar, self.seq, self.nm])
            

    def from_split_str(self, splits):
        """parse split string representation of chimeric half"""
        assert len(splits) == 7, (splits)
        i = 0
        self.flag = int(splits[i]) if splits[i] != 'None' else None
        i += 1
        self.rname = splits[i]
        i += 1
        self.pos = int(splits[i])
        i += 1
        self.aend = int(splits[i])
        i += 1
        self.mapq = int(splits[i]) if splits[i] != 'None' else None
        i += 1
        self.cigar = splits[i] if splits[i] != 'None' else None
        i += 1
        self.seq = splits[i] if splits[i] != 'None' else None


    def is_reverse(self):
        """return whether seq is mapped as reverse complement
        """
        return bool(self.flag & 0x10)


class Chimera(object):
    r"""Chimera object



    http://stackoverflow.com/questions/13106118/object-reuse-in-python-doctest

    >>> s = ['myqname', '0', 'abc', '0', '10', '60', '10M', 'AAAAAAAAAAA', '0', 'abc', '110', '115', '60', '5M', 'CCCCC', 'LR']
    >>> chim1 = Chimera('\t'.join(s))
    >>> chim1.dist2d()
    100
    >>> chim1.sanity_check()
    """


    def __init__(self, line=None, strict=True):
        """init"""
        self.qname = None# read name
        self.left = ChimericHalf()
        self.right = ChimericHalf()
        self.order = None
        self.extra = None# for eval chimera compat if strict not True. adding any after default elements in here
        self.strict = strict
        self.debug_stream = sys.stderr
        self.debug_on = False
        if line:
            self.from_line(line)


    def debug(self, msg):
        """debug printing
        """
        if self.debug_on:
            self.debug_stream.write("DEBUG: {}\n".format(msg))


    def __str__(self):
        """to string"""
        return "{}\t{}\t{}\t{}{}".format(
            self.qname, self.left, self.right, self.order, "\t{}".format(self.extra) if self.extra else "")


    def from_strings(self, strings):
        """parse chimera from values, e.g. already split line

        consumes strings
        """

        #import sys; sys.stderr.write("DEBUG strings={}\n".format(strings))
        self.qname = strings.pop(0)
        self.left.from_split_str([strings.pop(0) for _ in range(7)])
        self.right.from_split_str([strings.pop(0) for _ in range(7)])
        self.order = strings.pop(0)
        if len(strings) and not self.strict:
            self.extra = '\t'.join(strings)
            strings = []
        self.sanity_check()


    def from_line(self, line):
        """parse chimera from line
        """
        self.from_strings(line.rstrip().split("\t"))


    def sanity_check(self):
        """sanity check of values"""
        if self.order is not None:
            valid_order_values = ["LR", "RL", "invalid"]
            #assert self.order in ["LR", "RL"], ("meh {} {}".format(self.order, type(self.order)))
            if self.order not in valid_order_values:
                raise ValueError("Unknown order '{}' must be one of {}".format(self.order, valid_order_values))

        if self.left.rname == self.right.rname:
            assert self.left.is_reverse() == self.right.is_reverse()

        # assert sorting (reflected in order)
        if self.left.rname == self.right.rname:
            assert self.left.pos < self.right.pos
            assert self.left.aend <= self.right.pos
        else:
            assert self.left.rname < self.right.rname


    def type(self):
        """type of chimera: intra or inter
        """
        assert self.left.rname and self.right.rname
        if self.left.rname == self.right.rname:
            return "intra"
        else:
            return "inter"


    def sort_halves(self):
        """sort chimeric halves and adds order field which describes
        the read half orientation

        sorting allows to use the chimeric halves as identifiers even
        if their order was swapped due to ligation on different sites.
        since the latter is interesting for testing we still keep it.
        """

        assert self.order == None

        # first fix arbitray order defined by mapping/BAM so that
        # left/right are in order of read/cigar
        #
        left_first_match_pos = cigar_first_match_pos(self.left.cigar)
        right_first_match_pos = cigar_first_match_pos(self.right.cigar)
        if left_first_match_pos > right_first_match_pos:
            self.left, self.right = self.right, self.left

        # now apply genomic order and record in order field
        #
        # inter: arbitrary
        if self.left.rname > self.right.rname:
            self.left, self.right = self.right, self.left
            self.order = "RL"
        elif self.left.rname < self.right.rname:
            self.order = "LR"

        else:# intra LR/RL tag actually ignors whether both were
            # revcomp'ed or not. what matters is that both have same
            # orientation. if they don't, then chimera is invalid, yet
            # still okay for inter (since ref dependent)
            if not (self.left.is_reverse() == self.right.is_reverse()):
                self.order = "invalid"
            elif self.left.pos > self.right.pos:
                self.left, self.right = self.right, self.left
                self.order = "RL"
            else:
                self.order = "LR"


    def dist2d(self):
        """compute 2d i.e. sequence distance"""
        if self.left.rname != self.right.rname:
            return sys.maxint
        else:
            # assume it's sorted i.e. left is 5'
            assert self.order
            #if self.right.pos - self.left.aend < 0:
            #    sys.stderr.write("WARN: negative distance for {}\n".format(self))
            # NOTE if they slightly overlap we might get negative distances
            return self.right.pos - self.left.aend


def parse_chimeras(fh, strict=True):
    """Parse all chimeras from stream and yields them separately

    NOTE: load_chimeras() is a more sophisticated interface
    """

    for line in fh:
        if line.startswith('#'):
            continue
        yield Chimera(line, strict=strict)



def chim_is_rrna_rrna(chim):
    return all([x.rname.startswith("human-") for x in [chim.left, chim.right]])


def chim_gets_filtered(chim, min_mq=1, remove_compl=True):

    # mq filter
    if chim.left.mapq < min_mq or chim.right.mapq < min_mq:
        return True

    # remove if complementary
    if remove_compl:
        if any([half.flag & 0x10 for half in [chim.left, chim.right]]):
            return True

    return False


def load_chimeras(chimfile, min_mq=1,
                  remove_rrna_rrna=False, only_rrna_rrna=False,
                  remove_compl=True):
    """Reads chimeras listed in chimfile based on given filtering
    criteria.  Returned is dictionary with keys 'inter' and 'intra'
    each containing a list of chimeras of corresponding type

    """

    if remove_rrna_rrna:
        assert not only_rrna_rrna

    chimeras = {'intra': [],
                'inter': []}
    with open(chimfile) as fh:
        for line in fh:
            if line.startswith('#'):
                continue

            #chim = EvalrRNAChimera(line=line)
            chim = Chimera(line=line)

             # remove rrna rrna and assign to type.
             # works for inter and intra
            if remove_rrna_rrna:
                if chim_is_rrna_rrna(chim):
                    continue
            elif only_rrna_rrna:
                if not chim_is_rrna_rrna(chim):
                    continue

            if chim_gets_filtered(chim, min_mq, remove_compl):
                continue

            chimeras[chim.type()].append(chim)
    return chimeras


def write_chimeras(fname, chims):
    if fname == "-":
        fh = sys.stdout
    else:
        assert not os.path.exists(fname)
        fh = open(fname, 'w')
    for c in chims:
        fh.write("{}\n".format(c))
    if fh != sys.stdout:
        fh.close()


def get_midpoint_bin_key(chim, win=100):
    """return a bin key for a chimera

    Returns four tuple describing midpoint
    (actual region is returned region+100)
    """

    midpoint = dict()
    for name, half in (["L", chim.left],
                       ["R", chim.right]):
        reglen = half.aend - half.pos
        assert reglen > 0
        midpoint[name] = half.pos + reglen//2

    bin_reg = dict()
    for name in midpoint.keys():
        s = midpoint[name]//win*win
        bin_reg[name] = (s, s+win)

    return (chim.left.rname, bin_reg["L"],
            chim.right.rname, bin_reg["R"])


def midpoint_bin_key_to_str(key):
    """FIXME:add-doc"""
    return "{}\t{}\t{}\t{}\t{}\t{}".format(key[0], key[1][0], key[1][1], key[2], key[3][0], key[3][1])


def midpoint_binning(chims):
    """FIXME:add-doc"""
    bins = dict()
    for c in chims:
        key = get_midpoint_bin_key(c)
        if not key in bins:
            bins[key] = []
        bins[key].append(c)
    return bins
