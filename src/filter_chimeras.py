#!/mnt/software/unstowable/anaconda/bin/python
"""FIXME:add-doc
"""

import sys
import os
import argparse

from chimeras2 import parse_chimeras, chim_gets_filtered


def main():
    """main function
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--chim-in",
                        required=True,
                        help="Input chimera file")
    default = "-"
    parser.add_argument("-o", "--chim-out",
                        default=default,
                        help="Output chimera file")
    default = 1
    parser.add_argument("--min-mq",
                        default=default,
                        type=int,
                        help="Minimum mapping quality (0==ambigiously mapped; default = {})".format(default))
    parser.add_argument("--no-revcomp",
                        action="store_true",
                        help="Do not apply revcomp filtering")
    parser.add_argument("-t", "--type", choices=["intra", "inter"],
                        help="Type only")

    args = parser.parse_args()

    if not os.path.exists(args.chim_in) and  args.chim_in != "-":
        parser.error("Non-existant file {}".format(args.chim_in))
        sys.exit(1)

    if args.chim_out != '-':
        assert not os.path.exists(args.chim_out)
        fh_out = open(args.chim_out, 'w')
    else:
        fh_out = sys.stdout

    if args.chim_in == "-":
        fh = sys.stdin
    else:
        fh = open(args.chim_in)
    for c in parse_chimeras(fh):
        if not chim_gets_filtered(c, min_mq=args.min_mq, remove_compl=(not args.no_revcomp)):
            if args.type and c.type() != args.type:
                continue
            fh_out.write("{}\n".format(c))
    if args.chim_in != "-":
        fh.close()

    if fh_out != sys.stdout:
        fh_out.close()

if __name__ == "__main__":
    main()
