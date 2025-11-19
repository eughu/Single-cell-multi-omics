#-*- coding:utf-8 -*-
from __future__ import division
import re
import sys
import wWigIO
import functools
import numpy as np
import bw_bin
from optparse   import OptionParser

def prepare_optparser():
    usage ="""usage: %s [options] sample_input_information_xls

Get chromosomes lists from fai file and then sort it in ascii order.

Then merge it into one file, unzipped.

Example:
    python %s --ext_len 1000 --out_bin 100 /date/huboqiang/NORM_seq_PGC/human/02.SingleC/hPGC_Week13_rep2_1/singleC/all.GCA.GCC.GCT.NDR.bed /date/huboqiang/NORM_seq_PGC/human/02.SingleC/hPGC_Week13_rep2_1/singleC/all.GCA.GCC.GCT.bw

    """ % (sys.argv[0],sys.argv[0])

    description = "Merge singleC results."

    optparser = OptionParser(
        version="%s v0.1 2015.8" % (sys.argv[0]),
        description=description,
        usage=usage,
        add_help_option=False
    )

    optparser.add_option(
        "-e", "--ext_len", default=1000,
        help="\nLength of extension from the middle of NDR. [default: %default]"
    )
    
    optparser.add_option(
        "-b", "--out_bin", default=100,
        help="\nBin number for plot. [default: %default]"
    )
    return optparser

def main():
    prepare_optparser()
    (options,args) = prepare_optparser().parse_args()
    try:
        in_bed = args[0]
        in_bw = args[1]
        ext_len = int(options.ext_len)
        out_bin = int(options.out_bin)

    except IndexError:
        prepare_optparser().print_help()
        sys.exit(1)
    
    wWigIO.open(in_bw)
    with open(in_bed, "r") as f_inbed:
        for line in f_inbed:
            line = line.strip()
            f = line.split()
            chrom = f[0]
            begin = int(f[1])
            endin = int(f[2])
            beg_pos = int((begin+endin)/2) - ext_len
            end_pos = int((begin+endin)/2) + ext_len
            np_out = np.zeros([2, out_bin])
            np_bw = np.array(wWigIO.getIntervals(in_bw, chrom, beg_pos, end_pos))
            map(functools.partial(bw_bin.parse_array, beg_pos, end_pos, out_bin, np_out), np_bw)
            l_out = "\t".join( np.array(np_out[1,:]/np_out[0,:], dtype="string") )
            print "%s:%d-%d\t%s" % (chrom, beg_pos, end_pos, l_out)
    wWigIO.close()
    
if __name__ == '__main__':
    main()