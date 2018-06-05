#!/usr/bin/env python
"""
Divide the multiple alignment file in mln format into separate transcript files
"""

__date_ = "2014-05-23"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import os
import sys
import time
import itertools
import cPickle as cp
from argparse import ArgumentParser

parser = ArgumentParser(description=__doc__)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--mln",
                    dest="mln",
                    help="Alignment file in mln format")
parser.add_argument("--output-dir",
                    dest="output_dir",
                    help="Directory for the divided files")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def isa_group_separator(line):
    return line=='\n'

def read_mln_to_dict(filetoread):
    data = {}
    with open(filetoread) as f:
        for key, block in itertools.groupby(f, isa_group_separator):
            group = list(block)
            if not key:
                tmp = {}
                for i in range(0, len(group)-2, 3):
                    tmp[group[i].rstrip().split(" ")[-1]] = [group[i+1].rstrip(), group[i+2].rstrip()]
                data[group[i].rstrip().split(" ")[0][1:].replace('|', '__')] = tmp
    return data

def main(options):
    """Main logic of the script"""
    try:
        if options.verbose:
            syserr("Reading alignment file.\n")
        multiple_alignment_dict = read_mln_to_dict(options.mln)
    except IOError:
        raise IOError("Cannot open multiple alignment file %s" % options.mln)

    for name, value in multiple_alignment_dict.iteritems():
        if options.verbose:
            syserr("Saving %s to %s\r" % (name, options.output_dir))
        with open(os.path.join(options.output_dir, name), 'wb') as aln:
            cp.dump(value, aln)
    if options.verbose:
        syserr("\n")

if __name__ == '__main__':
    try:
        try:
            options = parser.parse_args()
        except Exception, e:
            parser.print_help()
            sys.exit()
        if options.verbose:
            start_time = time.time()
            start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
            syserr("############## Started script on %s ##############\n" % start_date)
        main(options)
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" % (time.time() - start_time, time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" % (time.time() - start_time))
        sys.exit(-1)
