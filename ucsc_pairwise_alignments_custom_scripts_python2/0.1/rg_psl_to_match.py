#!/usr/bin/env python
"""
Convert adjusted PSL file to the match.tab file defined by Jean Hausser
"""

__date__ = "2016-06-21"
__author__ = "Rafal Gumienny"
__email__ = "r.gumienny@unibas.ch"
__license__ = "GPL"

# imports
import sys
import time
import numpy as np
import pandas as pd
from contextlib import contextmanager
from argparse import ArgumentParser, RawTextHelpFormatter

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
parser.add_argument("--input",
                    dest="input",
                    default=sys.stdin,
                    help="Input file in adjusted PSL format. Defaults to sys.stdin.")
parser.add_argument("--output",
                    dest="output",
                    default=sys.stdout,
                    help="Output file in match.tab format. Defaults to sys.stdout.")


# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Main logic of the script"""
    all_tr = 0
    written_tr = 0
    with smart_open(options.input, 'r') as psl, smart_open(options.output, 'w') as out:
        for line in psl:
            all_tr += 1
            try:
                l = line.rstrip().split('\t')
                name = l[0]
                chr = l[1]
                strand = l[2]
                nblocks = int(l[3])
                domain = ''
            except KeyError:
                print name
                continue
            i = 1
            written_tr += 1
            for plen, trcoor, gcoor in zip(l[4].split(','), l[5].split(','), l[6].split(','))[:-1]:
                outtext = "%s%s.%i\t%s\t%i\t%i\t%s\t%i\t%i\n"%(name,
                                                               domain,
                                                               i,
                                                               chr,
                                                               int(gcoor)+1,
                                                               int(gcoor)+int(plen),
                                                               strand,
                                                               int(trcoor)+1,
                                                               int(trcoor)+int(plen))
                out.write(outtext)
                i+=1
    out.close()
    if options.verbose:
        syserr("Wrote %i sequences out of %i\n" % (written_tr, all_tr))


# this function is also defined in utils but I put it here to avoid
# unnecessary module import that might not be available everywhere as
# it is my own module
@contextmanager
def smart_open(filepath, mode='r'):
    """Open file intelligently depending on the source

    :param filepath: can be both path to file or sys.stdin or sys.stdout
    :param mode: mode can be read "r" or write "w". Defaults to "r"
    :yield: context manager for file handle

    """
    if mode == 'r':
        if filepath is not sys.stdin:
            fh = open(filepath, 'r')
        else:
            fh = filepath
        try:
            yield fh
        except IOError as e:
            if fh is not sys.stdin:
                fh.close()
            elif e.errno == errno.EPIPE:
                pass
        finally:
            if fh is not sys.stdin:
                fh.close()
    elif mode == 'w':
        if filepath is not sys.stdout:
            fh = open(filepath, 'w')
        else:
            fh = filepath
        try:
            yield fh
        finally:
            if fh is not sys.stdout:
                fh.close()
    else:
        raise NoSuchModeException("No mode %s for file" % mode)


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
            syserr("############## Started script on %s ##############\n" %
                   start_date)
        main(options)
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" %
                   (time.time() - start_time,
                    time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" %
               (time.time() - start_time))
        sys.exit(-1)
