#!/usr/bin/env python

__version__ = "0.1"
__author__ = "Foivos Gypas"
__contact__ = "foivos.gypas@unibas.ch"
__doc__ = "Split pairwise alignment in files based on the chromosome"

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# import needed (external) modules
# -----------------------------------------------------------------------------

import sys
import os
from argparse import ArgumentParser, RawTextHelpFormatter

# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Main function
# -----------------------------------------------------------------------------

def readAlignment(fp):

    """
    Read alignment generator:
    Consider lines staring with an integer and returns them together with the
    two following lines that represent the pairwise alignment.
    """

    name, alignment1, alignment2 = None, None, None

    for line in fp:
        if line[0].isdigit() and line is not os.linesep:
            name = line
        elif name is not None and alignment1 is None and line is not os.linesep:
            alignment1 = line
        elif name is not None and alignment1 is not None and alignment2 is None and line is not os.linesep:
            alignment2 = line
            yield (name, alignment1, alignment2)
        else:
            name, alignment1, alignment2 = None, None, None

def find_chromosome(l):

    """ Get name line and return chromosome """

    return l.strip().split(" ")[1]


def main():
    """ Main function """

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter
    )

    parser.add_argument(
        "--alignment",
        dest="alignment",
        help="Pairwise alignment file",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--out",
        dest="out",
        help="Output directory",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        dest="verbose",
        default=False,
        required=False,
        help="Verbose"
    )

    parser.add_argument(
        '--version',
        action='version',
        version=__version__
    )

    # _________________________________________________________________________
    # -------------------------------------------------------------------------
    # get the arguments
    # -------------------------------------------------------------------------
    try:
        options = parser.parse_args()
    except(Exception):
        parser.print_help()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # counter that keeps track of the numbering for each chromosome
    chromosomes_counter = {}
    # dictionaty of file handles for each chromosome
    chromosomes_fh = {}

    # parse the full alignment file
    with open(options.alignment) as fp:
        for name, alignment1, alignment2 in readAlignment(fp):

            # find the chomosome that it belongs
            chromosome = find_chromosome(name)

            # increase the counter for the correct chromosome
            if chromosome not in chromosomes_counter:
                chromosomes_counter[chromosome] = 0
            else:
                chromosomes_counter[chromosome] += 1

            # leave the file handle open for each chromosome
            # for faster writing to the file system
            if chromosome not in chromosomes_fh:
                chromosomes_fh[chromosome] = open(os.path.join(options.out, chromosome + '.axt'), 'w')

            # add the correct counting 
            name = name.split(" ") 
            name[0] = str(chromosomes_counter[chromosome])
            name = " ".join(name)

            # write output files
            chromosomes_fh[chromosome].write(name)
            chromosomes_fh[chromosome].write(alignment1)
            chromosomes_fh[chromosome].write(alignment2)
            chromosomes_fh[chromosome].write(os.linesep)

    # close all open file handles
    for key, value in chromosomes_fh.items():
        chromosomes_fh[key].close()

    # touch file that script is complete
    w = open(os.path.join(options.out, 'Done'), 'w')
    w.write("Done" + os.linesep)
    w.close()


# _____________________________________________________________________________
# -----------------------------------------------------------------------------
# Call the Main function and catch Keyboard interrups
# -----------------------------------------------------------------------------


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt!" + os.linesep)
        sys.exit(0)
