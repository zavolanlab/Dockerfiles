#!/usr/bin/env python

__version__ = "0.1"
__author__ = "Foivos Gypas"
__contact__ = "foivos.gypas@unibas.ch"
__doc__ = "Download pairwise alignments from ucsc"

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


def main():
    """ Main function """

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=RawTextHelpFormatter
    )

    parser.add_argument(
        "--organism",
        dest="organism",
        help="Reference organism (e.g. hg38, hg19, mm10)",
        required=True,
        metavar="FILE"
    )

    parser.add_argument(
        "--organism_to_download",
        dest="organism_to_download",
        help="Name of pairwise alignment organism",
        required=True
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

    start = "rsync://hgdownload.soe.ucsc.edu/goldenPath/" + options.organism + "/vs"
    f_s_l = options.organism_to_download
    f_s = f_s_l[0].upper() + f_s_l[1:]
    output_dir = os.path.join(options.out, options.organism + "_to_" + f_s_l)
    # file_name = options.organism + "." + f_s_l + ".net.axt.gz"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    download_link = start + f_s + "/axtNet/"

    try:
        if options.verbose:
            sys.stdout.write("Downloading " + download_link + os.linesep)
        os.system("rsync -a -P " + download_link + " " + output_dir)
    except:
        sys.stderr.write("Could not download " + download_link + os.linesep)

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
