#!/usr/bin/env python

"""blat_runner.py
Description:
    Takes a collection of fasta sequences, converts them to a 2Bit file, and 
    runs BLAT in an all-vs-all comparison.

Usage:
    blat_runner.py [options] <fasta>

Options:
    -h --help           show this message and exit
    -o, --out=FILE      specify output file (by default, input + .psl)
    -b, --blat=BLAT     specify the location of the BLAT binary
                            [default: ~/bin/blat]
    -f, --faTo2Bit=FA2  specify the location of the faToTwoBit binary 
                            [default: ~/bin/faToTwoBit]
    -m, --minScore=MS   specify the minimum BLAT score to keep [default: 50]
    -s, --stepSize=SS   specify the step size to use for BLAT [default: 20]
    -v, --verbose       be chatty
    
"""

from docopt import docopt
import subprocess
import os


# changes file extension
def ch_ext(fname, ext):
    fname = fname.split('.')
    fname.pop()
    return '.'.join(fname + [ext])


def main():
    args = docopt(__doc__)
    infile = args['<fasta>']
    blat = os.path.abspath(os.path.expanduser(args['--blat']))
    fa2b = os.path.abspath(os.path.expanduser(args['--faTo2Bit']))
    min_score = args['--minScore']
    step_size = args['--stepSize']
    verbose = args['--verbose']

    _2bit = ch_ext(infile, '2bit')
    cmd = [fa2b, infile, _2bit]
    if verbose:
        print("Converting {} to 2bit file...".format(infile))
        print(' '.join(cmd))
    pid = subprocess.Popen(cmd)
    pid.wait()

    psl = ch_ext(infile, 'psl')
    cmd = [blat, '-minScore={}'.format(min_score), '-repMatch=12321', 
           '-stepSize={}'.format(step_size), '-noHead', _2bit, infile, psl]
    if verbose: 
        print("Running all-vs-all BLAT...")
        print(' '.join(cmd))
    pid = subprocess.Popen(cmd)
    pid.wait()

    if verbose:
        print("Done.")
        print("PSL file: " + psl)


if __name__ == '__main__':
    main()