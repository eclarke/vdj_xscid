#!/usr/bin/env python

"""igblastn_split.py
Splits a list of sequences and runs igblastn against each chunk.

Usage:
    igblastn_split.py [options] <fasta>

Options:
    -h, --help               show this message and exit
    -v, --verbose            be verbose
    --chunks=CHUNKS          specify # of chunks (by default, number of cores on machine)
    --batch=BATCHSIZE        specify # of jobs launched simultaneously (by default, equal
                                 to the number of chunks)
    --igblast=IGBLAST        specify the location of the igblast folder [default: ~/igblast]
    --domain=DOMAIN          specify the domain system [default: kabat]

Author:
    Erik Clarke [ecl@mail.med.upenn.edu]
"""

from docopt import docopt
import multiprocessing
import os
import sys
from glob import glob

try:
    # the backport of subprocess from python3.2 contains bugfixes for long-
    # running processes such as igblast
    import subprocess32 as subprocess
except ImportError:
    import subprocess


def main():
    args = docopt(__doc__)
    verbose = args['--verbose']
    fasta = os.path.abspath(args['<fasta>'])
    igdir = os.path.abspath(os.path.expanduser(args['--igblast']))
    igdb = igdir + '/database'
    chunks = args['--chunks'] if args['--chunks'] else multiprocessing.cpu_count()
    batchsize = args['--batch'] if args['--batch'] else chunks
    igblast = ['bin/igblastn',
               '-germline_db_V', igdb+'/human_gl_V',
               '-germline_db_D', igdb+'/human_gl_D',
               '-germline_db_J', igdb+'/human_gl_J',
               '-organism', 'human',
               '-domain_system', args['--domain'],
               '-auxiliary_data', igdir+'/data/human_gl.aux',
               '-outfmt', '7',
               '-query']
    try:
        subprocess.check_call(['split', '-n', 'l/' + str(chunks), fasta,
                                 fasta + '.part.'])
    except subprocess.CalledProcessError:
        sys.stderr.write("Error in splitting file: does the file exist and " +
                         "can we write to the directory?")
        sys.exit(1)
    parts = glob(fasta + '.part.*')
    if len(parts) != chunks:
        # This could happen because either the split didn't work or
        # there's more than what we split in the directory (maybe left
        # over from another run?
        raise Exception("Did not find %i file(s) to process..." % chunks)
    b = 0
    processes = []
    for part in parts:
        b += 1
        if verbose:
            print("+ igblasting {part}\t[{i}/{total}]".format(part=part,
                                                              i=b,
                                                              total=len(parts)))
        processes.append(subprocess.Popen(igblast + [part, '-out', part+'.out'],
                                          cwd=igdir))
        if b % batchsize == 0:
            if verbose and b != len(parts):
                print("Batch limit reached, waiting to launch new processes...")
            for p in processes:
                p.wait()

    try:
        subprocess.check_call('cat {0}.part.*.out > {0}.igblast'.format(fasta),
                              shell=True)
        subprocess.check_call('rm -f {}.part.*'.format(fasta), shell=True)
        if verbose:
            print("Done. File outputted to {}.igblast".format(fasta))
    except subprocess.CalledProcessError:
        sys.stderr.write("Problem concatenating files back together, you'll " +
                         "have to do it manually.")
        sys.exit(0)

if __name__ == '__main__':
    main()
