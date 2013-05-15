#!/usr/bin/env python

"""PSL File Processor.
Description:
    Reads a .psl file from BLAT and aggregates the number of hits a subject
    sequence has above a set of percent identity and length thresholds.
    Outputs a tab-delimited file with the hits for each threshold next to the
    sequence's accession and unique index number.

Usage: 
    psl_processor.py [options] <input>

Options:
    -h --help                 show this help message and exit
    --max=MAX                 the maximum % identity threshold [default=0.99]
    --min=MIN                 the minimum % identity threshold [default=0.85]
    --step=STEP               spacing of intermediate % id thresholds [default=0.05]
    -o, --output FILE         specify output file [default <input>.processed]
    -v, --verbose             be verbose

Author:
    Erik Clarke [ecl@mail.med.upenn.edu]
    
"""

import csv
import sys
from docopt import docopt

# changes file extension
def ch_ext(fname, ext):
    _fname = fname.split('.')
    if len(_fname) > 1:
        _fname.pop()
        return '.'.join(_fname + [ext])
    else:
        return fname+'.'+ext

def frange(_min, _max, step):
    _min = int(_min * 100)
    _max = int(_max * 100)
    step = int(step * 100)
    while _max >= _min:
        yield _max / 100.0
        _max -= step

def parserow(row):
    row = row.split('\t')
    length = float(row[12]) - float(row[11])
    pct_id = float(row[0]) / length 
    idx = row[9]
    return (idx, pct_id, length)

def main():
    args = docopt(__doc__, version="1.0.0")
    infile = args["<input>"]
    outfile = args["--output"] if args["--output"] else ch_ext(args['<input>'], 'processed')
    _min = float(args['--min']) if args['--min'] else 0.84
    _max = float(args['--max']) if args['--max'] else 0.99
    step = float(args['--step']) if args['--step'] else 0.05
    verbose = args['--verbose']
    thresholds = [(threshold, int(60*threshold)) for threshold 
                  in frange(_min, _max, step)]
    
    min_id, min_al = thresholds[-1]
    if verbose:
        print("Thresholds:")
        [sys.stdout.write("{}|{}\t".format(_id, _al)) for _id, _al in thresholds]
        print("\n")
    if not args['--output'] and verbose:
        print("Your output file is {}\n".format(outfile))
    with open(outfile, 'w') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(["idx"] + ['id.'+str(int(t[0]*100)) for t in thresholds])
        with open(infile) as psl:
            oldidx = None
            matches = {threshold[0]:0 for threshold in thresholds}
            for i, row in enumerate(psl):
                try:
                    idx, pct_id, align_len = parserow(row)
                except (IndexError, ValueError) as e:  # probably due to a header
                    if verbose:
                        print(row, e)
                    continue
                if verbose:
                    print "{0}\r".format(i),
                if not oldidx:
                    oldidx = idx
                if idx != oldidx:
                    writer.writerow([oldidx] + [matches[_id] for _id, _al in thresholds])
                    oldidx = idx
                    matches = {threshold[0]:0 for threshold in thresholds}
                for _id, _al in thresholds:
                    if (pct_id < min_id or align_len < min_al):
                        break
                    if (pct_id >= _id and align_len >= _al):
                        matches[_id] += 1
                        break

if __name__ == '__main__':
    main()
