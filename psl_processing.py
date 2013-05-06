"""PSL File Processor.
Description:
    Aggregates the number matches above two thresholds

Usage: 
	psl_processor.py [options] <input>

Options:
  -h --help                 show this help message and exit
  --max=MAX                 the maximum % identity threshold [default=0.99]
  --min=MIN                 the minimum % identity threshold [default=0.85]
  --step=STEP               spacing of intermediate % id thresholds [default=0.05]
  -o, --output FILE         specify output file [default input + .out]
  -q --quiet                operate in quiet mode

"""

import csv, sys
from docopt import docopt

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
    accn_tmp = row[9].split('_')
    accn, idx = accn_tmp[0], accn_tmp[2]
    return (accn, idx, pct_id, length)

if __name__ == '__main__':
    args = docopt(__doc__, version="1.0.0")
    infile = args["<input>"]
    outfile = args["--output"] if args["--output"] else '.'.join([infile, "out"])
    _min = float(args['--min']) if args['--min'] else 0.84
    _max = float(args['--max']) if args['--max'] else 0.99
    step = float(args['--step']) if args['--step'] else 0.05
    thresholds = [(threshold, int(60*threshold)) for threshold 
                    in frange(_min, _max, step)]
    print thresholds

    print("Your output file is {}".format(outfile))
    print(outfile)
    with open(outfile, 'w') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(["accn", "idx"] + ['id.'+str(int(t[0]*100)) for t in thresholds])
        with open(infile) as psl:
            oldidx = None
            oldaccn = None
            matches = {threshold[0]:0 for threshold in thresholds}
            for i, row in enumerate(psl):
                accn, idx, pct_id, align_len = parserow(row)
                if not args["-q"]:
                    print "{0}\r".format(i),
                if not idx:
                    oldidx = idx
                    oldaccn = accn
                if idx != oldidx:
                    writer.writerow([oldaccn, oldidx, identical, similar])
                    oldidx = idx
                    oldaccn = accn
                    matches = {threshold[0]:0 for threshold in thresholds}
                for _id, _al in thresholds:
                    if (pct_id >= _id and align_len >= _al):
                        matches[_id] += 1
                        continue
