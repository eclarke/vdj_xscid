"""igblast_processing.py
For each igblastn query, outputs the number of V, D, and J genes found in the
sequence, as well as if all three were found (i.e., it was 'complete'), as a 
tab-delimited file for integration into R.

Usage:
    igblast_processing.py [options] <input>

Options:
    -h --help           show this message and exit
    -o, --out=FILE      specify output filename
    -v, --verbose       be verbose
    --threshold=ID      % identity threshold for considering a match valid
                            [default: 80.00]
    -x, --expanded      produce expanded output (experimental)

Author:
    Erik Clarke [ecl@mail.med.upenn.edu]
"""

from docopt import docopt
import csv, sys

"""
Writes the matches to the output. Adds a 'complete' marker if all three
genes were found.
"""
def write_matches(accn, idx, matches, writer):
    if accn and idx:
        complete = (1 if matches['V'] > 0 and matches['D'] > 0 and 
                    matches['J'] else 0)
        writer.writerow([accn, idx, matches['V'], matches['D'], matches['J'],
                         complete])

"""
Returns "GTSP0019_1_3001" as ('GTSP0019', '3001')
"""
def split_query(_query):
    query = _query.strip(' \n\t').split('_')
    accn, idx = query[0], query[2]
    return (accn, idx)

"""
The main loop here iterates through the file. Every time it encounters a
new igblastn result (demarcated by "# Query: GTSP..."), it writes the previous
query results to the file and resets the matches.
"""
def main():
    args        = docopt(__doc__)
    infile      = args['<input>']
    outfile     = args['--out'] if args['--out'] else infile + '.processed'
    threshold   = float(args['--threshold'])
    expanded    = args['--expanded']

    with open(outfile, 'w') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['accn', 'idx', 'v.genes', 'd.genes', 'j.genes',
                         'complete'])
        with open(infile) as results:
            oldidx  = None
            oldaccn = None
            matches = {"V":0, "D":0, "J":0, "complete":False}
            for i, row in enumerate(results):
                if "# Query: GTSP" in row:
                    write_matches(oldaccn, oldidx, matches, writer)
                    accn, idx = split_query(row[9:-1])
                    oldaccn, oldidx = accn, idx
                    matches = {"V":0, "D":0, "J":0}
                elif row.startswith("#"):
                    continue
                elif row.startswith(('V\t', 'D\t', 'J\t')):
                    row = row.split('\t')
                    gene_class = row[0]
                    pct_id = float(row[3])
                    if pct_id >= threshold:
                        matches[gene_class] += 1

if __name__ == '__main__':
    main()
