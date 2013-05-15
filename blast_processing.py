#!/usr/bin/env python
"""
Usage:
Enter `python blast_processing.py -h` to see usage description.

Background:
This script summarizes the results from all-vs-all megablast (see README for
more info). Specifically, for each sequence, it identifies other sequences
that were at or above certain thresholds for both percent_identity and
alignment_length scores. For example, we may define a sequence to be mostly
'identical' to another if it shares 99% identity and aligns across all 60
nucleotides. We may say it is 'similar' if it shares 95% identity and aligns
more than 40 nucleotides. These parameters are controllable through the command
line options.

Author:
Erik Clarke (ecl@mail.med.upenn.edu)
"""

import gzip
import csv
from argparse import ArgumentParser
parser = ArgumentParser()

parser.add_argument("input", help="tab-delimited, gzip'd blast results")
parser.add_argument("-o", "--out", dest="outfile",
	help="write results to given filename")
parser.add_argument("-id1", dest="pct_id_cutoff_1",
	default=99.00, type=float,
	help="The percent identity cutoff, that, with alignment length 1, " +
	"determines if the sequence is 'identical' to another.")
parser.add_argument("-al1", dest="align_len_cutoff_1",
	default=60, type=int,
	help="The alignment length cutoff, that, with percent identity 1, " +
	"determines if the sequence is 'identical' to another.")
parser.add_argument("-id2", dest="pct_id_cutoff_2",
	default=90.00, type=float,
	help="The percent identity cutoff, that, with alignment length 2, " +
	"determines if the sequence is 'similar' to another.")
parser.add_argument("-al2", dest="align_len_cutoff_2",
	default=57, type=int,
	help="The alignment length cutoff, that, with an percent identity 2, " +
	"determines if the sequence is 'similar' to another.")
parser.add_argument("-q", action="store_true", help="Be quiet.")

args = parser.parse_args()
outfile = args.outfile if args.outfile else args.input+'.out'
print "your input file is {}".format(args.input)
print "your out file is {}".format(outfile)

with open(outfile, 'w') as out:
	writer = csv.writer(out, delimiter = '\t')
	writer.writerow(["accn", "identical", "similar"])
	with gzip.open(args.input) as f:
		oldquery = None
		identical = 0
		similar = 0
		for i, line in enumerate(f):
			if not args.q:
				print "{0}\r".format(i),
			if line.startswith("#") or not line.startswith("GTSP"):
				continue
			query, subject, pct_id, align_len, mismatches, g_open, q_start, \
			q_end, s_start, s_end, e_val, bit_score = line.split('\t')
			pct_id = float(pct_id)
			align_len = int(align_len)
			if not oldquery:
				oldquery = query
			if query != oldquery:
				writer.writerow([oldquery, identical, similar])
				oldquery = query
				identical = 0
				similar = 0
			if ((pct_id >= args.pct_id_cutoff_1) and 
				(align_len >= args.align_len_cutoff_1)):
				identical += 1
			elif ((pct_id >= args.pct_id_cutoff_2) and 
				(align_len >= args.align_len_cutoff_2)):
				similar += 1

