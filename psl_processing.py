#!/usr/bin/env python

import csv
from argparse import ArgumentParser
parser = ArgumentParser()

parser.add_argument("input", help="psl blat results")
parser.add_argument("-o", "--out", dest="outfile",
	help="write results to given filename")
parser.add_argument("-id1", dest="pct_id_cutoff_1",
	default=0.99, type=float,
	help="The percent identity cutoff, that, with alignment length 1, " +
	"determines if the sequence is 'identical' to another.")
parser.add_argument("-al1", dest="align_len_cutoff_1",
	default=59, type=int,
	help="The alignment length cutoff, that, with percent identity 1, " +
	"determines if the sequence is 'identical' to another.")
parser.add_argument("-id2", dest="pct_id_cutoff_2",
	default=0.9, type=float,
	help="The percent identity cutoff, that, with alignment length 2, " +
	"determines if the sequence is 'similar' to another.")
parser.add_argument("-al2", dest="align_len_cutoff_2",
	default=50, type=int,
	help="The alignment length cutoff, that, with an percent identity 2, " +
	"determines if the sequence is 'similar' to another.")
parser.add_argument("-q", action="store_true", help="Be quiet.")

args = parser.parse_args()
outfile = args.outfile if args.outfile else args.input+'.out'
print "your input file is {}".format(args.input)
print "your out file is {}".format(outfile)

def parserow(row):
	row = row.split('\t')
	length = float(row[12]) - float(row[11])
	pct_id = float(row[0])/length
	accn_tmp = row[9].split('_')
	accn, idx = accn_tmp[0], accn_tmp[2]
	return (accn, idx, pct_id, length)

with open(outfile, 'w') as out:
	writer = csv.writer(out, delimiter='\t')
	writer.writerow(["accn", "idx", "identical", "similar"])
	with open(args.input) as psl:
		oldidx = None
		oldaccn = None
		identical = 0
		similar = 0
		for i, row in enumerate(psl):
			accn, idx, pct_id, align_len = parserow(row)
			if not args.q:
				print "{0}\r".format(i),
			if not idx:
				oldidx = idx
				oldaccn = accn
			if idx != oldidx:
				writer.writerow([oldaccn, oldidx, identical, similar])
				oldidx = idx
				oldaccn = accn
				identical = 0
				similar = 0
			if ((pct_id >= args.pct_id_cutoff_1) and 
				(align_len >= args.align_len_cutoff_1)):
				identical += 1
			elif ((pct_id >= args.pct_id_cutoff_2) and 
				(align_len >= args.align_len_cutoff_2)):
				similar += 1
			else:
				print "{} | {}".format(pct_id, align_len)
