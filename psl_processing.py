#!/usr/bin/env python

import csv
from argparse import ArgumentParser
parser = ArgumentParser()

parser.add_argument("input", help="psl blat results")
parser.add_argument("-o", "--out", dest="outfile",
	help="write results to given filename")
parser.add_argument("-id1", dest="pct_id_cutoff_1",
	default=0.99, type=float,
	help="The percent identity cutoff that determines if the sequence is "+
		    "'identical' to another. Default 0.99.")
# parser.add_argument("-al1", dest="align_len_cutoff_1",
# 	default=59, type=int,
# 	help="The alignment length cutoff, that, with percent identity 1, " +
# 	"determines if the sequence is 'identical' to another. Default 59.")
parser.add_argument("-id2", dest="pct_id_cutoff_2",
	default=0.9, type=float,
	help="The percent identity cutoff that determines if the sequence is "+
		    "'similar' to another. Default 0.9.")
# parser.add_argument("-al2", dest="align_len_cutoff_2",
# 	default=54, type=int,
# 	help="The alignment length cutoff, that, with an percent identity 2, " +
# 	"determines if the sequence is 'similar' to another. Default 54.")
parser.add_argument("-range", dest="pct_id_range",
		    type=float, nargs=3, help="The range of percent identities" +
		    " to test, in the format 'start' 'stop' 'step'. This " +
		    "calculates the appropriate alignment length for each out of 60. [NOT YET IMPLEMENTED]")
		    
parser.add_argument("-q", action="store_true", help="Be quiet.")

args = parser.parse_args()
outfile = args.outfile if args.outfile else args.input+'.out'
al1 = int(60 * args.pct_id_cutoff_1)
al2 = int(60 * args.pct_id_cutoff_2)
print "your input file is {}".format(args.input)
print "your out file is {}".format(outfile)

def frange(args):
	x, y, step = args
	while x > y:
		yield x
		x -= step

def parserow(row):
	row = row.split('\t')
	length = float(row[12]) - float(row[11])
	pct_id = float(row[0])/length
	accn_tmp = row[9].split('_')
	accn, idx = accn_tmp[0], accn_tmp[2]
	return (accn, idx, pct_id, length)

with open(outfile, 'w') as out:
	writer = csv.writer(out, delimiter='\t')
	if args.pct_id_range:
		writer.writerow(["accn", "idx"] + [pct for pct in frange(args.pct_id_range)])
	else:
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
			if args.pct_id_range:
				for pct_id_tmp in frange(args.pct_id_range):
					align_len_tmp = int(60*pct_id_tmp)
					if ((pct_id >= pct_id_tmp) and
					    (align_len >= align_len_tmp)):
						pass
			if ((pct_id >= args.pct_id_cutoff_1) and 
				(align_len >= al1)):
				identical += 1
			elif ((pct_id >= args.pct_id_cutoff_2) and 
				(align_len >= al2)):
				similar += 1
