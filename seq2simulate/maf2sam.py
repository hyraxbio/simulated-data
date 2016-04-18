from itertools import islice
import shutil


def parse_maf(ref, seq, ref_name="HXB2_pol"):
	""" The MAF format is:
		a
		s ref START_POSITION LENGTH ORIENTATION TOTAL_LENGTH SEQUENCE
		s READ_ID START_POSITION LENGTH ORIENTATION TOTAL_LENGTH SEQUENCE
		This function takes the second two lines above and returns a single 
		line of SAM format.

		Args:

			ref: reference sequence
			seq: query sequence
	"""
	ref_array = ref.split()
	seq_array = seq.split()

	if len(ref_array) < 7 or len(seq_array) < 7:
		return None

	name = seq_array[1]

	flag = 8
	if seq_array[4] == '-':
		flag += 16

	begin = int(ref_array[2])

	mate_ref = '='
	mate_location = 0
	ref_len = int(ref_array[3])
	query_len = int(seq_array[3])

	ungapped = ''.join([s for s in seq_array[6] if s != '-'])
	qual = 'I' * len(ungapped)

	cigar = ''
	count = 0
	current = 'X'
	for r, s in zip(ref_array[6], seq_array[6]):
		new = 'M'
		if r == '-':
			new = 'I'
		elif s == '-':
			new = 'D'

		if count > 0 and current != new:
			cigar += '%d%s' % (count, current)
			count = 1
		else:
			count += 1

		current = new

	return "\t".join([str(s) for s in [name, flag, ref_name, begin, mate_ref, 
			cigar, mate_location, 
			ref_len, query_len, ungapped, qual]])

def maf2sam(in_file, ref_file):
	""" Take in a MAF input file and a reference fasta and produce a .sam
		file.

		Args:
			in_file: the input file (must by in MAF format and end in .maf)
			ref_file: the fasta reference file.
	"""
	ref_name = ''
	ref_length = 0
	with open(ref_file, 'r') as ref:
		line = ref.readline()
		ref_name = line.strip()[1:]
		for line in ref:
			ref_length += len(line)

	sam_file = in_file.replace('.maf', '.sam')
	with open(sam_file, 'w') as sam:
		sam.write("@HD\tVN:1.5\tSO:coordinate\n")
		sam.write("@SQ\tSN: %s\tLN: %d\n" % (ref_name, ref_length))
	
		with open(in_file, 'r') as inf:
			while True:

				a = inf.readline()
				ref = inf.readline()
				seq = inf.readline()
				_ = inf.readline()
				if ref == '' and seq == '':
					break

				sam_line = parse_maf(ref, seq, ref_name)
				if sam_line is not None:
					sam.write(sam_line+"\n")














