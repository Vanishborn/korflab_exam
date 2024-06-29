#! /usr/bin/env python
# List found kmer and their locations

import argparse
import sys
import gzip


def read_fasta(filename):
	"""Iteratively read records from a FASTA file"""
	"""'Borrowed' from the MCB185 library"""
	try:
		if filename == '-':
			fp = sys.stdin
		elif filename.endswith('.gz'):
			fp = gzip.open(filename, 'rt')
		else:
			fp = open(filename)
	except Exception as e:
		sys.exit(f"Error opening file {filename}: {e}")

	name = None
	seqs = []
	while True:
		line = fp.readline()
		if line == '':
			break
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				yield (name, ''.join(seqs))
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield (name, ''.join(seqs))
	fp.close()


def anti_seq(seq):
	"""Get the reverse-complement of a sequence"""
	comp = str.maketrans('ACGTRYMKWSBDHVacgtrymkwsbdhv',
						 'TGCAYRKMWSVHDBtgcayrkmwsvhdb')
	anti = seq.translate(comp)[::-1]
	return anti


def find_kmers(seq, kmer_size):
	kmers = {}
	if len(seq) < kmer_size:
		return kmers

	kmer = seq[:kmer_size]
	if kmer not in kmers:
		kmers[kmer] = []
	kmers[kmer].append(str(1))

	for i in range(1, len(seq) - kmer_size + 1):
		kmer = kmer[1:] + seq[i + kmer_size - 1]
		if kmer not in kmers:
			kmers[kmer] = []
		kmers[kmer].append(str(i + 1))

	return kmers


def process_fasta_file(input_file, kmer_size, both_strands):
	for defline, seq in read_fasta(input_file):
		kmers = find_kmers(seq, kmer_size)
		if both_strands:
			rev_seq = anti_seq(seq)
			kmers_rev = find_kmers(rev_seq, kmer_size)
			for kmer, positions in kmers_rev.items():
				if kmer not in kmers:
					kmers[kmer] = []
				for p in positions:
					kmers[kmer].append(str(-int(p)))
		return kmers


def main():
	parser = argparse.ArgumentParser(
		description='Find k-mer locations in a DNA sequence.',
		formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('input_file', help='Input FASTA file')
	parser.add_argument('-k', "--kmer_size",
						type=int,
						default=3,
						help=('Length of k-mers\n'
							  'Default = 3'))
	parser.add_argument('-b', '--both_strands',
						action='store_true',
						help='Count k-mers on both strands')

	args = parser.parse_args()

	kmer_locations = process_fasta_file(
		args.input_file, args.kmer_size, args.both_strands)
	for kmer, positions in sorted(kmer_locations.items()):
		position_all = ' '.join(positions)
		print(f"{kmer} {position_all}")


if __name__ == '__main__':
	main()