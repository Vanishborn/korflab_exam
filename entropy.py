#! /usr/bin/env python
# Mask low entropy regions of nucleotide sequences
# Henry Li

import sys
import gzip
import math
import argparse
import os


def read_fasta(filename):
	"""Iteratively read records from a FASTA file"""
	"""Stolen from the MCB185 library"""
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


def entropy_of_seq(seq):
	"""Calculate entropy of a given DNA seq"""
	nts = 'ACGT'
	entropy = 0
	for nt in nts:
		p = seq.count(nt) / len(seq)
		if p > 0:
			entropy += -p * math.log2(p)
	return entropy


def mask_seq(seq, window_size, threshold, soft_mask=False):
	"""Perform masking of windows with entropy less than threshold"""
	seq_list = list(seq)
	for i in range(len(seq) - window_size + 1):
		window_seq = seq[i:i + window_size]
		entropy = entropy_of_seq(window_seq)
		if entropy < threshold:
			for j in range(i, i + window_size):
				if soft_mask:
					seq_list[j] = seq_list[j].lower()
				else:
					seq_list[j] = 'N'
	return ''.join(seq_list)


def process_fasta_file(input_file, output_file, window_size, threshold, soft_mask):
	"""Read in raw fasta file, output masked fasta file"""
	try:
		with open(output_file, 'w') as out_file:
			for defline, seq in read_fasta(input_file):
				masked_seq = mask_seq(seq, window_size, threshold, soft_mask)
				out_file.write(f'>{defline}\n')
				for i in range(0, len(masked_seq), 60):
					out_file.write(masked_seq[i:i+60] + '\n')
	except Exception as e:
		sys.exit(f"Error processing FASTA file {input_file}: {e}")


def main():
	"""argparse statements"""
	parser = argparse.ArgumentParser(
		description='Entropy filter for DNA sequences.',
		formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('input_file',
						help='Input FASTA file')
	parser.add_argument('-o', '--output',
						help=('Output FASTA file, file name length < 256\n'
							  'Default = [input_file_basename].masked.fasta'))
	parser.add_argument('-w', '--window_size',
						type=int,
						default=11,
						help=('Window size\n'
							  'Default = 11'))
	parser.add_argument('-t', '--threshold',
						type=float,
						default=1.4,
						help=('Entropy threshold\n'
							  'Default = 1.4'))
	parser.add_argument('-s', '--soft_mask',
						action='store_true',
						help=('Use soft masking\n'
							  'Masked to lowercase instead of N'))

	args = parser.parse_args()

	"""Input checks"""
	if not os.path.exists(args.input_file):
		sys.exit(f"Error: Input file {args.input_file} does not exist.")

	if not args.input_file.endswith('.fasta') and not args.input_file.endswith('.fasta.gz') and not args.input_file.endswith('.fa') and not args.input_file.endswith('.fa.gz'):
		sys.exit(
			"Error: Input file type error.\nFile type fasta/fa/fasta.gz/fa.gz expected.")

	"""Format outfile name"""
	if args.output:
		if len(args.output) >= 256:
			sys.exit(
				"Error: Output file name exceeds the maximum length of 255 characters.")
		output_file = args.output
	else:
		output_file = os.path.splitext(args.input_file)[0] + '.masked.fasta'

	"""Code body"""
	process_fasta_file(args.input_file, output_file,
					   args.window_size, args.threshold, args.soft_mask)


if __name__ == '__main__':
	main()
