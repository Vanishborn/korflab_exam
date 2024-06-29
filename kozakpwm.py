#! /usr/bin/env python
# Create a PWM of the Kozak consensus from a GBFF file

import argparse
import sys
import gzip
import json
import os


def read_gbff(filename):
	"""Read info and join seq from a GBFF file"""
	try:
		if filename == '-':
			fp = sys.stdin
		elif filename.endswith('.gz'):
			fp = gzip.open(filename, 'rt')
		else:
			fp = open(filename)
	except Exception as e:
		sys.exit(f"Error opening file {filename}: {e}")

	info = []
	seqs = []
	origin_found = False

	while True:
		line = fp.readline()
		line = line.rstrip()
		if 'ORIGIN' in line:
			origin_found = True
		if origin_found:
			if '//' in line:
				yield (info, ''.join(seqs))
				break
			seqs.append(''.join(line.split()[1:]))
		else:
			info.append(line)
	yield (info, ''.join(seqs))

	fp.close()


def anti_seq(seq):
	"""Get the reverse-complement of a sequence"""
	comp = str.maketrans('ACGTRYMKWSBDHVacgtrymkwsbdhv',
						 'TGCAYRKMWSVHDBtgcayrkmwsvhdb')
	anti = seq.translate(comp)[::-1]
	return anti


def create_pwm(info, seq, by_nucleotide):
	if by_nucleotide:
		pwm = {'A': [0] * 14, 'C': [0] * 14, 'G': [0] * 14, 'T': [0] * 14}
	else:
		pwm = []
		for _ in range(14):
			pwm.append({'A': 0, 'C': 0, 'G': 0, 'T': 0})

	for line in info:
		if line.startswith('     CDS'):
			if 'join' in line:
				continue

			cds_pos = line.split()[1]

			if 'complement' in cds_pos:
				end_pos = int(cds_pos.strip(')').split('..')[1])
				kozak = anti_seq(seq[end_pos-5:end_pos+9])
			else:
				start_pos = int(cds_pos.split('..')[0])
				kozak = seq[start_pos-10:start_pos+4]

			if len(kozak) != 14:
				continue

			for pos, nt in enumerate(kozak):
				if by_nucleotide:
					pwm[nt.upper()][pos] += 1
				else:
					pwm[pos][nt.upper()] += 1

	return pwm


def pwm_to_json(pwm, output_file):
	try:
		with open(output_file, 'w') as out_file:
			json.dump(pwm, out_file, indent=4)
	except IOError as e:
		sys.exit(f"Error writing to output file {output_file}: {e}")


def main():
	parser = argparse.ArgumentParser(
		description='Create a PWM of the Kozak consensus from a GBFF file.',
		formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('input_file',
						help=('Input GBFF file\n'
							  'Accept .gbff or .gbff.gz'))
	parser.add_argument('-o', '--output',
						help=('Output JSON file, file name length < 256\n'
							  'Default = [input_file_basename].pwm.json'))
	parser.add_argument('-n', '--by_nucleotide',
						action='store_true',
						help=('Output in nucleotide-by-nucleotide format\n'
							  'Default output in position-by-position format'))

	args = parser.parse_args()

	if args.output:
		if len(args.output) >= 256:
			sys.exit(
				"Error: Output file name exceeds the maximum length of 255 characters.")
		output_file = args.output
	else:
		output_file = os.path.splitext(args.input_file)[0] + ".pwm.json"

	if not os.path.exists(args.input_file):
		sys.exit(f"Error: Input file {args.input_file} does not exist.")

	if not args.input_file.endswith('.gbff') and not args.input_file.endswith('.gbff.gz'):
		sys.exit("Error: Input file type error.\nFile type gbff/gbff.gz expected.")

	for info, seq in read_gbff(args.input_file):
		pwm = create_pwm(info, seq, args.by_nucleotide)
		pwm_to_json(pwm, output_file)


if __name__ == '__main__':
	main()
