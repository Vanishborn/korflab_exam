#! /usr/bin/env python
# Turn output.txt into tsv for plotting
# Henry Li

import argparse


parser = argparse.ArgumentParser(
	description='Extract num of zones and time taken from output file.')
parser.add_argument('-i', '--input', required=True, help='Input text file.')
parser.add_argument('-o', '--output', required=True, help='Output TSV file.')
args = parser.parse_args()

with open(args.input, 'r') as f_in, open(args.output, 'w') as f_out:
	for line in f_in:
		if not line.startswith("###"):
			continue
		words = line.split()
		zones = words[4]
		seconds = words[8]
		f_out.write(f'{zones}\t{seconds}\n')
