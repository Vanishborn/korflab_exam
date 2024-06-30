#! /usr/bin/env python
# Find overlapping features between two GFF files
# Henry Li

import sys
import os
import gzip
import argparse


def read_gff(filename):
	"""Store features to dictionary by chromosome from a GFF file"""
	features = {}
	try:
		if filename.endswith('.gz'):
			fp = gzip.open(filename, 'rt')
		else:
			fp = open(filename)
	except Exception as e:
		sys.exit(f"Error opening file {filename}: {e}")

	while True:
		line = fp.readline()
		if line == '':
			break
		line = line.strip()
		if line.startswith('#'):
			continue
		cols = line.split()
		if len(cols) < 9:
			continue
		chrom = cols[0]
		start = int(cols[3])
		end = int(cols[4])
		feature_info = {
			"chrom": chrom,
			"source": cols[1],
			"feature_type": cols[2],
			"start": start,
			"end": end,
			"score": cols[5],
			"strand": cols[6],
			"frame": cols[7],
			"attribute": cols[8]
		}
		if chrom not in features:
			features[chrom] = []
		features[chrom].append(feature_info)

	fp.close()
	return features


def find_overlap(features, gff2_file):
	"""Overlaps each feature from gff2 with gff1 features"""
	"""Returns a list of dictionaries"""
	overlaps = []
	try:
		if gff2_file.endswith('.gz'):
			fp = gzip.open(gff2_file, 'rt')
		else:
			fp = open(gff2_file)
	except Exception as e:
		sys.exit(f"Error opening file {gff2_file}: {e}")

	while True:
		line = fp.readline()
		if line == '':
			break
		line = line.strip()
		if line.startswith('#'):
			continue
		fields = line.split()
		if len(fields) < 9:
			continue
		chrom = fields[0]
		start2 = int(fields[3])
		end2 = int(fields[4])
		feature2_info = {
			"chrom": chrom,
			"source": fields[1],
			"feature_type": fields[2],
			"start": start2,
			"end": end2,
			"score": fields[5],
			"strand": fields[6],
			"frame": fields[7],
			"attribute": fields[8]
		}
		if chrom in features:
			for feature1_info in features[chrom]:
				start1 = feature1_info["start"]
				end1 = feature1_info["end"]
				if start1 <= end2 and start2 <= end1:
					overlap_start = max(start1, start2)
					overlap_end = min(end1, end2)
					overlaps.append({
						"chrom": chrom,
						"overlap_start": overlap_start,
						"overlap_end": overlap_end,
						"feature1": feature1_info,
						"feature2": feature2_info
					})

	fp.close()
	return overlaps


def write_output(overlaps, output_file):
	"""Formats and writes to output file"""
	header = 'chr\toverlap_beg\toverlap_end\tbeg1\tend1\tbeg2\tend2\tsource1\tsource2\tfeature_type1\tfeature_type2\tscore1\tscore2\tstrand1\tstrand2\tframe1\tframe2\tattribute1\tattribute2'
	lines = [header]
	for overlap in overlaps:
		f1 = overlap['feature1']
		f2 = overlap['feature2']
		line = (f"{overlap['chrom']}\t{overlap['overlap_start']}\t{overlap['overlap_end']}\t"
				f"{f1['start']}\t{f1['end']}\t{f2['start']}\t{f2['end']}\t"
				f"{f1['source']}\t{f2['source']}\t{f1['feature_type']}\t{f2['feature_type']}\t"
				f"{f1['score']}\t{f2['score']}\t{f1['strand']}\t{f2['strand']}\t"
				f"{f1['frame']}\t{f2['frame']}\t{f1['attribute']}\t{f2['attribute']}")
		lines.append(line)

		output_content = '\n'.join(lines)

		with open(output_file, 'w') as f:
			f.write(output_content)


def main():
	"""argparse statements"""
	parser = argparse.ArgumentParser(
		description='Find overlapped features between two GFF files.',
		formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('gff1', 
					 	help='First input GFF file')
	parser.add_argument('gff2', 
					 	help='Second input GFF file')
	parser.add_argument('-o', '--output', 
						help=('Output TSV file, file name length < 256\n'
							  'Default = [gff1_basename].[gff2_basename].overlap.tsv'))

	args = parser.parse_args()

	"""Input checks"""
	if not os.path.exists(args.gff1):
		sys.exit(f"Error: Input gff1 {args.gff1} does not exist.")

	if not args.gff1.endswith('.gff') and not args.gff1.endswith('.gff.gz'):
		sys.exit("Error: Input gff1 type error.\nFile type gff/gff.gz expected.")

	if not os.path.exists(args.gff2):
		sys.exit(f"Error: Input gff2 {args.gff2} does not exist.")

	if not args.gff2.endswith('.gff') and not args.gff2.endswith('.gff.gz'):
		sys.exit("Error: Input gff2 type error.\nFile type gff/gff.gz expected.")

	"""Format outfile name"""
	if args.output:
		output_file = args.output
	else:
		base1 = os.path.splitext(args.gff1)[0]
		base2 = os.path.splitext(args.gff2)[0]
		output_file = f"{base1}.{base2}.overlap.tsv"

	"""Code body"""
	features = read_gff(args.gff1)
	overlaps = find_overlap(features, args.gff2)
	write_output(overlaps, output_file)


if __name__ == '__main__':
	main()
