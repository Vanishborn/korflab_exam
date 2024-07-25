#! /usr/bin/env python
# Find overlapping features between two GFF files using zone-based approach
# Henry Li

import argparse
import gzip
import time
import os
import sys


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


def chr_filter(features1, features2):
	features1_chr_filtered = {}
	features2_chr_filtered = {}
	for chr in features1:
		if chr in features2:
			features1_chr_filtered[chr] = features1[chr]
			features2_chr_filtered[chr] = features2[chr]
	return (features1_chr_filtered, features2_chr_filtered)


def find_zone_len_marks(features1, features2, num_zones):
	zone_len_marks_by_chr = {}
	for chr in features1:
		zone_len_marks_list = []
		unit_zone_len1 = (features1[chr][-1]["end"] // num_zones) + 1
		unit_zone_len2 = (features2[chr][-1]["end"] // num_zones) + 1
		if unit_zone_len1 >= unit_zone_len2:
			unit_zone_len = unit_zone_len1
		else:
			unit_zone_len = unit_zone_len2
		for i in range(1, num_zones + 1):
			zone_len_marks_list.append(unit_zone_len * i)
		zone_len_marks_by_chr[chr] = zone_len_marks_list
	return zone_len_marks_by_chr


def zoning(features, zone_len_marks, num_zones):
	"""Divide features into zones based on zone length marks and number of zones"""
	zoned_features_by_chr = {}
	for chr, feature_list in features.items():
		zones = {f"ZONE{i}": [] for i in range(1, num_zones + 1)}
		for feature in feature_list:
			start = feature['start']
			end = feature['end']
			for i in range(num_zones):
				zone_start = 1 if i == 0 else zone_len_marks[chr][i - 1] + 1
				zone_end = zone_len_marks[chr][i]
				if start > zone_end:
					continue
				if (start <= zone_end and end >= zone_start):
					zones[f"ZONE{i + 1}"].append(feature)
		zoned_features_by_chr[chr] = zones
	return zoned_features_by_chr


def find_overlap(zoned_features1, zoned_features2):
	"""Find overlapping features between two sets of zoned features"""
	overlaps = set()
	for chr in zoned_features1:
		for zone in zoned_features1[chr]:
			if zone in zoned_features2[chr]:
				for feature1 in zoned_features1[chr][zone]:
					for feature2 in zoned_features2[chr][zone]:
						if feature1['start'] <= feature2['end'] and feature2['start'] <= feature1['end']:
							overlap_start = max(
								feature1['start'], feature2['start'])
							overlap_end = min(feature1['end'], feature2['end'])
							overlap = {
								"chrom": feature1['chrom'],
								"overlap_start": overlap_start,
								"overlap_end": overlap_end,
								"feature1": feature1,
								"feature2": feature2
								}
							overlaps.add(str(overlap))
	return [eval(overlap) for overlap in overlaps]



def write_output(overlaps, output_file):
	"""Formats and writes overlaps to output file"""
	header = 'chr\toverlap_beg\toverlap_end\tbeg1\tend1\tbeg2\tend2\tsource1\tsource2\tfeature_type1\tfeature_type2\tscore1\tscore2\tstrand1\tstrand2\tframe1\tframe2\tattribute1\tattribute2'
	lines = [header]
	overlaps = sorted(overlaps, key=lambda x: (x['overlap_start'], x['overlap_end']))
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

	with open(output_file, 'w') as fp:
		fp.write(output_content)


"""argparse statements"""
parser = argparse.ArgumentParser(
	description='Find overlapped features between two GFF files using zone-based approach.',
	formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('gff1',
					help='First input GFF file')
parser.add_argument('gff2',
					help='Second input GFF file')
parser.add_argument('-z', '--zones', type=int, default=10,
					help=('Number of zones to divide each chromosome into\n'
						  'Default = 10'))
parser.add_argument('-o', '--output', type=str,
					help=('Output TSV file\n'
						  'Default = [gff1_basename].[gff2_basename].overlap.tsv'))

args = parser.parse_args()

if args.zones == 0:
	parser.error("The number of zones must be a non-zero positive integer.")

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
	base1 = os.path.splitext(os.path.basename(args.gff1))[0]
	base2 = os.path.splitext(os.path.basename(args.gff2))[0]
	output_file = f"{base1}.{base2}.overlap.tsv"

"""Code body"""
start_time = time.time()
features1 = read_gff(args.gff1)
features2 = read_gff(args.gff2)

features1, features2 = chr_filter(features1, features2)

zone_len_marks = find_zone_len_marks(features1, features2, args.zones)

zoned_features1 = zoning(features1, zone_len_marks, args.zones)
zoned_features2 = zoning(features2, zone_len_marks, args.zones)

overlaps = find_overlap(zoned_features1, zoned_features2)
write_output(overlaps, output_file)
end_time = time.time()

print(f"Overlap gff features with {args.zones} zones completed in {end_time - start_time} seconds")
