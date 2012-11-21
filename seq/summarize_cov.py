#!/usr/bin/env python

# summarize_cov.py
# Uses bedtools cov.out file (provided on command line), and summarizes coverage depth in ranges
# 
# Coverage depth ranges default to: 1,10,30,100 or can be provided as option (eg. -c 1,10,25,50,99)
#
# Created by Sue Grimes on 11/01/12.
# Copyright (c) 2012 Ji Research Group - Stanford Genome Technology Center. All rights reserved.

import sys, os, csv, re
from optparse import OptionParser
script_name = os.path.basename(__file__)

# Derive column headings for output file
def header_cols():
  depth_hdr_cols = []
  for i, dcol in enumerate(DEPTH_COLS):
    if dcol == 1:
      depth_hdr_cols.append('0s')
    else:
      depth_hdr_cols.append('< ' + str(dcol))
    if i == IDX_DEPTH_MAX-1:
      depth_hdr_cols.append(str(dcol) + '+')
  return ["Chr", "Loc Beg", "Loc End"] + depth_hdr_cols + ["Total Bases"] + depth_hdr_cols + ["Ave Cov"] 

def cov_range_callback(option, opt, value, parser):
  cov_array = [int(r_val) for r_val in value.split(',')]
  setattr(parser.values, option.dest, cov_array)

# Determine which range the specified depth is within, and return index to totals array accordingly
def categorize_depth(depth):
  for i, dcol in enumerate(DEPTH_COLS):
    if depth < dcol:
      return i

  return IDX_DEPTH_MAX

# Calculate proportion of bases covered at each depth range, and overall average coverage
def calc_averages(cov_totals, sum_depth):
  cov_averages = [0] * NR_COLS
  
  for i in range(IDX_DEPTH_MAX+1):
    cov_averages[i] = round_num((float(cov_totals[i]) / cov_totals[IDX_DEPTH_TOT]), 3)
  
  cov_averages[IDX_DEPTH_TOT] = round_num((float(sum_depth) / cov_totals[IDX_DEPTH_TOT]), 1)
  return cov_averages

# Rounding method 
def round_num(nr, precision=3):
  return "{0:.{1}f}".format(nr, precision)
 
# MAIN PROGRAM #
# Define and read command line options and arguments
usage = "usage: %prog <cov.out file>"
parser = OptionParser(usage)

parser.add_option("-c", "--covrange", dest="cov_range", type="string", default=[1,10,30,100],
                  action="callback", callback=cov_range_callback, help="Comma separated list for coverage range.  Default: 1,10,30,100")
(options, args) = parser.parse_args()

# Ensure that valid input file is provided as command line argument
if len(args) != 1:
  print "Incorrect number of arguments"
  parser.print_help()
  sys.exit(1)

elif args[0][-7:] != 'cov.out':
  print "File extension must be .cov.out"
  sys.exit(1)

try:
  in_file = args[0]
  cov_in = open(in_file, 'r')
except:
  print "Unable to open input file", in_file
  sys.exit(1)
  
try:
  out_file = in_file.replace('cov.out', 'cov_by_region.txt')
  cov_out = open(out_file, 'wb')
except:
  print "Unable to open output file ", out_file
  sys.exit(1)

# Initialize constants
# If DEPTH_COLS = [1,10,30,100]: array index values for cov_totals and cov_averages will reflect the following:
# 0 - 0 coverage
# 1 - 1-9 coverage
# 2 - 10-29 coverage
# 3 - 30-99 coverage
# 4 - 100+ coverage
# 5 - Total (or average) 
#
# IDX_DEPTH_MAX = 4 (index to cov_totals array value with maximum depth, eg 100+)
# IDX_DEPTH_TOT = 5 (index to cov_totals array value storing total coverage)
# NR_COLS       = 6 (total number of output columns for depths)
#
print "Coverage ranges: ", options.cov_range
DEPTH_COLS = options.cov_range
IDX_DEPTH_MAX = len(DEPTH_COLS)
IDX_DEPTH_TOT = len(DEPTH_COLS)+1
NR_COLS = len(DEPTH_COLS)+2

# Initialize other variables
cov_line = []

# cov_array will accumulate output in format: [chr, loc_beg, loc_end, cov_totals[], cov_averages[]]
# cov_totals, cov_averages and grand_tots are variable length arrays, depending on # depths specified in DEPTH_COLS
cov_array = []
cov_totals = [0] * NR_COLS
cov_averages = [0] * NR_COLS
grand_tots = [0] * NR_COLS
sum_depth = 0

# Start main loop
ln = 0
for cov_line in cov_in:
  cov_cols = cov_line.rstrip("\n").split("\t")

  chr_num = cov_cols[0]
  loc_beg = int(cov_cols[1])
  loc_end = int(cov_cols[2])
  rel_base = cov_cols[3]
  depth   = int(cov_cols[4])
  idx_depth = categorize_depth(depth)
  
  if ln == 0:
    cov_array.append([chr_num, loc_beg, loc_end, cov_totals, cov_averages])

  if chr_num == cov_array[-1][0] and loc_beg == cov_array[-1][1]:
    # Increment totals for current chromosome, start pos
    cov_totals[idx_depth] += 1
    cov_totals[IDX_DEPTH_TOT] += 1  
    sum_depth += depth

  else:
    # Write totals to array, and calculate averages (for previous chromosome, start pos)
    cov_array[-1][3] = cov_totals
    cov_array[-1][4] = calc_averages(cov_totals, sum_depth)
    
    # Reset variables to accumulate totals for new chromosome, position (from line already read)
    cov_totals = [0] * NR_COLS
    cov_totals[idx_depth] = 1
    cov_totals[IDX_DEPTH_TOT] = 1
    cov_averages = [0] * NR_COLS
    sum_depth = depth

    # Write new row in coverage array
    cov_array.append ([chr_num, loc_beg, loc_end, cov_totals, cov_averages])
  
  ln += 1

# EOF: Calculate averages and write totals to final row to coverage array
cov_array[-1][3] = cov_totals
cov_array[-1][4] = calc_averages(cov_totals, sum_depth)

# Write summarized coverage array to tab-delimited output file
cov_csv = csv.writer(cov_out, delimiter='\t')
cov_csv.writerow(header_cols())

for i, cov_line in enumerate(cov_array):
  cov_csv.writerow(cov_line[0:3] + cov_line[3] + cov_line[4])
  
  for i in range(NR_COLS):
    grand_tots[i] += cov_line[3][i]

cov_csv.writerow(['','','Total:'] + grand_tots)
cov_out.close
