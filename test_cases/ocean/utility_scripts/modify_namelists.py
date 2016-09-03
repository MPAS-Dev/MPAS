#!/usr/bin/env python
"""
This script modifies one or more namelist entires in one or more namelist
files.  The new namelist entries are supplied as a list of quoted strings
with the -n flag, for example:
-n "config_pio_num_iotasks = 8" "config_pio_stride = 16" \\
   "config_run_duration = '0001-00-00_00:00:00'"
Note that indentation should not be included in the namelist entries.
The namelist files are provided through the -f flag.
"""
import os
import glob
import argparse
from datetime import datetime

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("-n", "--new_entries", dest="new_entries", help="namelist variables and new values separated by equals signs", metavar="ENTRY", nargs="+", required=True)
parser.add_argument("-f", "--namelist_files", dest="namelist_files", help="The namelist files to be modified", metavar="NAMELIST", nargs="+", required=True)

args = parser.parse_args()

entries = {}
for entry in args.new_entries:
   var = entry.split('=')[0].strip()
   entries[var] = entry

for fileName in args.namelist_files:
    lines = []
    inFile = open(fileName, 'r')
    for line in inFile:
        lines.append(line)
    inFile.close()

    outFile = open(fileName, 'w')
    for line in lines:
        for var in entries:
            if var in line:
                line = "    %s\n"%entries[var]
        outFile.write(line)
    outFile.close()

