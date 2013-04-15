#!/usr/bin/python 

# File: process_yields.py
#  Cre: 2012-01-18
#  Mod: $Date: 2012/10/23 08:10:38 $ ($Revision: 1.1 $) 

# Transposes data files. Useful if you want to plot row against row. 

import sys
filename = sys.argv[1]
data = open(filename,"r")
i = 0
all_data = []
for line in data:
    ln = line.strip()
    row = ln.split()
    all_data.append(row)
transposed_data = map(list, zip(*all_data))    
for x in transposed_data:
    print "  ".join(['%6s'% x[n] for n in xrange(len(x))])

