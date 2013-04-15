#!/usr/bin/python

# File: rows.py
#  Cre: 2012-03
#  Mod: $Date: 2012/03/30 12:17:08 $

# Draws rows from the yields files.  
# Usage (example): ./rows.py data/yields.ww95bis.Z0

import sys 
datafile = open(sys.argv[1],"r")
data = datafile.readlines()
m = data[0].split()
c = data[6].split()
o = data[8].split()
fe = data[13].split()
for i in range(len(m)):
    print m[i], c[i], o[i], fe[i]



