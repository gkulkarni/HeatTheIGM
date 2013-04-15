#!/usr/bin/python 

# File: process_yields.py
#  Cre: 2012-01
#  Mod: $Date: 2012/10/23 08:10:38 $ ($Revision: 1.2 $) 

# Rewrite Emmanuel's yield tables in the format used by reion.

import sys 
from sets import Set 
filename = sys.argv[1] 
data = open(filename,"r")
i = 0 
for line in data: 
    ln = line.strip()
    if not ln.startswith("#"):
        i = i+1 
        if i not in Set([1,2,5,20,21,22,23]):
            if (i > 23): 
                break
            lnbrk = ln[15:].split()
            print "  ".join(['%10s'% lnbrk[n] for n in xrange(len(lnbrk))])
