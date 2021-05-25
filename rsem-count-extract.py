#!/usr/bin/env python
#coding: utf-8
import sys
from itertools import islice  
mydict = {}
for a in sys.argv[1:]:
    c = open(a, 'r')
    for line in islice(c, 1, None):
        a = line.strip().split()
        key = a[0]
        value = a[4]
        # 4:count; 5:TPM; 6:FPKM
        if key in mydict:
            mydict[key] = mydict[key] + '\t' + value
        else:
            mydict[key] = value
for key,value in mydict.items():
	print(key + '\t' + value)