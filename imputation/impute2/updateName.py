# -*- coding: utf-8 -*-

import sys

ped = sys.argv[1]
out = sys.argv[2]

f = open(ped)
fw = open(out,'w')

for r in f:
	r = r.strip()
	arr = r.split()
	fw.write(arr[1])
	fw.write('\t')
	if -1 != arr[1].find(':'):
		ar = arr[1].split(':')
		if ar[0].startswith('rs'):
			fw.write(ar[0])
		else:
			fw.write(ar[0]+':'+ar[1])
	else:
		fw.write(arr[1])
	fw.write('\n')
f.close()
fw.close()
