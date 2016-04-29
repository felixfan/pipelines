# -*- coding: utf-8 -*-

import sys

dat = sys.argv[1]
out = sys.argv[2]

f = open(dat)
fw = open(out,'w')

for r in f:
	r = r.strip()
	arr = r.split()
	fw.write(arr[0])
	for i in arr[1:6]:
		fw.write('\t')
		fw.write(i)
	for i in arr[6:]:
		fw.write('\t')
		if i == 'N':
			fw.write('0')	
		else:
			fw.write(i)
	fw.write('\n')
f.close()
fw.close()
