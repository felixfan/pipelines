import subprocess
import os

f = open("tsvFiles.txt")
f2 = open("id.txt","w")

for d in f:
	d = d.rstrip()
	words = d.split("/")
	cm = "cp " + d + " temp"
	os.system(cm)
	os.system("sed -i '1d' temp")
	os.system("awk '{print $2,$3,$4,$5,$6,$7}' temp >> optitype.hla.txt")
	f2.write(words[2])
	f2.write("\n")
f.close()
f2.close()


