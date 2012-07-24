#! /usr/bin/env python
import os
import sys
CASTORDIR = sys.argv[1]
FILTER = sys.argv[2]
os.system("rfdir "+CASTORDIR+" | grep default > tmp.list")
list = open("tmp.list")
mylist = open(FILTER+".list","w")
for line in list:
    sample = CASTORDIR+"/"+line.split(" ")[-1][:-1]
    os.system("rfdir "+sample+" | grep root >> tmp2.list")
    list2 = open("tmp2.list")
    for line2 in list2:
        mylist.write((sample+"/"+line2.split(" ")[-1]).replace("//","/"))
        continue
    list2.close()
    os.system("rm tmp2.list")
    continue
mylist.close()
list.close()
os.system("rm tmp.list")
if(FILTER == "NoFile"):
    os.system("wc "+FILTER+".list")
    os.system("rm "+FILTER+".list")
