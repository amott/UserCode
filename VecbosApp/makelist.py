#! /usr/bin/env python
import os
import sys
CASTORDIR = sys.argv[1]
FILTER = sys.argv[2]
os.system("rfdir "+CASTORDIR+" | grep "+FILTER+"> tmp.list")
list = open("tmp.list")
for line in list:
    sample = CASTORDIR+"/"+FILTER+line.split(FILTER)[-1][:-1]
    os.system("rfdir "+sample+" >> tmp2.list")
    list2 = open("tmp2.list")
    mylist = open(FILTER+line.split(FILTER)[-1][:-1]+".list","w")
    for line2 in list2:
        mylist.write((sample+"/de"+line2.split("de")[-1]).replace("//","/"))
        continue
    mylist.close()
    list2.close()
    os.system("rm tmp2.list")
    continue
list.close()
os.system("rm tmp.list")
