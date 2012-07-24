#! /usr/bin/env python
import os
import sys
reRecoJSON    = sys.argv[1]
promptRecoJSON = sys.argv[2]
pwd = os.environ['PWD']
lumifile = open("goodLumi.txt", "w")
for myFILE in os.popen("ls *json"):
    # ignore the Certification json files
    if myFILE.find("Cert_") != -1: continue
    if myFILE.find("json") == -1: continue
    myJSON = pwd+"/"+myFILE[:-1]
    #skip the file if not requested
    if myFILE.find("ReReco") != -1 and reRecoJSON == "none": continue
    if myFILE.find("ReReco") == -1 and promptRecoJSON == "none": continue
    # create the AND JSON
    if myFILE.find("ReReco") != -1: 
        os.system("cd /afs/cern.ch/user/m/mpierini/scratch0/CMSSW_4_2_0/src; eval `scramv1 run -sh`; compareJSON.py --and "+myJSON+" "+pwd+"/"+reRecoJSON+" "+myJSON+".good")
    else:
        os.system("cd /afs/cern.ch/user/m/mpierini/scratch0/CMSSW_4_2_0/src; eval `scramv1 run -sh`; compareJSON.py --and "+myJSON+" "+pwd+"/"+promptRecoJSON+" "+myJSON+".good")
    #run the LUMI
    os.system("cd /afs/cern.ch/user/m/mpierini/scratch0/CMSSW_4_2_0/src; eval `scramv1 run -sh`; lumiCalc.py -c frontier://LumiProd/CMS_LUMI_PROD -i "+myJSON+".good --nowarning overview > "+myJSON+".goodlumi")
    #get the lumi value
    mylumifile = open(myJSON+".goodlumi", 'r')
    lst = mylumifile.readlines()
    mylumi = float(lst[len(lst)-2].split("|")[4])
    mylumifile.close()
    lumifile.write(myJSON+"   "+str(int(mylumi/100000)/10.)+" pb-1\n")
    continue


os.system("rm *good")
os.system("rm *goodlumi")
