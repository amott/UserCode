#! /usr/bin/env python
import os
import sys
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################
if len(sys.argv) < 4:
    print "usage cmst3_submit_manyfilesperjob.py card njobs applicationName queue"
    print '      cmst3_submit_manyfilesperjob.py cmst3_52X/Data/Collisions8TeV_V14_52X/METRun2012A.list 50 1nh'
    sys.exit(1)
process = sys.argv[1].split("/")[1]
dataset = sys.argv[1].split("/")[-1].replace(".list","")

isData = False#False Means MC

##########Data##############
#isData = True
#JSON = sys.argv[4]
############################

inputlist =  sys.argv[1]
output = dataset
# choose among cmt3 8nm 1nh 8nh 1nd 1nw 
#queue = "cmst3"
#queue = "cms8nht3"
queue = sys.argv[3]
ijobmax = int(sys.argv[2])
# to write on the cmst3 cluster disks
#outputdir = "/castor/cern.ch/user/m/mpierini/CMST3/DiJet/"+dataset;
#outputdir = "/castor/cern.ch/user/m/mpierini/CMST3/RazorDiPhoton/"+dataset;
#outputdir = "/castor/cern.ch/user/m/mpierini/CMST3/RazorDiLepton/"+dataset;
#outputdir = "/castor/cern.ch/user/m/mpierini/VecbosApp/RazorHBB/"+dataset;
#os.system("rfmkdir "+outputdir)
################################################
os.system("mkdir -p "+process+"/"+output)
os.system("mkdir -p "+process+"/"+output+"/log/")
os.system("mkdir -p "+process+"/"+output+"/input/")
os.system("mkdir -p "+process+"/"+output+"/src/")
os.system("mkdir -p "+process+"/"+output+"/out/")
#look for the current directory
#######################################
pwd = os.environ['PWD']
#######################################
numfiles = reduce(lambda x,y: x+1, file(inputlist).xreadlines(), 0)
filesperjob = numfiles/ijobmax
extrafiles  = numfiles%ijobmax
input = open(inputlist)
######################################

for ijob in range(ijobmax):
    # prepare the list file
    inputfilename = pwd+"/"+process+"/"+output+"/input/input_"+str(ijob)+".list"
    inputfile = open(inputfilename,'w')
    # if it is a normal job get filesperjob lines
    if ijob != (ijobmax-1):
        for line in range(filesperjob):
            ntpfile = input.readline() 
            inputfile.write(ntpfile)
            continue
    else:
        # if it is the last job get ALL remaining lines
        ntpfile = input.readline()
        while ntpfile != '':
            inputfile.write(ntpfile)
            ntpfile = input.readline()
            continue
    inputfile.close()

    # prepare the script to run
    outputname = process+"/"+output+"/src/submit_"+str(ijob)+".src"
    print "OUTPUTFILE: ", outputname
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write("cd /afs/cern.ch/work/c/cpena/scratch_DM/CMSSW_5_2_3/src; eval `scramv1 run -sh`\n")
    outputfile.write('cd '+pwd+'\n')
    if isData:
        outputfile.write('./VecbosApp '+inputfilename+" "+process+"/"+output+"/out/"+output+"_"+str(ijob)+".root --isData -json="+JSON+"\n")
        print './VecbosApp '+inputfilename+" "+process+"/"+output+"/out/"+output+"_"+str(ijob)+".root --isData -json="+JSON+"\n"
    else:
        outputfile.write('./VecbosApp '+inputfilename+" "+process+"/"+output+"/out/"+output+"_testMC_"+str(ijob)+".root \n")
        print './VecbosApp '+inputfilename+" "+process+"/"+output+"/out/"+output+"_testMc_"+str(ijob)+".root \n"
    outputfile.close
    os.system("echo bsub -q "+queue+" -o /dev/null -e /dev/null source "+pwd+"/"+outputname)
    os.system("sleep .1; bsub -q "+queue+" -o /dev/null -e /dev/null source "+pwd+"/"+outputname)
    #os.system("sleep 1; bsub -q "+queue+" source "+pwd+"/"+outputname)
    ijob = ijob+1
    continue
