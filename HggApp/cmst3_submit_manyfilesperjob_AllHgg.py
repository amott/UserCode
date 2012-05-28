#! /usr/bin/env python
import os
import sys
import math
import re
doSub = True
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################
if len(sys.argv) < 7:
    print "usage cmst3_submit_manyfilesperjob.py dir dataset filePerJob queue reducerCfg selectorCfg [json]"
    print '      cmst3_submit_manyfilesperjob.py Summer11 TTJets_TuneZ2_7TeV-madgraph-tauola 21 1nh cfg cfg [json]'
    sys.exit(1)
dir = sys.argv[1]
dataset = sys.argv[2]
isData = False
inputlist = dir+"/"+dataset+".list"
output = dataset
# choose among cmt3 8nm 1nh 8nh 1nd 1nw 
#queue = "cmst3"
#queue = "cms8nht3"
queue = sys.argv[4]
filesPerJob = int(sys.argv[3])
redCfg    = sys.argv[5]
selCfg    = sys.argv[6]
json=''
if len(sys.argv)>8:
    json = sys.argv[8]
# to write on the cmst3 cluster disks
#outputdir = "/castor/cern.ch/user/a/amott/CMST3/DiJet/"+dataset;
#outputdir = "/castor/cern.ch/user/a/amott/CMST3/RazorDiPhoton/"+dataset;
#outputdir = "/castor/cern.ch/user/a/amott/CMST3/RazorDiLepton/"+dataset;
#outputdir = "/castor/cern.ch/user/a/amott/VecbosApp/RazorHBB/"+dataset;
#os.system("rfmkdir "+outputdir)
#outputdirLXCMS132 = "/data/amott/VecbosApp/RazorBoostedTop/MC/"
#outputdirLXCMS132 = "/data1/amott/VecbosApp/RazorBoostedTop/Data/JetPD/"
outputdirRed = "/castor/cern.ch/user/a/amott/CMST3/Hgg/2012/Reduced/"+dataset
outputdirSel = "/castor/cern.ch/user/a/amott/CMST3/Hgg/2012/Selected/"+dataset
os.system("rfmkdir -p "+outputdirRed)
os.system("rfmkdir -p "+outputdirSel)
################################################
os.system("mkdir -p "+dir+"/"+output)
os.system("mkdir -p "+dir+"/"+output+"/log/")
os.system("mkdir -p "+dir+"/"+output+"/input/")
os.system("mkdir -p "+dir+"/"+output+"/src/")
#os.system("ssh -o StrictHostKeyChecking=no lxcms132 'mkdir "+outputdirLXCMS132+"'")
#look for the current directory
#######################################
pwd = os.environ['PWD']
#######################################
numfiles = reduce(lambda x,y: x+1, file(inputlist).xreadlines(), 0)
nJob = math.ceil(numfiles/filesPerJob)
#filesperjob = numfiles/ijobmax
#extrafiles  = numfiles%ijobmax
input = open(inputlist)
######################################

for ijob in range(nJob):
    # prepare the list file
    inputfilename = pwd+"/"+dir+"/"+output+"/input/input_"+str(ijob)+".list"
    inputfile = open(inputfilename,'w')
    # if it is a normal job get filesperjob lines
    getMass=True
    Mass = 0
    if ijob != (nJob-1):
        for line in range(filesPerJob):
            ntpfile = input.readline() 
            if getMass:
                reg = re.search("M-[0-9]+",ntpfile)
                Mass = float(reg.group(0)[2:])
                getMass=False
            inputfile.write(ntpfile)
            continue
    else:
        # if it is the last job get ALL remaining lines
        ntpfile = input.readline()
        while ntpfile != '':
            if getMass:
                reg = re.search("M-[0-9]+",ntpfile)
                Mass = float(reg.group(0)[2:])
                getMass=False
            inputfile.write(ntpfile)
            ntpfile = input.readline()
            continue
    inputfile.close()

    if Mass==0:
        print "WARNING: COULD NOT FIGURE OUT MASS FOR THIS JOB!!!  using "+str(ijob)
        Mass = ijob
    # prepare the script to run
    outputname = dir+"/"+output+"/src/submit_"+str(ijob)+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export SCRAM_ARCH=slc5_amd64_gcc462\n')
    outputfile.write("cd /afs/cern.ch/user/a/amott/CMSSW_5_2_3/src/; eval `scramv1 run -sh`\n")
    outputfile.write('cd '+pwd+'\n')
    if isData: outputfile.write('./HggApp '+inputfilename+" /tmp/"+output+"_"+str(ijob)+".root "+redCfg+" --isData -json="+str(json)+"\n")
    else: outputfile.write('./HggApp '+inputfilename+" /tmp/"+output+"_M"+str(Mass)+".root "+redCfg+" \n")
    outputfile.write("rfcp /tmp/"+output+"_M"+str(Mass)+".root "+outputdirRed+"/\n") 
    outputfile.write("echo /tmp/"+output+"_M"+str(Mass)+".root > /tmp/sel_M"+str(Mass)+".list\n")
    outputfile.write("./HggSelectorApp /tmp/sel_M"+str(Mass)+".list /tmp/"+output+"_Sel_M"+str(Mass)+".root "+selCfg+"\n")
    outputfile.write("rfcp /tmp/"+output+"_Sel_M"+str(Mass)+".root "+outputdirSel+"/\n")     
    #outputfile.write("scp  -o StrictHostKeyChecking=no /tmp/"+output+"_"+str(ijob)+".root amott@lxcms132:"+outputdirLXCMS132+"/\n") 
    outputfile.write("rm /tmp/"+output+"_M"+str(Mass)+".root \n")
    outputfile.write("rm /tmp/"+output+"_Sel_M"+str(Mass)+".root \n")
    outputfile.close
    os.system("echo bsub -q "+queue+" -o /dev/null -e /dev/null source "+pwd+"/"+outputname)
    if doSub: os.system("sleep 1; bsub -q "+queue+" -o /dev/null -e /dev/null source "+pwd+"/"+outputname)
#    ijob = ijob+1
#    continue
