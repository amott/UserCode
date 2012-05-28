#! /usr/bin/env python
import os
import sys
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################
if len(sys.argv) < 7:
    print "usage cmst3_submit_manyfilesperjob.py process dataset njobs applicationName sample queue"
    print '      cmst3_submit_manyfilesperjob.py Summer11 TTJets_TuneZ2_7TeV-madgraph-tauola 50 VecbosApp TTJets 1nh'
    sys.exit(1)
process = sys.argv[1]
dataset = sys.argv[2]
isData = False
#isData = True
inputlist = "cmst3_42X/MC/"+process+"/"+dataset+".list"
if isData == True: inputlist = "cmst3_42X/"+process+"/"+dataset+".list"
output = dataset
# choose among cmt3 8nm 1nh 8nh 1nd 1nw 
#queue = "cmst3"
#queue = "cms8nht3"
queue = sys.argv[6]
ijobmax = int(sys.argv[3])
application = sys.argv[4]
sample = sys.argv[5]
# to write on the cmst3 cluster disks
#outputdir = "/castor/cern.ch/user/a/amott/CMST3/DiJet/"+dataset;
#outputdir = "/castor/cern.ch/user/a/amott/CMST3/RazorDiPhoton/"+dataset;
#outputdir = "/castor/cern.ch/user/a/amott/CMST3/RazorDiLepton/"+dataset;
#outputdir = "/castor/cern.ch/user/a/amott/VecbosApp/RazorHBB/"+dataset;
#os.system("rfmkdir "+outputdir)
#outputdirLXCMS132 = "/data/amott/VecbosApp/RazorBoostedTop/MC/"
#outputdirLXCMS132 = "/data1/amott/VecbosApp/RazorBoostedTop/MC/"
outputdirLXCMS132 = "/data1/amott/HiggsWW/VecbosApp/MC/"
################################################
os.system("mkdir -p "+process+"/"+output)
os.system("mkdir -p "+process+"/"+output+"/log/")
os.system("mkdir -p "+process+"/"+output+"/input/")
os.system("mkdir -p "+process+"/"+output+"/src/")
os.system("ssh -o StrictHostKeyChecking=no lxcms132 'mkdir "+outputdirLXCMS132+"'")
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
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export SCRAM_ARCH=slc5_amd64_gcc462')
    outputfile.write("cd /afs/cern.ch/user/a/amott/CMSSW_5_2_3/src; eval `scramv1 run -sh`\n")
    outputfile.write('cd '+pwd+'\n')
    if isData: outputfile.write('./VecbosApp '+inputfilename+" /tmp/"+output+"_"+str(ijob)+".root --isData\n")
    else: outputfile.write('./VecbosApp '+inputfilename+" /tmp/"+output+"_"+str(ijob)+".root \n")
    #outputfile.write("rfcp /tmp/"+output+"_"+str(ijob)+".root "+outputdir+"/\n") 
    outputfile.write("scp  -o StrictHostKeyChecking=no /tmp/"+output+"_"+str(ijob)+".root amott@lxcms132:"+outputdirLXCMS132+"/\n") 
    outputfile.write("rm /tmp/"+output+"_"+str(ijob)+".root \n")
    outputfile.close
    os.system("echo bsub -q "+queue+" -o /dev/null -e /dev/null source "+pwd+"/"+outputname)
    os.system("sleep 1; bsub -q "+queue+" -o /dev/null -e /dev/null source "+pwd+"/"+outputname)
    ijob = ijob+1
    continue
