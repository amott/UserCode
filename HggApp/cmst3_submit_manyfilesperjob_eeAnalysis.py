#! /usr/bin/env python
import os, sys, commands, time
# set parameters to use cmst3 batch 
#######################################
### usage  cmst3_submit_manyfilesperjob.py dataset njobs applicationName queue 
#######################################
if len(sys.argv) != 9:
    print "usage cmst3_submit_manyfilesperjob.py process dataset njobs applicationName sample queue dirname isMC"
    sys.exit(1)
process = sys.argv[1]
dataset = sys.argv[2]
#settingfile = "config/RSZZsettings.txt"
output = dataset
# choose among cmt3 8nm 1nh 8nh 1nd 1nw 
#queue = "cmst3"
#queue = "cms8nht3"
queue = sys.argv[6]
#ijobmax = 40
ijobmax = int(sys.argv[3])
#application = "VecbosApp"
application = sys.argv[4]
sample = sys.argv[5]
dirname = sys.argv[7]
isMC = int(sys.argv[8])
if isMC != 0:
    inputlist = "cmst3_39X/MC/PU_2010/"+process+"/"+dataset+".list"
else:
    inputlist = "cmst3_39X/"+process+"/"+dataset+".list"
# to write on the cmst3 cluster disks
################################################
#castordir = "/castor/cern.ch/user/m/mpierini/CMST3/Vecbos/output/"
#outputmain = castordir+output
# to write on local disks
################################################
castordir = "/castor/cern.ch/user/e/emanuele/VecBos7TeV/"+dirname+"/"
diskoutputdir = "/cmsrm/pc22_2/micheli/data/VecBos3.9.X/"+dirname+"/"
outputmain = castordir+"/"+process+"/"+output+"/"+sample
diskoutputmain = diskoutputdir+"/"+process+"/"+output+"/"+sample
# prepare job to write on the cmst3 cluster disks
################################################
os.system("rm -rf "+dirname+"/"+process+"/"+output+"/"+"/"+sample)
os.system("mkdir -p "+dirname+"/"+process+"/"+output+"/"+"/"+sample)
os.system("mkdir -p "+dirname+"/"+process+"/"+output+"/"+"/"+sample+"/log/")
os.system("mkdir -p "+dirname+"/"+process+"/"+output+"/"+"/"+sample+"/input/")
os.system("mkdir -p "+dirname+"/"+process+"/"+output+"/"+"/"+sample+"/src/")
outputroot = outputmain+"/root/"
if castordir != "none": 
#    os.system("rfrm -r "+outputmain)
    os.system("rfmkdir -p "+castordir)
    os.system("rfmkdir -p "+outputmain)
    os.system("rfmkdir -p "+outputroot)
    os.system("rfchmod 777 "+castordir)
    os.system("rfchmod 777 "+outputmain)
    os.system("rfchmod 777 "+outputroot)
else: os.system("mkdir -p "+outputroot)

if diskoutputdir != "none": 
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm21 rm -rf "+diskoutputmain)
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm21 mkdir -p "+diskoutputdir)
    os.system("ssh -o BatchMode=yes -o StrictHostKeyChecking=no pccmsrm21 mkdir -p "+diskoutputmain)

#look for the current directory
#######################################
pwd = os.environ['PWD']
#######################################
#print inputlist
inputListfile=open(inputlist)
inputfiles = inputListfile.readlines()
ijob=0

#copy the configuration in the actual run directory
os.system("cp -r config "+dirname)

while (len(inputfiles) > 0):
    inputfilename = pwd+"/"+dirname+"/"+process+"/"+output+"/"+"/"+sample+"/input/input_"+str(ijob)+".list"
    inputfile = open(inputfilename,'w')
    for line in range(min(ijobmax,len(inputfiles))):
        ntpfile = inputfiles.pop()
        if ntpfile != '':
            inputfile.write(ntpfile)


    inputfile.close()

    # prepare the script to run
    outputname = dirname+"/"+process+"/"+output+"/"+"/"+sample+"/src/submit_"+str(ijob)+".src"
    outputfile = open(outputname,'w')
    outputfile.write('#!/bin/bash\n')
    outputfile.write('export STAGE_HOST=castorcms\n')
    outputfile.write('export STAGE_SVCCLASS=cmst3\n')
    #    outputfile.write('cd '+pwd)
    outputfile.write('cp '+pwd+'/data/Z_calibFall08.root $WORKDIR\n')
    outputfile.write('cp -r '+pwd+'/data/ $WORKDIR\n')
    outputfile.write('cp '+pwd+'/pdfs_MC.root $WORKDIR\n')
    outputfile.write('cp -r '+pwd+"/"+dirname+'/config $WORKDIR\n')
    outputfile.write('cd $WORKDIR\n')
    if(isMC) outputfile.write(pwd+'/'+application+' '+inputfilename+" "+output+"_"+str(ijob)+"_ "+" -signal="+sample+"\n")
    else outputfile.write(pwd+'/'+application+' '+inputfilename+" "+output+"_"+str(ijob)+"_ "+" -signal="+sample+" --isData"+"\n")
#    if castordir != "none": outputfile.write('./VecbosApp '+inputfilename+" "+" rfio://"+outputroot+output+"_"+str(ijob)+".root\n")
#    else:  
    outputfile.write('ls *.root | grep -v Z_calibFall08 | xargs -i rfcp {} '+outputroot+'\n')
    outputfile.write('ls *.root | grep -v Z_calibFall08 | xargs -i scp -o BatchMode=yes -o StrictHostKeyChecking=no {} pccmsrm21:'+diskoutputmain+'/{}\n') 
    outputfile.close
    os.system("echo bsub -q "+queue+" -o "+output+"/"+sample+"/log/"+output+"_"+str(ijob)+".log source "+pwd+"/"+outputname)
    os.system("bsub -q "+queue+" -o "+dirname+"/"+process+"/"+output+"/"+"/"+sample+"/log/"+output+"_"+str(ijob)+".log source "+pwd+"/"+outputname+" -copyInput="+process+"_"+str(ijob))
    time.sleep(3)
    ijob = ijob+1
    continue
