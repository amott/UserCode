#! /usr/bin/env python
import os
################################################
# set parameters to use cmst3 batch 
################################################
inputlist = "cmst3/Chowder_PDmuon.list"
settingfile = "settings.txt"
output = "pdmuon"
# choose among cmt3 8nm 1nh 8nh 1nd 1nw 
queue = "cmst3"
ijobmin = 0
ijobmax = 10
################################################
# to write on the cmst3 cluster disks
################################################
#castordir = "/castor/cern.ch/user/m/mpierini/CMST3/Vecbos/output/"
#outputmain = castordir+output
################################################
# to write on local disks
################################################
castordir = "none"
outputmain = output
################################################
# prepare job to write on the cmst3 cluster disks
################################################
os.system("mkdir "+output)
os.system("mkdir "+output+"/log/")
os.system("mkdir "+output+"/input/")
os.system("mkdir "+output+"/src/")
outputroot = outputmain+"/root/"
if castordir != "none": 
    os.system("rfmkdir "+outputmain)
    os.system("rfmkdir "+outputroot)
    os.system("rfchmod 777 "+outputmain)
    os.system("rfchmod 777 "+outputroot)
else: os.system("mkdir "+outputroot)
################################################

#################################################
ijob = 0
input = open(inputlist)
for ntpfile in input:
    if ijob >= ijobmin: 
        if ijob < ijobmax:
            # prepare the list file
            inputfilename = output+"/input/input_"+str(ijob)+".list"
            inputfile = open(inputfilename,'w')
            inputfile.write(ntpfile)
            inputfile.close()
            # prepare the script to run
            outputname = output+"/src/submit_"+str(ijob)+".src"
            outputfile = open(outputname,'w')
            outputfile.write('/bin/tcsh\n')
            outputfile.write('setenv STAGE_HOST castorcms\n')
            outputfile.write('setenv STAGE_SVCCLASS cmst3\n')
            outputfile.write('cd '+pwd)
            if castordir != "none": outputfile.write('./VecbosApp '+inputfilename+" "+settingfile+" rfio://"+outputroot+output+"_"+str(ijob)+".root \n")
            else:  outputfile.write('./VecbosApp '+inputfilename+" "+settingfile+" "+outputroot+output+"_"+str(ijob)+".root\n") 
            outputfile.close
            os.system("bsub -q "+queue+" -o "+pwd[:-1]+"/"+output+"/log/"+output+"_"+str(ijob)+".log source "+pwd[:-1]+"/"+outputname)
            ijob = ijob+1
            continue
        continue
    ijob = ijob+1
    continue
##################################################
