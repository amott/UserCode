#!/usr/bin/env python

import os
import sys
import subprocess
import re

from optparse import OptionParser

ResubmitJobCodes = ['Cancelled','Aborted']
SafeResubmitErrorCodes = [60307,60317]

def getParser():
    parser = OptionParser()
    parser.add_option("-r","--register",  dest="register",default="", help="register crab directory for management")
    parser.add_option("-d","--deregister",dest="deregister",default="",help="de-register crab directory for management")
    parser.add_option("--deregisterAll",dest="deregisterAll",action="store_true",default="",help="de-register all crab directory for management")
    parser.add_option("-l","--list",dest="list",action="store_true",default=False,help="list the managed folders")
    parser.add_option("-c","--check",dest="check",action="store_true",default=False,help="check crab directories")
    parser.add_option("--resubmit",dest="resubmit",action="store_true",default=False,help="resubmit failed jobs with 'safe' error codes")
    parser.add_option("--resubmitAll",dest="resubmitAll",action="store_true",default=False,help="resubmit failed jobs with ALL error codes")
    parser.add_option("-v","--verbose",dest="verbose",action="store_true",default=False,help="run in verbose mode: list all output from crab -status")
    parser.add_option("--managerList",dest="managerList",default="~/.crabManager",help="specify the file listing the paths to manage [%default]")
    parser.add_option("-j","--job",dest="job",default=-1,type="int",help="Parser a specific job (specify the number from the -l command) [default=all]")
    parser.add_option("--command",dest="command",default="",help="issue a specific command (use -j to issue it to a specific job). specify crab commands without the leading '-'")
                      
    return parser

def addFolder(path,listFile):
    f = open(listFile,'a')
    f.write(path)
    f.write('\n')
    f.close()

def listFolders(listFile):
    f = open(listFile,'r')
    lines = f.readlines()
    f.close()
    return lines

def removeFolder(path,listFile):
    lines = listFolders(listFile)
    f=open(listFile,'w')
    for l in lines:
        if l != str(path)+'\n': f.write(l)
    f.close()
    
def crabManager():
    parser = getParser()
    (options,args)=parser.parse_args()

    managerList = os.path.abspath(os.path.expanduser(options.managerList))

    if not os.path.exists(options.managerList): os.system("touch %s" % options.managerList)

    if options.list:
        print "Managing: \n"
        for i,l in enumerate(listFolders(managerList)):
            print str(i)+" >> "+l.rstrip('\n')

    if options.register!="":
        addFolder(os.path.abspath(os.path.expanduser(options.register)),managerList)
        return

    if options.deregister!="":
        removeFolder(options.deregister,managerList)
        return

    if options.deregisterAll:
        for l in listFolders(managerList):
            print l
            removeFolder(l.rstrip('\n'),managerList)
        return


    if options.check or options.resubmit or options.resubmitAll or options.command!="":
        for i,l in enumerate(listFolders(managerList)):
            if options.job != -1 and options.job!=i: continue
            l=l.rstrip('\n')
            print "Processing Directory: "+str(l)
            [exitCodes,jobStatuses] = crab_check(l,options.verbose)
            if options.resubmit or options.resubmitAll:
                for status,jobs in jobStatuses.iteritems():
                    if status in ResubmitJobCodes:  # reprocess bad status codes
                        print "Status: "+status
                        crab_resubmit(l,jobs)
                for code,jobs in exitCodes.iteritems(): # resubmit jobs  that failed with error codes
                    if (int(code) in SafeResubmitErrorCodes) or (options.resubmitAll and int(code)!=0):
                        print "Code "+code
                        crab_resubmit(l,jobs)
            if options.command!="":
                crab_command(options.command,l)

def get_crab(option,folder,verbose):
    exe = ['crab',option,'-c',folder]
    return runProcess(exe,verbose)

def runProcess(exe,verbose):
    p = subprocess.Popen(exe, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    lines = []
    while(True):
        retcode = p.poll() #returns None while subprocess is running
        line = p.stdout.readline()
        lines.append(line)
        if verbose: print line.rstrip('\n')
        if(retcode is not None):
            break
    return lines

def crab_check(path,verbose):
    lines = get_crab('-status',path,verbose)

    ExitCode = re.compile(r"\s*>+\s*([0-9]+)[a-zA-Z\s]+:\s+([0-9]+)") # this regex matches the exit code line
    JobStatus = re.compile(r"\s*>+\s*([0-9]+)[a-zA-Z\s]+\s+([a-zA-Z]+)")# this regex matches the job status line

    JobList =re.compile("\s+[a-zA-Z\s]+:\s+([0-9,-]+)")
    exitCodes   = {}
    jobStatuses = {}
    isEX=-1
    isJS=""
    for l in lines:
        if isEX!=-1: # the list of jobs with a given exit code is always on the line after the exit code number
            m_codes = JobList.search(l)
            if m_codes is None:
                raise Exception("Invalid crab output format")
            exitCodes[isEX]=m_codes.group(1)
            isEX=-1
            continue

        if isJS!="":
            m_jobs = JobList.search(l)
            if m_jobs is None:
                continue ## the list of jobs isn't always on the line following the job status
            jobStatuses[isJS] = m_jobs.group(1)
            isJS=""
            continue

        ex = ExitCode.search(l)
        if ex is not None:
            print ex.group(0) # print the number of jobs with each exit code
            isEX=ex.group(2)  # record the exit code number for storage
            continue

        js = JobStatus.search(l)
        if js is not None:
            print js.group(0)
            isJS = js.group(2)
            continue

    return [exitCodes,jobStatuses]

def crab_command(command,path):
    cmd = "crab -%s -c %s" %(command,path)
    print cmd
    os.system(cmd)

def crab_resubmit(path,jobs):
    cmd = "crab -forceResubmit %s -c %s" %(jobs,path)
    print cmd
    os.system(cmd)

if __name__=='__main__':
    crabManager()
    
