#! /usr/bin/perl

use Getopt::Std;
getopts('q:w:f:n:C:d');

print "Starting...\n";

if(!$opt_w || !opt_f || $opt_h) { help(); }
if($opt_q){$queue = $opt_q;}
else{$queue="cmsslong";}
if($opt_w){$workdir=$opt_w;}
else{$workdir="/cmshome/dimarcoe/Vecbos/OfflineAnalysis/VecBosApp";}
if($opt_C){$CMSSWdir=$opt_C;}
else {$CMSSWdir="/cmshome/dimarcoe/Vecbos/CMSSW_RELEASES/CMSSW_2_0_11";}
if($opt_f){$file = $opt_f;}
else{$file="/cmshome/dimarcoe/Vecbos/OfflineAnalysis/VecBosApp/cmsrm-ui/ChowderPDelectron.list";}
if($opt_n){$njobs = $opt_n;}
else{$njobs = 1;}
if($opt_d){$donotsubmit=1;}
else{$donotsubmit=0;}

my $logdir = "$workdir/log";
my $rootdir = "$workdir/output";
my $scriptdir = "$workdir/batchscript";

# -- Create directories if not yet existent
if (-d "$logdir") {
    # do nothing
} else {
    system("/bin/mkdir $logdir"); 
    system("chmod 755 $logdir");
    if (-d "$logdir") {
        print " -> created $logdir\n";
    } else {
        die "run: cannot create $logdir\n";
    }
}
if (-d "$rootdir") {
    # do nothing
} else {
    system("/bin/mkdir $rootdir"); 
    system("chmod 755 $rootdir");
    if (-d "$rootdir") {
        print " -> created $rootdir\n";
    } else {
        die "run: cannot create $rootdir\n";
    }
}
if (-d "$scriptdir") {
    # clean the dir
    system("/bin/rm -f $scriptdir/*.list");
    system("/bin/rm -f $scriptdir/*.csh");
} else {
    system("/bin/mkdir $scriptdir"); 
    system("chmod 755 $scriptdir");
    if (-d "$scriptdir") {
        print " -> created $scriptdir\n";
    } else {
        die "run: cannot create $scriptdir\n";
    }
}


open(FULLINPUTLIST,"$file");
@datasets=<FULLINPUTLIST>;

$nfiles=$#datasets;
$stepsize=int($nfiles/$njobs)+1;
$job=1;

print "nfiles = $nfiles\n";
print "requested $njobs jobs\n";
print "number of datasets per job = $stepsize\n"; 

$firstline=0;
$lastline=0;
for($j=0; $j<$njobs; $j++) {
# create the split datasets
    $jobfile = "$scriptdir/dataset-$j.list";
    open(JOBFILE,">$jobfile");
    $firstline=$lastline;
    if(($stepsize)*($j+1)<$nfiles) {
	$lastline = ($stepsize)*($j+1)-1;
    } else {
	$lastline = $nfiles+1;
    }
    for($line=$firstline; $line<$lastline; $line++) {
	$lineToWrite = $datasets[$line];
	chop $lineToWrite;
	print JOBFILE "$lineToWrite\n";
    }

# create the scripts to submit
    $scriptfile = "$scriptdir/script-$j.csh";
    $prefix = "$rootdir/job-$j-";
    open(SCRIPTFILE,">$scriptfile");
    print SCRIPTFILE "\#\!/bin/tcsh\n\n";
    print SCRIPTFILE "cd $CMSSWdir\n";
    print SCRIPTFILE "eval `scramv1 runtime -csh`\n";
    print SCRIPTFILE "cd -\n";
    print SCRIPTFILE "$workdir/VecbosApp $jobfile $prefix\n";

#submit jobs 
    if(!$donotsubmit) {
	$logfile = "$logdir/vecbos-$j.log";
	$jobname = "VecbosApp-$j";
	system("bsub -q $queue -o $logfile -J $jobname < $scriptfile");
	print "submitting job n. $j to the queue $queue...\n";
    }
}

if(!$donotsubmit) {
    print "Submitted all the jobs. Enjoy the results of VecBosApp!\n";
} else {
    print "Scripts done, but you have given \"-d\" option. Jobs are not submitted.\n";
}


sub help() {
    print "Usage: ./submitBatch.pl  -w <workdir>  -f <filelist> [-q <queue>]  [-C <CMSSW_dir>]  [-n <numjobs>]  [-d] \n";
    print "Exemplum: ./submitBatch.pl -w /home/VecBosApp -f chowder.list -C /home/CMSSW_2_1_0 -n 10 \n";
    print "Options:\n";
    print "-w workdir:       the workdir where is the VecbosApp code\n";
    print "-f filelist:      the file with the list of all the ROOT vecbos ntuples \n";
    print "-q quesue:        the batch queue where to submit (default is cmslong) \n";
    print "-C CMSSW_dir:     the CMSSW directory where to do eval (to set up ROOT) \n";
    print "-n numjobs:       the number of jobs to submit \n";
    print "-d:               if given, create the scripts, the split lists, but do not submit the jobs \n";
    die "enjoy, vecbosapp is the spice of life...\n";
}
