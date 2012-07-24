#! /usr/bin/perl

use Getopt::Std;

$lumi = 100;
$chowder_lumi = 100;

$chowderFileList = "chowder_LogFiles.txt";
$ppEleX_pt14FileList = "ppEleX_pt14_LogFiles.txt";

print "creating tex table for lumi = $lumi pb-1....... \n";

# cross sections are taken from : https://twiki.cern.ch/twiki/bin/view/Main/AlpgenSummer07
# We use: Xsec x unweighting x matching efficiency
# ALPGEN SAMPLES
# $Xs_2jet_20_80 = 2.417E+8 * 0.11 * 0.34 ; 
# $Xs_2jet_80_140 = 9.951E+5 * 0.09 * 0.35 ;
# $Xs_2jet_140_5600 = 6.637E+4 * 0.09 * 0.24 ;

# $Xs_3jet_20_80 = 1.328E+7  * 0.13 * 0.25 ;
# $Xs_3jet_80_140 = 1.051E+6 * 0.15 * 0.23 ;
# $Xs_3jet_140_180 = 6.501E+4 * 0.24 * 0.17 ;
# $Xs_3jet_180_5600 = 2.861E+4 * 0.16 * 0.14 ;

# $Xs_4jet_20_100 = 1.835E+6 * 0.14 * 0.14 ;
# $Xs_4jet_100_160 = 1.625E+5 * 0.13 * 0.13 ;
# $Xs_4jet_160_200 = 1.912E+4 * 0.15 * 0.10 ;
# $Xs_4jet_200_250 = 7.283E+3 * 0.16 * 0.09 ;
# $Xs_4jet_250_400 = 3.498E+3 * 0.09 * 0.07 ;
# $Xs_4jet_400_5600 = 4.375E+2 * 0.14 * 0.07 ;

# $Xs_5jet_20_100 = 2.223E+5 * 0.13 * 0.10 ;
# $Xs_5jet_100_160 = 4.312E+4 * 0.10 * 0.08 ;
# $Xs_5jet_160_200 = 7.226E+3 * 0.11 * 0.07 ;
# $Xs_5jet_200_250 = 3.358E+3 * 0.10 * 0.065 ;
# $Xs_5jet_250_400 = 2.043E+3 * 0.09 * 0.06 ;
# $Xs_5jet_400_5600 = 2.371E+2 * 0.07 * 0.04 ;

# $Xs_6jet_20_100 = 3.196E+4 * 0.07 * 0.07 ; 
# $Xs_6jet_100_180 = 1.601E+4 * 0.11 * 0.08 ;
# $Xs_6jet_180_250 = 3.001E+3 * 0.14 * 0.07 ;
# $Xs_6jet_250_400 = 1.293E+3 * 0.07 * 0.07 ;
# $Xs_6jet_400_5600 = 2.315E+2 * 0.05 * 0.07 ;

# ppEleX Pythia
# https://cmsweb.cern.ch/dbs_discovery/getAppConfigs?dbsInst=cms_dbs_prod_global&appPath=*&procPath=/ppEleX_pt14/CMSSW_1_4_6-CSA07-3427/GEN-SIM&ajax=0&userMode=user
$Xs_ppEleX_pt14 = 55000000000 ; # pb
$filtereff_ppEleX_pt14 = 0.00004 ;

#init table
$textableW = "eff-W.tex";
open(TEXFILEW,">$textableW");
print TEXFILEW "\\begin\{sidewaystable\}\n";
print TEXFILEW "\\begin\{center\}\n";
print TEXFILEW "\\begin\{tabular\}\[t\]\{|l|c|c|c|c|\}\n";
print TEXFILEW "\\hline\n";
print TEXFILEW "\\hline\n";
print TEXFILEW "& & & & \\\\ \n";
print TEXFILEW "& \n";
print TEXFILEW "\$W+jets\$  \&\n";
print TEXFILEW "\$Z+jets\$\&\n";
print TEXFILEW "\$t \\bar\{t\}+jets\$ \&\n";
print TEXFILEW "\$QCD(pp \\rightarrow e\^{\\pm} X)\$\n";
print TEXFILEW "\\\\\n";
print TEXFILEW "\\hline\n";
print TEXFILEW "\\hline\n";

$textableZ = "eff-Z.tex";
open(TEXFILEZ,">$textableZ");
print TEXFILEZ "\\begin\{sidewaystable\}\n";
print TEXFILEZ "\\begin\{center\}\n";
print TEXFILEZ "\\begin\{tabular\}\[t\]\{|l|c|c|c|c|\}\n";
print TEXFILEZ "\\hline\n";
print TEXFILEZ "\\hline\n";
print TEXFILEZ "& & & & \\\\ \n";
print TEXFILEZ "& \n";
print TEXFILEZ "\$W+jets\$  \&\n";
print TEXFILEZ "\$Z+jets\$\&\n";
print TEXFILEZ "\$t \\bar\{t\}+jets\$ \&\n";
print TEXFILEZ "\$QCD(pp \\rightarrow e\^{\\pm} X)\$\n";
print TEXFILEZ "\\\\\n";
print TEXFILEZ "\\hline\n";
print TEXFILEZ "\\hline\n";





$Wj_event_Wsel = 0;
$Wj_trigger_Wsel = 0;
$Wj_nRecoEles_Wsel = 0;
$Wj_nAccEles_Wsel = 0;
$Wj_nIdEles_Wsel = 0;
$Wj_nTkIsolEles_Wsel = 0;
$Wj_nEcalIsolEles_Wsel = 0;
$Wj_nHcalIsolEles_Wsel = 0;
$Wj_nFinalEles_Wsel = 0;
$Wj_ZVeto_Wsel = 0;
$Wj_nPrimaryVertices_Wsel = 0;
$Wj_nDzVertex_Wsel = 0;
$Wj_nDxySigVertex_Wsel = 0;
$Wj_metCut_Wsel = 0;
$Wj_transvMassCut_Wsel = 0;

$Zj_event_Wsel = 0;
$Zj_trigger_Wsel = 0;
$Zj_nRecoEles_Wsel = 0;
$Zj_nAccEles_Wsel = 0;
$Zj_nIdEles_Wsel = 0;
$Zj_nTkIsolEles_Wsel = 0;
$Zj_nEcalIsolEles_Wsel = 0;
$Zj_nHcalIsolEles_Wsel = 0;
$Zj_nFinalEles_Wsel = 0;
$Zj_ZVeto_Wsel = 0;
$Zj_nPrimaryVertices_Wsel = 0;
$Zj_nDzVertex_Wsel = 0;
$Zj_nDxySigVertex_Wsel = 0;
$Zj_metCut_Wsel = 0;
$Zj_transvMassCut_Wsel = 0;

$ttbarj_event_Wsel = 0;
$ttbarj_trigger_Wsel = 0;
$ttbarj_nRecoEles_Wsel = 0;
$ttbarj_nAccEles_Wsel = 0;
$ttbarj_nIdEles_Wsel = 0;
$ttbarj_nTkIsolEles_Wsel = 0;
$ttbarj_nEcalIsolEles_Wsel = 0;
$ttbarj_nHcalIsolEles_Wsel = 0;
$ttbarj_nFinalEles_Wsel = 0;
$ttbarj_ZVeto_Wsel = 0;
$ttbarj_nPrimaryVertices_Wsel = 0;
$ttbarj_nDzVertex_Wsel = 0;
$ttbarj_nDxySigVertex_Wsel = 0;
$ttbarj_metCut_Wsel = 0;
$ttbarj_transvMassCut_Wsel = 0;

$Wj_event_Zsel = 0;
$Wj_trigger_Zsel = 0;
$Wj_nRecoEles_Zsel = 0;
$Wj_nAccEles_Zsel = 0;
$Wj_nIdEles_Zsel = 0;
$Wj_nTkIsolEles_Zsel = 0;
$Wj_nEcalIsolEles_Zsel = 0;
$Wj_nHcalIsolEles_Zsel = 0;
$Wj_nPrimaryVertices_Zsel = 0;
$Wj_nDzVertex_Zsel = 0;
$Wj_nDxySigVertex_Zsel = 0;
$Wj_meeCut_Zsel = 0;

$Zj_event_Zsel = 0;
$Zj_trigger_Zsel = 0;
$Zj_nRecoEles_Zsel = 0;
$Zj_nAccEles_Zsel = 0;
$Zj_nIdEles_Zsel = 0;
$Zj_nTkIsolEles_Zsel = 0;
$Zj_nEcalIsolEles_Zsel = 0;
$Zj_nHcalIsolEles_Zsel = 0;
$Zj_nPrimaryVertices_Zsel = 0;
$Zj_nDzVertex_Zsel = 0;
$Zj_nDxySigVertex_Zsel = 0;
$Zj_meeCut_Zsel = 0;

$ttbarj_event_Zsel = 0;
$ttbarj_trigger_Zsel = 0;
$ttbarj_nRecoEles_Zsel = 0;
$ttbarj_nAccEles_Zsel = 0;
$ttbarj_nIdEles_Zsel = 0;
$ttbarj_nTkIsolEles_Zsel = 0;
$ttbarj_nEcalIsolEles_Zsel = 0;
$ttbarj_nHcalIsolEles_Zsel = 0;
$ttbarj_nPrimaryVertices_Zsel = 0;
$ttbarj_nDzVertex_Zsel = 0;
$ttbarj_nDxySigVertex_Zsel = 0;
$ttbarj_meeCut_Zsel = 0;

$ppEleX_pt14_event_Wsel = 0;
$ppEleX_pt14_trigger_Wsel = 0;
$ppEleX_pt14_nRecoEles_Wsel = 0;
$ppEleX_pt14_nAccEles_Wsel = 0;
$ppEleX_pt14_nIdEles_Wsel = 0;
$ppEleX_pt14_nTkIsolEles_Wsel = 0;
$ppEleX_pt14_nEcalIsolEles_Wsel = 0;
$ppEleX_pt14_nHcalIsolEles_Wsel = 0;
$ppEleX_pt14_nFinalEles_Wsel = 0;
$ppEleX_pt14_ZVeto_Wsel = 0;
$ppEleX_pt14_nPrimaryVertices_Wsel = 0;
$ppEleX_pt14_nDzVertex_Wsel = 0;
$ppEleX_pt14_nDxySigVertex_Wsel = 0;
$ppEleX_pt14_metCut_Wsel = 0;
$ppEleX_pt14_transvMassCut_Wsel = 0;

$ppEleX_pt14_event_Zsel = 0;
$ppEleX_pt14_trigger_Zsel = 0;
$ppEleX_pt14_nRecoEles_Zsel = 0;
$ppEleX_pt14_nAccEles_Zsel = 0;
$ppEleX_pt14_nIdEles_Zsel = 0;
$ppEleX_pt14_nTkIsolEles_Zsel = 0;
$ppEleX_pt14_nEcalIsolEles_Zsel = 0;
$ppEleX_pt14_nHcalIsolEles_Zsel = 0;
$ppEleX_pt14_nPrimaryVertices_Zsel = 0;
$ppEleX_pt14_nDzVertex_Zsel = 0;
$ppEleX_pt14_nDxySigVertex_Zsel = 0;
$ppEleX_pt14_meeCut_Zsel = 0;



open(CHOWDERFILES,"$chowderFileList");
open(PPELEXFILES,"$ppEleX_pt14FileList");
@chowderfiles=<CHOWDERFILES>;
@ppelexfiles=<PPELEXFILES>;


foreach $Chowder_LogFile (@chowderfiles) {

    chop $Chowder_LogFile;
    print "now processing chowder file $Chowder_LogFile \n";

    $Wj_event_Wsel_tmp = 0;
    $Wj_trigger_Wsel_tmp = 0;
    $Wj_nRecoEles_Wsel_tmp = 0;
    $Wj_nAccEles_Wsel_tmp = 0;
    $Wj_nIdEles_Wsel_tmp = 0;
    $Wj_nTkIsolEles_Wsel_tmp = 0;
    $Wj_nEcalIsolEles_Wsel_tmp = 0;
    $Wj_nHcalIsolEles_Wsel_tmp = 0;
    $Wj_nFinalEles_Wsel_tmp = 0;
    $Wj_ZVeto_Wsel_tmp = 0;
    $Wj_nPrimaryVertices_Wsel_tmp = 0;
    $Wj_nDzVertex_Wsel_tmp = 0;
    $Wj_nDxySigVertex_Wsel_tmp = 0;
    $Wj_metCut_Wsel_tmp = 0;
    $Wj_transvMassCut_Wsel_tmp = 0;

    $Zj_event_Wsel_tmp = 0;
    $Zj_trigger_Wsel_tmp = 0;
    $Zj_nRecoEles_Wsel_tmp = 0;
    $Zj_nAccEles_Wsel_tmp = 0;
    $Zj_nIdEles_Wsel_tmp = 0;
    $Zj_nTkIsolEles_Wsel_tmp = 0;
    $Zj_nEcalIsolEles_Wsel_tmp = 0;
    $Zj_nHcalIsolEles_Wsel_tmp = 0;
    $Zj_nFinalEles_Wsel_tmp = 0;
    $Zj_ZVeto_Wsel_tmp = 0;
    $Zj_nPrimaryVertices_Wsel_tmp = 0;
    $Zj_nDzVertex_Wsel_tmp = 0;
    $Zj_nDxySigVertex_Wsel_tmp = 0;
    $Zj_metCut_Wsel_tmp = 0;
    $Zj_transvMassCut_Wsel_tmp = 0;

    $ttbarj_event_Wsel_tmp = 0;
    $ttbarj_trigger_Wsel_tmp = 0;
    $ttbarj_nRecoEles_Wsel_tmp = 0;
    $ttbarj_nAccEles_Wsel_tmp = 0;
    $ttbarj_nIdEles_Wsel_tmp = 0;
    $ttbarj_nTkIsolEles_Wsel_tmp = 0;
    $ttbarj_nEcalIsolEles_Wsel_tmp = 0;
    $ttbarj_nHcalIsolEles_Wsel_tmp = 0;
    $ttbarj_nFinalEles_Wsel_tmp = 0;
    $ttbarj_ZVeto_Wsel_tmp = 0;
    $ttbarj_nPrimaryVertices_Wsel_tmp = 0;
    $ttbarj_nDzVertex_Wsel_tmp = 0;
    $ttbarj_nDxySigVertex_Wsel_tmp = 0;
    $ttbarj_metCut_Wsel_tmp = 0;
    $ttbarj_transvMassCut_Wsel_tmp = 0;

    $Wj_event_Zsel_tmp = 0;
    $Wj_trigger_Zsel_tmp = 0;
    $Wj_nRecoEles_Zsel_tmp = 0;
    $Wj_nAccEles_Zsel_tmp = 0;
    $Wj_nIdEles_Zsel_tmp = 0;
    $Wj_nTkIsolEles_Zsel_tmp = 0;
    $Wj_nEcalIsolEles_Zsel_tmp = 0;
    $Wj_nHcalIsolEles_Zsel_tmp = 0;
    $Wj_nPrimaryVertices_Zsel_tmp = 0;
    $Wj_nDzVertex_Zsel_tmp = 0;
    $Wj_nDxySigVertex_Zsel_tmp = 0;
    $Wj_meeCut_Zsel_tmp = 0;

    $Zj_event_Zsel_tmp = 0;
    $Zj_trigger_Zsel_tmp = 0;
    $Zj_nRecoEles_Zsel_tmp = 0;
    $Zj_nAccEles_Zsel_tmp = 0;
    $Zj_nIdEles_Zsel_tmp = 0;
    $Zj_nTkIsolEles_Zsel_tmp = 0;
    $Zj_nEcalIsolEles_Zsel_tmp = 0;
    $Zj_nHcalIsolEles_Zsel_tmp = 0;
    $Zj_nPrimaryVertices_Zsel_tmp = 0;
    $Zj_nDzVertex_Zsel_tmp = 0;
    $Zj_nDxySigVertex_Zsel_tmp = 0;
    $Zj_meeCut_Zsel_tmp = 0;

    $ttbarj_event_Zsel_tmp = 0;
    $ttbarj_trigger_Zsel_tmp = 0;
    $ttbarj_nRecoEles_Zsel_tmp = 0;
    $ttbarj_nAccEles_Zsel_tmp = 0;
    $ttbarj_nIdEles_Zsel_tmp = 0;
    $ttbarj_nTkIsolEles_Zsel_tmp = 0;
    $ttbarj_nEcalIsolEles_Zsel_tmp = 0;
    $ttbarj_nHcalIsolEles_Zsel_tmp = 0;
    $ttbarj_nPrimaryVertices_Zsel_tmp = 0;
    $ttbarj_nDzVertex_Zsel_tmp = 0;
    $ttbarj_nDxySigVertex_Zsel_tmp = 0;
    $ttbarj_meeCut_Zsel_tmp = 0;


    open(CHOWDERLOGFILE,$Chowder_LogFile) or die "cannot open $Chowder_LogFile !\n";
    seek (CHOWDERLOGFILE, -30000, 2);
    $irow=1;
    while($row=<CHOWDERLOGFILE>) {
	if($row=~/Full\sZee\sselections\:\sZ\+jets\ssample\s\-\-\>/) {
	    last;
	}
	chop $row;

	if($row=~/\*\s+event\:\s+(\S+)$/) {
	    {$Wj_event_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+trigger\:\s+(\S+)$/) {
	    {$Wj_trigger_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nRecoEles\:\s+(\S+)$/) {
	    {$Wj_nRecoEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nAccEles\:\s+(\S+)$/) {
	    {$Wj_nAccEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nIdEles\:\s+(\S+)$/) {
	    {$Wj_nIdEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nTkIsolEles\:\s+(\S+)$/) {
	    {$Wj_nTkIsolEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nEcalIsolEles\:\s+(\S+)$/) {
	    {$Wj_nEcalIsolEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nHcalIsolEles\:\s+(\S+)$/) {
	    {$Wj_nHcalIsolEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nFinalEles\:\s+(\S+)$/) {
	    {$Wj_nFinalEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+ZVeto\:\s+(\S+)$/) {
	    {$Wj_ZVeto_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nPrimaryVertices\:\s+(\S+)$/) {
	    {$Wj_nPrimaryVertices_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nDzVertex\:\s+(\S+)$/) {
	    {$Wj_nDzVertex_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nDxySigVertex\:\s+(\S+)$/) {
	    {$Wj_nDxySigVertex_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+metCut\:\s+(\S+)$/) {
	    {$Wj_metCut_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+transvMassCut\:\s+(\S+)$/) {
	    {$Wj_transvMassCut_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

    }


    open(CHOWDERLOGFILE,$Chowder_LogFile) or die "cannot open $Chowder_LogFile !\n";
    seek (CHOWDERLOGFILE, -30000, 2);
    $irow=1;
    while($row=<CHOWDERLOGFILE>) {
	if($row=~/Full\sZee\sselections\:\sttbar\ssample\s\-\-\>/) {
	    last;
	}
	chop $row;

	if($row=~/\*\s+event\:\s+(\S+)$/) {
	    {$Zj_event_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+trigger\:\s+(\S+)$/) {
	    {$Zj_trigger_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nRecoEles\:\s+(\S+)$/) {
	    {$Zj_nRecoEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nAccEles\:\s+(\S+)$/) {
	    {$Zj_nAccEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nIdEles\:\s+(\S+)$/) {
	    {$Zj_nIdEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nTkIsolEles\:\s+(\S+)$/) {
	    {$Zj_nTkIsolEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nEcalIsolEles\:\s+(\S+)$/) {
	    {$Zj_nEcalIsolEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nHcalIsolEles\:\s+(\S+)$/) {
	    {$Zj_nHcalIsolEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nFinalEles\:\s+(\S+)$/) {
	    {$Zj_nFinalEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+ZVeto\:\s+(\S+)$/) {
	    {$Zj_ZVeto_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nPrimaryVertices\:\s+(\S+)$/) {
	    {$Zj_nPrimaryVertices_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nDzVertex\:\s+(\S+)$/) {
	    {$Zj_nDzVertex_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nDxySigVertex\:\s+(\S+)$/) {
	    {$Zj_nDxySigVertex_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+metCut\:\s+(\S+)$/) {
	    {$Zj_metCut_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+transvMassCut\:\s+(\S+)$/) {
	    {$Zj_transvMassCut_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

    }


    open(CHOWDERLOGFILE,$Chowder_LogFile) or die "cannot open $Chowder_LogFile !\n";
    seek (CHOWDERLOGFILE, -30000, 2);
    $irow=1;
    while($row=<CHOWDERLOGFILE>) {
	if($row=~/ElectronID\sselections\:/) {
	    last;
	}
	chop $row;

	if($row=~/\*\s+event\:\s+(\S+)$/) {
	    {$ttbarj_event_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+trigger\:\s+(\S+)$/) {
	    {$ttbarj_trigger_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nRecoEles\:\s+(\S+)$/) {
	    {$ttbarj_nRecoEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nAccEles\:\s+(\S+)$/) {
	    {$ttbarj_nAccEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nIdEles\:\s+(\S+)$/) {
	    {$ttbarj_nIdEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nTkIsolEles\:\s+(\S+)$/) {
	    {$ttbarj_nTkIsolEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nEcalIsolEles\:\s+(\S+)$/) {
	    {$ttbarj_nEcalIsolEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nHcalIsolEles\:\s+(\S+)$/) {
	    {$ttbarj_nHcalIsolEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nFinalEles\:\s+(\S+)$/) {
	    {$ttbarj_nFinalEles_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+ZVeto\:\s+(\S+)$/) {
	    {$ttbarj_ZVeto_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nPrimaryVertices\:\s+(\S+)$/) {
	    {$ttbarj_nPrimaryVertices_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nDzVertex\:\s+(\S+)$/) {
	    {$ttbarj_nDzVertex_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nDxySigVertex\:\s+(\S+)$/) {
	    {$ttbarj_nDxySigVertex_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+metCut\:\s+(\S+)$/) {
	    {$ttbarj_metCut_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+transvMassCut\:\s+(\S+)$/) {
	    {$ttbarj_transvMassCut_Wsel_tmp=$1*$lumi/$chowder_lumi;}
	}

    }

    open(CHOWDERLOGFILE,$Chowder_LogFile) or die "cannot open $Chowder_LogFile !\n";
    seek (CHOWDERLOGFILE, -30000, 2);
    $irow=1;
    while($row=<CHOWDERLOGFILE>) {
	if($row=~/Full\sWenu\sselections\:\sW\+jets\ssample\s\-\-\>/) {
	    last;
	}
	chop $row;

	if($row=~/\*\s+event\:\s+(\S+)$/) {
	    {$Wj_event_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+trigger\:\s+(\S+)$/) {
	    {$Wj_trigger_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nRecoEles\:\s+(\S+)$/) {
	    {$Wj_nRecoEles_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nAccEles\:\s+(\S+)$/) {
	    {$Wj_nAccEles_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nIdEles\:\s+(\S+)$/) {
	    {$Wj_nIdEles_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nTkIsolEles\:\s+(\S+)$/) {
	    {$Wj_nTkIsolEles_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nEcalIsolEles\:\s+(\S+)$/) {
	    {$Wj_nEcalIsolEles_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nHcalIsolEles\:\s+(\S+)$/) {
	    {$Wj_nHcalIsolEles_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nPrimaryVertices\:\s+(\S+)$/) {
	    {$Wj_nPrimaryVertices_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nDzVertex\:\s+(\S+)$/) {
	    {$Wj_nDzVertex_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nDxySigVertex\:\s+(\S+)$/) {
	    {$Wj_nDxySigVertex_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+meeCut\:\s+(\S+)$/) {
	    {$Wj_meeCut_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

    }


    open(CHOWDERLOGFILE,$Chowder_LogFile) or die "cannot open $Chowder_LogFile !\n";
    seek (CHOWDERLOGFILE, -30000, 2);
    $irow=1;
    while($row=<CHOWDERLOGFILE>) {
	if($row=~/Full\sWenu\sselections\:\sZ\+jets\ssample\s\-\-\>/) {
	    last;
	}
	chop $row;

	if($row=~/\*\s+event\:\s+(\S+)$/) {
	    {$Zj_event_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+trigger\:\s+(\S+)$/) {
	    {$Zj_trigger_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nRecoEles\:\s+(\S+)$/) {
	    {$Zj_nRecoEles_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nAccEles\:\s+(\S+)$/) {
	    {$Zj_nAccEles_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nIdEles\:\s+(\S+)$/) {
	    {$Zj_nIdEles_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nTkIsolEles\:\s+(\S+)$/) {
	    {$Zj_nTkIsolEles_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nEcalIsolEles\:\s+(\S+)$/) {
	    {$Zj_nEcalIsolEles_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nHcalIsolEles\:\s+(\S+)$/) {
	    {$Zj_nHcalIsolEles_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nPrimaryVertices\:\s+(\S+)$/) {
	    {$Zj_nPrimaryVertices_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nDzVertex\:\s+(\S+)$/) {
	    {$Zj_nDzVertex_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nDxySigVertex\:\s+(\S+)$/) {
	    {$Zj_nDxySigVertex_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+meeCut\:\s+(\S+)$/) {
	    {$Zj_meeCut_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

    }


    open(CHOWDERLOGFILE,$Chowder_LogFile) or die "cannot open $Chowder_LogFile !\n";
    seek (CHOWDERLOGFILE, -30000, 2);
    $irow=1;
    while($row=<CHOWDERLOGFILE>) {
	if($row=~/Full\sWenu\sselections\:\sttbar\ssample\s\-\-\>/) {
	    last;
	}
	chop $row;

	if($row=~/\*\s+event\:\s+(\S+)$/) {
	    {$ttbarj_event_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+trigger\:\s+(\S+)$/) {
	    {$ttbarj_trigger_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nRecoEles\:\s+(\S+)$/) {
	    {$ttbarj_nRecoEles_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nAccEles\:\s+(\S+)$/) {
	    {$ttbarj_nAccEles_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nIdEles\:\s+(\S+)$/) {
	    {$ttbarj_nIdEles_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nTkIsolEles\:\s+(\S+)$/) {
	    {$ttbarj_nTkIsolEles_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nEcalIsolEles\:\s+(\S+)$/) {
	    {$ttbarj_nEcalIsolEles_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nHcalIsolEles\:\s+(\S+)$/) {
	    {$ttbarj_nHcalIsolEles_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nPrimaryVertices\:\s+(\S+)$/) {
	    {$ttbarj_nPrimaryVertices_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nDzVertex\:\s+(\S+)$/) {
	    {$ttbarj_nDzVertex_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+nDxySigVertex\:\s+(\S+)$/) {
	    {$ttbarj_nDxySigVertex_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

	if($row=~/\*\s+meeCut\:\s+(\S+)$/) {
	    {$ttbarj_meeCut_Zsel_tmp=$1*$lumi/$chowder_lumi;}
	}

    }


    $Wj_event_Wsel += $Wj_event_Wsel_tmp;
    $Wj_trigger_Wsel += $Wj_trigger_Wsel_tmp;
    $Wj_nRecoEles_Wsel += $Wj_nRecoEles_Wsel_tmp;
    $Wj_nAccEles_Wsel += $Wj_nAccEles_Wsel_tmp;
    $Wj_nIdEles_Wsel += $Wj_nIdEles_Wsel_tmp;
    $Wj_nTkIsolEles_Wsel += $Wj_nTkIsolEles_Wsel_tmp;
    $Wj_nEcalIsolEles_Wsel += $Wj_nEcalIsolEles_Wsel_tmp;
    $Wj_nHcalIsolEles_Wsel += $Wj_nHcalIsolEles_Wsel_tmp;
    $Wj_nFinalEles_Wsel += $Wj_nFinalEles_Wsel_tmp;
    $Wj_ZVeto_Wsel += $Wj_ZVeto_Wsel_tmp;
    $Wj_nPrimaryVertices_Wsel += $Wj_nPrimaryVertices_Wsel_tmp;
    $Wj_nDzVertex_Wsel += $Wj_nDzVertex_Wsel_tmp;
    $Wj_nDxySigVertex_Wsel += $Wj_nDxySigVertex_Wsel_tmp;
    $Wj_metCut_Wsel += $Wj_metCut_Wsel_tmp;
    $Wj_transvMassCut_Wsel += $Wj_transvMassCut_Wsel_tmp;
    
    $Zj_event_Wsel += $Zj_event_Wsel_tmp;
    $Zj_trigger_Wsel += $Zj_trigger_Wsel_tmp;
    $Zj_nRecoEles_Wsel += $Zj_nRecoEles_Wsel_tmp;
    $Zj_nAccEles_Wsel += $Zj_nAccEles_Wsel_tmp;
    $Zj_nIdEles_Wsel += $Zj_nIdEles_Wsel_tmp;
    $Zj_nTkIsolEles_Wsel += $Zj_nTkIsolEles_Wsel_tmp;
    $Zj_nEcalIsolEles_Wsel += $Zj_nEcalIsolEles_Wsel_tmp;
    $Zj_nHcalIsolEles_Wsel += $Zj_nHcalIsolEles_Wsel_tmp;
    $Zj_nFinalEles_Wsel += $Zj_nFinalEles_Wsel_tmp;
    $Zj_ZVeto_Wsel += $Zj_ZVeto_Wsel_tmp;
    $Zj_nPrimaryVertices_Wsel += $Zj_nPrimaryVertices_Wsel_tmp;
    $Zj_nDzVertex_Wsel += $Zj_nDzVertex_Wsel_tmp;
    $Zj_nDxySigVertex_Wsel += $Zj_nDxySigVertex_Wsel_tmp;
    $Zj_metCut_Wsel += $Zj_metCut_Wsel_tmp;
    $Zj_transvMassCut_Wsel += $Zj_transvMassCut_Wsel_tmp;

    $ttbarj_event_Wsel += $ttbarj_event_Wsel_tmp;
    $ttbarj_trigger_Wsel += $ttbarj_trigger_Wsel_tmp;
    $ttbarj_nRecoEles_Wsel += $ttbarj_nRecoEles_Wsel_tmp;
    $ttbarj_nAccEles_Wsel += $ttbarj_nAccEles_Wsel_tmp;
    $ttbarj_nIdEles_Wsel += $ttbarj_nIdEles_Wsel_tmp;
    $ttbarj_nTkIsolEles_Wsel += $ttbarj_nTkIsolEles_Wsel_tmp;
    $ttbarj_nEcalIsolEles_Wsel += $ttbarj_nEcalIsolEles_Wsel_tmp;
    $ttbarj_nHcalIsolEles_Wsel += $ttbarj_nHcalIsolEles_Wsel_tmp;
    $ttbarj_nFinalEles_Wsel += $ttbarj_nFinalEles_Wsel_tmp;
    $ttbarj_ZVeto_Wsel += $ttbarj_ZVeto_Wsel_tmp;
    $ttbarj_nPrimaryVertices_Wsel += $ttbarj_nPrimaryVertices_Wsel_tmp;
    $ttbarj_nDzVertex_Wsel += $ttbarj_nDzVertex_Wsel_tmp;
    $ttbarj_nDxySigVertex_Wsel += $ttbarj_nDxySigVertex_Wsel_tmp;
    $ttbarj_metCut_Wsel += $ttbarj_metCut_Wsel_tmp;
    $ttbarj_transvMassCut_Wsel += $ttbarj_transvMassCut_Wsel_tmp;

    $Wj_event_Zsel += $Wj_event_Zsel_tmp;
    $Wj_trigger_Zsel += $Wj_trigger_Zsel_tmp;
    $Wj_nRecoEles_Zsel += $Wj_nRecoEles_Zsel_tmp;
    $Wj_nAccEles_Zsel += $Wj_nAccEles_Zsel_tmp;
    $Wj_nIdEles_Zsel += $Wj_nIdEles_Zsel_tmp;
    $Wj_nTkIsolEles_Zsel += $Wj_nTkIsolEles_Zsel_tmp;
    $Wj_nEcalIsolEles_Zsel += $Wj_nEcalIsolEles_Zsel_tmp;
    $Wj_nHcalIsolEles_Zsel += $Wj_nHcalIsolEles_Zsel_tmp;
    $Wj_nPrimaryVertices_Zsel += $Wj_nPrimaryVertices_Zsel_tmp;
    $Wj_nDzVertex_Zsel += $Wj_nDzVertex_Zsel_tmp;
    $Wj_nDxySigVertex_Zsel += $Wj_nDxySigVertex_Zsel_tmp;
    $Wj_meeCut_Zsel += $Wj_meeCut_Zsel_tmp;

    $Zj_event_Zsel += $Zj_event_Zsel_tmp;
    $Zj_trigger_Zsel += $Zj_trigger_Zsel_tmp;
    $Zj_nRecoEles_Zsel += $Zj_nRecoEles_Zsel_tmp;
    $Zj_nAccEles_Zsel += $Zj_nAccEles_Zsel_tmp;
    $Zj_nIdEles_Zsel += $Zj_nIdEles_Zsel_tmp;
    $Zj_nTkIsolEles_Zsel += $Zj_nTkIsolEles_Zsel_tmp;
    $Zj_nEcalIsolEles_Zsel += $Zj_nEcalIsolEles_Zsel_tmp;
    $Zj_nHcalIsolEles_Zsel += $Zj_nHcalIsolEles_Zsel_tmp;
    $Zj_nPrimaryVertices_Zsel += $Zj_nPrimaryVertices_Zsel_tmp;
    $Zj_nDzVertex_Zsel += $Zj_nDzVertex_Zsel_tmp;
    $Zj_nDxySigVertex_Zsel += $Zj_nDxySigVertex_Zsel_tmp;
    $Zj_meeCut_Zsel += $Zj_meeCut_Zsel_tmp;

    $ttbarj_event_Zsel += $ttbarj_event_Zsel_tmp;
    $ttbarj_trigger_Zsel += $ttbarj_trigger_Zsel_tmp;
    $ttbarj_nRecoEles_Zsel += $ttbarj_nRecoEles_Zsel_tmp;
    $ttbarj_nAccEles_Zsel += $ttbarj_nAccEles_Zsel_tmp;
    $ttbarj_nIdEles_Zsel += $ttbarj_nIdEles_Zsel_tmp;
    $ttbarj_nTkIsolEles_Zsel += $ttbarj_nTkIsolEles_Zsel_tmp;
    $ttbarj_nEcalIsolEles_Zsel += $ttbarj_nEcalIsolEles_Zsel_tmp;
    $ttbarj_nHcalIsolEles_Zsel += $ttbarj_nHcalIsolEles_Zsel_tmp;
    $ttbarj_nPrimaryVertices_Zsel += $ttbarj_nPrimaryVertices_Zsel_tmp;
    $ttbarj_nDzVertex_Zsel += $ttbarj_nDzVertex_Zsel_tmp;
    $ttbarj_nDxySigVertex_Zsel += $ttbarj_nDxySigVertex_Zsel_tmp;
    $ttbarj_meeCut_Zsel += $ttbarj_meeCut_Zsel_tmp;

}




foreach $ppEleX_pt14_LogFile (@ppelexfiles) {

    $ppEleX_pt14_event_Wsel_tmp = 0;
    $ppEleX_pt14_trigger_Wsel_tmp = 0;
    $ppEleX_pt14_nRecoEles_Wsel_tmp = 0;
    $ppEleX_pt14_nAccEles_Wsel_tmp = 0;
    $ppEleX_pt14_nIdEles_Wsel_tmp = 0;
    $ppEleX_pt14_nTkIsolEles_Wsel_tmp = 0;
    $ppEleX_pt14_nEcalIsolEles_Wsel_tmp = 0;
    $ppEleX_pt14_nHcalIsolEles_Wsel_tmp = 0;
    $ppEleX_pt14_nFinalEles_Wsel_tmp = 0;
    $ppEleX_pt14_ZVeto_Wsel_tmp = 0;
    $ppEleX_pt14_nPrimaryVertices_Wsel_tmp = 0;
    $ppEleX_pt14_nDzVertex_Wsel_tmp = 0;
    $ppEleX_pt14_nDxySigVertex_Wsel_tmp = 0;
    $ppEleX_pt14_metCut_Wsel_tmp = 0;
    $ppEleX_pt14_transvMassCut_Wsel_tmp = 0;

    $ppEleX_pt14_event_Zsel_tmp = 0;
    $ppEleX_pt14_trigger_Zsel_tmp = 0;
    $ppEleX_pt14_nRecoEles_Zsel_tmp = 0;
    $ppEleX_pt14_nAccEles_Zsel_tmp = 0;
    $ppEleX_pt14_nIdEles_Zsel_tmp = 0;
    $ppEleX_pt14_nTkIsolEles_Zsel_tmp = 0;
    $ppEleX_pt14_nEcalIsolEles_Zsel_tmp = 0;
    $ppEleX_pt14_nHcalIsolEles_Zsel_tmp = 0;
    $ppEleX_pt14_nPrimaryVertices_Zsel_tmp = 0;
    $ppEleX_pt14_nDzVertex_Zsel_tmp = 0;
    $ppEleX_pt14_nDxySigVertex_Zsel_tmp = 0;
    $ppEleX_pt14_meeCut_Zsel_tmp = 0;

    chop $ppEleX_pt14_LogFile;
    print "now processing ppEleX_pt14 file $ppEleX_pt14_LogFile \n";

    open(PPELEX_PT14LOGFILE,$ppEleX_pt14_LogFile) or die "cannot open $ppEleX_pt14_LogFile !\n";
    seek (PPELEX_PT14LOGFILE, -30000, 2);
    $irow=1;
    while($row=<PPELEX_PT14LOGFILE>) {
	if($row=~/ElectronID\sselections\:/) {
	    last;
	}
	chop $row;

	if($row=~/\*\s+event\:\s+(\S+)$/) {
	    {$ppEleX_pt14_event_Wsel_tmp=$1;}
	}

	if($row=~/\*\s+trigger\:\s+(\S+)$/) {
	    {$ppEleX_pt14_trigger_Wsel_tmp=$1;}
	}

	if($row=~/\*\s+nRecoEles\:\s+(\S+)$/) {
	    {$ppEleX_pt14_nRecoEles_Wsel_tmp=$1;}
	}

	if($row=~/\*\s+nAccEles\:\s+(\S+)$/) {
	    {$ppEleX_pt14_nAccEles_Wsel_tmp=$1;}
	}

	if($row=~/\*\s+nIdEles\:\s+(\S+)$/) {
	    {$ppEleX_pt14_nIdEles_Wsel_tmp=$1;}
	}

	if($row=~/\*\s+nTkIsolEles\:\s+(\S+)$/) {
	    {$ppEleX_pt14_nTkIsolEles_Wsel_tmp=$1;}
	}

	if($row=~/\*\s+nEcalIsolEles\:\s+(\S+)$/) {
	    {$ppEleX_pt14_nEcalIsolEles_Wsel_tmp=$1;}
	}

	if($row=~/\*\s+nHcalIsolEles\:\s+(\S+)$/) {
	    {$ppEleX_pt14_nHcalIsolEles_Wsel_tmp=$1;}
	}

	if($row=~/\*\s+nFinalEles\:\s+(\S+)$/) {
	    {$ppEleX_pt14_nFinalEles_Wsel_tmp=$1;}
	}

	if($row=~/\*\s+ZVeto\:\s+(\S+)$/) {
	    {$ppEleX_pt14_ZVeto_Wsel_tmp=$1;}
	}

	if($row=~/\*\s+nPrimaryVertices\:\s+(\S+)$/) {
	    {$ppEleX_pt14_nPrimaryVertices_Wsel_tmp=$1;}
	}

	if($row=~/\*\s+nDzVertex\:\s+(\S+)$/) {
	    {$ppEleX_pt14_nDzVertex_Wsel_tmp=$1;}
	}

	if($row=~/\*\s+nDxySigVertex\:\s+(\S+)$/) {
	    {$ppEleX_pt14_nDxySigVertex_Wsel_tmp=$1;}
	}

	if($row=~/\*\s+metCut\:\s+(\S+)$/) {
	    {$ppEleX_pt14_metCut_Wsel_tmp=$1;}
	}

	if($row=~/\*\s+transvMassCut\:\s+(\S+)$/) {
	    {$ppEleX_pt14_transvMassCut_Wsel_tmp=$1;}
	}

    }

    open(PPELEX_PT14LOGFILE,$ppEleX_pt14_LogFile) or die "cannot open $ppEleX_pt14_LogFile !\n";
    seek (PPELEX_PT14LOGFILE, -30000, 2);
    $irow=1;
    while($row=<PPELEX_PT14LOGFILE>) {
	if($row=~/Full\sWenu\sselections\:\s\-\-\>/) {
	    last;
	}
	chop $row;

	if($row=~/\*\s+event\:\s+(\S+)$/) {
	    {$ppEleX_pt14_event_Zsel_tmp=$1;}
	}

	if($row=~/\*\s+trigger\:\s+(\S+)$/) {
	    {$ppEleX_pt14_trigger_Zsel_tmp=$1;}
	}

	if($row=~/\*\s+nRecoEles\:\s+(\S+)$/) {
	    {$ppEleX_pt14_nRecoEles_Zsel_tmp=$1;}
	}

	if($row=~/\*\s+nAccEles\:\s+(\S+)$/) {
	    {$ppEleX_pt14_nAccEles_Zsel_tmp=$1;}
	}

	if($row=~/\*\s+nIdEles\:\s+(\S+)$/) {
	    {$ppEleX_pt14_nIdEles_Zsel_tmp=$1;}
	}

	if($row=~/\*\s+nTkIsolEles\:\s+(\S+)$/) {
	    {$ppEleX_pt14_nTkIsolEles_Zsel_tmp=$1;}
	}

	if($row=~/\*\s+nEcalIsolEles\:\s+(\S+)$/) {
	    {$ppEleX_pt14_nEcalIsolEles_Zsel_tmp=$1;}
	}

	if($row=~/\*\s+nHcalIsolEles\:\s+(\S+)$/) {
	    {$ppEleX_pt14_nHcalIsolEles_Zsel_tmp=$1;}
	}

	if($row=~/\*\s+nPrimaryVertices\:\s+(\S+)$/) {
	    {$ppEleX_pt14_nPrimaryVertices_Zsel_tmp=$1;}
	}

	if($row=~/\*\s+nDzVertex\:\s+(\S+)$/) {
	    {$ppEleX_pt14_nDzVertex_Zsel_tmp=$1;}
	}

	if($row=~/\*\s+nDxySigVertex\:\s+(\S+)$/) {
	    {$ppEleX_pt14_nDxySigVertex_Zsel_tmp=$1;}
	}

	if($row=~/\*\s+meeCut\:\s+(\S+)$/) {
	    {$ppEleX_pt14_meeCut_Zsel_tmp=$1;}
	}

    }

    $ppEleX_pt14_event_Wsel += $ppEleX_pt14_event_Wsel_tmp;
    $ppEleX_pt14_trigger_Wsel += $ppEleX_pt14_trigger_Wsel_tmp;
    $ppEleX_pt14_nRecoEles_Wsel += $ppEleX_pt14_nRecoEles_Wsel_tmp;
    $ppEleX_pt14_nAccEles_Wsel += $ppEleX_pt14_nAccEles_Wsel_tmp;
    $ppEleX_pt14_nIdEles_Wsel += $ppEleX_pt14_nIdEles_Wsel_tmp;
    $ppEleX_pt14_nTkIsolEles_Wsel += $ppEleX_pt14_nTkIsolEles_Wsel_tmp;
    $ppEleX_pt14_nEcalIsolEles_Wsel += $ppEleX_pt14_nEcalIsolEles_Wsel_tmp;
    $ppEleX_pt14_nHcalIsolEles_Wsel += $ppEleX_pt14_nHcalIsolEles_Wsel_tmp;
    $ppEleX_pt14_nFinalEles_Wsel += $ppEleX_pt14_nFinalEles_Wsel_tmp;
    $ppEleX_pt14_ZVeto_Wsel += $ppEleX_pt14_ZVeto_Wsel_tmp;
    $ppEleX_pt14_nPrimaryVertices_Wsel += $ppEleX_pt14_nPrimaryVertices_Wsel_tmp;
    $ppEleX_pt14_nDzVertex_Wsel += $ppEleX_pt14_nDzVertex_Wsel_tmp;
    $ppEleX_pt14_nDxySigVertex_Wsel += $ppEleX_pt14_nDxySigVertex_Wsel_tmp;
    $ppEleX_pt14_metCut_Wsel += $ppEleX_pt14_metCut_Wsel_tmp;
    $ppEleX_pt14_transvMassCut_Wsel += $ppEleX_pt14_transvMassCut_Wsel_tmp;

    $ppEleX_pt14_event_Zsel += $ppEleX_pt14_event_Zsel_tmp;
    $ppEleX_pt14_trigger_Zsel += $ppEleX_pt14_trigger_Zsel_tmp;
    $ppEleX_pt14_nRecoEles_Zsel += $ppEleX_pt14_nRecoEles_Zsel_tmp;
    $ppEleX_pt14_nAccEles_Zsel += $ppEleX_pt14_nAccEles_Zsel_tmp;
    $ppEleX_pt14_nIdEles_Zsel += $ppEleX_pt14_nIdEles_Zsel_tmp;
    $ppEleX_pt14_nTkIsolEles_Zsel += $ppEleX_pt14_nTkIsolEles_Zsel_tmp;
    $ppEleX_pt14_nEcalIsolEles_Zsel += $ppEleX_pt14_nEcalIsolEles_Zsel_tmp;
    $ppEleX_pt14_nHcalIsolEles_Zsel += $ppEleX_pt14_nHcalIsolEles_Zsel_tmp;
    $ppEleX_pt14_nPrimaryVertices_Zsel += $ppEleX_pt14_nPrimaryVertices_Zsel_tmp;
    $ppEleX_pt14_nDzVertex_Zsel += $ppEleX_pt14_nDzVertex_Zsel_tmp;
    $ppEleX_pt14_nDxySigVertex_Zsel += $ppEleX_pt14_nDxySigVertex_Zsel_tmp;
    $ppEleX_pt14_meeCut_Zsel += $ppEleX_pt14_meeCut_Zsel_tmp;

}


print "now writing tables...\n";

#--- HLT ---   
print TEXFILEW "HLT &\n";
$decimals = 0;
$n_Wj_trigger_Wsel = sprintf("%.4g", $Wj_trigger_Wsel);
$eff_Wj_trigger_Wsel = sprintf("%.0f", 100 * $Wj_trigger_Wsel/$Wj_event_Wsel);

$n_Zj_trigger_Wsel = sprintf("%.4g", $Zj_trigger_Wsel);
$eff_Zj_trigger_Wsel = sprintf("%.0f", 100 * $Zj_trigger_Wsel/$Zj_event_Wsel);

$n_ttbarj_trigger_Wsel = sprintf("%.4g", $ttbarj_trigger_Wsel);
$eff_ttbarj_trigger_Wsel = sprintf("%.0f", 100 * $ttbarj_trigger_Wsel/$ttbarj_event_Wsel);

$n_ppEleX_pt14_trigger_Wsel = sprintf("%.4g", $ppEleX_pt14_trigger_Wsel/$ppEleX_pt14_event_Wsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
$eff_ppEleX_pt14_trigger_Wsel = sprintf("%.0f", 100 * $ppEleX_pt14_trigger_Wsel/$ppEleX_pt14_event_Wsel);

print TEXFILEW "$n_Wj_trigger_Wsel ($eff_Wj_trigger_Wsel \\%)&\n";
print TEXFILEW "$n_Zj_trigger_Wsel ($eff_Zj_trigger_Wsel \\%)&\n";
print TEXFILEW "$n_ttbarj_trigger_Wsel ($eff_ttbarj_trigger_Wsel \\%)&\n";
print TEXFILEW "$n_ppEleX_pt14_trigger_Wsel ($eff_ppEleX_pt14_trigger_Wsel \\%) \\\\ \n";


print TEXFILEZ "HLT &\n";
$n_Wj_trigger_Zsel = sprintf("%.4g", $Wj_trigger_Zsel);
$eff_Wj_trigger_Zsel = sprintf("%.0f", 100 * $Wj_trigger_Zsel/$Wj_event_Zsel);

$n_Zj_trigger_Zsel = sprintf("%.4g", $Zj_trigger_Zsel);
$eff_Zj_trigger_Zsel = sprintf("%.0f", 100 * $Zj_trigger_Zsel/$Zj_event_Zsel);

$n_ttbarj_trigger_Zsel = sprintf("%.4g", $ttbarj_trigger_Zsel);
$eff_ttbarj_trigger_Zsel = sprintf("%.0f", 100 * $ttbarj_trigger_Zsel/$ttbarj_event_Zsel);

$n_ppEleX_pt14_trigger_Zsel = sprintf("%.4g", $ppEleX_pt14_trigger_Zsel/$ppEleX_pt14_event_Zsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
$eff_ppEleX_pt14_trigger_Zsel = sprintf("%.0f", 100 * $ppEleX_pt14_trigger_Zsel/$ppEleX_pt14_event_Zsel);

print TEXFILEZ "$n_Wj_trigger_Zsel ($eff_Wj_trigger_Zsel \\%)&\n";
print TEXFILEZ "$n_Zj_trigger_Zsel ($eff_Zj_trigger_Zsel \\%)&\n";
print TEXFILEZ "$n_ttbarj_trigger_Zsel ($eff_ttbarj_trigger_Zsel \\%)&\n";
print TEXFILEZ "$n_ppEleX_pt14_trigger_Zsel ($eff_ppEleX_pt14_trigger_Zsel \\%) \\\\ \n";



#--- electron reconstruction ---   
print TEXFILEW "\$e\^\{\\pm\}\$ reconstruction &\n";
$decimals = 0;
$n_Wj_nRecoEles_Wsel = sprintf("%.4g", $Wj_nRecoEles_Wsel);
$eff_Wj_nRecoEles_Wsel = sprintf("%.0f", 100 * $Wj_nRecoEles_Wsel/$Wj_trigger_Wsel);

$n_Zj_nRecoEles_Wsel = sprintf("%.4g", $Zj_nRecoEles_Wsel);
$eff_Zj_nRecoEles_Wsel = sprintf("%.0f", 100 * $Zj_nRecoEles_Wsel/$Zj_trigger_Wsel);

$n_ttbarj_nRecoEles_Wsel = sprintf("%.4g", $ttbarj_nRecoEles_Wsel);
$eff_ttbarj_nRecoEles_Wsel = sprintf("%.0f", 100 * $ttbarj_nRecoEles_Wsel/$ttbarj_trigger_Wsel);

$n_ppEleX_pt14_nRecoEles_Wsel = sprintf("%.4g", $ppEleX_pt14_nRecoEles_Wsel/$ppEleX_pt14_event_Wsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
$eff_ppEleX_pt14_nRecoEles_Wsel = sprintf("%.0f", 100 * $ppEleX_pt14_nRecoEles_Wsel/$ppEleX_pt14_trigger_Wsel);

print TEXFILEW "$n_Wj_nRecoEles_Wsel ($eff_Wj_nRecoEles_Wsel \\%)&\n";
print TEXFILEW "$n_Zj_nRecoEles_Wsel ($eff_Zj_nRecoEles_Wsel \\%)&\n";
print TEXFILEW "$n_ttbarj_nRecoEles_Wsel ($eff_ttbarj_nRecoEles_Wsel \\%)&\n";
print TEXFILEW "$n_ppEleX_pt14_nRecoEles_Wsel ($eff_ppEleX_pt14_nRecoEles_Wsel \\%) \\\\ \n";


print TEXFILEZ "\$e\^\{\\pm\}\$ reconstruction &\n";
$n_Wj_nRecoEles_Zsel = sprintf("%.4g", $Wj_nRecoEles_Zsel);
$eff_Wj_nRecoEles_Zsel = sprintf("%.0f", 100 * $Wj_nRecoEles_Zsel/$Wj_trigger_Zsel);

$n_Zj_nRecoEles_Zsel = sprintf("%.4g", $Zj_nRecoEles_Zsel);
$eff_Zj_nRecoEles_Zsel = sprintf("%.0f", 100 * $Zj_nRecoEles_Zsel/$Zj_trigger_Zsel);

$n_ttbarj_nRecoEles_Zsel = sprintf("%.4g", $ttbarj_nRecoEles_Zsel);
$eff_ttbarj_nRecoEles_Zsel = sprintf("%.0f", 100 * $ttbarj_nRecoEles_Zsel/$ttbarj_trigger_Zsel);

$n_ppEleX_pt14_nRecoEles_Zsel = sprintf("%.4g", $ppEleX_pt14_nRecoEles_Zsel/$ppEleX_pt14_event_Zsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
$eff_ppEleX_pt14_nRecoEles_Zsel = sprintf("%.0f", 100 * $ppEleX_pt14_nRecoEles_Zsel/$ppEleX_pt14_trigger_Zsel);

print TEXFILEZ "$n_Wj_nRecoEles_Zsel ($eff_Wj_nRecoEles_Zsel \\%)&\n";
print TEXFILEZ "$n_Zj_nRecoEles_Zsel ($eff_Zj_nRecoEles_Zsel \\%)&\n";
print TEXFILEZ "$n_ttbarj_nRecoEles_Zsel ($eff_ttbarj_nRecoEles_Zsel \\%)&\n";
print TEXFILEZ "$n_ppEleX_pt14_nRecoEles_Zsel ($eff_ppEleX_pt14_nRecoEles_Zsel \\%) \\\\ \n";



#--- electron acceptance and pT min ---   
print TEXFILEW "\$e\^\{\\pm\}\$ \$|\\eta|\<2.5\$ and \$p_T\>20\$ GeV/c &\n";
$decimals = 0;
$n_Wj_nAccEles_Wsel = sprintf("%.4g", $Wj_nAccEles_Wsel);
$eff_Wj_nAccEles_Wsel = sprintf("%.0f", 100 * $Wj_nAccEles_Wsel/$Wj_nRecoEles_Wsel);

$n_Zj_nAccEles_Wsel = sprintf("%.4g", $Zj_nAccEles_Wsel);
$eff_Zj_nAccEles_Wsel = sprintf("%.0f", 100 * $Zj_nAccEles_Wsel/$Zj_nRecoEles_Wsel);

$n_ttbarj_nAccEles_Wsel = sprintf("%.4g", $ttbarj_nAccEles_Wsel);
$eff_ttbarj_nAccEles_Wsel = sprintf("%.0f", 100 * $ttbarj_nAccEles_Wsel/$ttbarj_nRecoEles_Wsel);

$n_ppEleX_pt14_nAccEles_Wsel = sprintf("%.4g", $ppEleX_pt14_nAccEles_Wsel/$ppEleX_pt14_event_Wsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
$eff_ppEleX_pt14_nAccEles_Wsel = sprintf("%.0f", 100 * $ppEleX_pt14_nAccEles_Wsel/$ppEleX_pt14_nRecoEles_Wsel);

print TEXFILEW "$n_Wj_nAccEles_Wsel ($eff_Wj_nAccEles_Wsel \\%)&\n";
print TEXFILEW "$n_Zj_nAccEles_Wsel ($eff_Zj_nAccEles_Wsel \\%)&\n";
print TEXFILEW "$n_ttbarj_nAccEles_Wsel ($eff_ttbarj_nAccEles_Wsel \\%)&\n";
print TEXFILEW "$n_ppEleX_pt14_nAccEles_Wsel ($eff_ppEleX_pt14_nAccEles_Wsel \\%) \\\\ \n";


print TEXFILEZ "\$e\^\{\\pm\}\$ \$|\\eta|\<2.5\$ and \$p_T\>20\$ GeV/c &\n";
$n_Wj_nAccEles_Zsel = sprintf("%.4g", $Wj_nAccEles_Zsel);
$eff_Wj_nAccEles_Zsel = sprintf("%.0f", 100 * $Wj_nAccEles_Zsel/$Wj_nRecoEles_Zsel);

$n_Zj_nAccEles_Zsel = sprintf("%.4g", $Zj_nAccEles_Zsel);
$eff_Zj_nAccEles_Zsel = sprintf("%.0f", 100 * $Zj_nAccEles_Zsel/$Zj_nRecoEles_Zsel);

$n_ttbarj_nAccEles_Zsel = sprintf("%.4g", $ttbarj_nAccEles_Zsel);
$eff_ttbarj_nAccEles_Zsel = sprintf("%.0f", 100 * $ttbarj_nAccEles_Zsel/$ttbarj_nRecoEles_Zsel);

$n_ppEleX_pt14_nAccEles_Zsel = sprintf("%.4g", $ppEleX_pt14_nAccEles_Zsel/$ppEleX_pt14_event_Zsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
$eff_ppEleX_pt14_nAccEles_Zsel = sprintf("%.0f", 100 * $ppEleX_pt14_nAccEles_Zsel/$ppEleX_pt14_nRecoEles_Zsel);

print TEXFILEZ "$n_Wj_nAccEles_Zsel ($eff_Wj_nAccEles_Zsel \\%)&\n";
print TEXFILEZ "$n_Zj_nAccEles_Zsel ($eff_Zj_nAccEles_Zsel \\%)&\n";
print TEXFILEZ "$n_ttbarj_nAccEles_Zsel ($eff_ttbarj_nAccEles_Zsel \\%)&\n";
print TEXFILEZ "$n_ppEleX_pt14_nAccEles_Zsel ($eff_ppEleX_pt14_nAccEles_Zsel \\%) \\\\ \n";



#--- electron identification ---   
print TEXFILEW "electron ID &\n";
$decimals = 0;
$n_Wj_nIdEles_Wsel = sprintf("%.4g", $Wj_nIdEles_Wsel);
$eff_Wj_nIdEles_Wsel = sprintf("%.0f", 100 * $Wj_nIdEles_Wsel/$Wj_nAccEles_Wsel);

$n_Zj_nIdEles_Wsel = sprintf("%.4g", $Zj_nIdEles_Wsel);
$eff_Zj_nIdEles_Wsel = sprintf("%.0f", 100 * $Zj_nIdEles_Wsel/$Zj_nAccEles_Wsel);

$n_ttbarj_nIdEles_Wsel = sprintf("%.4g", $ttbarj_nIdEles_Wsel);
$eff_ttbarj_nIdEles_Wsel = sprintf("%.0f", 100 * $ttbarj_nIdEles_Wsel/$ttbarj_nAccEles_Wsel);

$n_ppEleX_pt14_nIdEles_Wsel = sprintf("%.4g", $ppEleX_pt14_nIdEles_Wsel/$ppEleX_pt14_event_Wsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
$eff_ppEleX_pt14_nIdEles_Wsel = sprintf("%.0f", 100 * $ppEleX_pt14_nIdEles_Wsel/$ppEleX_pt14_nAccEles_Wsel);

print TEXFILEW "$n_Wj_nIdEles_Wsel ($eff_Wj_nIdEles_Wsel \\%)&\n";
print TEXFILEW "$n_Zj_nIdEles_Wsel ($eff_Zj_nIdEles_Wsel \\%)&\n";
print TEXFILEW "$n_ttbarj_nIdEles_Wsel ($eff_ttbarj_nIdEles_Wsel \\%)&\n";
print TEXFILEW "$n_ppEleX_pt14_nIdEles_Wsel ($eff_ppEleX_pt14_nIdEles_Wsel \\%) \\\\ \n";


print TEXFILEZ "electron ID &\n";
$n_Wj_nIdEles_Zsel = sprintf("%.4g", $Wj_nIdEles_Zsel);
$eff_Wj_nIdEles_Zsel = sprintf("%.0f", 100 * $Wj_nIdEles_Zsel/$Wj_nAccEles_Zsel);

$n_Zj_nIdEles_Zsel = sprintf("%.4g", $Zj_nIdEles_Zsel);
$eff_Zj_nIdEles_Zsel = sprintf("%.0f", 100 * $Zj_nIdEles_Zsel/$Zj_nAccEles_Zsel);

$n_ttbarj_nIdEles_Zsel = sprintf("%.4g", $ttbarj_nIdEles_Zsel);
$eff_ttbarj_nIdEles_Zsel = sprintf("%.0f", 100 * $ttbarj_nIdEles_Zsel/$ttbarj_nAccEles_Zsel);

$n_ppEleX_pt14_nIdEles_Zsel = sprintf("%.4g", $ppEleX_pt14_nIdEles_Zsel/$ppEleX_pt14_event_Zsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
$eff_ppEleX_pt14_nIdEles_Zsel = sprintf("%.0f", 100 * $ppEleX_pt14_nIdEles_Zsel/$ppEleX_pt14_nAccEles_Zsel);

print TEXFILEZ "$n_Wj_nIdEles_Zsel ($eff_Wj_nIdEles_Zsel \\%)&\n";
print TEXFILEZ "$n_Zj_nIdEles_Zsel ($eff_Zj_nIdEles_Zsel \\%)&\n";
print TEXFILEZ "$n_ttbarj_nIdEles_Zsel ($eff_ttbarj_nIdEles_Zsel \\%)&\n";
print TEXFILEZ "$n_ppEleX_pt14_nIdEles_Zsel ($eff_ppEleX_pt14_nIdEles_Zsel \\%) \\\\ \n";




#--- electron tracker isolation ---   
print TEXFILEW "tracker isolation &\n";
$decimals = 0;
$n_Wj_nTkIsolEles_Wsel = sprintf("%.4g", $Wj_nTkIsolEles_Wsel);
$eff_Wj_nTkIsolEles_Wsel = sprintf("%.0f", 100 * $Wj_nTkIsolEles_Wsel/$Wj_nIdEles_Wsel);

$n_Zj_nTkIsolEles_Wsel = sprintf("%.4g", $Zj_nTkIsolEles_Wsel);
$eff_Zj_nTkIsolEles_Wsel = sprintf("%.0f", 100 * $Zj_nTkIsolEles_Wsel/$Zj_nIdEles_Wsel);

$n_ttbarj_nTkIsolEles_Wsel = sprintf("%.4g", $ttbarj_nTkIsolEles_Wsel);
$eff_ttbarj_nTkIsolEles_Wsel = sprintf("%.0f", 100 * $ttbarj_nTkIsolEles_Wsel/$ttbarj_nIdEles_Wsel);

$n_ppEleX_pt14_nTkIsolEles_Wsel = sprintf("%.4g", $ppEleX_pt14_nTkIsolEles_Wsel/$ppEleX_pt14_event_Wsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
$eff_ppEleX_pt14_nTkIsolEles_Wsel = sprintf("%.0f", 100 * $ppEleX_pt14_nTkIsolEles_Wsel/$ppEleX_pt14_nIdEles_Wsel);

print TEXFILEW "$n_Wj_nTkIsolEles_Wsel ($eff_Wj_nTkIsolEles_Wsel \\%)&\n";
print TEXFILEW "$n_Zj_nTkIsolEles_Wsel ($eff_Zj_nTkIsolEles_Wsel \\%)&\n";
print TEXFILEW "$n_ttbarj_nTkIsolEles_Wsel ($eff_ttbarj_nTkIsolEles_Wsel \\%)&\n";
print TEXFILEW "$n_ppEleX_pt14_nTkIsolEles_Wsel ($eff_ppEleX_pt14_nTkIsolEles_Wsel \\%) \\\\ \n";


print TEXFILEZ "tracker isolation &\n";
$n_Wj_nTkIsolEles_Zsel = sprintf("%.4g", $Wj_nTkIsolEles_Zsel);
$eff_Wj_nTkIsolEles_Zsel = sprintf("%.0f", 100 * $Wj_nTkIsolEles_Zsel/$Wj_nIdEles_Zsel);

$n_Zj_nTkIsolEles_Zsel = sprintf("%.4g", $Zj_nTkIsolEles_Zsel);
$eff_Zj_nTkIsolEles_Zsel = sprintf("%.0f", 100 * $Zj_nTkIsolEles_Zsel/$Zj_nIdEles_Zsel);

$n_ttbarj_nTkIsolEles_Zsel = sprintf("%.4g", $ttbarj_nTkIsolEles_Zsel);
$eff_ttbarj_nTkIsolEles_Zsel = sprintf("%.0f", 100 * $ttbarj_nTkIsolEles_Zsel/$ttbarj_nIdEles_Zsel);

$n_ppEleX_pt14_nTkIsolEles_Zsel = sprintf("%.4g", $ppEleX_pt14_nTkIsolEles_Zsel/$ppEleX_pt14_event_Zsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
$eff_ppEleX_pt14_nTkIsolEles_Zsel = sprintf("%.0f", 100 * $ppEleX_pt14_nTkIsolEles_Zsel/$ppEleX_pt14_nIdEles_Zsel);

print TEXFILEZ "$n_Wj_nTkIsolEles_Zsel ($eff_Wj_nTkIsolEles_Zsel \\%)&\n";
print TEXFILEZ "$n_Zj_nTkIsolEles_Zsel ($eff_Zj_nTkIsolEles_Zsel \\%)&\n";
print TEXFILEZ "$n_ttbarj_nTkIsolEles_Zsel ($eff_ttbarj_nTkIsolEles_Zsel \\%)&\n";
print TEXFILEZ "$n_ppEleX_pt14_nTkIsolEles_Zsel ($eff_ppEleX_pt14_nTkIsolEles_Zsel \\%) \\\\ \n";



#--- electron ECAL isolation ---   
print TEXFILEW "ECAL isolation &\n";
$decimals = 0;
$n_Wj_nEcalIsolEles_Wsel = sprintf("%.4g", $Wj_nEcalIsolEles_Wsel);
$eff_Wj_nEcalIsolEles_Wsel = sprintf("%.0f", 100 * $Wj_nEcalIsolEles_Wsel/$Wj_nTkIsolEles_Wsel);

$n_Zj_nEcalIsolEles_Wsel = sprintf("%.4g", $Zj_nEcalIsolEles_Wsel);
$eff_Zj_nEcalIsolEles_Wsel = sprintf("%.0f", 100 * $Zj_nEcalIsolEles_Wsel/$Zj_nTkIsolEles_Wsel);

$n_ttbarj_nEcalIsolEles_Wsel = sprintf("%.4g", $ttbarj_nEcalIsolEles_Wsel);
$eff_ttbarj_nEcalIsolEles_Wsel = sprintf("%.0f", 100 * $ttbarj_nEcalIsolEles_Wsel/$ttbarj_nTkIsolEles_Wsel);

$n_ppEleX_pt14_nEcalIsolEles_Wsel = sprintf("%.4g", $ppEleX_pt14_nEcalIsolEles_Wsel/$ppEleX_pt14_event_Wsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
$eff_ppEleX_pt14_nEcalIsolEles_Wsel = sprintf("%.0f", 100 * $ppEleX_pt14_nEcalIsolEles_Wsel/$ppEleX_pt14_nTkIsolEles_Wsel);

print TEXFILEW "$n_Wj_nEcalIsolEles_Wsel ($eff_Wj_nEcalIsolEles_Wsel \\%)&\n";
print TEXFILEW "$n_Zj_nEcalIsolEles_Wsel ($eff_Zj_nEcalIsolEles_Wsel \\%)&\n";
print TEXFILEW "$n_ttbarj_nEcalIsolEles_Wsel ($eff_ttbarj_nEcalIsolEles_Wsel \\%)&\n";
print TEXFILEW "$n_ppEleX_pt14_nEcalIsolEles_Wsel ($eff_ppEleX_pt14_nEcalIsolEles_Wsel \\%) \\\\ \n";


print TEXFILEZ "ECAL isolation &\n";
$n_Wj_nEcalIsolEles_Zsel = sprintf("%.4g", $Wj_nEcalIsolEles_Zsel);
$eff_Wj_nEcalIsolEles_Zsel = sprintf("%.0f", 100 * $Wj_nEcalIsolEles_Zsel/$Wj_nTkIsolEles_Zsel);

$n_Zj_nEcalIsolEles_Zsel = sprintf("%.4g", $Zj_nEcalIsolEles_Zsel);
$eff_Zj_nEcalIsolEles_Zsel = sprintf("%.0f", 100 * $Zj_nEcalIsolEles_Zsel/$Zj_nTkIsolEles_Zsel);

$n_ttbarj_nEcalIsolEles_Zsel = sprintf("%.4g", $ttbarj_nEcalIsolEles_Zsel);
$eff_ttbarj_nEcalIsolEles_Zsel = sprintf("%.0f", 100 * $ttbarj_nEcalIsolEles_Zsel/$ttbarj_nTkIsolEles_Zsel);

$n_ppEleX_pt14_nEcalIsolEles_Zsel = sprintf("%.4g", $ppEleX_pt14_nEcalIsolEles_Zsel/$ppEleX_pt14_event_Zsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
$eff_ppEleX_pt14_nEcalIsolEles_Zsel = sprintf("%.0f", 100 * $ppEleX_pt14_nEcalIsolEles_Zsel/$ppEleX_pt14_nTkIsolEles_Zsel);

print TEXFILEZ "$n_Wj_nEcalIsolEles_Zsel ($eff_Wj_nEcalIsolEles_Zsel \\%)&\n";
print TEXFILEZ "$n_Zj_nEcalIsolEles_Zsel ($eff_Zj_nEcalIsolEles_Zsel \\%)&\n";
print TEXFILEZ "$n_ttbarj_nEcalIsolEles_Zsel ($eff_ttbarj_nEcalIsolEles_Zsel \\%)&\n";
print TEXFILEZ "$n_ppEleX_pt14_nEcalIsolEles_Zsel ($eff_ppEleX_pt14_nEcalIsolEles_Zsel \\%) \\\\ \n";




#--- electron HCAL isolation ---   
print TEXFILEW "HCAL isolation &\n";
$decimals = 0;
$n_Wj_nHcalIsolEles_Wsel = sprintf("%.4g", $Wj_nHcalIsolEles_Wsel);
if($Wj_nEcalIsolEles_Wsel!=0) {
    $eff_Wj_nHcalIsolEles_Wsel = sprintf("%.0f", 100 * $Wj_nHcalIsolEles_Wsel/$Wj_nEcalIsolEles_Wsel);
} else {
    $eff_Wj_nHcalIsolEles_Wsel=0;
}

$n_Zj_nHcalIsolEles_Wsel = sprintf("%.4g", $Zj_nHcalIsolEles_Wsel);
if($Zj_nEcalIsolEles_Wsel!=0) {
    $eff_Zj_nHcalIsolEles_Wsel = sprintf("%.0f", 100 * $Zj_nHcalIsolEles_Wsel/$Zj_nEcalIsolEles_Wsel);
} else {
    $Zj_nHcalIsolEles_Wsel = 0;
}

$n_ttbarj_nHcalIsolEles_Wsel = sprintf("%.4g", $ttbarj_nHcalIsolEles_Wsel);
if($ttbarj_nEcalIsolEles_Wsel!=0) {
    $eff_ttbarj_nHcalIsolEles_Wsel = sprintf("%.0f", 100 * $ttbarj_nHcalIsolEles_Wsel/$ttbarj_nEcalIsolEles_Wsel);
} else {
    $ttbarj_nHcalIsolEles_Wsel = 0;
}

$n_ppEleX_pt14_nHcalIsolEles_Wsel = sprintf("%.4g", $ppEleX_pt14_nHcalIsolEles_Wsel/$ppEleX_pt14_event_Wsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
if($ppEleX_pt14_nEcalIsolEles_Wsel!=0) {
    $eff_ppEleX_pt14_nHcalIsolEles_Wsel = sprintf("%.0f", 100 * $ppEleX_pt14_nHcalIsolEles_Wsel/$ppEleX_pt14_nEcalIsolEles_Wsel);
} else {
    $ppEleX_pt14_nHcalIsolEles_Wsel = 0;
}
print TEXFILEW "$n_Wj_nHcalIsolEles_Wsel ($eff_Wj_nHcalIsolEles_Wsel \\%)&\n";
print TEXFILEW "$n_Zj_nHcalIsolEles_Wsel ($eff_Zj_nHcalIsolEles_Wsel \\%)&\n";
print TEXFILEW "$n_ttbarj_nHcalIsolEles_Wsel ($eff_ttbarj_nHcalIsolEles_Wsel \\%)&\n";
print TEXFILEW "$n_ppEleX_pt14_nHcalIsolEles_Wsel ($eff_ppEleX_pt14_nHcalIsolEles_Wsel \\%) \\\\ \n";


print TEXFILEZ "HCAL isolation &\n";
$n_Wj_nHcalIsolEles_Zsel = sprintf("%.4g", $Wj_nHcalIsolEles_Zsel);
if($Wj_nEcalIsolEles_Zsel!=0) {
    $eff_Wj_nHcalIsolEles_Zsel = sprintf("%.0f", 100 * $Wj_nHcalIsolEles_Zsel/$Wj_nEcalIsolEles_Zsel);
} else {
    $Wj_nHcalIsolEles_Zsel = 0;
}

$n_Zj_nHcalIsolEles_Zsel = sprintf("%.4g", $Zj_nHcalIsolEles_Zsel);
if($Zj_nEcalIsolEles_Zsel!=0) {
    $eff_Zj_nHcalIsolEles_Zsel = sprintf("%.0f", 100 * $Zj_nHcalIsolEles_Zsel/$Zj_nEcalIsolEles_Zsel);
} else {
    $Zj_nHcalIsolEles_Zsel = 0;
}

$n_ttbarj_nHcalIsolEles_Zsel = sprintf("%.4g", $ttbarj_nHcalIsolEles_Zsel);
if($ttbarj_nEcalIsolEles_Zsel!=0) {
    $eff_ttbarj_nHcalIsolEles_Zsel = sprintf("%.0f", 100 * $ttbarj_nHcalIsolEles_Zsel/$ttbarj_nEcalIsolEles_Zsel);
} else {
    $ttbarj_nHcalIsolEles_Zsel = 0;
}

$n_ppEleX_pt14_nHcalIsolEles_Zsel = sprintf("%.4g", $ppEleX_pt14_nHcalIsolEles_Zsel/$ppEleX_pt14_event_Zsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
if($ppEleX_pt14_nEcalIsolEles_Zsel!=0) {
    $eff_ppEleX_pt14_nHcalIsolEles_Zsel = sprintf("%.0f", 100 * $ppEleX_pt14_nHcalIsolEles_Zsel/$ppEleX_pt14_nEcalIsolEles_Zsel);
} else {
    $ppEleX_pt14_nHcalIsolEles_Zsel = 0;
}
print TEXFILEZ "$n_Wj_nHcalIsolEles_Zsel ($eff_Wj_nHcalIsolEles_Zsel \\%)&\n";
print TEXFILEZ "$n_Zj_nHcalIsolEles_Zsel ($eff_Zj_nHcalIsolEles_Zsel \\%)&\n";
print TEXFILEZ "$n_ttbarj_nHcalIsolEles_Zsel ($eff_ttbarj_nHcalIsolEles_Zsel \\%)&\n";
print TEXFILEZ "$n_ppEleX_pt14_nHcalIsolEles_Zsel ($eff_ppEleX_pt14_nHcalIsolEles_Zsel \\%) \\\\ \n";



#--- only for W selection: # of final electrons exactly 1 ---   
print TEXFILEW "exactly 1 good electron &\n";
$decimals = 0;
$n_Wj_nFinalEles_Wsel = sprintf("%.4g", $Wj_nFinalEles_Wsel);
if($Wj_nHcalIsolEles_Wsel!=0) {
    $eff_Wj_nFinalEles_Wsel = sprintf("%.0f", 100 * $Wj_nFinalEles_Wsel/$Wj_nHcalIsolEles_Wsel);
} else {
    $eff_Wj_nFinalEles_Wsel=0;
}

$n_Zj_nFinalEles_Wsel = sprintf("%.4g", $Zj_nFinalEles_Wsel);
if($Zj_nHcalIsolEles_Wsel!=0) {
    $eff_Zj_nFinalEles_Wsel = sprintf("%.0f", 100 * $Zj_nFinalEles_Wsel/$Zj_nHcalIsolEles_Wsel);
} else {
    $Zj_nFinalEles_Wsel = 0;
}

$n_ttbarj_nFinalEles_Wsel = sprintf("%.4g", $ttbarj_nFinalEles_Wsel);
if($ttbarj_nHcalIsolEles_Wsel!=0) {
    $eff_ttbarj_nFinalEles_Wsel = sprintf("%.0f", 100 * $ttbarj_nFinalEles_Wsel/$ttbarj_nHcalIsolEles_Wsel);
} else {
    $ttbarj_nFinalEles_Wsel = 0;
}

$n_ppEleX_pt14_nFinalEles_Wsel = sprintf("%.4g", $ppEleX_pt14_nFinalEles_Wsel/$ppEleX_pt14_event_Wsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
if($ppEleX_pt14_nHcalIsolEles_Wsel!=0) {
    $eff_ppEleX_pt14_nFinalEles_Wsel = sprintf("%.0f", 100 * $ppEleX_pt14_nFinalEles_Wsel/$ppEleX_pt14_nHcalIsolEles_Wsel);
} else {
    $ppEleX_pt14_nFinalEles_Wsel = 0;
}
print TEXFILEW "$n_Wj_nFinalEles_Wsel ($eff_Wj_nFinalEles_Wsel \\%)&\n";
print TEXFILEW "$n_Zj_nFinalEles_Wsel ($eff_Zj_nFinalEles_Wsel \\%)&\n";
print TEXFILEW "$n_ttbarj_nFinalEles_Wsel ($eff_ttbarj_nFinalEles_Wsel \\%)&\n";
print TEXFILEW "$n_ppEleX_pt14_nFinalEles_Wsel ($eff_ppEleX_pt14_nFinalEles_Wsel \\%) \\\\ \n";



#--- only for W selection: the good ele + one of the reco ele do not make Z mass ---   
print TEXFILEW "Z mass veto &\n";
$decimals = 0;
$n_Wj_ZVeto_Wsel = sprintf("%.4g", $Wj_ZVeto_Wsel);
if($Wj_nFinalEles_Wsel!=0) {
    $eff_Wj_ZVeto_Wsel = sprintf("%.0f", 100 * $Wj_ZVeto_Wsel/$Wj_nFinalEles_Wsel);
} else {
    $eff_Wj_ZVeto_Wsel=0;
}

$n_Zj_ZVeto_Wsel = sprintf("%.4g", $Zj_ZVeto_Wsel);
if($Zj_nFinalEles_Wsel!=0) {
    $eff_Zj_ZVeto_Wsel = sprintf("%.0f", 100 * $Zj_ZVeto_Wsel/$Zj_nFinalEles_Wsel);
} else {
    $Zj_ZVeto_Wsel = 0;
}

$n_ttbarj_ZVeto_Wsel = sprintf("%.4g", $ttbarj_ZVeto_Wsel);
if($ttbarj_nFinalEles_Wsel!=0) {
    $eff_ttbarj_ZVeto_Wsel = sprintf("%.0f", 100 * $ttbarj_ZVeto_Wsel/$ttbarj_nFinalEles_Wsel);
} else {
    $ttbarj_ZVeto_Wsel = 0;
}

$n_ppEleX_pt14_ZVeto_Wsel = sprintf("%.4g", $ppEleX_pt14_ZVeto_Wsel/$ppEleX_pt14_event_Wsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
if($ppEleX_pt14_nFinalEles_Wsel!=0) {
    $eff_ppEleX_pt14_ZVeto_Wsel = sprintf("%.0f", 100 * $ppEleX_pt14_ZVeto_Wsel/$ppEleX_pt14_nFinalEles_Wsel);
} else {
    $ppEleX_pt14_ZVeto_Wsel = 0;
}
print TEXFILEW "$n_Wj_ZVeto_Wsel ($eff_Wj_ZVeto_Wsel \\%)&\n";
print TEXFILEW "$n_Zj_ZVeto_Wsel ($eff_Zj_ZVeto_Wsel \\%)&\n";
print TEXFILEW "$n_ttbarj_ZVeto_Wsel ($eff_ttbarj_ZVeto_Wsel \\%)&\n";
print TEXFILEW "$n_ppEleX_pt14_ZVeto_Wsel ($eff_ppEleX_pt14_ZVeto_Wsel \\%) \\\\ \n";





#--- existence of PV ---   
print TEXFILEW "existence of PV &\n";
$decimals = 0;
$n_Wj_nPrimaryVertices_Wsel = sprintf("%.4g", $Wj_nPrimaryVertices_Wsel);
if($Wj_ZVeto_Wsel!=0) {
    $eff_Wj_nPrimaryVertices_Wsel = sprintf("%.0f", 100 * $Wj_nPrimaryVertices_Wsel/$Wj_ZVeto_Wsel);
} else {
    $eff_Wj_nPrimaryVertices_Wsel=0;
}

$n_Zj_nPrimaryVertices_Wsel = sprintf("%.4g", $Zj_nPrimaryVertices_Wsel);
if($Zj_ZVeto_Wsel!=0) {
    $eff_Zj_nPrimaryVertices_Wsel = sprintf("%.0f", 100 * $Zj_nPrimaryVertices_Wsel/$Zj_ZVeto_Wsel);
} else {
    $Zj_nPrimaryVertices_Wsel = 0;
}

$n_ttbarj_nPrimaryVertices_Wsel = sprintf("%.4g", $ttbarj_nPrimaryVertices_Wsel);
if($ttbarj_ZVeto_Wsel!=0) {
    $eff_ttbarj_nPrimaryVertices_Wsel = sprintf("%.0f", 100 * $ttbarj_nPrimaryVertices_Wsel/$ttbarj_ZVeto_Wsel);
} else {
    $ttbarj_nPrimaryVertices_Wsel = 0;
}

$n_ppEleX_pt14_nPrimaryVertices_Wsel = sprintf("%.4g", $ppEleX_pt14_nPrimaryVertices_Wsel/$ppEleX_pt14_event_Wsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
if($ppEleX_pt14_ZVeto_Wsel!=0) {
    $eff_ppEleX_pt14_nPrimaryVertices_Wsel = sprintf("%.0f", 100 * $ppEleX_pt14_nPrimaryVertices_Wsel/$ppEleX_pt14_ZVeto_Wsel);
} else {
    $ppEleX_pt14_nPrimaryVertices_Wsel = 0;
}
print TEXFILEW "$n_Wj_nPrimaryVertices_Wsel ($eff_Wj_nPrimaryVertices_Wsel \\%)&\n";
print TEXFILEW "$n_Zj_nPrimaryVertices_Wsel ($eff_Zj_nPrimaryVertices_Wsel \\%)&\n";
print TEXFILEW "$n_ttbarj_nPrimaryVertices_Wsel ($eff_ttbarj_nPrimaryVertices_Wsel \\%)&\n";
print TEXFILEW "$n_ppEleX_pt14_nPrimaryVertices_Wsel ($eff_ppEleX_pt14_nPrimaryVertices_Wsel \\%) \\\\ \n";


print TEXFILEZ "existence of PV &\n";
$n_Wj_nPrimaryVertices_Zsel = sprintf("%.4g", $Wj_nPrimaryVertices_Zsel);
if($Wj_nHcalIsolEles_Zsel!=0) {
    $eff_Wj_nPrimaryVertices_Zsel = sprintf("%.0f", 100 * $Wj_nPrimaryVertices_Zsel/$Wj_nHcalIsolEles_Zsel);
} else {
    $eff_Wj_nPrimaryVertices_Zsel = 0;
}

$n_Zj_nPrimaryVertices_Zsel = sprintf("%.4g", $Zj_nPrimaryVertices_Zsel);
if($Zj_nHcalIsolEles_Zsel!=0) {
    $eff_Zj_nPrimaryVertices_Zsel = sprintf("%.0f", 100 * $Zj_nPrimaryVertices_Zsel/$Zj_nHcalIsolEles_Zsel);
} else {
    $eff_Zj_nPrimaryVertices_Zsel = 0;
}

$n_ttbarj_nPrimaryVertices_Zsel = sprintf("%.4g", $ttbarj_nPrimaryVertices_Zsel);
if($ttbarj_nHcalIsolEles_Zsel!=0) {
    $eff_ttbarj_nPrimaryVertices_Zsel = sprintf("%.0f", 100 * $ttbarj_nPrimaryVertices_Zsel/$ttbarj_nHcalIsolEles_Zsel);
} else {
    $eff_ttbarj_nPrimaryVertices_Zsel = 0;
}

$n_ppEleX_pt14_nPrimaryVertices_Zsel = sprintf("%.4g", $ppEleX_pt14_nPrimaryVertices_Zsel/$ppEleX_pt14_event_Zsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
if($ppEleX_pt14_nHcalIsolEles_Zsel!=0) {
    $eff_ppEleX_pt14_nPrimaryVertices_Zsel = sprintf("%.0f", 100 * $ppEleX_pt14_nPrimaryVertices_Zsel/$ppEleX_pt14_nHcalIsolEles_Zsel);
} else {
    $eff_ppEleX_pt14_nPrimaryVertices_Zsel = 0;
}
print TEXFILEZ "$n_Wj_nPrimaryVertices_Zsel ($eff_Wj_nPrimaryVertices_Zsel \\%)&\n";
print TEXFILEZ "$n_Zj_nPrimaryVertices_Zsel ($eff_Zj_nPrimaryVertices_Zsel \\%)&\n";
print TEXFILEZ "$n_ttbarj_nPrimaryVertices_Zsel ($eff_ttbarj_nPrimaryVertices_Zsel \\%)&\n";
print TEXFILEZ "$n_ppEleX_pt14_nPrimaryVertices_Zsel ($eff_ppEleX_pt14_nPrimaryVertices_Zsel \\%) \\\\ \n";











#--- deltaz cut ---   
print TEXFILEW "\$\\vert z_\{e\} - z_\{PV\}\\vert < 0.02 cm\$&\n";
$decimals = 0;
$n_Wj_nDzVertex_Wsel = sprintf("%.4g", $Wj_nDzVertex_Wsel);
if($Wj_nPrimaryVertices_Wsel!=0) {
    $eff_Wj_nDzVertex_Wsel = sprintf("%.0f", 100 * $Wj_nDzVertex_Wsel/$Wj_nPrimaryVertices_Wsel);
} else {
    $eff_Wj_nDzVertex_Wsel=0;
}

$n_Zj_nDzVertex_Wsel = sprintf("%.4g", $Zj_nDzVertex_Wsel);
if($Zj_nDzVertex_Wsel!=0) {
    $eff_Zj_nPrimaryVertices_Wsel = sprintf("%.0f", 100 * $Zj_nDzVertex_Wsel/$Zj_nPrimaryVertices_Wsel);
} else {
    $Zj_nDzVertex_Wsel = 0;
}

$n_ttbarj_nDzVertex_Wsel = sprintf("%.4g", $ttbarj_nDzVertex_Wsel);
if($ttbarj_nPrimaryVertices_Wsel!=0) {
    $eff_ttbarj_nDzVertex_Wsel = sprintf("%.0f", 100 * $ttbarj_nDzVertex_Wsel/$ttbarj_nPrimaryVertices_Wsel);
} else {
    $ttbarj_nDzVertex_Wsel = 0;
}

$n_ppEleX_pt14_nDzVertex_Wsel = sprintf("%.4g", $ppEleX_pt14_nDzVertex_Wsel/$ppEleX_pt14_event_Wsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
if($ppEleX_pt14_nPrimaryVertices_Wsel!=0) {
    $eff_ppEleX_pt14_nDzVertex_Wsel = sprintf("%.0f", 100 * $ppEleX_pt14_nDzVertex_Wsel/$ppEleX_pt14_nPrimaryVertices_Wsel);
} else {
    $ppEleX_pt14_nDzVertex_Wsel = 0;
}
print TEXFILEW "$n_Wj_nDzVertex_Wsel ($eff_Wj_nDzVertex_Wsel \\%)&\n";
print TEXFILEW "$n_Zj_nDzVertex_Wsel ($eff_Zj_nDzVertex_Wsel \\%)&\n";
print TEXFILEW "$n_ttbarj_nDzVertex_Wsel ($eff_ttbarj_nDzVertex_Wsel \\%)&\n";
print TEXFILEW "$n_ppEleX_pt14_nDzVertex_Wsel ($eff_ppEleX_pt14_nDzVertex_Wsel \\%) \\\\ \n";


print TEXFILEZ "\$\\vert z_\{e\} - z_\{PV\}\\vert < 0.02 cm\$&\n";
$n_Wj_nDzVertex_Zsel = sprintf("%.4g", $Wj_nDzVertex_Zsel);
if($Wj_nPrimaryVertices_Zsel!=0) {
    $eff_Wj_nDzVertex_Zsel = sprintf("%.0f", 100 * $Wj_nDzVertex_Zsel/$Wj_nPrimaryVertices_Zsel);
} else {
    $eff_Wj_nDzVertex_Zsel = 0;
}

$n_Zj_nDzVertex_Zsel = sprintf("%.4g", $Zj_nDzVertex_Zsel);
if($Zj_nPrimaryVertices_Zsel!=0) {
    $eff_Zj_nDzVertex_Zsel = sprintf("%.0f", 100 * $Zj_nDzVertex_Zsel/$Zj_nPrimaryVertices_Zsel);
} else {
    $eff_Zj_nDzVertex_Zsel = 0;
}

$n_ttbarj_nDzVertex_Zsel = sprintf("%.4g", $ttbarj_nDzVertex_Zsel);
if($ttbarj_nPrimaryVertices_Zsel!=0) {
    $eff_ttbarj_nDzVertex_Zsel = sprintf("%.0f", 100 * $ttbarj_nDzVertex_Zsel/$ttbarj_nPrimaryVertices_Zsel);
} else {
    $eff_ttbarj_nDzVertex_Zsel = 0;
}

$n_ppEleX_pt14_nDzVertex_Zsel = sprintf("%.4g", $ppEleX_pt14_nDzVertex_Zsel/$ppEleX_pt14_event_Zsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
if($ppEleX_pt14_nPrimaryVertices_Zsel!=0) {
    $eff_ppEleX_pt14_nDzVertex_Zsel = sprintf("%.0f", 100 * $ppEleX_pt14_nDzVertex_Zsel/$ppEleX_pt14_nPrimaryVertices_Zsel);
} else {
    $eff_ppEleX_pt14_nDzVertex_Zsel = 0;
}
print TEXFILEZ "$n_Wj_nDzVertex_Zsel ($eff_Wj_nDzVertex_Zsel \\%)&\n";
print TEXFILEZ "$n_Zj_nDzVertex_Zsel ($eff_Zj_nDzVertex_Zsel \\%)&\n";
print TEXFILEZ "$n_ttbarj_nDzVertex_Zsel ($eff_ttbarj_nDzVertex_Zsel \\%)&\n";
print TEXFILEZ "$n_ppEleX_pt14_nDzVertex_Zsel ($eff_ppEleX_pt14_nDzVertex_Zsel \\%) \\\\ \n";




#--- dxy significance cut ---   
print TEXFILEW "\$\\vert \\delta_{xy}\\vert\/\\sigma(\\delta_{xy}) < 10 \$&\n";
$decimals = 0;
$n_Wj_nDxySigVertex_Wsel = sprintf("%.4g", $Wj_nDxySigVertex_Wsel);
if($Wj_nDzVertex_Wsel!=0) {
    $eff_Wj_nDxySigVertex_Wsel = sprintf("%.0f", 100 * $Wj_nDxySigVertex_Wsel/$Wj_nDzVertex_Wsel);
} else {
    $eff_Wj_nDxySigVertex_Wsel=0;
}

$n_Zj_nDxySigVertex_Wsel = sprintf("%.4g", $Zj_nDxySigVertex_Wsel);
if($Zj_nDzVertex_Wsel!=0) {
    $eff_Zj_nDxySigVertex_Wsel = sprintf("%.0f", 100 * $Zj_nDxySigVertex_Wsel/$Zj_nDzVertex_Wsel);
} else {
    $eff_Zj_nDxySigVertex_Wsel = 0;
}

$n_ttbarj_nDxySigVertex_Wsel = sprintf("%.4g", $ttbarj_nDxySigVertex_Wsel);
if($ttbarj_nDzVertex_Wsel!=0) {
    $eff_ttbarj_nDxySigVertex_Wsel = sprintf("%.0f", 100 * $ttbarj_nDxySigVertex_Wsel/$ttbarj_nDzVertex_Wsel);
} else {
    $eff_ttbarj_nDxySigVertex_Wsel = 0;
}

$n_ppEleX_pt14_nDxySigVertex_Wsel = sprintf("%.4g", $ppEleX_pt14_nDxySigVertex_Wsel/$ppEleX_pt14_event_Wsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
if($ppEleX_pt14_nDzVertex_Wsel!=0) {
    $eff_ppEleX_pt14_nDxySigVertex_Wsel = sprintf("%.0f", 100 * $ppEleX_pt14_nDxySigVertex_Wsel/$ppEleX_pt14_nDzVertex_Wsel);
} else {
    $eff_ppEleX_pt14_nDxySigVertex_Wsel = 0;
}
print TEXFILEW "$n_Wj_nDxySigVertex_Wsel ($eff_Wj_nDxySigVertex_Wsel \\%)&\n";
print TEXFILEW "$n_Zj_nDxySigVertex_Wsel ($eff_Zj_nDxySigVertex_Wsel \\%)&\n";
print TEXFILEW "$n_ttbarj_nDxySigVertex_Wsel ($eff_ttbarj_nDxySigVertex_Wsel \\%)&\n";
print TEXFILEW "$n_ppEleX_pt14_nDxySigVertex_Wsel ($eff_ppEleX_pt14_nDxySigVertex_Wsel \\%) \\\\ \n";


print TEXFILEZ "\$\\vert \\delta_{xy}\\vert\/\\sigma(\\delta_{xy}) < 10 \$&\n";
$n_Wj_nDxySigVertex_Zsel = sprintf("%.4g", $Wj_nDxySigVertex_Zsel);
if($Wj_nDzVertex_Zsel!=0) {
    $eff_Wj_nDxySigVertex_Zsel = sprintf("%.0f", 100 * $Wj_nDxySigVertex_Zsel/$Wj_nDzVertex_Zsel);
} else {
    $eff_Wj_nDxySigVertex_Zsel = 0;
}

$n_Zj_nDxySigVertex_Zsel = sprintf("%.4g", $Zj_nDxySigVertex_Zsel);
if($Zj_nDzVertex_Zsel!=0) {
    $eff_Zj_nDxySigVertex_Zsel = sprintf("%.0f", 100 * $Zj_nDxySigVertex_Zsel/$Zj_nDzVertex_Zsel);
} else {
    $eff_Zj_nDxySigVertex_Zsel = 0;
}

$n_ttbarj_nDxySigVertex_Zsel = sprintf("%.4g", $ttbarj_nDxySigVertex_Zsel);
if($ttbarj_nDzVertex_Zsel!=0) {
    $eff_ttbarj_nDxySigVertex_Zsel = sprintf("%.0f", 100 * $ttbarj_nDxySigVertex_Zsel/$ttbarj_nDzVertex_Zsel);
} else {
    $eff_ttbarj_nDxySigVertex_Zsel = 0;
}

$n_ppEleX_pt14_nDxySigVertex_Zsel = sprintf("%.4g", $ppEleX_pt14_nDxySigVertex_Zsel/$ppEleX_pt14_event_Zsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
if($ppEleX_pt14_nDzVertex_Zsel!=0) {
    $eff_ppEleX_pt14_nDxySigVertex_Zsel = sprintf("%.0f", 100 * $ppEleX_pt14_nDxySigVertex_Zsel/$ppEleX_pt14_nDzVertex_Zsel);
} else {
    $eff_ppEleX_pt14_nDxySigVertex_Zsel = 0;
}
print TEXFILEZ "$n_Wj_nDxySigVertex_Zsel ($eff_Wj_nDxySigVertex_Zsel \\%)&\n";
print TEXFILEZ "$n_Zj_nDxySigVertex_Zsel ($eff_Zj_nDxySigVertex_Zsel \\%)&\n";
print TEXFILEZ "$n_ttbarj_nDxySigVertex_Zsel ($eff_ttbarj_nDxySigVertex_Zsel \\%)&\n";
print TEXFILEZ "$n_ppEleX_pt14_nDxySigVertex_Zsel ($eff_ppEleX_pt14_nDxySigVertex_Zsel \\%) \\\\ \n";









#--- kinematics (MET || mee) ---   
print TEXFILEW "MET \$\>\$ 20 GeV &\n";
$decimals = 0;
$n_Wj_metCut_Wsel = sprintf("%.4g", $Wj_metCut_Wsel);
if($Wj_nDxySigVertex_Wsel!=0) {
    $eff_Wj_metCut_Wsel = sprintf("%.0f", 100 * $Wj_metCut_Wsel/$Wj_nDxySigVertex_Wsel);
} else {
    $eff_Wj_metCut_Wsel = 0;
}

$n_Zj_metCut_Wsel = sprintf("%.4g", $Zj_metCut_Wsel);
if($Zj_nDxySigVertex_Wsel!=0) {
    $eff_Zj_metCut_Wsel = sprintf("%.0f", 100 * $Zj_metCut_Wsel/$Zj_nDxySigVertex_Wsel);
} else {
    $eff_Zj_metCut_Wsel = 0;
}
$n_ttbarj_metCut_Wsel = sprintf("%.4g", $ttbarj_metCut_Wsel);
if($ttbarj_nDxySigVertex_Wsel!=0) {
    $eff_ttbarj_metCut_Wsel = sprintf("%.0f", 100 * $ttbarj_metCut_Wsel/$ttbarj_nDxySigVertex_Wsel);
} else {
    $eff_ttbarj_metCut_Wsel = 0;
}
$n_ppEleX_pt14_metCut_Wsel = sprintf("%.4g", $ppEleX_pt14_metCut_Wsel/$ppEleX_pt14_event_Wsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
if($ppEleX_pt14_nDxySigVertex_Wsel!=0) {
    $eff_ppEleX_pt14_metCut_Wsel = sprintf("%.0f", 100 * $ppEleX_pt14_metCut_Wsel/$ppEleX_pt14_nDxySigVertex_Wsel);
} else {
    $eff_ppEleX_pt14_metCut_Wsel = 0;
}

print TEXFILEW "$n_Wj_metCut_Wsel ($eff_Wj_metCut_Wsel \\%)&\n";
print TEXFILEW "$n_Zj_metCut_Wsel ($eff_Zj_metCut_Wsel \\%)&\n";
print TEXFILEW "$n_ttbarj_metCut_Wsel ($eff_ttbarj_metCut_Wsel \\%)&\n";
print TEXFILEW "$n_ppEleX_pt14_metCut_Wsel ($eff_ppEleX_pt14_metCut_Wsel \\%) \\\\ \n";

print TEXFILEW "\$W m_T \>\$ 30 GeV &\n";
$decimals = 0;
$n_Wj_transvMassCut_Wsel = sprintf("%.4g", $Wj_transvMassCut_Wsel);
if($Wj_metCut_Wsel!=0) {
    $eff_Wj_transvMassCut_Wsel = sprintf("%.0f", 100 * $Wj_transvMassCut_Wsel/$Wj_metCut_Wsel);
} else {
    $eff_Wj_transvMassCut_Wsel = 0;
}

$n_Zj_transvMassCut_Wsel = sprintf("%.4g", $Zj_transvMassCut_Wsel);
if($Zj_metCut_Wsel!=0) {
    $eff_Zj_transvMassCut_Wsel = sprintf("%.0f", 100 * $Zj_transvMassCut_Wsel/$Zj_metCut_Wsel);
} else {
    $eff_Zj_transvMassCut_Wsel = 0;
}
$n_ttbarj_transvMassCut_Wsel = sprintf("%.4g", $ttbarj_transvMassCut_Wsel);
if($ttbarj_metCut_Wsel!=0) {
    $eff_ttbarj_transvMassCut_Wsel = sprintf("%.0f", 100 * $ttbarj_transvMassCut_Wsel/$ttbarj_metCut_Wsel);
} else {
    $eff_ttbarj_transvMassCut_Wsel = 0;
}
$n_ppEleX_pt14_transvMassCut_Wsel = sprintf("%.4g", $ppEleX_pt14_transvMassCut_Wsel/$ppEleX_pt14_event_Wsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
if($ppEleX_pt14_metCut_Wsel!=0) {
    $eff_ppEleX_pt14_transvMassCut_Wsel = sprintf("%.0f", 100 * $ppEleX_pt14_transvMassCut_Wsel/$ppEleX_pt14_metCut_Wsel);
} else {
    $eff_ppEleX_pt14_transvMassCut_Wsel = 0;
}

print TEXFILEW "$n_Wj_transvMassCut_Wsel ($eff_Wj_transvMassCut_Wsel \\%)&\n";
print TEXFILEW "$n_Zj_transvMassCut_Wsel ($eff_Zj_transvMassCut_Wsel \\%)&\n";
print TEXFILEW "$n_ttbarj_transvMassCut_Wsel ($eff_ttbarj_transvMassCut_Wsel \\%)&\n";
print TEXFILEW "$n_ppEleX_pt14_transvMassCut_Wsel ($eff_ppEleX_pt14_transvMassCut_Wsel \\%) \\\\ \n";



print TEXFILEZ "70 \$\< m(e\^+e\^-)\<\$ 110 GeV/c\$\^2\$ &\n";
$decimals = 0;
$n_Wj_meeCut_Zsel = sprintf("%.4g", $Wj_meeCut_Zsel);
if($Wj_nDxySigVertex_Zsel!=0) {
    $eff_Wj_meeCut_Zsel = sprintf("%.0f", 100 * $Wj_meeCut_Zsel/$Wj_nDxySigVertex_Zsel);
} else {
    $eff_Wj_meeCut_Zsel = 0;
}

$n_Zj_meeCut_Zsel = sprintf("%.4g", $Zj_meeCut_Zsel);
if($Zj_nDxySigVertex_Zsel!=0) {
    $eff_Zj_meeCut_Zsel = sprintf("%.0f", 100 * $Zj_meeCut_Zsel/$Zj_nDxySigVertex_Zsel);
} else {
    $eff_Zj_meeCut_Zsel = 0;
}
$n_ttbarj_meeCut_Zsel = sprintf("%.4g", $ttbarj_meeCut_Zsel);
if($ttbarj_nDxySigVertex_Zsel!=0) {
    $eff_ttbarj_meeCut_Zsel = sprintf("%.0f", 100 * $ttbarj_meeCut_Zsel/$ttbarj_nDxySigVertex_Zsel);
} else {
    $eff_ttbarj_meeCut_Zsel = 0;
}
$n_ppEleX_pt14_meeCut_Zsel = sprintf("%.4g", $ppEleX_pt14_meeCut_Zsel/$ppEleX_pt14_event_Zsel * $Xs_ppEleX_pt14 * $filtereff_ppEleX_pt14 * $lumi );
if($ppEleX_pt14_nDxySigVertex_Zsel!=0) {
    $eff_ppEleX_pt14_meeCut_Zsel = sprintf("%.0f", 100 * $ppEleX_pt14_meeCut_Zsel/$ppEleX_pt14_nDxySigVertex_Zsel);
} else {
    $eff_ppEleX_pt14_meeCut_Zsel = 0;
}

print TEXFILEZ "$n_Wj_meeCut_Zsel ($eff_Wj_meeCut_Zsel \\%)&\n";
print TEXFILEZ "$n_Zj_meeCut_Zsel ($eff_Zj_meeCut_Zsel \\%)&\n";
print TEXFILEZ "$n_ttbarj_meeCut_Zsel ($eff_ttbarj_meeCut_Zsel \\%)&\n";
print TEXFILEZ "$n_ppEleX_pt14_meeCut_Zsel ($eff_ppEleX_pt14_meeCut_Zsel \\%) \\\\ \n";


# # final tables

print TEXFILEW "\\hline \n";
print TEXFILEW "\\hline \n";
print TEXFILEW "\\end\{tabular\}\n";
print TEXFILEW "\\caption\{\\emph Number of expected events for an integrated luminosity of $lumi pb\$^\{-1\}\$ after each cut for \$W(e \\nu) + jets\$ selection,
for the processes \$W+jets\$, \$Z+jets\$, \$t \\bar\{t\}+jets\$ and \$QCD\$. The relative efficiencies with respect to the previous cut are given in percent within the brackets.\}\n";
print TEXFILEW "\\label\{tab:effWsel\}\n";
print TEXFILEW "\\end\{center\}\n";
print TEXFILEW "\\end\{sidewaystable\}\n";

print TEXFILEZ "\\hline \n";
print TEXFILEZ "\\hline \n";
print TEXFILEZ "\\end\{tabular\}\n";
print TEXFILEZ "\\caption\{\\emph Number of expected events for an integrated luminosity of $lumi pb\$^\{-1\}\$ after each cut for \$Z(e\^+e\^-) + jets\$ selection,
for the processes \$W+jets\$, \$Z+jets\$, \$t \\bar\{t\}+jets\$ and \$QCD\$. The relative efficiencies with respect to the previous cut are given in percent within the brackets.\}\n";
print TEXFILEZ "\\label\{tab:effWsel\}\n";
print TEXFILEZ "\\end\{center\}\n";
print TEXFILEZ "\\end\{sidewaystable\}\n";


print "Done. Include in your tex and compile.\n";


