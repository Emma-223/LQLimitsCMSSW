#!/usr/bin/env python3
import os
import sys
from optparse import OptionParser
import subprocess
import shlex
import glob
from termcolor import colored
import time
import datetime
import numpy as np
import re
import math
import multiprocessing
import traceback
import json
import shutil
from rich.progress import Progress
from functools import partial
from pathlib import Path
from bisect import bisect
from tabulate import tabulate
from ROOT import TFile, TGraph, TSpline3

from BR_Sigma_EE_vsMass import BR_Sigma_EE_vsMass
from ComboPlotLQ1 import ComboPlot

from runLimits import (
    ErrorCallback,
    GetExecTimeStr,
    RunCommand,
#    GetEOSPathWithFullUrl,
    ConvertDatacardToWorkspace,
    FindCardWorkspace,
#    GenerateAsimovToyData,
#    FindAsimovToyData,
    FindFile,
    FindFileMultiple,
    GetFileList,
    GetFileListMultiple,
    GetFilenamesFromCondorOutFiles,
    MakeCondorDirs,
    ExtractLimitResult,
    ExtractAsymptoticLimitResultTxt,
    ReadXSecFile,
    CheckErrorFile,
    CheckForErrorAndResultFiles,
    ComputeBetaLimits,
    CreateComboArrays,
    get_value,
)

if not 'LQANA' in os.environ:
    raise RuntimeError('Need to define $LQANA environment variable to point to your git clone of the rootNtupleAnalyzerV2 repo.')
sys.path.append(os.getenv("LQANA").rstrip("/") + "/scripts/")
from combineCommon import SeparateDatacards, GetYearAndIntLumiFromDatacard
if not 'CMSSW_VERSION' in os.environ:
    raise RuntimeError('You must be in a CMSSW environment (release area + cmsenv)')
ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')

#NOTE: script was designed to work from CMS connect. Running on lxplus will likely not work without further development.
#On CMS connect, to use CMSSW you need to run in a container by running the command: cmssw-el9
#however, you can't submit to condor from within the container, so submitting the limit jobs runs in two steps.
#(1) write the condor.sh and condor.sub files (from inside the container)
#(2) loop over all directories and submit the jobs (from outside the container)

def GenerateAsimovToyDataBetaId(workspace,mass,toysDir,dictAsimovToysByBetaId, betaId, scaleFactor=1.0):
    cmd = 'combine -M GenerateOnly {} -t -1 --seed -1 --saveToys -n .asimov.betaId{} -m {} {}'.format(workspace, betaId, mass, commonCombineArgs.format(scaleFactor))
    RunCommand(cmd, toysDir, None, True)
    toyFile = Path(sorted(glob.glob(toysDir+'/higgsCombine.asimov.betaId{}.GenerateOnly.mH{}.*.root'.format(betaId, mass)), key=os.path.getmtime)[-1]).resolve()
    dictAsimovToysByBetaId[betaId] = toyFile
    return toyFile

def FindAsimovToyDataBetaId(mass, betaId, toysDir):
    listOfToyFiles = sorted(glob.glob(toysDir+'/higgsCombine.asimov.betaId{}.GenerateOnly.mH{}.*.root'.format(betaId, mass)), key=os.path.getmtime)
    if len(listOfToyFiles) > 0:
        return Path(listOfToyFiles[-1]).resolve()

def GetRMinAndRMax(mass, quant, betaId):
    rMed = rValuesByMassBetaAndQuantile[str(mass)][str(betaId)][str(quant)]
    if quant == 0.025 and mass > 800:
        if mass < 1200:
            rMax = 1.8
            rMin = 0.85
        else:
            rMax = 2.5
            rMin = 1.25
    elif quant == 0.16 and mass > 1000:
        rMax = 2.0
        rMin = 1.0
    elif quant == 0.975 and mass >= 1200:
        rMax = 0.9
        rMin = 0.3
    else:
        rMax = 1.3
        rMin = 0.75
    rMax *= rMed
    rMin *= rMed
    return rMin, rMax

def WriteCondorSubFile(condorDir,inputFiles,nJobs):
    filename = condorDir+"/condor.sub"
    with open(filename,'w') as f:
        f.write("Universe = vanilla\n")
        f.write("executable = "+filename.replace(".sub",".sh")+"\n")
        f.write("arguments = $(procId)\n")
        f.write("output                = "+filename.replace(".sub",".$(ClusterId).$(ProcId).out")+"\n")
        f.write("error                 = "+filename.replace(".sub",".$(ClusterId).$(ProcId).err")+"\n")
        f.write("log                   = "+filename.replace(".sub",".$(ClusterId).$(ProcId).log")+"\n")
        f.write("+REQUIRED_OS = \"rhel9\"\n\n")
        combineBin = shutil.which("combine")
        if combineBin is None:
            raise RuntimeError("No combine binary found in $PATH")
        cmssw_setup_script = "/home/emma-pearson/cmssw_setup_copyForDebug.sh" #Come up with a better way of doing this than hardcoding it
        inputFiles += [combineBin, sandbox, cmssw_setup_script]
        f.write("transfer_input_files = {}\n".format(",".join(inputFiles)))
        f.write("should_transfer_files = YES\n\n")
        f.write("max_materialize = 50\n")
        f.write("max_retries = 10\n")
        f.write("on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n")
        f.write("periodic_hold = (JobStatus == 2) && ((CurrentTime - EnteredCurrentStatus) > 24 * 60 * 60)\n")#runs for > 24hrs
        f.write("periodic_hold_reason = \"ran for more than 24 hrs\"\n")
        f.write("periodic_release =  (NumJobStarts < 10) && ((CurrentTime - EnteredCurrentStatus) > 300)\n\n")
        f.write("queue {}\n".format(nJobs))

def GetCommonCondorShFileLines(filesToMove):
    lines = []
    lines.append("#!/bin/sh\n")
    lines.append("unlimit -s unlimited\n")
    lines.append("set -e\n")
    lines.append("condorNodeBaseDir=$PWD\n")
    lines.append("source ./cmssw_setup_copyForDebug.sh\n")
    lines.append("export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch\n")
    lines.append("source $VO_CMS_SW_DIR/cmsset_default.sh\n")
    sandboxName = sandbox.split("/")[-1]
    lines.append("cmssw_setup {}\n".format(sandboxName))
    lines.append("cd $CMSSW_BASE\n")
    copyCmd = "cp"
    for f in filesToMove:
        copyCmd += " $condorNodeBaseDir/"+f
    copyCmd += " ."
    lines.append(copyCmd+"\n")
    return lines

def WriteCondorFilesLimits(cardWorkspace,mass,quant,betaList,condorDir,asimovToys,gridFiles):
    eosPath = eosDir + "/" + condorDir
    inputFiles = [cardWorkspace]+asimovToys
    nJobs = len(betaList)
    WriteCondorSubFile(condorDir, inputFiles, nJobs)
    print("INFO: wrote file {}".format(condorDir+"/condor.sub"))
    if not os.path.isfile(condorDir+"/condor.sub"):
        raise RuntimeError("something went wrong while writing condor.sub file for mass {}, quant {}, beta ids {} to {}".format(mass, quant, betaList[0], betaList[-1]))
    filesToMove = []
    for f in [cardWorkspace]+asimovToys:
        filesToMove.append(f.split("/")[-1])
    shFileCommonLines = GetCommonCondorShFileLines(filesToMove)
    shFileName = condorDir+"/condor.sh"
    with open(shFileName, 'w') as f:
        f.writelines(shFileCommonLines)
        f.write("\n")
        for i, betaId in enumerate(betaList):
            betaVal = betasToScan[betaId]
            if not options.disableSigRescaling:
                scaleFactor = sigRescaleByMassBetaAndQuantile[str(mass)][str(betaIdx)]['0.5'] 
            else:
                scaleFactor = 1.0
            scaleFactor *= betaVal**2
            rMin, rMax = GetRMinAndRMax(mass, quant,betaId)
            f.write("if [ $1 -eq {} ]; then\n".format(i))
            f.write("  sigSF={}\n".format(scaleFactor))
            f.write("  betaId={}\n".format(betaId))
            f.write("  rMin={}\n".format(rMin))
            f.write("  rMax={}\n".format(rMax))
            f.write("  toyFile={}\n".format(asimovToys[i].split("/")[-1]))
            f.write("  plotFile=limit_scan_m{}_quant{}_betaId{}.pdf\n".format(mass,str(quant).replace(".","p"),betaId))
            f.write("  combinedGridFile=higgsCombine.betaId{}.HybridNew.mH{}.{}.gridAll.root\n".format(betaId,mass,str(quant).replace(".","p")))
            gridFileListStr = " ".join(gridFiles[str(betaId)])
            f.write("  gridFiles=\"{}\"\n".format(gridFileListStr))
            f.write("fi\n\n")
        f.write("for gridFile in $gridFiles; do\n")
        f.write("  echo \"root://eoscms.cern.ch/$gridFile\"\n")
        f.write("  xrdcp -fs root://eoscms.cern.ch/$gridFile .\n")
        f.write("done\n")
        f.write("ls\n")
        f.write("hadd -fk207 $combinedGridFile higgsCombine.*.grid*.root\n\n")
        cmd = "combine {} -M HybridNew --LHCmode LHC-limits --rMin $rMin --rMax $rMax --grid=$combinedGridFile --rAbsAcc 0.00001 --rRelAcc 0.005 --plot $plotFile --readHybridResults -m {} -n .betaId$betaId".format(str(cardWorkspace).split("/")[-1],mass)
        if quant > 0:
            cmd += " --expectedFromGrid={}".format(quant)
        cmd += " -D $toyFile:toys/toy_asimov"
        f.write(cmd+"\n\n")
        f.write("ls higgsCombine*\n")
        f.write("limitOutputFile=higgsCombine.betaId$betaId.HybridNew.mH{}.quant{:.3f}.root\n".format(mass,quant))
        f.write("xrdcp -fs $plotFile {}/$plotFile\n".format(eosPath))
        f.write("xrdcp -fs $limitOutputFile {}/$limitOutputFile\n".format(eosPath))
        f.write("xrdcp -fs $combinedGridFile {}/$combinedGridFile\n".format(eosPath))
        f.write("limitFileEOS=$(xrdfs root://eoscms.cern.ch/ ls /{}/$limitOutputFile)\n".format(eosPath.split("//")[-1]))
        f.write("echo \"limit file: $limitFileEOS\"\n")

def WriteCondorFilesGridGen(cardWorkspace,mass,quant,betaList,condorDir,asimovToys):
    eosPath = eosDir + "/" + condorDir
    inputFiles = [cardWorkspace] + asimovToys
    nJobs = len(betaList)* int(50/nGridPointsPerSubFile)
    WriteCondorSubFile(condorDir, inputFiles, nJobs)
    filesToMove = []
    for f in [cardWorkspace]+asimovToys:
        filesToMove.append(f.split("/")[-1])
    shFileCommonLines = GetCommonCondorShFileLines(filesToMove)
    shFileName = condorDir+"/condor.sh"
    if quant == 0.025:
        nToys = 40000
    elif quant == 0.16:
        nToys = 30000
    else:
        nToys = 10000
    with open(shFileName,'w') as f:
        f.writelines(shFileCommonLines)
        f.write("\n")
        for i, betaIdx in enumerate(betaList):
            if not options.disableSigRescaling:
                scaleFactor = sigRescaleByMassBetaAndQuantile[str(mass)][str(betaIdx)]['0.5']
            else:
                scaleFactor = 1.0
            betaVal = betasToScan[betaIdx]
            scaleFactor *= betaVal**2
            jobIdMin = int(50/nGridPointsPerSubFile) * i
            jobIdMax = jobIdMin + int(50/nGridPointsPerSubFile) - 1
            rMin, rMax = GetRMinAndRMax(mass, quant,betaIdx)
         #   stepSize = (rMax - rMin) / nGridPointsPerSubFile
            f.write("if [[ $1 -ge {} ]] && [[ $1 -le {} ]]; then\n".format(jobIdMin, jobIdMax))
            f.write("  sigSF={}\n".format(scaleFactor))
            f.write("  betaId={}\n".format(betaIdx))
            f.write("  rMin={}\n".format(rMin))
            f.write("  rMax={}\n".format(rMax))
            f.write("  toyFile={}\n".format(asimovToys[i].split("/")[-1]))
            rVals = np.linspace(rMin,rMax,num=50)
            rValsList = []
            for r in rVals:
                rValsList.append(str(r))
            for i in range(int(50/nGridPointsPerSubFile)):
                f.write("  if [[ $1 -eq {} ]]; then\n".format(jobIdMin+i))
                rValsThisJob = rValsList[i*nGridPointsPerSubFile:(i+1)*nGridPointsPerSubFile]
                f.write("    rVals=\'{}\'\n".format(" ".join(rValsThisJob)))
                f.write("  fi\n")
            f.write("fi\n\n")

        #asimovGenCmd = "combine -M GenerateOnly {} -t -1 --seed -1 --saveToys -n .asimov.betaId$betaId -m {} ".format(cardWorkspace.split("/")[-1], mass) + commonCombineArgs.format("$sigSF")

        cmd = "combine -d {} -v1 -M HybridNew --LHCmode LHC-limits --saveHybridResult --saveToys --seed -1 --clsAcc 0 --rMin $rMin --rMax $rMax -T {} -m {} ".format(cardWorkspace.split("/")[-1], nToys, mass)
        if quant > 0:
            cmd += "  --expectedFromGrid {}".format(quant)
        cmd += " -D $toyFile:toys/toy_asimov"
        cmd += " "+commonCombineArgs.format("$sigSF")
        #f.write(asimovGenCmd+'\n')
        #f.write("toyFile=$(ls higgsCombine.asimov*.root)\n")
        #f.write("xrdcp -fs $toyFile {}/$toyFile\n".format(eosPath))
        #f.write("toyFileEOS=$(xrdfs root://eoscms.cern.ch ls /{}/$toyFile)\n".format(eosPath.split("//")[-1]))
        #f.write("echo \"asimov toy file: $toyFileEOS\"\n\n")

        f.write("for rval in $rVals; do\n")
        f.write("echo \"rValue $rval\"\n")
        f.write(cmd + " -n .betaId$betaId.POINT.$rval --singlePoint $rval\n")
        f.write("done\n\n")
        combinedOutputFile = "higgsCombine.betaId$betaId.HybridNew.mH{}.{}.grid$1.root".format(mass,str(quant).replace(".","p"))
        f.write("ls\n")
        f.write("hadd -fk207 {} higgsCombine.*POINT*.root\n".format(combinedOutputFile))
        f.write("xrdcp -fs {} {}/{}\n".format(combinedOutputFile, eosPath, combinedOutputFile))
        f.write("gridFileEOS=$(xrdfs root://eoscms.cern.ch ls /{}/{})\n".format(eosPath.split("//")[-1],combinedOutputFile))
        f.write("echo \"grid file: $gridFileEOS\"\n")

def WriteCondorFilesAsymptotic(cardWorkspace, betaList, mass, condorDir):
    inputFiles = [cardWorkspace]
    nJobs = len(betaList)
    WriteCondorSubFile(condorDir, inputFiles, nJobs)
    shFileCommonLines = GetCommonCondorShFileLines([cardWorkspace.split("/")[-1]])
    shFileName = condorDir+"/condor.sh"
    eosPath = eosDir + "/"+condorDir
    with open(shFileName,'w') as f:
        f.writelines(shFileCommonLines)
        f.write("\n")
        for i,betaIdx in enumerate(betaList):
            if not options.disableSigRescaling:
                scaleFactor = sigRescaleByMassBetaAndQuantile[str(mass)][str(betaIdx)]['0.5']
            else:
                scaleFactor = 1.0
            betaVal = betasToScan[betaIdx]
            cmd = "combine {} -v1 -M AsymptoticLimits --seed -1 --strict --rMax 200 --rMin 0 -m {} --rAbsAcc 0.0001 -n .betaId{}".format(cardWorkspace.split("/")[-1], mass, betaIdx)
            cmd += " "+commonCombineArgs.format(scaleFactor * betaVal**2)
            if blinded:
                cmd += " --run blind --expectSignal 0"
            f.write("if [ $1 -eq {} ]; then\n".format(i))
            f.write(cmd+"\n")
            f.write("fi\n\n")
        f.write("outputRootFile=$(ls higgsCombine*.root)\n")
        f.write("xrdcp -fs $outputRootFile {}/$outputRootFile\n".format(eosPath))
        f.write("outputEOSFile=$(xrdfs root://eoscms.cern.ch ls /{}/$outputRootFile)\n".format(eosPath.split("//")[-1]))
        f.write('echo "output eos file: $outputEOSFile "\n')

def SplitBetasBetweenCondorSubmissions(betas, nBetasPerFile):
    nListsToMake = int(math.ceil(len(betas)/nBetasPerFile))
    listOfBetaLists = []
    for i in range(nListsToMake):
        startIdx = i*nBetasPerFile
        if i==nListsToMake-1:
            betaList = betas[startIdx:]
        else:
            betaList = betas[startIdx:startIdx+nBetasPerFile]
        listOfBetaLists.append(betaList)
    return listOfBetaLists

def GetBetaRangeToAttempt(mass):
    minBetasPerMass = {300: 0,   400: 0,   500: 0, 600: 0, 700: 0, 800: 0, 900: 0, 1000: 0,  1100: 0.085, 1200: 0.095, 1300: 0.105, 1400: 0.115, 1500: 0.125, 1600: 0.125, 1700: 0.5, 1800: 0.7}
    maxBetasPerMass = {300: 0.075, 400: 0.1, 500: 0.15, 600: 0.2, 700: 0.25, 800: 0.35, 900: 0.45, 1000: 0.6, 1100: 0.7, 1200: 0.8, 1300: 0.9, 1400: 1.0, 1500: 1.0, 1600: 1.0, 1700: 1.0, 1800: 1.0}
    if mass in minBetasPerMass.keys():
        return minBetasPerMass[mass], maxBetasPerMass[mass]
    else:
        return 0.8,1.0
    #if mass > 1800:
    #    return 0.8,1.0
    #else:
    #    return minBetasPerMass[mass], maxBetasPerMass[mass]

def GetBetasToSubmit(mass):
    minBeta, maxBeta = GetBetaRangeToAttempt(int(mass))
    minId = bisect(betasToScan, minBeta)
    maxId = bisect(betasToScan, maxBeta)
    idList = list(range(minId, maxId))
    return idList

def CheckCondorJobsAndGetFilenames():
    foundAllFiles = False
    allFilesAreGood = False
    filesDict = {}
    missingFilesList = []
    badOutFiles = []
    with Progress() as progress:
        nFilesExpected = len(massList)*len(quantilesExpected)
        task_id = progress.add_task("[cyan]Checking output of condor jobs",total=nFilesExpected)
        for mass in massList:
            filesDict[str(mass)] = {}
            betas = GetBetasToSubmit(mass)
            betaLists = SplitBetasBetweenCondorSubmissions(betas, nBetasPerCondorSub)
            for quant in quantilesExpected:
                filesDict[str(mass)][str(quant)] = {}
                #print("INFO: check grid gen output files for mass {}, quantile {}".format(mass,quant))
                for betaList in betaLists:
                    jobId = 0
                    for b in betaList:
                        filesToCheck = []
                        if options.submitLimitJobsAfterGridGen:#if we're writing files for the limit jobs, we need to check the output of the gridGen jobs which have multiple out files in each directory per beta
                            for j in range(int(nGridPoints/nGridPointsPerSubFile)):
                                f = dirName+"/gridGen/hybridNewGridGen.M{}.{}/betaIds{}To{}/condor.*.{}.out".format(mass,str(quant).replace(".","p"),betaList[0],betaList[-1],jobId)
                                filesToCheck.append(f)
                                jobId+=1
                        else: #check limit jobs, 1 out file per directory per beta
                            filesToCheck = [dirName+"/limits/hybridNewLimits.M{}.{}/betaIds{}To{}/condor.*.{}.out".format(mass,str(quant).replace(".","p"),betaList[0],betaList[-1],jobId)]
                            jobId+=1
                        filesThisBeta = []
                        for f in filesToCheck:
                            try:
                                condorFile = FindFile(f)
                            except RuntimeError:
                                missingFilesList.append(f)
                                continue
                            condorFile = str(condorFile)
                            with open(condorFile,'r') as outFile:
                                lines = outFile.readlines()
                            if not len(lines)>0:
                                badOutFiles.append(f)
                                continue
                            lastLine = lines[-1]
                            if not ("grid file" in lastLine or "limit file" in lastLine):
                                badOutFiles.append(f)
                                continue
                            gridFile = lastLine.split()[-1].strip()
                            filesThisBeta.append(gridFile)
                        filesDict[str(mass)][str(quant)][str(b)] = filesThisBeta
                progress.advance(task_id)
    if options.submitLimitJobsAfterGridGen:
        outTxtFileName = dirName+"/incompleteGridJenJobs.txt"
    else:
        outTxtFileName = dirName+"/incompleteLimitJobs.txt"
    with open(outTxtFileName,'w') as f:
        if len(missingFilesList):
            f.write("# The following are missing condor.out files:\n")
            for missingFile in missingFilesList:
                f.write(missingFile+"\n")
        else:
            foundAllFiles = True
        if len(badOutFiles):
            f.write("# The following condor.out files do not have a grid root file as their last line\n")
            for badFile in badOutFiles:
                f.write(badFile+"\n")
        else:
            allFilesAreGood = True
    if not (allFilesAreGood and foundAllFiles):
        print("missing and incomplete condor.out files have been written to {}".format(outTxtFileName))
        print("run 'python3 scripts/submitJobsTxtFile.py {}' to resubmit these jobs to condor".format(outTxtFileName))

    return filesDict, (allFilesAreGood and foundAllFiles) 
                    
def ReadLimitJobResults(mass,quant,beta,rootFile):
    betaVal = betasToScan[beta]
    if options.disableSigRescaling:
        sigSF = 1.0
    else:
        sigSF = sigRescaleByMassBetaAndQuantile[str(mass)][str(betaIdx)]['0.5']
    limit, limitErr, qFromFile, sigSFFromFile = ExtractLimitResult(rootFile)
    rLimitsByMassBetaAndQuantile[str(quant)][str(beta)][str(mass)] = limit
    xsecLimitsByMassBetaAndQuantile[str(quant)][str(beta)][str(mass)] = limit * sigSF * xsThByMass[float(mass)]


#################################################################################
# Run
#################################################################################
if __name__ == "__main__":
    massList = list(range(300,2100,100))
    #massList = [1800]
    betasToScan = list(np.linspace(0.0, 1, 500))[:-1] + [0.9995]
    nBetasPerCondorSub = 10
    nGridPoints = 50
    nGridPointsPerSubFile = 5
    eosDir = "root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/eipearso"
    eosDirNoPrefix = eosDir[eosDir.rfind("//")+1:]
    quantilesExpected = [0.025, 0.16, 0.5, 0.84, 0.975]
    #quantilesExpected = [0.025]
    xsThFilename = "$LQANA/config/xsection_theory_13TeV_scalarPairLQ.txt"
    sigRescaleFile = "rValues_nominal_betaScan.json"
    commonCombineArgs = " --setParameters signalScaleParam={} --freezeParameters signalScaleParam --trackParameters signalScaleParam"
    blinded = True
    sandbox = os.getenv("HOME")+"/sandbox-CMSSW_14_1_0_pre4-8e58ceb.tar.bz2"
    ncores = 6
    parser = OptionParser(
        usage="%prog -d datacard -n name [options]"
    )
    parser.add_option(
        "-d",
        "--datacard",
        dest="datacard",
        help="combined datacard",
        metavar="datacard",
    )
    parser.add_option(
        "-n",
        "--name",
        dest="name",
        help="name of dir where results will be stored",
        metavar="dame",
    )
    parser.add_option(
        "-e",
        "--estimateRValueScanRange",
        dest="estimateRValueScanRange",
        help="run asymptotic limits to use for estimating r value scan ranges",
        metavar="estimateRValueScanRange",
        action="store_true",
        default=False,
    )
    parser.add_option(
        "-r",
        "--readResults",
        dest="readResults",
        help="read limit results. Pair with -e option to read asymptotic limit results",
        metavar="readResults",
        action="store_true",
        default=False,
    )
    parser.add_option(
        "-g",
        "--generateGrids",
        dest="generateGrids",
        help="generate CLs vs r grids for limit calculation",
        metavar="generateGrids",
        action="store_true",
        default=False,
    )
    parser.add_option(
        "-f",
        "--fileWithRValueRanges",
        dest="fileWithRValueRanges",
        help="file with the r values output by the asymptotic limit calculation",
        metavar="fileWithRValueRanges",
        default=None
    )
    parser.add_option(
        "-l",
        "--submitLimitJobsAfterGridGen",
        dest="submitLimitJobsAfterGridGen",
        help="write condor files for limit calcuation from grids",
        metavar="submitLimitJobsAfterGridGen",
        action="store_true",
        default=False,
    )
    parser.add_option(
        "--noSigSF",
        dest="disableSigRescaling",
        help="disable rescaling of signal yields to set r close to 1",
        metavar="disableSigRescaling",
        action="store_true",
        default=False,
    )
    (options, args) = parser.parse_args()
    combinedDatacard = options.datacard
    manager = multiprocessing.Manager()
    dictAsimovToysByBetaId = manager.dict()
    dirName = options.name
    eosDirName = eosDir+"/"+options.name
    eosDirNameNoPrefix = eosDirNoPrefix+"/"+options.name
    rValuesByMassBetaAndQuantile = {}
    xsThByMass, _, _ = ReadXSecFile(xsThFilename)

    if options.generateGrids or options.submitLimitJobsAfterGridGen:
        if options.fileWithRValueRanges is None:
            raise RuntimeError("Can't run grid gen or limits without approximate r values to estimate the scan ranges from")
        with open(options.fileWithRValueRanges) as f:
            rValuesByMassBetaAndQuantile = json.load(f)
    if options.submitLimitJobsAfterGridGen:
        if os.path.isfile(dirName+"/gridGenFilenames.json"):
            with open(dirName+"/gridGenFilenames.json",'r') as f:
                gridGenFilenamesByMassBetaAndQuantile = json.load(f)
        else:
            gridGenFilenamesByMassBetaAndQuantile, goodFiles = CheckCondorJobsAndGetFilenames()
            if not goodFiles:
                raise RuntimeError("Did not find all grid files needed for limit jobs")
            with open(dirName+"/gridGenFilenames.json",'w') as f:
                json.dump(gridGenFilenamesByMassBetaAndQuantile,f)
    if options.readResults and not options.estimateRValueScanRange:
        if os.path.isfile(dirName+"/limitFilenames.json"):
            with open(dirName+"/limitFilenames.json",'r') as f:
                limitFileNamesByMassBetaAndQuantile = json.load(f)
        else:
            limitFileNamesByMassBetaAndQuantile, goodFiles = CheckCondorJobsAndGetFilenames()
            if not goodFiles:
                raise RuntimeError("Did not find all limit job output files")
            with open(dirName+"/limitFilenames.json",'w') as f:
                json.dump(limitFileNamesByMassBetaAndQuantile,f)
                print("limit files written to "+dirName+"/limitFilenames.json")
                exit()
    if not options.disableSigRescaling:
        with open(sigRescaleFile,'r') as f:
            sigRescaleByMassBetaAndQuantile = json.load(f)

    if options.readResults:
        if options.estimateRValueScanRange:
            for mass in massList:
                betasToSubmit = GetBetasToSubmit(mass)
                betaLists = SplitBetasBetweenCondorSubmissions(betasToSubmit,nBetasPerCondorSub)
                rValuesByMassBetaAndQuantile[mass] = {}
                fileList = []
                errFiles = []
                for betaList in betaLists:
                    limitDirName = dirName+"/asymptoticLimits/asymptotic.M{}/betaIds{}To{}".format(mass,betaList[0],betaList[-1])
                    globString = limitDirName+"/condor.*.err"
                    errFiles = GetFileList(globString)
                    for f in errFiles:
                        CheckErrorFile(f,False)
                        outFile = str(FindFile(f.replace(".err",".out")))
                        jobId = outFile.split("/")[-1].split(".")[2]
                        betaId = betaList[int(jobId)]
                        limits, limitErrs, quantiles, signalScaleParam = ExtractAsymptoticLimitResultTxt(outFile)
                        rValuesByMassBetaAndQuantile[mass][betaId] = dict(zip(quantiles,limits))
            with open(dirName+"/asymptoticLimits/rValues.json","w") as f:
                json.dump(rValuesByMassBetaAndQuantile,f)
            exit()
        else:
            start = time.time()
            needToReadRootFiles = True
            if os.path.isfile(dirName+"/xsecLimits.json"):
                with open(dirName+"/xsecLimits.json",'r') as f:
                    xsecLimitsByMassBetaAndQuantile = json.load(f)
                needToReadRootFiles = False
            if os.path.isfile(dirName+"/rLimits.json"):
                with open(dirName+"/rLimits.json",'r') as f:
                    rLimitsByMassBetaAndQuantile = json.load(f)
            if needToReadRootFiles:
                rLimitsByMassBetaAndQuantile = manager.dict()
                xsecLimitsByMassBetaAndQuantile = manager.dict()
                for quant in quantilesExpected:
                    rLimitsByMassBetaAndQuantile[str(quant)] = manager.dict()
                    xsecLimitsByMassBetaAndQuantile[str(quant)] = manager.dict()
                for mass in massList:
                    betasToSubmit = GetBetasToSubmit(mass)
                    betaLists = SplitBetasBetweenCondorSubmissions(betasToSubmit,nBetasPerCondorSub)
                    for quant in quantilesExpected:
                        with Progress() as progress:
                            task_id = progress.add_task("[cyan]Reading batch results for mass {}, quantile {}".format(mass,quant),total=len(betasToSubmit))
                            with multiprocessing.Pool(ncores,maxtasksperchild=1) as pool:
                                partial_quit = partial(ErrorCallback,pool)
                                for beta in betasToSubmit:
                                    if not str(beta) in rLimitsByMassBetaAndQuantile[str(quant)].keys():
                                        rLimitsByMassBetaAndQuantile[str(quant)][str(beta)] = manager.dict()
                                    if not str(beta) in xsecLimitsByMassBetaAndQuantile[str(quant)].keys():
                                        xsecLimitsByMassBetaAndQuantile[str(quant)][str(beta)] = manager.dict()
                                    try:
                                        rootFile = limitFileNamesByMassBetaAndQuantile[str(mass)][str(quant)][str(beta)]
                                        rootFile = "root://eoscms.cern.ch/"+rootFile[0]
                                        #rootFile = "root://cmseos.fnal.gov/"+rootFile[0].replace("/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/","/store/user/")
                                        args = [mass,quant,beta,rootFile]
                                        pool.apply_async(ReadLimitJobResults,args, callback = lambda x: progress.advance(task_id), error_callback = partial_quit)
                                    except Exception as e:
                                        print("ERROR: caught exception reading limit job results for mass {}, quantile {}, beta id {}".format(mass, quantile, beta))
                                pool.close()
                                pool.join()
                with open(dirName+"/rLimits.json",'w') as f:
                    json.dump(get_value(rLimitsByMassBetaAndQuantile),f)
                with open(dirName+"/xsecLimits.json","w") as f:
                    json.dump(get_value(xsecLimitsByMassBetaAndQuantile),f)
            if not os.path.isfile(dirName+"/massLimits.json"):
                massLimits = ComputeBetaLimits(xsThByMass, xsecLimitsByMassBetaAndQuantile)
                with open(dirName+"/massLimits.json","w") as f:
                    json.dump(get_value(massLimits),f)
            else:
                with open(dirName+"/massLimits.json","r") as f:
                    massLimits = json.load(f)
            comboArraysByQuantile = CreateComboArrays(massLimits)
            with open(dirName+"/comboArrays.json","w") as f:
                json.dump(get_value(comboArraysByQuantile),f)
            end = time.time()
            tot = end - start
            print("Time to read results: {}".format(tot))
            years, intLumi = GetYearAndIntLumiFromDatacard(combinedDatacard)
            intLumi *= 0.001
            ComboPlot(dirName, intLumi, comboArraysByQuantile)
            exit()

    separateDatacardsDir = dirName+"/datacards"
    if not os.path.isdir(dirName):
        print("INFO: Making directory", dirName, flush=True)
        Path(dirName).mkdir(exist_ok=True)
    if not os.path.isdir(separateDatacardsDir):
        print("INFO: Making directory", separateDatacardsDir, flush=True)
        Path(separateDatacardsDir).mkdir(exist_ok=True)
    massListFromCards, cardFilesByMass, _ = SeparateDatacards(combinedDatacard, 0, separateDatacardsDir)
    for mass in massList:
        cardFile = cardFilesByMass[mass]
        cardWorkspace = FindCardWorkspace(cardFile, mass)
        if cardWorkspace is not None:
            print("INFO: Using previously-generated card workspace: {}".format(cardWorkspace), flush=True)
        else:
            cardWorkspace = ConvertDatacardToWorkspace(cardFile, mass)

        if options.estimateRValueScanRange and not options.readResults:
            #submit asymptotic limit jobs
            betasToSubmit = GetBetasToSubmit(mass)
            betaLists = SplitBetasBetweenCondorSubmissions(betasToSubmit,nBetasPerCondorSub)
            print("INFO: for mass {}, use beta lists".format(mass), betaLists)
            nCondorFiles = len(betaLists)
            if not os.path.isdir(dirName+"/asymptoticLimits"):
                print("INFO: Making directory", dirName+"/asymptoticLimits", flush=True)
                Path(dirName+"/asymptoticLimits").mkdir(exist_ok=True)
            if not os.path.isdir(dirName+"/asymptoticLimits/asymptotic.M{}".format(mass)):
                print("INFO: Making directory", dirName+"/asymptoticLimits/asymptotic.M{}".format(mass), flush=True)
                Path(dirName+"/asymptoticLimits/asymptotic.M{}".format(mass)).mkdir(exist_ok=True)
            with Progress() as progress:
                task_id = progress.add_task("[cyan]Writing {} condor.sub and .sh files...".format(nCondorFiles), total=nCondorFiles)
                with multiprocessing.Pool(ncores) as pool:
                    partial_quit = partial(ErrorCallback,pool)
                    for betaList in betaLists:
                        try:
                            condorDir = dirName+"/asymptoticLimits/asymptotic.M{}/betaIds{}To{}".format(mass,betaList[0],betaList[-1])
                            if not os.path.isdir(condorDir):
                                Path(condorDir).mkdir(exist_ok=True)
                            args = [str(cardWorkspace), betaList, mass, condorDir]
                            pool.apply_async(WriteCondorFilesAsymptotic,args,callback = lambda x: progress.advance(task_id), error_callback = partial_quit)
                        except Exception as e:
                            print("ERROR: caught exception while writing files for mass {}, betas {} to {}".format(mass, betaList[0], betaList[-1]))
                    pool.close()
                    pool.join()
        else:
            if not os.path.isdir(dirName+"/asimovToys"):
                print("INFO: Making directory", dirName+"/asimovToys")
                Path(dirName+"/asimovToys").mkdir(exist_ok=True)
            betasToSubmit = GetBetasToSubmit(mass)
            betaLists = SplitBetasBetweenCondorSubmissions(betasToSubmit,nBetasPerCondorSub)
            print("INFO: for mass {}, use beta lists".format(mass), betaLists)
            nCondorFiles = len(betaLists)
            with Progress() as progress:
                betaIdsToGen = []
                for betaId in betasToSubmit:
                    asimovToyFile = FindAsimovToyDataBetaId(mass, betaId, dirName+"/asimovToys")
                    if asimovToyFile is None:
                        betaIdsToGen.append(betaId)
                    else:
                        dictAsimovToysByBetaId[betaId] = asimovToyFile
                if len(betaIdsToGen):
                    task_id = progress.add_task("[cyan]Generating asimov toys for mass {}".format(mass),total = len(betaIdsToGen))
                    with multiprocessing.Pool(ncores) as pool:
                        partial_quit = partial(ErrorCallback,pool)
                        for betaId in betaIdsToGen:
                            betaVal = betasToScan[betaId]
                            signalScaleFactor = betaVal**2
                            try:
                                pool.apply_async(GenerateAsimovToyDataBetaId,[cardWorkspace,mass,dirName+"/asimovToys",dictAsimovToysByBetaId, betaId, signalScaleFactor],callback = lambda x: progress.advance(task_id), error_callback = partial_quit)
                            except Exception as e:
                                print("ERROR: caught exception in asimov toy job for mass {}, betaId {}".format(mass, betaId))
                                traceback.print_exc()
                                pool.terminate()
                                exit(-2)
                        pool.close()
                        pool.join()
            rValuesAtBeta = rValuesByMassBetaAndQuantile[str(mass)]
            if options.generateGrids:
                if not os.path.isdir(dirName+"/gridGen"):
                    print("INFO: Making directory", dirName+"/gridGen", flush=True)
                    Path(dirName+"/gridGen").mkdir(exist_ok=True)
                task_id = progress.add_task("[cyan]Writing gridGen files for mass {}".format(mass),total = len(betaLists)*len(quantilesExpected))
                with multiprocessing.Pool(ncores) as pool:
                    partial_quit = partial(ErrorCallback,pool)
                    for quant in quantilesExpected:
                        quantStr = str(quant).replace(".","p")
                        if not os.path.isdir(dirName+"/gridGen/hybridNewGridGen.M{}.{}".format(mass,quantStr)):
                            print("INFO: Making directory", dirName+"/gridGen/hybridNewGridGen.M{}.{}".format(mass,quantStr))
                            Path(dirName+"/gridGen/hybridNewGridGen.M{}.{}".format(mass,quantStr)).mkdir(exist_ok=True)
                        for betaList in betaLists:
                            condorDir = dirName+"/gridGen/hybridNewGridGen.M{}.{}/betaIds{}To{}".format(mass,quantStr,betaList[0],betaList[-1])
                            if not os.path.isdir(condorDir):
                                Path(condorDir).mkdir(exist_ok=True)
                            asimovToys = []
                            for betaId in betaList:
                                asimovToys.append(str(dictAsimovToysByBetaId[betaId]))
                            args = [str(cardWorkspace),mass,quant,betaList,condorDir,asimovToys]
                            try:
                                pool.apply_async(WriteCondorFilesGridGen,args,callback = lambda x: progress.advance(task_id), error_callback = partial_quit)
                            except Exception as e:
                                print("ERROR: caught exception writing gridGen files for mass {}, quantile {}, betaIds {} to {}".format(mass, quant, betaList[0], betaList[-1]))
                                traceback.print_exc()
                                pool.terminate()
                                exit(-2)
                    pool.close()
                    pool.join()
            elif options.submitLimitJobsAfterGridGen:
                if not os.path.isdir(dirName+"/limits"):
                    print("INFO: Making directory", dirName+"/limits", flush=True)
                    Path(dirName+"/limits").mkdir(exist_ok=True)
                with multiprocessing.Pool(ncores) as pool:
                    partial_quit = partial(ErrorCallback,pool)
                    for quant in quantilesExpected:
                        quantStr = str(quant).replace(".","p")
                        if not os.path.isdir(dirName+"/limits/hybridNewLimits.M{}.{}".format(mass,quantStr)):
                            print("INFO: Making directory", dirName+"/limits/hybridNewLimits.M{}.{}".format(mass,quantStr))
                            Path(dirName+"/limits/hybridNewLimits.M{}.{}".format(mass,quantStr)).mkdir(exist_ok=True)
                        for betaList in betaLists:
                            try:
                                condorDir = dirName+"/limits/hybridNewLimits.M{}.{}/betaIds{}To{}".format(mass,quantStr,betaList[0],betaList[-1])
                                if not os.path.isdir(condorDir):
                                    Path(condorDir).mkdir(exist_ok=True)
                                asimovToys = []
                                for betaId in betaList:
                                    asimovToys.append(str(dictAsimovToysByBetaId[betaId]))
                                args = [str(cardWorkspace),mass,quant,betaList,condorDir,asimovToys,gridGenFilenamesByMassBetaAndQuantile[str(mass)][str(quant)]]
                                pool.apply_async(WriteCondorFilesLimits,args,error_callback = partial_quit)
                            except Exception as e:
                                print("ERROR: caught exception writing limit job files for mass {}, quantile {}, betaIds {} to {}".format(mass, quant, betaList[0], betaList[-1]))
                                traceback.print_exc()
                                pool.terminate()
                                exit(-2)
                    pool.close()
                    pool.join()
