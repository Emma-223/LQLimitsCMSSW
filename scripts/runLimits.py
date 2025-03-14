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

if not 'LQANA' in os.environ:
    raise RuntimeError('Need to define $LQANA environment variable to point to your git clone of the rootNtupleAnalyzerV2 repo.')
sys.path.append(os.getenv("LQANA").rstrip("/") + "/scripts/")

from combineCommon import SeparateDatacards, GetYearAndIntLumiFromDatacard

if not 'CMSSW_VERSION' in os.environ:
    raise RuntimeError('You must be in a CMSSW environment (release area + cmsenv)')

ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')


def ErrorCallback(pool, e):
    print("Parallel job failed")
    traceback.print_exception(type(e), e, e.__traceback__)
    # raise e
    # raise RuntimeError("Parallel job failed with exception; exiting.")
    pool.terminate()


def GetExecTimeStr(startTime, stopTime):
    delta = datetime.timedelta(seconds=stopTime - startTime)
    hours, rem = divmod(delta.seconds, 3600)
    minutes, seconds = divmod(rem, 60)
    execTimeStr = "{}s".format(seconds)
    if minutes > 0:
        execTimeStr = "{}m".format(minutes) + execTimeStr
    if hours > 0:
        execTimeStr = "{}h".format(hours) + execTimeStr
    return execTimeStr


def RunCommand(cmd, workingDir=None, env=None, suppressOutput=False):
    origCmd = cmd
    if not suppressOutput:
        print(colored("\t{}".format(origCmd), "green"), flush=True)
    try:
        useShell = False
        # if env is None:
        #     useShell = True
        #     cmd = shlex.split(cmd)
        cmd = shlex.split(cmd)
        process = subprocess.run(cmd, shell=useShell, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=workingDir, env=env)
    except subprocess.CalledProcessError as e:
        result = ansi_escape.sub('', e.output.decode())
        if suppressOutput:
            print(colored("\t{}".format(origCmd), "green"), flush=True)
        print(colored("RunCommand had an error; output: {}".format(result), "red"), flush=True)
        raise e
    if not suppressOutput:
        print(process.stdout.decode(), flush=True)
        #print(process.stderr.decode())


def GetEOSPathWithFullURL(eosPath):
    if eosPath.startswith("/eos/user"):
        return "root://eosuser.cern.ch/" + eosPath
    elif eosPath.startswith("/eos/cms"):
        return "root://eoscms.cern.ch/" + eosPath
    raise RuntimeError("Don't know how to get full EOS path with URL from eosPath provided = '{}'. Please implement.".format(eosPath))


def ConvertDatacardToWorkspace(datacard, mass):
    workspaceFileName = datacard+".m{}.root".format(mass)
    cmd = 'text2workspace.py --channel-masks -m {} {} -o {} '.format(mass, datacard, workspaceFileName)
    RunCommand(cmd)
    return Path(workspaceFileName).resolve()


def FindCardWorkspace(cardFile, mass, exceptOnAbsence=False):
    workspaceFileName = cardFile+".m{}.root".format(mass)
    listOfWorkspaceFiles = sorted(glob.glob(workspaceFileName), key=os.path.getmtime)
    if exceptOnAbsence and len(listOfWorkspaceFiles) != 1:
        raise RuntimeError("Globbing for {} did not result in one file as expected, but {}".format(workspaceFileName, listOfWorkspaceFiles))
    if len(listOfWorkspaceFiles) > 0:
        return Path(listOfWorkspaceFiles[-1]).resolve()


def GenerateAsimovToyData(workspace, mass, toysDir, dictAsimovToysByScaleFactor, signalScaleFactor, expectSignal=False):
    #workspace, mass, toysDir, dictAsimovToysByScaleFactor, signalScaleFactor = args
    cmd = 'combine -M GenerateOnly {} -t -1 --seed -1 --saveToys -n .asimov.{} -m {} {}'.format(workspace, signalScaleFactor, mass, commonCombineArgs.format(signalScaleFactor))
    if expectSignal==True:
        cmd += " --expectSignal 1"
    RunCommand(cmd, toysDir, None, True)
    toyFile = Path(sorted(glob.glob(toysDir+'/higgsCombine.asimov.{}.GenerateOnly.mH{}.*.root'.format(signalScaleFactor, mass)), key=os.path.getmtime)[-1]).resolve()
    dictAsimovToysByScaleFactor[signalScaleFactor] = toyFile
    return toyFile


def FindAsimovToyData(mass, signalScaleFactor, toysDir):
    listOfToyFiles = sorted(glob.glob(toysDir+'/higgsCombine.asimov.{}.GenerateOnly.mH{}.*.root'.format(signalScaleFactor, mass)), key=os.path.getmtime)
    if len(listOfToyFiles) > 0:
        return Path(listOfToyFiles[-1]).resolve()


def FindFile(globString):
    # fileList = sorted(glob.glob(globString), key=os.path.getmtime)
    # if len(fileList) != 1:
    #     raise RuntimeError("Globbing for {} did not result in one file as expected, but {}".format(globString, fileList))
    # return Path(fileList[-1]).resolve()
    if not isinstance(globString, list):
        globString = [globString]
    return FindFileMultiple(globString)


def FindFileMultiple(globStringList):
    fileList = []
    globString = ""
    for globStr in globStringList:
        globString = globStr
        fileList = sorted(glob.glob(globString), key=os.path.getmtime)
        if len(fileList):
            break
    if len(fileList) <1: #!= 1:
        raise RuntimeError("Globbing for {} did not result in one file as expected; using globString={}, result is {}".format(globStringList, globString, fileList))
    return Path(fileList[-1]).resolve()


def GetFileList(globString):
    # fileList = sorted(glob.glob(globString), key=os.path.getmtime)
    # if len(fileList) < 1:
    #     raise RuntimeError("Globbing for {} did not result in any files.".format(globString))
    # return fileList
    if not isinstance(globString, list):
        globString = [globString]
    return GetFileListMultiple(globString)


def GetFileListMultiple(globStringList):
    fileList = []
    globString = ""
    for globStr in globStringList:
        globString = globStr
        fileList = sorted(glob.glob(globString), key=os.path.getmtime)
        if len(fileList):
            break
    if len(fileList) < 1:
        raise RuntimeError("Globbing for {} did not result in any files; using globString={}, result is {}".format(globStringList, globString, fileList))
    return fileList

def GetFilesInDirXrd(dirName):
    eosPathNoUrl = "/"+dirName.split("//")[-1]
    cmd = ["xrdfs","root://eoscms.cern.ch","ls","{}".format(eosPathNoUrl)] #returns an error if dir does not exist. Does not support wildcards.
    lsOutput = subprocess.check_output(cmd).decode('ASCII')
    fileList = lsOutput.split('\n')
    return fileList

def MakeCondorDirs(dirName):
    condorDir = dirName.strip("/")+"/condor"
    # condorSubfileDir = dirName.strip("/")+"/condor/subFiles"
    # condorOutDir = dirName.strip("/")+"/condor/out"
    # condorErrDir = dirName.strip("/")+"/condor/error"
    # condorLogDir = dirName.strip("/")+"/condor/log"
    dirsToHave = [condorDir]  # , condorSubfileDir, condorOutDir, condorErrDir, condorLogDir]
    for idir in dirsToHave:
        if not os.path.isdir(idir):
            print("INFO: Making directory", idir, flush=True)
            Path(idir).mkdir(exist_ok=True)

def GetFilenamesFromCondorOutFiles(dirName, nFilesExpected=1, clusterId = ''):
    if clusterId=='':
        condorOutFileList = GetFileList(dirName+"/condor.*.out")
    else:
        condorOutFileList = GetFileList(dirName+"/condor.{}*.out".format(clusterId))
    eosFileList = []
    filesMissingOutput = []
    for file in condorOutFileList:
        foundOutput = False
        with open(file,'r') as f:
            lines = f.readlines()
        for i in range(1,len(lines)):
            l = lines[len(lines)-i] #We expect it to be the last line, so start from the end
            if "output eos file" in l.lower():
                eosFileList.append("root://eoscms.cern.ch/"+l.split()[-1])
                foundOutput = True
                break
        if not foundOutput:
            filesMissingOutput.append(file)
    if len(filesMissingOutput):
        print("The following condor*.out files are missing the eos output file path: ",filesMissingOutput)
        raise RuntimeError("Can't find eos files in all condor*.out files")
    if len(eosFileList) < nFilesExpected:
        raise RuntimeError("In directory {}, expected {} files, found {}".format(dirName, nFilesExpected, len(eosFileList)))
    if nFilesExpected > 1:
        return eosFileList
    else:
        return eosFileList[0]

def CheckAllGridGenJobs(massList, quantileList, nFilesExpected):
    nMissingFiles = 0
    dirname = options.name+"/gridGen/hybridNewGridGen.M{}.{}"
    outputFiles = {}
    for mass in massList:
        for quant in quantileList:
            print("INFO: Check grid gen jobs for mass {} quantile {}".format(mass,quant))
            q = str(quant).replace(".","p")
            subdirList = os.listdir(dirname.format(mass,q))
            for d in subdirList:
                condorIds = []
                condorFiles = os.listdir(dirname.format(mass,q)+"/"+d)
                for f in condorFiles:
                    if ".sh" in f or ".sub" in f:
                        continue
                    clusterId = int(f.split(".")[1])
                    if not clusterId in condorIds:
                        condorIds.append(clusterId)
                if len(condorIds)>1:
                    idToUse = max(condorIds)
                    #print("WARNING: found multiple condor cluster Ids in directory {}. Using highest Id {}".format(dirname.format(mass,q)+"/"+d, idToUse))
                else:
                    idToUse = condorIds[0]
                try:
                    outputFiles[dirname.format(mass,q)+"/"+d] = GetFilenamesFromCondorOutFiles(dirname.format(mass,q)+"/"+d,nFilesExpected, idToUse)
                except Exception as e:
                    print(e)
                    nMissingFiles+=1
    if nMissingFiles > 0:
        raise RuntimeError("Some grid gen output is missing")
    return outputFiles

def WriteCondorSubFile(condorDir, filename, shFilename, nJobs, queue="workday", doGridGen=True, extractLimits=False, inputFiles=[]):
    # scanJobs = 1
    # if doGridGen:
    #     scanJobs = gridScanPoints / gridPointsPerJob
    maxMaterialize = 1
    if options.cmsConnectMode:
        basename = str(shFilename).replace("condor_", "").replace(".sh", "")
    else:
        basename = shFilename.name.replace("condor_", "").replace(".sh", "")
    parentDir = str(shFilename.parent.resolve().name)
    if options.doSignificance:
        parentDir+=".significance"
    baseDir = str(Path(condorDir).parent.parent.resolve().name)
    with open(filename, "w") as subfile:
        if options.cmsConnectMode:
            subfile.write("Universe = vanilla\n")
            subfile.write("executable = " + str(shFilename) + "\n")
        else:
            subfile.write("executable = " + shFilename.name + "\n")
        subfile.write("arguments = $(ProcId)\n")
        subfile.write("output                = " + basename + ".$(ClusterId).$(ProcId).out\n")
        subfile.write("error                 = " + basename + ".$(ClusterId).$(ProcId).err\n")
        subfile.write("log                   = " + basename + ".$(ClusterId).$(ProcId).log\n")
        if options.cmsConnectMode:
            subfile.write('+REQUIRED_OS = "rhel9"\n')
            subfile.write('periodic_hold = ((time() - EnteredCurrentStatus) > 28800)\n') #28800 seconds = 8 hours
        subfile.write("\n")
        if options.cmsConnectMode:
            fullPathInputFiles = [
                options.cmsswSandbox,
                "/home/emma-pearson/cmssw_setup_copyForDebug.sh",
            ]
        else:
            fullPathInputFiles = []
        if len(inputFiles):
            for thisFile in inputFiles:
                if thisFile.startswith("/eos"):
                    thisFile = GetEOSPathWithFullURL(thisFile)
                fullPathInputFiles.append(thisFile)
        if len(fullPathInputFiles):
            subfile.write("transfer_input_files = {}\n".format(",".join(fullPathInputFiles)))
            subfile.write("\n")
        # subfile.write('transfer_output_remaps = "combine_logger.out=out/' + basename + '.$(ClusterId).$(ProcId).combine_logger.out"\n')
        subfile.write("\n")
        subfile.write("should_transfer_files = YES\n")
        if not options.cmsConnectMode:
            subfile.write("when_to_transfer_output = ON_EXIT\n")
            subfile.write("output_destination = "+eosDir+"/"+baseDir+"/"+parentDir+"/\n")
            subfile.write("MY.XRDCP_CREATE_DIR = True\n")
            subfile.write("\n")
        subfile.write("# Send the job to Held state on failure.\n")
        subfile.write("on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n")
        subfile.write("\n")
        subfile.write("# Periodically retry the jobs every 10 minutes, up to a maximum of 5 retries.\n")
        subfile.write("periodic_release =  (NumJobStarts < 10) && ((CurrentTime - EnteredCurrentStatus) > 300)\n")
        subfile.write("\n")
        # if doGridGen:
        #     subfile.write("max_materialize={}\n".format(maxMaterialize))
        if not options.cmsConnectMode:
            subfile.write("+JobFlavour=\"{}\"\n".format(queue))
            subfile.write("MY.WantOS = \"el9\"\n")
        # subfile.write("queue {}\n".format(int(scanJobs)))
        subfile.write("queue {}\n".format(nJobs))


def WriteCondorShFile(condorDir, filename, mass, combineCmds=[], limitCmds=[], signalScaleFactor=1.0, quantiles=[], betaIndex=-1):
    if not isinstance(quantiles, list):
        quantiles = [quantiles]
    jobIdx = 0
    with open(filename, "w") as shfile:
        shfile.write("#!/bin/sh\n")
        shfile.write("ulimit -s unlimited\n")
        shfile.write("set -e\n")
        if options.cmsConnectMode:
            shfile.write("condorNodeBaseDir=$PWD\n")
            shfile.write("source ./cmssw_setup_copyForDebug.sh\n")
            shfile.write("export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch\n")
            shfile.write("source $VO_CMS_SW_DIR/cmsset_default.sh\n")
            sandbox = Path(options.cmsswSandbox).resolve().name
            shfile.write("cmssw_setup {}\n".format(sandbox))
            shfile.write("cd $CMSSW_BASE\n")
            for i,c in enumerate(combineCmds):
                c = c.replace("./combine","combine")
                combineCmds[i] = c
            for i,c in enumerate(limitCmds):
                c = c.replace("./combine","combine")
                limitCmds[i] = c

        else:
            shfile.write("cd {}\n".format(condorDir))
            shfile.write("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
            shfile.write("eval `scramv1 runtime -sh`\n")
            shfile.write("cd -\n")
        shfile.write("\n")
        rootFileString = ""
        sigSFRound = round(signalScaleFactor, 6)
        if betaIndex != -1:
            rootFileString = ".betaIndex"+str(betaIndex)
        combinedOutputFiles = []
        for idx, quantile in enumerate(quantiles):
            if len(combineCmds):
                rMin = -1
                rMax = -1
                if quantile != -1 and not options.doSignificance:
                   rMin, rMax = GetRMinAndRMax(mass, quantile, signalScaleFactor)
                if quantile == -2: #Observed limit
                    rMin, rMax = GetRMinAndRMax(mass, 0.5, signalScaleFactor)
                if rMin > 0 and rMax > 0:
                    shfile.write("if [ $1 -eq {} ]; then\n".format(jobIdx))
                    if quantile < 0.2 and quantile > 0:
                        pointsPerJob = gridPointsPerJobLowQuantiles
                    else:
                        pointsPerJob = gridPointsPerJob
                    #print("Info: using {} grid points per job for quantile {}".format(pointsPerJob, quantile))
                    stepSize = (rMax-rMin)/gridScanPoints
                    if gridScanPoints < pointsPerJob:
                        pointsPerJob = gridScanPoints
                    scanPoints = np.linspace(rMin, rMax, num=gridScanPoints)
                    # print("INFO: Writing commands to compute grid of {} limits in the range r=[{}, {}] by steps of {} for quantile {}".format(gridScanPoints, rMin, rMax, stepSize, quantile))
                    for rValIdx, rVal in enumerate(scanPoints):
                        thisStepCmd = combineCmds[idx]
                        thisStepCmd += ' -n .signalScaleFactor{}.POINT.{}'.format(sigSFRound, rVal)
                        if len(rootFileString):
                            thisStepCmd += "." + rootFileString
                        thisStepCmd += ' --singlePoint {}'.format(rVal)
                        shfile.write("  {}\n".format(thisStepCmd))
                        if (rValIdx+1) % pointsPerJob == 0:
                            quant = "quant{}.".format(quantile) if quantile > 0 else ""
                            combinedOutputFile = "higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}grid{}.root".format(sigSFRound, mass, quant, jobIdx)
                            shfile.write("  hadd -fk207 {} higgsCombine.*POINT.*.root\n".format(combinedOutputFile))
                            # shfile.write("  ls -ltrph\n")
                            shfile.write("  rm higgsCombine.*POINT.*root\n")  # remove individual point files so that condor doesn't transfer them back
                            if rValIdx+1 < len(scanPoints):
                                jobIdx += 1
                                shfile.write("fi\n")
                                shfile.write("if [ $1 -eq {} ]; then\n".format(jobIdx))
                            combinedOutputFiles.append(combinedOutputFile)
                    if len(limitCmds):
                        shfile.write(limitCmds[idx].format(combinedOutputFile)+"\n")
                    shfile.write("fi\n")
                    jobIdx += 1
                else:
                    if options.doSignificance:
                        njobs = 10
                        for ijob in range(njobs): #10 jobs of 10000 toys each
                            shfile.write("if [ $1 -eq {} ]; then\n".format(ijob))
                            cmd = combineCmds[idx]
                            cmd += ' -n .jobId{}'.format(ijob)
                            shfile.write(cmd+"\n")
                            shfile.write("fi\n")
                            jobIdx+=1
                    else:
                        if len(rootFileString):
                            combineCmds[idx] += ' -n {}'.format(rootFileString)
                        # print("WriteCondorShFile for mass {}: write cmd='{}'".format(mass, combineCmd))
                        shfile.write("  {}\n".format(combineCmds[idx]))
            if len(limitCmds) and (not len(combineCmds) or not len(combineCmds[idx])):
                quant = "quant{}.".format(quantile) if quantile > 0 else ""
                if options.cmsConnectMode:
                    nFiles = int(math.ceil(gridScanPoints/gridPointsPerJob))
                    #filenameBase = eosDir+"/"+options.name+"/hybridNewGridGen.M{}.{}/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}grid{}.root"
                    #filesToCheck = []
                    #for i in range(nFiles):
                    #    f = filenameBase.format(mass,str(quantile).replace(".","p"),sigSFRound,mass,quant,i)
                    #    shfile.write("xrdcp -fs {} {}\n".format(f, f.split("/")[-1]))
                    #    filesToCheck.append(f)
                    gridGenCondorDir = options.name+"/gridGen/hybridNewGridGen.M{}.{}/signalScaleFactor{}".format(mass,str(quantile).replace(".","p"),sigSFRound)
                    eosFileList = gridGenFilesByCondorDir[gridGenCondorDir]#GetFilenamesFromCondorOutFiles(gridGenCondorDir,nFiles)
                    for f in eosFileList:
                        shfile.write("xrdcp -fs {} {}\n".format(f, f.split("/")[-1]))

                combinedOutputFile = "higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}gridAll.root".format(sigSFRound, mass, quant)
                shfile.write("hadd -fk207 {} higgsCombine.*.root\n".format(combinedOutputFile))
                shfile.write(limitCmds[idx].format(combinedOutputFile)+"\n")
                combinedOutputFiles.append(combinedOutputFile)
        if options.cmsConnectMode:
            shfile.write("\n")
            parentDir = str(Path(filename).parent.resolve().name)
            baseDir = str(Path(filename).parent.parent.resolve().name)
            eosPath = eosDir+"/"+options.name+"/"+baseDir#+"/"+parentDir
            if len(limitCmds):
                if "0.5" in quant:
                    quantWithTrailing0s = quant.replace('0.5','0.500')
                elif "0.16" in quant:
                    quantWithTrailing0s = quant.replace('0.16','0.160')
                elif "0.84" in quant:
                    quantWithTrailing0s = quant.replace('0.84','0.840')
                else:
                    quantWithTrailing0s = quant
                files = [
                    "higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}gridAll.root".format(sigSFRound, mass, quant),
                    "higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}root".format(sigSFRound, mass, quantWithTrailing0s),
                    "limit_scan_signalScaleFactor{}_m{}_{}pdf".format(sigSFRound, mass, quant),
                ]
                for f in files:
                    shfile.write("xrdcp -fs {} {}/{}\n".format(f,eosPath,f))
                shfile.write("outputEOSFile=$(xrdfs root://eoscms.cern.ch ls /{}/{})\n".format(eosPath.split("//")[-1],files[1]))
                shfile.write('echo "output eos file: $outputEOSFile "')
            else:
                shfile.write("ls higgsCombine*.root\n")
                shfile.write("outputRootFile=$(ls higgsCombine*.root)\n")
                shfile.write("xrdcp -fs $outputRootFile {}/$outputRootFile\n".format(eosPath))
                shfile.write("outputEOSFile=$(xrdfs root://eoscms.cern.ch ls /{}/$outputRootFile)\n".format(eosPath.split("//")[-1]))
                shfile.write('echo "output eos file: $outputEOSFile "')
    return combinedOutputFiles, jobIdx if jobIdx > 0 else 1


def GetHybridNewCommandArgs(workspace, mass, dirName, quantiles, genAsimovToyFile, signalScaleFactor, batch):
    if not isinstance(quantiles, list):
        quantiles = [quantiles]
    rAbsAcc = 0.0001
    clsAcc = 0 if batch else 0.001
    toysDefault = 10000
    # toys = 500  # reduced for shape-based limits as test, but still took hours
    cmds = []
    for quantile in quantiles:
        if not options.doSignificance:
            if quantile>0:
                rMin, rMax = GetRMinAndRMax(mass, quantile, signalScaleFactor)
            else:
                rMin, rMax = GetRMinAndRMax(mass, 0.5, signalScaleFactor)
            rMin*=0.9
            rMax*=1.1
        else:
            rMin = -1
            rMax = -1
        toys = toysDefault
        if quantile == 0.025:
            toys = 40000
        elif quantile == 0.16:
            toys = 30000
        if options.cmsConnectMode:
            cmd = '-d $condorNodeBaseDir/{}'.format(str(workspace).split("/")[-1])
        else:
            cmd = '-d {}'.format(str(workspace).split("/")[-1] if doBatch else workspace)
        cmd += ' -v1'
        cmd += ' -M HybridNew'
        if options.doSignificance:
            cmd += ' --LHCmode LHC-significance'
        else:
            cmd += ' --LHCmode LHC-limits'
        cmd += ' --saveHybridResult'
        cmd += ' --saveToys'
        cmd += ' --seed -1'
        if not batch:
            cmd += ' --iterations 2'
            cmd += ' --rAbsAcc {}'.format(rAbsAcc)
        if not options.doSignificance:
            cmd += ' --clsAcc {}'.format(clsAcc)
        if rMin >= 0 and rMax >= 0:
            cmd += ' --rMin {}'.format(rMin)
            cmd += ' --rMax {}'.format(rMax)
        cmd += ' -T {}'.format(toys)
        cmd += ' -m {}'.format(mass)
        if not batch:
            cmd += ' -n .signalScaleFactor{}'.format(round(signalScaleFactor, 6))
        # cmd += ' -H AsymptoticLimits'
        # cmd += ' --fork 4'
        cmd += commonCombineArgs.format(signalScaleFactor)
        if quantile > 0:
            cmd += ' --expectedFromGrid {}'.format(quantile)
            #print("INFO: Run with --expectedFromGrid {}".format(quantile))
        if genAsimovToyFile != "":
            if options.cmsConnectMode:
                cmd += ' -D $condorNodeBaseDir/{}:toys/toy_asimov'.format(str(genAsimovToyFile).split("/")[-1])
            else:
                cmd += ' -D {}:toys/toy_asimov'.format(str(genAsimovToyFile).split("/")[-1] if doBatch else genAsimovToyFile)
            if quantile < 0:
                print("INFO: Found quantile {}: Running observed limit with Asimov toys".format(quantile))
        cmds.append(cmd)
    return cmds


def GetRMinAndRMax(mass, quantileExp, signalScaleFactor=1.0, betaId = -1):
    #print("starting rvalue = {} for quantile {} and mass{}".format(rValuesByMassAndQuantile[str(mass)][str(quantileExp)], quantileExp, mass))
    if betaId < 0:# or rescaleSignal:
        rValuesByQuantile = rValuesByMassAndQuantile[str(mass)]
        if quantileExp == 0.025 and mass > 800 and mass <=2000: # adjust scan range upwards for lowest quantile and higher masses
            rMax = rValuesByQuantile[str(quantileExp)]*2.5
            rMin = rValuesByQuantile[str(quantileExp)]*1.25
        #elif quantileExp == 0.025 and mass > 1000 and mass < 2000:
         #   rMax = rValuesByQuantile[str(quantileExp)]*3.5
          #  rMin = rValuesByQuantile[str(quantileExp)]*2.25
        elif quantileExp == 0.16 and mass > 1000 and mass <2000:
            rMax = rValuesByQuantile[str(quantileExp)]*2.0
            rMin = rValuesByQuantile[str(quantileExp)]*1.0
        elif quantileExp == 0.975 and mass > 1500:  # adjust scan range downwards here
             rMax = rValuesByQuantile[str(quantileExp)]*1.0
             rMin = rValuesByQuantile[str(quantileExp)]*0.45
        elif (quantileExp == 0.025 or quantileExp==0.16) and mass >=2000 and mass <2500:
            rMax = rValuesByQuantile[str(quantileExp)]*5.0
            rMin = rValuesByQuantile[str(quantileExp)]*1.0
        elif quantileExp==0.84 and mass >= 2000:
            rMax = rValuesByQuantile[str(quantileExp)]*1.0
            rMin = rValuesByQuantile[str(quantileExp)]*0.45
        elif (quantileExp==0.025 or quantileExp==0.16) and mass >= 2500:
            rMax = rValuesByQuantile[str(quantileExp)]*5.0
            rMin = rValuesByQuantile[str(quantileExp)]*1.0
        elif quantileExp==0.5 and mass >= 2800:
            rMax = rValuesByQuantile[str(quantileExp)]*2
            rMin = rValuesByQuantile[str(quantileExp)]*1
        else:
            rMax = rValuesByQuantile[str(quantileExp)]*1.3
            rMin = rValuesByQuantile[str(quantileExp)]*0.75
    else:
        #beta = math.sqrt(signalScaleFactor)
        if quantileExp == 0.025 and mass > 800 and mass <1200:  # adjust scan range upwards for lowest quantile and higher masses
            rMax = rValuesAtBeta[str(betaId)][str(quantileExp)]*1.8
            rMin = rValuesAtBeta[str(betaId)][str(quantileExp)]*0.85
        elif quantileExp == 0.025 and mass >= 1200:
            rMax = rValuesAtBeta[str(betaId)][str(quantileExp)]*2.5
            rMin = rValuesAtBeta[str(betaId)][str(quantileExp)]*1.25
        elif quantileExp == 0.16 and mass > 1000:
            rMax = rValuesAtBeta[str(betaId)][str(quantileExp)]*2.0
            rMin = rValuesAtBeta[str(betaId)][str(quantileExp)]*1.0
        elif quantileExp == 0.975 and mass >= 1200:  # adjust scan range downwards here
            rMax = rValuesAtBeta[str(beta)][str(quantileExp)]*0.9
            rMin = rValuesAtBeta[str(beta)][str(quantileExp)]*0.3
        else:
            rMax = rValuesAtBeta[str(betaId)][str(quantileExp)]*1.3
            rMin = rValuesAtBeta[str(betaId)][str(quantileExp)]*0.75
    #print("rValue range for quantile {}, mass {} = {} - {}".format(quantileExp, mass, rMin, rMax))
    return rMin, rMax


def GetAsymptoticCommandArgs(workspace, mass, dirName, blinded, signalScaleFactor, rMin, rMax):
    rAbsAcc = 0.0001
    # rMax = 150
    if options.cmsConnectMode:
        cmd = ' $condorNodeBaseDir/{}'.format(str(workspace).split("/")[-1])
    else:
        cmd = ' {}'.format(str(workspace).split("/")[-1] if doBatch else workspace)
    cmd += ' -v1'
    cmd += ' -M AsymptoticLimits'
    cmd += ' --seed -1'
    cmd += ' --strict'
    cmd += ' --rMax {}'.format(rMax)
    if rMin >= 0:
        cmd += ' --rMin {}'.format(rMin)
    cmd += ' -m {}'.format(mass)
    cmd += ' --rAbsAcc {}'.format(rAbsAcc)
    # cmd += ' -n .signalScaleFactor{}'.format(round(signalScaleFactor, 6))
    cmd += commonCombineArgs.format(signalScaleFactor)
    if blinded:
        cmd += ' --run blind --expectSignal 0'
    return cmd


def SubmitLimitJobsCrab(mass, dirName, listFailedCommands, quantile, signalScaleFactor, rMin, rMax, taskName, cmdArgs, limitCmd, queue="workday", inputFiles=[]):
    if len(limitCmd):
        raise RuntimeError("Currently, we don't support extracting the limits with crab.")
    subfileDir = ""
    if "gridgen" in taskName.lower():
        subfileDir = dirName.strip("/") + "/gridGen/" + taskName
    elif "limits" in taskName.lower():
        subfileDir = dirName.strip("/") + "/limits/" + taskName
    Path(subfileDir).mkdir(parents=True)
    betaIndex = -1
    if "betaIndex" in taskName:
        betaIndex = taskName[taskName.find("betaIndex")+9:]
    customCrabFile = str(Path("scripts/custom_crab.py").resolve())
    cmd = "combineTool.py {} --singlePoint {}:{}:{} --job-mode crab3 --task-name {} --custom-crab {}".format(cmdArgs, rMin, rMax, gridScanPoints, taskName, customCrabFile)
    if options.dryRun:
        cmd += " --dry-run"
    if len(inputFiles):
        cmd += " --crab-extra-files {}".format(" ".join(inputFiles))
    # runCommandArgs = [cmd, subfileDir, None, True]
    runCommandArgs = [cmd, subfileDir]
    trial = 0
    maxTries = 5
    while trial < maxTries:
        try:
            RunCommand(*runCommandArgs)
        except subprocess.CalledProcessError as e:
            if trial < maxTries-1:
                print(colored("Caught exception running combineTool.py; retry", "yellow"), flush=True)
            trial += 1
            continue
        break
    if trial >= maxTries:
        listFailedCommands.append(" ".join(runCommandArgs[0:1]))
        return " ".join(runCommandArgs[0:1])


def SubmitLimitJobsBatch(mass, dirName, listFailedCommands, quantiles, signalScaleFactor, taskName, cmdArgs, limitCmd, queue="workday", inputFiles=[]):
    # print("SubmitLimitJobsBatch() called with cmdArgs='{}'".format(cmdArgs))
    # condorDir = dirName.strip("/")+"/condor"
    condorDir = dirName.strip("/")
    # condorSubfileDir = dirName.strip("/")+"/condor/subFiles"
    doGridGen = False
    # condorSubfileDir = ""
    condorSubfileDir = condorDir + "/" + taskName
    if "gridgen" in taskName.lower():
        if options.doSignificance:
            condorSubfileDir = condorDir + "/significance/" + taskName
        else:
            condorSubfileDir = condorDir + "/gridGen/" + taskName
        doGridGen = True
    elif "limits" in taskName.lower():
        condorSubfileDir = condorDir + "/limits/" + taskName
    elif "asymptotic" in taskName.lower():
        condorSubfileDir = condorDir + "/asymptoticLimits/" + taskName
    elif "significance" in taskName.lower():
        condorSubfileDir = condorDir + "/significance/" + taskName

    if options.doBetaScan:
        condorSubfileDir = condorSubfileDir.replace(".signalScaleFactor","/signalScaleFactor")
    Path(condorSubfileDir).mkdir(parents=True, exist_ok=True)
    cmds = []
    if len(cmdArgs):
        cmds = ["./combine " + args for args in cmdArgs] if isinstance(cmdArgs, list) else ["./combine " + cmdArgs]
    shFilename = condorSubfileDir + "/condor.sh"
    condorFile = shFilename.replace(".sh", ".sub")
    betaIndex = -1
    if "betaIndex" in taskName:
        betaIndex = taskName[taskName.find("betaIndex")+9:]
    _, nJobs = WriteCondorShFile(Path(condorDir).resolve(), shFilename, mass, cmds, limitCmd, signalScaleFactor, quantiles, betaIndex)
    WriteCondorSubFile(Path(condorSubfileDir).resolve(), condorFile, Path(shFilename), nJobs, queue, doGridGen, len(limitCmd) > 0, inputFiles)
    if not options.dryRun:
        cmd = "condor_submit {}".format(Path(condorFile).name)
        runCommandArgs = [cmd, condorSubfileDir, None, True]
        trial = 0
        maxTries = 5
        while trial < maxTries:
            try:
                RunCommand(*runCommandArgs)
            except subprocess.CalledProcessError as e:
                if trial < maxTries-1:
                    print(colored("Caught exception running condor_submit; retry", "yellow"), flush=True)
                trial += 1
                continue
            break
        if trial >= maxTries:
            listFailedCommands.append(" ".join(runCommandArgs[0:1]))
            return " ".join(runCommandArgs[0:1])
    #else:
     #   print("INFO: dry run enabled; not submitting {} to condor".format(str(Path(condorFile))))


def SubmitHybridNewBatch(args):
    workspace, mass, dirName, listFailedCommands, quantiles, genAsimovToyFile, signalScaleFactor = args
    quantileStr = "."+str(quantiles).replace(".", "p") if not isinstance(quantiles, list) else ""
    if options.submitLimitJobsAfterGridGen:
        queue = "longlunch"
        #FIXME: handle multiple quantiles here if needed
        quant = ""
        if not isinstance(quantiles, list) and quantiles > 0:
            quant = "quant{}.".format(quantiles)
        sigSFRound = round(signalScaleFactor, 6)
        if options.cmsConnectMode:
            inputFiles = []
            baseDirName = eosDirName.rstrip("/")+'/hybridNewGridGen.M{}{}/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}grid{}.root'
            nFiles = int(math.ceil(gridScanPoints/gridPointsPerJob))
            for i in range(nFiles):
                f = baseDirName.format(mass,quantileStr,sigSFRound,mass,quant,i)
                missingFiles = []
                if not f.split("//")[-1] in inputFileList:
                    missingFiles.append(f)
            #if len(missingFiles):
            #    print("did not find files: ",missingFiles)
            #    raise RuntimeError("Some input files are missing")
        else:
            if options.doSignificance:
                condorDir = eosDirName.rstrip("/")+"/hybridNewGridGen.M{}{}.significance".format(mass, quantileStr)
                globStrs = [str(Path(condorDir).resolve())+"/higgsCombine.*.HybridNew.mH{}.*.root".format(mass)]
            else:
                condorDir = dirName.strip("/")+"/gridGen/hybridNewGridGen.M{}{}".format(mass, quantileStr)
                if signalScaleFactor != 1.0:
                    condorDir += '.signalScaleFactor{}'.format(sigSFRound)
                globStrs = [str(Path(condorDir).resolve())+"/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}grid*.root".format(sigSFRound, mass, quant)]
                if signalScaleFactor != 1.0:
                    globStrs.append(eosDirName.rstrip("/")+'/hybridNewGridGen.M{}{}.signalScaleFactor{}/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}grid*.root'.format(mass, quantileStr, sigSFRound, sigSFRound, mass, quant))
                else:
                    globStrs.append(eosDirName.rstrip("/")+'/hybridNewGridGen.M{}{}/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}grid*.root'.format(mass, quantileStr, sigSFRound, mass, quant))
            inputFiles = GetFileList(globStrs)
        # of course, now the condorDir part has to be removed
        # inputFiles = [f[f.rfind("/")+1:] for f in inputFiles]
        inputFiles.append(str(workspace))
        combinedOutputFile = "higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}gridAll.root".format(sigSFRound, mass, quant)
        limitCmd = GetComputeLimitsFromGridCommand(workspace, mass, combinedOutputFile, quantiles, sigSFRound, toyFile=genAsimovToyFile)
        if genAsimovToyFile != "":
            inputFiles.append(str(genAsimovToyFile))
        if options.doSignificance:
            taskName = 'hybridNewSignificance.M{}{}'.format(mass, quantileStr)
        else:
            taskName = 'hybridNewLimits.M{}{}'.format(mass, quantileStr)
        if signalScaleFactor != 1.0:
            taskName += '.signalScaleFactor{}'.format(round(signalScaleFactor, 6))
        if options.crabMode:
            SubmitLimitJobsCrab(mass, dirName, listFailedCommands, quantiles, signalScaleFactor, taskName, "", limitCmd, queue, inputFiles)
        else:
            combineBin = shutil.which("combine")
            if combineBin is None:
                raise RuntimeError("No combine binary found in $PATH. Did you compile combine?")
            inputFiles.append(combineBin)
            SubmitLimitJobsBatch(mass, dirName, listFailedCommands, quantiles, signalScaleFactor, taskName, "", limitCmd, queue, inputFiles)
    else:
        cmdArgs = GetHybridNewCommandArgs(workspace, mass, dirName, quantiles, genAsimovToyFile, signalScaleFactor, True)
        queue = "tomorrow"
        if options.generateGridsOnly or gridScanPoints / gridPointsPerJob > 1 or options.doSignificance:
            taskName = 'hybridNewGridGen.M{}{}'.format(mass, quantileStr)
            limitCmds = []
            if options.doSignificance:
                queue = "longlunch"
        else:
            taskName = 'hybridNewGridGenAndLimits.M{}{}'.format(mass, quantileStr)
            sigSFRound = round(signalScaleFactor, 6)
            limitCmds = GetComputeLimitsFromGridCommand(workspace, mass, "", quantiles, sigSFRound)
        if signalScaleFactor != 1.0:
            taskName += '.signalScaleFactor{}'.format(round(signalScaleFactor, 6))
        inputFiles = [str(workspace)]
        if len(str(genAsimovToyFile)):
            inputFiles.append(str(genAsimovToyFile))
        if options.crabMode:
            SubmitLimitJobsCrab(mass, dirName, listFailedCommands, quantiles, signalScaleFactor, taskName, cmdArgs, limitCmds, queue, inputFiles)
        else:
            combineBin = shutil.which("combine")
            if combineBin is None:
                raise RuntimeError("No combine binary found in $PATH. Did you compile combine?")
            inputFiles.append(combineBin)
            SubmitLimitJobsBatch(mass, dirName, listFailedCommands, quantiles, signalScaleFactor, taskName, cmdArgs, limitCmds, queue, inputFiles)


def SubmitAsymptoticBatch(args):
    workspace, mass, dirName, listFailedCommands, blinded, signalScaleFactor, index = args
    rMin = 0
    rMax = 300
    cmdArgs = GetAsymptoticCommandArgs(workspace, mass, dirName, blinded, signalScaleFactor, rMin, rMax)
    # cmd = "combine " + cmdArgs
    taskName = 'asymptotic.M{}'.format(mass)
    if signalScaleFactor != 1.0:
        taskName += '.signalScaleFactor{}'.format(round(signalScaleFactor, 6))
        if index != -1:
            taskName += '.betaIndex{}'.format(index)
    queue = "espresso"
    inputFiles = [str(workspace)]
    combineBin = shutil.which("combine")
    if combineBin is None:
        raise RuntimeError("No combine binary found in $PATH. Did you compile combine?")
    inputFiles.append(combineBin)
    SubmitLimitJobsBatch(mass, dirName, listFailedCommands, -1, signalScaleFactor, taskName, cmdArgs, "", queue, inputFiles)


def RunHybridNewInteractive(workspace, mass, dirName, quantile=-1, genAsimovToyFile="", signalScaleFactor=1.0, rMin=0, rMax=200):
    cmds = GetHybridNewCommandArgs(workspace, mass, dirName, quantile, genAsimovToyFile, signalScaleFactor, False)
    resultFiles = []
    for cmd in cmds:
        cmd = 'combine ' + cmd
        RunCommand(cmd, dirName)
        quantileString = "quant{:.3f}".format(quantile) if quantile > 0 else "."
        globString = dirName+'/higgsCombineTest.HybridNew.mH{}.*.{}.root'.format(mass, quantileString)
        resultFiles.append(FindFile(globString))
    return resultFiles


def RunAsymptoticInteractive(workspace, mass, dirName, batch, blinded=True, signalScaleFactor=1.0, rMin=-1, rMax=150):
    cmd = GetAsymptoticCommandArgs(workspace, mass, dirName, blinded, signalScaleFactor, rMin, rMax)
    cmd = 'combine ' + cmd
    RunCommand(cmd, dirName)
    return Path(sorted(glob.glob(dirName+'/higgsCombineTest.AsymptoticLimits.mH{}.*.root'.format(mass)), key=os.path.getmtime)[-1]).resolve()


def GetComputeLimitsFromGridCommand(workspace, mass, filename, quantiles, signalScaleFact, rMin=0, rMax=1000, toyFile=""):
    rAbsAcc = 0.00001
    rRelAcc = 0.005
    if signalScaleFact != -1:
        sigSFStr = "signalScaleFactor"+str(signalScaleFact)+"_"
    else:
        sigSFStr = ""
    cmds = []
    if not isinstance(quantiles, list):
        quantiles = [quantiles]
    for quantile in quantiles:
        # if rMin==-1 and rMax==-1:
        #     rMin, rMax = GetRMinAndRMax(mass, quantile, signalScaleFact)
        plotFilename = "limit_scan_{}m{}_quant{}.pdf".format("signalScaleFactor"+str(signalScaleFact)+"_" if signalScaleFact != -1 else "", mass, quantile)
        cmd = './combine ' if doBatch else 'combine '
        if options.cmsConnectMode:
            cmd += ' $condorNodeBaseDir/{}'.format(str(workspace).split("/")[-1])
        else:
            cmd += '{}'.format(str(workspace).split("/")[-1] if doBatch else workspace)
        # cmd += ' -v1'
        cmd += ' -M HybridNew'
        if options.doSignificance:
            cmd += ' --LHCmode LHC-significance'
            if len(filename):
                cmd += ' --toysFile={}'.format(filename)
            else:
                cmd += ' --toysFile={}'
            cmd += ' --readHybridResult'
        else:
            cmd += ' --rMin {}'.format(rMin)
            cmd += ' --rMax {}'.format(rMax)
            cmd += ' --LHCmode LHC-limits'
            if len(filename):
                cmd += ' --grid={}'.format(filename)
            else:
                cmd += ' --grid={}'
            cmd += ' --rAbsAcc {}'.format(rAbsAcc)
            cmd += ' --rRelAcc {}'.format(rRelAcc)
            cmd += ' --plot={}'.format(plotFilename)
            cmd += ' --readHybridResults'
        if quantile > 0:
            cmd += ' --expectedFromGrid={}'.format(quantile)
        if toyFile!="":
            if options.cmsConnectMode:
                cmd += " -D $condorNodeBaseDir/{}:toys/toy_asimov".format(str(toyFile).split("/")[-1])
            else:
                cmd += " -D {}:toys/toy_asimov".format(str(toyFile).split("/")[-1] if doBatch else toyFile)
        cmd += ' -m {}'.format(mass)
        #cmd += ' --rAbsAcc {}'.format(rAbsAcc)
        #cmd += ' --rRelAcc {}'.format(rRelAcc)
        #cmd += ' --plot={}'.format(plotFilename)
        if signalScaleFact != -1 and not options.doSignificance:
            cmd += ' -n .signalScaleFactor{}'.format(signalScaleFact)
        cmd += commonCombineArgs.format(signalScaleFact)
        cmds.append(cmd)
    return cmds


def ComputeLimitsFromGrid(workspace, mass, dirName, filename, quantile, signalScaleFact=-1):
    cmd = GetComputeLimitsFromGridCommand(workspace, mass, filename, quantile, signalScaleFact)
    RunCommand(cmd, dirName)
    if signalScaleFact != -1:
        return Path(sorted(glob.glob(dirName+'/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}.root'.format(signalScaleFact, mass, "quant{:.3f}".format(quantile))), key=os.path.getmtime)[-1]).resolve()
    else:
        return Path(sorted(glob.glob(dirName+'/higgsCombineTest.HybridNew.mH{}.{}.root'.format(mass, "quant{:.3f}".format(quantile))), key=os.path.getmtime)[-1]).resolve()


def ExtractLimitResult(rootFile):
    # higgsCombineTest.MethodName.mH$MASS.[word$WORD].root
    # higgsCombineTest.HybridNew.mH120.quant0.500.root
    if not os.path.isfile(rootFile) and not "eoscms.cern.ch" in rootFile and not "fnal.gov" in rootFile:
        raise RuntimeError("ERROR: Did not find the root file {}. Exiting.".format(rootFile))
    tfile = TFile.Open(str(rootFile))
    limitTree = tfile.Get("limit")
    bytesRead = limitTree.GetEntry(0)
    if bytesRead <= 0:
        raise RuntimeError("ERROR: Something went wrong: read {} bytes from 'limit' tree in file {}. Exiting.".format(bytesRead, tfile))
    limit = limitTree.limit
    limitErr = limitTree.limitErr
    quantile = limitTree.quantileExpected
    signalScaleParam = 0
    if "trackedParam_signalScaleParam" in limitTree.GetListOfBranches():
        # print("Scale original limit r={}+/-{} by signalScaleParam={}; final limit={}+/-{}".format(limit, limitErr, limitTree.trackedParam_signalScaleParam, limit/limitTree.trackedParam_signalScaleParam, limitErr/limitTree.trackedParam_signalScaleParam))
        # limit /= limitTree.trackedParam_signalScaleParam
        # limitErr /= limitTree.trackedParam_signalScaleParam
        signalScaleParam = limitTree.trackedParam_signalScaleParam
    tfile.Close()
    return limit, limitErr, quantile, signalScaleParam


def ExtractAsymptoticLimitResultRoot(rootFile):
    if not os.path.isfile(rootFile):
        raise RuntimeError("ERROR: Did not find the root file {}. Exiting.".format(rootFile))
    tfile = TFile.Open(str(rootFile))
    limitTree = tfile.Get("limit")
    bytesRead = limitTree.GetEntry(0)
    if bytesRead <= 0:
        raise RuntimeError("ERROR: Something went wrong: read {} bytes from 'limit' tree in file {}. Exiting.".format(bytesRead, tfile))
    if "trackedParam_signalScaleParam" in limitTree.GetListOfBranches():
        # print("Scale original limit r={}+/-{} by signalScaleParam={}; final limit={}+/-{}".format(limit, limitErr, limitTree.trackedParam_signalScaleParam, limit/limitTree.trackedParam_signalScaleParam, limitErr/limitTree.trackedParam_signalScaleParam))
        signalScaleParam = limitTree.trackedParam_signalScaleParam
    else:
        signalScaleParam = 1.0
    limits = []
    limitErrs = []
    quantiles = []
    for iEntry in range(0, limitTree.GetEntries()):
        limitTree.GetEntry(iEntry)
        limits.append(limitTree.limit)
        limitErrs.append(limitTree.limitErr)
        quantiles.append(round(limitTree.quantileExpected, 3))
    tfile.Close()
    return limits, limitErrs, quantiles, signalScaleParam


def ExtractAsymptoticLimitResultTxt(txtFile):
    # NB: limit errors are never in the output txt, so cannot extract them here
    limits = []
    quantiles = []
    with open(txtFile, "r") as theFile:
        for line in theFile:
            if "Set Default Value of Parameter signalScaleParam To :" in line:
                signalScaleParam = float(line.split(":")[-1].strip())
            elif "Expected  2.5%:" in line:
                quantiles.append(0.025)
                limits.append(float(line.split("<")[-1].strip()))
            elif "Expected 16.0%:" in line:
                quantiles.append(0.16)
                limits.append(float(line.split("<")[-1].strip()))
            elif "Expected 50.0%:" in line:
                quantiles.append(0.5)
                limits.append(float(line.split("<")[-1].strip()))
            elif "Expected 84.0%:" in line:
                quantiles.append(0.84)
                limits.append(float(line.split("<")[-1].strip()))
            elif "Expected 97.5%:" in line:
                quantiles.append(0.975)
                limits.append(float(line.split("<")[-1].strip()))
    return limits, [], quantiles, signalScaleParam


def ReadXSecFile(filename):
    xsThByMass = {}
    yPDFupByMass = {}
    yPDFdownByMass = {}
    with open(os.path.expandvars(filename), "r") as xsecFile:
        for line in xsecFile:
            line = line.strip()
            if line.startswith("#"):
                continue
            split = line.split()
            mass = float(split[0])
            xs = float(split[1])
            xsThByMass[mass] =  xs
            # these are in any case in a different order and unit for the vector files vs. the scalar files,
            # and for the stop pair prod, we don't have these.
            # are these used anyway?
            #
            # if len(split)== 7:
            #     yPDFupByMass[mass] = xs*(1+float(split[5])/100.)
            #     yPDFdownByMass[mass] = xs*(1-float(split[6])/100.)
    return xsThByMass, yPDFupByMass, yPDFdownByMass


def CreateArraysForPlotting(xsecLimitsByMassAndQuantile):
    massList = list(xsecLimitsByMassAndQuantile.keys())
    if type(massList[0]) is str:
        for i,m in enumerate(massList):
            massList[i] = int(m)
        massList.sort()
        for i,m in enumerate(massList):
            massList[i] = str(m)
    else:
        massList.sort()
    shadeMassList = []
    xs_medExpList = []
    xs_oneSigmaExpList = []
    xs_twoSigmaExpList = []
    xs_obsList = []
    if str(-1) in xsecLimitsByMassAndQuantile[list(xsecLimitsByMassAndQuantile.keys())[0]].keys():
        hasObserved = True
    else:
        hasObserved = False
    for mass in massList:
        try:
            xs_medExpList.append(xsecLimitsByMassAndQuantile[mass][str(0.5)])
        except KeyError:
            xs_medExpList.append(xsecLimitsByMassAndQuantile[mass]["0.500"])
        if hasObserved:
            xs_obsList.append(xsecLimitsByMassAndQuantile[mass][str(-1)])
        try:
            xs_oneSigmaExpList.append(xsecLimitsByMassAndQuantile[mass][str(0.16)])
        except KeyError:
            xs_oneSigmaExpList.append(xsecLimitsByMassAndQuantile[mass]["0.160"])
        xs_twoSigmaExpList.append(xsecLimitsByMassAndQuantile[mass][str(0.025)])
    for mass in reversed(massList):
        try:
            xs_oneSigmaExpList.append(xsecLimitsByMassAndQuantile[mass][str(0.84)])
        except KeyError:
            xs_oneSigmaExpList.append(xsecLimitsByMassAndQuantile[mass]["0.840"])
        xs_twoSigmaExpList.append(xsecLimitsByMassAndQuantile[mass][str(0.975)])
    masses = np.array(massList, dtype="f")
    shadeMasses = np.concatenate([masses, np.flip(masses)])
    xsMedExp = np.array(xs_medExpList, dtype="f")
    xsObs = np.array(xs_obsList, dtype="f")
    xsOneSigmaExp = np.array(xs_oneSigmaExpList, dtype="f")
    xsTwoSigmaExp = np.array(xs_twoSigmaExpList, dtype="f")
    return masses, shadeMasses, xsMedExp, xsOneSigmaExp, xsTwoSigmaExp, xsObs


def GetBetaRangeToAttempt(mass):
    # get a restricted beta range to test for each mass point
    minBetasPerMass = {300: 0,   400: 0,   500: 0.025, 600: 0.025, 700: 0.035, 800: 0.04, 900: 0.06, 1000: 0.075,  1100: 0.085, 1200: 0.095, 1300: 0.105, 1400: 0.115, 1500: 0.125, 1600: 0.125, 1700: 0.5, 1800: 0.7}
    maxBetasPerMass = {300: 0.075, 400: 0.1, 500: 0.15, 600: 0.2, 700: 0.25, 800: 0.35, 900: 0.45, 1000: 0.6, 1100: 0.7, 1200: 0.8, 1300: 0.9, 1400: 1.0, 1500: 1.0, 1600: 1.0, 1700: 1.0}
    if mass > 1700:
        return 0.9, 1.0
    else:
        return minBetasPerMass[mass], maxBetasPerMass[mass]


def GetBetasToSubmit(mass):
    minBeta, maxBeta = GetBetaRangeToAttempt(int(mass))
    minBetaPosition = bisect(betasToScan, minBeta)
    maxBetaPosition = bisect(betasToScan, maxBeta)
    print("INFO: for mass {}, scan beta range {} to {} or indices {} to {}".format(mass, minBeta, maxBeta, minBetaPosition, maxBetaPosition))
    return betasToScan[minBetaPosition : maxBetaPosition]


def CheckErrorFile(errFileName, throwException=True):
    print("\tChecking error file {}...".format(errFileName), end = "", flush=True)
    messagesToIgnore = ["[WARNING] Minimization finished with status 1 (covariance matrix forced positive definite), this could indicate a problem with the minimum!"]
    messagesToIgnore.append("Info in <TCanvas::Print>: pdf file limit_scan")
    messagesToIgnore.append("WARNING: Environment variable XrdSecGSISRVNAMES already has value")
    messagesToIgnore.append("WARNING: While bind mounting '/etc/hosts:/etc/hosts': destination is already in the mount point list")
    messagesToIgnore.append("WARNING: In non-interactive mode release checks e.g. deprecated releases, production architectures are disabled.")
    messagesToIgnore.append("write error: Broken pipe")
    messagesToIgnore.append("INFO")
    messagesToIgnore.append("DEBUG")
    typesOfIgnoredMessages = ["minimization status 1 warning", "limit scan pdf created","CMS connect env var warning","CMS connect mount warning", "CMS connect deprecated release warning","CMS connect broken pipe error","INFO message ended up in .err file","DEBUG message"]
    showIgnoreMessages = [False, False, False, False, False, False, False, False]
    instancesOfIgnoredMessages = [0]*len(messagesToIgnore)
    with open(errFileName, "r") as errFile:
        for line in errFile:
            # if any(msg in line for msg in messagesToIgnore):
            if len(line):
                ignoreLine = False
                for idx, msg in enumerate(messagesToIgnore):
                    if msg in line:
                        instancesOfIgnoredMessages[idx] += 1
                        ignoreLine = True
                if not ignoreLine:
                    #print(colored("Found unexpected content '{}' in {}".format(line, errFileName), "red"))
                    if throwException:
                        raise RuntimeError("Found unexpected content '{}' in {}".format(line, errFileName))
    for idx, count in enumerate(instancesOfIgnoredMessages):
        if count > 0 and showIgnoreMessages[idx]:
            print(colored("Found {} error message(s) of type {} in {}".format(count, typesOfIgnoredMessages[idx], errFileName), "red"))
    # print("OK")


def CheckForErrorAndResultFiles(massList, condorDir):
    # eosDirNoPrefix = eosDir[eosDir.rfind("//")+1:]
    # eosDirName = eosDirNoPrefix.rstrip("/") + "/" + dirName
    errFilesByMassAndQuantile = {}
    rootFilesByMassAndQuantile = {}
    for mass in massList:
        errFilesByMassAndQuantile[mass] = {}
        rootFilesByMassAndQuantile[mass] = {}
        for quantileExp in quantilesExpected:
            signalScaleFactor = 1.0
            if rescaleSignal:
                signalScaleFactor = round(xsRescaleByMass[mass]/xsThByMass[mass], 6)  # we assume that the signal yields as given have been normalized with the scalar pair-prod xsec
            quantile = "quant{:.3f}".format(quantileExp) if quantileExp > 0 else "."
            quantileStr = str(quantileExp).replace(".", "p")
            #globStrings = [condorDir+'/error/hybridNewLimits.M{}.{}.*.0.err'.format(mass, quantileStr)]
            #globStrings.append(condorDir+'/gridGen/hybridNewGridGenAndLimits.M{}.{}/condor.*.0.err'.format(mass, quantileStr))
            #globStrings.append(condorDir+'/limits/hybridNewLimits.M{}.{}.signalScaleFactor{}/condor.*.0.err'.format(mass, quantileStr,signalScaleFactor))
            #if signalScaleFactor != 1.0:
            #    globStrings.append(eosDirName+'/hybridNewLimits.M{}.{}.signalScaleFactor{}/condor.*.0.err'.format(mass, quantileStr, signalScaleFactor))
            #else:
            #    globStrings.append(eosDirName+'/hybridNewLimits.M{}.{}/condor.*.0.err'.format(mass, quantileStr))
            globStrings = [eosDirName+'/hybridNewLimits.M{}.{}.signalScaleFactor*/condor.*.0.err'.format(mass, quantileStr)]
            errFileName = FindFile(globStrings)
            errFilesByMassAndQuantile[mass][quantileExp] = errFileName
            globStrings = [condorDir+'/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}.root'.format(signalScaleFactor, mass, quantile)]
            globStrings.append(condorDir+'/gridGen/hybridNewGridGenAndLimits.M{}.{}/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}.root'.format(mass, quantileStr,signalScaleFactor, mass, quantile))
            if signalScaleFactor != 1.0:
                if not quantile == ".":
                    globStrings.append(eosDirName+'/hybridNewLimits.M{}.{}.signalScaleFactor{}/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}.root'.format(mass, quantileStr, signalScaleFactor, signalScaleFactor, mass, quantile))
                else:
                    globStrings.append(eosDirName+"/hybridNewLimits.M{}.{}.signalScaleFactor{}/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.root".format(mass, quantileStr, signalScaleFactor, signalScaleFactor, mass))
            else:
                globStrings.append(eosDirName+'/hybridNewLimits.M{}.{}/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}.root'.format(mass, quantileStr, signalScaleFactor, mass, quantile))
            globStrings.append(eosDirName+'/hybridNewLimits.M{}.{}.signalScaleFactor*/higgsCombine.signalScaleFactor*.HybridNew.mH{}.{}.root'.format(mass,quantileStr,mass, quantile))
            rootFileName = FindFile(globStrings)
            rootFilesByMassAndQuantile[mass][quantileExp] = rootFileName
    return errFilesByMassAndQuantile, rootFilesByMassAndQuantile


def ReadBatchResults(massList, condorDir):
    errFilesByMassAndQuantile, rootFilesByMassAndQuantile = CheckForErrorAndResultFiles(massList, condorDir)
    for mass in massList:
        workspaceGlobString = dirName + "/datacards/*.m{}.root".format(mass)
        listOfWorkspaceFiles = sorted(glob.glob(workspaceGlobString), key=os.path.getmtime)
        if len(listOfWorkspaceFiles) != 1:
            raise RuntimeError("Globbing for {} did not result in one file as expected, but {}".format(workspaceGlobString, listOfWorkspaceFiles))
        rLimitsByMassAndQuantile[mass] = manager.dict()
        xsecLimitsByMassAndQuantile[mass] = manager.dict()
        signalScaleFactorsByMassAndQuantile[mass] = manager.dict()
        cardWorkspace = Path(listOfWorkspaceFiles[-1]).resolve()
        for quantileExp in quantilesExpected:
            signalScaleFactorsByMassAndQuantile[mass][str(quantileExp)] = {}
            # quantile = "quant{:.3f}".format(quantileExp) if quantileExp > 0 else "."
            # quantileStr = str(quantileExp).replace(".", "p")
            # globString = condorDir+'/error/hybridNewLimits.M{}.{}.*.0.err'.format(mass, quantileStr)
            # errFileName = None
            # try:
            #     errFileName = FindFile(globString)
            # except RuntimeError as e:
            #     try:
            #         globString2 = condorDir+'/gridGen/hybridNewGridGenAndLimits.M{}.{}/condor.*.0.err'.format(mass, quantileStr)
            #         errFileName = FindFile(globString2)
            #     except RuntimeError as e2:
            #         print("Caught exceptions looking for error files: '{}' and '{}'".format(e, e2))
            errFileName = errFilesByMassAndQuantile[mass][quantileExp]
            CheckErrorFile(errFileName, False)  # don't throw exception
            # # now run HybridNew again to get limits from grid
            # globString = condorDir+'/higgsCombine.signalScaleFactor*.HybridNew.mH{}.quant{}.grid.root'.format(mass, quantileExp)
            # rootFileName = FindFile(globString)
            # resultFileName = ComputeLimitsFromGrid(cardWorkspace, mass, condorDir, rootFileName, quantileExp)
            # # now extract limits from the root file produced above
            # limit, limitErr, quantileFromFile, signalScaleFactor = ExtractLimitResult(resultFileName)
            # signalScaleFactor = 1.0
            # if rescaleSignal:
            #     signalScaleFactor = round(xsRescaleByMass[mass]/xsThByMass[mass], 6)  # we assume that the signal yields as given have been normalized with the scalar pair-prod xsec
            # globString = condorDir+'/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}.root'.format(signalScaleFactor, mass, quantile)
            # rootFileName = None
            # try:
            #     rootFileName = FindFile(globString)
            # except RuntimeError as e:
            #     try:
            #         globString2 = condorDir+'/gridGen/hybridNewGridGenAndLimits.M{}.{}/higgsCombine.signalScaleFactor1.0.HybridNew.mH{}.{}.root'.format(mass, quantileStr, mass ,quantile)
            #         rootFileName = FindFile(globString2)
            #     except RuntimeError as e2:
            #         print("Caught exceptions looking for root files: '{}' and '{}'".format(e, e2))
            rootFileName = rootFilesByMassAndQuantile[mass][quantileExp]
            limit, limitErr, quantileFromFile, signalScaleFactor = ExtractLimitResult(rootFileName)
            rLimitsByMassAndQuantile[mass][str(quantileExp)] = limit
            #print(xsThByMass[float(mass)])
            #print(xsRescaleByMass[mass])
            xsecLimitsByMassAndQuantile[mass][str(quantileExp)] = limit * signalScaleFactor * (xsThByMass[float(mass)])# if signalScaleFactor == 1.0 else xsRescaleByMass[mass])
            signalScaleFactorsByMassAndQuantile[mass][str(quantileExp)] = signalScaleFactor
    return xsecLimitsByMassAndQuantile, signalScaleFactorsByMassAndQuantile


def ReadBatchResultsBetaScanSerial(massList, condorDir):
    for quantileExp in quantilesExpected:
        rLimitsByMassAndQuantile[str(quantileExp)] = {}
        xsecLimitsByMassAndQuantile[str(quantileExp)] = {}
        signalScaleFactorsByMassAndQuantile[str(quantileExp)] = {}
        for mass in massList:
            listOfWorkspaceFiles = sorted(glob.glob(dirName + "/datacards/*.m{}.root".format(mass)), key=os.path.getmtime)
            if len(listOfWorkspaceFiles) != 1:
                raise RuntimeError("Globbing for {} did not result in one file as expected, but {}".format(dirName + "/datacards/*.m{}.root".format(mass), listOfWorkspaceFiles))
            cardWorkspace = Path(listOfWorkspaceFiles[-1]).resolve()
            quantile = "quant{:.3f}".format(quantileExp) if quantileExp > 0 else "."
            quantile = str(quantile).rstrip("0")
            quantileStr = str(quantileExp).replace(".", "p").rstrip("0")
            globString = condorDir+'/error/hybridNewLimits.M{}.{}.*.0.err'.format(mass, quantileStr)
            errorFiles = GetFileList(globString)
            for errFile in errorFiles:
                CheckErrorFile(errFile)
                lastPart = errFile.split("{}.".format(quantileStr))[-1]
                lastPartCut = lastPart.rstrip(".0.err")
                sigScaleFact = lastPartCut[0:lastPartCut.rfind(".")].strip("signalScaleFactor")
                sigScaleFactRound = round(float(sigScaleFact), 6)
                rootGlobString = condorDir+'/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}.grid.root'.format(sigScaleFactRound, mass, quantile)
                rootFileName = FindFile(rootGlobString)
                resultFileName = ComputeLimitsFromGrid(cardWorkspace, mass, condorDir, rootFileName, quantileExp, sigScaleFact)
                limit, limitErr, quantileFromFile, signalScaleFactor = ExtractLimitResult(resultFileName)
                betaVal = math.sqrt(signalScaleFactor)
                betaValRound = str(round(betaVal, 6))
                if not betaValRound in rLimitsByMassAndQuantile[str(quantileExp)].keys():
                    rLimitsByMassAndQuantile[str(quantileExp)][betaValRound] = {}
                    xsecLimitsByMassAndQuantile[str(quantileExp)][betaValRound] = {}
                    signalScaleFactorsByMassAndQuantile[str(quantileExp)][betaValRound] = {}
                rLimitsByMassAndQuantile[str(quantileExp)][betaValRound][mass] = limit
                xsecLimitsByMassAndQuantile[str(quantileExp)][betaValRound][mass] = limit * xsThByMass[float(mass)]
                signalScaleFactorsByMassAndQuantile[str(quantileExp)][betaValRound][mass] = signalScaleFactor
                # print("\tINFO: ReadBatchResultsBetaScan() - fill xsecLimitsByMassAndQuantile[{}][{}][{}]".format(quantileExp, betaValRound, mass))
    return xsecLimitsByMassAndQuantile


def get_value(d):
    return {
        key: get_value(sub_d)
        if isinstance(sub_d, multiprocessing.managers.DictProxy) else sub_d
        for key, sub_d in d.items()
            }


def ReadBatchResultParallel(args):
    mass, errFile, rootFileName, betaValRound, quantileExp = args
    print("ReadBatchResult for mass {}, beta {}".format(mass, betaValRound))
    # quantile = "quant{:.3f}".format(quantileExp) if quantileExp > 0 else "."
    # # quantile = str(quantile).rstrip("0")
    # quantileStr = str(quantileExp).replace(".", "p").rstrip("0")
    # CheckErrorFile(errFile)
    CheckErrorFile(errFile, False)
    # lastPart = errFile.split("{}.".format(quantileStr))[-1]
    # lastPartSplit = lastPart.split("/condor")
    # sigScaleFact = lastPartSplit[0].split("signalScaleFactor")[-1]
    # signalScaleFactor = float(sigScaleFact) # FIXME: here, like previously, we should probably use the beta index instead to keep track of the unrounded scale factor
    # # rootGlobString = condorDir+'/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}.root'.format(sigScaleFactRound, mass, quantile)
    # rootGlobString = str(Path(errFile).parent)+'/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}.root'.format(sigScaleFact, mass, quantile)
    # rootFileName = FindFile(rootGlobString)
    # resultFileName = ComputeLimitsFromGrid(cardWorkspace, mass, condorDir, rootFileName, quantileExp, sigScaleFact)
    limit, limitErr, quantileFromFile, signalScaleFactor = ExtractLimitResult(rootFileName)
    print("limit = {} for mass {}, beta {}".format(limit, mass, betaValRound))
    # betaVal = math.sqrt(signalScaleFactor)
    # betaVal = math.sqrt(float(sigScaleFact))
    # betaValRound = str(round(betaVal, 6))
    # print()
    # print("INFO: [3] ReadBatchResultParallel(): store limit*xs[{}]={}*{}={} in rLimitsByMassAndQuantile={} for key={} at betaValRound={} from rootFileName={}".format(
    #     mass, limit, xsThByMass[float(mass)], limit * xsThByMass[float(mass)], get_value(rLimitsByMassAndQuantile), str(quantileExp), betaValRound, rootFileName), flush=True)
    rLimitsByMassAndQuantile[str(quantileExp)][betaValRound][mass] = limit
    xsecLimitsByMassAndQuantile[str(quantileExp)][betaValRound][mass] = limit * xsThByMass[float(mass)]
    signalScaleFactorsByMassAndQuantile[str(quantileExp)][betaValRound][mass] = signalScaleFactor
    # print("INFO: [5] ReadBatchResultParallel(): DONE storing xsecLimitsByMassAndQuantile[{}]={}".format(str(quantileExp), get_value(xsecLimitsByMassAndQuantile)), flush=True)


def ReadBatchResultsBetaScan(massList, condorDir):
    with Progress() as progress:
        nJobs = 0
        for quantileExp in quantilesExpected:
            for mass in massList:
                quantileStr = str(quantileExp).replace(".", "p").rstrip("0")
                globString = condorDir+'/limits/hybridNewLimits.M{}.{}/signalScaleFactor*/condor*.0.err'
                #globString = eosDirName.rstrip("/")+"/hybridNewLimits.M{}.{}.*/condor*.0.err"
                try:
                    errorFiles = GetFileList(globString.format(mass, quantileStr))
                    nJobs += len(errorFiles)
                except RuntimeError:
                    globString2 = condorDir+'/gridGen/hybridNewGridGenAndLimits.M{}*/condor*.0.err'
                    try:
                        errorFiles = GetFileList(globString2.format(mass))
                        nJobs += 5 * len(errorFiles)
                        globString = globString2
                    except RuntimeError as e:
                        raise RuntimeError("Caught exception: "+str(e)+"; this is after already trying and failing to glob for {}".format(globString))
                # print("Got {} error files for mass {}, quantile {}".format(len(errorFiles), mass, quantileExp))
        # print("INFO: reading {} jobs".format(nJobs), flush=True)
        task_id = progress.add_task("[cyan]Reading batch results...", total=nJobs)
        with multiprocessing.Pool(ncores) as pool:
            partial_quit = partial(ErrorCallback, pool)
            for quantileExp in quantilesExpected:
                rLimitsByMassAndQuantile[str(quantileExp)] = manager.dict()
                xsecLimitsByMassAndQuantile[str(quantileExp)] = manager.dict()
                signalScaleFactorsByMassAndQuantile[str(quantileExp)] = manager.dict()
                for mass in massList:
                    listOfWorkspaceFiles = sorted(glob.glob(dirName + "/datacards/*.m{}.root".format(mass)), key=os.path.getmtime)
                    if len(listOfWorkspaceFiles) != 1:
                        raise RuntimeError("Globbing for {} did not result in one file as expected, but {}".format(workspaceFileName, listOfWorkspaceFiles))
                    cardWorkspace = Path(listOfWorkspaceFiles[-1]).resolve()
                    quantileStr = str(quantileExp).replace(".", "p").rstrip("0")
                    quantile = "quant{:.3f}".format(quantileExp) if quantileExp > 0 else "."
                    errorFiles = GetFileList(globString.format(mass, quantileStr))
                    for errFile in errorFiles:
                        if "GridGenAndLimits" in errFile:
                            taskName = Path(errFile).parent.name
                            sigScaleFact = taskName.split("signalScaleFactor")[-1]
                            signalScaleFactor = float(sigScaleFact) # FIXME: here, like previously, we should probably use the beta index instead to keep track of the unrounded scale factor
                            eosPath = eosDir[eosDir.rfind("//")+1:]  # remove root://.../ prefix
                            rootGlobString = eosPath+"/"+dirName+"/"+taskName+'/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}.root'.format(sigScaleFact, mass, quantile)
                            rootFileName = FindFile(rootGlobString)
                        else:
                            lastPart = errFile.split("{}.".format(quantileStr))[-1]
                            lastPartSplit = lastPart.split("/condor")
                            sigScaleFact = lastPartSplit[0].split("signalScaleFactor")[-1]
                            signalScaleFactor = float(sigScaleFact) # FIXME: here, like previously, we should probably use the beta index instead to keep track of the unrounded scale factor
                            #rootGlobString = str(Path(errFile).parent)+'/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}.root'.format(sigScaleFact, mass, quantile)
                            outFile = errFile.split("/condor")[0]
                            #rootFileName = FindFile(rootGlobString)
                            rootFileName = eosDir + "/" + options.name + "/hybridNewLimits.M{}.{}/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}.root".format(mass,quantileStr,sigScaleFact,mass,quantile)#GetFilenamesFromCondorOutFiles(outFile)
                        betaVal = math.sqrt(signalScaleFactor)
                        betaValRound = str(round(betaVal, 6))
                        if not betaValRound in rLimitsByMassAndQuantile[str(quantileExp)].keys():
                            rLimitsByMassAndQuantile[str(quantileExp)][betaValRound] = manager.dict()
                        if not betaValRound in xsecLimitsByMassAndQuantile[str(quantileExp)].keys():
                            xsecLimitsByMassAndQuantile[str(quantileExp)][betaValRound] = manager.dict()
                        if not betaValRound in signalScaleFactorsByMassAndQuantile[str(quantileExp)].keys():
                            signalScaleFactorsByMassAndQuantile[str(quantileExp)][betaValRound] = manager.dict()
                        # print()
                        # print("INFO: [3] ReadBatchResultsBetaScan(): store limit*xs[{}] in rLimitsByMassAndQuantile={} for key={} at betaValRound={} from errFile={}".format(
                        #     mass, get_value(rLimitsByMassAndQuantile), str(quantileExp), betaValRound, errFile), flush=True)
                        # rLimitsByMassAndQuantile[str(quantileExp)][betaValRound][mass] = limit
                        # xsecLimitsByMassAndQuantile[str(quantileExp)][betaValRound][mass] = limit * xsThByMass[float(mass)]
                        # signalScaleFactorsByMassAndQuantile[str(quantileExp)][betaValRound][mass] = signalScaleFactor
                        # print("INFO: [5] ReadBatchResultsBetaScan(): DONE storing xsecLimitsByMassAndQuantile[{}]={}".format(str(quantileExp), get_value(xsecLimitsByMassAndQuantile)), flush=True)
                        try:
                            pool.apply_async(ReadBatchResultParallel, [[mass, errFile, rootFileName, betaValRound, quantileExp]], callback = lambda x: progress.advance(task_id), error_callback = partial_quit)
                            # nJobs += 1
                        except KeyboardInterrupt:
                            print("\n\nCtrl-C detected: Bailing.")
                            pool.terminate()
                            sys.exit(1)
                        except Exception as e:
                            print("ERROR: caught exception while queuing ReadBatchResultParallel job with err file {}".format(errFile))
                            # traceback.print_exc()
                            # pool.terminate()
                            # exit(-2)
            # now close the pool and wait for jobs to finish
            pool.close()
            # if nJobs > 0:
            #     print("Waiting for {} limit calculations to finish for all beta values...".format(nJobs))
            pool.join()


def sigmaval(spline, mval):
    return spline.Eval(mval)
    
 
def mval(spline, sigma):
    testm = 150
    oldtestm = 1500
    inc = 50
    dif = 55
    olddif = 0
    while abs(oldtestm - testm) > 0.01:
        testsigma = sigmaval(spline, testm)
        olddif = dif
        dif = testsigma - sigma
        if testm > 1500:
            break
        if dif*olddif <= 0.0:
            inc = -inc/2.3
        oldtestm = testm
        #print '**' + str(testm) + '  ' + str(testsigma) +'  ' +str(dif) + '   ' + str(dif*olddif)
        testm = testm + inc
    return testm


def loggraph(inputarrayX, inputarrayY):
    logarray = []
    for j in inputarrayY:
        logarray.append(math.log10(j))
        # logarray.append(math.log(j))
    x = np.array(inputarrayX, dtype="f")
    y = np.array(logarray, dtype="f")
    g = TGraph(len(x), x, y)
    return g


def logspline(inputarrayX,inputarrayY):
    logarray = []
    for j in inputarrayY:
        logarray.append(math.log(j))
    x = np.array("d",inputarrayX)
    y = np.array("d",logarray)
    g = TGraph(len(x),x,y)
    outspline = TSpline3("",g)
    return outspline


def get_simple_intersection(graph1, graph2, xmin, xmax, verbose=False):
    num = (xmax-xmin)*10
    inc = (xmax - xmin)/(1.0*num)
    dif = []
    sdif = []
    x = xmin +0.1
    xvals = []
    xx = []
    yy = []
    xvals = []
    while x < (xmax-.1):
        # print("\tINFO: get_simple_intersection() - x=", str(x) + '   xmax-.1='+ str(xmax-.1))
        # print("\tINFO: get_simple_intersection() - graph1.Eval(x) =", graph1.Eval(x))
        # print("\tINFO: get_simple_intersection() - math.pow(graph1.Eval(x), 10) =", math.pow(graph1.Eval(x), 10))
        # print("\tINFO: get_simple_intersection() - graph2.Eval(x) =", graph2.Eval(x))
        # print("\tINFO: get_simple_intersection() - math.pow(graph2.Eval(x), 10) =", math.pow(graph2.Eval(x), 10))
        # print("\tINFO: get_simple_intersection() - math.pow(graph1.Eval(x), 10) - math.pow(graph2.Eval(x), 10) =", math.pow(graph1.Eval(x), 10) - math.pow(graph2.Eval(x), 10))
        #
        # thisdif = (math.exp(graph1.Eval(x)) - math.exp(graph2.Eval(x)))
        # xx.append(math.exp(graph1.Eval(x)))
        # yy.append(math.exp(graph2.Eval(x)))
        #XXX SIC TEST Jun. 26 2024: use 10^x instead of exp(x)
        thisdif = (math.pow(10, graph1.Eval(x)) - math.pow(10, graph2.Eval(x)))
        xx.append(math.pow(10, graph1.Eval(x)))
        yy.append(math.pow(10, graph2.Eval(x)))
        #
        sdif.append(thisdif)
        dif.append(abs(thisdif))
        xvals.append(x)
        if verbose:
            print("\tINFO: get_simple_intersection() [1] - mass=", str(x) + '   th=' +str(math.pow(10, graph1.Eval(x))) + '    lim='+str(math.pow(10, graph2.Eval(x))) + '    th-lim=' + str(thisdif))
        x = x+inc
	#print 'Done Looping for Difs'
    mindif = min(dif)
    bestmass = 0	
    
    for x in range(len(dif)-2):
       a = sdif[x]
       b = sdif[x+1]
       if verbose:
           print("\tINFO: get_simple_intersection() [2] -", str(xvals[x+1]) +'    '+str(a)  + '     ' +str(b))
       if ((a/abs(a))*(b/abs(b))) < 0.0 and a > 0.0 :
           if verbose:
               print('\tINFO: get_simple_intersection() [3] - Limit found at: '+ (str(xvals[x])))
           bestmass = xvals[x]
           break
    return [bestmass, mindif]


def ComputeBetaLimits(xsThByMass, xsecLimitsByMassAndQuantile):
    mTh = np.array(list(xsThByMass.keys()), dtype="f")
    xsTh = np.array(list(xsThByMass.items()), dtype="f")
    print("INFO: Theory masses = ",mTh)
    print("INFO: Theory xs = ",xsTh)
    g = TGraph(len(mTh), mTh, xsTh);
    # spline = TSpline3("xsection", g)
    logtheory = loggraph(list(xsThByMass.keys()), list(xsThByMass.values()))  # make graph with log10 of y values
    massLimitsByQuantileAndBetaVal = {}
    # print("ComputeBetaLimits: xsecLimitsByMassAndQuantile[0.5]={}".format(get_value(xsecLimitsByMassAndQuantile["0.5"])))
    for quantile in xsecLimitsByMassAndQuantile.keys():
        with Progress() as progress:
            verbose = True
            if quantile == "0.5":
                verbose = True
            if verbose:
                print("INFO: ComputeBetaLimits() - examine quantile={}".format(quantile))
            else:
                task_id = progress.add_task("[cyan]Computing beta limits for quantile={}...".format(quantile), total=len(xsecLimitsByMassAndQuantile[quantile].keys()))
            massLimitsByQuantileAndBetaVal[quantile] = {}
            betaIds = list(xsecLimitsByMassAndQuantile[quantile].keys())
            for i,b in enumerate(betaIds):
                betaIds[i] = int(b)
            for betaId in sorted(betaIds):
                betasToScan = list(np.linspace(0.0, 1, 500))[:-1] + [0.9995]
                betaVal = betasToScan[int(betaId)]
                limit_set = list(xsecLimitsByMassAndQuantile[quantile][str(betaId)].values())
                massList = [int(mass) for mass in xsecLimitsByMassAndQuantile[quantile][str(betaId)].keys()]
                if verbose:
                    print("\tINFO: examine betaId={}, value={}".format(betaId,betaVal))
                    print("\tINFO: examine massList={}, limit_set={}".format(massList, limit_set))
                    if len(massList) < 2:
                        print("\tWARN: Skipping beta value={} as we only have one mass point tested here({})! Need to adjust beta scan range.".format(betaVal, massList))
                        continue
                # print("\tINFO: make loggraph for fitted_limits")
                try:
                    fitted_limits = loggraph(massList, limit_set)
                except TypeError as e:
                    print("ERROR: Caught TypeError in loggraph() for quantile={}, betaVal={}, with massList={}, limit_set={}: {}".format(quantile, betaVal, massList, limit_set, e))
                    raise e
                #FIXME: get_simple_intersection seems to assume both graphs are in log base e rather than log10, as they actually are
                if verbose:
                    print("\tINFO: ComputeBetaLimits() - get_simple_intersection")
                goodm = get_simple_intersection(logtheory, fitted_limits, min(mTh), max(mTh)) # , verbose and float(betaVal) > 0.05 and float(betaVal) < 0.06)
                if verbose:
                    print("\tINFO: for betaVal={}: bestmass, mindif={}".format(betaVal, goodm))
                if goodm[0] < min(massList) or goodm[0] > max(massList):
                    if verbose:
                        print("\tWARN: For beta value={}, intersection mass {} is outside of massList range {}! Need to adjust beta scan range.".format(betaVal, goodm[0],  massList))
                else:
                    massLimitsByQuantileAndBetaVal[quantile][str(betaId)] = round(goodm[0], 3)
                if not verbose:
                    progress.advance(task_id)
            # print("ComputeBetaLimits: finished quantile={}".format(quantile), flush=True)
    return massLimitsByQuantileAndBetaVal



def CreateComboArrays(xsecLimitsByMassAndQuantile):
    betasToScan = list(np.linspace(0.0, 1, 500))[:-1] + [0.9995]
    retVal = {}
    for quantile in xsecLimitsByMassAndQuantile.keys():
        retVal[quantile] = {}
        retVal[quantile]["betas"] = []
        retVal[quantile]["massLimits"] = []
        sortedBetas = sorted([int(beta) for beta in xsecLimitsByMassAndQuantile[quantile].keys()])
        # for betaVal, mass in xsecLimitsByMassAndQuantile[quantile].items():
        for betaId in sortedBetas:
            betaVal = betasToScan[int(betaId)]
            retVal[quantile]["betas"].append(betaVal)
            if str(betaId) not in xsecLimitsByMassAndQuantile[quantile].keys():
                print("Oops! now trying to get out {} from array which does not have it in its keys: {}".format(str(betaVal), xsecLimitsByMassAndQuantile[quantile].keys()))
            retVal[quantile]["massLimits"].append(float(xsecLimitsByMassAndQuantile[quantile][str(betaId)]))
    return retVal


def MakeResultTable(masses, xsMedExp, xsObs, xsOneSigmaExp, xsTwoSigmaExp):
    doObs = False
    if len(xsObs):
        doObs = True
    table = []
    for idx, mass in enumerate(list(masses)):
        tableLine = [int(mass), xsMedExp[idx], xsOneSigmaExp[idx], xsTwoSigmaExp[idx]]
        if doObs:
            tableLine.append(xsObs[idx])
        table.append(tableLine)
    columnNames = ["MLQ", "xs exp [pb]", "xs 1sigma exp [pb]", "xs 2sigma exp [pb]"]
    if doObs:
        columnNames.append("xs obs [pb]")
    return table, columnNames


def MakeResultTableBetaScan(masses, dictByMassAndQuantile, quantileExp="0.5"):
    table = []
    headers = ["MLQ / beta="]
    betaValsSorted = sorted(dictByMassAndQuantile[quantileExp].keys())
    for mIdx, mass in enumerate(masses):
        tableRow = [mass]
        for betaVal in betaValsSorted:
            if mIdx == 0:
                headers.append(str(betaVal))
            if str(mass) in dictByMassAndQuantile[quantileExp][betaVal].keys():
                xsecLimit = dictByMassAndQuantile[quantileExp][betaVal][str(mass)]
            elif mass in dictByMassAndQuantile[quantileExp][betaVal].keys():
                xsecLimit = dictByMassAndQuantile[quantileExp][betaVal][mass]
            else:
                xsecLimit = -1  # not all masses have all beta values
            tableRow.append(xsecLimit)
        table.append(tableRow)
    return table, headers


def InitCardsAndWorkspaces(dirName):
    scaleFactorsFile = dirName+'/signalScaleFactors.json'
    if os.path.isfile(scaleFactorsFile):
        print("INFO: Reading scale factors from {}...".format(scaleFactorsFile), flush=True)
        with open(scaleFactorsFile, 'r') as f:
            signalScaleFactorsByMassAndQuantile = json.load(f)
    else:
        print("Did not find scale factors JSON file expected: '{}'; will assume 1.0 for all signal scale factors.".format(scaleFactorsFile))
        signalScaleFactorsByMassAndQuantile = {}
        for mass in massList:
            signalScaleFactorsByMassAndQuantile[str(mass)] = {}
            signalScaleFactorsByMassAndQuantile[str(mass)][str(quantilesExpected[0])] = 1.0
    needInit = True
    for mass in massList:
        datacardDir = dirName.strip("/")+"/datacards"
        cardWorkspace = FindCardWorkspace(datacardDir + "/*", mass)
        if cardWorkspace is None:
            if not doShapeBasedLimits and needInit:
                separateDatacardsDir = dirName+"/datacards"
                if not os.path.isdir(dirName):
                    print("INFO: Making directory", dirName, flush=True)
                    Path(dirName).mkdir(exist_ok=True)
                if not os.path.isdir(separateDatacardsDir):
                    print("INFO: Making directory", separateDatacardsDir, flush=True)
                    Path(separateDatacardsDir).mkdir(exist_ok=True)
                massListFromCards, cardFilesByMass, _ = SeparateDatacards(combinedDatacard, 0, separateDatacardsDir)
                needInit = False
            cardFile = combinedDatacard if doShapeBasedLimits else cardFilesByMass[mass]
            cardWorkspace = ConvertDatacardToWorkspace(cardFile, mass)
    return signalScaleFactorsByMassAndQuantile


# http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#nuisance-parameter-impacts
def MakeImpacts(workspace, mass, dirName, signalScaleFactor, asimovData=False, signal=False):
    impactsDir = dirName.strip("/")+"/impacts"
    if not os.path.isdir(impactsDir):
        print("INFO: Making directory", impactsDir, flush=True)
        Path(impactsDir).mkdir()

    if asimovData:
        if signal:
            print("Using Asimov data with signal + bkg")
            impactsDir+="/sigPlusBkgToys"
        else:
            print("Using Asimov data with bkg only")
            impactsDir+="/bkgOnlyToys"
    else:
        print("Using data")
        impactsDir+="/data"
        
    if not os.path.isdir(impactsDir):
        print("INFO: Making directory", impactsDir, flush=True)
        Path(impactsDir).mkdir()
    with Progress() as progress:
        #signalScaleFactor = signalScaleFactor * 10
        #if signalScaleFactor <=1.0:
        #    signalScaleFactor = 1.0
        task_id = progress.add_task("[cyan]Making impact plots for LQ{}".format(mass), total=5)
        rMin = -10
        if mass==3000:
            rMin = -5
            #signalScaleFactor *= 3
        rMax = 5
        print("INFO: Use rMin = {}".format(rMin))
        stepSize = 0.2
        combToolCmd = "combineTool.py -v3"
        cmd = combToolCmd + " -M Impacts -d {} -m {} --doInitialFit --robustFit 1 --rMin {}".format(workspace, mass, rMin, rMax)
        #cmd += " --cminDefaultMinimizerType GSLMultiMin --cminDefaultMinimizerAlgo SteepestDescent "
        cmd += " --cminDefaultMinimizerStrategy 1 --cminFallbackAlgo Minuit2,Migrad,0:0.1 --cminFallbackAlgo Minuit2,Migrad,1:1.0 --cminFallbackAlgo Minuit2,Migrad,0:1.0 --X-rtd MINIMIZER_MaxCalls=999999999 --X-rtd MINIMIZER_analytic --X-rtd FAST_VERTICAL_MORPH"
        if asimovData:
            cmd += " -t -1"
            if signal:
                cmd += " --expectSignal 1"
        cmd += commonCombineArgs.format(signalScaleFactor)
        print(cmd)
        RunCommand(cmd, impactsDir, None, True)
        if os.path.isfile(impactsDir+"/combine_logger.out"):
            shutil.copy(impactsDir+"/combine_logger.out", impactsDir+"/combine_logger_initialFit{}.out".format(mass))
        progress.update(task_id, advance=1)
        cmd = combToolCmd + " -M Impacts -d {} -m {} --robustFit 1 --doFits --parallel 4 --rMin {}".format(workspace, mass, rMin, rMax)
        #cmd += " --cminDefaultMinimizerType GSLMultiMin --cminDefaultMinimizerAlgo SteepestDescent "
        cmd+=" --rMax 3"
        cmd += " --cminDefaultMinimizerStrategy 1 --cminFallbackAlgo Minuit2,Migrad,0:0.1 --cminFallbackAlgo Minuit2,Migrad,1:1.0 --cminFallbackAlgo Minuit2,Migrad,0:1.0 --X-rtd MINIMIZER_MaxCalls=999999999 --X-rtd MINIMIZER_analytic --X-rtd FAST_VERTICAL_MORPH"
        if asimovData:
            cmd += " -t -1"
            if signal:
                cmd +=" --expectSignal 1"
        cmd += commonCombineArgs.format(signalScaleFactor)
        RunCommand(cmd, impactsDir, None, True)
        if os.path.isfile(impactsDir+"/combine_logger.out"):
            shutil.copy(impactsDir+"/combine_logger.out", impactsDir+"/combine_logger_doFits{}.out".format(mass))
        progress.update(task_id, advance=1)
        cmd = combToolCmd + " -M Impacts -d {} -m {} -o impacts.m{}.json --rMin {}".format(workspace, mass, mass, rMin)
        if asimovData:
            cmd += " -t -1"
            if signal:
                cmd += " --expectSignal 1"
        cmd += commonCombineArgs.format(signalScaleFactor)
        RunCommand(cmd, impactsDir, None, True)
        if os.path.isfile(impactsDir+"/combine_logger.out"):
            shutil.copy(impactsDir+"/combine_logger.out", impactsDir+"/combine_logger_makeJson{}.out".format(mass))
        progress.update(task_id, advance=1)
        cmd = "plotImpacts.py --summary -i impacts.m{}.json -o impacts.m{}".format(mass, mass)
        if not asimovData:
            cmd += " --blind"
        RunCommand(cmd, impactsDir, None, True)
        progress.update(task_id, advance=1)

        cmd = "combine -M MultiDimFit {} --algo grid --saveFitResult --rMin {} --points=1000 --mass {} -v2".format(workspace, rMin, mass)
        cmd += " --cminDefaultMinimizerStrategy 1 --cminFallbackAlgo Minuit2,Migrad,0:0.1 --cminFallbackAlgo Minuit2,Migrad,1:1.0 --cminFallbackAlgo Minuit2,Migrad,0:1.0 --X-rtd MINIMIZER_MaxCalls=999999999 --X-rtd MINIMIZER_analytic --X-rtd FAST_VERTICAL_MORPH"
        if asimovData:
            cmd += " -t -1"
            if signal:
                cmd += " --expectSignal 1"
        cmd += commonCombineArgs.format(signalScaleFactor)
        RunCommand(cmd, impactsDir, None, True)

        cmd = "/afs/cern.ch/user/e/eipearso/public/leptoquark_analysis/CMSSW_14_1_0_pre4/src/HiggsAnalysis/CombinedLimit/scripts/plot1DScan.py higgsCombineTest.MultiDimFit.mH{}.root -o likelihoodScan_LQ{}".format(mass,mass)
        RunCommand(cmd, impactsDir, None, True)
        progress.update(task_id, advance=1)
        if os.path.isfile(impactsDir+"/combine_logger.out"):
            shutil.copy(impactsDir+"/combine_logger.out", impactsDir+"/combine_logger_likelihoodScan{}.out".format(mass))

    
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/HiggsPAGPreapprovalChecks
# but this implementation is not finished
def DoPreapprovalChecks(workspace, mass, dirName, signalScaleFactor):
    preappDir = dirName.strip("/")+"/preappChecks"
    if not os.path.isdir(preappDir):
        print("INFO: Making directory", preappDir, flush=True)
        Path(preappDir).mkdir()
    cmd = "combine -M FitDiagnostics -d {} -m {} -n .m{} -t -1 --expectSignal 0".format(workspace, mass, mass)
    cmd += commonCombineArgs.format(signalScaleFactor)
    RunCommand(cmd, preappDir, None, True)
    # is the result of the fit r=0 ?
    diffNuisances = os.getenv("CMSSW_BASE") + "src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py"
    cmd = "python3 " + diffNuisances + " -a fitDiagnostics.m{}.root -g plots_m{}.root".format(mass, mass)

    
####################################################################################################
# Run
####################################################################################################
if __name__ == "__main__":
    doBatch = True
    gridPointsPerJob = 10
    gridPointsPerJobLowQuantiles = 5
    gridScanPoints = 50
    doShapeBasedLimits = False
    doAsymptoticLimits = False
    blinded = True
    doObservedLimit = False
    useAsimovData = True #do or don't use asimov data in place of real data
    #NOTE: if doObservedLimit==True and useAsimovData==True, combine will calculate the observed limit using asimov toys in place of real data
    ncores = 6
    #massList = list(range(1400, 3100, 100))
    #massList = list(range(300,2100, 100))
    massList = [1500,1600,1700] #[1100,1200,1300,1400,1500,1600,1700,1800,1900,2000]
    betasToScan = list(np.linspace(0.0, 1, 500))[:-1] + [0.9995]
    eosDir = "root://eoscms.cern.ch//eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/eipearso"
    eosDirNoPrefix = eosDir[eosDir.rfind("//")+1:]
    
    quantilesExpected = [0.025, 0.16, 0.5, 0.84, 0.975]
    #quantilesExpected = [0.5]
    if doObservedLimit == True:
        quantilesExpected.append(-2)
    xsThFilename = "$LQANA/config/xsection_theory_13TeV_scalarPairLQ.txt"
    #xsThFilename = "xsection_theory_13TeV_atlas.txt"
    sigRescaleFile = "rValues_nominal.json"
    #sigRescaleFile = "rValues_YM.json"
    #sigRescaleFile = "rValues_minimal.json"
    #xsThFilename = "$LQANA/config/xsection_theory_13TeV_stopPair.txt"
    commonCombineArgs = " --setParameters signalScaleParam={} --freezeParameters signalScaleParam --trackParameters signalScaleParam"
    
    parser = OptionParser(
        usage="%prog -d datacard -n name [commands/options]",
    )
    parser.add_option(
        "-d",
        "--datacard",
        dest="datacard",
        help="combined datacard",
        metavar="datacard",
    )
    parser.add_option(
        "-r",
        "--readResults",
        action="store_true",
        help="read limit results from batch",
        metavar="readResults",
        default=False,
    )
    parser.add_option(
        "-n",
        "--name",
        dest="name",
        help="name of limit calculation (for bookkeeping)",
        metavar="name",
    )
    parser.add_option(
        "-t",
        "--dryRun",
        dest="dryRun",
        help="don't submit jobs to condor, just write the submission files (dry run)",
        metavar="dryRun",
        action="store_true",
        default=False,
    )
    parser.add_option(
        "-i",
        "--doImpacts",
        dest="doImpacts",
        help="produce impact plots",
        metavar="doImpacts",
        action="store_true",
        default=False,
    )
    parser.add_option(
        "-b",
        "--doBetaScan",
        dest="doBetaScan",
        help="do scan of beta values",
        metavar="doBetaScan",
        action="store_true",
        default=False,
    )
    parser.add_option(
        "-e",
        "--estimateRValueScanRange",
        dest="estimateRValueScanRange",
        help="asymptotic limits to estimate r-value scan ranges only",
        metavar="estimateRValueScanRange",
        action="store_true",
        default=False,
    )
    parser.add_option(
        "-l",
        "--submitLimitsAfterGridGen",
        dest="submitLimitJobsAfterGridGen",
        help="submit the limit jobs (after grid generation)",
        metavar="submitLimitJobsAfterGridGen",
        action="store_true",
        default=False,
    )
    parser.add_option(
        "--crab",
        dest="crabMode",
        help="submit jobs with crab3 (crab mode)",
        metavar="crabMode",
        action="store_true",
        default=False,
    )
    parser.add_option(
        "-c",
        "--cmsConnect",
        dest="cmsConnectMode",
        help="submit jobs with CMS COnnect",
        metavar="cmsConnectMode",
        action="store_true",
        default=False,
    )
    parser.add_option(
        "-s",
        "--sandbox",
        dest="cmsswSandbox",
        help="CMSSW sandbox for CMS Connect",
        metavar="cmsswSandbox",
        default=None,
    )
    parser.add_option(
        "-f",
        "--fileWithRValueRanges",
        dest="fileWithRValueRanges",
        help="file containing r-value ranges",
        metavar="fileWithRValueRanges",
        default=None,
    )
    parser.add_option(
        "-x",
        "--xsecFileForSignalRescaling",
        dest="xsecFileForSignalRescaling",
        help="file containing cross sections for rescaling signal samples",
        metavar="xsecFileForSignalRescaling",
        default=None,
    )
    parser.add_option(
        "-g",
        "--generateGridsOnly",
        dest="generateGridsOnly",
        help="Do grid generation only, without computing the limits (need to use -l afterwards to do that in a separate step)",
        metavar="generateGridsOnly",
        action="store_true",
        default=False,
    )
    parser.add_option(
        "-p",
        "--doHiggsPreapprovalChecks",
        dest="doHiggsPreapprovalChecks",
        help="do Higgs PAG preapproval check",
        metavar="doHiggsPreapprovalChecks",
        action="store_true",
        default=False,
    )
    parser.add_option(
        "--signif",
        dest="doSignificance",
        help="do significance calculation",
        metavar="doSignificance",
        action="store_true",
        default=False,
    )
    parser.add_option(
        "--noSigSF",
        dest="disableSigRescaling",
        help="By default, signal is rescaled in the limit calc to set the median expected limit near r=1. Option diables this.",
        metavar="disableSigRescaling",
        action="store_true",
        default=False,
    )
    (options, args) = parser.parse_args()
    if options.datacard is None and options.readResults is None:
        raise RuntimeError("Need either option d to specify datacard, or option r to specify reading limit results from batch")
    if options.name is None:
        raise RuntimeError("Option n to specify name of limit results dir is required")
    if options.doBetaScan and not doBatch:
        raise RuntimeError("Won't do beta scan without batch submission enabled (see doBatch parameter inside script).")
    if options.cmsConnectMode and options.cmsswSandbox is None:
        raise RuntimeError("Can't do CMS Connect submission without specifying a cmssw sandbox using -s. (CMSSW sandbox can be created with the cmssw-sandbox command)")
    #if options.submitLimitJobsAfterGridGen:
    #    gridGenFilesByCondorDir = CheckAllGridGenJobs(massList, quantilesExpected, math.ceil(gridScanPoints/gridPointsPerJob))
    if doBatch and options.crabMode:
        RunCommand("crab createmyproxy")
    manager = multiprocessing.Manager()
    dictAsimovToysByScaleFactor = manager.dict()
    rLimitsByMassAndQuantile = manager.dict()
    xsecLimitsByMassAndQuantile = manager.dict()
    signalScaleFactorsByMassAndQuantile = manager.dict()
    listFailedCommands = manager.list()
    dirName = options.name
    #eosDir += "/"+dirName
    eosDirName = eosDirNoPrefix.rstrip("/") + "/" + dirName
    xsThByMass, _, _ = ReadXSecFile(xsThFilename)
    if not options.disableSigRescaling:
        with open(sigRescaleFile, 'r') as f:
            sigRescaleFactors = json.load(f) #By mass and quantile, but they're the same for all quantiles for a given mass
    rescaleSignal = True
    if options.xsecFileForSignalRescaling is not None:
        xsRescaleByMass, _, _ = ReadXSecFile(options.xsecFileForSignalRescaling)
        if not options.disableSigRescaling:
            for mass in massList:
                xsRescaleByMass[mass] *= sigRescaleFactors[str(mass)]['0.5']
    else:
        if options.disableSigRescaling:
            rescaleSignal = False
        else:
            xsRescaleByMass = {}
            for mass in massList:
                xsRescaleByMass[mass] = xsThByMass[mass] * sigRescaleFactors[str(mass)]['0.5']
    failedBatchCommands = []
    combinedDatacard = options.datacard

    if options.doImpacts:
        print("INFO: Making nuisance parameter impact plots...", flush=True, end="")
        #signalScaleFactorsByMassAndQuantile = InitCardsAndWorkspaces(dirName)
        InitCardsAndWorkspaces(dirName)
        datacardDir = dirName.strip("/")+"/datacards"
        for mass in massList:
            sigSF = sigRescaleFactors[str(mass)]['0.5']
            print("INFO: For mass {} use scale factor {}".format(mass,sigSF))
            cardWorkspace = FindCardWorkspace(datacardDir + "/*", mass)
            MakeImpacts(cardWorkspace, mass, dirName, sigSF)#signalScaleFactorsByMassAndQuantile[str(mass)][str(quantilesExpected[0])]) #Data
            #MakeImpacts(cardWorkspace, mass, dirName, sigSF, True, True)#signalScaleFactorsByMassAndQuantile[str(mass)][str(quantilesExpected[0])], True, True) #Toys, signal + bkg
            #MakeImpacts(cardWorkspace, mass, dirName, sigSF, True, False) #signalScaleFactorsByMassAndQuantile[str(mass)][str(quantilesExpected[0])],True, False) #Toys, bkg only
        print("DONE", flush=True)

    elif options.doHiggsPreapprovalChecks:
        # see: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsWG/HiggsPAGPreapprovalChecks
        print("INFO: Doing steps for Higgs PAG preapproval checks...", flush=True, end="")
        InitCardsAndWorkspaces(dirName)
        datacardDir = dirName.strip("/")+"/datacards"
        for mass in massList:
            cardWorkspace = FindCardWorkspace(datacardDir + "/*", mass)
            DoPreapprovalChecks(cardWorkspace, mass, dirName, signalScaleFactorsByMassAndQuantile[str(mass)][str(quantilesExpected[0])])
        print("DONE", flush=True)

    elif options.doSignificance:
        signalScaleFactor = 1.0
        if not options.readResults:
            separateDatacardsDir = dirName+"/datacards"
            asimovToysDir = dirName + "/asimovDataWithSignal"
            if not os.path.isdir(dirName):
                print("INFO: Making directory",dirName,flush=True)
                Path(dirName).mkdir(exist_ok=True)
            if not os.path.isdir(separateDatacardsDir):
                print("INFO: Making directory",separateDatacardsDir,flush=True)
                Path(separateDatacardsDir).mkdir(exist_ok=True)
            if not os.path.isdir(asimovToysDir):
                print("INFO: Making directory",asimovToysDir,flush=True)
                Path(asimovToysDir).mkdir(exist_ok=True)
            massListFromCards, cardFilesByMass, _ = SeparateDatacards(combinedDatacard, 0, separateDatacardsDir)
            for mass in massList:
                cardFile = cardFilesByMass[mass]
                cardWorkspace = FindCardWorkspace(cardFile, mass)
                if cardWorkspace is not None:
                    print("INFO: Using previously-generated card workspace: {}".format(cardWorkspace), flush=True)
                else:
                    cardWorkspace = ConvertDatacardToWorkspace(cardFile, mass)
                
                asimovToyFile = FindAsimovToyData(mass, signalScaleFactor, asimovToysDir)
                if asimovToyFile is not None:
                    print("INFO: Using previously-generated Asimov toy file: {}".format(asimovToyFile), flush=True)
                else:
                    print("INFO: Generating Asimov toy file", flush=True)
                    asimovToyFile = GenerateAsimovToyData(cardWorkspace, mass, asimovToysDir, dictAsimovToysByScaleFactor, signalScaleFactor, True)
                print("Submitting HybridNew jobs to batch for mass {}...".format(mass))
                SubmitHybridNewBatch((cardWorkspace, mass, dirName, listFailedCommands, 0.5, asimovToyFile, signalScaleFactor)) 

    elif options.readResults and options.estimateRValueScanRange:
        if not options.doBetaScan:
            rValuesByMassAndQuantile = {}
            scanDirName = dirName+"/asymptoticLimits/"
            eosDirNoPrefix = eosDir[eosDir.rfind("//")+1:]
            eosScanDirName = eosDirNoPrefix.rstrip("/") + "/" + dirName
            for mass in massList:
                taskName = 'asymptotic.M{}'
                if rescaleSignal:
                    taskName += ".signalScaleFactor*"
                globString = scanDirName+'/'+taskName.format(mass)+'/condor.*.0.err'
                globString2 = eosScanDirName+'/'+taskName.format(mass)+'/condor.*.0.err'
                errFiles = GetFileList([globString, globString2])
                for f in errFiles:
                    CheckErrorFile(f)
                    globString = scanDirName+"/asymptotic.M{}*/higgsCombineTest.AsymptoticLimits.mH{}.*.root".format(mass, mass)
                    globString2 = eosScanDirName+"/asymptotic.M{}*/higgsCombineTest.AsymptoticLimits.mH{}.*.root".format(mass, mass)
                    rootFile = FindFile([globString, globString2])
                    limits, limitErrs, quantiles, _ = ExtractAsymptoticLimitResultRoot(rootFile)
                    rValuesByMassAndQuantile[mass] = dict(zip(quantiles, limits))
            with open(scanDirName+'/rValues.json', 'w') as f:
                json.dump(rValuesByMassAndQuantile, f)
        else:
            rValuesByMassBetaAndQuantile = {}
            for mass in massList:
                scanDirName = dirName+"/asymptoticLimits/asymptotic.M{}".format(mass)
                betasToSubmit = GetBetasToSubmit(mass)
                rValuesByMassBetaAndQuantile[mass] = {}
                globString = scanDirName+'/signalScaleFactor*.betaIndex*/condor*.0.err'.format(mass)
                print(globString)
                errFiles = GetFileList(globString)
                fileList = []
                for f in errFiles:
                    CheckErrorFile(f, False)
                    betaValIndex = int(f[f.rfind("betaIndex")+9:f.find("/", f.rfind("betaIndex")+9)])
                    if not options.cmsConnectMode:
                        #glob does not seem to be able to find files in the exo eos space when running on cms connect
                        #you can list/copy them with xrd from the exo space to cms connect, but you have to know the exact filename
                        #Since the asymptotic jobs have the seed in the filename, and xrd doesn't support wildcards
                        #(as far as I know), when we run from cms connect we have to get the limit result from the condor.out file
                        #-Emma, 15 jan 2025
                        fileList.append(str(FindFile(eosDir+"/"+dirName+"/asymptotic.M{}/higgsCombine.betaIndex{}.AsymptoticLimits.mH{}.*.root".format(mass, betaValIndex, mass))))
                    else:
                        fileList.append(str(FindFile(f.replace(".err",".out"))))
                for f in fileList:
                    betaValIndex = int(f[f.rfind("betaIndex")+9:f.find("/", f.rfind("betaIndex")+9)])
                    if ".root" in f:
                        limits, limitErrs, quantiles, signalScaleParam = ExtractAsymptoticLimitResultRoot(f)
                    elif ".out" in f or ".txt" in f:
                        limits, limitErrs, quantiles, signalScaleParam = ExtractAsymptoticLimitResultTxt(f)
                    else:
                        print("Error: Don't know how to extract limit from file {}".format(f))
                    print(betaValIndex)
                    rValuesByMassBetaAndQuantile[mass][betasToScan[betaValIndex]] = dict(zip(quantiles, limits))
            if not os.path.isdir(dirName+"/asymptoticLimits"):
                os.mkdir(dirName+"/asymptoticLimits")
            with open(dirName+"/asymptoticLimits"+'/rValues.json', 'w') as f:
                json.dump(rValuesByMassBetaAndQuantile, f)

    elif not doBatch or options.readResults and not options.estimateRValueScanRange:
        years, intLumi = GetYearAndIntLumiFromDatacard(combinedDatacard)
        intLumi = intLumi/1000.0
        if intLumi < 100:
            intLumi = round(intLumi, 1)
        else:
            intLumi = round(intLumi, 0)
        print("Found total int lumi = {}/fb for years = {}".format(intLumi, years))
        if not options.doBetaScan:
            if options.readResults:
                if os.path.isfile(dirName+'/xsecLimits.json') and os.path.isfile(dirName+'/signalScaleFactors.json'):
                    print("INFO: read xsection limits from "+dirName+'/xsecLimits.json')
                    with open(dirName+'/xsecLimits.json','r') as f:
                        xsecLimitsByMassAndQuantile = json.load(f)
                    with open(dirName+'/signalScaleFactors.json','r') as f:
                        signalScaleFactorsByMassAndQuantile = json.load(f)
                else:
                    condorDir = dirName.rstrip("/") # +"/condor"
                    print("INFO: Reading results from batch from {}...".format(condorDir), flush=True)
                    xsecLimitsByMassAndQuantile, signalScaleFactorsByMassAndQuantile = ReadBatchResults(massList, condorDir)
                    if len(xsecLimitsByMassAndQuantile) and len(signalScaleFactorsByMassAndQuantile):
                        print("DONE", flush=True)
                    else:
                        raise RuntimeError("Got empty dicts back from ReadBatchResults; something wrong happened there.")
                    with open(dirName+'/xsecLimits.json', 'w') as f:
                        json.dump(get_value(xsecLimitsByMassAndQuantile), f)
                    with open(dirName+'/signalScaleFactors.json', 'w') as f:
                        json.dump(get_value(signalScaleFactorsByMassAndQuantile), f)
            masses, shadeMasses, xsMedExp, xsOneSigmaExp, xsTwoSigmaExp, xsObs = CreateArraysForPlotting(xsecLimitsByMassAndQuantile)
            print("mData =", list(masses))
            # print("x_shademasses =", list(shadeMasses))
            # print("xsUp_expected =", list(xsMedExp))
            # print("xsUp_observed =", list(xsObs))
            # print("y_1sigma =", list(xsOneSigmaExp))
            # print("y_2sigma =", list(xsTwoSigmaExp))
            table, columnNames = MakeResultTable(masses, xsMedExp, xsObs, xsOneSigmaExp, xsTwoSigmaExp)
            print(tabulate(table, headers=columnNames, tablefmt="fancy_grid", floatfmt=".6f"))
            stopTime = time.time()
            if not doBatch:
                execTimeStr = GetExecTimeStr(startTime, stopTime)
                print("Total limit calculation execution time:", execTimeStr)
            print("Make plot and calculate mass limit")
            if options.xsecFileForSignalRescaling is not None:
                thXsecToPlot = options.xsecFileForSignalRescaling
            else:
                thXsecToPlot = xsThFilename
            BR_Sigma_EE_vsMass(dirName, intLumi, masses, shadeMasses, xsMedExp, xsObs, xsOneSigmaExp, xsTwoSigmaExp, thXsecToPlot)
        else:
            mainDirName = dirName
            betaDirName = dirName #+"/betaScan"
            condorDir = betaDirName.strip("/")#+"/condor"
            if os.path.isfile(betaDirName+'/xsecLimits.json'):
                print("INFO: Reading beta scan results from json files in {}...".format(betaDirName), flush=True)
                with open(betaDirName+'/xsecLimits.json', 'r') as f:
                    xsecLimitsByMassAndQuantile = json.load(f)
                with open(betaDirName+'/rLimits.json', 'r') as f:
                    rLimitsByMassAndQuantile = json.load(f)
                with open(betaDirName+'/massLimits.json', 'r') as f:
                    massLimitsByQuantileAndBetaVal = json.load(f)
                with open(betaDirName+'/comboArrays.json', 'r') as f:
                    comboArraysByQuantile = json.load(f)
            else:
                print("INFO: Reading beta scan results from batch from {}...".format(condorDir), flush=True)
                ReadBatchResultsBetaScan(massList, condorDir)
                # print("DONE", flush=True)
                # print("xsecLimitsByMassAndQuantile={}".format(dict(xsecLimitsByMassAndQuantile)))
                # print("xsecLimitsByMassAndQuantile at quantile 0.5={}".format(get_value(xsecLimitsByMassAndQuantile["0.5"])))
                with open(betaDirName+'/xsecLimits.json', 'w') as f:
                    json.dump(get_value(xsecLimitsByMassAndQuantile), f)
                with open(betaDirName+'/rLimits.json', 'w') as f:
                    json.dump(get_value(rLimitsByMassAndQuantile), f)
                print("INFO: Compute beta limits...", flush=True)
                massLimitsByQuantileAndBetaVal = ComputeBetaLimits(xsThByMass, xsecLimitsByMassAndQuantile)
                with open(betaDirName+'/massLimits.json', 'w') as f:
                    json.dump(get_value(massLimitsByQuantileAndBetaVal), f)
                # print("massLimitsByQuantileAndBetaVal={}".format(massLimitsByQuantileAndBetaVal))
                # print("DONE", flush=True)
                comboArraysByQuantile = CreateComboArrays(massLimitsByQuantileAndBetaVal)
                with open(betaDirName+'/comboArrays.json', 'w') as f:
                    json.dump(comboArraysByQuantile, f)
            # print("comboArraysByQuantile['0.5']['betas']={}".format(comboArraysByQuantile["0.5"]["betas"]))
            # print("comboArraysByQuantile['0.5']['massLimits']={}".format(comboArraysByQuantile["0.5"]["massLimits"]))
            ComboPlot(mainDirName, intLumi, comboArraysByQuantile) # , m_observed_lljj)
            for quantile in quantilesExpected:
                table, columnNames = MakeResultTableBetaScan(massList, rLimitsByMassAndQuantile, str(quantile))
                print("rLimits beta scan for quantile={}:".format(quantile))
                print(tabulate(table, headers=columnNames, tablefmt="fancy_grid", floatfmt=".6f"))
            # table, columnNames = MakeResultTableBetaScan(massList, xsecLimitsByMassAndQuantile)
            # print("xsecLimits beta scan for quantile=0.5:")
            # print(tabulate(table, headers=columnNames, tablefmt="fancy_grid", floatfmt=".6f"))
            headers = ["MLQ / beta="]
            headers.extend(comboArraysByQuantile["0.5"]["betas"])
            tableRow = ["mass"]
            tableRow.extend(comboArraysByQuantile["0.5"]["massLimits"])
            table = [tableRow]
            print("mass limits beta scan for quantile=0.5:")
            print(tabulate(table, headers=headers, tablefmt="fancy_grid", floatfmt=".6f"))

    else:
        # submit limit jobs
        startTime = time.time()
        #if options.submitLimitJobsAfterGridGen:
         #   CheckAllGridGenJobs(massList, quantilesExpected, math.ceil(gridScanPoints/gridPointsPerJob))
        #TODO: check for previous datacards?
        separateDatacardsDir = dirName+"/datacards"
        asimovToysDir = dirName+"/asimovData"
        if not os.path.isdir(dirName):
            print("INFO: Making directory", dirName, flush=True)
            Path(dirName).mkdir(exist_ok=True)
        if not os.path.isdir(separateDatacardsDir):
            print("INFO: Making directory", separateDatacardsDir, flush=True)
            Path(separateDatacardsDir).mkdir(exist_ok=True)
        if not os.path.isdir(asimovToysDir) and useAsimovData:
            print("INFO: Making directory", asimovToysDir, flush=True)
            Path(asimovToysDir).mkdir(exist_ok=True)
        
        if not doShapeBasedLimits:
            massListFromCards, cardFilesByMass, _ = SeparateDatacards(combinedDatacard, 0, separateDatacardsDir)
        for mass in massList:
            cardFile = combinedDatacard if doShapeBasedLimits else cardFilesByMass[mass]
            if not options.estimateRValueScanRange and not options.doImpacts:
                print("INFO: Computing limits for mass {}".format(mass), flush=True)
            else:
                print("INFO: Estimate r-value scan range for mass {}".format(mass), flush=True)
            rLimitsByMassAndQuantile[mass] = {}
            xsecLimitsByMassAndQuantile[mass] = {}
            signalScaleFactorsByMassAndQuantile[mass] = {}
            cardWorkspace = FindCardWorkspace(cardFile, mass)
            if cardWorkspace is not None:
                print("INFO: Using previously-generated card workspace: {}".format(cardWorkspace), flush=True)
            else:
                cardWorkspace = ConvertDatacardToWorkspace(cardFile, mass)
            if not options.estimateRValueScanRange and doAsymptoticLimits:
                print("INFO: Doing interactive AsymptoticLimits for mass {}".format(mass), flush=True)
                asStartTime = time.time()
                rootFileName = RunAsymptotic(cardWorkspace, mass, dirName, blinded, 1.0)
                #TODO: should we also run again, scaling the signal yields, as per below? it might not matter.
                asStopTime = time.time()
                execTimeStr = GetExecTimeStr(asStartTime, asStopTime)
                print("Asymptotic calculation execution time:", execTimeStr, flush=True)
                limits, limitErrs, quantiles = ExtractAsymptoticLimitResultRoot(rootFileName)
                for index, quantileExp in enumerate(quantiles):
                    limit = limits[index]
                    rLimitsByMassAndQuantile[mass][str(quantileExp)] = limit
                    xsecLimitsByMassAndQuantile[mass][str(quantileExp)] = limit * xsThByMass[float(mass)]
            else:
                if options.estimateRValueScanRange:
                    if not options.doBetaScan:
                        rValuesByMassAndQuantile = {}
                        if doBatch:
                            # scanDirName = dirName+"/asymptoticLimits"
                            # if not os.path.isdir(scanDirName):
                            #     print("INFO: Making directory", scanDirName, flush=True)
                            #     Path(scanDirName).mkdir(exist_ok=True)
                            # MakeCondorDirs(scanDirName)
                            signalScaleFactor = 1.0
                            if rescaleSignal:
                                signalScaleFactor = xsRescaleByMass[mass]/xsThByMass[mass]  # we assume that the signal yields as given have been normalized with the scalar pair-prod xsec
                            SubmitAsymptoticBatch((cardWorkspace, mass, dirName, listFailedCommands, blinded, signalScaleFactor, -1))
                            continue
                        else:
                            print("INFO: Doing interactive AsymptoticLimits for mass {} to estimate scan range".format(mass), flush=True)
                            asymptoticRootFileName = RunAsymptoticInteractive(cardWorkspace, mass, dirName, blinded, 1.0)
                            limits, limitErrs, quantiles = ExtractAsymptoticLimitResultRoot(asymptoticRootFileName)
                            rValuesByQuantile = dict(zip(quantiles, limits))
                            #FIXME: add writing of r-value scan file here
                            rValuesByMassAndQuantile[mass] = dict(zip(quantiles, limits))
                    else:
                        rValuesByMassBetaAndQuantile = {}
                        betaDirName = dirName #+"/betaScan/"
                        #if not os.path.isdir(betaDirName):
                        #    print("INFO: Making directory", betaDirName, flush=True)
                        #    Path(betaDirName).mkdir(exist_ok=True, parents=True)
                        #MakeCondorDirs(betaDirName)
                        betasToSubmit = GetBetasToSubmit(mass)
                        rValuesAtBeta = {}
                        if doBatch:
                            with Progress() as progress:
                                asymptJobs = len(betasToSubmit)
                                task_id = progress.add_task("[cyan]Submitting asymptotic limit jobs to batch...", total=asymptJobs)
                                with multiprocessing.Pool(ncores) as pool:
                                    partial_quit = partial(ErrorCallback, pool)
                                    for beta in betasToSubmit:
                                        idx = betasToScan.index(beta)
                                        signalScaleFactor = beta*beta
                                        try:
                                            pool.apply_async(SubmitAsymptoticBatch, [[cardWorkspace, mass, betaDirName, listFailedCommands, blinded, signalScaleFactor, idx]], callback = lambda x: progress.advance(task_id), error_callback = partial_quit)
                                            # asimovJobs += 1
                                        except KeyboardInterrupt:
                                            print("\n\nCtrl-C detected: Bailing.")
                                            pool.terminate()
                                            sys.exit(1)
                                        except Exception as e:
                                            print("ERROR: caught exception in asimov job submission for scale factor: {}".format(signalScaleFactor))
                                            # traceback.print_exc()
                                            # pool.terminate()
                                            # exit(-2)
                                    # now close the pool and wait for jobs to finish
                                    pool.close()
                                    # if asimovJobs > 0:
                                    #     print("Waiting for {} AsymptoticLimits job submission to finish for all beta values for mass={}...".format(asimovJobs, mass))
                                    pool.join()
                                continue
                        else:
                            print("INFO: Doing interactive AsymptoticLimits for mass {}, beta {} to estimate scan range".format(mass, beta), flush=True)
                            asymptoticRootFileName = RunAsymptoticInteractive(cardWorkspace, mass, betaDirName, blinded, signalScaleFactor)
                            limits, limitErrs, quantiles = ExtractAsymptoticLimitResultRoot(asymptoticRootFileName)
                            rValuesByQuantile = dict(zip(quantiles, limits))
                            rValuesAtBeta[beta] = rValuesByQuantile
                            rValuesByMassBetaAndQuantile[mass] = rValuesAtBeta
                elif options.fileWithRValueRanges is not None:
                    with open(options.fileWithRValueRanges) as f:
                        if options.doBetaScan:
                            rValuesByMassBetaAndQuantile = json.load(f)
                        else:
                            rValuesByMassAndQuantile = json.load(f)
                else:
                    raise RuntimeError("Either the r-value scan range must be estimated first with option '-e', or there must be a file of scan ranges given with '-f' in order to submit limit jobs")

                # now submit jobs for limits
                signalScaleFactor = 1.0
                if rescaleSignal:
                    signalScaleFactor = xsRescaleByMass[mass]/xsThByMass[mass]  # we assume that the signal yields as given have been normalized with the scalar pair-prod xsec
                if not options.doBetaScan:
                    printedAsimovToyInfo = False
                    with Progress() as progress:
                            if useAsimovData:
                                asimovToyFile = FindAsimovToyData(mass, signalScaleFactor, asimovToysDir)
                                if asimovToyFile is not None:
                                    if not printedAsimovToyInfo:
                                        print("INFO: Using previously-generated Asimov toy file: {}".format(asimovToyFile), flush=True)
                                        printedAsimovToyInfo = True
                                else:
                                    print("INFO: Generating Asimov toy file", flush=True)
                                    asimovToyFile = GenerateAsimovToyData(cardWorkspace, mass, asimovToysDir, dictAsimovToysByScaleFactor, signalScaleFactor)
                                    printedAsimovToyInfo = True
                            else:
                                print("INFO: not using Asimov toys")
                                asimovToyFile = ""
                            if doBatch:
                                # MakeCondorDirs(dirName)
                                # rValuesByQuantile = rValuesByMassAndQuantile[str(mass)]
                                # failedCmd = SubmitHybridNewBatch((cardWorkspace, mass, dirName, listFailedCommands, quantileExp, asimovToyFile, signalScaleFactor))
                                # if failedCmd is not None:
                                #     failedBatchCommands.append(failedCmd)
                                jobs = len(quantilesExpected)
                                task_id = progress.add_task("[cyan]Submitting HybridNewjobs to batch...", total=jobs)
                                with multiprocessing.Pool(ncores) as pool:
                                    partial_quit = partial(ErrorCallback, pool)
                                    for quantileExp in quantilesExpected:
                                        try:
                                            pool.apply_async(SubmitHybridNewBatch, [[cardWorkspace, mass, dirName, listFailedCommands, quantileExp, asimovToyFile, signalScaleFactor]], callback = lambda x: progress.advance(task_id), error_callback = partial_quit)
                                        except KeyboardInterrupt:
                                            print("\n\nCtrl-C detected: Bailing.")
                                            pool.terminate()
                                            sys.exit(1)
                                        except Exception as e:
                                            print("ERROR: caught exception in asimov job submission for scale factor: {}".format(signalScaleFactor))
                                    # now close the pool and wait for jobs to finish
                                    pool.close()
                                    pool.join()
                                continue
                            else:
                                for quantileExp in quantilesExpected:
                                    hnStartTime = time.time()
                                    rootFileName = RunHybridNewInteractive(cardWorkspace, mass, dirName, quantile=quantileExp, genAsimovToyFile=asimovToyFile, signalScaleFactor=signalScaleFactor)
                                    hnStopTime = time.time()
                                    execTimeStr = GetExecTimeStr(hnStartTime, hnStopTime)
                                    print("HybridNew calculation execution time:", execTimeStr, flush=True)
                                    limit, limitErr, quantile, signalScaleFactor = ExtractLimitResult(rootFileName)
                                    rLimitsByMassAndQuantile[mass][str(quantileExp)] = limit
                                    xsecLimitsByMassAndQuantile[mass][str(quantileExp)] = limit * xsThByMass[float(mass)]
                                    signalScaleFactorsByMassAndQuantile[mass][str(quantileExp)] = signalScaleFactor

                else:
                    # beta scan jobs
                    betaDirName = dirName # +"/betaScan"
                    asimovToysBetaDir = betaDirName+"/asimovToys"
                    if not os.path.isdir(asimovToysBetaDir):
                        print("INFO: beta scan - Making directory", asimovToysBetaDir, flush=True)
                        Path(asimovToysBetaDir).mkdir(parents=True)
                    #MakeCondorDirs(betaDirName)
                    rValuesAtBeta = rValuesByMassBetaAndQuantile[str(mass)]
                    with Progress() as progress:
                        rValuesAtBeta = rValuesByMassBetaAndQuantile[str(mass)]
                        asimovScaleFactorsToGen = []
                        for idx, beta in enumerate(rValuesAtBeta.keys()):
                            signalScaleFactor = float(beta)*float(beta)
                            asimovToyFile = FindAsimovToyData(mass, signalScaleFactor, asimovToysBetaDir)
                            #print(asimovToyFile)
                            if asimovToyFile is None:
                                asimovScaleFactorsToGen.append(signalScaleFactor)
                        if len(asimovScaleFactorsToGen):
                            task_id = progress.add_task("[cyan]Running asimov toy jobs for mass {}".format(mass), total=len(asimovScaleFactorsToGen))
                            with multiprocessing.Pool(ncores) as pool:
                                partial_quit = partial(ErrorCallback, pool)
                                for signalScaleFactor in asimovScaleFactorsToGen:
                                    try:
                                        #print("generate asimov toy for {}".format(signalScaleFactor))
                                        pool.apply_async(GenerateAsimovToyData, [cardWorkspace, mass, asimovToysBetaDir, dictAsimovToysByScaleFactor, signalScaleFactor], callback = lambda x: progress.advance(task_id), error_callback = partial_quit)
                                        # asimovJobs += 1
                                    except KeyboardInterrupt:
                                        print("\n\nCtrl-C detected: Bailing.")
                                        pool.terminate()
                                        sys.exit(1)
                                    except Exception as e:
                                        print("ERROR: caught exception in asimov toy job for scale factor: {}; exiting".format(signalScaleFactor))
                                        traceback.print_exc()
                                        pool.terminate()
                                        exit(-2)
                                # now close the pool and wait for jobs to finish
                                pool.close()
                                pool.join()
                    # HybridNew
                    with Progress() as progress:
                        rValuesAtBeta = rValuesByMassBetaAndQuantile[str(mass)]
                        hybridNewJobs = len(rValuesAtBeta.keys())
                        task_id = progress.add_task("[cyan]Submitting HybridNew jobs to batch for mass {}".format(mass), total=hybridNewJobs)
                        with multiprocessing.Pool(ncores,maxtasksperchild=1) as pool:
                            partial_quit = partial(ErrorCallback, pool)
                            for qExp in quantilesExpected:
                                #if options.cmsConnectMode and options.submitLimitJobsAfterGridGen:
                                #    inputFileList = GetFilesInDirXrd(eosDir+"/"+options.name+"/hybridNewGridGen.M{}.{}".format(mass,str(qExp).replace(".","p")))
                                #else:
                                inputFileList = []
                                for idx, beta in enumerate(rValuesAtBeta.keys()):
                                    signalScaleFactor = float(beta)*float(beta)
                                    try:
                                        asimovToyFile = FindAsimovToyData(mass, signalScaleFactor, asimovToysBetaDir)
                                        pool.apply_async(SubmitHybridNewBatch, [[cardWorkspace, mass, betaDirName, listFailedCommands, qExp, asimovToyFile, signalScaleFactor,inputFileList]], callback = lambda x: progress.advance(task_id), error_callback = partial_quit)
                                    # hybridNewJobs += 1
                                    except KeyboardInterrupt:
                                        print("\n\nCtrl-C detected: Bailing.")
                                        pool.terminate()
                                        sys.exit(1)
                                    except Exception as e:
                                        print("ERROR: caught exception in hybridnew job for scale factor: {}; exiting".format(signalScaleFactor))
                                        traceback.print_exc()
                                        pool.terminate()
                                        exit(-2)
                            # now close the pool
                            pool.close()
                            # print("Waiting for submission of {} HybridNew jobs for mass={}, quantile={}...".format(hybridNewJobs, mass, quantileExp))
                            pool.join()
                            # if failedCmd is not None:
                            #     failedBatchCommands.append(failedCmd)
        
    if len(failedBatchCommands):
        print(colored("batch commands failed:", "red"))
        for item in failedBatchCommands:
            print(item)
    if len(listFailedCommands):
        print(colored("batch commands failed:", "red"))
        for item in listFailedCommands:
            print(item)

