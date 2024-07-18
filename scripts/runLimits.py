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
from pathlib import Path
from bisect import bisect
from tabulate import tabulate
from ROOT import TFile, TGraph, TSpline3

if not 'LQANA' in os.environ:
    raise RuntimeError('Need to define $LQANA environment variable to point to your git clone of the rootNtupleAnalyzerV2 repo.')
sys.path.append(os.getenv("LQANA").rstrip("/") + "/scripts/")

from combineCommon import SeparateDatacards, GetYearAndIntLumiFromDatacard
from BR_Sigma_EE_vsMass import BR_Sigma_EE_vsMass
from ComboPlotLQ1 import ComboPlot

if not 'CMSSW_VERSION' in os.environ:
    raise RuntimeError('You must be in a CMSSW environment (release area + cmsenv)')

ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')

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
    print(colored("\t{}".format(cmd), "green"), flush=True)
    try:
        useShell = False
        # if env is None:
        #     useShell = True
        #     cmd = shlex.split(cmd)
        cmd = shlex.split(cmd)
        process = subprocess.run(cmd, shell=useShell, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=workingDir, env=env)
    except subprocess.CalledProcessError as e:
        result = ansi_escape.sub('', e.output.decode())
        print(colored("RunCommand had an error; output: {}".format(result), "red"), flush=True)
        raise e
    if not suppressOutput:
        print(process.stdout.decode(), flush=True)
        #print(process.stderr.decode())


def ConvertDatacardToWorkspace(datacard, mass):
    workspaceFileName = datacard+".m{}.root".format(mass)
    cmd = 'text2workspace.py -m {} {} -o {} '.format(mass, datacard, workspaceFileName)
    RunCommand(cmd)
    return Path(workspaceFileName).resolve()


def FindCardWorkspace(cardFile, mass, exceptOnAbsence=False):
    workspaceFileName = cardFile+".m{}.root".format(mass)
    listOfWorkspaceFiles = sorted(glob.glob(workspaceFileName), key=os.path.getmtime)
    if exceptOnAbsence and len(listOfWorkspaceFiles) != 1:
        raise RuntimeError("Globbing for {} did not result in one file as expected, but {}".format(workspaceFileName, listOfWorkspaceFiles))
    if len(listOfWorkspaceFiles) > 0:
        return Path(listOfWorkspaceFiles[-1]).resolve()


def GenerateAsimovToyData(args):
    workspace, mass, toysDir, dictAsimovToysByScaleFactor, signalScaleFactor = args
    cmd = 'combine -M GenerateOnly {} -t -1 --seed -1 --saveToys -n .asimov.{} -m {} {}'.format(workspace, signalScaleFactor, mass, commonCombineArgs.format(signalScaleFactor))
    RunCommand(cmd, toysDir)
    toyFile = Path(sorted(glob.glob(toysDir+'/higgsCombine.asimov.{}.GenerateOnly.mH{}.*.root'.format(signalScaleFactor, mass)), key=os.path.getmtime)[-1]).resolve()
    dictAsimovToysByScaleFactor[signalScaleFactor] = toyFile
    return toyFile


def FindAsimovToyData(mass, signalScaleFactor, toysDir):
    listOfToyFiles = sorted(glob.glob(toysDir+'/higgsCombine.asimov.{}.GenerateOnly.mH{}.*.root'.format(signalScaleFactor, mass)), key=os.path.getmtime)
    if len(listOfToyFiles) > 0:
        return Path(listOfToyFiles[-1]).resolve()


def FindFile(globString):
    fileList = sorted(glob.glob(globString), key=os.path.getmtime)
    if len(fileList) != 1:
        raise RuntimeError("Globbing for {} did not result in one file as expected, but {}".format(globString, fileList))
    return Path(fileList[-1]).resolve()


def GetFileList(globString):
    fileList = sorted(glob.glob(globString), key=os.path.getmtime)
    if len(fileList) < 1:
        raise RuntimeError("Globbing for {} did not result in any files.".format(globString))
    return fileList


def MakeCondorDirs(dirName):
    condorDir = dirName.strip("/")+"/condor"
    condorSubfileDir = dirName.strip("/")+"/condor/subFiles"
    condorOutDir = dirName.strip("/")+"/condor/out"
    condorErrDir = dirName.strip("/")+"/condor/error"
    condorLogDir = dirName.strip("/")+"/condor/log"
    dirsToHave = [condorDir, condorSubfileDir, condorOutDir, condorErrDir, condorLogDir]
    for idir in dirsToHave:
        if not os.path.isdir(idir):
            print("INFO: Making directory", idir, flush=True)
            Path(idir).mkdir(exist_ok=True)


def WriteCondorSubFile(filename, shFilename, rMin, rMax, queue="workday", extractLimits=False, inputFiles=[]):
    scanJobs = 1
    if rMin > 0 and rMax > 0 and not extractLimits:
        scanJobs = gridScanPoints / gridPointsPerJob
    basename = shFilename.replace("condor_", "").replace(".sh", "")
    with open(filename, "w") as subfile:
        subfile.write("executable = subFiles/" + shFilename + "\n")
        subfile.write("arguments = $(ProcId)\n")
        subfile.write("output                = out/"   + basename + ".$(ClusterId).$(ProcId).out\n")
        subfile.write("error                 = error/" + basename + ".$(ClusterId).$(ProcId).err\n")
        subfile.write("log                   = log/"   + basename + ".$(ClusterId).$(ProcId).log\n")
        subfile.write("\n")
        if len(inputFiles):
            subfile.write("transfer_input_files = {}\n".format(",".join(inputFiles)))
            subfile.write("\n")
        subfile.write('transfer_output_remaps = "combine_logger.out=out/' + basename + '.$(ClusterId).$(ProcId).combine_logger.out"\n')
        subfile.write("\n")
        subfile.write("should_transfer_files = YES\n")
        subfile.write("when_to_transfer_output = ON_EXIT\n")
        subfile.write("\n")
        subfile.write("# Send the job to Held state on failure.\n")
        subfile.write("on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n")
        subfile.write("\n")
        subfile.write("# Periodically retry the jobs every 10 minutes, up to a maximum of 5 retries.\n")
        subfile.write("periodic_release =  (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)\n")
        subfile.write("\n")
        subfile.write("+JobFlavour=\"{}\"\n".format(queue))
        subfile.write("MY.WantOS = \"el9\"\n")
        subfile.write("queue {}\n".format(int(scanJobs)))


def WriteCondorShFile(condorDir, filename, mass, combineCmd, limitCmd="", signalScaleFactor=1.0, quantile=-1, rMin=-1, rMax=-1, betaIndex=-1):
    with open(filename, "w") as shfile:
        shfile.write("#!/bin/sh\n")
        shfile.write("ulimit -s unlimited\n")
        shfile.write("set -e\n")
        shfile.write("cd {}\n".format(condorDir))
        shfile.write("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
        shfile.write("eval `scramv1 runtime -sh`\n")
        shfile.write("cd -\n")
        shfile.write("\n")
        rootFileString = ""
        sigSFRound = round(signalScaleFactor, 6)
        if betaIndex != -1:
            rootFileString = "betaIndex"+betaIndex
        if len(limitCmd):
            if len(combineCmd):
                raise RuntimeError("Can't handle case when we have both a combine command '{}' and a limit command '{}' specified.". format(combineCmd, limitCmd))
            quant = "quant{}.".format(quantile) if quantile > 0 else ""
            combinedOutputFile = "higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}gridAll.root".format(sigSFRound, mass, quant)
            shfile.write("hadd -fk207 {} higgsCombine.*.root\n".format(combinedOutputFile))
            shfile.write(limitCmd.format(combinedOutputFile)+"\n")
        elif rMin > 0 and rMax > 0:
            stepSize = (rMax-rMin)/gridScanPoints
            pointsPerJob = gridPointsPerJob
            if gridScanPoints < pointsPerJob:
                pointsPerJob = gridScanPoints
            print("INFO: Writing commands to compute grid of {} limits in the range r=[{}, {}] by steps of {}".format(gridScanPoints, rMin, rMax, stepSize))
            shfile.write("if [ $1 -eq 0 ]; then\n")
            jobIdx = 0
            scanPoints = np.linspace(rMin, rMax, num=gridScanPoints)
            for rValIdx, rVal in enumerate(scanPoints):
                thisStepCmd = combineCmd
                thisStepCmd += ' -n .signalScaleFactor{}.POINT.{}'.format(sigSFRound, rVal)
                if len(rootFileString):
                    thisStepCmd += "." + rootFileString
                thisStepCmd += ' --singlePoint {}'.format(rVal)
                shfile.write("  {}\n".format(thisStepCmd))
                if (rValIdx+1) % pointsPerJob == 0:
                    quant = "quant{}.".format(quantile) if quantile > 0 else ""
                    combinedOutputFile = "higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}grid{}.root".format(sigSFRound, mass, quant, jobIdx)
                    shfile.write("  hadd -fk207 {} higgsCombine.*.root\n".format(combinedOutputFile))
                    shfile.write("  ls -ltrph\n")
                    shfile.write("  rm higgsCombine.*POINT.*root\n")  # remove individual point files so that condor doesn't transfer them back
                    shfile.write("fi\n")
                    jobIdx += 1
                    if rValIdx+1 < len(scanPoints):
                        shfile.write("if [ $1 -eq {} ]; then\n".format(jobIdx))
        else:
            if len(rootFileString):
                combineCmd += ' -n {}'.format(rootFileString)
            # print("WriteCondorShFile for mass {}: write cmd='{}'".format(mass, combineCmd))
            shfile.write("  {}\n".format(combineCmd))
            combinedOutputFile = ""
    return combinedOutputFile


def GetHybridNewCommandArgs(workspace, mass, dirName, quantile, genAsimovToyFile, signalScaleFactor, batch, rMin, rMax):
    rAbsAcc = 0.0001
    clsAcc = 0 if batch else 0.001
    toys = 10000
    # toys = 500  # reduced for shape-based limits as test, but still took hours
    cmd = '{}'.format(str(workspace).split("/")[-1] if doBatch else workspace)
    cmd += ' -v1'
    cmd += ' -M HybridNew'
    cmd += ' --LHCmode LHC-limits'
    cmd += ' --saveHybridResult'
    cmd += ' --saveToys'
    cmd += ' --seed -1'
    if not batch:
        cmd += ' --iterations 2'
        cmd += ' --rAbsAcc {}'.format(rAbsAcc)
    cmd += ' --clsAcc {}'.format(clsAcc)
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
        if genAsimovToyFile != "":
            cmd += ' -D {}:toys/toy_asimov'.format(str(genAsimovToyFile).split("/")[-1] if doBatch else genAsimovToyFile)
    return cmd


def GetAsymptoticCommandArgs(workspace, mass, dirName, blinded, signalScaleFactor, rMin, rMax):
    rAbsAcc = 0.0001
    # rMax = 150
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


def SubmitLimitJobsBatch(mass, dirName, listFailedCommands, quantile, signalScaleFactor, rMin, rMax, taskName, cmdArgs, limitCmd, queue="workday", inputFiles=[]):
    # print("SubmitLimitJobsBatch() called with cmdArgs='{}'".format(cmdArgs))
    condorDir = dirName.strip("/")+"/condor"
    condorSubfileDir = dirName.strip("/")+"/condor/subFiles"
    cmd = ""
    if len(cmdArgs):
        cmd = "combine " + cmdArgs
    shFilename = condorSubfileDir + "/condor_" + taskName + ".sh"
    condorFile = shFilename.replace(".sh", ".sub")
    if "betaIndex" in taskName:
        betaIndex = taskName[taskName.find("betaIndex")+9:]
    else:
        betaIndex = -1
    WriteCondorShFile(Path(condorDir).resolve(), shFilename, mass, cmd, limitCmd, signalScaleFactor, quantile, rMin, rMax, betaIndex)
    WriteCondorSubFile(condorFile, Path(shFilename).name, rMin, rMax, queue, len(limitCmd), inputFiles)
    if not options.dryRun:
        cmd = "condor_submit subFiles/{}".format(Path(condorFile).name)
        runCommandArgs = [cmd, condorDir, None, True]
        try:
            RunCommand(*runCommandArgs)
        except subprocess.CalledProcessError as e:
            try:
                print(colored("Caught exception running condor_submit; retry", "yellow"), flush=True)
                RunCommand(*runCommandArgs)
            except subprocess.CalledProcessError as e:
                listFailedCommands.append( runCommandArgs[0:1])
                return runCommandArgs[0:1]
    else:
        print("INFO: dry run enabled; not submitting subFiles/{} to condor".format(Path(condorFile).name))


def SubmitHybridNewBatch(args):
    workspace, mass, dirName, listFailedCommands, quantile, genAsimovToyFile, signalScaleFactor, rMin, rMax = args
    cmdArgs = GetHybridNewCommandArgs(workspace, mass, dirName, quantile, genAsimovToyFile, signalScaleFactor, True, rMin*0.9, rMax*1.1)
    quantileStr = str(quantile).replace(".", "p")
    if options.submitLimitJobsAfterGridGen:
        queue = "espresso"
        taskName = 'hybridNewLimits.M{}.{}'.format(mass, quantileStr)
        if signalScaleFactor != 1.0:
            taskName += '.signalScaleFactor{}'.format(round(signalScaleFactor, 6))
        quant = "quant{}.".format(quantile) if quantile > 0 else ""
        sigSFRound = round(signalScaleFactor, 6)
        condorDir = dirName.strip("/")+"/condor"
        globStr = condorDir+"/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}grid*.root".format(sigSFRound, mass, quant)
        inputFiles = GetFileList(globStr)
        # of course, now the condorDir part has to be removed
        inputFiles = [f[f.rfind("/")+1:] for f in inputFiles]
        inputFiles.append(str(workspace))
        combinedOutputFile = "higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}gridAll.root".format(sigSFRound, mass, quant)
        limitCmd = GetComputeLimitsFromGridCommand(workspace, mass, combinedOutputFile, quantile, sigSFRound, rMin, rMax)
        SubmitLimitJobsBatch(mass, dirName, listFailedCommands, quantile, signalScaleFactor, rMin, rMax, taskName, "", limitCmd, queue, inputFiles)
    else:
        queue = "longlunch"
        taskName = 'hybridNewGridGen.M{}.{}'.format(mass, quantileStr)
        if signalScaleFactor != 1.0:
            taskName += '.signalScaleFactor{}'.format(round(signalScaleFactor, 6))
        inputFiles = [str(workspace)]
        if len(str(genAsimovToyFile)):
            inputFiles.append(str(genAsimovToyFile))
        SubmitLimitJobsBatch(mass, dirName, listFailedCommands, quantile, signalScaleFactor, rMin, rMax, taskName, cmdArgs, "", queue, inputFiles)


def SubmitAsymptoticBatch(args):
    workspace, mass, dirName, listFailedCommands, blinded, signalScaleFactor, index, rMin, rMax = args
    cmdArgs = GetAsymptoticCommandArgs(workspace, mass, dirName, blinded, signalScaleFactor, rMin, rMax)
    # cmd = "combine " + cmdArgs
    taskName = 'asymptotic.M{}'.format(mass)
    if signalScaleFactor != 1.0:
        taskName += '.signalScaleFactor{}'.format(round(signalScaleFactor, 6))
        taskName += '.betaIndex{}'.format(index)
    queue = "espresso"
    inputFiles = [str(workspace)]
    SubmitLimitJobsBatch(mass, dirName, listFailedCommands, -1, signalScaleFactor, rMin, rMax, taskName, cmdArgs, "", queue, inputFiles)


# http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#nuisance-parameter-impacts
def MakeImpacts(workspace, mass, dirName, signalScaleFactor):
    impactsDir = dirName.strip("/")+"/impacts"
    if not os.path.isdir(impactsDir):
        print("INFO: Making directory", impactsDir, flush=True)
        Path(impactsDir).mkdir(exist_ok=True)
    cmd += "cd {} && combineTool.py ".format(cwd + "/" + impactsDir)
    cmd += " -M Impacts -d {} -m {} --doInitialFit --robustFit 1".format(workspace, mass)
    cmd += commonCombineArgs.format(signalScaleFactor)
    RunCommand(cmd, dirName, None, True)
    cmd += "cd {} && combineTool.py ".format(cwd + "/" + impactsDir)
    cmd += " -M Impacts -d {} -m {} --robustFit 1 --doFits".format(workspace, mass)
    cmd += commonCombineArgs.format(signalScaleFactor)
    RunCommand(cmd, dirName, None, True)
    cmd += "cd {} && combineTool.py ".format(cwd + "/" + impactsDir)
    cmd += " -M Impacts -d {} -m {} -o impacts.m{}.json".format(workspace, mass, mass)
    cmd += commonCombineArgs.format(signalScaleFactor)
    RunCommand(cmd, dirName, None, True)
    cmd += "cd {} && plotImpacts.py ".format(cwd + "/" + impactsDir)
    cmd += " -i impacts.m{}.json -o impacts.m{}".format(mass, mass)
    RunCommand(cmd, dirName, None, True)

    
def RunHybridNewInteractive(workspace, mass, dirName, quantile=-1, genAsimovToyFile="", signalScaleFactor=1.0):
    cmd = GetHybridNewCommandArgs(workspace, mass, dirName, quantile, genAsimovToyFile, signalScaleFactor)
    cmd = 'combine ' + cmd
    RunCommand(cmd, dirName)
    quantileString = "quant{:.3f}".format(quantile) if quantile > 0 else "."
    globString = dirName+'/higgsCombineTest.HybridNew.mH{}.*.{}.root'.format(mass, quantileString)
    return FindFile(globString)


def RunAsymptoticInteractive(workspace, mass, dirName, batch, blinded=True, signalScaleFactor=1.0, rMin=-1, rMax=150):
    cmd = GetAsymptoticCommandArgs(workspace, mass, dirName, blinded, signalScaleFactor, rMin, rMax)
    cmd = 'combine ' + cmd
    RunCommand(cmd, dirName)
    return Path(sorted(glob.glob(dirName+'/higgsCombineTest.AsymptoticLimits.mH{}.*.root'.format(mass)), key=os.path.getmtime)[-1]).resolve()


def GetComputeLimitsFromGridCommand(workspace, mass, filename, quantile, signalScaleFact, rMin=1, rMax=200):
    rAbsAcc = 0.00001
    rRelAcc = 0.005
    if signalScaleFact != -1:
        sigSFStr = "signalScaleFactor"+str(signalScaleFact)+"_"
    else:
        sigSFStr = ""
    plotFilename = "limit_scan_{}m{}_quant{}.pdf".format("signalScaleFactor"+str(signalScaleFact)+"_" if signalScaleFact != -1 else "", mass, quantile)
    cmd = 'combine'
    cmd += ' {}'.format(workspace)
    # cmd += ' -v1'
    cmd += ' -M HybridNew'
    if rMin > 0:
        cmd += ' --rMin {}'.format(rMin)
    cmd += ' --rMax {}'.format(rMax)
    cmd += ' --LHCmode LHC-limits'
    cmd += ' --readHybridResults'
    if len(filename):
        cmd += ' --grid={}'.format(filename)
    else:
        cmd += ' --grid={}'
    cmd += ' --expectedFromGrid {}'.format(quantile)
    cmd += ' -m {}'.format(mass)
    cmd += ' --rAbsAcc {}'.format(rAbsAcc)
    cmd += ' --rRelAcc {}'.format(rRelAcc)
    cmd += ' --plot={}'.format(plotFilename)
    if signalScaleFact != -1:
        cmd += ' -n .signalScaleFactor{}'.format(signalScaleFact)
    cmd += commonCombineArgs.format(signalScaleFact)
    return cmd


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
    if not os.path.isfile(rootFile):
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
    massList = sorted(list(xsecLimitsByMassAndQuantile.keys()))
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
    minBetasPerMass = {300: 0, 400: 0, 500: 0.05, 600: 0.05, 700: 0.1, 800: 0.12, 900: 0.15, 1000: 0.2, 1100: 0.25, 1200: 0.3, 1300: 0.45, 1400: 0.55, 1500: 0.65}
    maxBetasPerMass = {300: 0.1, 400: 0.1, 500: 0.15, 600: 0.2, 700: 0.2, 800: 0.25, 900: 0.3, 1000: 0.35, 1100: 0.45, 1200: 0.55, 1300: 0.75, 1400: 1.0, 1500: 1.0}
    if mass > 1500:
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
    # print("\tChecking error file {}...".format(errFileName), end = "", flush=True)
    messagesToIgnore = ["[WARNING] Minimization finished with status 1 (covariance matrix forced positive definite), this could indicate a problem with the minimum!"]
    messagesToIgnore.append("Info in <TCanvas::Print>: pdf file limit_scan")
    typesOfIgnoredMessages = ["minimization status 1 warning", "limit scan pdf created"]
    showIgnoreMessages = [True, False]
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
                    print(colored("Found unexpected content '{}' in {}".format(line, errFileName), "red"))
                    if throwException:
                        raise RuntimeError("Found unexpected content '{}' in {}".format(line, errFileName))
    for idx, count in enumerate(instancesOfIgnoredMessages):
        if count > 0 and showIgnoreMessages[idx]:
            print(colored("Found {} error message(s) of type {} in {}".format(count, typesOfIgnoredMessages[idx], errFileName), "red"))
    # print("OK")


def CheckForErrorAndResultFiles(massList, condorDir):
    allFilesPresent = True
    for mass in massList:
        for quantileExp in quantilesExpected:
            quantile = "quant{:.3f}".format(quantileExp) if quantileExp > 0 else "."
            quantileStr = str(quantileExp).replace(".", "p")
            globString = condorDir+'/error/hybridNewLimits.M{}.{}.*.0.err'.format(mass, quantileStr)
            try:
                errFileName = FindFile(globString)
            except RuntimeError as e:
                print("Caught exception trying to find error file: '{}'".format(e))
                allFilesPresent = False
                continue
            globString = condorDir+'/higgsCombine.signalScaleFactor1.0.HybridNew.mH{}.{}.root'.format(mass, quantile)
            try:
                rootFileName = FindFile(globString)
            except RuntimeError as e:
                print("Caught exception: trying to find root file '{}'".format(e))
                allFilesPresent = False
    return allFilesPresent


def ReadBatchResults(massList, condorDir):
    if not CheckForErrorAndResultFiles(massList, condorDir):
        return {}, {}
    for mass in massList:
        listOfWorkspaceFiles = sorted(glob.glob(dirName + "/datacards/*.m{}.root".format(mass)), key=os.path.getmtime)
        if len(listOfWorkspaceFiles) != 1:
            raise RuntimeError("Globbing for {} did not result in one file as expected, but {}".format(workspaceFileName, listOfWorkspaceFiles))
        rLimitsByMassAndQuantile[mass] = manager.dict()
        xsecLimitsByMassAndQuantile[mass] = manager.dict()
        signalScaleFactorsByMassAndQuantile[mass] = manager.dict()
        cardWorkspace = Path(listOfWorkspaceFiles[-1]).resolve()
        for quantileExp in quantilesExpected:
            signalScaleFactorsByMassAndQuantile[mass][str(quantileExp)] = {}
            quantile = "quant{:.3f}".format(quantileExp) if quantileExp > 0 else "."
            quantileStr = str(quantileExp).replace(".", "p")
            globString = condorDir+'/error/hybridNewLimits.M{}.{}.*.0.err'.format(mass, quantileStr)
            try:
                errFileName = FindFile(globString)
            except RuntimeError as e:
                print("Caught exception: {}".format(e.what()))
            CheckErrorFile(errFileName, False)  # don't throw exception
            # # now run HybridNew again to get limits from grid
            # globString = condorDir+'/higgsCombine.signalScaleFactor*.HybridNew.mH{}.quant{}.grid.root'.format(mass, quantileExp)
            # rootFileName = FindFile(globString)
            # resultFileName = ComputeLimitsFromGrid(cardWorkspace, mass, condorDir, rootFileName, quantileExp)
            # # now extract limits from the root file produced above
            # limit, limitErr, quantileFromFile, signalScaleFactor = ExtractLimitResult(resultFileName)
            globString = condorDir+'/higgsCombine.signalScaleFactor1.0.HybridNew.mH{}.{}.root'.format(mass, quantile)
            rootFileName = FindFile(globString)
            limit, limitErr, quantileFromFile, signalScaleFactor = ExtractLimitResult(rootFileName)
            rLimitsByMassAndQuantile[mass][str(quantileExp)] = limit
            xsecLimitsByMassAndQuantile[mass][str(quantileExp)] = limit * xsThByMass[float(mass)]
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
                raise RuntimeError("Globbing for {} did not result in one file as expected, but {}".format(workspaceFileName, listOfWorkspaceFiles))
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
                betaValRound = round(betaVal, 6)
                if not betaValRound in rLimitsByMassAndQuantile[str(quantileExp)].keys():
                    rLimitsByMassAndQuantile[str(quantileExp)][betaValRound] = {}
                    xsecLimitsByMassAndQuantile[str(quantileExp)][betaValRound] = {}
                    signalScaleFactorsByMassAndQuantile[str(quantileExp)][betaValRound] = {}
                rLimitsByMassAndQuantile[str(quantileExp)][betaValRound][mass] = limit
                xsecLimitsByMassAndQuantile[str(quantileExp)][betaValRound][mass] = limit * xsThByMass[float(mass)]
                signalScaleFactorsByMassAndQuantile[str(quantileExp)][betaValRound][mass] = signalScaleFactor
                # print("\tINFO: ReadBatchResultsBetaScan() - fill xsecLimitsByMassAndQuantile[{}][{}][{}]".format(quantileExp, betaValRound, mass))
    return xsecLimitsByMassAndQuantile


def ReadBatchResultParallel(args):
    cardWorkspace, mass, condorDir, errFile, quantileExp = args
    quantile = "quant{:.3f}".format(quantileExp) if quantileExp > 0 else "."
    quantile = str(quantile).rstrip("0")
    quantileStr = str(quantileExp).replace(".", "p").rstrip("0")
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
    betaValRound = round(betaVal, 6)
    print("INFO: [3] ReadBatchResultParallel(): store limit*xs[{}]={}*{}={} in rLimitsByMassAndQuantile={} for key={} at betaValRound={} from resultFileName={}".format(
        mass, limit, xsThByMass[float(mass)], limit * xsThByMass[float(mass)], dict(rLimitsByMassAndQuantile), str(quantileExp), betaValRound, resultFileName), flush=True)
    rLimitsByMassAndQuantile[str(quantileExp)][betaValRound][mass] = limit
    xsecLimitsByMassAndQuantile[str(quantileExp)][betaValRound][mass] = limit * xsThByMass[float(mass)]
    signalScaleFactorsByMassAndQuantile[str(quantileExp)][betaValRound][mass] = signalScaleFactor
    print("INFO: [5] ReadBatchResultParallel(): DONE storing xsecLimitsByMassAndQuantile[{}]={}".format(str(quantileExp), dict(xsecLimitsByMassAndQuantile)), flush=True)


def ReadBatchResultsBetaScan(massList, condorDir):
    ncores = 8
    with multiprocessing.Pool(ncores) as pool:
        nJobs = 0
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
                quantile = str(quantile).rstrip("0")
                globString = condorDir+'/error/hybridNewLimits.M{}.{}.*.0.err'.format(mass, quantileStr)
                errorFiles = GetFileList(globString)
                for errFile in errorFiles:
                    lastPart = errFile.split("{}.".format(quantileStr))[-1]
                    lastPartCut = lastPart.rstrip(".0.err")
                    sigScaleFact = lastPartCut[0:lastPartCut.rfind(".")].strip("signalScaleFactor")
                    sigScaleFactRound = round(float(sigScaleFact), 6)
                    rootGlobString = condorDir+'/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}.grid.root'.format(sigScaleFactRound, mass, quantile)
                    rootFileName = FindFile(rootGlobString)
                    limit, limitErr, quantileFromFile, signalScaleFactor = ExtractLimitResult(rootFileName)
                    betaVal = math.sqrt(signalScaleFactor)
                    betaValRound = round(betaVal, 6)
                    if not betaValRound in rLimitsByMassAndQuantile[str(quantileExp)].keys():
                        rLimitsByMassAndQuantile[str(quantileExp)][betaValRound] = manager.dict()
                    if not betaValRound in xsecLimitsByMassAndQuantile[str(quantileExp)].keys():
                        xsecLimitsByMassAndQuantile[str(quantileExp)][betaValRound] = manager.dict()
                    if not betaValRound in signalScaleFactorsByMassAndQuantile[str(quantileExp)].keys():
                        signalScaleFactorsByMassAndQuantile[str(quantileExp)][betaValRound] = manager.dict()
                    try:
                        pool.apply_async(ReadBatchResultParallel, [[cardWorkspace, mass, condorDir, errFile, quantileExp]])
                        nJobs += 1
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
        if nJobs > 0:
            print("Waiting for {} limit calculations to finish for all beta values...".format(nJobs))
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
    # print("\tINFO: loggraph() - x={}, y={}".format(x, y))
    g = TGraph(len(x), x, y)
    return g


def logspline(inputarrayX,inputarrayY):
    logarray = []
    for j in inputarrayY:
        logarray.append(math.log(j))
    x = array("d",inputarrayX)
    y = array("d",logarray)
    g = TGraph(len(x),x,y)
    outspline = TSpline3("",g)
    return outspline


def get_simple_intersection(graph1, graph2, xmin, xmax):
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
        # print("\tINFO: get_simple_intersection() - math.exp(graph1.Eval(x)) =", math.exp(graph1.Eval(x)))
        # print("\tINFO: get_simple_intersection() - graph2.Eval(x) =", graph2.Eval(x))
        # print("\tINFO: get_simple_intersection() - math.exp(graph2.Eval(x)) =", math.exp(graph2.Eval(x)))
        # print("\tINFO: get_simple_intersection() - math.exp(graph1.Eval(x)) - math.exp(graph2.Eval(x)) =", math.exp(graph1.Eval(x)) - math.exp(graph2.Eval(x)))
        #
        # thisdif = (math.exp(graph1.Eval(x)) - math.exp(graph2.Eval(x)))
        # xx.append(math.exp(graph1.Eval(x)))
        # yy.append(math.exp(graph2.Eval(x)))
        #XXX SIC TEST Jun. 26 2024: use 10^x instead of exp(x)
        thisdif = (math.pow(graph1.Eval(x), 10) - math.pow(graph2.Eval(x), 10))
        xx.append(math.pow(graph1.Eval(x), 10))
        yy.append(math.pow(graph2.Eval(x), 10))
        #
        sdif.append(thisdif)
        dif.append(abs(thisdif))
        xvals.append(x)
        # print("\tINFO: get_simple_intersection() [1] -", str(x) + '   ' +str(math.exp(graph1.Eval(x))) + '    '+str(math.exp(graph2.Eval(x))) + '    ' + str(thisdif))
        x = x+inc
	#print 'Done Looping for Difs'
    mindif = min(dif)
    bestmass = 0	
    
    for x in range(len(dif)-2):
       a = sdif[x]
       b = sdif[x+1]
       # print("\tINFO: get_simple_intersection() [2] -", str(xvals[x+1]) +'    '+str(a)  + '     ' +str(b))
       if ((a/abs(a))*(b/abs(b))) < 0.0 and a > 0.0 :
           # print('\tINFO: get_simple_intersection() - Limit found at: '+ (str(xvals[x])))
           bestmass = xvals[x]
           break
    return [bestmass, mindif]


def ComputeBetaLimits(xsThByMass, xsecLimitsByMassAndQuantile):
    mTh = np.array(list(xsThByMass.keys()), dtype="f")
    xsTh = np.array(list(xsThByMass.items()), dtype="f")
    g = TGraph(len(mTh), mTh, xsTh);
    # spline = TSpline3("xsection", g)
    logtheory = loggraph(list(xsThByMass.keys()), list(xsThByMass.values()))  # make graph with log10 of y values
    massLimitsByQuantileAndBetaVal = {}
    for quantile in xsecLimitsByMassAndQuantile.keys():
        print("INFO: ComputeBetaLimits() - examine quantile={}".format(quantile))
        massLimitsByQuantileAndBetaVal[quantile] = {}
        for betaVal in xsecLimitsByMassAndQuantile[quantile].keys():
            print("\tINFO: examine betaVal={}".format(betaVal))
            limit_set = list(xsecLimitsByMassAndQuantile[quantile][betaVal].values())
            massList = list(xsecLimitsByMassAndQuantile[quantile][betaVal].keys())
            if len(massList) < 2:
                print("\tWARN: Skipping beta value={} as we only have one mass point tested here! Need to adjust beta scan range.".format(betaVal))
                continue
            print("\tINFO: examine massList={}, limit_set={}".format(massList, limit_set))
            # print("\tINFO: make loggraph for fitted_limits")
            fitted_limits = loggraph(massList, limit_set)
            #FIXME: get_simple_intersection seems to assume both graphs are in log base e rather than log10, as they actually are
            # print("\tINFO: ComputeBetaLimits() - get_simple_intersection")
            goodm = get_simple_intersection(logtheory, fitted_limits, min(mTh), max(mTh))
            massLimitsByQuantileAndBetaVal[quantile][betaVal] = round(goodm[0], 3)
            print("\tINFO: for betaVal={}: bestmass, mindif={}".format(betaVal, goodm))
            if goodm[0] < min(massList) or goodm[0] > max(massList):
                print("\tWARN: For beta value={}, intersection mass is outside of massList range {}! Need to adjust beta scan range.".format(betaVal, massList))
    return massLimitsByQuantileAndBetaVal



def CreateComboArrays(xsecLimitsByMassAndQuantile):
    retVal = {}
    for quantile in xsecLimitsByMassAndQuantile.keys():
        retVal[quantile] = {}
        retVal[quantile]["betas"] = []
        retVal[quantile]["massLimits"] = []
        sortedBetas = list(sorted(xsecLimitsByMassAndQuantile[quantile].keys()))
        # for betaVal, mass in xsecLimitsByMassAndQuantile[quantile].items():
        for betaVal in sortedBetas:
            retVal[quantile]["betas"].append(betaVal)
            retVal[quantile]["massLimits"].append(xsecLimitsByMassAndQuantile[quantile][betaVal])
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


def MakeResultTableBetaScan(masses, xsecLimitsByMassAndQuantile, quantileExp="0.5"):
    table = []
    headers = ["MLQ"]
    betaValsSorted = sorted(xsecLimitsByMassAndQuantile[quantileExp].keys())
    for mIdx, mass in enumerate(masses):
        tableRow = [mass]
        for betaVal in betaValsSorted:
            if mIdx == 0:
                headers.append(str(betaVal))
            if mass not in xsecLimitsByMassAndQuantile[quantileExp][betaVal].keys():
                continue  # not all masses have all beta values
            xsecLimit = xsecLimitsByMassAndQuantile[quantileExp][betaVal][mass]
            tableRow.append(xsecLimit)
        table.append(tableRow)
    return table, headers


####################################################################################################
# Run
####################################################################################################
if __name__ == "__main__":
    doBatch = True
    gridPointsPerJob = 10
    gridScanPoints = 50
    doShapeBasedLimits = False
    doAsymptoticLimits = False
    blinded = True
    massList = list(range(300, 3100, 100))
    # massList = list(range(300, 2200, 100))
    betasToScan = list(np.linspace(0.002, 1, 500))[:-1] + [0.9995]
    
    quantilesExpected = [0.025, 0.16, 0.5, 0.84, 0.975]
    xsThFilename = "$LQANA/config/xsection_theory_13TeV_scalarPairLQ.txt"
    # xsThFilename = "$LQANA/config/xsection_theory_13TeV_stopPair.txt"
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
        "-s",
        "--skipStandardLimits",
        dest="skipStandardLimits",
        help="don't submit standard limits",
        metavar="skipStandardLimits",
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
        "-f",
        "--fileWithRValueRanges",
        dest="fileWithRValueRanges",
        help="file containing r-value ranges",
        metavar="fileWithRValueRanges",
        default=None,
    )
    parser.add_option(
        "-j",
        "--fileWithRValueRangesBetaScan",
        dest="fileWithRValueRangesBetaScan",
        help="file containing r-value ranges for beta scan",
        metavar="fileWithRValueRangesBetaScan",
        default=None,
    )
    (options, args) = parser.parse_args()
    if options.datacard is None and options.readResults is None:
        raise RuntimeError("Need either option d to specify datacard, or option r to specify reading limit results from batch")
    if options.name is None:
        raise RuntimeError("Option n to specify name of limit results dir is required")
    if options.doBetaScan and not doBatch:
        raise RuntimeError("Won't do beta scan without batch submission enabled (see doBatch parameter inside script).")
    
    manager = multiprocessing.Manager()
    dictAsimovToysByScaleFactor = manager.dict()
    rLimitsByMassAndQuantile = manager.dict()
    xsecLimitsByMassAndQuantile = manager.dict()
    signalScaleFactorsByMassAndQuantile = manager.dict()
    listFailedCommands = manager.list()
    dirName = options.name
    xsThByMass, _, _ = ReadXSecFile(xsThFilename)
    failedBatchCommands = []
    combinedDatacard = options.datacard
    if not options.readResults:
        startTime = time.time()
        #TODO: check for previous datacards?
        separateDatacardsDir = dirName+"/datacards"
        asimovToysDir = dirName+"/asimovData"
        if not os.path.isdir(dirName):
            print("INFO: Making directory", dirName, flush=True)
            Path(dirName).mkdir(exist_ok=True)
        if not os.path.isdir(separateDatacardsDir):
            print("INFO: Making directory", separateDatacardsDir, flush=True)
            Path(separateDatacardsDir).mkdir(exist_ok=True)
        if not os.path.isdir(asimovToysDir):
            print("INFO: Making directory", asimovToysDir, flush=True)
            Path(asimovToysDir).mkdir(exist_ok=True)
        
        if not doShapeBasedLimits:
            massListFromCards, cardFilesByMass, _ = SeparateDatacards(combinedDatacard, 0, separateDatacardsDir)
        for mass in massList:
            cardFile = combinedDatacard if doShapeBasedLimits else cardFilesByMass[mass]
            if not options.estimateRValueScanRange:
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
                    if not options.skipStandardLimits:
                        rValuesByMassAndQuantile = {}
                        if doBatch:
                            scanDirName = dirName+"/asymptoticLimits"
                            if not os.path.isdir(scanDirName):
                                print("INFO: Making directory", scanDirName, flush=True)
                                Path(scanDirName).mkdir(exist_ok=True)
                            MakeCondorDirs(scanDirName)
                            SubmitAsymptoticBatch((cardWorkspace, mass, scanDirName, listFailedCommands, blinded, 1.0, -1, -1, 150))
                            continue
                        else:
                            print("INFO: Doing interactive AsymptoticLimits for mass {} to estimate scan range".format(mass), flush=True)
                            asymptoticRootFileName = RunAsymptoticInteractive(cardWorkspace, mass, dirName, blinded, 1.0)
                            limits, limitErrs, quantiles = ExtractAsymptoticLimitResultRoot(asymptoticRootFileName)
                            rValuesByQuantile = dict(zip(quantiles, limits))
                            #FIXME: add writing of r-value scan file here
                            rValuesByMassAndQuantile[mass] = dict(zip(quantiles, limits))
                    if options.doBetaScan:
                        rValuesByMassBetaAndQuantile = {}
                        betaDirName = dirName+"/betaScan/asymptoticLimits"
                        if not os.path.isdir(betaDirName):
                            print("INFO: Making directory", betaDirName, flush=True)
                            Path(betaDirName).mkdir(exist_ok=True, parents=True)
                        MakeCondorDirs(betaDirName)
                        betasToSubmit = GetBetasToSubmit(mass)
                        rValuesAtBeta = {}
                        if doBatch:
                            ncores = 8
                            with multiprocessing.Pool(ncores) as pool:
                                asimovJobs = 0
                                for idx, beta in enumerate(betasToSubmit):
                                    signalScaleFactor = beta*beta
                                    try:
                                        pool.apply_async(SubmitAsymptoticBatch, [[cardWorkspace, mass, betaDirName, listFailedCommands, blinded, signalScaleFactor, idx, -1, 150]])
                                        asimovJobs += 1
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
                                if asimovJobs > 0:
                                    print("Waiting for {} AsymptoticLimits job submission to finish for all beta values for mass={}...".format(asimovJobs, mass))
                                pool.join()
                            continue
                        else:
                            print("INFO: Doing interactive AsymptoticLimits for mass {}, beta {} to estimate scan range".format(mass, beta), flush=True)
                            asymptoticRootFileName = RunAsymptoticInteractive(cardWorkspace, mass, betaDirName, blinded, signalScaleFactor)
                            limits, limitErrs, quantiles = ExtractAsymptoticLimitResultRoot(asymptoticRootFileName)
                            rValuesByQuantile = dict(zip(quantiles, limits))
                            rValuesAtBeta[beta] = rValuesByQuantile
                            rValuesByMassBetaAndQuantile[mass] = rValuesAtBeta
                elif not options.skipStandardLimits:
                    if options.fileWithRValueRanges is not None:
                        with open(options.fileWithRValueRanges) as f:
                            rValuesByMassAndQuantile = json.load(f)
                    else:
                        raise RuntimeError("Standard limits are requested, but no r-value scan range file has been passed with the '-f' option.")
                elif options.doBetaScan:
                    if options.fileWithRValueRangesBetaScan is not None:
                        with open(options.fileWithRValueRangesBetaScan) as f:
                            rValuesByMassBetaAndQuantile = json.load(f)
                    else:
                        raise RuntimeError("Beta scan limit requested ('-b'), but no r-value scan range file has been passed with the '-j' option.")
                else:
                    raise RuntimeError("Either the r-value scan range must be estimated first with option '-e', or there must be a file of scan ranges given with '-f' in order to submit limit jobs")

                # now submit jobs for limits
                for quantileExp in quantilesExpected:
                    signalScaleFactor = 1.0
                    if not options.skipStandardLimits:
                        asimovToyFile = FindAsimovToyData(mass, signalScaleFactor, asimovToysDir)
                        if asimovToyFile is not None:
                            print("INFO: Using previously-generated Asimov toy file: {}".format(asimovToyFile), flush=True)
                        else:
                            asimovToyFile = GenerateAsimovToyData((cardWorkspace, mass, asimovToysDir, dictAsimovToysByScaleFactor, signalScaleFactor))
                    if doBatch:
                        MakeCondorDirs(dirName)
                        if not options.skipStandardLimits:
                            rValuesByQuantile = rValuesByMassAndQuantile[str(mass)]
                            if quantileExp == 0.025 and mass > 800:  # adjust scan range upwards for lowest quantile and higher masses
                                rMax = rValuesByQuantile[str(quantileExp)]*1.8
                                rMin = rValuesByQuantile[str(quantileExp)]*0.8
                            elif quantileExp == 0.016 and mass > 1000:
                                rMax = rValuesByQuantile[str(quantileExp)]*1.9
                                rMin = rValuesByQuantile[str(quantileExp)]*0.8
                            # elif quantileExp == 0.975 and mass > 1500:  # adjust scan range downwards here
                            #     rMax = rValuesByQuantile[str(quantileExp)]*1.0
                            #     rMin = rValuesByQuantile[str(quantileExp)]*0.45
                            else:
                                rMax = rValuesByQuantile[str(quantileExp)]*1.3
                                rMin = rValuesByQuantile[str(quantileExp)]*0.75
                            failedCmd = SubmitHybridNewBatch((cardWorkspace, mass, dirName, listFailedCommands, quantileExp, asimovToyFile, signalScaleFactor, rMin, rMax))
                            if failedCmd is not None:
                                failedBatchCommands.append(failedCmd)
                        # beta scan jobs
                        if options.doBetaScan:
                            betaDirName = dirName+"/betaScan"
                            MakeCondorDirs(betaDirName)
                            rValuesAtBeta = rValuesByMassBetaAndQuantile[str(mass)]
                            ncores = 8
                            with multiprocessing.Pool(ncores) as pool:
                                asimovJobs = 0
                                for idx, beta in enumerate(rValuesAtBeta.keys()):
                                    signalScaleFactor = float(beta)*float(beta)
                                    asimovToysBetaDir = betaDirName+"/asimovToys"
                                    if not os.path.isdir(asimovToysBetaDir):
                                        print("INFO: beta scan - Making directory", asimovToysBetaDir, flush=True)
                                        Path(asimovToysBetaDir).mkdir(exist_ok=True)
                                    asimovToyFile = FindAsimovToyData(mass, signalScaleFactor, asimovToysBetaDir)
                                    if asimovToyFile is not None:
                                        # print("INFO: beta scan - Using previously-generated Asimov toy file: {}".format(asimovToyFile), flush=True)
                                        pass
                                    else:
                                        try:
                                            pool.apply_async(GenerateAsimovToyData, [[cardWorkspace, mass, asimovToysBetaDir, dictAsimovToysByScaleFactor, signalScaleFactor]])
                                            asimovJobs += 1
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
                                if asimovJobs > 0:
                                    print("Waiting for {} Asimov toy generation jobs to finish for all beta values for mass={}, quantile={}...".format(asimovJobs, mass, quantileExp))
                                pool.join()
                            # HybridNew
                            with multiprocessing.Pool(ncores) as pool:
                                hybridNewJobs = 0
                                rValuesAtBeta = rValuesByMassBetaAndQuantile[str(mass)]
                                for idx, beta in enumerate(rValuesAtBeta.keys()):
                                    signalScaleFactor = float(beta)*float(beta)
                                    if quantileExp == 0.025 and mass > 800:  # adjust scan range upwards for lowest quantile and higher masses
                                        rMax = rValuesAtBeta[str(beta)][str(quantileExp)]*1.8
                                        rMin = rValuesAtBeta[str(beta)][str(quantileExp)]*0.8
                                    elif quantileExp == 0.016 and mass > 1000:
                                        rMax = rValuesAtBeta[str(beta)][str(quantileExp)]*1.9
                                        rMin = rValuesAtBeta[str(beta)][str(quantileExp)]*0.8
                                    # elif quantileExp == 0.975 and mass > 1500:  # adjust scan range downwards here
                                    #     rMax = rValuesAtBeta[str(beta)][str(quantileExp)]*1.0
                                    #     rMin = rValuesAtBeta[str(beta)][str(quantileExp)]*0.45
                                    else:
                                        rMax = rValuesAtBeta[str(beta)][str(quantileExp)]*1.3
                                        rMin = rValuesAtBeta[str(beta)][str(quantileExp)]*0.75
                                    try:
                                        asimovToyFile = FindAsimovToyData(mass, signalScaleFactor, asimovToysBetaDir)
                                        pool.apply_async(SubmitHybridNewBatch, [[cardWorkspace, mass, betaDirName, listFailedCommands, quantileExp, asimovToyFile, signalScaleFactor, rMin, rMax]])
                                        hybridNewJobs += 1
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
                                print("Waiting for submission of {} HybridNew jobs for mass={}, quantile={}...".format(hybridNewJobs, mass, quantileExp))
                                pool.join()
                                # if failedCmd is not None:
                                #     failedBatchCommands.append(failedCmd)
                    else:
                        hnStartTime = time.time()
                        rootFileName = RunHybridNewInteractive(cardWorkspace, mass, dirName, quantile=quantileExp, genAsimovToyFile=asimovToyFile, signalScaleFactor=signalScaleFactor)
                        hnStopTime = time.time()
                        execTimeStr = GetExecTimeStr(hnStartTime, hnStopTime)
                        print("HybridNew calculation execution time:", execTimeStr, flush=True)
                        limit, limitErr, quantile, signalScaleFactor = ExtractLimitResult(rootFileName)
                        rLimitsByMassAndQuantile[mass][str(quantileExp)] = limit
                        xsecLimitsByMassAndQuantile[mass][str(quantileExp)] = limit * xsThByMass[float(mass)]
                        signalScaleFactorsByMassAndQuantile[mass][str(quantileExp)] = signalScaleFactor
        
    if len(failedBatchCommands):
        print(colored("batch commands failed:", "red"))
        for item in failedBatchCommands:
            print(item)
    if len(listFailedCommands):
        print(colored("batch commands failed:", "red"))
        for item in listFailedCommands:
            print(item)

    if options.readResults and options.estimateRValueScanRange:
        if not options.skipStandardLimits:
            rValuesByMassAndQuantile = {}
            scanDirName = dirName+"/asymptoticLimits/condor"
            for mass in massList:
                globString = scanDirName+'/error/asymptotic.M{}.*.0.err'.format(mass)
                errFiles = GetFileList(globString)
                for f in errFiles:
                    CheckErrorFile(f)
                    rootFile = FindFile(scanDirName+"/higgsCombineTest.AsymptoticLimits.mH{}.*.root".format(mass))
                    limits, limitErrs, quantiles, _ = ExtractAsymptoticLimitResultRoot(rootFile)
                    rValuesByMassAndQuantile[mass] = dict(zip(quantiles, limits))
            with open(scanDirName+'/rValues.json', 'w') as f:
                json.dump(rValuesByMassAndQuantile, f)
        if options.doBetaScan:
            rValuesByMassBetaAndQuantile = {}
            scanDirName = dirName+"/betaScan/asymptoticLimits/condor"
            for mass in massList:
                betasToSubmit = GetBetasToSubmit(mass)
                rValuesByMassBetaAndQuantile[mass] = {}
                globString = scanDirName+'/error/asymptotic.M{}.*.0.err'.format(mass)
                errFiles = GetFileList(globString)
                for f in errFiles:
                    CheckErrorFile(f)
                fileList = GetFileList(scanDirName+"/higgsCombine*.AsymptoticLimits.mH{}.*.root".format(mass))
                for f in fileList:
                    betaValIndex = int(f[f.find("betaIndex")+9:f.find(".", f.find("betaIndex")+9)])
                    limits, limitErrs, quantiles, signalScaleParam = ExtractAsymptoticLimitResultRoot(f)
                    rValuesByMassBetaAndQuantile[mass][betasToSubmit[betaValIndex]] = dict(zip(quantiles, limits))
            with open(scanDirName+'/rValues.json', 'w') as f:
                json.dump(rValuesByMassBetaAndQuantile, f)
    elif not doBatch or options.readResults and not options.estimateRValueScanRange:
        years, intLumi = GetYearAndIntLumiFromDatacard(combinedDatacard)
        intLumi = intLumi/1000.0
        if intLumi < 100:
            intLumi = round(intLumi, 1)
        else:
            intLumi = round(intLumi, 0)
        print("Found total int lumi = {}/fb for years = {}".format(intLumi, years))
        if not options.skipStandardLimits:
            if options.readResults:
                condorDir = dirName.rstrip("/")+"/condor"
                print("INFO: Reading results from batch from {}...".format(condorDir), flush=True)
                xsecLimitsByMassAndQuantile, signalScaleFactorsByMassAndQuantile = ReadBatchResults(massList, condorDir)
                if len(xsecLimitsByMassAndQuantile) and len(signalScaleFactorsByMassAndQuantile):
                    print("DONE", flush=True)
                else:
                    raise RuntimeError("Got empty dicts back from ReadBatchResults; something wrong happened there.")
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
            BR_Sigma_EE_vsMass(dirName, intLumi, masses, shadeMasses, xsMedExp, xsObs, xsOneSigmaExp, xsTwoSigmaExp)
        if options.doBetaScan:
            mainDirName = dirName
            betaDirName = dirName+"/betaScan"
            condorDir = betaDirName.strip("/")+"/condor"
            print("INFO: Reading results from batch from {}...".format(condorDir), flush=True, end="")
            ReadBatchResultsBetaScan(massList, condorDir)
            print("DONE", flush=True)
            table, columnNames = MakeResultTableBetaScan(massList, xsecLimitsByMassAndQuantile)
            print("xsecLimits beta scan for quantile=0.5:")
            print(tabulate(table, headers=columnNames, tablefmt="fancy_grid", floatfmt=".6f"))
            print("xsecLimitsByMassAndQuantile={}".format(dict(xsecLimitsByMassAndQuantile)))
            print("xsecLimitsByMassAndQuantile at quantile 0.5={}".format(dict(xsecLimitsByMassAndQuantile["0.5"])))
            print("INFO: Compute beta limits...", flush=True, end="")
            massLimitsByQuantileAndBetaVal = ComputeBetaLimits(xsThByMass, xsecLimitsByMassAndQuantile)
            print("massLimitsByQuantileAndBetaVal={}".format(massLimitsByQuantileAndBetaVal))
            print("DONE", flush=True)
            comboArraysByQuantile = CreateComboArrays(massLimitsByQuantileAndBetaVal)
            print("comboArraysByQuantile['0.5']['betas']={}".format(comboArraysByQuantile["0.5"]["betas"]))
            print("comboArraysByQuantile['0.5']['massLimits']={}".format(comboArraysByQuantile["0.5"]["massLimits"]))
            ComboPlot(mainDirName, comboArraysByQuantile["0.5"]["betas"], comboArraysByQuantile["0.5"]["massLimits"]) # , m_observed_lljj)
        if options.doImpacts:
            print("INFO: Making nuisance parameter impact plots...", flush=True, end="")
            for mass in massList:
                datacardDir = dirName.strip("/")+"/datacards"
                cardWorkspace = FindCardWorkspace(datacardDir + "/*", mass, True)
                MakeImpacts(cardWorkspace, mass, dirName, signalScaleFactorsByMassAndQuantile[mass][str(quantilesExpected[0])])
            print("DONE", flush=True)
