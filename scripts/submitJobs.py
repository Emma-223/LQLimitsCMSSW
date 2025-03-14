import sys
import os
import subprocess
import multiprocessing
#from rich.progress import Progress
from pathlib import Path
import time
from random import randint

def GetOkCondorSites():
    siteList = ""
    with open("/etc/ciconnect/config.ini", "r") as siteConfFile:
        for line in siteConfFile:
            if "DefaultSites" in line:
                siteList = line.split("=")[-1].strip().split(",")
    sitesToExclude = ["T2_TR_METU", "T2_UA_KIPT", "T1_ES_PIC", "T2_TW_NCHC", "T2_US_MIT", "T2_UA_KIPT", "T2_ES_IFCA", "T2_EE_Estonia", "T2_ES_CIEMAT"]
    for s in sitesToExclude:
        if s in siteList:
            siteList.remove(s)
    okSitesCmd = "CONDOR_DEFAULT_DESIRED_SITES=\""+",".join(siteList)
    okSitesCmd = okSitesCmd.strip(",")+"\""
    print(okSitesCmd)
    return okSitesCmd

def progress_bar(i, tot, prefix='', suffix = '', length = 50, fill = '*'):
    percent = ("{0:.1f}").format(100* (i/float(tot)))
    filled = int(length * i //tot)
    bar = fill * filled + '-' * (length - filled -1)
    sys.stdout.write(f'\r{prefix} |{bar}|{percent}%{suffix} ')
    sys.stdout.flush()

def SubmitJob(f,dirname):
    #print("submit file "+parentDir+"/"+dirname+"/"+f+"/condor.sub")
    filename = parentDir+"/"+dirname+"/"+f+"/condor.sub"
    #time.sleep(randint(2,5)) # for troubleshooting the progress bar
    #subprocess.check_output([okSitesCmd,"condor_submit",filename])
    exitCode = os.WEXITSTATUS(os.system(okSitesCmd+" condor_submit "+filename+" > {}/{}/{}/condor_submit_log.txt".format(parentDir,dirname,f)))
    finishedFilesList.append(f)
    progress_bar(len(finishedFilesList), nFiles,dirname)

if __name__ == "__main__":

    parentDir = sys.argv[1]

    okSitesCmd = GetOkCondorSites()

    dirList = os.listdir(parentDir)
    for d in dirList:
        if not os.path.isdir(parentDir+"/"+d):
            continue
            #pool.apply_async(SubmitJobsInDir,[d])
        if not ("1800" in d):
            continue
        filesToSubmit = os.listdir(parentDir+"/"+d)
        filesToSubmit.sort()
        nFiles = len(filesToSubmit)
        print("\nSubmit {} condor.sub files for {}".format(nFiles,d))
        manager = multiprocessing.Manager()
        finishedFilesList = manager.list()
        with multiprocessing.Pool(8) as pool:
            for i,f in enumerate(filesToSubmit):
                if not "betaIds3" in f:
                    continue
                #SubmitJob(f,d)
            #progress_bar(i+1, nFiles, d)
                pool.apply_async(SubmitJob, [f,d])#, callback = progress_bar(i+1, nFiles, d))
            pool.close()
            pool.join()
    print("\nDone submitting jobs")
