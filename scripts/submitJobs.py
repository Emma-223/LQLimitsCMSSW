import sys
import os
import subprocess
import multiprocessing
#from rich.progress import Progress
from pathlib import Path
import time
from random import randint

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
    subprocess.check_output(["condor_submit",filename])
    finishedFilesList.append(f)
    progress_bar(len(finishedFilesList), nFiles,dirname)

if __name__ == "__main__":

    parentDir = sys.argv[1]

    dirList = os.listdir(parentDir)
    for d in dirList:
        if not os.path.isdir(parentDir+"/"+d):
            continue
            #pool.apply_async(SubmitJobsInDir,[d])
        if "0p5" in d:
            continue
        filesToSubmit = os.listdir(parentDir+"/"+d)
        nFiles = len(filesToSubmit)
        print("\nSubmit {} condor.sub files for {}".format(nFiles,d))
        manager = multiprocessing.Manager()
        finishedFilesList = manager.list()
        with multiprocessing.Pool(4) as pool:
            for i,f in enumerate(filesToSubmit):
                #SubmitJob(f,d)
            #progress_bar(i+1, nFiles, d)
                pool.apply_async(SubmitJob, [f,d])#, callback = progress_bar(i+1, nFiles, d))
            pool.close()
            pool.join()
    print("\nDone submitting jobs")
