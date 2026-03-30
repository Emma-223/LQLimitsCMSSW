import sys
import os
import subprocess
import glob
import multiprocessing

from submitJobs import GetOkCondorSites

def GetDirList(fileLines):
    dirList = []
    for filepath in fileLines:
        if filepath.startswith("#"):
            continue
        filepath = filepath.strip()
        filename = filepath.split("/")[-1]
        directory = filepath.replace(filename,"")
        if not directory in dirList:
            dirList.append(directory)
    return dirList

def ResubmitJob(directory):
    outFiles = glob.glob(directory+"condor.*.out")
    errFiles = glob.glob(directory+"condor.*.err")
    filesToRemove = outFiles+errFiles
    if len(filesToRemove):
        print("remove {} .out and {} .err files from ".format(len(outFiles),len(errFiles))+directory)
        cmd = ["rm"] + filesToRemove
        subprocess.run(cmd)
    subfile = directory+"condor.sub"
    nSubFiles = len(glob.glob(subfile))
    if not nSubFiles==1:
        print("ERROR: found {} condor.sub files in ".format(nSubFiles)+directory)
    #    exit()
    print("submit file "+subfile)
    subprocess.run(["condor_submit",subfile])
    #print("test")
    #print(okSitesCmd+" condor_submit "+subfile)
    #print("test")
    #exitCode = os.WEXITSTATUS(os.system(okSitesCmd+" condor_submit "+subfile))

if __name__=="__main__":
    manager = multiprocessing.Manager()
    txtfile = sys.argv[1]
    filesSubmitted = manager.list()
    with open(txtfile) as f:
        lines = f.readlines()

    okSitesCmd = GetOkCondorSites()

    dirList = GetDirList(lines)
    print("INFO: {} condor.sub files to submit".format(len(dirList)))
    with multiprocessing.Pool(6) as pool:
        for d in dirList:
            pool.apply_async(ResubmitJob,[d])
        pool.close()
        pool.join()
