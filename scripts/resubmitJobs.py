import sys
import os
import subprocess
import glob

jobIds = []

txtfile = "heldJobs.txt"
#print("get held jobs")
#subprocess.run(["condor_q","-hold",">",txtfile])

with open(txtfile) as f:
    lines = f.readlines()

for l in lines:
    if "CPU usage" in l:
        jobIds.append(l.split()[0])

print("Ids to resubmit: ",jobIds)
condorFile = "condor.{}.log"

for jobId in jobIds:
    string = "*/*/*/*/condor.{}.log".format(jobId)
    path = glob.glob(string)
    path = path[0]
    directory = path.replace("/condor.{}.log".format(jobId),"")
    outFiles = glob.glob(directory+"/condor.*.out")
    errFiles = glob.glob(directory+"/condor.*.err")
    filesToRemove = outFiles+errFiles
    if len(filesToRemove):
        print("remove {} .out and {} .err files from ".format(len(outFiles),len(errFiles))+directory)
        cmd = ["rm"] + filesToRemove
        subprocess.run(cmd)
    subfile = directory+"/condor.sub"
    nSubFiles = len(glob.glob(subfile))
    if not nSubFiles==1:
        print("ERROR: found {} condor.sub files in ".format(nSubFiles)+directory)
        exit()
    print("submit file "+subfile)
    subprocess.run(["condor_submit",subfile])
    print("remove held job {}".format(jobId))
    subprocess.run(["condor_rm",str(jobId)])
