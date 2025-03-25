#!/usr/bin.env python3 
import os
import sys
import subprocess
import ROOT
from ROOT import TFile, TCanvas, RDataFrame, TKey

if __name__ == "__main__":
    txtFile = "limitDirs_6mar25.txt"
    limitCalcName = "limits_6mar25"
    with open(txtFile,'r') as f:
        eosDirList = f.readlines()

    limits = {}
    for m in range(300, 3100, 100):
        limits[str(m)] = {}
    filenameList = []
    cBlank = TCanvas()

    for d in eosDirList:
        if not "Limits" in d:
            continue
        d = d.strip()
        mass = d.split(".")[1].replace("M","")
        if int(mass) >= 2000: # or int(mass) < 1500:
            continue
        quantile = d.split(".")[2].replace("p",".")
        #if not "025" in quantile:
         #   continue
        sigSF = d.split("signalScaleFactor")[1]
        filename = "higgsCombine.signalScaleFactor{}.HybridNew.mH{}.quant{}.gridAll.root".format(sigSF,mass,quantile)
        filename = "root://eoscms/"+d+"/"+filename

        if "0.16" in quantile:
            quantLong = quantile.replace("0.16","0.160")
        elif "0.5" in quantile:
            quantLong = quantile.replace("0.5","0.500")
        elif "0.84" in quantile:
            quantLong = quantile.replace("0.84","0.840")
        else:
            quantLong = quantile
        #tfileName = "root://eoscms/"+d+"/"+"higgsCombine.signalScaleFactor{}.HybridNew.mH{}.quant{}.root".format(sigSF,mass,quantLong)
        tfileName = d+"/"+"higgsCombine.signalScaleFactor{}.HybridNew.mH{}.quant{}.root".format(sigSF,mass,quantLong)
        tfile = TFile.Open(tfileName)
        limitTree = tfile.Get("limit")
        df = RDataFrame(limitTree)
        limit = df.Sum("limit").GetValue() #Not the correct way to do this, but faster than figuring how how to do it correctly in pyroot and works since there's only one entry
        if limit < 1:
            limits[str(mass)][quantile] = round(limit,6)
        else:
            limits[str(mass)][quantile] = round(limit,5)


        commandArgs = ["python3","/afs/cern.ch/user/e/eipearso/public/leptoquark_analysis/CMSSW_14_1_0_pre4/src/HiggsAnalysis/CombinedLimit/test/plotTestStatCLs.py","--input",filename,"--poi","r","--val","all","--mass",mass,"--expected","-q",quantile]
        subprocess.run(commandArgs)

        canvasesFile = "test_stat_distributions.root"
        tfile = TFile.Open(canvasesFile)
        names = tfile.GetListOfKeys()
        smallestDiff = 999
        for iname in names:
            canvas = iname.ReadObj()
            cname = canvas.GetName()
            cLimit = float(cname.split("_")[2])
            diff = abs(limit - cLimit)
            if diff < smallestDiff:
                print("{} - {} = {}".format(limit, cLimit, abs(limit - cLimit)))
                smallestDiff = diff
                cToGet = canvas
        print("limit: {}".format(limit))
        print(cToGet.GetName())
        cToGet.Print("qPlot_{}_{}.pdf".format(mass,quantile))
        cernboxDir = os.getenv("LQDATAEOS")+"/limits/{}/testStatistics".format(limitCalcName)
        if not os.path.isdir(cernboxDir):
            os.mkdir(cernboxDir)
        if not os.path.isdir(cernboxDir+"/"+mass):
            os.mkdir(cernboxDir+"/"+mass)
        #cernboxDir = "root://eosuser.cern.ch/"+cernboxDir

        subprocess.run(["cp",canvasesFile, cernboxDir+"/"+mass+"/"+canvasesFile.replace(".root","{}.root".format(quantile))])
        subprocess.run(["cp","qPlot_{}_{}.pdf".format(mass,quantile),cernboxDir+"/"+mass])

        os.remove(canvasesFile)
        os.remove("qPlot_{}_{}.pdf".format(mass,quantile))
