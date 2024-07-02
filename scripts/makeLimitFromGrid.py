#!/usr/bin/env python3
import os
from optparse import OptionParser
import subprocess
import shlex
import glob
from termcolor import colored
import time
import datetime
import numpy as np
import re
from pathlib import Path
from ROOT import TFile
from combineCommon import SeparateDatacards
from BR_Sigma_EE_vsMass import BR_Sigma_EE_vsMass
from runLimits import ExtractLimitResult, ReadXSecFile, CreateArraysForPlotting


quantilesExpected = ["0.025", "0.160", "0.500", "0.840", "0.975"]
xsThFilename = "$LQANA/config/xsection_theory_13TeV_scalarPairLQ.txt"
xsThByMass, yPDFUpByMass, yPDFDownByMass = ReadXSecFile(xsThFilename)
rLimitsByMassAndQuantile = {}
xsecLimitsByMassAndQuantile = {}

massList = list(range(1300, 1900, 100))
print("INFO: considering LQ masses {}".format(massList))
for mass in massList:
    rLimitsByMassAndQuantile[mass] = {}
    xsecLimitsByMassAndQuantile[mass] = {}
    # print("INFO: considering LQ mass {}".format(mass))
    for quantileExp in quantilesExpected:
        # LQ1300/higgsCombineTest.HybridNew.mH1300.quant0.
        # print("INFO: considering quantileExp={}".format(quantileExp))
        rootFileName = Path(sorted(glob.glob('LQ{}/higgsCombineTest.HybridNew.mH{}.quant{}.root'.format(mass, mass, quantileExp)), key=os.path.getmtime)[-1]).resolve()
        limit, limitErr, quantile = ExtractLimitResult(rootFileName)
        rLimitsByMassAndQuantile[mass][str(quantileExp)] = limit
        xsecLimitsByMassAndQuantile[mass][str(quantileExp)] = limit * xsThByMass[float(mass)]

masses, shadeMasses, xsMedExp, xsOneSigmaExp, xsTwoSigmaExp, xsObs = CreateArraysForPlotting(xsecLimitsByMassAndQuantile)
print("mData =", list(masses))
print("x_shademasses =", list(shadeMasses))
print("xsUp_expected =", list(xsMedExp))
print("xsUp_observed =", list(xsObs))
print("y_1sigma =", list(xsOneSigmaExp))
print("y_2sigma =", list(xsTwoSigmaExp))
print("Make plot and calculate mass limit")
BR_Sigma_EE_vsMass(os.getcwd(), masses, shadeMasses, xsMedExp, xsObs, xsOneSigmaExp, xsTwoSigmaExp)
