#!/usr/bin/env python3

from ROOT import gROOT, gStyle, gPad, TCanvas, TH2F, TGraph, kBlue, kCyan, kGreen, kYellow, kBlack, TPolyLine, TLegend, TColor
import math
import numpy
import os
import cmsstyle as CMS


def ReadXSecFile(filename):
    masses = []
    xsTh = []
    yPDF_up = []
    yPDF_down = []
    with open(os.path.expandvars(filename), "r") as xsecFile:
        for line in xsecFile:
            line = line.strip()
            if line.startswith("#"):
                continue
            split = line.split()
            if len(split)==7:
                #print("length of this line is not 7; don't know how to handle it. Quitting.  Line looks like '"+line+"'")
                #exit(-1)
                masses.append(float(split[0]))
                xs = float(split[1])
                xsTh.append(xs)
                if "scalar" in filename:
                    yPDF_up.append(xs*(1+float(split[5])/100.))
                    yPDF_down.append(xs*(1-float(split[6])/100.))
                elif "vector" in filename:
                    yPDF_up.append(float(split[6]))
                    yPDF_down.append(float(split[5]))
                else:
                    print("Do not know how to assign yPDF_up/down for this xsec file. Quitting.")
                    exit(-1)
            elif len(split)==3: #stopPair is like this
                masses.append(float(split[0]))
                xs = float(split[1])
                xsTh.append(xs)
                yPDF_up.append(xs*(1+float(split[2])/100.))
                yPDF_down.append(xs*(1+float(split[2])/100.))
            elif len(split)==2: #ATLAS xsec file
                masses.append(float(split[0]))
                xs = float(split[1])
                xsTh.append(xs)
                yPDF_up.append(0)
                yPDF_down.append(0)
            else:
                print("length of this line is not 7 or 3 or 2; don't know how to handle it. Quitting.  Line looks like '"+line+"'")
                exit(-1)
            #yPDF_up.append(float(split[6]))
            #yPDF_down.append(float(split[5]))
    return masses, xsTh, yPDF_up, yPDF_down


def BR_Sigma_EE_vsMass(dirName=".", intLumi="35.9", mData = [], x_shademasses = [], xsUp_expected = [], xsUp_observed = [], y_1sigma = [], y_2sigma = [], xsThFilename="xsection_theory_13TeV_scalarPairLQ.txt"):
    #xsThFilename = "$LQANA/config/xsection_theory_13TeV_scalarPairLQ.txt"
    #xsThFilename = "$LQANA/config/xsection_theory_13TeV_stopPair.txt"
    #xsThFilename = "xsection_theory_13TeV_vectorSingletU1_YM.txt"
    #xsThFilename = "xsection_theory_13TeV_atlas.txt"
    mTh, xsTh, y_pdf, yPDFDown = ReadXSecFile(xsThFilename)
    print(xsTh)
    nTH = len(mTh)
    massPoints = len(mData)
    x_pdf = mTh
    x_pdf.extend(list(reversed(mTh)))
    y_pdf.extend(list(reversed(yPDFDown)))

    # filename for the final plot (NB: changing the name extension changes the file format)
    fileName2 = dirName+"/BR_Sigma_EE.pdf"
    fileName3 = dirName+"/BR_Sigma_EE.png"
    fileName1 = dirName+"/BR_Sigma_EE.eps"

    # integrated luminosity
    CMS.SetLumi(intLumi)
    CMS.SetExtraText("Preliminary")
    CMS.ResetAdditionalInfo()
    # CMS.cmsGrid(True)

    #mData = [1400.0 , 1500.0 , 1600.0 , 1700.0 , ]
    #x_shademasses = [1400.0 , 1500.0 , 1600.0 , 1700.0 , 1700.0 , 1600.0 , 1500.0 , 1400.0 , ]
    #xsUp_expected = [0.000171843 , 0.00015757172 , 0.00014484 , 0.000138962523 , ]
    #xsUp_observed = []
    #y_1sigma = [0.00011263875 , 0.00010142096 , 9.211185e-05 , 8.8379716e-05 , 0.000230907651 , 0.00024067935 , 0.000256793 , 0.00027595575 , ]
    #y_2sigma = [8.0556e-05 , 7.201988e-05 , 6.44964e-05 , 6.1881149e-05 , 0.000373152242 , 0.00038581755 , 0.00040444652 , 0.00042817425 , ]

    doObserved = False
    if len(xsUp_observed) > 0:
        doObserved = True

    gROOT.SetBatch(True)

    plotLow = 0.00001
    plotHigh = 30.

    c = CMS.cmsCanvas('c', 300, 3000, plotLow, plotHigh, 'M_{LQ} (GeV)', '#sigma #times #beta^{2} (pb)', True, 0, extraSpace=0.025)
    # c.SetBottomMargin(0.13)
    # c.SetLeftMargin(0.14)
    c.SetRightMargin(0.06)
    c.SetLogy()
    # c.cd()
    # c.SetGridx()
    # c.SetGridy()
    # CMS.cmsStyle.SetPadGridX(True)
    # CMS.cmsStyle.SetPadGridY(True)
    # CMS.UpdatePad()
    gPad.SetGridx()
    gPad.SetGridy()
    gPad.RedrawAxis("g")

    CMS.GetcmsCanvasHist(c).SetLabelSize(0.04, "X")
    CMS.GetcmsCanvasHist(c).SetLabelSize(0.04, "Y")
    # CMS.UpdatePad()
    # c.RedrawAxis()
    # c.GetFrame().Draw()
    # bg = TH2F("bg", title, 800, 200., 1800., 500, plotLow, plotHigh)
    # # bg = TH2F("bg", title, 900, 0., 1800., 500, plotLow, plotHigh)
    # bg.GetXaxis().CenterTitle()
    # bg.GetYaxis().CenterTitle()
    # bg.SetStats(False)
    # bg.SetTitleOffset(1., "X")
    # bg.SetTitleOffset(1.25, "Y")
    # bg.SetTitleSize(0.05, "X")
    # bg.SetTitleSize(0.05, "Y")
    # bg.SetLabelSize(0.04, "X")
    # bg.SetLabelSize(0.04, "Y")
    # bg.Draw()

    xsTh_vs_m = TGraph(nTH, numpy.array(mTh, dtype="f"), numpy.array(xsTh, dtype="f"))
    # xsTh_vs_m.SetLineWidth(2)
    # xsTh_vs_m.SetLineColor(kBlue)
    # xsTh_vs_m.SetFillColor(kCyan-6)
    # xsTh_vs_m.SetMarkerSize(0.00001)
    # xsTh_vs_m.SetMarkerStyle(22)
    # xsTh_vs_m.SetMarkerColor(kBlue)

    xsData_vs_m_expected = TGraph(massPoints, numpy.array(mData, dtype="f"), numpy.array(xsUp_expected, dtype="f"))
    # xsData_vs_m_expected.SetMarkerStyle(0)
    # xsData_vs_m_expected.SetMarkerColor(kBlack)
    # xsData_vs_m_expected.SetLineColor(kBlack)
    # xsData_vs_m_expected.SetLineWidth(2)
    # xsData_vs_m_expected.SetLineStyle(7)
    # xsData_vs_m_expected.SetMarkerSize(0.001)

    xsData_vs_m_observed = TGraph()
    if doObserved:
        xsData_vs_m_observed = TGraph(massPoints, mData, xsUp_observed)
        xsData_vs_m_observed.SetMarkerStyle(21)
        xsData_vs_m_observed.SetMarkerColor(kBlack)
        xsData_vs_m_observed.SetLineColor(kBlack)
        xsData_vs_m_observed.SetLineWidth(2)
        xsData_vs_m_observed.SetLineStyle(1)
        xsData_vs_m_observed.SetMarkerSize(1)

    xsUp_observed_logY = []
    xsUp_expected_logY = []
    xsTh_logY = []
    if doObserved:
        for ii in range(0, massPoints):
            xsUp_observed_logY.append(math.log(xsUp_observed[ii]))
        xsData_vs_m_observed_log = TGraph(massPoints, numpy.array(mData, dtype="f"), numpy.array(xsUp_observed_logY, dtype="f"))
    for ii in range(0, massPoints):
        xsUp_expected_logY.append(math.log(xsUp_expected[ii]))
    for ii in range(0, nTH):
        xsTh_logY.append(math.log(xsTh[ii]))
    xsTh_vs_m_log = TGraph(nTH, numpy.array(mTh, dtype="f"), numpy.array(xsTh_logY, dtype="f"))
    xsData_vs_m_expected_log = TGraph(massPoints, numpy.array(mData, dtype="f"), numpy.array(xsUp_expected_logY, dtype="f"))

    # xsTh_vs_m_tgraph = TGraph(221, numpy.array(mTh, dtype="f"), numpy.array(xsTh, dtype="f"))

    obslim = 0.0
    exlim = 0.0
    for mtest in numpy.linspace(1300.0, 3000.0, 17000, endpoint=False):
        if pow(10.0, xsData_vs_m_expected_log.Eval(mtest))/pow(10.0, xsTh_vs_m_log.Eval(mtest)) < 1.0 and pow(10.0, xsData_vs_m_expected_log.Eval(mtest+0.1))/pow(10.0, xsTh_vs_m_log.Eval(mtest+0.10)) > 1.0:
            exlim = mtest
        if doObserved:
            if pow(10.0, xsData_vs_m_observed_log.Eval(mtest))/pow(10.0, xsTh_vs_m_log.Eval(mtest)) < 1.0 and pow(10.0, xsData_vs_m_observed_log.Eval(mtest+0.1))/pow(10.0, xsTh_vs_m_log.Eval(mtest+0.10)) > 1.0:
                obslim = mtest
    print("## LLJJ expected limit:", exlim, "GeV")
    if doObserved:
        print("## LLJJ observed limit:", obslim, "GeV")


   #  for 1-sigma and 2-sigma expected for long-lived 2D
   # exlim = 0.0
   # for (Double_t mtest=mData[0]+.10 mtest<mData[massPoints-1]-.10 mtest = mtest+0.10){
   #   if(( xsy_2sigma_2->Eval(mtest)/xsTh_vs_m_tgraph->Eval(mtest) ) < 1.0 && ( xsy_2sigma_2->Eval(mtest+0.1)/xsTh_vs_m_tgraph->Eval(mtest+0.10) ) > 1.0) exlim = mtest 
   #  }
   #  std::cout<<"## LLJJ expected limit -2 sigma: "<<exlim<<" GeV"<<std::endl

   # exlim = 0.0
   # for (Double_t mtest=mData[0]+.10 mtest<mData[massPoints-1]-.10 mtest = mtest+0.10){
   #   if(( xsy_1sigma_2->Eval(mtest)/xsTh_vs_m_tgraph->Eval(mtest) ) < 1.0 && ( xsy_1sigma_2->Eval(mtest+0.1)/xsTh_vs_m_tgraph->Eval(mtest+0.10) ) > 1.0) exlim = mtest 
   #  }
   #  std::cout<<"## LLJJ expected limit -1 sigma: "<<exlim<<" GeV"<<std::endl

   # exlim = 0.0
   # for (Double_t mtest=mData[0]+.10 mtest<mData[massPoints-1]-.10 mtest = mtest+0.10){
   #   if(( xsy_1sigma_1->Eval(mtest)/xsTh_vs_m_tgraph->Eval(mtest) ) < 1.0 && ( xsy_1sigma_1->Eval(mtest+0.1)/xsTh_vs_m_tgraph->Eval(mtest+0.10) ) > 1.0) exlim = mtest 
   #  }
   #  std::cout<<"## LLJJ expected limit +1 sigma: "<<exlim<<" GeV"<<std::endl

   # exlim = 0.0
   # for (Double_t mtest=mData[0]+.10 mtest<mData[massPoints-1]-.10 mtest = mtest+0.10){
   #   if(( xsy_2sigma_1->Eval(mtest)/xsTh_vs_m_tgraph->Eval(mtest) ) < 1.0 && ( xsy_2sigma_1->Eval(mtest+0.1)/xsTh_vs_m_tgraph->Eval(mtest+0.10) ) > 1.0) exlim = mtest 
   #  }
   #  std::cout<<"## LLJJ expected limit +2 sigma: "<<exlim<<" GeV"<<std::endl

    # #  region excluded by Tevatron limits
    # # x_shaded[5] = {1000,1080,1080,1000,1000}// CHANGED FOR LQ2
    # x_shaded = [200, 1080, 1080, 200, 200]  # CHANGED FOR LQ2
    # y_shaded = [plotLow, plotLow, plotHigh, plotHigh, plotLow]  # CHANGED FOR LQ2

    # x_shaded2 = [200, 1000, 1000, 200, 200]  # CHANGED FOR LQ2
    # y_shaded2 = [plotLow, plotLow, plotHigh, plotHigh, plotLow]  # CHANGED FOR LQ2

    # # x_shaded3 = [840, obslim, obslim, 840, 840]  # CHANGED FOR LQ2
    # x_shaded3 = [1080, obslim, obslim, 1080, 1080]  # CHANGED FOR LQ2
    # y_shaded3 = [plotLow, plotLow, plotHigh, plotHigh, plotLow]  # CHANGED FOR LQ2

    # p2 = TPolyLine(5, numpy.array(x_shaded2, dtype="f"), numpy.array(y_shaded2, dtype="f"), "")
    # #  pl.SetFillStyle(3001)
    # p2.SetFillColor(8)
    # p2.SetFillStyle(3345)
    # p2.SetLineColor(8)   # CHANGED FOR LQ2
    # # p2.Draw()
    # # p2.Draw("F")

    # pl = TPolyLine(5, numpy.array(x_shaded, dtype="f"), numpy.array(y_shaded, dtype="f"), "")
    # CMS.ResetAdditionalInfo()
    # #  pl.SetFillStyle(3001)
    # pl.SetLineColor(14)
    # pl.SetFillColor(14)
    # # pl.SetFillStyle(3344)
    # pl.SetFillStyle(3354)
    # # pl.Draw()
    # # pl.Draw("F")

    # p3 = TPolyLine(5, numpy.array(x_shaded3, dtype="f"), numpy.array(y_shaded3, dtype="f"), "")
    # p3.SetLineColor(46)
    # p3.SetFillColor(46)
    # # p3.SetFillStyle(3354)
    # p3.SetFillStyle(3345)
    # # p3.Draw()
    # # p3.Draw("F")

    exshade1 = TGraph(2*massPoints, numpy.array(x_shademasses, dtype="f"), numpy.array(y_1sigma, dtype="f"))
    #exshade1.SetFillColor(kGreen)
    # exshade1.SetFillColor(TColor.GetColor("#607641"))
    exshade2 = TGraph(2*massPoints, numpy.array(x_shademasses, dtype="f"), numpy.array(y_2sigma, dtype="f"))
    #exshade2.SetFillColor(kYellow)
    # exshade2.SetFillColor(TColor.GetColor("#F5BB54"))

    # exshade2.Draw("f")
    # exshade1.Draw("f")
    CMS.cmsDraw(exshade2, "f", fcolor=TColor.GetColor("#F5BB54"))
    CMS.cmsDraw(exshade1, "f", fcolor=TColor.GetColor("#607641"))

    # gPad.RedrawAxis()

    grshade = TGraph(2*nTH, numpy.array(x_pdf, dtype="f"), numpy.array(y_pdf, dtype="f"))
    # grshade.SetFillColor(kCyan-6)
    # grshade.SetFillStyle(3001)
    # grshade.Draw("f")
    CMS.cmsDraw(grshade, "f", fcolor=kCyan-6, fstyle=3001)

    # xsTh_vs_m.Draw("L")
    CMS.cmsDraw(xsTh_vs_m, "l", fcolor=kCyan-6, lwidth=2, lcolor=kBlue, msize=0.00001, marker=22, mcolor=kBlue)
    # xsData_vs_m_expected.Draw("LP")
    CMS.cmsDraw(xsData_vs_m_expected, "lp", marker=0, mcolor=kBlack, lcolor=kBlack, lwidth=2, lstyle=7, msize=0.001)
    if doObserved:
        xsData_vs_m_observed.Draw("LP")

    # grshade.SetFillStyle(1001)

    legXSize = 0.5
    legYSize = 0.27
    legXStart = 0.4
    legYStart = 0.62
    # legend = TLegend(legXStart, legYStart, legXStart+legXSize, legYStart+legYSize)
    legend = CMS.cmsLeg(legXStart, legYStart, legXStart+legXSize, legYStart+legYSize, textSize=0.036, textFont=42)
    legend.SetBorderSize(1)
    legend.SetFillColor(0)
    legend.SetFillStyle(1001)
    # legend.SetTextSize(.036)
    # legend.SetTextFont(42)
    legend.SetMargin(0.15)
    if "vector" in xsThFilename:
        legend.SetHeader("Vector LQ #bar{LQ} #rightarrow eejj")
    else:
        legend.SetHeader("Scalar LQ #bar{LQ} #rightarrow eejj")
    # legend.AddEntry(p2,"ATLAS exclusion (20 fb^{-1}, 8TeV)","f")
    # legend.AddEntry(pl,"CMS exclusion (19.7 fb^{-1}, 8 TeV)","f")
    # legend.AddEntry(p3,"CMS exclusion (2.7 fb^{-1}, 13 TeV)","f")

    legend.AddEntry(xsTh_vs_m, "#sigma_{theory}#times #beta^{2}  with unc. ( #beta=1)", "lf")
    # legend.AddEntry(xsTh_vs_m,"#sigma_{theory}#times#beta^{2}  with unc. (#beta=0.5)","lf")
    # legend.AddEntry(xsTh_vs_m,"#sigma_{theory}#times2#beta(1-#beta)  with unc. (#beta=0.5)","lf")
    legend.AddEntry(xsData_vs_m_expected,  "Expected 95% CL upper limit", "lp")
    if doObserved:
        legend.AddEntry(xsData_vs_m_observed,  "Observed 95% CL upper limit", "lp")
    # legend.Draw()

    # c.Modified()
    # c.Update()
    # # draw the lumi text on the canvas
    # CMS_lumi.lumi_13TeV = lint
    # CMS_lumi.writeExtraText = 1
    # CMS_lumi.extraText = "Preliminary"
    # CMS_lumi.lumiTextSize = 0.7
    # CMS_lumi.relPosX = 0.1  # control position of extraText
    # CMS_lumi.hOffset = 0.0
    # iPos = 0
    # CMS_lumi.CMS_lumi(c, 4, iPos)

    # c.SetGridx()
    # c.SetGridy()
    # c.RedrawAxis()
    # legend.Draw()
    # gStyle.SetPadGridX(True)
    # gStyle.SetPadGridY(True)

    c.SaveAs(fileName1)
    c.SaveAs(fileName2)
    c.SaveAs(fileName3)


if __name__ == "__main__":
    BR_Sigma_EE_vsMass()
