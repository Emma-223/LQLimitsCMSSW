#!/usr/bin/env python3

from ROOT import gROOT, gStyle, gPad, TCanvas, TH2F, TGraph, kBlue, kCyan, kGreen, kYellow, kBlack, kRed, TPolyLine, TLegend, TLatex
import CMS_lumi
from tdrstyle import setTDRStyle
import math
import numpy
import os


def ComboPlot(dirName=".", beta_vals=[], m_expected_lljj=[], m_observed_lljj=[]):
    doObserved = True if len(m_observed_lljj) else False
    numbetas = len(beta_vals)
    beta_valsdub = []
    ind = numbetas
    for ii in range(0, numbetas*2):
        if ii < numbetas:
            beta_valsdub.append(beta_vals[ii])
        else:
            ind -= 1
            beta_valsdub.append(beta_vals[ind])
    
    for ii in range(0, numbetas):
        # if beta_vals[ii] == 0.5:
        #     print("## ---------------- B = 1/2 Limits ----------------------")
        #     print("##   lvjj Expected :", m_expected_lvjj[ii])
        #     print("##   lvjj Observed :", m_observed_lvjj[ii])
        #     print("##  combo Expected :", m_expected_combo[ii])
        #     print("##  combo Observed :", m_observed_combo[ii])
        #     print("##   ")
        if beta_vals[ii] > 0.999:
            print("## ----------------   B = 1 Limits ----------------------")
            print("##   lljj Expected :", m_expected_lljj[ii])
            if doObserved:
                print("##   lljj Observed :", m_observed_lljj[ii])
            # print("##  combo Expected :", m_expected_combo[ii])
            # print("##  combo Observed :", m_observed_combo[ii])
            print("##   ")
    
    
    # filename for the final plot (NB: changing the name extension changes the file format)
    fileNamePdf = dirName+"/LQ1BetaScan.pdf"
    fileNamePng = dirName+"/LQ1BetaScan.png"
    
    # axes labels for the final plot
    title = ";M_{LQ} [GeV];#beta"
    
    # integrated luminosity
    lint = "35.9 fb^{-1}"
    sqrts = "#sqrt{s} = 13 TeV"
     
    n_exclude = 204
    
    beta_CMSATLAS = [0.018,0.018,0.02,0.022,0.024,0.026,0.028,0.03,0.032,0.034,0.036,0.038,0.04,0.042,0.044,0.046,0.048,0.05,0.052,0.054,0.056,0.058,0.06,0.062,0.064,0.066,0.068,0.07,0.072,0.074,0.076,0.078,0.08,0.082,0.084,0.086,0.088,0.09,0.092,0.094,0.096,0.098,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.8,0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.9,0.902,0.904,0.906,0.908,0.91,0.912,0.914,0.916,0.918,0.92,0.922,0.924,0.926,0.928,0.93,0.932,0.934,0.936,0.938,0.94,0.942,0.944,0.946,0.948,0.95,0.952,0.954,0.956,0.958,0.96,0.962,0.964,0.966,0.968,0.97,0.972,0.974,0.976,0.978,0.98,0.982,0.984,0.986,0.988,0.99,0.992,0.994,0.996,0.998,0.9995,0.9995, 0.9459430964276001, 0.9113315376333271, 0.880788245243444, 0.8298887187568689, 0.7789832314403973, 0.7321489909432641, 0.6934572440848077, 0.6506942504070129, 0.552954522593408, 0.5977471788509756, 0.5020579765217812, 0.4572384965296794, 0.4246416982404375, 0.39408648419076153, 0.3676025169604238, 0.3044564654514025, 0.23723022587819786, 0.34112153014503444, 0.26778842034282224, 0.18626214984781253, 0.20463342758895597, 0.16789087210666898, 0.15155521777519476, 0.13521956344372055, 0.11888390911224611, 0.10254527436582372, 0.09028086685368775, 0.07185892205842359, 0.08412332957055957, 0.06366576136562607, 0.05140135385349032]
    mass_CMSATLAS = [300.,317.5,342.2,358.0,369.5,379.8,389.5,398.6,407.1,414.6,422.1,429.1,435.6,441.5,447.3,452.8,457.1,461.1,465.1,468.7,472.3,475.7,479.3,482.9,485.9,488.7,491.5,494.1,496.6,499.9,501.4,502.8,504.1,505.4,506.8,508.0,509.6,510.8,511.9,513.0,514.0,515.6,516.6,521.5,525.9,530.2,533.7,537.3,540.3,543.1,545.7,548.1,552.0,561.4,569.0,576.5,583.9,596.0,600.1,603.2,605.6,607.2,609.4,611.8,613.3,615.7,617.6,619.6,621.8,623.1,625.2,626.6,628.6,630.0,632.2,633.5,635.5,636.8,638.8,640.4,642.5,643.8,645.1,646.9,648.7,651.9,655.3,661.7,666.6,671.8,678.4,685.3,693.6,704.1,716.6,729.2,746.9,751.7,753.7,756.1,758.6,760.7,763.3,765.5,767.7,770.0,772.3,774.7,776.9,779.3,781.7,783.9,786.5,788.6,790.9,793.4,795.5,797.7,800.0,802.5,805.2,807.2,807.8,808.6,808.9,809.5,809.9,810.5,810.4,811.3,811.7,812.2,812.9,813.2,813.9,814.4,814.6,815.3,815.8,816.1,816.6,817.2,817.7,818.1,818.4,818.9,819.6,820.2,820.4,821.2,821.7,822.0,822.3,823.0,823.5,823.7,824.2,824.7,825.0,825.8,826.2,826.8,827.1,827.8,828.3,828.6,829.0,829.9,830.3,830.8,831.1,831.3, 659.1508052708638, 654.4655929721814, 652.1229868228403, 648.6090775988287, 645.0951683748169, 639.2386530014641, 633.3821376281112, 627.5256222547583, 621.6691068814055, 609.9560761346999, 613.4699853587115, 607.6134699853586, 593.5578330893118, 583.0161054172767, 574.8169838945827, 566.6178623718887, 550.2196193265007, 530.3074670571009, 559.5900439238652, 539.6778916544655, 499.8535871156661, 519.7657393850658, 479.94143484626636, 460.02928257686676, 440.11713030746705, 420.2049780380672, 399.121522693997, 379.2093704245973, 339.38506588579787, 359.2972181551976, 319.47291361639816, 299.56076134699845]
     
    mass_CMS = [300.0000,319.8520,331.7728,339.8246,351.0135,359.7973,370.1496,379.8745,389.8086,399.8472,409.4675,420.0290,429.8584,439.4788,449.3082,460.0788,469.4900,478.3784,486.8485,494.2728,498.7693,501.8018,508.4942,515.2912,522.0882,527.5257,532.4405,536.7278,541.8517,545.7207,549.0669,551.5766,557.4324,564.9614,573.7452,580.5421,586.1889,592.5676,597.2732,600.0965,602.0833,603.2336,605.7432,608.1483,610.8671,613.0631,615.0499,617.0367,618.7098,621.1149,623.7291,625.9250,627.4936,630.2124,633.2449,635.5454,638.2642,640.7738,642.9698,645.8977,647.6754,649.4530,650.0804,653.4266,656.4591,658.1322,660.8510,665.2429,669.6348,673.9221,677.7912,682.0785,685.5293,689.3983,694.5222,699.1232,704.2471,709.3710,713.7629,719.3050,723.5923,728.2979,730.9122,733.5264,736.8726,741.8919,746.8066,750.7802,755.5904,759.5640,761.1326,764.5833,768.8707,772.1123,773.3671,777.8636,782.9875,786.4382,788.9479,793.4443,798.4636,801.7053,806.2017,810.2799,814.7764,817.1815,823.3510,827.6384,832.6577,836.2130,838.1998,840.9186,841.9643,842.2780,300.0000,300.0000]
    beta_CMS = [0.0333,0.0337,0.0345,0.0350,0.0357,0.0376,0.0411,0.0445,0.0478,0.0516,0.0550,0.0593,0.0636,0.0681,0.0724,0.0767,0.0804,0.0841,0.0881,0.0914,0.0934,0.0963,0.1043,0.1135,0.1225,0.1303,0.1382,0.1453,0.1541,0.1613,0.1678,0.1710,0.1756,0.1824,0.1903,0.1966,0.2027,0.2098,0.2154,0.2197,0.2271,0.2328,0.2440,0.2547,0.2679,0.2781,0.2896,0.2992,0.3090,0.3218,0.3360,0.3485,0.3591,0.3764,0.3936,0.4092,0.4252,0.4411,0.4563,0.4722,0.4854,0.4974,0.5005,0.5085,0.5167,0.5215,0.5281,0.5378,0.5489,0.5597,0.5697,0.5786,0.5876,0.5956,0.6076,0.6185,0.6279,0.6383,0.6478,0.6586,0.6668,0.6757,0.6811,0.6874,0.6929,0.7014,0.7093,0.7205,0.7330,0.7441,0.7503,0.7589,0.7727,0.7815,0.7886,0.7994,0.8154,0.8281,0.8350,0.8470,0.8626,0.8724,0.8855,0.8974,0.9094,0.9200,0.9356,0.9498,0.9632,0.9746,0.9823,0.9897,0.9961,1.0000,1.0000,0.0333]
    
    gROOT.SetBatch(True)
    setTDRStyle()
    gStyle.SetPadLeftMargin(0.14)
    gROOT.ForceStyle()
    
    bg = TH2F("bg",title, 20, 300., 2000., 100, 0., 1.)
    bg.SetStats(False)
    bg.GetXaxis().CenterTitle()
    bg.GetYaxis().CenterTitle()
    bg.SetTitleOffset(1.,"X")
    bg.SetTitleOffset(1.05,"Y")
    #  bg.GetXaxis().SetNdivisions(505)
    bg.GetXaxis().SetLabelSize(0.03)
    bg.GetYaxis().SetLabelSize(0.03)
    
    c = TCanvas("c","",800,800)
    c.cd()
    
    bg.Draw()
    
    gr_excl_CMS = TPolyLine(116, numpy.array(mass_CMS, dtype="f"), numpy.array(beta_CMS, dtype="f"), "")
    # gr_excl_CMS.SetLineWidth(3)
    gr_excl_CMS.SetLineColor(14)
    gr_excl_CMS.SetFillColor(14)
    gr_excl_CMS.SetFillStyle(3344)
    
    # gr_excl_CMSborder = TPolyLine(n_exclude, numpy.array(mass_CMSATLAS, dtype="f"), numpy.array(beta_CMSATLAS, dtype="f"), "")
    # # gr_excl_CMS.SetLineWidth(3)
    # gr_excl_CMSborder.SetLineColor(14)
    # gr_excl_CMSborder.SetFillColor(14)
    # gr_excl_CMSborder.SetFillStyle(3344)
    
    #  ------------------------------------ LLJJ CURVES  -------------------------------------
    if doObserved:
        beta_vs_m_lljj_observed = TPolyLine(numbetas, numpy.array(m_observed_lljj, dtype="f"), numpy.array(beta_vals, dtype="f"))
    beta_vs_m_lljj_expected = TPolyLine(numbetas, numpy.array(m_expected_lljj, dtype="f"), numpy.array(beta_vals, dtype="f"))
    # beta_vs_m_lljj_observed_shade = TPolyLine(numbetas*2, m_2sigma_lljj, beta_valsdub)
    if doObserved:
        beta_vs_m_lljj_observed.SetLineWidth(2)
        beta_vs_m_lljj_observed.SetLineColor(kBlue)
        beta_vs_m_lljj_observed.SetLineStyle(1)
    beta_vs_m_lljj_expected.SetLineWidth(2)
    beta_vs_m_lljj_expected.SetLineColor(kBlue)
    beta_vs_m_lljj_expected.SetLineStyle(2)
    #beta_vs_m_lljj_observed_shade.SetLineWidth(3)
    #beta_vs_m_lljj_observed_shade.SetLineColor(kGreen+2)
    #beta_vs_m_lljj_observed_shade.SetLineStyle(4)
    #beta_vs_m_lljj_observed_shade.SetFillStyle(3425)
    #beta_vs_m_lljj_observed_shade.SetFillColor(kGreen+2)
    #beta_vs_m_lljj_observed_shade.Draw("f")
    
    # ------------------------------------ LVJJ CURVES  -------------------------------------
    # beta_vs_m_lvjj_observed = TPolyLine(numbetas, numpy.array(m_observed_lvjj, dtype="f"), numpy.array(beta_vals, dtype="f"))
    # beta_vs_m_lvjj_expected = TPolyLine(numbetas, numpy.array(m_expected_lvjj, dtype="f"), numpy.array(beta_vals, dtype="f"))
    # # beta_vs_m_lvjj_observed_shade = TPolyLine(numbetas*2, m_2sigma_lvjj, beta_valsdub)
    # beta_vs_m_lvjj_observed.SetLineWidth(2)
    # beta_vs_m_lvjj_observed.SetLineColor(kRed)
    # beta_vs_m_lvjj_observed.SetLineStyle(1)
    # beta_vs_m_lvjj_expected.SetLineWidth(2)
    # beta_vs_m_lvjj_expected.SetLineColor(kRed)
    # beta_vs_m_lvjj_expected.SetLineStyle(2)
    # #beta_vs_m_lvjj_observed_shade.SetFillStyle(3452)
    # #beta_vs_m_lvjj_observed_shade.SetFillColor(kRed)
    # #beta_vs_m_lvjj_observed_shade.Draw("f")
    
    # ------------------------------------ COMBO CURVES  -------------------------------------
    # beta_vs_m_comb_observed = TPolyLine(numbetas, numpy.array(m_observed_combo, dtype="f"), numpy.array(beta_vals, dtype="f"))
    # beta_vs_m_comb_expected = TPolyLine(numbetas, numpy.array(m_expected_combo, dtype="f"), numpy.array(beta_vals, dtype="f"))
    # beta_vs_m_comb_expected_shade2 = TGraph(numbetas*2, numpy.array(m_2sigma_combo, dtype="f"), numpy.array(beta_valsdub, dtype="f"))
    # beta_vs_m_comb_expected_shade1 = TGraph(numbetas*2, numpy.array(m_1sigma_combo, dtype="f"), numpy.array(beta_valsdub, dtype="f"))
    # beta_vs_m_comb_observed.SetLineWidth(4)
    # beta_vs_m_comb_observed.SetLineColor(kBlack)
    # beta_vs_m_comb_expected.SetLineStyle(1)
    # # beta_vs_m_comb_expected_shade.SetLineWidth(1)
    # beta_vs_m_comb_expected_shade1.SetFillStyle(1001)
    # beta_vs_m_comb_expected_shade1.SetFillColor(kGreen)
    # beta_vs_m_comb_expected_shade2.SetFillStyle(1001)
    # beta_vs_m_comb_expected_shade2.SetFillColor(kYellow)
    # beta_vs_m_comb_expected.SetLineWidth(4)
    # beta_vs_m_comb_expected.SetLineColor(kBlack)
    # beta_vs_m_comb_expected.SetLineStyle(2)
    
    # ------------------------------------ DRAW  -------------------------------------
    # beta_vs_m_comb_expected_shade2.Draw("f")
    # beta_vs_m_comb_expected_shade1.Draw("f")
    
    # gr_excl_CMSborder.Draw("L")
    # gr_excl_CMSborder.Draw("F")
    
    beta_vs_m_lljj_expected.Draw("C")
    # beta_vs_m_lvjj_expected.Draw("C")
    if doObserved:
        beta_vs_m_lljj_observed.Draw("C")
    # beta_vs_m_lvjj_observed.Draw("C")
    # beta_vs_m_comb_expected.Draw("C")
    # beta_vs_m_comb_observed.Draw("C")
    
    gPad.RedrawAxis()
    
    legend = TLegend(.63,.2,.93,.44)
    legend.SetBorderSize(1)
    # legend.SetFillColor(0)
    # legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetMargin(0.15)
    # legend.SetHeader("95% CL limits")
    legend.SetHeader("95% CL limits, scalar LQ#bar{LQ}")
    # legend.AddEntry(beta_vs_m_comb_observed,"CMS eejj + e#nujj (Obs.)","l")
    # legend.AddEntry(beta_vs_m_comb_expected,"CMS eejj + e#nujj (Exp.)","l")
    # 
    # legend.AddEntry(beta_vs_m_lvjj_observed,"CMS e#nujj (Obs.)","l")
    # legend.AddEntry(beta_vs_m_lvjj_expected,"CMS e#nujj (Exp.)","l")
    
    if doObserved:
        legend.AddEntry(beta_vs_m_lljj_observed,"CMS eejj (Obs.)","l")
    legend.AddEntry(beta_vs_m_lljj_expected,"CMS eejj (Exp.)","l")
    
    # legend.AddEntry(gr_excl_CMSborder,"CMS, 7 TeV, 5.0 fb^{-1}","f")
    
    # legend.AddEntry(gr_excl_D0,"D#oslash (Obs.), 1 fb^{-1}","l")
    legend.Draw()
    
    # draw the lumi text on the canvas
    CMS_lumi.lumi_13TeV = lint
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumiTextSize = 0.7
    CMS_lumi.relPosX = 0.1  # control position of extraText
    CMS_lumi.hOffset = 0.0
    iPos = 0
    CMS_lumi.CMS_lumi(c, 4, iPos)
    
    # c.SetGridx()
    # c.SetGridy()
    c.RedrawAxis()
    
    l3 = TLatex()
    l3.SetTextAlign(12)
    l3.SetTextFont(42)
    l3.SetNDC()
    l3.SetTextSize(0.06)
    l3.SetTextColor(kBlue)
    l3.SetTextAngle(30)
    l3.DrawLatex(0.24,0.35,"eejj")
    # l3.SetTextColor(kRed)
    # l3.SetTextAngle(-20)
    # l3.DrawLatex(0.25,0.84,"e#nujj")
    # l3.SetTextColor(kBlack)
    # l3.SetTextAngle(79)
    # l3.DrawLatex(0.655,0.50,"eejj + e#nujj")
    
    # c.SetGridx()
    # c.SetGridy()
    c.SetLogy(0)
    c.SaveAs(fileNamePdf)
    c.SaveAs(fileNamePng)


if __name__ == "__main__":
    ComboPlot()
