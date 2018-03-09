#!/usr/bin/env python

from ROOT import *
import sys
import numpy as np


systematics =  ["_JER_do","_JER_up",
                "_JES_do","_JES_up",
                "_gluonSplitting_do","_gluonSplitting_up",
                "_bFragmentation_do","_bFragmentation_up",
                "_cFragmentation_do","_cFragmentation_up",
                "_cdFragmentation_do","_cdFragmentation_up",
                "_v0_do","_v0_up",
                "_bSV_do","_bSV_up",
                "_cSV_do","_cSV_up",
                "_lSV_do","_lSV_up",
                ]


def compute(filein, dirin, match, value) :
   result={}    
   thedir = filein.Get(dirin) 
   thedir.cd()
   nextiter =TIter(thedir.GetListOfKeys())
   while (True):
      key = nextiter();
      if key == None: break;
      cl = gROOT.GetClass(key.GetClassName())
      if not cl.InheritsFrom("TH1"): continue
      h = key.ReadObj();
      if (h.Integral()>0 and h.GetName().endswith(match) and "SVmass" in h.GetName() and "NoSV" not in h.GetName() and "_SV" not in h.GetName()):
        isCentral=True
        print h.GetName()
        for ext in systematics:
          if ext in h.GetName():
            print "this is a systematic variation", ext
            isCentral=False
            break
        if isCentral:
          svEff = (h.Integral()-h.GetBinContent(1))/h.Integral()
          svEffUp = min(svEff*value, 1.)
          svEffDo = svEff/value
          aUp=np.array([[h.GetBinContent(1), h.Integral()-h.GetBinContent(1)],[-svEffUp*h.GetBinContent(1), (1-svEffUp)*(h.Integral()-h.GetBinContent(1))]])
          aDo=np.array([[h.GetBinContent(1), h.Integral()-h.GetBinContent(1)],[-svEffDo*h.GetBinContent(1), (1-svEffDo)*(h.Integral()-h.GetBinContent(1))]])
          b=np.array([h.Integral(), 0.])
          xUp=np.linalg.solve(aUp, b)
          xDo=np.linalg.solve(aDo, b)
          result[h.GetName()]={}
          result[h.GetName()]["up"] = xUp
          result[h.GetName()]["do"] = xDo

          #print svEffUp, xUp
          #print svEffDo, xDo
          #hUp = h.Clone()
          #hUp.SetNameTitle(h.GetName()+match+"_up",h.GetName()+match+"_up")
          #hDo = h.Clone()
          #hDo.SetNameTitle(h.GetName()+match+"_do",h.GetName()+match+"_do")
          ##hUp.SetBinContent(1, h.GetBinContent(1)*xUp[0])
          ##hDo.SetBinContent(1, h.GetBinContent(1)*xDo[0])
          #for ibin in range(2, h.GetNbinsX()+1):
          #  hUp.SetBinContent(ibin, h.GetBinContent(ibin)*value) 
          #  hDo.SetBinContent(ibin, h.GetBinContent(ibin)/value)
          #hUp.Scale(h.Integral()/hUp.Integral())  
          #hDo.Scale(h.Integral()/hDo.Integral())  
          ##hUp.Write()
          ##hDo.Write()
          #plots.append(hUp)
          #plots.append(hDo)
    
   #for plot in plots:
   # plot.Write()
   return result

def apply(filein, dirin, result, suffix) :
   plots=[] 
   thedir = filein.Get(dirin)
   thedir.cd()
   nextiter =TIter(thedir.GetListOfKeys())
   while (True):
      key = nextiter();
      if key == None: break;
      cl = gROOT.GetClass(key.GetClassName())
      if not cl.InheritsFrom("TH1"): continue
      h = key.ReadObj();
      if (h.Integral()>0 and "SVmass" in h.GetName() and "NoSV" not in h.GetName() and "_SV" not in h.GetName()):
        isCentral=True
        #print h.GetName()
        for ext in systematics:
          if ext in h.GetName():
            #print "this is a systematic variation", ext
            isCentral=False
            break
        if isCentral:
          for basename in result.keys():
            if basename in h.GetName():
              print "applying correction", basename, "to", h.GetName()
              hUp = h.Clone()
              hUp.SetNameTitle(h.GetName()+"_"+suffix+"_up",h.GetName()+"_"+suffix+"_up")
              hDo = h.Clone()
              hDo.SetNameTitle(h.GetName()+"_"+suffix+"_do",h.GetName()+"_"+suffix+"_do")
              hUp.SetBinContent(1, h.GetBinContent(1)*result[basename]['up'][0])
              hDo.SetBinContent(1, h.GetBinContent(1)*result[basename]['do'][0])
              for ibin in range(2, h.GetNbinsX()+1):
                hUp.SetBinContent(ibin, h.GetBinContent(ibin)*result[basename]['up'][1]) 
                hDo.SetBinContent(ibin, h.GetBinContent(ibin)*result[basename]['do'][1])   
              plots.append(hUp)
              plots.append(hDo)
   for plot in plots:
     plot.Write()
            
   
filein=TFile(sys.argv[1], "update")
result_b_20_30 = compute(filein, "mc", "20.0-30.0_b", 1.2)
result_b_30_50 = compute(filein, "mc", "30.0-50.0_b", 1.1)
result_b_50_70 = compute(filein, "mc", "50.0-70.0_b", 1.1)
result_b_70_100 = compute(filein, "mc", "70.0-100.0_b", 1.05)
result_b_100_140 = compute(filein, "mc", "100.0-140.0_b", 1.05)
result_b_140_200 = compute(filein, "mc", "140.0-200.0_b", 1.05)
result_b_200_300 = compute(filein, "mc", "200.0-300.0_b", 1.05)
result_b_300_670 = compute(filein, "mc", "300.0-670.0_b", 1.05)
result_b_670_1000 = compute(filein, "mc", "670.0-1000.0_b", 1.05)
result_b_dummy = compute(filein, "mc", "_b", 1.)
result_c_20_30 = compute(filein, "mc", "20.0-30.0_c", 1.2)
result_c_30_50 = compute(filein, "mc", "30.0-50.0_c", 1.1)
result_c_50_70 = compute(filein, "mc", "50.0-70.0_c", 1.1)
result_c_70_100 = compute(filein, "mc", "70.0-100.0_c", 1.05)
result_c_100_140 = compute(filein, "mc", "100.0-140.0_c", 1.05)
result_c_140_200 = compute(filein, "mc", "140.0-200.0_c", 1.05)
result_c_200_300 = compute(filein, "mc", "200.0-300.0_c", 1.05)
result_c_300_670 = compute(filein, "mc", "300.0-670.0_c", 1.05)
result_c_670_1000 = compute(filein, "mc", "670.0-1000.0_c", 1.05)
result_c_dummy = compute(filein, "mc", "_c", 1.)
result_l_20_30 = compute(filein, "mc", "20.0-30.0_l", 1.2)
result_l_30_50 = compute(filein, "mc", "30.0-50.0_l", 1.1)
result_l_50_70 = compute(filein, "mc", "50.0-70.0_l", 1.1)
result_l_70_100 = compute(filein, "mc", "70.0-100.0_l", 1.05)
result_l_100_140 = compute(filein, "mc", "100.0-140.0_l", 1.05)
result_l_140_200 = compute(filein, "mc", "140.0-200.0_l", 1.05)
result_l_200_300 = compute(filein, "mc", "200.0-300.0_l", 1.05)
result_l_300_670 = compute(filein, "mc", "300.0-670.0_l", 1.05)
result_l_670_1000 = compute(filein, "mc", "670.0-1000.0_l", 1.05)
result_l_dummy = compute(filein, "mc", "_l", 1.)


apply(filein, "mc", result_b_20_30, "bSV")
apply(filein, "mc", result_b_30_50, "bSV")
apply(filein, "mc", result_b_50_70, "bSV")
apply(filein, "mc", result_b_70_100, "bSV")
apply(filein, "mc", result_b_100_140, "bSV")
apply(filein, "mc", result_b_140_200, "bSV")
apply(filein, "mc", result_b_200_300, "bSV")
apply(filein, "mc", result_b_300_670, "bSV")
apply(filein, "mc", result_b_670_1000, "bSV")
apply(filein, "mc", result_b_dummy, "cSV")
apply(filein, "mc", result_b_dummy, "lSV")

apply(filein, "mc", result_c_20_30, "cSV")
apply(filein, "mc", result_c_30_50, "cSV")
apply(filein, "mc", result_c_50_70, "cSV")
apply(filein, "mc", result_c_70_100, "cSV")
apply(filein, "mc", result_c_100_140, "cSV")
apply(filein, "mc", result_c_140_200, "cSV")
apply(filein, "mc", result_c_200_300, "cSV")
apply(filein, "mc", result_c_300_670, "cSV")
apply(filein, "mc", result_c_670_1000, "cSV")
apply(filein, "mc", result_c_dummy, "bSV")
apply(filein, "mc", result_c_dummy, "lSV")

apply(filein, "mc", result_l_20_30, "lSV")
apply(filein, "mc", result_l_30_50, "lSV")
apply(filein, "mc", result_l_50_70, "lSV")
apply(filein, "mc", result_l_70_100, "lSV")
apply(filein, "mc", result_l_100_140, "lSV")
apply(filein, "mc", result_l_140_200, "lSV")
apply(filein, "mc", result_l_200_300, "lSV")
apply(filein, "mc", result_l_300_670, "lSV")
apply(filein, "mc", result_l_670_1000, "lSV")
apply(filein, "mc", result_l_dummy, "bSV")
apply(filein, "mc", result_l_dummy, "cSV")
