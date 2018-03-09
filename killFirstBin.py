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


def kill(filein, dirin, match) :
   result={}    
   thedir = filein.Get(dirin) 
   thedir.cd()
   nextiter =TIter(thedir.GetListOfKeys())
   print thedir.GetListOfKeys().LastIndex()
   i=0
   plotsToUpdate=[]
   while (True):
      key = nextiter();
      if key == None: break;
      cl = gROOT.GetClass(key.GetClassName())
      if not cl.InheritsFrom("TH1"): continue
      h = key.ReadObj();
      if match in h.GetName():
        #print h.GetName()
        h.SetBinContent(1, 0.)
        h.SetBinError(1, 0.)
        #h.Write(h.GetName(), TObject.kOverwrite)
      plotsToUpdate.append(h)
      i+=1
      if (i%1000==0): print "processed",i,"from",thedir.GetListOfKeys().LastIndex()
   for plot in plotsToUpdate:
     plot.Write(plot.GetName(), TObject.kOverwrite)
    
filein=TFile(sys.argv[1], "update")

print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>MC"
kill(filein, "mc", "SVmass_")
print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>DATA"
kill(filein, "data", "SVmass_")
print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>DATA_INVERTED"
kill(filein, "data_inverted", "SVmass_")
