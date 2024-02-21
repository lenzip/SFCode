#!/usr/bin/env python

from ROOT import *

import sys

systematics =  ["_JER_do","_JER_up",
                "_JES_do","_JES_up",
                "_gluonSplitting_do","_gluonSplitting_up",
                "_bFragmentation_do","_bFragmentation_up",
                "_cFragmentation_do","_cFragmentation_up",
                "_cdFragmentation_do","_cdFragmentation_up",
                "_v0_do","_v0_up"]


def doMove (filein, fileout, dirout, normalize) :
   directory = fileout.mkdir(dirout) 
   nextiter =TIter(filein.GetListOfKeys())
   while (True):
      key = nextiter(); 
      if key == None: break;
      cl = gROOT.GetClass(key.GetClassName())
      if not cl.InheritsFrom("TH1"): continue
      h = key.ReadObj();
      if (normalize and h.Integral()>0):
        print(h.GetName())
        for ext in systematics:
          if ext in h.GetName():
            print("this is a systematic variation", ext)
            # get the central
            hcentral = filein.Get(h.GetName().replace(ext, ""))
            h.Scale(hcentral.Integral()/h.Integral());
            break
        else:
          print("this is the central")
      directory.cd();
      h.Write();

def merge(filenamein, filenameout, dirname, normalize) :
  fileout = TFile(filenameout, "update");
  filein = TFile(filenamein);

  doMove(filein, fileout, dirname, normalize);

  filein.Close();  
  fileout.Close();


merge(sys.argv[1], sys.argv[2], sys.argv[3], bool(sys.argv[4]))
