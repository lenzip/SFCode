#!/bin/bash

from ROOT import *
import sys
import math
from array import *

filein=TFile(sys.argv[1])
directory=sys.argv[2]
flavor=sys.argv[3]

output=TFile("correlations2.root", "update")

bins = [ "20.0-30.0",
         "30.0-50.0",
         "50.0-70.0",
         "70.0-100.0",
         "100.0-140.0",
         "140.0-200.0",
         "200.0-300.0",
         "300.0-670.0",
         "670.0-1000.0" ]

wps = ["L", "M", "T"]

x  = array("d", [25., 40., 60., 85., 120., 160., 250., 485., 835.])
ex = array("d", [5. , 10., 10., 15.,  20.,  30.,  50., 185., 165.])

for wp in wps:
  correlations = []
  ey = []
  for b in bins:
    print b
    nametotal  = directory+"/SVmass__ptbin_"+b+"_DeepCSVBDiscL"
    nametagged = directory+"/SVmass__ptbin_"+b+"_DeepCSVBDisc"+wp
    if flavor != "":
      nametotal  = directory+"/SVmass__ptbin_"+b+"_"+flavor+"_DeepCSVBDiscL"
      nametagged = directory+"/SVmass__ptbin_"+b+"_"+flavor+"_DeepCSVBDisc"+wp
    print nametotal, nametagged  
    total  = filein.Get(nametotal)
    tagged = filein.Get(nametagged)
    hasSV = total.Integral() - total.GetBinContent(1) 
    isTagged = tagged.Integral()
    isTaggedAndHasSV = tagged.Integral() - tagged.GetBinContent(1)

    correlation = isTaggedAndHasSV/math.sqrt(isTagged*hasSV)
    correlations.append(correlation)
    ey.append(0.)
  
    print "correlation WP", wp, "bin", b, correlation

  graph = TGraphErrors(9, x, array("d", correlations), ex, array("d", ey))
  graph.SetNameTitle("correlation"+wp+directory+flavor, "correlation"+wp+directory+flavor)
  output.cd()
  graph.Write()
  c=TCanvas()
  c.cd()
  graph.Draw("AP")
  c.SetLogx()
  c.SaveAs("correlation"+wp+directory+".png")
