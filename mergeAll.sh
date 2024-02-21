#!/bin/bash

rm total_2022EE.root
python3 merge.py ../BTagAnalysis_RDF/Data/BTagMu_Run2022EE_rdf_PNetBDisc.root total_2022EE.root data True
python3 merge.py ../BTagAnalysis_RDF/Data_inverted/BTagMu_Run2022EE_rdf_PNetBDisc.root total_2022EE.root data_inverted True
python3 merge.py ../BTagAnalysis_RDF/mc/QCD_MuEnriched_rdf_PNetBDisc.root total_2022EE.root mc True
