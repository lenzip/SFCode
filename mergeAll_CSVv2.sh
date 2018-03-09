#!/bin/bash

rm total_CSVv2_2017.root
python merge.py ../../BTagAnalysis/plots_data_CSVv2/Run2017.root total_CSVv2_2017.root data True
python merge.py ../../BTagAnalysis/plots_data_inverted_CSVv2/Run2017_MCJP.root total_CSVv2_2017.root data_inverted True
python merge.py ../../BTagAnalysis/plots_mc_runBCDEF_CSVv2/QCD_MuEnriched.root total_CSVv2_2017.root mc True

rm total_CSVv2_2017B.root
python merge.py ../../BTagAnalysis/plots_data_CSVv2/Run2017B.root total_CSVv2_2017B.root data True
python merge.py ../../BTagAnalysis/plots_data_inverted_CSVv2/Run2017B_MCJP.root total_CSVv2_2017B.root data_inverted True
python merge.py ../../BTagAnalysis/plots_mc_runBCDEF_CSVv2/QCD_MuEnriched.root total_CSVv2_2017B.root mc True

rm total_CSVv2_2017CDElo.root
python merge.py ../../BTagAnalysis/plots_data_CSVv2/Run2017CDElo.root total_CSVv2_2017CDElo.root data True
python merge.py ../../BTagAnalysis/plots_data_inverted_CSVv2/Run2017CDElo_MCJP.root total_CSVv2_2017CDElo.root data_inverted True
python merge.py ../../BTagAnalysis/plots_mc_runBCDEF_CSVv2/QCD_MuEnriched.root total_CSVv2_2017CDElo.root mc True

rm total_CSVv2_2017EhiF.root
python merge.py ../../BTagAnalysis/plots_data_CSVv2/Run2017EhiF.root total_CSVv2_2017EhiF.root data True
python merge.py ../../BTagAnalysis/plots_data_inverted_CSVv2/Run2017EhiF_MCJP.root total_CSVv2_2017EhiF.root data_inverted True
python merge.py ../../BTagAnalysis/plots_mc_runBCDEF_CSVv2/QCD_MuEnriched.root total_CSVv2_2017EhiF.root mc True



