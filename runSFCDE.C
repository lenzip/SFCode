{
   //gROOT->SetStyle("Plain");
   //gROOT->ProcessLine(".L ~/tdrStyle.C");
   //setTDRStyle();
   gROOT->Macro("load.C");
   //gROOT->Macro("computeSFs.C(\"total_2017.root\", \"L\", \"Full2017\")");
   //gROOT->Macro("computeSFs.C(\"total_2017.root\", \"M\", \"Full2017\")");
   //gROOT->Macro("computeSFs.C(\"total_2017.root\", \"T\", \"Full2017\")");
   //gROOT->Macro("computeSFs.C(\"total_2017B.root\", \"L\", \"Full2017B\")");
   //gROOT->Macro("computeSFs.C(\"total_2017B.root\", \"M\", \"Full2017B\")");
   //gROOT->Macro("computeSFs.C(\"total_2017B.root\", \"T\", \"Full2017B\")");
   gROOT->Macro("computeSFs_CSVv2.C(\"total_CSVv2_2017CDElo.root\", \"L\", \"Full2017CDElo_CSVv2\")");
   gROOT->Macro("computeSFs_CSVv2.C(\"total_CSVv2_2017CDElo.root\", \"M\", \"Full2017CDElo_CSVv2\")");
   gROOT->Macro("computeSFs_CSVv2.C(\"total_CSVv2_2017CDElo.root\", \"T\", \"Full2017CDElo_CSVv2\")");
   //gROOT->Macro("computeSFs.C(\"total_2017EhiF.root\", \"L\", \"Full2017EhiF\")");
   //gROOT->Macro("computeSFs.C(\"total_2017EhiF.root\", \"M\", \"Full2017EhiF\")");
   //gROOT->Macro("computeSFs.C(\"total_2017EhiF.root\", \"T\", \"Full2017EhiF\")");
}
