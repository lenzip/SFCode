#include "../CFIT/cfit.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TSystem.h"
#include "TMath.h"
#include <string>
using namespace std;

void getResults(CFIT::cfit *cf,float *par,float *err);


double computeSF(CFIT::cfit *cfSV, CFIT::cfit *cfJP, CFIT::cfit *cfJPtagged, bool debug=false){
   float parSV[100];
   float errSV[100];
   float par_tagSV[100];
   float err_tagSV[100];

   float parJP[100];
   float errJP[100];
   float par_tagJP[100];
   float err_tagJP[100]; 
 
   float parJPtagged[100];
   float errJPtagged[100];
   float par_tagJPtagged[100];
   float err_tagJPtagged[100];  
 
   cfJP->Run();   
   getResults(cfJP,parJP,errJP);
   float nmc1 = cfJP->GetNTemplate("b");
   float ndata1 = nmc1*parJP[0];

   if (debug) {
     cout << "JP total fit, used to get the SV finding efficiency" << endl;
     cout << "p (b-jet) " << parJP[0] << "\\pm" << errJP[0] << endl;
     cout << "p (c-jet) " << parJP[1] << "\\pm" << errJP[1] << endl;
     cout << "p (l-jet) " << parJP[2] << "\\pm" << errJP[2] << endl;
     cout << "n B jets in mc   " << nmc1 << endl;
     cout << "n B jets in data " << ndata1 << endl;
     cout << "N_{dof} " << cfJP->GetNDOF() << endl;
     cout << "\\Chi^2 " << cfJP->GetChisq() << endl;
   } 

   cfJP->Run("tag");
   getResults(cfJP,par_tagJP,err_tagJP);
   float nmc1_SV = cfJP->GetNTemplate("b");
   float ndata1_SV = nmc1_SV*par_tagJP[0]; 
   float effMC_SV = nmc1_SV/nmc1;
   float effDATA_SV = ndata1_SV/ndata1;
   if (debug) {
     cout << "JP tag fit, used to get the SV finding efficiency" << endl;
     cout << "p (b-jet) " << parJP[0] << "\\pm" << errJP[0] << endl;
     cout << "p (c-jet) " << parJP[1] << "\\pm" << errJP[1] << endl;
     cout << "p (l-jet) " << parJP[2] << "\\pm" << errJP[2] << endl;
     cout << "n B jets in mc   with a SV " << nmc1_SV << endl;
     cout << "n B jets in data with a SV " << ndata1_SV << endl;
     cout << "N_{dof} " << cfJP->GetNDOF() << endl;
     cout << "\\Chi^2 " << cfJP->GetChisq() << endl;
     cout << "SV finding efficiency in data " << effDATA_SV << endl;
     cout << "SV finding efficiency in MC   " << effMC_SV << endl;
   } 


   cfJPtagged->Run();
   getResults(cfJPtagged,parJPtagged,errJPtagged);
   float nmc1_tagged = cfJPtagged->GetNTemplate("b");
   float ndata1_tagged = nmc1_tagged*parJPtagged[0];

   if (debug) {
     cout << "JPtagged total fit, used to get the SV finding efficiency in tagged jets" << endl;
     cout << "p (b-jet) " << parJPtagged[0] << "\\pm" << errJPtagged[0] << endl;
     cout << "p (c-jet) " << parJPtagged[1] << "\\pm" << errJPtagged[1] << endl;
     cout << "p (l-jet) " << parJPtagged[2] << "\\pm" << errJPtagged[2] << endl;
     cout << "n B jets in mc tagged   " << nmc1_tagged << endl;
     cout << "n B jets in data tagged " << ndata1_tagged << endl;
     cout << "N_{dof} " << cfJPtagged->GetNDOF() << endl;
     cout << "\\Chi^2 " << cfJPtagged->GetChisq() << endl;
   }

   cfJPtagged->Run("tag");
   getResults(cfJPtagged,par_tagJPtagged,err_tagJPtagged);
   float nmc1_tagged_SV = cfJPtagged->GetNTemplate("b");
   float ndata1_tagged_SV = nmc1_tagged_SV*par_tagJPtagged[0];
   float effMC_SV_tagged = nmc1_tagged_SV/nmc1_tagged;
   float effDATA_SV_tagged = ndata1_tagged_SV/ndata1_tagged;
   if (debug) {
     cout << "JPtagged tag fit, used to get the SV finding efficiency in tagged jets" << endl;
     cout << "p (b-jet) " << parJPtagged[0] << "\\pm" << errJPtagged[0] << endl;
     cout << "p (c-jet) " << parJPtagged[1] << "\\pm" << errJPtagged[1] << endl;
     cout << "p (l-jet) " << parJPtagged[2] << "\\pm" << errJPtagged[2] << endl;
     cout << "n B jets in mc   with a SV " << nmc1_tagged_SV << endl;
     cout << "n B jets in data with a SV " << ndata1_tagged_SV << endl;
     cout << "N_{dof} " << cfJPtagged->GetNDOF() << endl;
     cout << "\\Chi^2 " << cfJPtagged->GetChisq() << endl;
     cout << "SV finding efficiency in data for tagged jets " << effDATA_SV_tagged << endl;
     cout << "SV finding efficiency in MC   for tagged jets " << effMC_SV_tagged << endl;
   } 


   cfSV->Run();
   getResults(cfSV,parSV,errSV);
   float nmc1_SV_SVfit = cfSV->GetNTemplate("b");
   float ndata1_SV_SVfit    = nmc1_SV_SVfit*parSV[0]; 
   if (debug) {
     cout << "SV total fit, used to find the total number of events with a SV" << endl;
     cout << "p (b-jet) " << parSV[0] << "\\pm" << errSV[0] << endl;
     cout << "p (c-jet) " << parSV[1] << "\\pm" << errSV[1] << endl;
     cout << "p (l-jet) " << parSV[2] << "\\pm" << errSV[2] << endl;
     cout << "n B jets in mc   with a SV " << nmc1_SV_SVfit << endl;
     cout << "n B jets in data with a SV " << ndata1_SV_SVfit << endl;
     cout << "N_{dof} " << cfSV->GetNDOF() << endl;
     cout << "\\Chi^2 " << cfSV->GetChisq() << endl;
   }
  


   cfSV->Run("tag");
   getResults(cfSV,par_tagSV,err_tagSV);
   float nmc1_tagged_SV_SVfit = cfSV->GetNTemplate("b");
   float ndata1_tagged_SV_SVfit = nmc1_tagged_SV_SVfit*par_tagSV[0];
   if (debug) {
     cout << "SV tag fit,  used to find the total number of events with a SV that are also tagged" << endl;
     cout << "p (b-jet) tag " << par_tagSV[0] << "\\pm" << err_tagSV[0] << endl;
     cout << "p (c-jet) tag " << par_tagSV[1] << "\\pm" << err_tagSV[1] << endl;
     cout << "p (l-jet) tag " << par_tagSV[2] << "\\pm" << err_tagSV[2] << endl;
     cout << "n B jets in mc   with a SV and tagged " << nmc1_tagged_SV_SVfit << endl;
     cout << "n B jets in data with a SV and tagged " << ndata1_tagged_SV_SVfit << endl;
     cout << "N_{dof} " << cfSV->GetNDOF() << endl;
     cout << "\\Chi^2 " << cfSV->GetChisq() << endl;
   }

   float nmc1_tag_estim = nmc1_tagged_SV_SVfit/effMC_SV_tagged;
   float nmc1_estim     = nmc1_SV_SVfit/effMC_SV;  
   float effMC = nmc1_tag_estim/nmc1_estim;

   float ndata1_tag_estim = ndata1_tagged_SV_SVfit/effDATA_SV_tagged;
   float ndata1_estim     = ndata1_SV_SVfit/effDATA_SV; 
   float effDATA = ndata1_tag_estim/ndata1_estim;
  
   if (debug) { 
     std::cout << "effMC = " << effMC << std::endl;
     std::cout << "effDATA = " << effDATA << std::endl;
     std::cout << "sf = " << effDATA/effMC << std::endl;
   } 
  
   return effDATA/effMC; 

}




void computeSFsWithJPandSV(const char* filein, const char* wp, const char* outdir)
{
   gROOT->SetBatch();
    
    
   gSystem->Load("../CFIT/libCFIT.so");
   //gROOT->SetStyle("Plain");
   //gROOT->ProcessLine(".L ~/tdrStyle.C");
   //setTDRStyle();


   float par[100];
   float err[100];
   float par_tag[100];
   float err_tag[100];

   string bins[9] = { "20.0-30.0",
                      "30.0-50.0", 
                      "50.0-70.0",
                      "70.0-100.0",
                      "100.0-140.0",
                      "140.0-200.0",
                      "200.0-300.0",
                      "300.0-670.0",
                      "670.0-1000.0"
                      };
   //double correlations[9] = {0.86011, 0.921278, 0.9391, 0.943887, 0.939105, 0.928875, 0.911364, 0.887669, 0.851008};                   
   float correlations[9] = {0.84289,0.888427,0.900466,0.903794,0.897022,0.881543,0.845106,0.782661,0.612723};                   
   // array of SF values per bin                   
   double v_sf[9] = {1.,1.,1.,1.,1.,1.,1.,1.,1.};
   // array of stat errors per bin
   double v_sfstat[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
   // array of total errors up per bin
   double v_sftotup[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
   //array of total errors down per bin
   double v_sftotdo[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};

   // map of syst errors per bin
   double v_sf_syst[23][9];

   // arrays of x values
   double v_pt[9]={25., 40., 60., 85., 120., 170., 250., 485., 835.};
   double v_pterr[9] = {5, 10, 10, 15, 20, 30, 50, 185, 165};

   string systematics[23] = {  "_JER_do","_JER_up",
                               "_JES_do","_JES_up",
                               "_gluonSplitting_do","_gluonSplitting_up",
                               "_bFragmentation_do","_bFragmentation_up",
                               "_cFragmentation_do","_cFragmentation_up",
                               "_cdFragmentation_do","_cdFragmentation_up",
                               "_v0_do","_v0_up",
                               "_cSVmass_do", "_cSVmass_up",
                               "_morph1",
                               "_morph2",
                               "_morph3",
                               "_morph4",
                               "_JPMC",
                               "_Cb_up",
                               "_Cb_do"};


   for (unsigned ibin = 0; ibin < 9; ++ibin){
     std::cout << ">>>>>>>>>>>>>>>> Running bin " << bins[ibin] << std::endl;
  
     CFIT::cfit *cfSV = new CFIT::cfit("SV");
     cfSV->SetPlotDirName("picsSV");
     CFIT::cfit *cfJP = new CFIT::cfit("JP");
     cfJP->SetPlotDirName("picsJP");
     CFIT::cfit *cfJPtagged = new CFIT::cfit("JPtagged");
     cfJPtagged->SetPlotDirName("picsJPtagged");
     
     cfSV->SetOptimization(OPT_MORPH_SGN_SIGMA);
     cfSV->SetMorphing(OPTMORPH_CUTOFF,0.5);
     cfSV->ProducePlots(1);
     cfJP->SetOptimization(OPT_MORPH_SGN_SIGMA);
     cfJP->SetMorphing(OPTMORPH_CUTOFF,0.5);
     cfJP->ProducePlots(1);
     cfJPtagged->SetOptimization(OPT_MORPH_SGN_SIGMA);
     cfJPtagged->SetMorphing(OPTMORPH_CUTOFF,0.5);
     cfJPtagged->ProducePlots(1);  
     
     //cf->SetInputFile("total_syst.root");
     cfSV->SetInputFile(string(filein));
     cfJP->SetInputFile(string(filein));
     cfJPtagged->SetInputFile(string(filein));

     cfSV->AddSys("JER","_JER_do","_JER_up");
     cfSV->AddSys("JES","_JES_do","_JES_up");
     cfSV->AddSys("gluonSplitting","_gluonSplitting_do","_gluonSplitting_up");
     cfSV->AddSys("cdFragmentation","_cdFragmentation_do","_cdFragmentation_up");
     cfSV->AddSys("cFragmentation","_cFragmentation_do","_cFragmentation_up");
     cfSV->AddSys("bFragmentation","_bFragmentation_do","_bFragmentation_up");
     cfSV->AddSys("v0","_v0_do","_v0_up");
     cfSV->AddSys("cSVmass", "_cSVmass_do", "_cSVmass_up");
  
     cfJP->AddSys("JER","_JER_do","_JER_up");
     cfJP->AddSys("JES","_JES_do","_JES_up");
     cfJP->AddSys("gluonSplitting","_gluonSplitting_do","_gluonSplitting_up");
     cfJP->AddSys("cdFragmentation","_cdFragmentation_do","_cdFragmentation_up");
     cfJP->AddSys("cFragmentation","_cFragmentation_do","_cFragmentation_up");
     cfJP->AddSys("bFragmentation","_bFragmentation_do","_bFragmentation_up");
     cfJP->AddSys("v0","_v0_do","_v0_up");
     cfJP->AddSys("cSVmass", "_cSVmass_do", "_cSVmass_up");
 
     cfJPtagged->AddSys("JER","_JER_do","_JER_up");
     cfJPtagged->AddSys("JES","_JES_do","_JES_up");
     cfJPtagged->AddSys("gluonSplitting","_gluonSplitting_do","_gluonSplitting_up");
     cfJPtagged->AddSys("cdFragmentation","_cdFragmentation_do","_cdFragmentation_up");
     cfJPtagged->AddSys("cFragmentation","_cFragmentation_do","_cFragmentation_up");
     cfJPtagged->AddSys("bFragmentation","_bFragmentation_do","_bFragmentation_up");
     cfJPtagged->AddSys("v0","_v0_do","_v0_up");
     cfJPtagged->AddSys("cSVmass", "_cSVmass_do", "_cSVmass_up");
     
     const char* thebin = bins[ibin].c_str(); 

     cfSV->SetMatrixName("matrixSV");
     cfSV->SetMatrixOption("WRITE");   

     cfSV->SetData(Form("data/SVmass__ptbin_%s", thebin));
     cfSV->SetDataTag(Form("data/SVmass__ptbin_%s_DeepFlavour%s", thebin, wp));
     cfSV->SetDataUntag(Form("data/SVmass__ptbin_%s_DeepFlavour%s_Fail", thebin, wp));

     cfSV->AddTemplate("b",Form("mc/SVmass__ptbin_%s_b", thebin),2);
     cfSV->AddTemplate("c",Form("mc/SVmass__ptbin_%s_c", thebin),3);
     cfSV->AddTemplate("l",Form("mc/SVmass__ptbin_%s_l", thebin),4);

     cfSV->AddTemplateTag("b",Form("mc/SVmass__ptbin_%s_b_DeepFlavour%s", thebin, wp),2);
     cfSV->AddTemplateTag("c",Form("mc/SVmass__ptbin_%s_c_DeepFlavour%s", thebin, wp),3);
     cfSV->AddTemplateTag("l",Form("mc/SVmass__ptbin_%s_l_DeepFlavour%s", thebin, wp),4);

     cfSV->AddTemplateUntag("b",Form("mc/SVmass__ptbin_%s_b_DeepFlavour%s_Fail", thebin, wp),2);
     cfSV->AddTemplateUntag("c",Form("mc/SVmass__ptbin_%s_c_DeepFlavour%s_Fail", thebin, wp),3);
     cfSV->AddTemplateUntag("l",Form("mc/SVmass__ptbin_%s_l_DeepFlavour%s_Fail", thebin, wp),4);   


     cfJP->SetMatrixName("matrixJP");
     cfJP->SetMatrixOption("WRITE");

     cfJP->SetData(Form("data/JP__ptbin_%s", thebin));
     cfJP->SetDataTag(Form("data/JP__ptbin_%s_SV", thebin));
     cfJP->SetDataUntag(Form("data/JP__ptbin_%s_NoSV", thebin));
     //cfJP->SetDataTag(Form("data/JP__ptbin_%s_DeepFlavour%s", thebin, wp));
     //cfJP->SetDataUntag(Form("data/JP__ptbin_%s_DeepFlavour%s_Fail", thebin, wp));

     cfJP->AddTemplate("b",Form("mc/JP__ptbin_%s_b", thebin),2);
     cfJP->AddTemplate("c",Form("mc/JP__ptbin_%s_c", thebin),3);
     cfJP->AddTemplate("l",Form("mc/JP__ptbin_%s_l", thebin),4);

     //cfJP->AddTemplateTag("b",Form("mc/JP__ptbin_%s_b_DeepFlavour%s", thebin, wp),2);
     //cfJP->AddTemplateTag("c",Form("mc/JP__ptbin_%s_c_DeepFlavour%s", thebin, wp),3);
     //cfJP->AddTemplateTag("l",Form("mc/JP__ptbin_%s_l_DeepFlavour%s", thebin, wp),4);
     cfJP->AddTemplateTag("b",Form("mc/JP__ptbin_%s_b_SV", thebin),2);
     cfJP->AddTemplateTag("c",Form("mc/JP__ptbin_%s_c_SV", thebin),3);
     cfJP->AddTemplateTag("l",Form("mc/JP__ptbin_%s_l_SV", thebin),4);

     cfJP->AddTemplateUntag("b",Form("mc/JP__ptbin_%s_b_NoSV", thebin),2);
     cfJP->AddTemplateUntag("c",Form("mc/JP__ptbin_%s_c_NoSV", thebin),3);
     cfJP->AddTemplateUntag("l",Form("mc/JP__ptbin_%s_l_NoSV", thebin),4);


     cfJPtagged->SetMatrixName("matrixJPtagged");
     cfJPtagged->SetMatrixOption("WRITE");

     cfJPtagged->SetData(Form("data/JP__ptbin_%s_DeepFlavour%s", thebin, wp));
     cfJPtagged->SetDataTag(Form("data/JP__ptbin_%s_DeepFlavour%sSV", thebin,wp));
     cfJPtagged->SetDataUntag(Form("data/JP__ptbin_%s_DeepFlavour%sNoSV", thebin,wp));

     cfJPtagged->AddTemplate("b",Form("mc/JP__ptbin_%s_b_DeepFlavour%s", thebin, wp),2);
     cfJPtagged->AddTemplate("c",Form("mc/JP__ptbin_%s_c_DeepFlavour%s", thebin, wp),3);
     cfJPtagged->AddTemplate("l",Form("mc/JP__ptbin_%s_l_DeepFlavour%s", thebin, wp),4);

     cfJPtagged->AddTemplateTag("b",Form("mc/JP__ptbin_%s_b_DeepFlavour%s_SV", thebin, wp),2);
     cfJPtagged->AddTemplateTag("c",Form("mc/JP__ptbin_%s_c_DeepFlavour%s_SV", thebin, wp),3);
     cfJPtagged->AddTemplateTag("l",Form("mc/JP__ptbin_%s_l_DeepFlavour%s_SV", thebin, wp),4);

     cfJPtagged->AddTemplateUntag("b",Form("mc/JP__ptbin_%s_b_DeepFlavour%s_NoSV", thebin, wp),2);
     cfJPtagged->AddTemplateUntag("c",Form("mc/JP__ptbin_%s_c_DeepFlavour%s_NoSV", thebin, wp),3);
     cfJPtagged->AddTemplateUntag("l",Form("mc/JP__ptbin_%s_l_DeepFlavour%s_NoSV", thebin, wp),4);

     double sfcentral = computeSF(cfSV, cfJP, cfJPtagged, true); 
     system(Form("mv picsSV picsSV_%s", thebin));
     system(Form("mv picsJP picsJP_%s", thebin));
     system(Form("mv picsJPtagged picsJPtagged_%s", thebin));

     // these are the 19 systematic variations for this bin
     double sfvars[23]={}; 

     cout << "sf_central " <<  sfcentral << endl;
     // perform statistical variation
    
     double statSum = 0;
     double statSqSum = 0; 
     int nVariations = 10;
     for(int is=0;is<nVariations;is++)
     { 
        cout << "doing stat variation: " << is << endl; 
        int isys = 666+is;
        cfJP->SetMatrixOption("READ");
        cfJP->ProducePlots(0);
        cfJP->SetStatVariation(isys);
        cfJPtagged->SetMatrixOption("READ");
        cfJPtagged->ProducePlots(0);
        cfJPtagged->SetStatVariation(isys);
        cfSV->SetMatrixOption("READ");
        cfSV->ProducePlots(0);
        cfSV->SetStatVariation(isys);
        double sfstat = computeSF(cfSV, cfJP, cfJPtagged);
        statSum += sfstat;
        statSqSum += sfstat*sfstat;
  
     }  
     double averageSF = statSum/nVariations;
     double rmsstat =  sqrt(statSqSum/nVariations - averageSF*averageSF);
    

     // perform systematic variation
     //
     for (unsigned int i = 0; i < 16; ++i){
      cout << "doing systematic variation " << systematics[i] << endl;
      cfJP->SetMatrixOption("READ");
      cfJP->SetSysVariation(systematics[i]);
      cfJPtagged->SetMatrixOption("READ");
      cfJPtagged->SetSysVariation(systematics[i]);
      cfSV->SetMatrixOption("READ");
      cfSV->SetSysVariation(systematics[i]);
      double sf = computeSF(cfSV, cfJP, cfJPtagged);
      sfvars[i] = sf;
      v_sf_syst[i][ibin] = sf;
      //cout << "sf" << systematics[i] << " "<< sf << endl;
     } 
   
     double errdosq = 0;
     double errupsq = 0;
     for (unsigned int i = 0; i < 16; ++i){
       //if (i%2 == 0)
       if ((sfvars[i] - sfcentral) < 0.)
         errdosq += (sfvars[i] - sfcentral)*(sfvars[i] - sfcentral);
       else
         errupsq += (sfvars[i] - sfcentral)*(sfvars[i] - sfcentral);
     }
     for (unsigned int ifit = 0; ifit < 7; ++ifit){
        cout << "doing systematic variation " << systematics[16+ifit] << endl;
       cfJP->SetMatrixOption("WRITE");
       cfJP->ProducePlots(0);
       cfJP->SetOptimization(OPT_MORPH_SGN_SIGMA);
       cfJP->SetMorphing(OPTMORPH_CUTOFF,0.5);
       cfJP->SetData(Form("data/JP__ptbin_%s", thebin));
       cfJP->SetDataTag(Form("data/JP__ptbin_%s_SV", thebin));
       cfJP->SetDataUntag(Form("data/JP__ptbin_%s_NoSV", thebin));
       cfJPtagged->SetData(Form("data/JP__ptbin_%s_DeepFlavour%s", thebin, wp));
       cfJPtagged->SetDataTag(Form("data/JP__ptbin_%s_DeepFlavour%sSV", thebin,wp));
       cfJPtagged->SetDataUntag(Form("data/JP__ptbin_%s_DeepFlavour%sNoSV", thebin,wp)); 
       cfSV->SetMatrixOption("WRITE");
       cfSV->ProducePlots(0);
       cfSV->SetOptimization(OPT_MORPH_SGN_SIGMA);
       cfSV->SetMorphing(OPTMORPH_CUTOFF,0.5);
       cfSV->SetData(Form("data/SVmass__ptbin_%s", thebin));
       cfSV->SetDataTag(Form("data/SVmass__ptbin_%s_DeepFlavour%s", thebin, wp));
       cfSV->SetDataUntag(Form("data/SVmass__ptbin_%s_DeepFlavour%s_Fail", thebin, wp)); 
       if (ifit == 0) {cfSV->SetOptimization(OPT_NOCORR); cfJP->SetOptimization(OPT_NOCORR); cfJPtagged->SetOptimization(OPT_NOCORR);}
       if (ifit == 1) {cfSV->SetMorphing(OPTMORPH_CUTOFF,0.25); cfJP->SetMorphing(OPTMORPH_CUTOFF,0.25); cfJPtagged->SetMorphing(OPTMORPH_CUTOFF,0.25);}
       if (ifit == 2) {cfSV->SetMorphing(OPTMORPH_CUTOFF,0.75); cfJP->SetMorphing(OPTMORPH_CUTOFF,0.75); cfJPtagged->SetMorphing(OPTMORPH_CUTOFF,0.75);}
       if (ifit == 3) {cfSV->SetMorphing(OPTMORPH_GEOMETRIC); cfJP->SetMorphing(OPTMORPH_GEOMETRIC); cfJPtagged->SetMorphing(OPTMORPH_GEOMETRIC);}
       if (ifit == 4) { //use MC JP calib
        cfJP->SetData(Form("data_inverted/JP__ptbin_%s", thebin));
        cfJP->SetDataTag(Form("data_inverted/JP__ptbin_%s_SV", thebin));
        cfJP->SetDataUntag(Form("data_inverted/JP__ptbin_%s_NoSV", thebin));
        cfJPtagged->SetData(Form("data_inverted/JP__ptbin_%s_DeepFlavour%s", thebin, wp));
        cfJPtagged->SetDataTag(Form("data_inverted/JP__ptbin_%s_DeepFlavour%sSV", thebin,wp));
        cfJPtagged->SetDataUntag(Form("data_inverted/JP__ptbin_%s_DeepFlavour%sNoSV", thebin,wp));
        cfSV->SetData(Form("data_inverted/SVmass__ptbin_%s", thebin));
        cfSV->SetDataTag(Form("data_inverted/SVmass__ptbin_%s_DeepFlavour%s", thebin, wp));
        cfSV->SetDataUntag(Form("data_inverted/SVmass__ptbin_%s_DeepFlavour%s_Fail", thebin, wp));
       } 
       if (ifit < 5){ 
         if (ifit == 4){
           system("mkdir picsSV");
           system("mkdir picsJP");
           system("mkdir picsJPtagged");
         }  
         double  sf = computeSF(cfSV, cfJP, cfJPtagged);
         if (ifit == 4){
            system(Form("mv picsSV picsSV_inverted_%s", thebin));
            system(Form("mv picsJP picsJP_inverted_%s", thebin));
            system(Form("mv picsJPtagged picsJPtagged_inverted_%s", thebin));
         }   
         sfvars[ifit+16] = sf;
         v_sf_syst[ifit+16][ibin] = sf;
         errdosq += (sf - sfcentral)*(sf - sfcentral);
         errupsq += (sf - sfcentral)*(sf - sfcentral);
       } else { // Cb variation
         sfvars[ifit+16] = sfcentral;
         TFile fileinraw(filein);
         double nbmc = ((TH1*)fileinraw.Get(Form("mc/SVmass__ptbin_%s_b", thebin)))->Integral(); 
         double nbmc_jp = nbmc - ((TH1*)fileinraw.Get(Form("mc/SVmass__ptbin_%s_b_JP0", thebin)))->Integral();
         double nbmctag = ((TH1*)fileinraw.Get(Form("mc/SVmass__ptbin_%s_b_DeepFlavour%s", thebin, wp)))->Integral() ;
         double nbmctag_jp = nbmctag - ((TH1*)fileinraw.Get(Form("mc/SVmass__ptbin_%s_b_DeepFlavour%s_JP0", thebin, wp)))->Integral();
         cout << "nbmc " << nbmc << " nbmc_jp " << nbmc_jp << " nbmctag " << nbmctag << "nbmctag_jp " << nbmctag_jp << endl;
         double Cb=(nbmctag/nbmctag_jp)*(nbmc_jp/nbmc);
         double Cb_error = fabs(1-Cb)/2;
         cout << "Cb " << Cb << " Cb_error " << Cb_error << endl;
         if (ifit == 5){
           v_sf_syst[ifit+16][ibin] = sfcentral*(1+Cb_error);
           errupsq += (Cb_error*sfcentral)*(Cb_error*sfcentral);
         } 
         if (ifit == 6){
           v_sf_syst[ifit+16][ibin] = sfcentral*(1-Cb_error);
           errdosq += (Cb_error*sfcentral)*(Cb_error*sfcentral);
         }
       }
     }

     cout << "sf " << sfcentral << " +"<< sqrt(errupsq) << " -" <<  sqrt(errdosq) << "(syst) +/-" << rmsstat << "(stat)" << endl; 
     v_sf[ibin] = sfcentral;
     v_sfstat[ibin] = rmsstat;
     v_sftotup[ibin] = sqrt(v_sfstat[ibin]*v_sfstat[ibin] + errupsq);
     v_sftotdo[ibin] = sqrt(v_sfstat[ibin]*v_sfstat[ibin] + errdosq);

     delete cfSV;
     delete cfJP;
   }
  
   TFile output(Form("DeepFlavour%s.root", wp), "recreate"); 

   TGraphErrors gcentral(9, v_pt, v_sf, v_pterr, v_sfstat);
   gcentral.SetNameTitle("statistical", "statistical");
   TGraphAsymmErrors gtotal(9, v_pt, v_sf, v_pterr, v_pterr, v_sftotdo, v_sftotup);
   gtotal.SetNameTitle("total", "total");
   gtotal.SetFillColor(kYellow);
   gtotal.SetFillStyle(3001);

   TCanvas * c = new TCanvas();

   c->SetLogx();

   gtotal.Draw("AE2");
   gcentral.Draw("P");
   c->SetGridx();
   c->SetGridy();
   gPad->RedrawAxis("g");

   c->Write();

   gcentral.Write();
   gtotal.Write();
   for (unsigned int i = 0; i < 23; ++i){
     TGraphErrors gsyst(9, v_pt, v_sf_syst[i], v_pterr, 0);
     gsyst.SetNameTitle(Form("variation%s", systematics[i].c_str()), Form("variation%s", systematics[i].c_str()));
     gsyst.Write();
   }


   system(Form("mkdir -p %s/%s", outdir, wp));
   system(Form("mv pics* %s/%s", outdir, wp));
   system(Form("mv DeepFlavour%s.root %s/%s", wp, outdir, wp));

   //gApplication->Terminate();
}

void getResults(CFIT::cfit *cf,float *par,float *err)
{   
   float chi2 = cf->GetChisq();
   int nPar = cf->GetNPar();
   for(int i=0;i<nPar;i++)
     {
	par[i] = cf->GetPar(i);
	err[i] = cf->GetParErr(i);	
     }
}
