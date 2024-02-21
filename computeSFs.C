#include "../CFIT/cfit.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TSystem.h"
#include <string>
#include <map>

using namespace std;

void getResults(CFIT::cfit *cf,float *par,float *err);


std::pair<double, double> computeSF(CFIT::cfit *cf, bool debug=true){
   float par[100];
   float err[100];
   float par_tag[100];
   float err_tag[100];
  
   cf->Run();   
   getResults(cf,par,err);
   float chi2 = cf->GetChisq();
   float ndata = cf->GetNData();
   float nmc1 = cf->GetNTemplate("b");
   float nmc = cf->GetNTemplate("b")+cf->GetNTemplate("c")+cf->GetNTemplate("l");

   if (debug) {
     cout << "p (b-jet) " << par[0] << "\\pm" << err[0] << endl;
     cout << "p (c-jet) " << par[1] << "\\pm" << err[1] << endl;
     cout << "p (l-jet) " << par[2] << "\\pm" << err[2] << endl;
     cout << "N_{dof} " << cf->GetNDOF() << endl;
     cout << "\\Chi^2 " << cf->GetChisq() << endl;
   } 

   cf->Run("tag");
   getResults(cf,par_tag,err_tag);
   float chi2_tag = cf->GetChisq();
   float ndata_tag = cf->GetNData();
   float nmc1_tag = cf->GetNTemplate("b");
   float nmc_tag = cf->GetNTemplate("b")+cf->GetNTemplate("c")+cf->GetNTemplate("l");

   float fr = nmc1/nmc;
   float fr_tag = nmc1_tag/nmc_tag;

   float effMC = nmc1_tag/nmc1;
   float effDATA = nmc_tag/nmc*par_tag[0]/par[0]*fr_tag/fr;

   if (debug) {
     cout << "p (b-jet) tag " << par_tag[0] << "\\pm" << err_tag[0] << endl;
     cout << "p (c-jet) tag " << par_tag[1] << "\\pm" << err_tag[1] << endl;
     cout << "p (l-jet) tag " << par_tag[2] << "\\pm" << err_tag[2] << endl;
     cout << "N_{dof} " << cf->GetNDOF() << endl;
     cout << "\\Chi^2 " << cf->GetChisq() << endl;
     std::cout << "effMC = " << effMC << std::endl;
     std::cout << "effDATA = " << effDATA << std::endl;
     std::cout << "sf = " << effDATA/effMC << std::endl;
   }
   double sf=effDATA/effMC;
   double rel_err = sqrt(std::pow(err[0]/par[0],2)+std::pow(err_tag[0]/par_tag[0],2));
  
   return std::make_pair(sf, rel_err*sf); 

}




void computeSFs(const char* filein, const char* wp, const char* outdir)
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

   string bins[9] = { "20_30",
                      "30_50", 
                      "50_70",
                      "70_100",
                      "100_140",
                      "140_200",
                      "200_300",
                      "300_600",
                      "600_1000"
                      };
   // array of SF values per bin                   
   double v_sf[9] = {1.,1.,1.,1.,1.,1.,1.,1.,1.};
   // array of stat errors per bin
   double v_sfstat[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
   // array of total errors up per bin
   double v_sftotup[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};
   //array of total errors down per bin
   double v_sftotdo[9] = {0.,0.,0.,0.,0.,0.,0.,0.,0.};

   // map of syst errors per bin
   double v_sf_syst[21][9];

   // arrays of x values
   double v_pt[9]={25., 40., 60., 85., 120., 170., 250., 450., 800.};
   double v_pterr[9] = {5, 10, 10, 15, 20, 30, 50, 150, 200};

   string systematics[21] = {  "_JER_do","_JER_up",
                               "_JES_do","_JES_up",
                               "_gluonSplitting_do","_gluonSplitting_up",
                               "_bFragmentation_do","_bFragmentation_up",
                               "_cFragmentation_do","_cFragmentation_up",
                               "_cdFragmentation_do","_cdFragmentation_up",
                               "_v0_do","_v0_up",
                               //"_cSVmass_do", "_cSVmass_up",
                               "_morph1",
                               "_morph2",
                               "_morph3",
                               "_morph4",
                               "_JPMC",
                               "_Cb_up",
                               "_Cb_do"};


   for (unsigned ibin = 0; ibin < 9; ++ibin){
     std::cout << ">>>>>>>>>>>>>>>> Running bin " << bins[ibin] << std::endl;
  
     CFIT::cfit *cf = new CFIT::cfit("JP");
     
     cf->SetOptimization(OPT_MORPH_SGN_SIGMA);
     cf->SetMorphing(OPTMORPH_CUTOFF,0.5);
     cf->ProducePlots(1);
     
     //cf->SetInputFile("total_syst.root");
     cf->SetInputFile(string(filein));

     cf->AddSys("JER","_JER_do","_JER_up");
     cf->AddSys("JES","_JES_do","_JES_up");
     cf->AddSys("gluonSplitting","_gluonSplitting_do","_gluonSplitting_up");
     cf->AddSys("cdFragmentation","_cdFragmentation_do","_cdFragmentation_up");
     cf->AddSys("cFragmentation","_cFragmentation_do","_cFragmentation_up");
     cf->AddSys("bFragmentation","_bFragmentation_do","_bFragmentation_up");
     cf->AddSys("v0","_v0_do","_v0_up");
     //cfSV->AddSys("cSVmass", "_cSVmass_do", "_cSVmass_up"); 
        
     const char* thebin = bins[ibin].c_str(); 


     cf->SetMatrixOption("WRITE");   

     cf->SetData(Form("data/JP__ptbin_%s", thebin));
     cf->SetDataTag(Form("data/JP__ptbin_%s_PNetBDisc%s", thebin, wp));
     cf->SetDataUntag(Form("data/JP__ptbin_%s_PNetBDisc%s_Fail", thebin, wp));

     cf->AddTemplate("b",Form("mc/JP__ptbin_%s_b", thebin),2);
     cf->AddTemplate("c",Form("mc/JP__ptbin_%s_c", thebin),3);
     cf->AddTemplate("l",Form("mc/JP__ptbin_%s_l", thebin),4);

     cf->AddTemplateTag("b",Form("mc/JP__ptbin_%s_b_PNetBDisc%s", thebin, wp),2);
     cf->AddTemplateTag("c",Form("mc/JP__ptbin_%s_c_PNetBDisc%s", thebin, wp),3);
     cf->AddTemplateTag("l",Form("mc/JP__ptbin_%s_l_PNetBDisc%s", thebin, wp),4);

     cf->AddTemplateUntag("b",Form("mc/JP__ptbin_%s_b_PNetBDisc%s_Fail", thebin, wp),2);
     cf->AddTemplateUntag("c",Form("mc/JP__ptbin_%s_c_PNetBDisc%s_Fail", thebin, wp),3);
     cf->AddTemplateUntag("l",Form("mc/JP__ptbin_%s_l_PNetBDisc%s_Fail", thebin, wp),4);   

     //double sfcentral = computeSF(cf, true); 
     std::pair<double, double> sfcentral = computeSF(cf, true);
     system(Form("mv pics pics_%s", thebin));

     // these are the 19 systematic variations for this bin
     double sfvars[19]={}; 

     //cout << "sf_central " <<  sfcentral << endl;
     // perform statistical variation
    
    
     double statSum = 0;
     double statSqSum = 0; 
     int nVariations = 100;
     for(int is=0;is<nVariations;is++)
     {  
        int isys = 666+is;
        cf->SetMatrixOption("READ");
        cf->ProducePlots(0);
        cf->SetStatVariation(isys);
        //double sfstat = computeSF(cf, false);
        std::pair<double, double> sfstat = computeSF(cf, false);
        statSum += sfstat.first;
        statSqSum += sfstat.first*sfstat.first;
  
     }  
     double averageSF = statSum/nVariations;
     double rmsstat =  sqrt(statSqSum/nVariations - averageSF*averageSF);
     // perform systematic variation
     //
     for (unsigned int i = 0; i < 14; ++i){ 
      cf->SetMatrixOption("READ");
      cf->SetSysVariation(systematics[i]);
      //double sf = computeSF(cf, false);
      std::pair<double,double> sf = computeSF(cf, false);
      sfvars[i] = sf.first;
      v_sf_syst[i][ibin] = sf.first;
      cout << "sf" << systematics[i] << " "<< sf.first << endl;
     } 
   
     double errdosq = 0;
     double errupsq = 0;
     for (unsigned int i = 0; i < 14; ++i){
       //if (i%2 == 0)
       if ((sfvars[i] - sfcentral.first) < 0.)
         errdosq += (sfvars[i] - sfcentral.first)*(sfvars[i] - sfcentral.first);
       else
         errupsq += (sfvars[i] - sfcentral.first)*(sfvars[i] - sfcentral.first);
     }
     for (unsigned int ifit = 0; ifit < 7; ++ifit){
       cf->SetMatrixOption("WRITE");
       cf->ProducePlots(0);
       cf->SetOptimization(OPT_MORPH_SGN_SIGMA);
       cf->SetMorphing(OPTMORPH_CUTOFF,0.5);
       cf->SetData(Form("data/JP__ptbin_%s", thebin));
       cf->SetDataTag(Form("data/JP__ptbin_%s_PNetBDisc%s", thebin, wp));
       cf->SetDataUntag(Form("data/JP__ptbin_%s_PNetBDisc%s_Fail", thebin, wp));
       if (ifit == 0) cf->SetOptimization(OPT_NOCORR);
       if (ifit == 1) cf->SetMorphing(OPTMORPH_CUTOFF,0.25);
       if (ifit == 2) cf->SetMorphing(OPTMORPH_CUTOFF,0.75);
       if (ifit == 3) cf->SetMorphing(OPTMORPH_GEOMETRIC);
       if (ifit == 4) { //use MC JP calib
        cf->SetData(Form("data_inverted/JP__ptbin_%s", thebin));
        cf->SetDataTag(Form("data_inverted/JP__ptbin_%s_PNetBDisc%s", thebin, wp));
        cf->SetDataUntag(Form("data_inverted/JP__ptbin_%s_PNetBDisc%s_Fail", thebin, wp));
       } 
       if (ifit < 5){ 
         if (ifit == 4)
           system("mkdir pics");
         //double  sf = computeSF(cf, false);
         std::pair<double,double> sf = computeSF(cf, false);
         if (ifit == 4){
            system(Form("mv pics pics_inverted_%s", thebin));
            cout << "sf JP inv" << sf.first << endl;
         }   
         sfvars[ifit+14] = sf.first;
         v_sf_syst[ifit+14][ibin] = sf.first;
         errdosq += (sf.first - sfcentral.first)*(sf.first - sfcentral.first);
         errupsq += (sf.first - sfcentral.first)*(sf.first - sfcentral.first);
       } else { // Cb variation
         sfvars[ifit+14] = sfcentral.first;
         TFile fileinraw(filein);
         double nbmc = ((TH1*)fileinraw.Get(Form("mc/JP__ptbin_%s_b", thebin)))->Integral(); 
         double nbmc_jp = nbmc - ((TH1*)fileinraw.Get(Form("mc/JP__ptbin_%s_b_JP0", thebin)))->Integral();
         double nbmctag = ((TH1*)fileinraw.Get(Form("mc/SVmass__ptbin_%s_b_PNetBDisc%s", thebin, wp)))->Integral() ;
         double nbmctag_jp = nbmctag - ((TH1*)fileinraw.Get(Form("mc/JP__ptbin_%s_b_PNetBDisc%s_JP0", thebin, wp)))->Integral();
         cout << "nbmc " << nbmc << " nbmc_jp " << nbmc_jp << " nbmctag " << nbmctag << "nbmctag_jp " << nbmctag_jp << endl;
         double Cb=(nbmctag/nbmctag_jp)*(nbmc_jp/nbmc);
         double Cb_error = fabs(1-Cb)/2;
         cout << "Cb " << Cb << " Cb_error " << Cb_error << endl;
         if (ifit == 5){
           v_sf_syst[ifit+14][ibin] = sfcentral.first*(1+Cb_error);
           errupsq += (Cb_error*sfcentral.first)*(Cb_error*sfcentral.first);
         } 
         if (ifit == 6){
           v_sf_syst[ifit+14][ibin] = sfcentral.first*(1-Cb_error);
           errdosq += (Cb_error*sfcentral.first)*(Cb_error*sfcentral.first);
         }
       }
     }

     cout << "sf " << sfcentral.first << " +"<< sqrt(errupsq) << " -" <<  sqrt(errdosq) << "(syst) +/-" << rmsstat << "(stat)" << endl; 
     v_sf[ibin] = sfcentral.first;
     //v_sfstat[ibin] = sqrt(std::pow(sfcentral.second, 2)+std::pow(rmsstat,2));
     v_sfstat[ibin] = rmsstat;
     v_sftotup[ibin] = sqrt(v_sfstat[ibin]*v_sfstat[ibin] + errupsq);
     v_sftotdo[ibin] = sqrt(v_sfstat[ibin]*v_sfstat[ibin] + errdosq);

     delete cf;
   }
  
   TFile output(Form("PNetBDisc%s.root", wp), "recreate"); 

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
   for (unsigned int i = 0; i < 21; ++i){
     TGraphErrors gsyst(9, v_pt, v_sf_syst[i], v_pterr, 0);
     gsyst.SetNameTitle(Form("variation%s", systematics[i].c_str()), Form("variation%s", systematics[i].c_str()));
     gsyst.Write();
   }


   system(Form("mkdir -p %s/%s", outdir, wp));
   system(Form("mv pics* %s/%s", outdir, wp));
   system(Form("mv PNetBDisc%s.root %s/%s", wp, outdir, wp));

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
