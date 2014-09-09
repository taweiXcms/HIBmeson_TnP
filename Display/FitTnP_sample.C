#include "utilities.h"

#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooPlot.h>
#include <RooGaussian.h>
#include <RooExponential.h>
#include <RooChebychev.h>
#include <RooAddPdf.h>
#include <RooDecay.h>
#include <RooGaussModel.h>
#include <RooProdPdf.h>
#include <RooAddModel.h>
#include <RooCBShape.h>
#include <RooCategory.h>
#include <RooSimultaneous.h>

bool isDataInput = true;
//bool isDataInput = false;

void FitTnP_sample(){

  TFile*finput=new TFile("../../TnPfiles/foutputData.root","read");
  TString plotfolder = "PlotsBinFitData";
  if(!isDataInput){
    finput=new TFile("../../TnPfiles/foutputMC.root","read");
    plotfolder = "PlotsBinFitMC";
  }
  finput->cd();

  using namespace std;
  using namespace RooFit;

  int numCPU = 6;
  bool quiet = false;
  //bool quiet = true;

  const int nMuPtBin = 7;
  const double MuPtBin[nMuPtBin+1] = {0.0,1.5,3.0,4.5,6.0,9.0,20.0,30.0};
  //const int nMuEtaBin = 5;

  double mumulowMuId[nMuPtBin] = {2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6};
  double mumuhighMuId[nMuPtBin] = {3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5};
  double mumulowtrg[nMuPtBin] = {2.6, 2.6, 2.6, 2.6, 2.6, 2.6, 2.6};
  double mumuhightrg[nMuPtBin] = {3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5};
  double mumutrklow[nMuPtBin] = {2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
  double mumutrkhigh[nMuPtBin] = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
  
  //Crystall Ball parameters
  ////for trigger efficiency
  double CB_trg_Pt_mean_min[nMuPtBin] = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
  double CB_trg_Pt_mean_max[nMuPtBin] = {3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2};
  double CB_trg_Pt_alpha_mean[nMuPtBin] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  double CB_trg_Pt_alpha_min[nMuPtBin] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double CB_trg_Pt_alpha_max[nMuPtBin] = {10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
  double CB_trg_Pt_n_mean[nMuPtBin] ={1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  double CB_trg_Pt_n_min[nMuPtBin] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
  double CB_trg_Pt_n_max[nMuPtBin] = {100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0};
  double CB_trg_Pt_sigma[nMuPtBin] = {0.00, 0.005, 0.005, 0.05, 0.05, 0.05, 0.05};
  double CB_trg_Pt_sigma_min[nMuPtBin] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double CB_trg_Pt_sigma_max[nMuPtBin] = {1, 1, 1, 1, 1, 1, 1};

  ////for Muon ID efficiency
  double CB_MuId_Pt_mean_min[nMuPtBin] = {3.07, 3.07, 3.05, 3.0, 3.0, 3.0, 3.0};
  double CB_MuId_Pt_mean_max[nMuPtBin] = {3.12, 3.12, 3.15, 3.2, 3.2, 3.2, 3.2};
  double CB_MuId_Pt_alpha_mean[nMuPtBin] = {4.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  double CB_MuId_Pt_alpha_min[nMuPtBin] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double CB_MuId_Pt_alpha_max[nMuPtBin] = {10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
  double CB_MuId_Pt_n_mean[nMuPtBin] ={10.0, 5.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  double CB_MuId_Pt_n_min[nMuPtBin] = {5.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
  double CB_MuId_Pt_n_max[nMuPtBin] = {50.0, 50.0, 100.0, 200.0, 100.0, 100.0, 100.0};
  double CB_MuId_Pt_sigma[nMuPtBin] = {0.044, 0.001, 0.01, 0.05, 0.05, 0.05, 0.05};
  //double CB_MuId_Pt_sigma_min[nMuPtBin] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  //double CB_MuId_Pt_sigma_max[nMuPtBin] = {0.1, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1};

  //Gaussian parameters
  double G_Pt_min[nMuPtBin] = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0};
  double G_Pt_max[nMuPtBin] = {3.2, 3.2, 3.2, 3.2, 3.2, 3.2, 3.2};
  double G1_Pt_sigma_mean[nMuPtBin] = {0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15};
  double G1_Pt_sigma_min[nMuPtBin] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
  double G1_Pt_sigma_max[nMuPtBin] = {0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25};
  double G2_Pt_sigma_mean[nMuPtBin] = {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02};
  double G2_Pt_sigma_min[nMuPtBin] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
  double G2_Pt_sigma_max[nMuPtBin] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
  double G_Pt_frac_mean[nMuPtBin] = {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2};
  double G_Pt_frac_min[nMuPtBin] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double G_Pt_frac_max[nMuPtBin] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  //Chebychev parameters
  ////for Muon ID efficiency
  double Cheb_MuId_Pt_p1_Pass_mean[nMuPtBin] = {-0.80, 0., 0., 0., 0., 0., 0.};
  double Cheb_MuId_Pt_p1_Pass_min[nMuPtBin] = {-1., -1., -1., -1., -1., -1., -1.};
  double Cheb_MuId_Pt_p1_Pass_max[nMuPtBin] = {+1., +1., +1., +1., +1., +1., +1.};
  double Cheb_MuId_Pt_p1_Fail_mean[nMuPtBin] = {0.05, 0., 0., 0., 0., 0., 0.};
  double Cheb_MuId_Pt_p1_Fail_min[nMuPtBin] = {-1., -1., -1., -1., -1., -1., -1.};
  double Cheb_MuId_Pt_p1_Fail_max[nMuPtBin] = {+1., +1., +1., +1., +1., +1., +1.};

  double Cheb_MuId_Pt_p2_Pass_mean[nMuPtBin] = {-0.002, 0., 0., 0., 0., 0., 0.};
  double Cheb_MuId_Pt_p2_Pass_min[nMuPtBin] = {-1., -1., -1., -1., -1., -1., -1.};
  double Cheb_MuId_Pt_p2_Pass_max[nMuPtBin] = {+1., +1., +1., +1., +1., +1., +1.};
  double Cheb_MuId_Pt_p2_Fail_mean[nMuPtBin] = {0.002, 0., 0., 0., 0., 0., 0.};
  double Cheb_MuId_Pt_p2_Fail_min[nMuPtBin] = {-1., -1., -1., -1., -1., -1., -1.};
  double Cheb_MuId_Pt_p2_Fail_max[nMuPtBin] = {+1., +1., +1., +1., +1., +1., +1.};

  ////for tracking efficinecy
  double Cheb_trk_Pt_p1_Pass_mean[nMuPtBin] = {0., 0., 0., 0., 0., 0., 0.};
  double Cheb_trk_Pt_p1_Pass_min[nMuPtBin] = {-1., -1., -1., -1., -1., -1., -1.};
  double Cheb_trk_Pt_p1_Pass_max[nMuPtBin] = {+1., +1., +1., +1., +1., +1., +1.};
  double Cheb_trk_Pt_p1_Fail_mean[nMuPtBin] = {0., 0., 0., 0., 0., 0., 0.};
  double Cheb_trk_Pt_p1_Fail_min[nMuPtBin] = {-1., -1., -1., -1., -1., -1., -1.};
  double Cheb_trk_Pt_p1_Fail_max[nMuPtBin] = {+1., +1., +1., +1., +1., +1., +1.};

  double Cheb_trk_Pt_p2_Pass_mean[nMuPtBin] = {0., 0., 0., 0., 0., 0., 0.};
  double Cheb_trk_Pt_p2_Pass_min[nMuPtBin] = {-1., -1., -1., -1., -1., -1., -1.};
  double Cheb_trk_Pt_p2_Pass_max[nMuPtBin] = {+1., +1., +1., +1., +1., +1., +1.};
  double Cheb_trk_Pt_p2_Fail_mean[nMuPtBin] = {0., 0., 0., 0., 0., 0., 0.};
  double Cheb_trk_Pt_p2_Fail_min[nMuPtBin] = {-1., -1., -1., -1., -1., -1., -1.};
  double Cheb_trk_Pt_p2_Fail_max[nMuPtBin] = {+1., +1., +1., +1., +1., +1., +1.};

  //exponential parameters
  double lp_Pt_mean[nMuPtBin] = {-5., 0., 0., 0., 0., 0., 0.};
  double lp_Pt_min[nMuPtBin] = {-10., -10., -5., -5., -5., -5., -5.};
  double lp_Pt_max[nMuPtBin] = {10., 10., 5., 5., 5., 5., 5.};
  double lf_Pt_mean[nMuPtBin] = {0., 0., 0., 0., 0., 0., 0.};
  double lf_Pt_min[nMuPtBin] = {-10., -10., -5., -5., -5., -5., -5.};
  double lf_Pt_max[nMuPtBin] = {10., 10., 5., 5., 5., 5., 5.};

  //make canvas
  ////tracking efficiency
  TCanvas* ctrkpass = new TCanvas("ctrkpass", "", 0, 20, 600, 600);
  TCanvas* ctrkfail = new TCanvas("ctrkfail", "", 0, 40, 600, 600);
  TCanvas* ctrkall = new TCanvas("ctrkall", "", 0, 60, 600, 600);
  ////Muon ID efficiency
  TCanvas* cMuIdpass = new TCanvas("cMuIdpass", "", 10, 20, 600, 600);
  TCanvas* cMuIdfail = new TCanvas("cMuIdfail", "", 10, 40, 600, 600);
  TCanvas* cMuIdall = new TCanvas("cMuIdall", "", 10, 60, 600, 600);
  ////trigger efficiency
  TCanvas* ctrgpass = new TCanvas("ctrgpass", "", 20, 20, 600, 600);
  TCanvas* ctrgfail = new TCanvas("ctrgfail", "", 20, 40, 600, 600);
  TCanvas* ctrgall = new TCanvas("ctrgall", "", 20, 60, 600, 600);

  //define pt histogram
  ////tracking efficiency
  TH1D* hTrkPtPass[nMuPtBin];
  TH1D* hTrkPtFail[nMuPtBin];
  ////Muon Id efficinecy
  TH1D* hMuIdPtPass[nMuPtBin];
  TH1D* hMuIdPtFail[nMuPtBin];
  ////trigger efficiency
  TH1D* hTrgPtPass[nMuPtBin];
  TH1D* hTrgPtFail[nMuPtBin];
/*
  TH1D* hTrkEtaPass[nMuEtaBin];
  TH1D* hTrkEtaFail[nMuEtaBin];
  TH1D* hTrkEtaAll[nMuEtaBin];
*/

  TH1D* EffTrig = myTH1D("EffTrig", "Mu Pt (GeV)", "efficiency", 2, nMuPtBin, MuPtBin);
  TH1D* EffTrk = myTH1D("EffTrk", "Mu Pt (GeV)", "efficiency", 3, nMuPtBin, MuPtBin);
  TH1D* EffMuId = myTH1D("EffMuId", "Mu Pt (GeV)", "efficiency", 4, nMuPtBin, MuPtBin);

  //make categrory
  ////tracking efficiency
  RooCategory sampleTrk("sampleTrk", "sampleTrk");
  sampleTrk.defineType("Pass");
  sampleTrk.defineType("Fail");
  ////Muon ID efficiency
  RooCategory sampleMuId("sampleMuId", "sampleMuId");
  sampleMuId.defineType("Pass");
  sampleMuId.defineType("Fail");
  ////trigger efficiency
  RooCategory sampleTrg("sampleTrg", "sampleTrg");
  sampleTrg.defineType("Pass");
  sampleTrg.defineType("Fail");

  //get pt histogram
  for(int i = 0; i < nMuPtBin; i++)
//  for(int i = 1; i < 2; i++)
  {
    RooRealVar masstrk("masstrk", "masstrk", mumutrklow[i], mumutrkhigh[i]);
    masstrk.setBins(100);
    RooRealVar massMuId("massMuId", "massMuId", mumulowMuId[i], mumuhighMuId[i]);
    massMuId.setBins(100);
    RooRealVar masstrg("masstrg", "masstrg", mumulowtrg[i], mumuhightrg[i]);
    masstrg.setBins(100);

    //tracking efficiency
    hTrkPtPass[i] = (TH1D*)finput->Get(Form("hTrkPtPass%d", i));
    hTrkPtFail[i] = (TH1D*)finput->Get(Form("hTrkPtFail%d", i));
/*
    double pbinnumtrk = hTrkPtPass[i]->GetNbinsX();
    double plowbintrk = hTrkPtPass[i]->GetBinLowEdge(1);
    double phighbintrk = hTrkPtPass[i]->GetBinLowEdge(100) + hTrkPtPass[i]->GetBinWidth(100);
    cout << "Passtrk number of bin: " << pbinnumtrk << endl;
    cout << "Passtrk low bin: " << plowbintrk << endl;
    cout << "Passtrk high bin: " << phighbintrk << endl;

    double fbinnumtrk = hTrkPtFail[i]->GetNbinsX();
    double flowbintrk = hTrkPtFail[i]->GetBinLowEdge(1);
    double fhighbintrk = hTrkPtFail[i]->GetBinLowEdge(100) + hTrkPtFail[i]->GetBinWidth(100);
    cout << "Failtrk number of bin: " << fbinnumtrk << endl;
    cout << "Failtrk low bin: " << flowbintrk << endl;
    cout << "Failtrk high bin: " << fhighbintrk << endl;
*/
    //Muon Id efficiency
    hMuIdPtPass[i] = (TH1D*)finput->Get(Form("hMuIdPtPass%d", i));
    hMuIdPtFail[i] = (TH1D*)finput->Get(Form("hMuIdPtFail%d", i));
/*
    double pbinnumMuId = hMuIdPtPass[i]->GetNbinsX();
    double plowbinMuId = hMuIdPtPass[i]->GetBinLowEdge(0);
    double phighbinMuId = hMuIdPtPass[i]->GetBinLowEdge(99);
    cout << "PassMuId number of bin: " << pbinnumMuId << endl;
    cout << "PassMuId low bin: " << plowbinMuId << endl;
    cout << "PassMuId high bin: " << phighbinMuId << endl;

    double fbinnumMuId = hMuIdPtFail[i]->GetNbinsX();
    double flowbinMuId = hMuIdPtFail[i]->GetBinLowEdge(0);
    double fhighbinMuId = hMuIdPtFail[i]->GetBinLowEdge(99);
    cout << "FailMuId number of bin: " << fbinnumMuId << endl;
    cout << "FailMuId low bin: " << flowbinMuId << endl;
    cout << "FailMuId high bin: " << fhighbinMuId << endl;
*/
    //trigger efficiency
    hTrgPtPass[i] = (TH1D*)finput->Get(Form("hTrigPtPass%d", i));
    hTrgPtFail[i] = (TH1D*)finput->Get(Form("hTrigPtFail%d", i));
/*
    double pbinnumtrg = hTrgPtPass[i]->GetNbinsX();
    double plowbintrg = hTrgPtPass[i]->GetBinLowEdge(0);
    double phighbintrg = hTrgPtPass[i]->GetBinLowEdge(99);
    cout << "Passtrg number of bin: " << pbinnumtrg << endl;
    cout << "Passtrg low bin: " << plowbintrg << endl;
    cout << "Passtrg high bin: " << phighbintrg << endl;

    double fbinnumtrg = hTrkPtFail[i]->GetNbinsX();
    double flowbintrg = hTrkPtFail[i]->GetBinLowEdge(0);
    double fhighbintrg = hTrkPtFail[i]->GetBinLowEdge(99);
    cout << "Failtrg number of bin: " << fbinnumtrg << endl;
    cout << "Failtrg low bin: " << flowbintrg << endl;
    cout << "Failtrg high bin: " << fhighbintrg << endl;
*/

    //make CB function for trigger efficiency
    RooRealVar mean_CB_trg("mean_CB_trg", "mean_CB_trg", 3.1, CB_trg_Pt_mean_min[i], CB_trg_Pt_mean_max[i]);
    //RooRealVar sigma_CB_trg("sigma_CB_trg", "sigma_CB_trg", CB_trg_Pt_sigma[i], CB_trg_Pt_sigma_min[i], CB_trg_Pt_sigma_max[i]);
    RooRealVar sigma_CB_trg("sigma_CB_trg", "sigma_CB_trg", CB_trg_Pt_sigma[i]);
    sigma_CB_trg.setConstant(kFALSE);
    RooRealVar alpha_CB_trg("alpha_CB_trg", "alpha_CB_trg", CB_trg_Pt_alpha_mean[i], CB_trg_Pt_alpha_min[i], CB_trg_Pt_alpha_max[i]);
    RooRealVar n_CB_trg("n_CB_trg", "n_CB_trg", CB_trg_Pt_n_mean[i], CB_trg_Pt_n_min[i], CB_trg_Pt_n_max[i]);

    RooCBShape signal_CB_trg_Pass("signal_CB_trg_Pass", "signal_CB_trg_Pass", masstrg, mean_CB_trg, sigma_CB_trg, alpha_CB_trg, n_CB_trg);
    RooCBShape signal_CB_trg_Fail("signal_CB_trg_Fail", "signal_CB_trg_Fail", masstrg, mean_CB_trg, sigma_CB_trg, alpha_CB_trg, n_CB_trg);

    //make CB function for Muon ID efficiency
    RooRealVar mean_CB_MuId("mean_CB_MuId", "mean_CB_MuId", 3.1, CB_MuId_Pt_mean_min[i], CB_MuId_Pt_mean_max[i]);
    //RooRealVar sigma_CB_MuId("sigma_CB_MuId", "sigma_CB_MuId", CB_MuId_Pt_sigma[i], CB_MuId_Pt_sigma_min[i], CB_MuId_Pt_sigma_max[i]);
    RooRealVar sigma_CB_MuId("sigma_CB_MuId", "sigma_CB_MuId", CB_MuId_Pt_sigma[i]);
    sigma_CB_MuId.setConstant(kFALSE);
    RooRealVar alpha_CB_MuId("alpha_CB_MuId", "alpha_CB_MuId", CB_MuId_Pt_alpha_mean[i], CB_MuId_Pt_alpha_min[i], CB_MuId_Pt_alpha_max[i]);
    RooRealVar n_CB_MuId("n_CB_MuId", "n_CB_MuId", CB_MuId_Pt_n_mean[i], CB_MuId_Pt_n_min[i], CB_MuId_Pt_n_max[i]);

    RooCBShape signal_CB_MuId_Pass("signal_CB_MuId_Pass", "signal_CB_MuId_Pass", massMuId, mean_CB_MuId, sigma_CB_MuId, alpha_CB_MuId, n_CB_MuId);
    RooCBShape signal_CB_MuId_Fail("signal_CB_MuId_Fail", "signal_CB_MuId_Fail", massMuId, mean_CB_MuId, sigma_CB_MuId, alpha_CB_MuId, n_CB_MuId);

    //make 2gauss function
    RooRealVar mean_gaus("mean_gaus", "mean_gaus", 3.1, G_Pt_min[i], G_Pt_max[i]);
    RooRealVar sigma1_gaus("sigma1_gaus", "sigma1_gaus", G1_Pt_sigma_mean[i], G1_Pt_sigma_min[i], G1_Pt_sigma_max[i]);
    RooRealVar sigma2_gaus("sigma2_gaus", "sigma2_gaus", G2_Pt_sigma_mean[i], G2_Pt_sigma_min[i], G2_Pt_sigma_max[i]);
    RooGaussian gauss1("gauss1", "gauss1", masstrk, mean_gaus, sigma1_gaus);
    RooGaussian gauss2("gauss2", "gauss2", masstrk, mean_gaus, sigma2_gaus);
    RooRealVar mfrac_gaus("mfrac_gaus", "mfrac_gaus", G_Pt_frac_mean[i], G_Pt_frac_min[i], G_Pt_frac_max[i]);

    RooAddPdf signal_gauss_Pass("signal_gauss_Pass", "signal_gauss_Pass", RooArgList(gauss1, gauss2), RooArgList(mfrac_gaus));
    RooAddPdf signal_gauss_Fail("signal_gauss_Fail", "signal_gauss_Fail", RooArgList(gauss1, gauss2), RooArgList(mfrac_gaus));

    //make Chebyschev function for Muon ID efficiency
    RooRealVar cheb_MuId_p1Pass("cheb_MuId_p1Pass", "cheb_MuId_p1Pass", Cheb_MuId_Pt_p1_Pass_mean[i], Cheb_MuId_Pt_p1_Pass_min[i], Cheb_MuId_Pt_p1_Pass_max[i]);
    RooRealVar cheb_MuId_p2Pass("cheb_MuId_p2Pass", "cheb_MuId_p2Pass", Cheb_MuId_Pt_p2_Pass_mean[i], Cheb_MuId_Pt_p2_Pass_min[i], Cheb_MuId_Pt_p2_Pass_max[i]);
    RooChebychev background_cheb_MuId_Pass("background_cheb_MuId_Pass", "background_cheb_MuId_Pass", massMuId, RooArgList(cheb_MuId_p1Pass, cheb_MuId_p2Pass));
    RooRealVar cheb_MuId_p1Fail("cheb_MuId_p1Fail", "cheb_MuId_p1Fail", Cheb_MuId_Pt_p1_Fail_mean[i], Cheb_MuId_Pt_p1_Fail_min[i], Cheb_MuId_Pt_p1_Fail_max[i]);
    RooRealVar cheb_MuId_p2Fail("cheb_MuId_p2Fail", "cheb_MuId_p2Fail", Cheb_MuId_Pt_p2_Fail_mean[i], Cheb_MuId_Pt_p2_Fail_min[i], Cheb_MuId_Pt_p2_Fail_max[i]);
    RooChebychev background_cheb_MuId_Fail("background_cheb_MuId_Fail", "background_cheb_MuId_Fail", massMuId, RooArgList(cheb_MuId_p1Fail, cheb_MuId_p2Fail));

    //make Chebyschev function for tracking efficiency
    RooRealVar cheb_trk_p1Pass("cheb_trk_p1Pass", "cheb_trk_p1Pass", Cheb_trk_Pt_p1_Pass_mean[i], Cheb_trk_Pt_p1_Pass_min[i], Cheb_trk_Pt_p1_Pass_max[i]);
    RooRealVar cheb_trk_p2Pass("cheb_trk_p2Pass", "cheb_trk_p2Pass", Cheb_trk_Pt_p2_Pass_mean[i], Cheb_trk_Pt_p2_Pass_min[i], Cheb_trk_Pt_p2_Pass_max[i]);
    RooChebychev background_cheb_trk_Pass("background_cheb_trk_Pass", "background_cheb_trk_Pass", masstrk, RooArgList(cheb_trk_p1Pass, cheb_trk_p2Pass));
    RooRealVar cheb_trk_p1Fail("cheb_trk_p1Fail", "cheb_trk_p1Fail", Cheb_trk_Pt_p1_Fail_mean[i], Cheb_trk_Pt_p1_Fail_min[i], Cheb_trk_Pt_p1_Fail_max[i]);
    RooRealVar cheb_trk_p2Fail("cheb_trk_p2Fail", "cheb_trk_p2Fail", Cheb_trk_Pt_p2_Fail_mean[i], Cheb_trk_Pt_p2_Fail_min[i], Cheb_trk_Pt_p2_Fail_max[i]);
    RooChebychev background_cheb_trk_Fail("background_cheb_trk_Fail", "background_cheb_trk_Fail", masstrk, RooArgList(cheb_trk_p1Fail, cheb_trk_p2Fail));

    //make exponential function
    RooRealVar lp("lp", "lp", lp_Pt_mean[i], lp_Pt_min[i], lp_Pt_max[i]);
    RooRealVar lf("lf", "lf", lf_Pt_mean[i], lf_Pt_min[i], lf_Pt_max[i]);
    RooExponential background_exp_Pass("background_exp_Pass", "background_exp_Pass", masstrg, lp);
    RooExponential background_exp_Fail("background_exp_Fail", "background_exp_Fail", masstrg, lf);

    //set axis variable
    RooRealVar weight("weight", "weight", 1., 0, 1E9);

    //put hist to dataset
    ////tracking efficinecy
    RooDataHist dataTrkPass("dataTrkPass", "dataTrkPass", masstrk, hTrkPtPass[i]);
    RooDataHist dataTrkFail("dataTrkFail", "dataTrkFail", masstrk, hTrkPtFail[i]);

    RooDataSet dataTrkPassset("dataTrkPassset", "masstrk", RooArgSet(masstrk, weight), WeightVar("weight"));
    RooDataSet dataTrkFailset("dataTrkFailset", "masstrk", RooArgSet(masstrk, weight), WeightVar("weight"));

    ////Muon ID efficiency
    RooDataHist dataMuIdPass("dataMuIdPass", "dataMuIdPass", massMuId, hMuIdPtPass[i]);
    RooDataHist dataMuIdFail("dataMuIdFail", "dataMuIdFail", massMuId, hMuIdPtFail[i]);

    RooDataSet dataMuIdPassset("dataMuIdPassset", "massMuId", RooArgSet(massMuId, weight), WeightVar("weight"));
    RooDataSet dataMuIdFailset("dataMuIdFailset", "massMuId", RooArgSet(massMuId, weight), WeightVar("weight"));

    ////trigger efficiency
    RooDataHist dataTrgPass("dataTrgPass", "dataTrgPass", masstrg, hTrgPtPass[i]);
    RooDataHist dataTrgFail("dataTrgFail", "dataTrgFail", masstrg, hTrgPtFail[i]);

    RooDataSet dataTrgPassset("dataTrgPassset", "masstrg", RooArgSet(masstrg, weight), WeightVar("weight"));
    RooDataSet dataTrgFailset("dataTrgFailset", "masstrg", RooArgSet(masstrg, weight), WeightVar("weight"));

    //get entries
    for(int j = 0; j < dataTrkPass.numEntries(); j++)
    {
      ////tracking efficinecy
      dataTrkPassset.add(*dataTrkPass.get(j), dataTrkPass.weight());
      dataTrkFailset.add(*dataTrkFail.get(j), dataTrkFail.weight());
    }    
    for(int j = 0; j < dataMuIdPass.numEntries(); j++)
    {
      ////Muon ID efficinecy
      dataMuIdPassset.add(*dataMuIdPass.get(j), dataMuIdPass.weight());
      dataMuIdFailset.add(*dataMuIdFail.get(j), dataMuIdFail.weight());
    }
    for(int j = 0; j < dataTrgPass.numEntries(); j++)
    {
      ////trigger efficinecy
      dataTrgPassset.add(*dataTrgPass.get(j), dataTrgPass.weight());
      dataTrgFailset.add(*dataTrgFail.get(j), dataTrgFail.weight());
    }

    //make data set by adding passing and failing
    RooDataSet dataTrkAllset("dataTrkAllset", "masstrk", RooArgSet(masstrk, weight), Index(sampleTrk), Import("Pass", dataTrkPassset), Import("Fail", dataTrkFailset), WeightVar(weight));
    RooDataSet dataMuIdAllset("dataMuIdAllset", "massMuId", RooArgSet(massMuId, weight), Index(sampleMuId), Import("Pass", dataMuIdPassset), Import("Fail", dataMuIdFailset), WeightVar(weight));
    RooDataSet dataTrgAllset("dataTrgAllset", "masstrg", RooArgSet(masstrg, weight), Index(sampleTrg), Import("Pass", dataTrgPassset), Import("Fail", dataTrgFailset), WeightVar(weight));

    //set efficiency and number of signal
    //tracking efficiency
    RooRealVar nsigtrk("nsigtrk", "nsigtrk", 1, 0, 1E9);
    RooRealVar efficiencytrk("efficiencytrk", "efficiencytrk", 0.9, 0., 1.);
    RooFormulaVar nsigPasstrk("nsigPasstrk", "nsigPasstrk", "nsigtrk*efficiencytrk", RooArgList(nsigtrk, efficiencytrk));
    RooFormulaVar nsigFailtrk("nsigFailtrk", "nsigFailtrk", "nsigtrk*(1-efficiencytrk)", RooArgList(nsigtrk, efficiencytrk));

    RooRealVar nbkgPasstrk("nbkgPasstrk", "nbkgPasstrk", 1, 0, 1E9);
    RooRealVar nbkgFailtrk("nbkgFailtrk", "nbkgFailtrk", 1, 0, 1E9);

    //Muon ID efficiency
    RooRealVar nsigMuId("nsigMuId", "nsigMuId", 1, 0, 1E9);
    RooRealVar efficiencyMuId("efficiencyMuId", "efficiencyMuId", 0.9, 0., 1.);
    RooFormulaVar nsigPassMuId("nsigPassMuId", "nsigPassMuId", "nsigMuId*efficiencyMuId", RooArgList(nsigMuId, efficiencyMuId));
    RooFormulaVar nsigFailMuId("nsigFailMuId", "nsigFailMuId", "nsigMuId*(1-efficiencyMuId)", RooArgList(nsigMuId, efficiencyMuId));

    RooRealVar nbkgPassMuId("nbkgPassMuId", "nbkgPassMuId", 1, 0, 1E9);
    RooRealVar nbkgFailMuId("nbkgFailMuId", "nbkgFailMuId", 1, 0, 1E9);

    //tigger efficiency
    RooRealVar nsigtrg("nsigtrg", "nsigtrg", 1, 0, 1E9);
    RooRealVar efficiencytrg("efficiencytrg", "efficiencytrg", 0.9, 0., 1.);
    RooFormulaVar nsigPasstrg("nsigPasstrg", "nsigPasstrg", "nsigtrg*efficiencytrg", RooArgList(nsigtrg, efficiencytrg));
    RooFormulaVar nsigFailtrg("nsigFailtrg", "nsigFailtrg", "nsigtrg*(1-efficiencytrg)", RooArgList(nsigtrg, efficiencytrg));

    RooRealVar nbkgPasstrg("nbkgPasstrg", "nbkgPasstrg", 1, 0, 1E9);
    RooRealVar nbkgFailtrg("nbkgFailtrg", "nbkgFailtrg", 1, 0, 1E9);

    //define PDFs

    ////CB plus Chebyschev
    RooAddPdf model_CB_cheb_Pass("model_CB_cheb_Pass", "model_CB_cheb_Pass", RooArgList(signal_CB_MuId_Pass, background_cheb_MuId_Pass), RooArgList(nsigPassMuId, nbkgPassMuId));
    RooAddPdf model_CB_cheb_Fail("model_CB_cheb_Fail", "model_CB_cheb_Fail", RooArgList(signal_CB_MuId_Fail, background_cheb_MuId_Fail), RooArgList(nsigFailMuId, nbkgFailMuId));

    ////CB plus exponential
    RooAddPdf model_CB_exp_Pass("model_CB_exp_Pass", "model_CB_exp_Pass", RooArgList(signal_CB_trg_Pass, background_exp_Pass), RooArgList(nsigPasstrg, nbkgPasstrg));
    RooAddPdf model_CB_exp_Fail("model_CB_exp_Fail", "model_CB_exp_Fail", RooArgList(signal_CB_trg_Fail, background_exp_Fail), RooArgList(nsigFailtrg, nbkgFailtrg));

    ////Gaussian plus Chebyschev
    RooAddPdf model_gaus_cheb_Pass("model_gaus_cheb_Pass", "model_gaus_cheb_Pass", RooArgList(signal_gauss_Pass, background_cheb_trk_Pass), RooArgList(nsigPasstrk, nbkgPasstrk));
    RooAddPdf model_gaus_cheb_Fail("model_gaus_cheb_Fail", "model_gaus_cheb_Fail", RooArgList(signal_gauss_Fail, background_cheb_trk_Fail), RooArgList(nsigFailtrk, nbkgFailtrk));

    //fit
    ////tracking efficiency
    RooSimultaneous simTrkPdf("simTrkPdf", "simTrkPdf", sampleTrk);
    simTrkPdf.addPdf(model_gaus_cheb_Pass, "Pass");
    simTrkPdf.addPdf(model_gaus_cheb_Fail, "Fail");
    ////Muon ID efficiency
    RooSimultaneous simMuIdPdf("simMuIdPdf", "simMuIdPdf", sampleMuId);
    simMuIdPdf.addPdf(model_CB_cheb_Pass, "Pass");
    simMuIdPdf.addPdf(model_CB_cheb_Fail, "Fail");
    ////trigger efficiency
    RooSimultaneous simTrgPdf("simTrgPdf", "simTrgPdf", sampleTrg);
    simTrgPdf.addPdf(model_CB_exp_Pass, "Pass");
    simTrgPdf.addPdf(model_CB_exp_Fail, "Fail");

//    simTrkPdf.fitTo(dataTrkAllset, Extended(true), NumCPU(numCPU), PrintLevel(quiet?-1:1), Warnings(!quiet));
//    simMuIdPdf.fitTo(dataMuIdAllset, Extended(true), NumCPU(numCPU), PrintLevel(quiet?-1:1), Warnings(!quiet));
    simTrgPdf.fitTo(dataTrgAllset, Extended(true), NumCPU(numCPU), PrintLevel(quiet?-1:1), Warnings(!quiet));

    //draw tracking efficiency mass distribution
    ////pass
    RooPlot *frameTrkPass = masstrk.frame();
    dataTrkAllset.plotOn(frameTrkPass, Cut("sampleTrk==sampleTrk::Pass"));
    simTrkPdf.plotOn(frameTrkPass, Slice(sampleTrk, "Pass"), ProjWData(sampleTrk,dataTrkAllset),LineStyle(kDashed), LineColor(kBlue));
    simTrkPdf.plotOn(frameTrkPass, Slice(sampleTrk, "Pass"), Components("signal_gauss_Pass"),ProjWData(sampleTrk,dataTrkAllset),LineStyle(kDashed), LineColor(kGreen));
    simTrkPdf.plotOn(frameTrkPass, Slice(sampleTrk, "Pass"), Components("background_cheb_trk_Pass"),ProjWData(sampleTrk,dataTrkAllset),LineStyle(kDashed), LineColor(kRed));

    ctrkpass->cd();
    frameTrkPass->Draw();
    ctrkpass->SaveAs(Form("%s/Trk_Pt_bin%d_Pass.pdf", plotfolder.Data(), i));

    ////fail
    RooPlot *frameTrkFail = masstrk.frame();
    dataTrkAllset.plotOn(frameTrkFail, Cut("sampleTrk==sampleTrk::Fail"));

    simTrkPdf.plotOn(frameTrkFail, Slice(sampleTrk, "Fail"), ProjWData(sampleTrk,dataTrkAllset),LineStyle(kDashed), LineColor(kBlue));
    simTrkPdf.plotOn(frameTrkFail, Slice(sampleTrk, "Fail"), Components("signal_gauss_Fail"),ProjWData(sampleTrk,dataTrkAllset),LineStyle(kDashed), LineColor(kGreen));
    simTrkPdf.plotOn(frameTrkFail, Slice(sampleTrk, "Fail"), Components("background_cheb_trk_Fail"),ProjWData(sampleTrk,dataTrkAllset),LineStyle(kDashed), LineColor(kRed));

    ctrkfail->cd();
    frameTrkFail->Draw();
    ctrkfail->SaveAs(Form("%s/Trk_Pt_bin%d_Fail.pdf", plotfolder.Data(), i));

    ////all
    RooPlot *frameTrkAll = masstrk.frame();
    dataTrkAllset.plotOn(frameTrkAll);

    simTrkPdf.plotOn(frameTrkAll, ProjWData(sampleTrk, dataTrkAllset), LineColor(kBlue));
    simTrkPdf.plotOn(frameTrkAll, Components("signal_gauss_Pass, signal_gauss_Fail"), ProjWData(sampleTrk,dataTrkAllset),LineStyle(kDashed), LineColor(kGreen));
    simTrkPdf.plotOn(frameTrkAll, Components("background_cheb_trk_Pass, background_cheb_trk_Fail"), ProjWData(sampleTrk,dataTrkAllset),LineStyle(kDashed), LineColor(kRed));

    ctrkall->cd();
    frameTrkAll->Draw();
    ctrkall->SaveAs(Form("%s/Trk_Pt_bin%d_All.pdf", plotfolder.Data(), i));

    //draw MuonId efficiency mass distribution
    ////pass
    RooPlot *frameMuIdPass = massMuId.frame();
    dataMuIdAllset.plotOn(frameMuIdPass, Cut("sampleMuId==sampleMuId::Pass"));

    simMuIdPdf.plotOn(frameMuIdPass, Slice(sampleMuId, "Pass"), ProjWData(sampleMuId,dataMuIdAllset),LineStyle(kDashed), LineColor(kBlue));
    simMuIdPdf.plotOn(frameMuIdPass, Slice(sampleMuId, "Pass"), Components("signal_CB_MuId_Pass"),ProjWData(sampleMuId,dataMuIdAllset),LineStyle(kDashed), LineColor(kGreen));
    simMuIdPdf.plotOn(frameMuIdPass, Slice(sampleMuId, "Pass"), Components("background_cheb_MuId_Pass"),ProjWData(sampleMuId,dataMuIdAllset),LineStyle(kDashed), LineColor(kRed));

    cMuIdpass->cd();
    frameMuIdPass->Draw();
    cMuIdpass->SaveAs(Form("%s/MuId_Pt_bin%d_Pass.pdf", plotfolder.Data(), i));

    ////fail
    RooPlot *frameMuIdFail = massMuId.frame();
    dataMuIdAllset.plotOn(frameMuIdFail, Cut("sampleMuId==sampleMuId::Fail"));

    simMuIdPdf.plotOn(frameMuIdFail, Slice(sampleMuId, "Fail"), ProjWData(sampleMuId,dataMuIdAllset),LineStyle(kDashed), LineColor(kBlue));
    simMuIdPdf.plotOn(frameMuIdFail, Slice(sampleMuId, "Fail"), Components("signal_CB_MuId_Fail"),ProjWData(sampleMuId,dataMuIdAllset),LineStyle(kDashed), LineColor(kGreen));
    simMuIdPdf.plotOn(frameMuIdFail, Slice(sampleMuId, "Fail"), Components("background_cheb_MuId_Fail"),ProjWData(sampleMuId,dataMuIdAllset),LineStyle(kDashed), LineColor(kRed));

    cMuIdfail->cd();
    frameMuIdFail->Draw();
    cMuIdfail->SaveAs(Form("%s/MuId_Pt_bin%d_Fail.pdf", plotfolder.Data(), i));

    ////all
    RooPlot *frameMuIdAll = massMuId.frame();
    dataMuIdAllset.plotOn(frameMuIdAll);

    simMuIdPdf.plotOn(frameMuIdAll, ProjWData(sampleMuId, dataMuIdAllset), LineColor(kBlue));
    simMuIdPdf.plotOn(frameMuIdAll, Components("signal_CB_MuId_Pass, signal_CB_MuId_Fail"), ProjWData(sampleMuId,dataMuIdAllset),LineStyle(kDashed), LineColor(kGreen));
    simMuIdPdf.plotOn(frameMuIdAll, Components("background_cheb_MuId_Pass, background_cheb_MuId_Fail"), ProjWData(sampleMuId,dataMuIdAllset),LineStyle(kDashed), LineColor(kRed));

    cMuIdall->cd();
    frameMuIdAll->Draw();
    cMuIdall->SaveAs(Form("%s/MuId_Pt_bin%d_All.pdf", plotfolder.Data(), i));

    //draw trigger efficiency mass distribution
    ////pass
    RooPlot *frameTrgPass = masstrg.frame();
    dataTrgAllset.plotOn(frameTrgPass, Cut("sampleTrg==sampleTrg::Pass"));

    simTrgPdf.plotOn(frameTrgPass, Slice(sampleTrg, "Pass"), ProjWData(sampleTrg,dataTrgAllset),LineStyle(kDashed), LineColor(kBlue));
    simTrgPdf.plotOn(frameTrgPass, Slice(sampleTrg, "Pass"), Components("signal_CB_trg_Pass"),ProjWData(sampleTrg,dataTrgAllset),LineStyle(kDashed), LineColor(kGreen));
    simTrgPdf.plotOn(frameTrgPass, Slice(sampleTrg, "Pass"), Components("background_exp_Pass"),ProjWData(sampleTrg,dataTrgAllset),LineStyle(kDashed), LineColor(kRed));

    ctrgpass->cd();
    frameTrgPass->Draw();
    ctrgpass->SaveAs(Form("%s/Trg_Pt_bin%d_Pass.pdf", plotfolder.Data(), i));

    ////fail
    RooPlot *frameTrgFail = masstrg.frame();
    dataTrgAllset.plotOn(frameTrgFail, Cut("sampleTrg==sampleTrg::Fail"));

    simTrgPdf.plotOn(frameTrgFail, Slice(sampleTrg, "Fail"), ProjWData(sampleTrg,dataTrgAllset),LineStyle(kDashed), LineColor(kBlue));
    simTrgPdf.plotOn(frameTrgFail, Slice(sampleTrg, "Fail"), Components("signal_CB_trg_Fail"),ProjWData(sampleTrg,dataTrgAllset),LineStyle(kDashed), LineColor(kGreen));
    simTrgPdf.plotOn(frameTrgFail, Slice(sampleTrg, "Fail"), Components("background_exp_Fail"),ProjWData(sampleTrg,dataTrgAllset),LineStyle(kDashed), LineColor(kRed));

    ctrgfail->cd();
    frameTrgFail->Draw();
    ctrgfail->SaveAs(Form("%s/Trg_Pt_bin%d_Fail.pdf", plotfolder.Data(), i));

    ////all
    RooPlot *frameTrgAll = masstrg.frame();
    dataTrgAllset.plotOn(frameTrgAll);

    simTrgPdf.plotOn(frameTrgAll, ProjWData(sampleTrg, dataTrgAllset), LineColor(kBlue));
    simTrgPdf.plotOn(frameTrgAll, Components("signal_CB_trg_Pass, signal_CB_trg_Fail"), ProjWData(sampleTrg,dataTrgAllset),LineStyle(kDashed), LineColor(kGreen));
    simTrgPdf.plotOn(frameTrgAll, Components("background_exp_Pass, background_exp_Fail"), ProjWData(sampleTrg,dataTrgAllset),LineStyle(kDashed), LineColor(kRed));

    ctrgall->cd();
    frameTrgAll->Draw();
    ctrgall->SaveAs(Form("%s/Trg_Pt_bin%d_All.pdf", plotfolder.Data(), i));

    //Get efficiency
    EffTrig->SetBinContent(i+1, efficiencytrg.getVal());
    EffTrig->SetBinError(i+1, 0.00001);
    EffTrk->SetBinContent(i+1, efficiencytrk.getVal());
    EffTrk->SetBinError(i+1, 0.00001);
    EffMuId->SetBinContent(i+1, efficiencyMuId.getVal());
    EffMuId->SetBinError(i+1, 0.00001);
  }
  //Save efficiency plots
  TString outputfile;
  TFile*foutput;
  if(isDataInput)foutput=new TFile(Form("%s/foutputData.root",plotfolder.Data()),"recreate");
  else foutput=new TFile(Form("%s/foutputMC.root",plotfolder.Data()),"recreate");
  TCanvas *cforSave= new TCanvas("cforSave","",600,600);
  cforSave->cd();
  EffTrig->Draw("pe");
  EffTrk->Draw("pe same");
  EffMuId->Draw("pe same");
  TLegend *leg = myLegend(0.6778523,0.399651,0.8775168,0.6073298);
  leg->AddEntry(EffTrig, "Trig", "p");
  leg->AddEntry(EffTrk, "Trk", "p");
  leg->AddEntry(EffMuId, "ID", "p");
  leg->Draw("same");
  EffTrig->Write();
  EffTrk->Write();
  EffMuId->Write();
  cforSave->SaveAs(Form("%s/eff.pdf",plotfolder.Data()));
  cforSave->Write();
  EffTrig->Draw("pe");
  cforSave->SaveAs(Form("%s/EffTrig.pdf",plotfolder.Data()));
  EffTrk->Draw("pe");
  cforSave->SaveAs(Form("%s/EffTrk.pdf",plotfolder.Data()));
  EffMuId->Draw("pe");
  cforSave->SaveAs(Form("%s/EffMuId.pdf",plotfolder.Data()));

/*
  for(int i = 0; i < nMuEtaBin; i++)
  {
    RooDataHist dataTrkPass("dataTrkPass", "dataTrkPass", mass, hTrkEtaPass[i]);
    RooDataHist dataTrkFail("dataTrkFail", "dataTrkFail", mass, hTrkEtaFail[i]);
    RooDataHist dataTrkAll("dataTrkAll", "dataTrkAll", mass, hTrkEtaAll[i]);
  }
*/
}


