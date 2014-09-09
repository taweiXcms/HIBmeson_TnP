#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <TCanvas.h>
#include "TH1D.h"

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

#define NUM_BX 9000


void myfittino(){

  using namespace std;
  using namespace RooFit;

  int numCPU = 6;
  bool quiet = false;

  double mumulow=2.6;
  double mumuhigh=3.5;
  double mumutrklow=1.;
  double mumutrkhigh=5.;
  
  int rebin=1;
  
  RooRealVar mass("mass", "mass", mumutrklow, mumutrkhigh);
  mass.setBins(100);

  TH1D* hTrkPtPass;
  TH1D* hTrkPtFail;
  TH1D* hTrkPtAll;

  TFile*finput=new TFile("/afs/cern.ch/user/g/ginnocen/public/TnPResults/foutputData.root","read");
  finput->cd();

  RooCategory sampleTrk("sampleTrk", "sampleTrk");
  sampleTrk.defineType("Pass");
  sampleTrk.defineType("Fail");

  TCanvas* cpass = new TCanvas("cpass", "", 0, 20, 600, 600);
  TCanvas* cfail = new TCanvas("cfail", "", 0, 40, 600, 600);
  TCanvas* call = new TCanvas("call", "", 0, 60, 600, 600);

  hTrkPtPass = (TH1D*)finput->Get("hTrkPtPass0");
  hTrkPtFail = (TH1D*)finput->Get("hTrkPtFail0");
  hTrkPtAll = (TH1D*)finput->Get("hTrkPtAll0");
  
  hTrkPtPass->Rebin(rebin);
  hTrkPtFail->Rebin(rebin);
  hTrkPtAll->Rebin(rebin);
  

  //make CB function
  RooRealVar mean_CB("mean_CB", "mean_CB", 3.1, 3.0, 3.2);
  RooRealVar sigma_CB("sigma_CB", "sigma_CB", 0.05);
  sigma_CB.setConstant(kFALSE);
  RooRealVar alpha_CB("alpha_CB", "alpha_CB", 2.0, 1.0, 5.0);
  RooRealVar n_CB("n_CB", "n_CB", 1, 0.5, 100.0);

  RooCBShape signal_CB_Pass("signal_CB_Pass", "signal_CB_Pass", mass, mean_CB, sigma_CB, alpha_CB, n_CB);
  RooCBShape signal_CB_Fail("signal_CB_Fail", "signal_CB_Fail", mass, mean_CB, sigma_CB, alpha_CB, n_CB);

  //make 2gauss function
  RooRealVar mean_gaus("mean_gaus", "mean_gaus", 3.1, 3.0, 3.2);
  RooRealVar sigma1_gaus("sigma1_gaus", "sigma1_gaus", 0.15, 0.05, 0.25);
  RooRealVar sigma2_gaus("sigma2_gaus", "sigma2_gaus", 0.02, 0.01, 0.1);
  RooGaussian gauss1("gauss1", "gauss1", mass, mean_gaus, sigma1_gaus);
  RooGaussian gauss2("gauss2", "gauss2", mass, mean_gaus, sigma2_gaus);
  RooRealVar mfrac_gaus("mfrac_gaus", "mfrac_gaus", 0.2, 0.0, 1.0);

  RooAddPdf signal_gauss_Pass("signal_gauss_Pass", "signal_gauss_Pass", RooArgList(gauss1, gauss2), RooArgList(mfrac_gaus));
  RooAddPdf signal_gauss_Fail("signal_gauss_Fail", "signal_gauss_Fail", RooArgList(gauss1, gauss2), RooArgList(mfrac_gaus));

  //make Chebyschev function
  RooRealVar cheb_p1Pass("cheb_paPass", "cheb_p1Pass", 0., -1., +1.);
  RooRealVar cheb_p2Pass("cheb_p2Pass", "cheb_p2Pass", 0., -1., +1.);
  RooChebychev background_cheb_Pass("background_cheb_Pass", "background_cheb_Pass", mass, RooArgList(cheb_p1Pass, cheb_p2Pass));
  RooRealVar cheb_p1Fail("cheb_p1Fail", "cheb_p1Fail", 0., -1., +1.);
  RooRealVar cheb_p2Fail("cheb_p2Fail", "cheb_p2Fail", 0., -1., +1.);
  RooChebychev background_cheb_Fail("background_cheb_Fail", "background_cheb_Fail", mass, RooArgList(cheb_p1Fail, cheb_p2Fail));

  //make exponential function
  RooRealVar lp("lp", "lp", 0.0, -5., 5.);
  RooRealVar lf("lf", "lf", 0.0, -5., 5.);
  RooExponential background_exp_Pass("background_exp_Pass", "background_exp_Pass", mass, lp);
  RooExponential background_exp_Fail("background_exp_Fail", "background_exp_Fail", mass, lf);

  //set axis variable
  RooRealVar weight("weight", "weight", 1., 0, 1E9);

  RooDataHist dataTrkPass("dataTrkPass", "dataTrkPass", mass, hTrkPtPass);
  RooDataHist dataTrkFail("dataTrkFail", "dataTrkFail", mass, hTrkPtFail);
  RooDataHist dataTrkAll("dataTrkAll", "dataTrkFAll", mass, hTrkPtAll);

  RooDataSet dataTrkPassset("dataTrkPassset", "mass", RooArgSet(mass, weight), WeightVar("weight"));
  RooDataSet dataTrkFailset("dataTrkFailset", "mass", RooArgSet(mass, weight), WeightVar("weight"));

  for(int i = 0; i < dataTrkPass.numEntries(); i++)
  {
    dataTrkPassset.add(*dataTrkPass.get(i), dataTrkPass.weight());
    dataTrkFailset.add(*dataTrkFail.get(i), dataTrkFail.weight());
  }
  RooDataSet dataTrkAllset("dataTrkAllset", "mass", RooArgSet(mass, weight), Index(sampleTrk), Import("Pass", dataTrkPassset), Import("Fail", dataTrkFailset), WeightVar(weight));

  RooRealVar nsig("nsig", "nsig", 1, 0, 1E9);
  RooRealVar efficiency("efficiency", "efficiency", 0.9, 0., 1.);
  RooFormulaVar nsigPass("nsigPass", "nsigPass", "nsig*efficiency", RooArgList(nsig, efficiency));
  RooFormulaVar nsigFail("nsigFail", "nsigFail", "nsig*(1-efficiency)", RooArgList(nsig, efficiency));

  RooRealVar nbkgPass("nbkgPass", "nbkgPass", 1, 0, 1E9);
  RooRealVar nbkgFail("nbkgFail", "nbkgFail", 1, 0, 1E9);

  RooAddPdf model_gaus_cheb_Pass("model_gaus_cheb_Pass", "model_gaus_cheb_Pass", RooArgList(signal_gauss_Pass, background_cheb_Pass), RooArgList(nsigPass, nbkgPass));
  RooAddPdf model_gaus_cheb_Fail("model_gaus_cheb_Fail", "model_gaus_cheb_Fail", RooArgList(signal_gauss_Fail, background_cheb_Fail), RooArgList(nsigFail, nbkgFail));

  //fit
  RooSimultaneous simPdf("simPdf", "simPdf", sampleTrk);
  simPdf.addPdf(model_gaus_cheb_Pass, "Pass");
  simPdf.addPdf(model_gaus_cheb_Fail, "Fail");

  simPdf.fitTo(dataTrkAllset, Extended(true), NumCPU(numCPU), PrintLevel(quiet?-1:1), Warnings(!quiet));

  RooPlot *frame1 = mass.frame();
  dataTrkAllset.plotOn(frame1, Cut("sampleTrk==sampleTrk::Pass"));

  simPdf.plotOn(frame1, Slice(sampleTrk, "Pass"), ProjWData(sampleTrk,dataTrkAllset),LineStyle(kDashed), LineColor(kBlue));
  simPdf.plotOn(frame1, Slice(sampleTrk, "Pass"), Components("signal_gauss_Pass"),ProjWData(sampleTrk,dataTrkAllset),LineStyle(kDashed), LineColor(kGreen));
  simPdf.plotOn(frame1, Slice(sampleTrk, "Pass"), Components("background_cheb_Pass"),ProjWData(sampleTrk,dataTrkAllset),LineStyle(kDashed), LineColor(kRed));

  cpass->cd();
  frame1->Draw();

  RooPlot *frame2 = mass.frame();
  dataTrkAllset.plotOn(frame2, Cut("sampleTrk==sampleTrk::Fail"));

  simPdf.plotOn(frame2, Slice(sampleTrk, "Fail"), ProjWData(sampleTrk,dataTrkAllset),LineStyle(kDashed), LineColor(kBlue));
  simPdf.plotOn(frame2, Slice(sampleTrk, "Fail"), Components("signal_gauss_Fail"),ProjWData(sampleTrk,dataTrkAllset),LineStyle(kDashed), LineColor(kGreen));
  simPdf.plotOn(frame2, Slice(sampleTrk, "Fail"), Components("background_cheb_Fail"),ProjWData(sampleTrk,dataTrkAllset),LineStyle(kDashed), LineColor(kRed));

  cfail->cd();
  frame2->Draw();

  RooPlot *frame3 = mass.frame();
  dataTrkAllset.plotOn(frame3);

  simPdf.plotOn(frame3, ProjWData(sampleTrk, dataTrkAllset), LineColor(kGreen));
  simPdf.plotOn(frame3, Components("signal_gauss_Pass, signal_gauss_Fail"), ProjWData(sampleTrk,dataTrkAllset),LineStyle(kDashed), LineColor(kGreen));
//  simPdf.plotOn(frame3, ProjWData(sampleTrk,dataTrkAllset),LineStyle(kDashed), LineColor(kRed));

  call->cd();
  frame3->Draw();

}


