#include "format.h"
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
  
using namespace std;
using namespace RooFit;

//const int nMuPtBin = 1;
const int nMuPtBin = 7;
//const double MuPtBin[nMuPtBin+1] = {0,1.5};
const double MuPtBin[nMuPtBin+1] = {0.0,1.5,3.0,4.5,6.0,9.0,20.0,30.0};
int numCPU = 6;
bool quiet = false;

TString ninf = "/afs/cern.ch/user/g/ginnocen/public/foutputData.root";
TString noutf = "ResultsData/TnPRooFit.root";
TString plotfolder = "PlotsRooFitData";

//TString ninf = "/afs/cern.ch/user/g/ginnocen/public/foutputMC.root";
//TString noutf = "ResultsMC/TnPRooFit.root";
//TString plotfolder = "PlotsRooFitMC";

void doRooFit(double &npass, double &nfail, double &npass_err, double &nfail_err,
TCanvas* cpass, TCanvas* cfail, TCanvas* call, double ptlow, double pthigh, TString tree_type,
int pdfset
){
  TFile *f = new TFile(ninf); 
  TTree* TnPtree;
  TnPtree = (TTree*)f->Get("TnPtree"+tree_type);

  //RooRealVar mass("mass","mass [GeV]", 2.6, 3.5);
  RooRealVar mass("mass","mass [GeV]", 2., 5.);
  RooRealVar pass("pass","pass", 0, 1);
  RooRealVar pt("pt","pt",0,100);	
  RooRealVar eta("eta","eta",-5,5);	
  RooDataSet* data_set;
  RooDataSet* data_setPass;
  RooDataSet* data_setFail;
  data_set = new RooDataSet("data_set","data_set",TnPtree,RooArgSet(mass, pt, eta, pass));
  TString CutPass, CutFail;
  CutPass = TString::Format("pt > %f && pt < %f && pass == 1", ptlow, pthigh);
  data_setPass = (RooDataSet*)data_set->reduce(CutPass);
  CutFail = TString::Format("pt > %f && pt < %f && pass == 0", ptlow, pthigh);
  data_setFail = (RooDataSet*)data_set->reduce(CutFail);
  RooCategory sample("sample","sample");
  sample.defineType("Pass");
  sample.defineType("Fail");
  RooDataSet combData("combData","combData",mass,Index(sample),Import("Pass",*data_setPass),Import("Fail",*data_setFail));
cout<<TnPtree->GetEntries()<<endl;
cout<<data_set->sumEntries()<<endl;
cout<<data_setPass->sumEntries()<<endl;
cout<<data_setFail->sumEntries()<<endl;

  RooRealVar nsigPass("nsigPass","nsigPass",data_setPass->sumEntries()*0.9,0.,1E10);
  RooRealVar nsigFail("nsigFail","nsigFail",data_setFail->sumEntries()*0.9,0.,1E10);
  RooRealVar nbkgPass("nbkgPass","nbkgPass",data_setPass->sumEntries()*0.1,0.,1E10);
  RooRealVar nbkgFail("nbkgFail","nbkgFail",data_setFail->sumEntries()*0.1,0.,1E10);

  RooRealVar mean_CB("mean_CB","mean_CB",3.1, 3.0, 3.2);
  RooRealVar sigma_CB("sigma_CB","sigma_CB", 0.05);
  sigma_CB.setConstant(kFALSE);
  RooRealVar alpha_CB("alpha_CB","alpha_CB", 2.0, 1.0, 5.0);
  RooRealVar n_CB("n_CB","n_CB", 1, 0.5, 100.0);
  RooCBShape signal_CB_Pass("signal_CB_Pass", "signal_CB_Pass", mass, mean_CB, sigma_CB, alpha_CB, n_CB);
  RooCBShape signal_CB_Fail("signal_CB_Fail", "signal_CB_Fail", mass, mean_CB, sigma_CB, alpha_CB, n_CB);

  RooRealVar cheb_p1Pass("cheb_p1Pass","cheb_p1Pass",0.,-1.,+1.);
  RooRealVar cheb_p2Pass("cheb_p2Pass","cheb_p2Pass",0.,-1.,+1.);
  RooChebychev background_cheb_Pass("background_cheb_Pass","background_cheb_Pass",mass,RooArgList(cheb_p1Pass,cheb_p2Pass));
  RooRealVar cheb_p1Fail("cheb_p1Fail","cheb_p1Fail",0.,-1.,+1.);
  RooRealVar cheb_p2Fail("cheb_p2Fail","cheb_p2Fail",0.,-1.,+1.);
  RooChebychev background_cheb_Fail("background_cheb_Fail","background_cheb_Fail",mass,RooArgList(cheb_p1Fail,cheb_p2Fail));

  RooRealVar lp("lp","lp",0.0, -5., 5.);
  RooRealVar lf("lt","lt",0.0, -5., 5.);
  RooExponential background_exp_Pass("background_exp_Pass","background_exp_Pass",mass, lp);
  RooExponential background_exp_Fail("background_exp_Fail","background_exp_Fail",mass, lf);
    
  RooRealVar mean_gaus("mean_gaus","mean_gaus",3.1, 3.0, 3.2);
  RooRealVar sigma1_gaus("sigma1_gaus","sigma1_gaus",0.15,0.05,0.25);
  RooRealVar sigma2_gaus("sigma2_gaus","sigma2_gaus",0.02 ,0.01,0.1);
  RooGaussian gauss1("gauss1","gauss1",mass,mean_gaus,sigma1_gaus);
  RooGaussian gauss2("gauss2","gauss2",mass,mean_gaus,sigma2_gaus);
  RooRealVar mfrac_gaus("mfrac_gaus","mfrac_gaus",0.2,0.0,1.0);
  RooAddPdf signal_gaus_Pass("signal_gaus_Pass","signal_gaus_Pass",RooArgList(gauss1,gauss2),RooArgList(mfrac_gaus));
  RooAddPdf signal_gaus_Fail("signal_gaus_Fail","signal_gaus_Fail",RooArgList(gauss1,gauss2),RooArgList(mfrac_gaus));

  //ID CB + cheb
  //Trg CB + exp
  //Trk 2-Gaus + cheb
  RooAddPdf model_CB_cheb_Pass("model_CB_cheb_Pass","model_CB_cheb_Pass",RooArgList(signal_CB_Pass, background_cheb_Pass),RooArgList(nsigPass, nbkgPass));
  RooAddPdf model_CB_cheb_Fail("model_CB_cheb_Fail","model_CB_cheb_Fail",RooArgList(signal_CB_Fail, background_cheb_Fail),RooArgList(nsigFail, nbkgFail));
  RooAddPdf model_CB_exp_Pass("model_CB_exp_Pass","model_CB_exp_Pass",RooArgList(signal_CB_Pass, background_exp_Pass),RooArgList(nsigPass, nbkgPass));
  RooAddPdf model_CB_exp_Fail("model_CB_exp_Fail","model_CB_exp_Fail",RooArgList(signal_CB_Fail, background_exp_Fail),RooArgList(nsigFail, nbkgFail));
  RooAddPdf model_gaus_cheb_Pass("model_gaus_cheb_Pass","model_gaus_cheb_Pass",RooArgList(signal_gaus_Pass, background_cheb_Pass),RooArgList(nsigPass, nbkgPass));
  RooAddPdf model_gaus_cheb_Fail("model_gaus_cheb_Fail","model_gaus_cheb_Fail",RooArgList(signal_gaus_Fail, background_cheb_Fail),RooArgList(nsigFail, nbkgFail));

  RooSimultaneous simPdf("simPdf","simpdf",sample);
  //pdf set 1, ID
  if(pdfset == 1){
    simPdf.addPdf(model_CB_cheb_Pass,"Pass");
    simPdf.addPdf(model_CB_cheb_Fail,"Fail");
  }
  if(pdfset == 2){
    simPdf.addPdf(model_CB_exp_Pass,"Pass");
    simPdf.addPdf(model_CB_exp_Fail,"Fail");
  }
  if(pdfset == 3){
    simPdf.addPdf(model_gaus_cheb_Pass,"Pass");
    simPdf.addPdf(model_gaus_cheb_Fail,"Fail");
  }
  simPdf.fitTo(combData, Extended(true), NumCPU(numCPU), PrintLevel(quiet?-1:1), PrintEvalErrors(quiet?-1:1), Warnings(!quiet));

  RooPlot *frame1 = mass.frame(); 
  combData.plotOn(frame1, Binning(100), Cut("sample==sample::Pass"));
  simPdf.plotOn(frame1,Slice(sample,"Pass"),ProjWData(sample,combData));
  if(pdfset == 1){
    simPdf.plotOn(frame1,Slice(sample,"Pass"),Components("signal_CB_Pass"),ProjWData(sample,combData),LineStyle(kDashed), LineColor(kGreen));
    simPdf.plotOn(frame1,Slice(sample,"Pass"),Components("background_cheb_Pass"),ProjWData(sample,combData),LineStyle(kDashed), LineColor(kRed));
  }
  if(pdfset == 2){
    simPdf.plotOn(frame1,Slice(sample,"Pass"),Components("signal_CB_Pass"),ProjWData(sample,combData),LineStyle(kDashed), LineColor(kGreen));
    simPdf.plotOn(frame1,Slice(sample,"Pass"),Components("background_exp_Pass"),ProjWData(sample,combData),LineStyle(kDashed), LineColor(kRed));
  }
  if(pdfset == 3){
   simPdf.plotOn(frame1,Slice(sample,"Pass"),Components("signal_gaus_Pass"),ProjWData(sample,combData),LineStyle(kDashed), LineColor(kGreen));
    simPdf.plotOn(frame1,Slice(sample,"Pass"),Components("background_cheb_Pass"),ProjWData(sample,combData),LineStyle(kDashed), LineColor(kRed));
  }
  cpass->cd();
  frame1->SetYTitle("");
  frame1->Draw();
  TLatex t1 = myLatex();
  t1.DrawLatex(0.17, 0.85, "Passing probes");
  t1.DrawLatex(0.17, 0.80, "efficiency = "+tree_type);
  t1.DrawLatex(0.17, 0.75, Form("p_{T} = %.1f ~ %.1f GeV", ptlow, pthigh));
  npass = nsigPass.getVal();
  npass_err = nsigPass.getError();

  RooPlot *frame2 = mass.frame();
  combData.plotOn(frame2, Binning(100), Cut("sample==sample::Fail"));
  simPdf.plotOn(frame2,Slice(sample,"Fail"),ProjWData(sample,combData));
  if(pdfset == 1){
    simPdf.plotOn(frame2,Slice(sample,"Fail"),Components("signal_CB_Fail"),ProjWData(sample,combData),LineStyle(kDashed), LineColor(kGreen));
    simPdf.plotOn(frame2,Slice(sample,"Fail"),Components("background_cheb_Fail"),ProjWData(sample,combData),LineStyle(kDashed), LineColor(kRed));
  }
  if(pdfset == 2){
    simPdf.plotOn(frame2,Slice(sample,"Fail"),Components("signal_CB_Fail"),ProjWData(sample,combData),LineStyle(kDashed), LineColor(kGreen));
    simPdf.plotOn(frame2,Slice(sample,"Fail"),Components("background_exp_Fail"),ProjWData(sample,combData),LineStyle(kDashed), LineColor(kRed));
  }
  if(pdfset == 3){
    simPdf.plotOn(frame2,Slice(sample,"Fail"),Components("signal_gaus_Fail"),ProjWData(sample,combData),LineStyle(kDashed), LineColor(kGreen));
    simPdf.plotOn(frame2,Slice(sample,"Fail"),Components("background_cheb_Fail"),ProjWData(sample,combData),LineStyle(kDashed), LineColor(kRed));
  }
  cfail->cd();
  frame2->Draw();
  frame2->SetYTitle("");
  t1.DrawLatex(0.17, 0.85, "Failing probes");
  t1.DrawLatex(0.17, 0.80, "efficiency = "+tree_type);
  t1.DrawLatex(0.17, 0.75, Form("p_{T} = %.1f ~ %.1f GeV", ptlow, pthigh));
  nfail = nsigFail.getVal();
  nfail_err = nsigFail.getError();

  RooPlot *frame3 = mass.frame();
  combData.plotOn(frame3, Binning(100));
  simPdf.plotOn(frame3,ProjWData(sample,combData));
  if(pdfset == 1){
    simPdf.plotOn(frame3,Components("signal_CB_Pass, signal_CB_Fail"),ProjWData(sample,combData),LineStyle(kDashed), LineColor(kGreen));
    simPdf.plotOn(frame3,Components("background_cheb_Pass, background_cheb_Fail"),ProjWData(sample,combData),LineStyle(kDashed), LineColor(kRed));
  }
  if(pdfset == 2){
    simPdf.plotOn(frame3,Components("signal_CB_Pass, signal_CB_Fail"),ProjWData(sample,combData),LineStyle(kDashed), LineColor(kGreen));
    simPdf.plotOn(frame3,Components("background_exp_Fail, background_exp_Fail"),ProjWData(sample,combData),LineStyle(kDashed), LineColor(kRed));
  }
  if(pdfset == 3){
    simPdf.plotOn(frame3,Components("signal_gaus_Pass, signal_gaus_Fail"),ProjWData(sample,combData),LineStyle(kDashed), LineColor(kGreen));
    simPdf.plotOn(frame3,Components("background_cheb_Pass, background_cheb_Fail"),ProjWData(sample,combData),LineStyle(kDashed), LineColor(kRed));
  }
  call->cd();
  frame3->Draw();
  frame3->SetYTitle("");
  t1.DrawLatex(0.17, 0.85, "All probes");
  t1.DrawLatex(0.17, 0.80, "efficiency = "+tree_type);
  t1.DrawLatex(0.17, 0.75, Form("p_{T} = %.1f ~ %.1f GeV", ptlow, pthigh));
}

void doTnPRooFit(){
  TCanvas *croofit_trg_pass[nMuPtBin];
  TCanvas *croofit_trk_pass[nMuPtBin];
  TCanvas *croofit_id_pass[nMuPtBin];
  TCanvas *croofit_trg_fail[nMuPtBin];
  TCanvas *croofit_trk_fail[nMuPtBin];
  TCanvas *croofit_id_fail[nMuPtBin];
  TCanvas *croofit_trg_all[nMuPtBin];
  TCanvas *croofit_trk_all[nMuPtBin];
  TCanvas *croofit_id_all[nMuPtBin];
  TH1D* eff_trg = myTH1D("eff_trg", "Mu Pt (GeV)", "efficiency", 2, nMuPtBin, MuPtBin);
  TH1D* eff_trk = myTH1D("eff_trk", "Mu Pt (GeV)", "efficiency", 3, nMuPtBin, MuPtBin);
  TH1D* eff_id = myTH1D("eff_id", "Mu Pt (GeV)", "efficiency", 4, nMuPtBin, MuPtBin);
  double npass = 0; double nfail = 0;
  double npass_err = 0; double nfail_err = 0;
  for(int i = 0; i < nMuPtBin; i++){
    croofit_trg_pass[i] = new TCanvas(Form("croofit_trg_pass%d", i), "", 600, 600);
    croofit_trg_fail[i] = new TCanvas(Form("croofit_trg_fail%d", i), "", 600, 600);
    croofit_trg_all[i] = new TCanvas(Form("croofit_trg_all%d", i), "", 600, 600);
    croofit_trk_pass[i] = new TCanvas(Form("croofit_trk_pass%d", i), "", 600, 600);
    croofit_trk_fail[i] = new TCanvas(Form("croofit_trk_fail%d", i), "", 600, 600);
    croofit_trk_all[i] = new TCanvas(Form("croofit_trk_all%d", i), "", 600, 600);
    croofit_id_pass[i] = new TCanvas(Form("croofit_id_pass%d", i), "", 600, 600);
    croofit_id_fail[i] = new TCanvas(Form("croofit_id_fail%d", i), "", 600, 600);
    croofit_id_all[i] = new TCanvas(Form("croofit_id_all%d", i), "", 600, 600);
    doRooFit(npass, nfail, npass_err, nfail_err, croofit_trg_pass[i], croofit_trg_fail[i], croofit_trg_all[i], MuPtBin[i], MuPtBin[i+1], "Trg", 2); 
    eff_trg->SetBinContent(i+1, npass/(npass+nfail));
    eff_trg->SetBinError(i+1, 0.00001);
    doRooFit(npass, nfail, npass_err, nfail_err, croofit_trk_pass[i], croofit_trk_fail[i], croofit_trk_all[i], MuPtBin[i], MuPtBin[i+1], "Trk", 3); 
    eff_trk->SetBinContent(i+1, npass/(npass+nfail));
    eff_trk->SetBinError(i+1, 0.00001);
    doRooFit(npass, nfail, npass_err, nfail_err, croofit_id_pass[i], croofit_id_fail[i], croofit_id_all[i], MuPtBin[i], MuPtBin[i+1], "ID", 1); 
//    doRooFit(npass, nfail, npass_err, nfail_err, croofit_id_pass[i], croofit_id_fail[i], MuPtBin[i], MuPtBin[i+1], "ID", 2); 
    eff_id->SetBinContent(i+1, npass/(npass+nfail));
    eff_id->SetBinError(i+1, 0.00001);
  }
  TCanvas *c= new TCanvas("c","",600,600);
  c->SetTopMargin(0.07504363);
  c->cd();
  eff_trg->Draw("pe");
  eff_trk->Draw("pe same");
  eff_id->Draw("pe same");
  TLegend *leg = myLegend(0.6778523,0.399651,0.8775168,0.6073298);
  leg->AddEntry(eff_trg, "Trg", "p");
  leg->AddEntry(eff_trk, "Trk", "p");
  leg->AddEntry(eff_id, "ID", "p");
  leg->Draw("same");

  TFile *outf = new TFile(noutf,"recreate");
  outf->cd();
  eff_trg->Write();
  eff_trk->Write();
  eff_id->Write();
  c->SaveAs(Form("%s/eff.pdf",plotfolder.Data()));
  c->Write();

  eff_trg->Draw("pe");
  c->SaveAs(Form("%s/eff_trg.pdf",plotfolder.Data()));
  eff_trk->Draw("pe");
  c->SaveAs(Form("%s/eff_trk.pdf",plotfolder.Data()));
  eff_id->Draw("pe");
  c->SaveAs(Form("%s/eff_id.pdf",plotfolder.Data()));
  for(int i = 0; i < nMuPtBin; i++){
    croofit_trg_pass[i]->Write();
    croofit_trk_pass[i]->Write();
    croofit_id_pass[i]->Write();
    croofit_trg_fail[i]->Write();
    croofit_trk_fail[i]->Write();
    croofit_id_fail[i]->Write();
    croofit_trg_all[i]->Write();
    croofit_trk_all[i]->Write();
    croofit_id_all[i]->Write();
    croofit_trg_pass[i]->SaveAs(Form("%s/croofit_trg_pass_%d.pdf",plotfolder.Data(), i));
    croofit_trk_pass[i]->SaveAs(Form("%s/croofit_trk_pass_%d.pdf",plotfolder.Data(), i));
    croofit_id_pass[i]->SaveAs(Form("%s/croofit_id_pass_%d.pdf",plotfolder.Data(), i));
    croofit_trg_fail[i]->SaveAs(Form("%s/croofit_trg_fail_%d.pdf",plotfolder.Data(), i));
    croofit_trk_fail[i]->SaveAs(Form("%s/croofit_trk_fail_%d.pdf",plotfolder.Data(), i));
    croofit_id_fail[i]->SaveAs(Form("%s/croofit_id_fail_%d.pdf",plotfolder.Data(), i));
    croofit_trg_all[i]->SaveAs(Form("%s/croofit_trg_all_%d.pdf",plotfolder.Data(), i));
    croofit_trk_all[i]->SaveAs(Form("%s/croofit_trk_all_%d.pdf",plotfolder.Data(), i));
    croofit_id_all[i]->SaveAs(Form("%s/croofit_id_all_%d.pdf",plotfolder.Data(), i));
  }
  outf->Close();
}
