#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
using namespace std;
string effname = "trg";

void DrawSF(){
  gStyle->SetOptStat("0");
//  TFile *infMC = new TFile("../TnPOutMC.root");
//  TFile *inf = new TFile("../TnPOut.root");
  TFile *infMC = new TFile("TnPRooFit.root");
  TFile *inf = new TFile("TnPRooFitMC.root");
  TH1D* effMC = (TH1D*)infMC->Get(Form("eff_%s",effname.c_str()));
  TH1D* eff = (TH1D*)inf->Get(Form("eff_%s",effname.c_str()));
  TH1D* sf = (TH1D*)eff->Clone("sf");
  sf->Divide(effMC);
  TCanvas *c1 = new TCanvas("c1", "c1",0,44,600,600);
  c1->SetLeftMargin(0.15);
//  sf->SetMaximum(1.5);
  sf->SetMinimum(0);
  sf->SetYTitle(Form("%s scaling factor", effname.c_str()));
  sf->GetYaxis()->SetTitleOffset(1.2);
  sf->Draw("p");
  c1->SaveAs(Form("SF_%s.pdf", effname.c_str()));
}
