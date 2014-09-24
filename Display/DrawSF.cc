#include <TH1.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include "utilities.h"
using namespace std;
string effname = "";
struct name_struct{string _name;};
struct name_struct _type[3] = {{"Trig"}, {"Trk"}, {"MuId"}};
TCanvas *c[3];
void DrawSF(){
  gStyle->SetOptStat("0");
  TFile *infMC = new TFile("0913_DefaultMuId_PlotsBinFitMC/foutputMC.root");
  TFile *inf = new TFile("0913_DefaultMuId_PlotsBinFitData/foutputData.root");
//  for(int i = 0; i < 1; i++){
  for(int i = 0; i < 3; i++){
    c[i] = new TCanvas("c1", "c1",0,44,600,600);
    c[i]->cd();
    effname = _type[i]._name;
    TH1D* effMC = (TH1D*)infMC->Get(Form("Eff%s",effname.c_str()));
    TH1D* eff = (TH1D*)inf->Get(Form("Eff%s",effname.c_str()));
	removeError(effMC);
	removeError(eff);
    TH1D* sf = (TH1D*)eff->Clone("sf");
    sf->Divide(effMC);
    c[i]->SetLeftMargin(0.15);
//    sf->SetMaximum(1.5);
    sf->SetMinimum(0);
    sf->SetYTitle(Form("%s scaling factor", effname.c_str()));
    sf->GetYaxis()->SetTitleOffset(1.2);
    sf->Draw("pe");
    c[i]->SaveAs(Form("0913ScalingFactor/SF%s.pdf", effname.c_str()));
  }
}
