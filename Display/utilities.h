#include <iostream>
#include <utility>
#include <math.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCut.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TChain.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TString.h>
#include "format.h"

#define MAX_MUON 64
#define MUON_MASS   0.10565837
#define JPSI_MASS   3.096916

TLatex myLatex(){
  TLatex t1;
  t1.SetNDC();
  t1.SetTextFont(22);
  t1.SetTextColor(2);
  t1.SetTextSize(0.05);
  return t1;
}


TH1D* myTH1D(std::string th1name, std::string xtitle, std::string ytitle, int _color, int nMuPtBin, const double MuPtBin[]){
  TH1D* eff = new TH1D(th1name.c_str(), "", nMuPtBin, MuPtBin);
  eff->SetMinimum(0);
  eff->SetMaximum(1.2);
  eff->SetLineWidth(2);
  eff->SetMarkerSize(1.2);
  eff->SetMarkerColor(_color);
  eff->SetXTitle(xtitle.c_str());
  eff->SetYTitle(ytitle.c_str());
  return eff;
};

TF1 *myTF1(std::string tf1name, double blow, double bhigh){
  TF1 *mytf1 = new TF1(tf1name.c_str(),"[0]*([4]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[4])*Gaus(x,[1],[3])/(sqrt(2*3.14159)*[3]))+[5]+[6]*x");
//  TF1 *mytf1 = new TF1(tf1name.c_str(),"[0]*([4]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[4])*Gaus(x,[1],[3])/(sqrt(2*3.14159)*[3]))");
  mytf1->SetParameter(1,3.097);
  mytf1->SetParameter(2,0.01);
  mytf1->SetParameter(3,0.03);
  mytf1->SetParameter(4,0.5);
  mytf1->SetParameter(5,0.);
  mytf1->SetParameter(6,0.);
  mytf1->SetParLimits(2,0,0.15);
  mytf1->SetParLimits(3,0,0.15);
  mytf1->SetParLimits(4,0.,1.);
//  mytf1->SetParLimits(5,0,9999);
//  mytf1->SetParLimits(6,-9999,0);
  mytf1->SetRange(blow, bhigh);
  return mytf1;
};

TLegend *myLegend(double x1,double y1,double x2, double y2){
  TLegend *leg = new TLegend(x1,y1,x2,y2);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  return leg;
};

bool SoftMuSel(MuonInfoBranches MuonInfo, int mu1){
  if(
     (MuonInfo.isTrackerMuon[mu1] || MuonInfo.isGlobalMuon[mu1]) &&
     (MuonInfo.muqual[mu1] & 1<<12) &&
     (MuonInfo.dzPV[mu1] < 30.) &&
     (MuonInfo.dxyPV[mu1] < 3.) &&
     (MuonInfo.normchi2[mu1] < 1.8) &&
     (MuonInfo.i_nPixelLayer[mu1] >= 1) &&
     ((MuonInfo.i_nPixelLayer[mu1]+MuonInfo.i_nStripLayer[mu1]) >= 6)
  ) return true;
  else return false;
};
bool KisooTrackSel(MuonInfoBranches MuonInfo, int mu1){
  if(
     (MuonInfo.dzPV[mu1] < 30.) &&
     (MuonInfo.dxyPV[mu1] < 3.) &&
     (MuonInfo.normchi2[mu1] < 1.8) &&
     (MuonInfo.i_nPixelLayer[mu1] > 1) &&//kisoo's sel
     ((MuonInfo.i_nPixelLayer[mu1]+MuonInfo.i_nStripLayer[mu1]) >= 6)
  ) return true;
  else return false;
};
bool KisooGlobalSel(MuonInfoBranches MuonInfo, int mu1){
  if(
    (MuonInfo.muqual[mu1] & 1<<12) &&
    (MuonInfo.muqual[mu1] & 1<<4)
  ) return true;
  else return false;
};
bool Acc(MuonInfoBranches MuonInfo, int mu1){
  if(fabs(MuonInfo.eta[mu1]) <= 1.3 && MuonInfo.pt[mu1] > 3.3
  ) return true;
  else if(fabs(MuonInfo.eta[mu1]) > 1.3 && fabs(MuonInfo.eta[mu1]) <= 2.2 && MuonInfo.pt[mu1] > 2.9
  ) return true;
  if(fabs(MuonInfo.eta[mu1]) > 2.2 && fabs(MuonInfo.eta[mu1]) <= 2.4 && MuonInfo.pt[mu1] > 0.8
  ) return true;
  else return false;
};
