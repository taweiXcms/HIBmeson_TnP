#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "../Display/utilities.h"
#include "TH1D.h"
#define NUM_BX 9000


//void fitJpsi(bool isDataInput=true){
void fitJpsi(bool isDataInput=false){

  TString plotfolder = "PlotsData";
  if(!isDataInput) plotfolder = "PlotsMC";

  int nBin = 400;
  double mumulow=2.6;
  double mumuhigh=3.5;
  double mumuTrklow=2.0;
  double mumuTrkhigh=5.0;

  TString outputfile;
  if(isDataInput) outputfile="../../TnPfiles/foutputData.root";
  else outputfile="../../TnPfiles/foutputMC.root";
  TFile*foutput=new TFile(outputfile.Data(),"recreate");
  double massTrg, massTrk, massID;
  double passTrg, passTrk, passID;
  double ptTrg, ptTrk, ptID;
  double etaTrg, etaTrk, etaID;
  double genTrg, genTrk, genID;

  TTree* TnPtreeTrg = new TTree("TnPtreeTrg","");
  TTree* TnPtreeTrk = new TTree("TnPtreeTrk","");
  TTree* TnPtreeID = new TTree("TnPtreeID","");
  
  TnPtreeTrg->Branch("mass",&massTrg);
  TnPtreeTrg->Branch("pt",&ptTrg);
  TnPtreeTrg->Branch("eta",&etaTrg);
  TnPtreeTrg->Branch("pass",&passTrg);
  TnPtreeTrg->Branch("gen",&genTrg);

  TnPtreeTrk->Branch("mass",&massTrk);
  TnPtreeTrk->Branch("pt",&ptTrk);
  TnPtreeTrk->Branch("eta",&etaTrk);
  TnPtreeTrk->Branch("pass",&passTrk);
  TnPtreeTrk->Branch("gen",&genTrk);

  TnPtreeID->Branch("mass",&massID);
  TnPtreeID->Branch("pt",&ptID);
  TnPtreeID->Branch("eta",&etaID);
  TnPtreeID->Branch("pass",&passID);
  TnPtreeID->Branch("gen",&genID);

  const int nMuPtBin = 7;
//  const int nMuPtBin = 10;
//  const int nMuPtBin = 18;
  const double MuPtBin[nMuPtBin+1] = {0.0,1.5,3.0,4.5,6.0,9.0,20.0,30.0};
//  const double MuPtBin[nMuPtBin+1] = {0.0, 1.5, 2.9, 3.3, 3.5, 4.0, 4.5, 6.0, 9.0, 20.0, 30.0};
//  const double MuPtBin[nMuPtBin+1] = {0.0, 1.5, 2.5, 3.0, 3.1, 3.125, 3.15, 3.175, 3.2, 3.225, 3.25, 3.3, 3.5, 4.0, 4.5, 6.0, 9.0, 20.0, 30.0};

  const int nMuEtaBin = 5;
  const double MuEtaBin[nMuEtaBin+1] = {-2.4,-1.5,-0.5,0.5,1.5,2.4};
    
  Bool_t IsMuonInAcceptance(Float_t,Float_t,Float_t);
  Bool_t IsTag(Bool_t, Int_t, Bool_t, Bool_t, Bool_t);
  
  TString infname;
  
  if(isDataInput) infname="/net/hisrv0001/home/tawei/tawei/TnPBntuple/TnPnt_20140802_PAMuon_HIRun2013_28Sep2013_v1.root";
  else infname="/net/hisrv0001/home/tawei/tawei/TnPBntuple/TnPnt_BoostedMC_20140806_Kp.root";
//  else infname="/net/hisrv0001/home/tawei/tawei/TnPBntuple/TnPnt_BoostedMC_20140905_Kp.root";
//  else infname="/net/hisrv0001/home/tawei/tawei/TnPBntuple/TnPnt_BoostedMC_20140906_Kp.root";
  TFile *inf = new TFile(infname.Data());
  
  TTree *ntuple = (TTree*) inf->Get("ntJpsi");
  TTree *nt_mcgen = (TTree*)inf->Get("ntGen");
  ntuple->AddFriend(nt_mcgen);
  
  //Pt probes
  //trigger efficiency
  TH1D* hTrigPtPass[nMuPtBin];
  TH1D* hTrigPtFail[nMuPtBin];
  TH1D* hTrigPtAll[nMuPtBin];
  TH1D* hEtaTrigPtPass[nMuPtBin];
  TH1D* hEtaTrigPtFail[nMuPtBin];
  TH1D* hEtaTrigPtAll[nMuPtBin];
  
  //tracking efficiency
  TH1D* hTrkPtPass[nMuPtBin];
  TH1D* hTrkPtFail[nMuPtBin];
  TH1D* hTrkPtAll[nMuPtBin];
  TH1D* hEtaTrkPtPass[nMuPtBin];
  TH1D* hEtaTrkPtFail[nMuPtBin];
  TH1D* hEtaTrkPtAll[nMuPtBin];
  
  //Muon ID efficiency
  TH1D* hMuIdPtPass[nMuPtBin];
  TH1D* hMuIdPtFail[nMuPtBin];
  TH1D* hMuIdPtAll[nMuPtBin];
  TH1D* hEtaMuIdPtPass[nMuPtBin];
  TH1D* hEtaMuIdPtFail[nMuPtBin];
  TH1D* hEtaMuIdPtAll[nMuPtBin];

  //Eta probes
  //trigger efficiency
  TH1D* hTrigEtaPass[nMuEtaBin];
  TH1D* hTrigEtaFail[nMuEtaBin];
  TH1D* hTrigEtaAll[nMuEtaBin];
  
  //tracking efficiency
  TH1D* hTrkEtaPass[nMuEtaBin];
  TH1D* hTrkEtaFail[nMuEtaBin];
  TH1D* hTrkEtaAll[nMuEtaBin];
  
  //Muon ID efficinecy
  TH1D* hMuIdEtaPass[nMuEtaBin];
  TH1D* hMuIdEtaFail[nMuEtaBin];
  TH1D* hMuIdEtaAll[nMuEtaBin];

  //Use to store quick sideband method plots
  TH1D* EffTrig = myTH1D("EffTrig", "Mu Pt (GeV)", "efficiency", 2, nMuPtBin, MuPtBin);
  TH1D* EffTrk = myTH1D("EffTrk", "Mu Pt (GeV)", "efficiency", 3, nMuPtBin, MuPtBin);
  TH1D* EffMuId = myTH1D("EffMuId", "Mu Pt (GeV)", "efficiency", 4, nMuPtBin, MuPtBin);

  for(int i = 0; i < nMuPtBin; i++){
    hTrigPtPass[i] = new TH1D(Form("hTrigPtPass%d",i),"",nBin,mumulow,mumuhigh);
    hTrigPtFail[i] = new TH1D(Form("hTrigPtFail%d",i),"",nBin,mumulow,mumuhigh);
    hTrigPtAll[i] = new TH1D(Form("hTrigPtAll%d",i),"",nBin,mumulow,mumuhigh);
    hTrkPtPass[i] = new TH1D(Form("hTrkPtPass%d",i),"",nBin,mumuTrklow,mumuTrkhigh);
    hTrkPtFail[i] = new TH1D(Form("hTrkPtFail%d",i),"",nBin,mumuTrklow,mumuTrkhigh);
    hTrkPtAll[i] = new TH1D(Form("hTrkPtAll%d",i),"",nBin,mumuTrklow,mumuTrkhigh);
    hMuIdPtPass[i] = new TH1D(Form("hMuIdPtPass%d",i),"",nBin,mumulow,mumuhigh);
    hMuIdPtFail[i] = new TH1D(Form("hMuIdPtFail%d",i),"",nBin,mumulow,mumuhigh);
    hMuIdPtAll[i] = new TH1D(Form("hMuIdPtAll%d",i),"",nBin,mumulow,mumuhigh);

    hEtaTrigPtPass[i] = new TH1D(Form("hEtaTrigPtPass%d",i),"",nBin,-3, 3);
    hEtaTrigPtFail[i] = new TH1D(Form("hEtaTrigPtFail%d",i),"",nBin,-3, 3);
    hEtaTrigPtAll[i]  = new TH1D(Form("hEtaTrigPtAll%d",i),"" ,nBin,-3, 3);
    hEtaTrkPtPass[i]  = new TH1D(Form("hEtaTrkPtPass%d",i),"" ,nBin,-3, 3);
    hEtaTrkPtFail[i]  = new TH1D(Form("hEtaTrkPtFail%d",i),"" ,nBin,-3, 3);
    hEtaTrkPtAll[i]   = new TH1D(Form("hEtaTrkPtAll%d",i),""  ,nBin,-3, 3);
    hEtaMuIdPtPass[i] = new TH1D(Form("hEtaMuIdPtPass%d",i),"",nBin,-3, 3);
    hEtaMuIdPtFail[i] = new TH1D(Form("hEtaMuIdPtFail%d",i),"",nBin,-3, 3);
    hEtaMuIdPtAll[i]  = new TH1D(Form("hEtaMuIdPtAll%d",i),"" ,nBin,-3, 3);
  }
  
  for(int i = 0; i < nMuEtaBin; i++){
    hTrigEtaPass[i] = new TH1D(Form("hTrigEtaPass%d",i),"",nBin,mumulow,mumuhigh);
    hTrigEtaFail[i] = new TH1D(Form("hTrigEtaFail%d",i),"",nBin,mumulow,mumuhigh);
    hTrigEtaAll[i] = new TH1D(Form("hTrigEtaAll%d",i),"",nBin,mumulow,mumuhigh);
    hTrkEtaPass[i] = new TH1D(Form("hTrkEtaPass%d",i),"",nBin,mumuTrklow,mumuTrkhigh);
    hTrkEtaFail[i] = new TH1D(Form("hTrkEtaFail%d",i),"",nBin,mumuTrklow,mumuTrkhigh);
    hTrkEtaAll[i] = new TH1D(Form("hTrkEtaAll%d",i),"",nBin,mumuTrklow,mumuTrkhigh);
    hMuIdEtaPass[i] = new TH1D(Form("hMuIdEtaPass%d",i),"",nBin,mumulow,mumuhigh);
    hMuIdEtaFail[i] = new TH1D(Form("hMuIdEtaFail%d",i),"",nBin,mumulow,mumuhigh);
    hMuIdEtaAll[i] = new TH1D(Form("hMuIdEtaAll%d",i),"",nBin,mumulow,mumuhigh);
  }

  int Run,size,Event;
  int isTriggered1[NUM_BX];
  int isTriggered2[NUM_BX];
  Float_t pt[NUM_BX];
  Float_t mass[NUM_BX];
  Float_t genpt[NUM_BX];
  Float_t eta[NUM_BX];
  Float_t y[NUM_BX];
  Float_t phi[NUM_BX];
  int isTracker1[NUM_BX];
  int isTracker2[NUM_BX];
  Float_t pt1[NUM_BX];
  Float_t pt2[NUM_BX];
  Float_t eta1[NUM_BX];
  Float_t eta2[NUM_BX];
  Float_t p1[NUM_BX];
  Float_t p2[NUM_BX];
  Float_t phi1[NUM_BX];
  Float_t phi2[NUM_BX];
  int id1[NUM_BX];
  int id2[NUM_BX];
  Bool_t outerTrackisNonnull1[NUM_BX];
  Bool_t outerTrackisNonnull2[NUM_BX];
  int gen[NUM_BX];
  Bool_t isTrackerMuArbitrated1[NUM_BX];
  Bool_t isTrackerMuArbitrated2[NUM_BX];
  Bool_t isTMOneStationTight1[NUM_BX];
  Bool_t isTMOneStationTight2[NUM_BX];
  Int_t isCalo1[NUM_BX];
  Int_t isCalo2[NUM_BX];
  

  ntuple->SetBranchAddress("Run",&Run);
  ntuple->SetBranchAddress("Event",&Event);
  ntuple->SetBranchAddress("size",&size);
  ntuple->SetBranchAddress("mass",mass);
  ntuple->SetBranchAddress("pt",pt);
  ntuple->SetBranchAddress("genpt",genpt);
  ntuple->SetBranchAddress("eta",eta);
  ntuple->SetBranchAddress("y",y);
  ntuple->SetBranchAddress("phi",phi);
  ntuple->SetBranchAddress("isTracker1",isTracker1);
  ntuple->SetBranchAddress("isTracker2",isTracker2);
  ntuple->SetBranchAddress("pt1",pt1);
  ntuple->SetBranchAddress("pt2",pt2);
  ntuple->SetBranchAddress("eta1",eta1);
  ntuple->SetBranchAddress("eta2",eta2);
  ntuple->SetBranchAddress("p1",p1);
  ntuple->SetBranchAddress("p2",p2);
  ntuple->SetBranchAddress("phi1",phi1);
  ntuple->SetBranchAddress("phi2",phi2);
  ntuple->SetBranchAddress("id1",id1);
  ntuple->SetBranchAddress("id2",id2);
  ntuple->SetBranchAddress("outerTrackisNonnull1",outerTrackisNonnull1);
  ntuple->SetBranchAddress("outerTrackisNonnull2",outerTrackisNonnull2);
  ntuple->SetBranchAddress("isTrackerMuArbitrated1",isTrackerMuArbitrated1);
  ntuple->SetBranchAddress("isTrackerMuArbitrated2",isTrackerMuArbitrated2);
  ntuple->SetBranchAddress("isTMOneStationTight1",isTMOneStationTight1);
  ntuple->SetBranchAddress("isTMOneStationTight2",isTMOneStationTight2);
  ntuple->SetBranchAddress("gen",gen);
  ntuple->SetBranchAddress("isTriggered1",isTriggered1);
  ntuple->SetBranchAddress("isTriggered2",isTriggered2);
  ntuple->SetBranchAddress("isCalo1",isCalo1);
  ntuple->SetBranchAddress("isCalo2",isCalo2);

  Bool_t qualitycut1;
  Bool_t qualitycut2;
  Bool_t glb_cut1;
  Bool_t glb_cut2;
  Bool_t isTag1;
  Bool_t isTag2;
  Bool_t isacceptance1;
  Bool_t isacceptance2;

  Int_t entries = (Int_t)ntuple->GetEntries();
  
  for (int i=0; i<entries; i++){
    if (i%10000==0) cout <<i<<" / "<<entries<<endl;
    ntuple->GetEntry(i);

    for(int j=0;j<size;j++){
    
      qualitycut1 = false;
      qualitycut2 = false;
      glb_cut1 = false;
      glb_cut2 = false;
      isTag1=false;
      isTag2=false;
      isacceptance1=false;
      isacceptance2=false;

      if(id1[j]) qualitycut1 = true;
      if(id2[j]) qualitycut2 = true;
      if(isTrackerMuArbitrated1[j]&&isTMOneStationTight1[j]) glb_cut1 = true;
      if(isTrackerMuArbitrated2[j]&&isTMOneStationTight2[j]) glb_cut2 = true;
      isacceptance1=IsMuonInAcceptance(pt1[j],p1[j],eta1[j]);
      isacceptance2=IsMuonInAcceptance(pt2[j],p2[j],eta2[j]);
      
      isTag1=IsTag(isacceptance1,isTriggered1[j],qualitycut1,glb_cut1,isTracker1[j]);
      isTag2=IsTag(isacceptance2,isTriggered2[j],qualitycut2,glb_cut2,isTracker2[j]);

      if(isTag1||isTag2){
        if(isTag1){
          
          for(int m = 0; m < nMuPtBin; m++){
            if(pt2[j] > MuPtBin[m] && pt2[j] < MuPtBin[m+1]){
              //tracking efficiency
//              if(outerTrackisNonnull2[j]){
              if(outerTrackisNonnull2[j] && isacceptance2){
                hTrkPtAll[m]->Fill(mass[j]);
                hEtaTrkPtAll[m]->Fill(eta2[j]);
                massTrk=mass[j]; ptTrk = pt2[j]; etaTrk = eta2[j]; passTrk = 0; genTrk = gen[j];
                if(qualitycut2&&isTracker2[j]) {hTrkPtPass[m]->Fill(mass[j]); passTrk=1; hEtaTrkPtPass[m]->Fill(eta2[j]);}
                else { hTrkPtFail[m]->Fill(mass[j]); hTrkPtFail[m]->Fill(eta2[j]);}
                TnPtreeTrk->Fill();
              }//tracking efficinecy over

              //muon ID efficiency
//              if(qualitycut2&&isCalo2[j]&&isacceptance2){
              if(qualitycut2&&isTracker2[j]&&isacceptance2){
                hMuIdPtAll[m]->Fill(mass[j]);
                hEtaMuIdPtAll[m]->Fill(eta2[j]);
                massID=mass[j]; ptID = pt2[j]; etaID = eta2[j]; passID = 0; genID = gen[j];
//                if(isTracker2[j]&&isacceptance2){ hMuIdPtPass[m]->Fill(mass[j]); passID=1; hEtaMuIdPtPass[m]->Fill(eta2[j]);}
                if(glb_cut2){ hMuIdPtPass[m]->Fill(mass[j]); passID=1; hEtaMuIdPtPass[m]->Fill(eta2[j]);}
                else {hMuIdPtFail[m]->Fill(mass[j]); hEtaMuIdPtFail[m]->Fill(eta2[j]);}
                TnPtreeID->Fill();
              }//muon ID efficiency over

              //trigger efficiency
              if(qualitycut2&&glb_cut2&&isTracker2[j]&&isacceptance2){
                hTrigPtAll[m]->Fill(mass[j]);
                hEtaTrigPtAll[m]->Fill(eta2[j]);
                massTrg=mass[j]; ptTrg = pt2[j]; etaTrg = eta2[j]; passTrg = 0; genTrg = gen[j];
                if(isTriggered2[j]) {hTrigPtPass[m]->Fill(mass[j]); passTrg = 1; hEtaTrigPtPass[m]->Fill(eta2[j]);}
                else { hTrigPtFail[m]->Fill(mass[j]); hEtaTrigPtFail[m]->Fill(eta2[j]);}
                TnPtreeTrg->Fill();
              }//trigger efficiency over
            }//if proper pt bin
          }//loop over pt bins
        }//if isTag1
        else{
          
          for(int m = 0; m < nMuPtBin; m++){
            if(pt1[j] > MuPtBin[m] && pt1[j] < MuPtBin[m+1]){
              //tracking efficiency
//              if(outerTrackisNonnull1[j]){
              if(outerTrackisNonnull1[j] && isacceptance1){
                massTrk=mass[j]; ptTrk = pt1[j]; etaTrk = eta1[j]; passTrk = 0; genTrk = gen[j];
                hTrkPtAll[m]->Fill(mass[j]);
                hEtaTrkPtAll[m]->Fill(eta1[j]);
                if(qualitycut1&&isTracker1[j]) {hTrkPtPass[m]->Fill(mass[j]); passTrk = 1; hEtaTrkPtPass[m]->Fill(eta1[j]);}
                else { hTrkPtFail[m]->Fill(mass[j]); hEtaTrkPtFail[m]->Fill(eta1[j]);}
                TnPtreeTrk->Fill();
              }//tracking efficinecy over

              //muon ID efficiency
//              if(qualitycut1&&isCalo1[j]&&isacceptance1){
              if(qualitycut1&&isTracker1[j]&&isacceptance1){
                hMuIdPtAll[m]->Fill(mass[j]);
                hEtaMuIdPtAll[m]->Fill(eta1[j]);
                massID=mass[j]; ptID = pt1[j]; etaID = eta1[j]; passID = 0; genID = gen[j];
//                if(isTracker1[j]&&isacceptance1) {hMuIdPtPass[m]->Fill(mass[j]); passID = 1; hEtaMuIdPtPass[m]->Fill(eta1[j]);}
                if(glb_cut1) {hMuIdPtPass[m]->Fill(mass[j]); passID = 1; hEtaMuIdPtPass[m]->Fill(eta1[j]);}
                else {hMuIdPtFail[m]->Fill(mass[j]); hEtaMuIdPtFail[m]->Fill(eta1[j]);}
                TnPtreeID->Fill();
              }//muon ID efficiency over

              //trigger efficiency
              if(qualitycut1&&glb_cut1&&isTracker1[j]&&isacceptance1){
                hTrigPtAll[m]->Fill(mass[j]);
                hEtaTrigPtAll[m]->Fill(eta1[j]);
                massTrg=mass[j]; ptTrg = pt1[j]; etaTrg = eta1[j]; passTrg = 0; genTrg = gen[j];
                if(isTriggered1[j]) {hTrigPtPass[m]->Fill(mass[j]); passTrg = 1; hEtaTrigPtPass[m]->Fill(eta1[j]);}
                else { hTrigPtFail[m]->Fill(mass[j]); hEtaTrigPtFail[m]->Fill(eta1[j]);}
                TnPtreeTrg->Fill();
              }//trigger efficiency over
            }//if proper pt bin
          }//loop over pt bins
        }//if isTag2
      }//if isTag1||isTag2
      
      if(isTag1||isTag2){
        if(isTag1){
          
          for(int m = 0; m < nMuEtaBin; m++){
            if(eta2[j] > MuEtaBin[m] && eta2[j] < MuEtaBin[m+1]){
              //tracking efficiency
//              if(outerTrackisNonnull2[j]){
              if(outerTrackisNonnull2[j] && isacceptance2){
                hTrkEtaAll[m]->Fill(mass[j]);
                massTrk=mass[j]; etaTrk = eta2[j]; etaTrk = eta2[j]; passTrk = 0; genTrk = gen[j];
                if(qualitycut2&&isTracker2[j]) {hTrkEtaPass[m]->Fill(mass[j]); passTrk=1;}
                else hTrkEtaFail[m]->Fill(mass[j]);
                TnPtreeTrk->Fill();
              }//tracking efficinecy over

              //muon ID efficiency
//              if(qualitycut2&&isCalo2[j]){
              if(qualitycut2&&isTracker2[j]&&isacceptance2){
                hMuIdEtaAll[m]->Fill(mass[j]);
                massID=mass[j]; etaID = eta2[j]; etaID = eta2[j]; passID = 0; genID = gen[j];
//                if(isTracker2[j]&&isacceptance2){ hMuIdEtaPass[m]->Fill(mass[j]); passID=1; }
                if(glb_cut2){ hMuIdEtaPass[m]->Fill(mass[j]); passID=1; }
                else hMuIdEtaFail[m]->Fill(mass[j]);
                TnPtreeID->Fill();
              }//muon ID efficiency over

              //trigger efficiency
              if(qualitycut2&&glb_cut2&&isTracker2[j]&&isacceptance2){
                hTrigEtaAll[m]->Fill(mass[j]);
                massTrg=mass[j]; etaTrg = eta2[j]; etaTrg = eta2[j]; passTrg = 0; genTrg = gen[j];
                if(isTriggered2[j]) {hTrigEtaPass[m]->Fill(mass[j]); passTrg = 1;}
                else hTrigEtaFail[m]->Fill(mass[j]);
                TnPtreeTrg->Fill();
              }//trigger efficiency over
            }//if proper eta bin
          }//loop over eta bins
        }//if isTag1
        else{
          
          for(int m = 0; m < nMuEtaBin; m++){
            if(eta1[j] > MuEtaBin[m] && eta1[j] < MuEtaBin[m+1]){
              //tracking efficiency
 //             if(outerTrackisNonnull1[j]){
              if(outerTrackisNonnull1[j] && isacceptance1){
                massTrk=mass[j]; etaTrk = eta1[j]; etaTrk = eta1[j]; passTrk = 0; genTrk = gen[j];
                hTrkEtaAll[m]->Fill(mass[j]);
                if(qualitycut1&&isTracker1[j]) {hTrkEtaPass[m]->Fill(mass[j]); passTrk = 1;}
                else hTrkEtaFail[m]->Fill(mass[j]);
                TnPtreeTrk->Fill();
              }//tracking efficinecy over

              //muon ID efficiency
//              if(qualitycut1&&isCalo1[j]){
              if(qualitycut1&&isTracker1[j]&&isacceptance1){
                hMuIdEtaAll[m]->Fill(mass[j]);
                massID=mass[j]; etaID = eta1[j]; etaID = eta1[j]; passID = 0; genID = gen[j];
//                if(isTracker1[j]&&isacceptance1) {hMuIdEtaPass[m]->Fill(mass[j]); passID = 1;}
                if(glb_cut1) {hMuIdEtaPass[m]->Fill(mass[j]); passID = 1;}
                else hMuIdEtaFail[m]->Fill(mass[j]);
                TnPtreeID->Fill();
              }//muon ID efficiency over

              //trigger efficiency
              if(qualitycut1&&glb_cut1&&isTracker1[j]&&isacceptance1){
                hTrigEtaAll[m]->Fill(mass[j]);
                massTrg=mass[j]; etaTrg = eta1[j]; etaTrg = eta1[j]; passTrg = 0; genTrg = gen[j];
                if(isTriggered1[j]) {hTrigEtaPass[m]->Fill(mass[j]); passTrg = 1;}
                else hTrigEtaFail[m]->Fill(mass[j]);
                TnPtreeTrg->Fill();
              }//trigger efficiency over
            }//if proper eta bin
          }//loop over eta bins
        }//if isTag2
      }//if isTag1||isTag2

    }//loop over candidates  
  }// loop over events
  
  foutput->cd();
  TCanvas *cforSave= new TCanvas("cforSave","",600,600);
  cforSave->cd();
  TLegend *legforSave = myLegend(0.77,0.50,0.90,0.90);
  for(int i = 0; i < nMuPtBin; i++){
    hTrigPtPass[i]->Write();
    hTrigPtFail[i]->Write();
    hTrigPtAll[i]->Write();
    hTrkPtPass[i]->Write();
    hTrkPtFail[i]->Write();
    hTrkPtAll[i]->Write();
    hMuIdPtPass[i]->Write();
    hMuIdPtFail[i]->Write();
    hMuIdPtAll[i]->Write();

    hEtaTrigPtPass[i]->Write();
    hEtaTrigPtFail[i]->Write();
    hEtaTrigPtAll[i]->Write();
    hEtaTrkPtPass[i]->Write();
    hEtaTrkPtFail[i]->Write();
    hEtaTrkPtAll[i]->Write();
    hEtaMuIdPtPass[i]->Write();
    hEtaMuIdPtFail[i]->Write();
    hEtaMuIdPtAll[i]->Write();

    hEtaTrigPtAll[i]->Draw();
    hEtaTrigPtPass[i]->SetLineColor(3);
    hEtaTrigPtPass[i]->Draw("same");
    hEtaTrigPtFail[i]->SetLineColor(2);
    hEtaTrigPtFail[i]->Draw("same");
    legforSave->Clear(); legforSave->SetTextSize(0.03); legforSave->SetHeader(Form("pt: %.1f ~ %.1f", MuPtBin[i], MuPtBin[i+1]));
    legforSave->AddEntry(hEtaTrigPtAll[i],  "all", "l");
    legforSave->AddEntry(hEtaTrigPtPass[i], "pass", "l");
    legforSave->AddEntry(hEtaTrigPtFail[i], "fail", "l");
    legforSave->Draw("same");
    cforSave->SaveAs(Form("%s/hEtaTrigPtAll_%d.pdf",plotfolder.Data(), i));

    hEtaTrkPtAll[i]->Draw();
    hEtaTrkPtPass[i]->SetLineColor(3);
    hEtaTrkPtPass[i]->Draw("same");
    hEtaTrkPtFail[i]->SetLineColor(2);
    hEtaTrkPtFail[i]->Draw("same");
    legforSave->Clear(); legforSave->SetTextSize(0.03); legforSave->SetHeader(Form("pt: %.1f ~ %.1f", MuPtBin[i], MuPtBin[i+1]));
    legforSave->AddEntry(hEtaTrkPtAll[i],  "all", "l");
    legforSave->AddEntry(hEtaTrkPtPass[i], "pass", "l");
    legforSave->AddEntry(hEtaTrkPtFail[i], "fail", "l");
    legforSave->Draw("same");
    cforSave->SaveAs(Form("%s/hEtaTrkPtAll_%d.pdf",plotfolder.Data(), i));

    hEtaMuIdPtAll[i]->Draw();
    hEtaMuIdPtPass[i]->SetLineColor(3);
    hEtaMuIdPtPass[i]->Draw("same");
    hEtaMuIdPtFail[i]->SetLineColor(2);
    hEtaMuIdPtFail[i]->Draw("same");
    legforSave->Clear(); legforSave->SetTextSize(0.03); legforSave->SetHeader(Form("pt: %.1f ~ %.1f", MuPtBin[i], MuPtBin[i+1]));
    legforSave->AddEntry(hEtaMuIdPtAll[i],  "all", "l");
    legforSave->AddEntry(hEtaMuIdPtPass[i], "pass", "l");
    legforSave->AddEntry(hEtaMuIdPtFail[i], "fail", "l");
    legforSave->Draw("same");
    cforSave->SaveAs(Form("%s/hEtaMuIdPtAll_%d.pdf",plotfolder.Data(), i));

	//Do simple sideband calcualtion for fast check
	double nall, npass;
//	double sideband_width = 0.2;
	double sideband_width = 0.15;
    int siglow = -1;
    int sighigh = -1;
    for(int b = 0; b < nBin; b++){
     if(siglow == -1 && hTrigPtAll[i]->GetBinCenter(b+1) > JPSI_MASS-sideband_width) siglow = b+1;
     if(hTrigPtAll[i]->GetBinCenter(b+1) < JPSI_MASS+sideband_width) sighigh = b+1;
    }
    int sigwidth = int((sighigh-siglow+1)/2);
    int sblow = siglow-sigwidth;
    int sbhigh = sighigh+sigwidth;
    if(sblow < 0) sblow = 0;
    if(sbhigh > nBin) sbhigh = nBin;
    nall = hTrigPtAll[i]->Integral(siglow,sighigh)-(hTrigPtAll[i]->Integral(siglow-sigwidth, siglow-1)+hTrigPtAll[i]->Integral(sighigh+1,sighigh+sigwidth));
    npass = hTrigPtPass[i]->Integral(siglow,sighigh)-(hTrigPtPass[i]->Integral(siglow-sigwidth, siglow-1)+hTrigPtPass[i]->Integral(sighigh+1,sighigh+sigwidth));
    EffTrig->SetBinContent(i+1, npass/nall);
    EffTrig->SetBinError(i+1, 0.00001);
//    cout<<"his trg all : "<<hTrigPtAll[i]->Integral()<<endl;
//    cout<<"his trg pass: "<<hTrigPtPass[i]->Integral()<<endl;
//    cout<<"trg all : "<<nall<<endl;
//    cout<<"trg pass: "<<npass<<endl;
//    cout<<"eff trg: "<<npass/nall<<endl;

    //Using sideband subtraction
    nall = hTrkPtAll[i]->Integral(siglow,sighigh)-(hTrkPtAll[i]->Integral(siglow-sigwidth, siglow-1)+hTrkPtAll[i]->Integral(sighigh+1,sighigh+sigwidth));
    npass = hTrkPtPass[i]->Integral(siglow,sighigh)-(hTrkPtPass[i]->Integral(siglow-sigwidth, siglow-1)+hTrkPtPass[i]->Integral(sighigh+1,sighigh+sigwidth));
    EffTrk->SetBinContent(i+1, npass/nall);
    EffTrk->SetBinError(i+1, 0.00001);
//    cout<<"his trk all : "<<hTrkPtAll[i]->Integral()<<endl;
//    cout<<"his trk pass: "<<hTrkPtPass[i]->Integral()<<endl;
//    cout<<"trk all : "<<nall<<endl;
//    cout<<"trk pass: "<<npass<<endl;
//    cout<<"eff trk: "<<npass/nall<<endl;

    //Using sideband subtraction
    nall = hMuIdPtAll[i]->Integral(siglow,sighigh)-(hMuIdPtAll[i]->Integral(siglow-sigwidth, siglow-1)+hMuIdPtAll[i]->Integral(sighigh+1,sighigh+sigwidth));
    npass = hMuIdPtPass[i]->Integral(siglow,sighigh)-(hMuIdPtPass[i]->Integral(siglow-sigwidth, siglow-1)+hMuIdPtPass[i]->Integral(sighigh+1,sighigh+sigwidth));
    EffMuId->SetBinContent(i+1, npass/nall);
    EffMuId->SetBinError(i+1, 0.00001);
//    cout<<"his id all : "<<hMuIdPtAll[i]->Integral()<<endl;
//    cout<<"his id pass: "<<hMuIdPtPass[i]->Integral()<<endl;
//    cout<<"id all : "<<nall<<endl;
//    cout<<"id pass: "<<npass<<endl;
//    cout<<"eff id: "<<npass/nall<<endl;
  }
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


  for(int i = 0; i < nMuEtaBin; i++){
    hTrigEtaPass[i]->Write();
    hTrigEtaFail[i]->Write();
    hTrigEtaAll[i]->Write();
    hTrkEtaPass[i]->Write();
    hTrkEtaFail[i]->Write();
    hTrkEtaAll[i]->Write();
    hMuIdEtaPass[i]->Write();
    hMuIdEtaFail[i]->Write();
    hMuIdEtaAll[i]->Write();
  }
  TnPtreeTrg->Write();
  TnPtreeTrk->Write();
  TnPtreeID->Write();

  foutput->Close();
  delete foutput;

}

Bool_t IsMuonInAcceptance(Float_t pt,Float_t p,Float_t eta){
  Bool_t isselected=false;
  isselected=(abs(eta)<1.3&&pt>3.3)||(abs(eta)>1.3&&abs(eta)<2.2&&p>2.9)||(abs(eta)>2.2&&abs(eta)<2.4&&pt>0.8);
  return isselected;
}
Bool_t IsTag(Bool_t isacceptance, Int_t istrigger, Bool_t qualitycut, Bool_t glb_cut, Bool_t istracker)
{
  Bool_t isTag=false;
  isTag=(isacceptance && istrigger && qualitycut && glb_cut && istracker);
  return isTag;
}
