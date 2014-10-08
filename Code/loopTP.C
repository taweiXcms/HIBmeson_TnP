//T&P version of loop.C
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include <TMath.h>
#include <iostream>
#include <TNtuple.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <cmath>
#include "loopTP.h"

#define MUON_MASS   0.10565837
#define PION_MASS   0.13957018
#define KAON_MASS   0.493677
#define KSHORT_MASS 0.497614
#define KSTAR_MASS  0.89594
#define PHI_MASS    1.019455
#define JPSI_MASS   3.096916
static const unsigned int CaloMuon =  1<<4;

int signalGen(int Btype, int j)
{
  float BId,MId,tk1Id,tk2Id;
  int twoTks;
  //tk1:positive, tk2:negtive
  if(Btype==1)
    {
      BId = 521;//B+-
      MId = -1;
      tk1Id = 321;//K+-
      tk2Id = -1;
      twoTks = 0;
    }
  if(Btype==2)
    {
      BId = 521;//B+-
      MId = -1;
      tk1Id = 211;//pi+-
      tk2Id = -1;
      twoTks = 0;
    }
  if(Btype==3)
    {
      BId = 511;//B0
      MId = 310;//Ks
      tk1Id = 211;//pi+
      tk2Id = -211;//pi-
      twoTks = 1;
    }
  if(Btype==4)
    {
      BId = 511;//B0
      MId = 313;//K*0
      tk1Id = 321;//K+
      tk2Id = -211;//pi-
      twoTks = 1;
    }
  if(Btype==5)
    {
      BId = 511;//B0
      MId = 313;//K*0
      tk1Id = -321;//pi+
      tk2Id = 211;//K-
      twoTks = 1;
    }
  if(Btype==6)
    {
      BId = 531;//Bs
      MId = 333;//phi
      tk1Id = 321;//K+
      tk2Id = -321;//K-
      twoTks = 1;
    }

  int flag=0;
  if (abs(GenInfo_pdgId[j])==BId&&GenInfo_nDa[j]==2&&GenInfo_da1[j]!=-1&&GenInfo_da2[j]!=-1)
    {
      if (abs(GenInfo_pdgId[GenInfo_da1[j]]==443))//jpsi
	{
	  if(GenInfo_da1[GenInfo_da1[j]]!=-1&&GenInfo_da2[GenInfo_da1[j]]!=-1)
	    {
	      if(abs(GenInfo_pdgId[GenInfo_da1[GenInfo_da1[j]]])==13&&abs(GenInfo_pdgId[GenInfo_da2[GenInfo_da1[j]]])==13)
		{
		  if(!twoTks)
		    {
		      if(abs(GenInfo_pdgId[GenInfo_da2[j]])==tk1Id) flag++;
		    }
		  else
		    {
		      if (abs(GenInfo_pdgId[GenInfo_da2[j]])==MId) 
			{
			  if(GenInfo_da1[GenInfo_da2[j]]!=-1 && GenInfo_da2[GenInfo_da2[j]]!=-1)
			    {
			      if(GenInfo_pdgId[GenInfo_da1[GenInfo_da2[j]]]==tk1Id && GenInfo_pdgId[GenInfo_da2[GenInfo_da2[j]]]==tk2Id) flag++;
			    }
			}
		    }
		}
	    }
	}
    }
  return flag;
}

void loopTP(string infile="", string outfile="", bool REAL=1,bool PbpMC=0,int startEntries=0,int nEntries=0){
//void loopTP(string infile="", string outfile="", bool REAL=0,bool PbpMC=0,int startEntries=0,int nEntries=0){
//////////////////////////////////////////////////////////Phi
//   This file has been automatically generated 
//     (Thu Nov 21 13:34:42 2013 by ROOT version5.27/06b)
//   from TTree root/root
//   found on file: merged_pPbData_20131114.root
//////////////////////////////////////////////////////////

  const char* infname;
  const char* outfname;
  if(REAL){
//    infile="/mnt/hadoop/cms/store/user/tawei/Bfinder/Bfinder_20140930_PAMuon_HIRun2013_28Sep2013_v1/Bfinder_all_*.root";
//    outfile="/export/d00/scratch/jwang/tawei/TnPBntuple/TnPnt_20140930_PAMuon_HIRun2013_28Sep2013_v1.rootT";
    infile="/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_*.root";
//    outfile="/export/d00/scratch/jwang/tawei/TnPBntuple/TnPnt_20140930_PAMuon_HIRun2013_PromptReco_v1_1.rootT";
//    outfile="/export/d00/scratch/jwang/tawei/TnPBntuple/TnPnt_20140930_PAMuon_HIRun2013_PromptReco_v1_2.rootT";
    outfile="/export/d00/scratch/jwang/tawei/TnPBntuple/TnPnt_20140930_PAMuon_HIRun2013_PromptReco_v1_3.rootT";
//    outfile="/export/d00/scratch/jwang/tawei/TnPBntuple/TnPnt_20140930_PAMuon_HIRun2013_PromptReco_v1_4.rootT";
//    outfile="/export/d00/scratch/jwang/tawei/TnPBntuple/TnPnt_20140930_PAMuon_HIRun2013_PromptReco_v1_5.rootT";
//    outfile="/export/d00/scratch/jwang/tawei/TnPBntuple/TnPnt_20140930_PAMuon_HIRun2013_PromptReco_v1_6.rootT";
//    outfile="test.root";
  }
  else {
    infile="/mnt/hadoop/cms/store/user/tawei/Bfinder/BfinderBoostedMC_20140930_Kp/*.root";
    outfile="/export/d00/scratch/jwang/tawei/TnPBntuple/TnPnt_BoostedMC_20140930_Kp.rootT";
  }

  if(REAL) cout<<"--- REAL DATA ---"<<endl;
  else {
     cout<<"--- MC ---"<<endl;
     if(PbpMC) cout<<"--- Pbp ---"<<endl;
     else cout<<"--- pPb ---"<<endl;
  }
  infname = infile.c_str();
  outfname = outfile.c_str();

  //TFile *f = new TFile(infname);
  //TTree *root = (TTree*)f->Get("demo/root");
  //TTree *hlt = (TTree*)f->Get("hltanalysis/HltTree");
  TChain *root = new TChain("demo/root");
  TChain *hlt = new TChain("hltanalysis/HltTree");
//  root->Add(infile.c_str());
//  hlt->Add(infile.c_str());
//////
//  root->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_?_?_???.root");
//  root->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_??_?_???.root");
//  root->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_1??_?_???.root");
//  root->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_2??_?_???.root");
//  root->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_3??_?_???.root");
//  hlt->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_?_?_???.root");
//  hlt->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_??_?_???.root");
//  hlt->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_1??_?_???.root");
//  hlt->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_2??_?_???.root");
//  hlt->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_3??_?_???.root");
//////
//  root->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_4??_?_???.root");
//  root->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_5??_?_???.root");
//  root->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_6??_?_???.root");
//  hlt->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_4??_?_???.root");
//  hlt->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_5??_?_???.root");
//  hlt->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_6??_?_???.root");
//////
  root->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_7??_?_???.root");
  root->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_8??_?_???.root");
  root->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_9??_?_???.root");
  hlt->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_7??_?_???.root");
  hlt->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_8??_?_???.root");
  hlt->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_9??_?_???.root");
//////
//  root->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_10??_?_???.root");
//  root->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_11??_?_???.root");
//  root->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_12??_?_???.root");
//  hlt->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_10??_?_???.root");
//  hlt->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_11??_?_???.root");
//  hlt->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_12??_?_???.root");
//////
//  root->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_13??_?_???.root");
//  root->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_14??_?_???.root");
//  root->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_15??_?_???.root");
//  hlt->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_13??_?_???.root");
//  hlt->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_14??_?_???.root");
//  hlt->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_15??_?_???.root");
//////
//  root->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_16??_?_???.root");
//  root->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_17??_?_???.root");
//  hlt->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_16??_?_???.root");
//  hlt->Add("/mnt/hadoop/cms/store/user/twang/HI_BfinderNtuple/20140930_PAMuon_HIRun2013_PromptReco_v1/Bfinder_all_17??_?_???.root");
//////
  if (root->GetEntries()!=hlt->GetEntries()) {
     cout <<"Inconsistent number of entries!!! "<<infile<<endl;
     cout <<"HLT tree: "<<hlt->GetEntries()<<endl;
     cout <<"Bfinder tree: "<<root->GetEntries()<<endl;
  }
  
  TFile *outf = new TFile(outfname,"recreate");
  setBranch(root);
  setHltBranch(hlt);
  TTree* nt0 = new TTree("ntJpsi","");
  buildBranch(nt0);
  TTree* ntGen = new TTree("ntGen","");
  buildGenBranch(ntGen);
  cout<<"--- Tree building finished ---"<<endl;
  
  Long64_t nentries = root->GetEntries();
  Long64_t nbytes = 0;
  TLorentzVector* b4P = new TLorentzVector;
  //TLorentzVector* b4Pout = new TLorentzVector;
  TLorentzVector bGen;
  int flagEvt=0;  
  int offsetHltTree=0;
  if (nEntries!=0) nentries=nEntries;
  TLorentzVector* mu1P4 = new TLorentzVector;
  TLorentzVector* mu2P4 = new TLorentzVector;
  TLorentzVector* JpsiP4 = new TLorentzVector;
//  startEntries=7772040;
//  startEntries=17191;
//  startEntries=0;
//  nentries = 0;
//cout<<nentries<<endl; return;
  for (Long64_t i=startEntries; i<nentries;i++) {
    nbytes += root->GetEntry(i);
    flagEvt=0;
    if (i%10000==0) cout <<i<<" / "<<nentries<<"   offset HLT:"<<offsetHltTree<<endl;
    size=0;
    for (int muIt1=0;muIt1<MuonInfo_size;muIt1++) {
      //get the second muon
      for (int muIt2=muIt1+1;muIt2<MuonInfo_size;muIt2++) {
        if(size >= MAX_XB) {
		  cout<<"Fatal!!! exceed allowed total jpsi cand. size."<<endl; 
          cout<<"EvtInfo_RunNo: "<<EvtInfo_RunNo<<" / EvtInfo_EvtNo: "<<EvtInfo_EvtNo<<endl;
          cout<<"MuonInfo_size: "<<MuonInfo_size<<endl;
		  break;}
        //if(MuonInfo_charge[muIt1]==MuonInfo_charge[muIt2]) continue;
        //if( (MuonInfo_charge[muIt1]*MuonInfo_charge[muIt2]!=-1) && (MuonInfo_StandAloneMuon_charge[muIt1]*MuonInfo_StandAloneMuon_charge[muIt2]!=-1)) continue;
        //if( (MuonInfo_charge[muIt1]*MuonInfo_charge[muIt2]!=-1) && (MuonInfo_StandAloneMuon_charge[muIt1] == 0) && (MuonInfo_StandAloneMuon_charge[muIt2] == 0)) continue;
        if( (MuonInfo_charge[muIt1]*MuonInfo_charge[muIt2]!=-1) && 
			(MuonInfo_charge[muIt1]*MuonInfo_StandAloneMuon_charge[muIt2]!=-1) && 
			(MuonInfo_charge[muIt2]*MuonInfo_StandAloneMuon_charge[muIt1]!=-1)) continue;
        int mu1,mu2;	 
        //mu1 is charged +1
		if(MuonInfo_charge[muIt1]>0) {
		  mu1=muIt1;
		  mu2=muIt2;
		} 
        else {
		  mu1=muIt2;
		  mu2=muIt1;
		}	 
        //get the first muon
	    //store jpsi info
        mu1P4->SetPtEtaPhiM(MuonInfo_pt[mu1],MuonInfo_eta[mu1],MuonInfo_phi[mu1],MUON_MASS);
        mu2P4->SetPtEtaPhiM(MuonInfo_pt[mu2],MuonInfo_eta[mu2],MuonInfo_phi[mu2],MUON_MASS);
        JpsiP4 ->SetPxPyPzE(mu1P4->Px()+mu2P4->Px(),mu1P4->Py()+mu2P4->Py(),mu1P4->Pz()+mu2P4->Pz(),mu1P4->E()+mu2P4->E());
        mass[size][0]=JpsiP4->Mag();
        pt[size][0]  =JpsiP4->Pt();
        eta[size][0] =JpsiP4->Eta();
        y[size][0]   =JpsiP4->Y();
        phi[size][0] =JpsiP4->Phi();
        //reco + StandAloneMuon jpsi
        if((MuonInfo_charge[mu1]*MuonInfo_StandAloneMuon_charge[mu2]==-1)){
          mu1P4->SetPtEtaPhiM(MuonInfo_pt[mu1],MuonInfo_eta[mu1],MuonInfo_phi[mu1],MUON_MASS);
          mu2P4->SetPtEtaPhiM(MuonInfo_StandAloneMuon_pt[mu2],MuonInfo_StandAloneMuon_eta[mu2],MuonInfo_StandAloneMuon_phi[mu2],MUON_MASS);
          JpsiP4 ->SetPxPyPzE(mu1P4->Px()+mu2P4->Px(),mu1P4->Py()+mu2P4->Py(),mu1P4->Pz()+mu2P4->Pz(),mu1P4->E()+mu2P4->E());
          mass[size][1]=JpsiP4->Mag();
          pt[size][1]  =JpsiP4->Pt();
          eta[size][1] =JpsiP4->Eta();
          y[size][1]   =JpsiP4->Y();
          phi[size][1] =JpsiP4->Phi();
		}
        //get the first muon
        else{
          mass[size][1]= -1;
          pt[size][1]  = -1;
          eta[size][1] = -1;
          y[size][1]   = -1;
          phi[size][1] = -1;
        }
        if((MuonInfo_charge[mu2]*MuonInfo_StandAloneMuon_charge[mu1]==-1)){
          mu1P4->SetPtEtaPhiM(MuonInfo_StandAloneMuon_pt[mu1],MuonInfo_StandAloneMuon_eta[mu1],MuonInfo_StandAloneMuon_phi[mu1],MUON_MASS);
          mu2P4->SetPtEtaPhiM(MuonInfo_pt[mu2],MuonInfo_eta[mu2],MuonInfo_phi[mu2],MUON_MASS);
          JpsiP4 ->SetPxPyPzE(mu1P4->Px()+mu2P4->Px(),mu1P4->Py()+mu2P4->Py(),mu1P4->Pz()+mu2P4->Pz(),mu1P4->E()+mu2P4->E());
          mass[size][2]=JpsiP4->Mag();
          pt[size][2]  =JpsiP4->Pt();
          eta[size][2] =JpsiP4->Eta();
          y[size][2]   =JpsiP4->Y();
          phi[size][2] =JpsiP4->Phi();
		}
       //get the first muon
        else{
          mass[size][2]= -1;
          pt[size][2]  = -1;
          eta[size][2] = -1;
          y[size][2]   = -1;
          phi[size][2] = -1;
        }
        //if (JpsiP4->Mag()>5) continue;
        if(mass[size][0]>5 && mass[size][1]>5 && mass[size][2]>5) continue;
        if(mass[size][0]<2 && mass[size][1]<2 && mass[size][2]<2) continue;
        

        //store StandAloneMuon info
	    isStandAloneMuon1[size] = MuonInfo_isStandAloneMuon[mu1];
  	    StandAloneMuon_charge1[size] = MuonInfo_StandAloneMuon_charge[mu1];
    	StandAloneMuon_pt1[size] = MuonInfo_StandAloneMuon_pt[mu1];
	    StandAloneMuon_eta1[size] = MuonInfo_StandAloneMuon_eta[mu1];
	    StandAloneMuon_phi1[size] = MuonInfo_StandAloneMuon_phi[mu1];
	    isStandAloneMuon2[size] = MuonInfo_isStandAloneMuon[mu2];
  	    StandAloneMuon_charge2[size] = MuonInfo_StandAloneMuon_charge[mu2];
    	StandAloneMuon_pt2[size] = MuonInfo_StandAloneMuon_pt[mu2];
	    StandAloneMuon_eta2[size] = MuonInfo_StandAloneMuon_eta[mu2];
	    StandAloneMuon_phi2[size] = MuonInfo_StandAloneMuon_phi[mu2];
	    //store muon info
	    charge1[size]=MuonInfo_charge[mu1];
	    charge2[size]=MuonInfo_charge[mu2];
	    isTracker1[size]=MuonInfo_isTrackerMuon[mu1];
	    isTracker2[size]=MuonInfo_isTrackerMuon[mu2];
	    isGlobal1[size]=MuonInfo_isGlobalMuon[mu1];
	    isGlobal2[size]=MuonInfo_isGlobalMuon[mu2];
	    isCalo1[size]=MuonInfo_type[mu1]&CaloMuon;
	    isCalo2[size]=MuonInfo_type[mu2]&CaloMuon;
	    isTriggered1[size]=MuonInfo_isTriggered[mu1];
	    isTriggered2[size]=MuonInfo_isTriggered[mu2];
	    pt1[size]=MuonInfo_pt[mu1];
	    pt2[size]=MuonInfo_pt[mu2];
	    eta1[size]=MuonInfo_eta[mu1];
	    eta2[size]=MuonInfo_eta[mu2];
	    phi1[size]=MuonInfo_phi[mu1];
	    phi2[size]=MuonInfo_phi[mu2];	  
        //get muon p
	    float mu1px,mu1py,mu1pz;
        mu1px = MuonInfo_pt[mu1]*cos(MuonInfo_phi[mu1]);
        mu1py = MuonInfo_pt[mu1]*sin(MuonInfo_phi[mu1]);
        mu1pz = MuonInfo_pt[mu1]*sinh(MuonInfo_phi[mu1]);
        b4P->SetXYZM(mu1px,mu1py,mu1pz,MUON_MASS);
        p1[size] = b4P->P();
        float mu2px,mu2py,mu2pz;
        mu2px = MuonInfo_pt[mu2]*cos(MuonInfo_phi[mu2]);
        mu2py = MuonInfo_pt[mu2]*sin(MuonInfo_phi[mu2]);
        mu2pz = MuonInfo_pt[mu2]*sinh(MuonInfo_phi[mu2]);
        b4P->SetXYZM(mu2px,mu2py,mu2pz,MUON_MASS);
        p2[size] = b4P->P();	  
        //StandAloneMuon p
        b4P->SetXYZM(MuonInfo_StandAloneMuon_pt[mu1]*cos(MuonInfo_StandAloneMuon_phi[mu1]), 
					 MuonInfo_StandAloneMuon_pt[mu1]*sin(MuonInfo_StandAloneMuon_phi[mu1]), 
					 MuonInfo_StandAloneMuon_pt[mu1]*sinh(MuonInfo_StandAloneMuon_phi[mu1]), 
					 MUON_MASS);
        StandAloneMuon_p1[size] = b4P->P();
        b4P->SetXYZM(MuonInfo_StandAloneMuon_pt[mu2]*cos(MuonInfo_StandAloneMuon_phi[mu2]), 
					 MuonInfo_StandAloneMuon_pt[mu2]*sin(MuonInfo_StandAloneMuon_phi[mu2]), 
					 MuonInfo_StandAloneMuon_pt[mu2]*sinh(MuonInfo_StandAloneMuon_phi[mu2]), 
					 MUON_MASS);
        StandAloneMuon_p2[size] = b4P->P();

	    outerTrackisNonnull1[size]=MuonInfo_outerTrackisNonnull[mu1];
	    outerTrackisNonnull2[size]=MuonInfo_outerTrackisNonnull[mu2];
        nPixel1[size] = MuonInfo_i_nPixelLayer[mu1];
        nPixel2[size] = MuonInfo_i_nPixelLayer[mu2];
        nTracker1[size] = MuonInfo_i_nStripLayer[mu1] + MuonInfo_i_nPixelLayer[mu1];
        nTracker2[size] = MuonInfo_i_nStripLayer[mu2] + MuonInfo_i_nPixelLayer[mu2];
        dxy1[size] = MuonInfo_dxyPV[mu1];
        dxy2[size] = MuonInfo_dxyPV[mu2];
        dz1[size] = MuonInfo_dzPV[mu1];
        dz2[size] = MuonInfo_dzPV[mu2];
        chisq1[size] = MuonInfo_normchi2[mu1];
	    chisq2[size] = MuonInfo_normchi2[mu2];
		innerTrackQuality1[size] = MuonInfo_innerTrackQuality[mu1];
		innerTrackQuality2[size] = MuonInfo_innerTrackQuality[mu2];

	    StandAloneMuon_dzPV1[size] = MuonInfo_StandAloneMuon_dzPV[mu1];
	    StandAloneMuon_dxyPV1[size] = MuonInfo_StandAloneMuon_dxyPV[mu1];
	    StandAloneMuon_dzPV2[size] = MuonInfo_StandAloneMuon_dzPV[mu2];
	    StandAloneMuon_dxyPV2[size] = MuonInfo_StandAloneMuon_dxyPV[mu2];
	    
        if(MuonInfo_muqual[mu1]&16) isTrackerMuArbitrated1[size] = 1;
        else isTrackerMuArbitrated1[size] = 0;
        if(MuonInfo_muqual[mu1]&4096) isTMOneStationTight1[size] = 1;
        else isTMOneStationTight1[mu1] = 0;
 
        if(MuonInfo_muqual[mu2]&16) isTrackerMuArbitrated2[size] = 1;
        else isTrackerMuArbitrated2[size] = 0;
        if(MuonInfo_muqual[mu2]&4096) isTMOneStationTight2[size] = 1;
        else isTMOneStationTight2[mu2] = 0;
	
        id1[size]=1;
  	    if(abs(MuonInfo_dxyPV[mu1])>=3. || abs(MuonInfo_dzPV[mu1])>=30.) id1[size]=0;
	    if(MuonInfo_i_nPixelLayer[mu1]<1.) id1[size]=0;
	    if(MuonInfo_normchi2[mu1]>1.8) id1[size]=0;
	    if((MuonInfo_i_nStripLayer[mu1]+MuonInfo_i_nPixelLayer[mu1])<6.) id1[size]=0;
        id2[size]=1;
  	    if(abs(MuonInfo_dxyPV[mu2])>=3. || abs(MuonInfo_dzPV[mu2])>=30.) id2[size]=0;
	    if(MuonInfo_i_nPixelLayer[mu2]<1.) id2[size]=0;
	    if(MuonInfo_normchi2[mu2]>1.8) id2[size]=0;
	    if((MuonInfo_i_nStripLayer[mu2]+MuonInfo_i_nPixelLayer[mu2])<6.) id2[size]=0;

        StandAloneMuon_id1[size]=1;
		if(!MuonInfo_isStandAloneMuon[mu1]) StandAloneMuon_id1[size]=0;
  	    if(abs(MuonInfo_StandAloneMuon_dxyPV[mu1])>=3. || abs(MuonInfo_StandAloneMuon_dzPV[mu1])>=30.) StandAloneMuon_id1[size]=0;
	    if(MuonInfo_i_nPixelLayer[mu1]<1.) StandAloneMuon_id1[size]=0;
	    if(MuonInfo_normchi2[mu1]>1.8) StandAloneMuon_id1[size]=0;
	    if((MuonInfo_i_nStripLayer[mu1]+MuonInfo_i_nPixelLayer[mu1])<6.) StandAloneMuon_id1[size]=0;
        StandAloneMuon_id2[size]=1;
		if(!MuonInfo_isStandAloneMuon[mu2]) StandAloneMuon_id2[size]=0;
  	    if(abs(MuonInfo_StandAloneMuon_dxyPV[mu2])>=3. || abs(MuonInfo_StandAloneMuon_dzPV[mu2])>=30.) StandAloneMuon_id2[size]=0;
	    if(MuonInfo_i_nPixelLayer[mu2]<1.) StandAloneMuon_id2[size]=0;
	    if(MuonInfo_normchi2[mu2]>1.8) StandAloneMuon_id2[size]=0;
	    if((MuonInfo_i_nStripLayer[mu2]+MuonInfo_i_nPixelLayer[mu2])<6.) StandAloneMuon_id2[size]=0;

        gen[size]=0;
    	genpt[size]=-1;
        int level1=0;
  	    int level2=0;
        if(MuonInfo_geninfo_index[mu1]>-1) {
	      if(abs(GenInfo_pdgId[MuonInfo_geninfo_index[mu1]])==13) {
	        level1=1;
	        if(GenInfo_mo1[MuonInfo_geninfo_index[mu1]]>-1) {
	          if(GenInfo_pdgId[GenInfo_mo1[MuonInfo_geninfo_index[mu1]]]==443) {
	            level1=2;
	          }
	        }
	        gen[size]+=(level1*1);
	      }
        }
	  
        if(MuonInfo_geninfo_index[mu2]>-1) {
	      if(abs(GenInfo_pdgId[MuonInfo_geninfo_index[mu2]])==13) {
	        level2=1;
	        if(GenInfo_mo1[MuonInfo_geninfo_index[mu2]]>-1) {
		      if(GenInfo_pdgId[GenInfo_mo1[MuonInfo_geninfo_index[mu2]]]==443) {
		        level2=2;
	  	      }
	        }
	        gen[size]+=(level2*10);
	      }
        }  
	  
	    if (level1==2&&level2==2) {
	      if (GenInfo_mo1[MuonInfo_geninfo_index[mu2]]==GenInfo_mo1[MuonInfo_geninfo_index[mu1]]) {
	        gen[size]+=100;
		    genpt[size]=GenInfo_pt[GenInfo_mo1[MuonInfo_geninfo_index[mu2]]];
	      }
	    }
        size++;
      }//mu2
    }//mu1
    nt0->Fill();

    if(!REAL){
	  Gensize = 0;
	  for (int j=0;j<GenInfo_size;j++){
	    if (GenInfo_pdgId[j]!=443) continue;
	    bGen.SetPtEtaPhiM(GenInfo_pt[j],GenInfo_eta[j],GenInfo_phi[j],GenInfo_mass[j]);
	    Geny[j] = bGen.Rapidity();
	    Geneta[j] = bGen.Eta();
	    Genphi[j] = bGen.Phi();
	    Genpt[j] = bGen.Pt();
	    GenpdgId[j] = GenInfo_pdgId[j];
	    Gensize++;
	  }
	  ntGen->Fill();
	}
  }//event
  outf->Write();
  outf->Close();
}//void
