#include "format.h"
#include "utilities.h"

using namespace std;

//const int nMuPtBin = 1;
const int nMuPtBin = 7;
//const int nMuPtBin = 10;
//const int nMuPtBin = 18;
//const double MuPtBin[nMuPtBin+1] = {0,30};
const double MuPtBin[nMuPtBin+1] = {0.0,1.5,3.0,4.5,6.0,9.0,20.0,30.0};
//const double MuPtBin[nMuPtBin+1] = {0.0, 1.5, 2.9, 3.3, 3.5, 4.0, 4.5, 6.0, 9.0, 20.0, 30.0};
//const double MuPtBin[nMuPtBin+1] = {0.0, 1.5, 2.5, 3.0, 3.1, 3.125, 3.15, 3.175, 3.2, 3.225, 3.25, 3.3, 3.5, 4.0, 4.5, 6.0, 9.0, 20.0, 30.0};
////TH1 setting of mumu
int nBin = 100;
double blow = 2.5;
double bhigh = 3.5;
double measure = (bhigh-blow)/nBin;
//width 2.9~3.3
double sideband_width = 0.2;

TString ninf = "../BfinderBoostedMC_20140701_inclBtoPsiMuMu_pa_STARTHI53_V27-v1_00000_reduce.root";
//TString ninf = "../BfinderBoostedMC_20140701_inclBtoPsiMuMu_pa_STARTHI53_V27-v1_20000_reduce.root";
//TString ninf = "../CaloAcc_BfinderBoostedMC_20140701_inclBtoPsiMuMu_pa_STARTHI53_V27-v1_00000_reduce.root";
//TString ninf = "../BfinderBoostedMC_20140708_inclBtoPsiMuMu_pa_STARTHI53_V27-v1_00000_reduce.root";
//TString ninf = "../BfinderBoostedMC_20140708_hckim-HIJINGemb_inclBtoPsiMuMu_reduced.root";
//TString ninf = "../20140707_PAMuon_HIRun2013_28Sep2013_v1_reduced.root";
TString noutf = "../TnPOutMC.root";
//TString noutf = "../TnPOut.root";
TString plotfolder = "PlotsSideBandMC";
//TString plotfolder = "PlotsSideBand";

//TTree variable
double massTrg, massTrk, massID;
double passTrg, passTrk, passID;
double ptTrg, ptTrk, ptID;
double etaTrg, etaTrk, etaID;

void doTnP(TH1D* hmumu_trg_all[], TH1D* hmumu_trg_pass[], TH1D* hmumu_trk_all[], TH1D* hmumu_trk_pass[], TH1D* hmumu_id_all[], TH1D* hmumu_id_pass[], 
TH1D* heta_trg_all[], TH1D* heta_trk_all[], TH1D* heta_id_all[], 
TH1D* heta_trg_pass[], TH1D* heta_trg_fail[], TH1D* heta_trk_pass[], TH1D* heta_trk_fail[], TH1D* heta_id_pass[], TH1D* heta_id_fail[], 
TTree* TnPtreeTrg, TTree* TnPtreeTrk, TTree* TnPtreeID){
  TChain *nt = new TChain("demo/root");
  //TChain *hlt = new TChain("hltanalysis/HltTree");
  TChain *hlt = new TChain("demo/HltTree");
  nt->Add(ninf);
  hlt->Add(ninf);

  nt->SetBranchStatus("*" ,0);
  nt->SetBranchStatus("MuonInfo*" ,1);
  hlt->SetBranchStatus("*" ,0);
  hlt->SetBranchStatus("HLT_PAMu3_v1" ,1);              
  nt->AddFriend(hlt);

  MuonInfoBranches MuonInfo;
  MuonInfo.setbranchadd(nt);

  Int_t HLT_PAMu3_v1; 
  nt->SetBranchAddress("HLT_PAMu3_v1",&HLT_PAMu3_v1);                         

  TLorentzVector mu1Vec;
  TLorentzVector mu2Vec;
  TLorentzVector jpsiVec;
  int nevents_total = nt->GetEntries();                                                  
//  nevents_total = 10000;
  for(int entry=0; entry<nevents_total; entry++){
    if ((entry%10000) == 0) printf("Loading event #%d of %d.\n",entry,nevents_total);
    nt->GetEntry(entry);
    if(MuonInfo.size<2) continue;
    if(!HLT_PAMu3_v1) continue;
	for(int mu1 = 0; mu1 < MuonInfo.size; mu1++){
	  for(int mu2 = mu1+1; mu2 < MuonInfo.size; mu2++){
		if(MuonInfo.charge[mu1]*MuonInfo.charge[mu2] > 0) continue;
        mu1Vec.SetPtEtaPhiM(MuonInfo.pt[mu1], MuonInfo.eta[mu1], MuonInfo.phi[mu1], MUON_MASS);
        mu2Vec.SetPtEtaPhiM(MuonInfo.pt[mu2], MuonInfo.eta[mu2], MuonInfo.phi[mu2], MUON_MASS);
		jpsiVec = mu1Vec+mu2Vec;
		////Trigger
		bool mu1Trg = false;
		bool mu2Trg = false;
		//std::cout<<MuTrgMatchTrgObjE->at(mu1).at(0)<<std::endl;
		if(MuonInfo.MuTrgMatchTrgObjE->at(mu1).at(0) > 0) mu1Trg = true;
		if(MuonInfo.MuTrgMatchTrgObjE->at(mu2).at(0) > 0) mu2Trg = true;
		////
		////Trigger efficiency
		int tagchg = 0;//1 = tag positive, -1 = tag negative, other = tag both
		//if(int(MuonInfo.pt[mu1]*10)%2==1) pickmu1 = true;
		if(MuonInfo.isTrackerMuon[mu1] && KisooTrackSel(MuonInfo, mu1) && KisooGlobalSel(MuonInfo, mu1) && Acc(MuonInfo, mu1) && mu1Trg && MuonInfo.charge[mu1] != -tagchg &&
		   MuonInfo.isTrackerMuon[mu2] && KisooTrackSel(MuonInfo, mu2) && KisooGlobalSel(MuonInfo, mu2) && Acc(MuonInfo, mu2) 
		){
		  //Fill Tree
		  massTrg = jpsiVec.Mag(); ptTrg = MuonInfo.pt[mu2]; etaTrg = MuonInfo.eta[mu2]; passTrg = 0;
          if(mu2Trg) {passTrg = 1;}
		  TnPtreeTrg->Fill();
		  //Fill Histogram
		  for(int i = 0; i < nMuPtBin; i++){
            if(MuonInfo.pt[mu2] > MuPtBin[i] && MuonInfo.pt[mu2] < MuPtBin[i+1]){
		      hmumu_trg_all[i]->Fill(jpsiVec.Mag());
			    if(mu2Trg){
                hmumu_trg_pass[i]->Fill(jpsiVec.Mag());
			    }
              break;
            }
          }
        }
		if(MuonInfo.isTrackerMuon[mu2] && KisooTrackSel(MuonInfo, mu2) && KisooGlobalSel(MuonInfo, mu2) && Acc(MuonInfo, mu2) && mu2Trg && MuonInfo.charge[mu2] != -tagchg &&
		   MuonInfo.isTrackerMuon[mu1] && KisooTrackSel(MuonInfo, mu1) && KisooGlobalSel(MuonInfo, mu1) && Acc(MuonInfo, mu1) 
		){
		  //Fill Tree
		  massTrg = jpsiVec.Mag(); ptTrg = MuonInfo.pt[mu1]; etaTrg = MuonInfo.eta[mu1]; passTrg = 0;
          if(mu1Trg) {passTrg = 1;}
		  TnPtreeTrg->Fill();
		  //Fill Histogram
		  for(int i = 0; i < nMuPtBin; i++){
            if(MuonInfo.pt[mu1] > MuPtBin[i] && MuonInfo.pt[mu1] < MuPtBin[i+1]){
		      hmumu_trg_all[i]->Fill(jpsiVec.Mag());
			    if(mu1Trg){
                hmumu_trg_pass[i]->Fill(jpsiVec.Mag());
			    }
              break;
            }
          }
        }
        ////
		////Tracking efficiency
		if(MuonInfo.isTrackerMuon[mu1] && KisooTrackSel(MuonInfo, mu1) && KisooGlobalSel(MuonInfo, mu1) && Acc(MuonInfo, mu1) && mu1Trg	&& MuonInfo.charge[mu1] != -tagchg &&
           MuonInfo.outerTrackisNonnull[mu2]
		){
		  //Fill Tree
		  massTrk = jpsiVec.Mag(); ptTrk = MuonInfo.pt[mu2]; etaTrk = MuonInfo.eta[mu2]; passTrk = 0;
          if(MuonInfo.isTrackerMuon[mu2] && KisooTrackSel(MuonInfo, mu2)) {passTrk = 1;}
		  TnPtreeTrk->Fill();
		  //Fill Histogram
		  for(int i = 0; i < nMuPtBin; i++){
            if(MuonInfo.pt[mu2] > MuPtBin[i] && MuonInfo.pt[mu2] < MuPtBin[i+1]){
		      hmumu_trk_all[i]->Fill(jpsiVec.Mag());
			    if(MuonInfo.isTrackerMuon[mu2] && KisooTrackSel(MuonInfo, mu2)){
                hmumu_trk_pass[i]->Fill(jpsiVec.Mag());
			    }
              break;
            }
          }
        }
		if(MuonInfo.isTrackerMuon[mu2] && KisooTrackSel(MuonInfo, mu2) && KisooGlobalSel(MuonInfo, mu2) && Acc(MuonInfo, mu2) && mu2Trg	&& MuonInfo.charge[mu2] != -tagchg &&
           MuonInfo.outerTrackisNonnull[mu1]  
		){
		  //Fill Tree
		  massTrk = jpsiVec.Mag(); ptTrk = MuonInfo.pt[mu1]; etaTrk = MuonInfo.eta[mu1]; passTrk = 0;
          if(MuonInfo.isTrackerMuon[mu1] && KisooTrackSel(MuonInfo, mu1)) {passTrk = 1;}
		  TnPtreeTrk->Fill();
		  //Fill Histogram
		  for(int i = 0; i < nMuPtBin; i++){
            if(MuonInfo.pt[mu1] > MuPtBin[i] && MuonInfo.pt[mu1] < MuPtBin[i+1]){
		      hmumu_trk_all[i]->Fill(jpsiVec.Mag());
			    if(MuonInfo.isTrackerMuon[mu1] && KisooTrackSel(MuonInfo, mu1)){
                hmumu_trk_pass[i]->Fill(jpsiVec.Mag());
			    }
              break;
            }
          }
        }
		////
		////ID efficiency
        if(MuonInfo.isTrackerMuon[mu1] && KisooTrackSel(MuonInfo, mu1) && KisooGlobalSel(MuonInfo, mu1) && Acc(MuonInfo, mu1) && mu1Trg && MuonInfo.charge[mu1] != -tagchg &&
           (MuonInfo.type[mu2] & (1<<4)) && KisooTrackSel(MuonInfo, mu2) && Acc(MuonInfo, mu2) 
//           (MuonInfo.type[mu2] & (1<<4)) && KisooTrackSel(MuonInfo, mu2) && Acc(MuonInfo, mu2) && jpsiVec.Pt()>3.0
//           (MuonInfo.type[mu2] & (1<<4)) && KisooTrackSel(MuonInfo, mu2) && Acc(MuonInfo, mu2) && MuonInfo.geninfo_index[mu2] != -1 
//           (MuonInfo.type[mu2] & (1<<4)) && KisooTrackSel(MuonInfo, mu2)
        ){ 
		  //Fill Tree
		  massID = jpsiVec.Mag(); ptID = MuonInfo.pt[mu2]; etaID = MuonInfo.eta[mu2]; passID = 0;
          if(MuonInfo.isTrackerMuon[mu2] && KisooTrackSel(MuonInfo, mu2) && KisooGlobalSel(MuonInfo, mu2) && Acc(MuonInfo, mu2)) {passID = 1;}
		  TnPtreeID->Fill();
		  //Fill Histogram
		  for(int i = 0; i < nMuPtBin; i++){
            if(MuonInfo.pt[mu2] > MuPtBin[i] && MuonInfo.pt[mu2] < MuPtBin[i+1]){
		      hmumu_id_all[i]->Fill(jpsiVec.Mag());
		      heta_id_all[i]->Fill(MuonInfo.eta[mu2]);
			    if(MuonInfo.isTrackerMuon[mu2] && KisooTrackSel(MuonInfo, mu2) && KisooGlobalSel(MuonInfo, mu2) && Acc(MuonInfo, mu2)){
//			    if(MuonInfo.isTrackerMuon[mu2] && KisooTrackSel(MuonInfo, mu2) && KisooGlobalSel(MuonInfo, mu2) && Acc(MuonInfo, mu2) && mu2Trg){
                  hmumu_id_pass[i]->Fill(jpsiVec.Mag());
                  heta_id_pass[i]->Fill(MuonInfo.eta[mu2]);
			    }
				else{
				  heta_id_fail[i]->Fill(MuonInfo.eta[mu2]);
				}
              break;
            }
          }
        }
        if(MuonInfo.isTrackerMuon[mu2] && KisooTrackSel(MuonInfo, mu2) && KisooGlobalSel(MuonInfo, mu2) && Acc(MuonInfo, mu2) && mu2Trg && MuonInfo.charge[mu2] != -tagchg &&
           (MuonInfo.type[mu1] & (1<<4)) && KisooTrackSel(MuonInfo, mu1) && Acc(MuonInfo, mu1) 
//           (MuonInfo.type[mu1] & (1<<4)) && KisooTrackSel(MuonInfo, mu1) && Acc(MuonInfo, mu1) && jpsiVec.Pt() > 3.0
//           (MuonInfo.type[mu1] & (1<<4)) && KisooTrackSel(MuonInfo, mu1) && Acc(MuonInfo, mu1) && MuonInfo.geninfo_index[mu1] != -1
//           (MuonInfo.type[mu1] & (1<<4)) && KisooTrackSel(MuonInfo, mu1)
        ){  
		  //Fill Tree
		  massID = jpsiVec.Mag(); ptID = MuonInfo.pt[mu1]; etaID = MuonInfo.eta[mu1]; passID = 0;
          if(MuonInfo.isTrackerMuon[mu1] && KisooTrackSel(MuonInfo, mu1) && KisooGlobalSel(MuonInfo, mu1) && Acc(MuonInfo, mu1)) {passID = 1;}
		  TnPtreeID->Fill();
		  //Fill Histogram
		  for(int i = 0; i < nMuPtBin; i++){
            if(MuonInfo.pt[mu1] > MuPtBin[i] && MuonInfo.pt[mu1] < MuPtBin[i+1]){
		      hmumu_id_all[i]->Fill(jpsiVec.Mag());
		      heta_id_all[i]->Fill(MuonInfo.eta[mu1]);
			    if(MuonInfo.isTrackerMuon[mu1] && KisooTrackSel(MuonInfo, mu1) && KisooGlobalSel(MuonInfo, mu1) && Acc(MuonInfo, mu1)){
//			    if(MuonInfo.isTrackerMuon[mu1] && KisooTrackSel(MuonInfo, mu1) && KisooGlobalSel(MuonInfo, mu1) && Acc(MuonInfo, mu1) && mu1Trg){
                  hmumu_id_pass[i]->Fill(jpsiVec.Mag());
                  heta_id_pass[i]->Fill(MuonInfo.eta[mu1]);
			    }
				else{
				  heta_id_fail[i]->Fill(MuonInfo.eta[mu1]);
				}
              break;
            }
          }
        }
	    ////
      }
    }
  }
}

int TnPeff(){
  gStyle->SetOptStat("0");
  TTree* TnPtreeTrg = new TTree("TnPtreeTrg","");
  TTree* TnPtreeTrk = new TTree("TnPtreeTrk","");
  TTree* TnPtreeID = new TTree("TnPtreeID","");
  TnPtreeTrg->Branch("mass",&massTrg);
  TnPtreeTrg->Branch("pt",&ptTrg);
  TnPtreeTrg->Branch("eta",&etaTrg);
  TnPtreeTrg->Branch("pass",&passTrg);

  TnPtreeTrk->Branch("mass",&massTrk);
  TnPtreeTrk->Branch("pt",&ptTrk);
  TnPtreeTrk->Branch("eta",&etaTrk);
  TnPtreeTrk->Branch("pass",&passTrk);

  TnPtreeID->Branch("mass",&massID);
  TnPtreeID->Branch("pt",&ptID);
  TnPtreeID->Branch("eta",&etaID);
  TnPtreeID->Branch("pass",&passID);
  
  TH1D* hmumu_trg_all[nMuPtBin];
  TH1D* hmumu_trg_pass[nMuPtBin];
  TH1D* hmumu_trg_fail[nMuPtBin];
  TH1D* hmumu_trk_all[nMuPtBin];
  TH1D* hmumu_trk_pass[nMuPtBin];
  TH1D* hmumu_trk_fail[nMuPtBin];
  TH1D* hmumu_id_all[nMuPtBin];
  TH1D* hmumu_id_pass[nMuPtBin];
  TH1D* hmumu_id_fail[nMuPtBin];

  TH1D* heta_trg_all[nMuPtBin];
  TH1D* heta_trk_all[nMuPtBin];
  TH1D* heta_id_all[nMuPtBin];
  TH1D* heta_trg_pass[nMuPtBin];
  TH1D* heta_trg_fail[nMuPtBin];
  TH1D* heta_trk_pass[nMuPtBin];
  TH1D* heta_trk_fail[nMuPtBin];
  TH1D* heta_id_pass[nMuPtBin];
  TH1D* heta_id_fail[nMuPtBin];

  TF1* tf1mumu_trg_all[nMuPtBin];
  TF1* tf1mumu_trg_pass[nMuPtBin];
  TF1* tf1mumu_trk_all[nMuPtBin];
  TF1* tf1mumu_trk_pass[nMuPtBin];
  TF1* tf1mumu_id_all[nMuPtBin];
  TF1* tf1mumu_id_pass[nMuPtBin];
  TH1D* eff_trg = myTH1D("eff_trg", "Mu Pt (GeV)", "efficiency", 2, nMuPtBin, MuPtBin);
  TH1D* eff_trk = myTH1D("eff_trk", "Mu Pt (GeV)", "efficiency", 3, nMuPtBin, MuPtBin);
  TH1D* eff_id = myTH1D("eff_id", "Mu Pt (GeV)", "efficiency", 4, nMuPtBin, MuPtBin);
  for(int i = 0; i < nMuPtBin; i++){
    hmumu_trg_all[i] = new TH1D(Form("hmumu_trg_all%d",i),"",nBin,blow,bhigh);
    hmumu_trg_pass[i] = new TH1D(Form("hmumu_trg_pass%d",i),"",nBin,blow,bhigh);
    hmumu_trg_fail[i] = new TH1D(Form("hmumu_trg_fail%d",i),"",nBin,blow,bhigh);
    hmumu_trk_all[i] = new TH1D(Form("hmumu_trk_all%d",i),"",nBin,blow,bhigh);
    hmumu_trk_pass[i] = new TH1D(Form("hmumu_trk_pass%d",i),"",nBin,blow,bhigh);
    hmumu_trk_fail[i] = new TH1D(Form("hmumu_trk_fail%d",i),"",nBin,blow,bhigh);
    hmumu_id_all[i] = new TH1D(Form("hmumu_id_all%d",i),"",nBin,blow,bhigh);
    hmumu_id_pass[i] = new TH1D(Form("hmumu_id_pass%d",i),"",nBin,blow,bhigh);
    hmumu_id_fail[i] = new TH1D(Form("hmumu_id_fail%d",i),"",nBin,blow,bhigh);

    heta_trg_all[i] = new TH1D(Form("heta_trg_all%d",i),"",nBin,-4,4);
    heta_trk_all[i] = new TH1D(Form("heta_trk_all%d",i),"",nBin,-4,4);
    heta_id_all[i]  = new TH1D(Form("heta_id_all%d", i),"",nBin,-4,4);
    heta_trg_pass[i] = new TH1D(Form("heta_trg_pass%d",i),"",nBin,-4,4);
    heta_trg_fail[i] = new TH1D(Form("heta_trg_fail%d",i),"",nBin,-4,4);
    heta_trk_pass[i] = new TH1D(Form("heta_trk_pass%d",i),"",nBin,-4,4);
    heta_trk_fail[i] = new TH1D(Form("heta_trk_fail%d",i),"",nBin,-4,4);
    heta_id_pass[i]  = new TH1D(Form("heta_id_pass%d", i),"",nBin,-4,4);
    heta_id_fail[i]  = new TH1D(Form("heta_id_fail%d", i),"",nBin,-4,4);

	tf1mumu_trg_all[i] = myTF1(Form("tf1_trg_all%d",i), blow, bhigh);
	tf1mumu_trg_pass[i] = myTF1(Form("tf1_trg_pass%d",i), blow, bhigh);
	tf1mumu_trk_all[i] = myTF1(Form("tf1_trk_all%d",i), blow, bhigh);
	tf1mumu_trk_pass[i] = myTF1(Form("tf1_trk_pass%d",i), blow, bhigh);
	tf1mumu_id_all[i] = myTF1(Form("tf1_id_all%d",i), blow, bhigh);
	tf1mumu_id_pass[i] = myTF1(Form("tf1_id_pass%d",i), blow, bhigh);
  }
  TCanvas *c= new TCanvas("c","",600,600);
  c->SetTopMargin(0.07504363);
  c->cd();
  doTnP(hmumu_trg_all, hmumu_trg_pass, hmumu_trk_all, hmumu_trk_pass, hmumu_id_all, hmumu_id_pass, 
        heta_trg_all, heta_trk_all, heta_id_all,
        heta_trg_pass, heta_trg_fail, heta_trk_pass, heta_trk_fail, heta_id_pass, heta_id_fail, 
        TnPtreeTrg, TnPtreeTrk, TnPtreeID);
  TF1 *signal = new TF1("signal","[0]*([4]*Gaus(x,[1],[2])/(sqrt(2*3.14159)*[2])+(1-[4])*Gaus(x,[1],[3])/(sqrt(2*3.14159)*[3]))");
  signal->SetRange(blow, bhigh);
  for(int i = 0; i < nMuPtBin; i++){
	cout<<"Mu Pt Bin(GeV): " <<MuPtBin[i]<<endl;
	////side band method
	int siglow = -1;
	int sighigh = -1;
	for(int b = 0; b < nBin; b++){
	 if(siglow == -1 && hmumu_trg_all[i]->GetBinCenter(b+1) > JPSI_MASS-sideband_width) siglow = b+1; 
	 if(hmumu_trg_all[i]->GetBinCenter(b+1) < JPSI_MASS+sideband_width) sighigh = b+1; 
	}
	int sigwidth = int((sighigh-siglow+1)/2);
	int sblow = siglow-sigwidth;
	int sbhigh = sighigh+sigwidth;
	if(sblow < 0) sblow = 0;
	if(sbhigh > nBin) sbhigh = nBin;
//	cout<<siglow-sigwidth<<endl;cout<<siglow-1<<endl;cout<<sighigh+1<<endl;cout<<sighigh+sigwidth<<endl;

    //Get failing probe histogram
    hmumu_trg_fail[i] = (TH1D*)hmumu_trg_all[i]->Clone();
    hmumu_trg_fail[i]->Add(hmumu_trg_pass[i], -1);
    hmumu_trk_fail[i] = (TH1D*)hmumu_trk_all[i]->Clone();
    hmumu_trk_fail[i]->Add(hmumu_trk_pass[i], -1);
    hmumu_id_fail[i] = (TH1D*)hmumu_id_all[i]->Clone();
    hmumu_id_fail[i]->Add(hmumu_id_pass[i], -1);

	double nall, npass;
	//Using rootfit
//	tf1mumu_trg_all[i]->SetParameter(0,hmumu_trg_all[i]->Integral());
//	hmumu_trg_all[i]->Fit(Form("tf1_trg_all%d",i), "L q", "", blow, bhigh);
//	hmumu_trg_all[i]->Fit(Form("tf1_trg_all%d",i), "L m q", "", blow, bhigh);
//	for(int k = 0; k < 5; k++){
//	  signal->SetParameter(k, tf1mumu_trg_all[i]->GetParameter(k));}
//	nall = signal->Integral(blow, bhigh)/measure;
//    tf1mumu_trg_pass[i]->SetParLimits(0,0,nall);
//	hmumu_trg_pass[i]->Fit(Form("tf1_trg_pass%d",i), "L q", "", blow, bhigh);
//	hmumu_trg_pass[i]->Fit(Form("tf1_trg_pass%d",i), "L m q", "", blow, bhigh);
//	for(int k = 0; k < 5; k++){
//	  signal->SetParameter(k, tf1mumu_trg_pass[i]->GetParameter(k));}
//	npass = signal->Integral(blow, bhigh)/measure;

	//Using sideband subtraction
	nall = hmumu_trg_all[i]->Integral(siglow,sighigh)-(hmumu_trg_all[i]->Integral(siglow-sigwidth, siglow-1)+hmumu_trg_all[i]->Integral(sighigh+1,sighigh+sigwidth));
	npass = hmumu_trg_pass[i]->Integral(siglow,sighigh)-(hmumu_trg_pass[i]->Integral(siglow-sigwidth, siglow-1)+hmumu_trg_pass[i]->Integral(sighigh+1,sighigh+sigwidth));
	eff_trg->SetBinContent(i+1, npass/nall);
	eff_trg->SetBinError(i+1, 0.00001);
//	  cout<<"his trg all : "<<hmumu_trg_all[i]->Integral()<<endl;
//    cout<<"his trg pass: "<<hmumu_trg_pass[i]->Integral()<<endl;
//    cout<<"trg all : "<<nall<<endl;
//    cout<<"trg pass: "<<npass<<endl;
//	  cout<<"eff trg: "<<npass/nall<<endl;

	//Using sideband subtraction
	nall = hmumu_trk_all[i]->Integral(siglow,sighigh)-(hmumu_trk_all[i]->Integral(siglow-sigwidth, siglow-1)+hmumu_trk_all[i]->Integral(sighigh+1,sighigh+sigwidth));
	npass = hmumu_trk_pass[i]->Integral(siglow,sighigh)-(hmumu_trk_pass[i]->Integral(siglow-sigwidth, siglow-1)+hmumu_trk_pass[i]->Integral(sighigh+1,sighigh+sigwidth));
	eff_trk->SetBinContent(i+1, npass/nall);
	eff_trk->SetBinError(i+1, 0.00001);
//	  cout<<"his trk all : "<<hmumu_trk_all[i]->Integral()<<endl;
//	  cout<<"his trk pass: "<<hmumu_trk_pass[i]->Integral()<<endl;
//    cout<<"trk all : "<<nall<<endl;
//    cout<<"trk pass: "<<npass<<endl;
//	  cout<<"eff trk: "<<npass/nall<<endl;

	//Using sideband subtraction
	nall = hmumu_id_all[i]->Integral(siglow,sighigh)-(hmumu_id_all[i]->Integral(siglow-sigwidth, siglow-1)+hmumu_id_all[i]->Integral(sighigh+1,sighigh+sigwidth));
	npass = hmumu_id_pass[i]->Integral(siglow,sighigh)-(hmumu_id_pass[i]->Integral(siglow-sigwidth, siglow-1)+hmumu_id_pass[i]->Integral(sighigh+1,sighigh+sigwidth));
	eff_id->SetBinContent(i+1, npass/nall);
	eff_id->SetBinError(i+1, 0.00001);
    cout<<"his id all : "<<hmumu_id_all[i]->Integral()<<endl;
    cout<<"his id pass: "<<hmumu_id_pass[i]->Integral()<<endl;
    cout<<"id all : "<<nall<<endl;
    cout<<"id pass: "<<npass<<endl;
    cout<<"eff id: "<<npass/nall<<endl;
  }
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
  TnPtreeTrg->Write();
  TnPtreeTrk->Write();
  TnPtreeID->Write();
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

  TCanvas *cforSave= new TCanvas("cforSave","",600,600);
  TLegend *legforSave = myLegend(0.77,0.50,0.90,0.90);
  cforSave->cd();
  for(int i = 0; i < nMuPtBin; i++){
    hmumu_trg_all[i]->Write();
    hmumu_trg_pass[i]->Write();
    hmumu_trg_fail[i]->Write();
    hmumu_trk_all[i]->Write();
    hmumu_trk_pass[i]->Write();
    hmumu_trk_fail[i]->Write();
    hmumu_id_all[i]->Write();
    hmumu_id_pass[i]->Write();
    hmumu_id_fail[i]->Write();

    hmumu_trg_all[i]->Draw();
    hmumu_trg_pass[i]->Draw("same");
    hmumu_trg_pass[i]->SetLineColor(3);
    hmumu_trg_fail[i]->Draw("same");
    hmumu_trg_fail[i]->SetLineColor(2);
	legforSave->Clear(); legforSave->SetTextSize(0.03); legforSave->SetHeader(Form("pt: %.1f ~ %.1f", MuPtBin[i], MuPtBin[i+1]));
    legforSave->AddEntry(hmumu_trg_all[i],  "all", "l");
    legforSave->AddEntry(hmumu_trg_pass[i], "pass", "l");
    legforSave->AddEntry(hmumu_trg_fail[i], "fail", "l");
    legforSave->Draw("same");
    cforSave->SaveAs(Form("%s/hmumu_trg_all_%d.pdf",plotfolder.Data(), i));

    hmumu_trk_all[i]->Draw();
    hmumu_trk_pass[i]->Draw("same");
    hmumu_trk_pass[i]->SetLineColor(3);
    hmumu_trk_fail[i]->Draw("same");
    hmumu_trk_fail[i]->SetLineColor(2);
	legforSave->Clear(); legforSave->SetTextSize(0.03); legforSave->SetHeader(Form("pt: %.1f ~ %.1f", MuPtBin[i], MuPtBin[i+1]));
    legforSave->AddEntry(hmumu_trk_all[i],  "all", "l");
    legforSave->AddEntry(hmumu_trk_pass[i], "pass", "l");
    legforSave->AddEntry(hmumu_trk_fail[i], "fail", "l");
    legforSave->Draw("same");
    cforSave->SaveAs(Form("%s/hmumu_trk_all_%d.pdf",plotfolder.Data(), i));

    hmumu_id_all[i]->Draw();
    hmumu_id_pass[i]->Draw("same");
    hmumu_id_pass[i]->SetLineColor(3);
    hmumu_id_fail[i]->Draw("same");
    hmumu_id_fail[i]->SetLineColor(2);
	legforSave->Clear(); legforSave->SetTextSize(0.03); legforSave->SetHeader(Form("pt: %.1f ~ %.1f", MuPtBin[i], MuPtBin[i+1]));
    legforSave->AddEntry(hmumu_id_all[i],  "all", "l");
    legforSave->AddEntry(hmumu_id_pass[i], "pass", "l");
    legforSave->AddEntry(hmumu_id_fail[i], "fail", "l");
    legforSave->Draw("same");
    cforSave->SaveAs(Form("%s/hmumu_id_all_%d.pdf",plotfolder.Data(), i));

    heta_trg_all[i]->Draw("p");
    heta_trg_pass[i]->SetLineColor(3);
    heta_trg_pass[i]->Draw("same");
    heta_trg_fail[i]->SetLineColor(2);
    heta_trg_fail[i]->Draw("same");
    cforSave->SaveAs(Form("%s/heta_trg_all_%d.pdf",plotfolder.Data(), i));

    heta_trk_all[i]->Draw();
    heta_trk_fail[i]->SetLineColor(3);
    heta_trk_pass[i]->Draw();
    heta_trk_fail[i]->SetLineColor(2);
    heta_trk_fail[i]->Draw();
    cforSave->SaveAs(Form("%s/heta_trk_all_%d.pdf",plotfolder.Data(), i));

    heta_id_all[i]->Draw();
    heta_id_pass[i]->SetLineColor(3);
    heta_id_pass[i]->Draw("same");
    heta_id_fail[i]->SetLineColor(2);
    heta_id_fail[i]->Draw("same");
	legforSave->Clear(); legforSave->SetTextSize(0.03); legforSave->SetHeader(Form("pt: %.1f ~ %.1f", MuPtBin[i], MuPtBin[i+1]));
    legforSave->AddEntry(heta_id_all[i],  "all", "l");
    legforSave->AddEntry(heta_id_pass[i], "pass", "l");
    legforSave->AddEntry(heta_id_fail[i], "fail", "l");
    legforSave->Draw("same");
    cforSave->SaveAs(Form("%s/heta_id_all_%d.pdf",plotfolder.Data(), i));
  }
  outf->Close();
  return 0;
}

