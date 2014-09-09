#include "TFile.h"                                                                                                                                                                                          
#include "TChain.h"
#include "TTree.h"
#include "format.h"
#include "iostream"
using namespace std;
void reduce(){
//	string inf = "/net/hisrv0001/home/tawei/tawei/Bfinder/BfinderBoostedMC_20140701_inclBtoPsiMuMu_pa_STARTHI53_V27-v1_00000/*.root";
//	string inf = "/net/hisrv0001/home/tawei/tawei/Bfinder/BfinderBoostedMC_20140701_inclBtoPsiMuMu_pa_STARTHI53_V27-v1_20000/*.root";
//	string inf = "/net/hisrv0001/home/tawei/tawei/Bfinder/Calo_BfinderBoostedMC_20140701_inclBtoPsiMuMu_pa_STARTHI53_V27-v1_00000/*.root";
//	string inf = "/mnt/hadoop/cms/store/user/twang/HI_Btuple/20140707_PAMuon_HIRun2013_28Sep2013_v1/*.root";
//	string inf = "/net/hisrv0001/home/tawei/tawei/Bfinder/BfinderBoostedMC_20140708_inclBtoPsiMuMu_pa_STARTHI53_V27-v1_00000/*.root";
	string inf = "/mnt/hadoop/cms/store/user/twang/HI_Btuple/hckim-HIJINGemb_inclBtoPsiMuMu/*.root";
    TChain *nt = new TChain("demo/root");
    TChain *hlt = new TChain("hltanalysis/HltTree");
	nt->Add(inf.c_str());
	hlt->Add(inf.c_str());
    nt->SetBranchStatus("*" ,0);
    nt->SetBranchStatus("MuonInfo*" ,1);
    hlt->SetBranchStatus("*" ,0);
    hlt->SetBranchStatus("HLT_PAMu3_v1" ,1);
    //nt->AddFriend(hlt);
	MuonInfoBranches MuonInfo;
    MuonInfo.setbranchadd(nt);
    Int_t HLT_PAMu3_v1;
    hlt->SetBranchAddress("HLT_PAMu3_v1",&HLT_PAMu3_v1);
    TFile *file = new TFile ("test.root","recreate");
    file->mkdir("demo");
    file->cd("demo");
    TTree *tree = nt->CloneTree(0);
    TTree *hlttree = hlt->CloneTree(0);

    int nevents_total = nt->GetEntries();
    int nevents = 0;
    for(int entry=0;entry<nevents_total;entry++) {
    if ((entry%10000) == 0) printf("Loading event #%d of %d.\n",entry,nevents_total);
      nt->GetEntry(entry);
      hlt->GetEntry(entry);
      if(MuonInfo.size<2) continue;
	  if(!HLT_PAMu3_v1) continue;
	  tree->Fill(); // book this event to the new tree
      hlttree->Fill();
      nevents++;
    }
    printf("Finished - %d events saved from %d events.\n",nevents,nevents_total);
    tree->Write();
    hlttree->Write();
}
