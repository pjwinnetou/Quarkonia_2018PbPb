#include <ctime>

#include <TLorentzVector.h>
#include "commonUtility.h"
#include "HiEvtPlaneList.h"
#include "cutsAndBinChiB.h"
#include "TreeSettingOtaVer.h"

static const long MAXTREESIZE = 1000000000000;

bool isAcceptance(double eta, double pt)
{
  if (fabs(eta) > 2.4) return false;
  if (fabs(eta) < 0.3 && pt < 3.4) return false;
  if (fabs(eta) >= 0.3 && fabs(eta) < 1.1 && pt < 3.3) return false;
  if (fabs(eta) >= 1.1 && fabs(eta) < 2.1 && (pt < 5.5 - 2.0*fabs(eta))) return false;
  if (fabs(eta) >= 2.1 && pt<1.3) return false;
  return true;
}


void SkimTree_OtaVer(int nevt=-1, bool isMC = false, bool SignQQ = true, bool dirPb = true) 
{

  using namespace std;

  //Read Input File
  TString fnameData_PbpRun = "/home/samba.old/CMS_Files/ChiAna/skimmedFiles/OtaAna/chi_Pbp_tree_merged_v2.root";
  TString fnameData_pPbRun = "/home/samba.old/CMS_Files/ChiAna/skimmedFiles/OtaAna/chi_pPb_tree_merged_v2.root";
  TString fnameMC = "/home/deathold/work/CMS/analysis/chiBAna/chiBAna/oniaTree_Gen_v4_merged.root";

  TChain *mytree = new TChain("ChiRootuple/event_tree");
  if(!isMC){
    if(dirPb) mytree->Add(fnameData_PbpRun.Data());
    else if(!dirPb) mytree->Add(fnameData_pPbRun.Data());
  }
  else if(isMC){
    mytree->Add(fnameMC.Data());
  }

  //SetBranchAddress
  SetTree settree_;
  settree_.TreeSetting(mytree,isMC);
//  SetTree::TreeSetting(mytree,isMC,1);
  

  
  //For New Writing File
  TFile* newfile;
  if(isMC) newfile = new TFile(Form("skimmedFiles/ChiSkim_isMC%d_sameSign%d.root",isMC,SignQQ),"recreate");
  else if(!isMC) newfile = new TFile(Form("skimmedFiles/ChiSkim_isMC%d_sameSign%d_RunPb%d.root",isMC,SignQQ,dirPb),"recreate");

  //For New Tree
  const static int nMaxDimu = 1000;
  int evt;
  int runN;
  int lumi;
  int cBin;
  int nCand;
  double vz;
  float Dimu_mass[nMaxDimu];
  float Dimu_pt[nMaxDimu];
  float Dimu_y[nMaxDimu];
  float Dimu_phi[nMaxDimu];
  float Dimu_eta[nMaxDimu];
  float Simu_eta1[nMaxDimu];
  float Simu_eta2[nMaxDimu];
  float Simu_phi1[nMaxDimu];
  float Simu_phi2[nMaxDimu];
  float Simu_pt1[nMaxDimu];
  float Simu_pt2[nMaxDimu];
  int recoQQsign[nMaxDimu];
  
  float pho_mass[nMaxDimu];
  float pho_eta[nMaxDimu];
  float pho_y[nMaxDimu];
  float pho_pt[nMaxDimu];
  float pho1_pt[nMaxDimu];
  float pho2_pt[nMaxDimu];
  float pho1_eta[nMaxDimu];
  float pho2_eta[nMaxDimu];

  float chi_mass[nMaxDimu];
  float chi_pt[nMaxDimu];
  float chi_y[nMaxDimu];

  float chi_dimu_mass_diff[nMaxDimu];
  float chi_b_pdg[nMaxDimu];
  float chi_c_pdg[nMaxDimu];


  bool isTrackerHitConv[nMaxDimu];
  bool isReco_conv_chi2[nMaxDimu];
  bool isReco_conv_chi2_probability[nMaxDimu];
  bool isReco_conv_dphitrksatvtx[nMaxDimu];
  bool isReco_conv_distofminapproach[nMaxDimu];
  bool isReco_conv_tk2_d0[nMaxDimu];
  bool isReco_conv_zofprimvtxfromtrks[nMaxDimu];
  bool isdzofdimuon[nMaxDimu];
  int nTracks;


  TTree* mmevttree = new TTree("mmepevt","dimuonAndEventPlanes in event based");
  mmevttree->SetMaxTreeSize(MAXTREESIZE);
  mmevttree->Branch("event",&evt,"event/I");
//  mmevttree->Branch("runN",&runN,"runN/I");
//  mmevttree->Branch("lumi",&lumi,"lumi/I");
  mmevttree->Branch("nCand",&nCand,"nCand/I");
  mmevttree->Branch("Dimu_mass",Dimu_mass,"Dimu_mass[nCand]/F");
  mmevttree->Branch("Dimu_y",Dimu_y,"Dimu_y[nCand]/F");
  mmevttree->Branch("Dimu_pt",Dimu_pt,"Dimu_pt[nCand]/F");
  mmevttree->Branch("Dimu_eta",Dimu_eta,"Dimu_eta[nCand]/F");
  mmevttree->Branch("Simu_pt1",Simu_pt1,"Simu_pt1[nCand]/F");
  mmevttree->Branch("Simu_pt2",Simu_pt2,"Simu_pt2[nCand]/F");
  mmevttree->Branch("Simu_eta1",Simu_eta1,"Simu_eta1[nCand]/F");
  mmevttree->Branch("Simu_eta2",Simu_eta2,"Simu_eta2[nCand]/F");
  mmevttree->Branch("recoQQsign",recoQQsign,"recoQQsign[nCand]/I");

  mmevttree->Branch("pho_mass",pho_mass,"pho_mass[nCand]/F");
  mmevttree->Branch("pho_eta",pho_eta,"pho_eta[nCand]/F");
  mmevttree->Branch("pho_y",pho_y,"pho_y[nCand]/F");
  mmevttree->Branch("pho_pt",pho_pt,"pho_pt[nCand]/F");
//  mmevttree->Branch("pho1_pt",pho1_pt,"pho1_pt[nCand]/F");
//  mmevttree->Branch("pho2_pt",pho2_pt,"pho2_pt[nCand]/F");
//  mmevttree->Branch("pho1_eta",pho1_eta,"pho1_eta[nCand]/F");
//  mmevttree->Branch("pho2_eta",pho2_eta,"pho2_eta[nCand]/F");

  mmevttree->Branch("chi_mass",chi_mass,"chi_mass[nCand]/F");
  mmevttree->Branch("chi_pt",chi_pt,"chi_pt[nCand]/F");
  mmevttree->Branch("chi_y",chi_y,"chi_y[nCand]/F");
  
  mmevttree->Branch("chi_dimu_mass_diff",chi_dimu_mass_diff,"chi_dimu_mass_diff[nCand]/F");
  mmevttree->Branch("chi_b_pdg",chi_b_pdg,"chi_b_pdg[nCand]/F");
  mmevttree->Branch("chi_c_pdg",chi_c_pdg,"chi_c_pdg[nCand]/F");
  mmevttree->Branch("nTracks",&nTracks,"nTracks/I");
 

  float gen_dimu_mass[nMaxDimu];
  float gen_dimu_pt[nMaxDimu];
  float gen_dimu_y[nMaxDimu];
  float gen_chi_mass[nMaxDimu];
  float gen_chi_pt[nMaxDimu];
  float gen_chi_y[nMaxDimu];
  if(isMC){

    mmevttree->Branch("Gen_dimu_mass",gen_dimu_mass,"Gen_dimu_mass[nCand]/F");
    mmevttree->Branch("Gen_dimu_pt",gen_dimu_pt,"Gen_dimu_pt[nCand]/F");
    mmevttree->Branch("Gen_dimu_y",gen_dimu_y,"Gen_dimu_y[nCand]/F");
    
    mmevttree->Branch("Gen_chi_mass",gen_chi_mass,"Gen_chi_mass[nCand]/F");
    mmevttree->Branch("Gen_chi_pt",gen_chi_pt,"Gen_chi_pt[nCand]/F");
    mmevttree->Branch("Gen_chi_y",gen_chi_y,"Gen_chi_y[nCand]/F");
  }
      


  ////////////////////////////////////////////////////////////////////////
  ////////////////// TLorentzVector dummies & conversion photon
  ////////////////////////////////////////////////////////////////////////
  TLorentzVector* JP_Reco     = new TLorentzVector;
  TLorentzVector* mupl_Reco   = new TLorentzVector;
  TLorentzVector* mumi_Reco   = new TLorentzVector;
  
  TLorentzVector* conv_4mom   = new TLorentzVector;
  TLorentzVector* conv_4mom_  = new TLorentzVector;
  TLorentzVector* chi_4mom   = new TLorentzVector;



  // event loop start
  if(nevt == -1) nevt = mytree->GetEntries();
  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);

    evt = (Int_t) eventNumber; 
    nTracks = ntracks_inEvent;

    nCand=0;
    for(int ir=0; ir < chiCandPerEvent; ir++){
      int irqq = chi_daughterJpsi_position->at(ir);
      if(dimuon_vtxProb->at(irqq) <= 0.01) continue;

      bool isOS = (dimuon_charge->at(irqq)==0) ? true : false;
      if(SignQQ){if(!isOS) continue;}
      else if(!SignQQ){if(isOS) continue;}

      int irsq1 = dimuon_muon1_position->at(irqq);
      int irsq2 = dimuon_muon1_position->at(irqq);
      if(!settree_.SoftMuIdCut(irsq1)) continue;
      if(!settree_.SoftMuIdCut(irsq2)) continue;
      if(!settree_.TriggerL1DBMuOpenCut(irsq1)) continue;
      if(!settree_.TriggerL1DBMuOpenCut(irsq2)) continue;

      int ipho = chi_daughterConv_position->at(ir);
      if(!settree_.PhoHighPurityIdCut(ipho)) continue;
      if(abs(conv_eta->at(ipho))>2.4 || conv_pt->at(ipho)<0.2) continue;

      JP_Reco = (TLorentzVector*) dimuon_p4->At(irqq);
      mupl_Reco = (TLorentzVector*) muon_p4->At(irsq1);
      mumi_Reco = (TLorentzVector*) muon_p4->At(irsq2);

      if(!isAcceptance(mupl_Reco->Eta(), mupl_Reco->Pt())) continue;
      if(!isAcceptance(mumi_Reco->Eta(), mumi_Reco->Pt())) continue;

      conv_4mom = (TLorentzVector*) conv_p4->At(ipho);
      chi_4mom  = (TLorentzVector*) chi_p4->At(ir);

      double Qval = chi_4mom->M()-JP_Reco->M();
      recoQQsign[nCand] = dimuon_charge->at(irqq); 
      Dimu_mass[nCand] = JP_Reco->M();
      Dimu_phi[nCand] = JP_Reco->Phi();
      Simu_phi1[nCand] = mupl_Reco->Phi();
      Simu_phi2[nCand] = mumi_Reco->Phi();
      Dimu_eta[nCand] = JP_Reco->Eta();
      Dimu_y[nCand] = JP_Reco->Rapidity();
      Dimu_pt[nCand] = JP_Reco->Pt();
      Simu_pt1[nCand] = mupl_Reco->Pt();
      Simu_pt2[nCand] = mumi_Reco->Pt();
      Simu_eta1[nCand] = mupl_Reco->Eta();
      Simu_eta2[nCand] = mumi_Reco->Eta();

      pho_mass[nCand] = conv_4mom->M(); 
      pho_pt[nCand] = conv_4mom->Pt(); 
      pho_y[nCand] = conv_4mom->Rapidity(); 
      pho_eta[nCand] = conv_4mom->Eta(); 

      chi_mass[nCand] = chi_4mom->M();
      chi_y[nCand] = chi_4mom->Rapidity();
      chi_pt[nCand] = chi_4mom->Pt();

      chi_dimu_mass_diff[nCand] = Qval;
      chi_b_pdg[nCand] = pdgMass.Y1S + Qval;
      chi_c_pdg[nCand] = pdgMass.JPsi + Qval;

      if(dirPb){
        Dimu_eta[nCand] = -JP_Reco->Eta();
        Dimu_y[nCand] = -JP_Reco->Rapidity();
        Simu_eta1[nCand] = -mupl_Reco->Eta();
        Simu_eta2[nCand] = -mumi_Reco->Eta();
        pho_eta[nCand] = -conv_4mom->Eta();
        pho_y[nCand] = -conv_4mom->Rapidity();
        chi_y[nCand] = -chi_4mom->Rapidity();
      }
      nCand++;
    } 

    if(nCand>0) mmevttree->Fill();

  } //end of event loop

  newfile->cd();
  mmevttree->Write();
  newfile->Close();
  
} 
