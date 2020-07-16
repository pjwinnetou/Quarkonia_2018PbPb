#include <ctime>

#include <TLorentzVector.h>
#include "commonUtility.h"
#include "HiEvtPlaneList.h"
#include "cutsAndBinChiB.h"
#include "TreeSetting.h"

static const long MAXTREESIZE = 1000000000000;

void SkimTree(int nevt=-1, bool isMC = true, int kTrigSel = kL1DoubleMuOpen, bool dirPb = false) 
{

  using namespace std;

  //Read Input File
  TString fnameData_PbpRun = "/home/samba.old/UpsilonAnalysis/DataFiles/Onia2018/Conversion/oniaTreeConversion_Pbp_GTV19_20200630_merged.root";
  TString fnameData_pPbRun = "/home/samba.old/UpsilonAnalysis/DataFiles/Onia2018/Conversion/oniaTreeConversion_pPb_GTV19_20200630_merged.root";
  TString fnameMC = "/home/deathold/work/CMS/analysis/chiBAna/chiBAna/oniaTree_ChiBMC_merged.root";

  TChain *mytree = new TChain("hionia/myTree");
  if(!isMC){
    if(dirPb) mytree->Add(fnameData_PbpRun.Data());
    else if(!dirPb) mytree->Add(fnameData_pPbRun.Data());
  }
  else if(isMC){
    mytree->Add(fnameMC.Data());
  }

  //SetBranchAddress
  SetTree settree_;
  settree_.TreeSetting(mytree,isMC,1);
//  SetTree::TreeSetting(mytree,isMC,1);
  

  //TriggerIndex
  int trigIndx=0;
  if(kTrigSel == kL1DoubleMuOpen) trigIndx=0;
  else if(kTrigSel == kTrigUps) trigIndx=1;
  else if(kTrigSel == kTrigL1DBOS40100) trigIndx=2;
  else if(kTrigSel == kTrigL1DB50100) trigIndx=3;
  
  //For New Writing File
  TFile* newfile;
  if(isMC) newfile = new TFile(Form("ChiBSkim_%sTrig_isMC%d.root",fTrigName[trigIndx].Data(),isMC),"recreate");
  else if(!isMC) newfile = new TFile(Form("ChiBSkim_%sTrig_isMC%d_RunPb%d.root",fTrigName[trigIndx].Data(),isMC,dirPb),"recreate");

  //For New Tree
  const static int nMaxDimu = 1000;
  int evt;
  int runN;
  int lumi;
  int cBin;
  int nCand;
  float vz;
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

  TTree* mmevttree = new TTree("mmepevt","dimuonAndEventPlanes in event based");
  mmevttree->SetMaxTreeSize(MAXTREESIZE);
  mmevttree->Branch("event",&evt,"event/I");
  mmevttree->Branch("runN",&runN,"runN/I");
  mmevttree->Branch("lumi",&lumi,"lumi/I");
  mmevttree->Branch("vz",&vz,"vz/F");
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


      


  ////////////////////////////////////////////////////////////////////////
  ////////////////// TLorentzVector dummies & conversion photon
  ////////////////////////////////////////////////////////////////////////
  TLorentzVector* JP_Reco = new TLorentzVector;
  TLorentzVector* mupl_Reco = new TLorentzVector;
  TLorentzVector* mumi_Reco = new TLorentzVector;
  
  TLorentzVector* conv_4mom = new TLorentzVector;
  TLorentzVector* conv_4mom_ = new TLorentzVector;
  TLorentzVector* chib_4mom = new TLorentzVector;


  // event loop start
  if(nevt == -1) nevt = mytree->GetEntries();
  for(int iev=0; iev<nevt ; ++iev)
  {
    if(iev%100000==0) cout << ">>>>> EVENT " << iev << " / " << mytree->GetEntries() <<  " ("<<(int)(100.*iev/mytree->GetEntries()) << "%)" << endl;

    mytree->GetEntry(iev);

    //Trigger Matching
//    if(!isMC){
      if(!( (HLTriggers&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;
//    }

    //Dimuon Loop
    nCand = 0;
    for (Int_t irqq=0; irqq<Reco_QQ_size; ++irqq) 
    {
      //Dimuon Filter Matching
//      if(!isMC){
        if(!( (Reco_QQ_trig[irqq]&((ULong64_t)pow(2, kTrigSel))) == ((ULong64_t)pow(2, kTrigSel)) ) ) continue;
//      }

      runN = runNb;
      evt = eventNb;
      lumi = LS;
      vz = zVtx;

      JP_Reco = (TLorentzVector*) Reco_QQ_4mom->At(irqq);
      mupl_Reco = (TLorentzVector*) Reco_QQ_mupl_4mom->At(irqq);
      mumi_Reco = (TLorentzVector*) Reco_QQ_mumi_4mom->At(irqq);
      
   
      //Soft Muon Id cut
      if(!settree_.SoftMuIdCut(irqq)) continue;
     
      //Dimuon Vertex probability cut 
      if ( Reco_QQ_VtxProb[irqq]  < 0.01 ) continue;

      //Acceptance cut --> Based on Jpsi trigger but just for skimming at this stage
      if(! (IsAcceptanceQQ(mupl_Reco->Pt(), mupl_Reco->Eta()) && IsAcceptanceQQ(mumi_Reco->Pt(), mumi_Reco->Eta()) ) ) continue;


      //Loop for Converted photon 

      for(int ipho =0; ipho<Reco_conv_size; ipho++)
      {
        conv_4mom = (TLorentzVector*) Reco_conv_4mom->At(ipho);
        conv_4mom_ -> SetPtEtaPhiM( conv_4mom->Pt(), conv_4mom->Eta(), conv_4mom->Phi(), 0 );
        
        //oposite charge
        if(Reco_conv_tk1_charge[ipho]+Reco_conv_tk2_charge[ipho]!=0) continue;

        //4 and 3 hits in tracker
        bool isTrackerHitConv = false;
        if(Reco_conv_tk1_nhit[ipho]>4 && Reco_conv_tk2_nhit[ipho]>3) isTrackerHitConv=true;
        if(Reco_conv_tk2_nhit[ipho]>4 && Reco_conv_tk1_nhit[ipho]>3) isTrackerHitConv=true;
//        if(isTrackerHitConv==false) continue;

        //difference in z coordinate less than 5 cm 

        //reduced chi2 less than 10
//        if(Reco_conv_chi2[ipho]>=10) continue;

        //chi2 probability cut
//        if(Reco_conv_chi2_probability[ipho]<=5*10e-4) continue;

        //delta phi cot theta cut
        if(Reco_conv_dphitrksatvtx[ipho]>=0.2 || Reco_conv_paircotthetasep[ipho]>=0.1) continue;

        //distance of minimum approach cut 
//        if(! (Reco_conv_distofminapproach[ipho]<1 && Reco_conv_distofminapproach[ipho]>-0.25) ) continue;

        //impact parameter d0 times charge 
//        if(Reco_conv_tk2_d0[ipho]*Reco_conv_tk2_charge[ipho]<0) continue;
        
        //Dalitz decay pi0 suppression
//        if(abs(Reco_conv_zofprimvtxfromtrks[ipho])<1.5) continue;
        
        //photon eta cut
//        if(abs(conv_4mom->Eta())>1.4) continue;
      
        //make chi_b candidate and Q-value cut
        *chib_4mom = *JP_Reco+*conv_4mom_;
        double Qval = chib_4mom->M()-JP_Reco->M();

//        if(abs(Qval)>=2) continue;
//        cout << "nothing1??" << endl;

        //dz of photon w.r.t dimuon less than 1mm 
        float dzofdimuon = sqrt(Reco_QQ_ctau3D[irqq]*Reco_QQ_ctau3D[irqq] - Reco_QQ_ctau[irqq]*Reco_QQ_ctau[irqq]);
//        if(abs( (Reco_conv_zofprimvtxfromtrks[ipho]-dzofdimuon) ) >=0.1) continue;
//        cout << "nothing2??" << endl;
      

        recoQQsign[nCand] = Reco_QQ_sign[irqq];     
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

        chi_mass[nCand] = chib_4mom->M();
        chi_y[nCand] = chib_4mom->Rapidity();
        chi_pt[nCand] = chib_4mom->Pt();
        
        chi_dimu_mass_diff[nCand] = Qval;
        nCand++;
        
      } // end of photon loop


    } // end of dimuon loop

    if(nCand>0) mmevttree->Fill();
    
  } //end of event loop

  newfile->cd();
  mmevttree->Write();
  newfile->Close();
  
} 
