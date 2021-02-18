#ifndef TreeSettingOtaVer_h
#define TreeSettingOtaVer_h

#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TROOT.h"
#include "TClonesArray.h"
#include "TChain.h"

using namespace std;

const int maxBranchSize = 1000;

// import the tree to the RooDataSet
const int PythCode_chic0 = 10441; //Pythia codes
const int PythCode_chic1 = 20443;
const int PythCode_chic2 = 445;


int runNumber, eventNumber;
TLorentzVector* LVchic, *LVdimuon, *LVconv, *LVmuon1, *LVmuon2;
TLorentzVector* LVchic_rot, *LVchic_rotGamma, *LVdimuon_rot, *LVconv_rot, *LVmuon1_rot, *LVmuon2_rot;
TLorentzVector LVaux;


int ntracks_inEvent;
double hfTowerSum_inEvent;
int Trig_Event_HLTDoubleMuOpen;
int nPrimVertices;
int chiCandPerEvent;


//
vector <double>* pvtx_z = 0;
vector <double>* pvtx_zError = 0;
vector <double>* pvtx_x = 0;
vector <double>* pvtx_y = 0;
vector <double>* pvtx_nTracks = 0;
vector <bool>* pvtx_isFake = 0;




//muon info
TClonesArray*  muon_p4;
vector <bool>* muonIsHLTDoubleMuOpen = 0;
vector <bool>* muonIsGlobal = 0;
vector <bool>* muonIsTracker = 0;
vector <bool>* muonIsPF = 0;
vector <bool>* muonIsSoft = 0;
vector <bool>* muonIsTight = 0;
vector <int>* muonTrackerLayersWithMeasurement = 0;
vector <double>* muon_eta = 0;
vector <double>* muon_pt = 0;


//dimuon
TClonesArray*  dimuon_p4;
vector <double>* dimuon_eta = 0;
vector <double>* dimuon_pt = 0;
vector <double>* dimuon_charge = 0;
vector <int>*  dimuon_pvtx_index = 0;
vector <double>* dimuon_dz_dimuonvtx_pvtx = 0;
vector <double>* dimuon_vtxProb = 0;
vector <int>* dimuon_muon1_position = 0; //stores position of first muon in muon collection (no specific order)
vector <int>* dimuon_muon2_position = 0; //stores position of second muon in muon collection (no specific order)
vector <double>* dimuon_ctpv = 0;
vector <double>* dimuon_ctpvError = 0;





// Chi
TClonesArray*  chi_p4 ;

vector <double>* chi_eta = 0;
vector <double>* chi_pt = 0;
vector <int>* chi_daughterJpsi_position = 0; //stores position of daughter Jpsi in dimuon collection
vector <int>* chi_daughterConv_position = 0; //stores position of daughter photon (conversion)
vector <double>* chi_dzPhotToDimuonVtx = 0; //z distance of photon to dimuon vertex when dxy is minimal
vector <double>* chi_dxyPhotToDimuonVtx = 0; //dxy distance of photon to dimuon vertex when dz is 0 - probably not too good for very midrapidity conversions




//// Conversions  /////


TClonesArray*  conv_p4; 
vector <bool>* convQuality_isHighPurity = 0;
vector <bool>* convQuality_isGeneralTracksOnly = 0;
vector <double>* conv_vertexPositionRho = 0;
vector <double>* conv_sigmaTkVtx1 = 0;
vector <double>* conv_sigmaTkVtx2 = 0;
vector <bool>* conv_tkVtxCompatibilityOK = 0;


vector <int>* conv_compatibleInnerHitsOK = 0; //-1: less than 2 tracks, 0: not compatible, 1: yes
vector <double>* conv_vertexChi2Prob = 0;
vector <double>* conv_zOfPriVtx = 0;
vector <double>* conv_zOfPriVtxFromTracks = 0;
vector <double>* conv_dzToClosestPriVtx = 0;
vector <double>* conv_dxyPriVtx_Tr1 = 0;
vector <double>* conv_dxyPriVtx_Tr2 = 0;
vector <double>* conv_dxyPriVtxTimesCharge_Tr1 = 0;
vector <double>* conv_dxyPriVtxTimesCharge_Tr2 = 0;
vector <double>* conv_dxyError_Tr1 = 0;
vector <double>* conv_dxyError_Tr2 = 0;


vector <int>* conv_tk1NumOfDOF = 0;
vector <int>* conv_tk2NumOfDOF = 0;
vector <double>* conv_track1Chi2 = 0;
vector <double>* conv_track2Chi2 = 0;
vector <double>* conv_minDistanceOfApproach = 0;
vector <double>* conv_eta = 0;
vector <double>* conv_pt = 0;


//MC
vector <int>* gen_Jpsi_matchPosition = 0;
vector <int>* gen_conv_matchPosition = 0;
vector <int>* convGen_motherCode;

class SetTree
{
  public:
    SetTree(){};

    virtual ~SetTree();
    virtual void TreeSetting(TTree* event_tree, bool isMC);
    Bool_t SoftMuIdCut(int irsq);
    Bool_t TriggerL1DBMuOpenCut(int irsq);
    Bool_t PhoHighPurityIdCut(int ipho);
};

SetTree::~SetTree()
{
}


void SetTree::TreeSetting(TTree* event_tree, bool isMC)
{
  chi_p4=0;
  dimuon_p4=0;
  conv_p4=0;


  event_tree->SetBranchAddress("eventNumber", &eventNumber);
  event_tree->SetBranchAddress("nPrimVertices", &nPrimVertices);
  event_tree->SetBranchAddress("ntracks_inEvent", &ntracks_inEvent);
  event_tree->SetBranchAddress("chiCandPerEvent", &chiCandPerEvent);


  event_tree->SetBranchAddress("pvtx_z", &pvtx_z);
  event_tree->SetBranchAddress("pvtx_zError", &pvtx_zError);
  event_tree->SetBranchAddress("pvtx_x", &pvtx_x);
  event_tree->SetBranchAddress("pvtx_y", &pvtx_y);
  event_tree->SetBranchAddress("pvtx_nTracks", &pvtx_nTracks);
  event_tree->SetBranchAddress("pvtx_isFake", &pvtx_isFake);


  event_tree->SetBranchAddress("muon_p4", &muon_p4);
  event_tree->SetBranchAddress("muonIsHLTDoubleMuOpen", &muonIsHLTDoubleMuOpen);
  event_tree->SetBranchAddress("muonIsGlobal", &muonIsGlobal);
  event_tree->SetBranchAddress("muonIsTracker", &muonIsTracker);
  event_tree->SetBranchAddress("muonIsPF", &muonIsPF);
  event_tree->SetBranchAddress("muonIsSoft", &muonIsSoft);
  event_tree->SetBranchAddress("muonIsTight", &muonIsTight);
  event_tree->SetBranchAddress("muonTrackerLayersWithMeasurement", &muonTrackerLayersWithMeasurement);
  event_tree->SetBranchAddress("muon_eta", &muon_eta);
  event_tree->SetBranchAddress("muon_pt", &muon_pt);


  event_tree->SetBranchAddress("dimuon_p4", &dimuon_p4);
  event_tree->SetBranchAddress("dimuon_eta", &dimuon_eta);
  event_tree->SetBranchAddress("dimuon_pt", &dimuon_pt);
  event_tree->SetBranchAddress("dimuon_charge", &dimuon_charge);
  event_tree->SetBranchAddress("dimuon_pvtx_index", &dimuon_pvtx_index);
  event_tree->SetBranchAddress("dimuon_dz_dimuonvtx_pvtx", &dimuon_dz_dimuonvtx_pvtx);
  event_tree->SetBranchAddress("dimuon_vtxProb", &dimuon_vtxProb);
  event_tree->SetBranchAddress("dimuon_muon1_position", &dimuon_muon1_position);
  event_tree->SetBranchAddress("dimuon_muon2_position", &dimuon_muon2_position);
  event_tree->SetBranchAddress("dimuon_ctpv", &dimuon_ctpv);
  event_tree->SetBranchAddress("dimuon_ctpvError", &dimuon_ctpvError);


  event_tree->SetBranchAddress("chi_p4", &chi_p4);
  event_tree->SetBranchAddress("chi_eta", &chi_eta);
  event_tree->SetBranchAddress("chi_pt", &chi_pt);
  event_tree->SetBranchAddress("chi_daughterJpsi_position", &chi_daughterJpsi_position);
  event_tree->SetBranchAddress("chi_daughterConv_position", &chi_daughterConv_position);
  event_tree->SetBranchAddress("chi_dzPhotToDimuonVtx", &chi_dzPhotToDimuonVtx);
  event_tree->SetBranchAddress("chi_dxyPhotToDimuonVtx", &chi_dxyPhotToDimuonVtx);


  //Conversions
  event_tree->SetBranchAddress("conv_p4", &conv_p4);
  event_tree->SetBranchAddress("convQuality_isHighPurity", &convQuality_isHighPurity);
  event_tree->SetBranchAddress("convQuality_isGeneralTracksOnly", &convQuality_isGeneralTracksOnly);
  event_tree->SetBranchAddress("conv_vertexPositionRho", &conv_vertexPositionRho);
  event_tree->SetBranchAddress("conv_sigmaTkVtx1", &conv_sigmaTkVtx1);
  event_tree->SetBranchAddress("conv_sigmaTkVtx2", &conv_sigmaTkVtx2);
  event_tree->SetBranchAddress("conv_tkVtxCompatibilityOK", &conv_tkVtxCompatibilityOK);


  event_tree->SetBranchAddress("conv_compatibleInnerHitsOK", &conv_compatibleInnerHitsOK);
  event_tree->SetBranchAddress("conv_vertexChi2Prob", &conv_vertexChi2Prob);
  event_tree->SetBranchAddress("conv_zOfPriVtx", &conv_zOfPriVtx);
  event_tree->SetBranchAddress("conv_zOfPriVtxFromTracks", &conv_zOfPriVtxFromTracks);
  event_tree->SetBranchAddress("conv_dzToClosestPriVtx", &conv_dzToClosestPriVtx);
  event_tree->SetBranchAddress("conv_dxyPriVtx_Tr1", &conv_dxyPriVtx_Tr1);
  event_tree->SetBranchAddress("conv_dxyPriVtx_Tr2", &conv_dxyPriVtx_Tr2);
  event_tree->SetBranchAddress("conv_dxyPriVtxTimesCharge_Tr1", &conv_dxyPriVtxTimesCharge_Tr1);
  event_tree->SetBranchAddress("conv_dxyPriVtxTimesCharge_Tr2", &conv_dxyPriVtxTimesCharge_Tr2);
  event_tree->SetBranchAddress("conv_dxyError_Tr1", &conv_dxyError_Tr1);
  event_tree->SetBranchAddress("conv_dxyError_Tr2", &conv_dxyError_Tr2);


  event_tree->SetBranchAddress("conv_tk1NumOfDOF", &conv_tk1NumOfDOF);
  event_tree->SetBranchAddress("conv_tk2NumOfDOF", &conv_tk2NumOfDOF);
  event_tree->SetBranchAddress("conv_track1Chi2", &conv_track1Chi2);
  event_tree->SetBranchAddress("conv_track2Chi2", &conv_track2Chi2);
  event_tree->SetBranchAddress("conv_minDistanceOfApproach", &conv_minDistanceOfApproach);
  event_tree->SetBranchAddress("conv_eta", &conv_eta);
  event_tree->SetBranchAddress("conv_pt", &conv_pt);


};

Bool_t SetTree::SoftMuIdCut(int irsq)
{
  bool isSoftMuon = false;
  if(muonIsSoft->at(irsq)) isSoftMuon=true;
  return isSoftMuon;
}  

Bool_t SetTree::TriggerL1DBMuOpenCut(int irsq)
{
  bool isTriggerMuon = false;
  if(muonIsHLTDoubleMuOpen->at(irsq)) isTriggerMuon=true;
  return isTriggerMuon;
}  

Bool_t SetTree::PhoHighPurityIdCut(int ipho)  //uses variables loaded in main function
{
  bool isHighPurityPhoton = true;
  if (convQuality_isHighPurity->at(ipho) != 1) isHighPurityPhoton = false;
  if (convQuality_isGeneralTracksOnly->at(ipho) != 1) isHighPurityPhoton = false;
  if (conv_vertexPositionRho->at(ipho) <= 1.5) isHighPurityPhoton = false;
  if (conv_sigmaTkVtx1->at(ipho) > 10) isHighPurityPhoton = false;
  if (conv_sigmaTkVtx2->at(ipho) > 10) isHighPurityPhoton = false;
  if (conv_tkVtxCompatibilityOK->at(ipho) != 1) isHighPurityPhoton=false;
  if (conv_compatibleInnerHitsOK->at(ipho) != 1) isHighPurityPhoton=false;
  if (conv_vertexChi2Prob->at(ipho) <= 0.01) isHighPurityPhoton=false;
  if (fabs(conv_zOfPriVtx->at(ipho)) >= 20) isHighPurityPhoton=false;
  if (fabs(conv_dzToClosestPriVtx->at(ipho)) >= 10) isHighPurityPhoton=false;
  if (conv_tk1NumOfDOF->at(ipho) < 2.5) isHighPurityPhoton=false;
  if (conv_tk2NumOfDOF->at(ipho) < 2.5) isHighPurityPhoton=false;
  if (conv_track1Chi2->at(ipho) >= 10) isHighPurityPhoton=false;
  if (conv_track2Chi2->at(ipho) >= 10) isHighPurityPhoton=false;
  if (conv_minDistanceOfApproach->at(ipho) <= -0.25) isHighPurityPhoton=false;
  if (conv_minDistanceOfApproach->at(ipho) >= 1.00) isHighPurityPhoton=false;

  return isHighPurityPhoton;
}
#endif 
