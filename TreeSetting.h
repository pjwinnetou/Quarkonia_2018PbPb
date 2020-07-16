#ifndef TreeSetting_h
#define TreeSetting_h

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
UInt_t          runNb;
UInt_t          eventNb, LS;
float           zVtx;
Int_t           Centrality;
ULong64_t       HLTriggers;
Int_t           Reco_QQ_size;
TClonesArray    *Reco_QQ_4mom;
TClonesArray    *Reco_QQ_mupl_4mom;
TClonesArray    *Reco_QQ_mumi_4mom;
ULong64_t       Reco_QQ_trig[maxBranchSize];   //[Reco_QQ_size]
Float_t         Reco_QQ_VtxProb[maxBranchSize];   //[Reco_QQ_size]
TBranch        *b_runNb;   //!
TBranch        *b_eventNb;   //!
TBranch        *b_LS;
TBranch        *b_zVtx;   //!
TBranch        *b_Centrality;   //!
TBranch        *b_HLTriggers;   //!
TBranch        *b_Reco_QQ_size;   //!
TBranch        *b_Reco_QQ_4mom;   //!
TBranch        *b_Reco_QQ_mupl_4mom;   //!
TBranch        *b_Reco_QQ_mumi_4mom;   //!
TBranch        *b_Reco_QQ_trig;   //!
TBranch        *b_Reco_QQ_VtxProb;   //!

Bool_t          Reco_QQ_mupl_highPurity[maxBranchSize];   //[Reco_QQ_size]
Bool_t          Reco_QQ_mumi_highPurity[maxBranchSize];   //[Reco_QQ_size]
TBranch        *b_Reco_QQ_mupl_highPurity;   //!
TBranch        *b_Reco_QQ_mumi_highPurity;   //!

//  mytree->GetBranch("Reco_QQ_mupl_4mom")->SetAutoDelete(kFALSE);
//  mytree->GetBranch("Reco_QQ_mumi_4mom")->SetAutoDelete(kFALSE);

//  muon id 
Int_t           Reco_QQ_mupl_nTrkHits[maxBranchSize];   //[Reco_QQ_size]
Int_t           Reco_QQ_mumi_nTrkHits[maxBranchSize];   //[Reco_QQ_size]
Int_t           Reco_mu_nTrkHits[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_QQ_mupl_nTrkHits;   //!
TBranch        *b_Reco_QQ_mumi_nTrkHits;   //!
TBranch        *b_Reco_mu_nTrkHits;   //!
Float_t         Reco_QQ_mupl_normChi2_global[maxBranchSize];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_normChi2_global[maxBranchSize];   //[Reco_QQ_size]
Float_t         Reco_mu_normChi2_global[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_QQ_mupl_normChi2_global;   //!
TBranch        *b_Reco_QQ_mumi_normChi2_global;   //!
TBranch        *b_Reco_mu_normChi2_global;   //!
Int_t           Reco_QQ_mupl_nMuValHits[maxBranchSize];   //[Reco_QQ_size]
Int_t           Reco_QQ_mumi_nMuValHits[maxBranchSize];   //[Reco_QQ_size]
Int_t           Reco_mu_nMuValHits[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_QQ_mupl_nMuValHits;   //!
TBranch        *b_Reco_QQ_mumi_nMuValHits;   //!
TBranch        *b_Reco_mu_nMuValHits;   //!
Int_t           Reco_QQ_mupl_StationsMatched[maxBranchSize];   //[Reco_QQ_size]
Int_t           Reco_QQ_mumi_StationsMatched[maxBranchSize];   //[Reco_QQ_size]
Int_t           Reco_mu_StationsMatched[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_QQ_mupl_StationsMatched;   //!
TBranch        *b_Reco_QQ_mumi_StationsMatched;   //!
TBranch        *b_Reco_mu_StationsMatched;   //!
Float_t         Reco_QQ_mupl_dxy[maxBranchSize];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_dxy[maxBranchSize];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_dxyErr[maxBranchSize];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_dxyErr[maxBranchSize];   //[Reco_QQ_size]
Float_t         Reco_mu_dxy[maxBranchSize];   //[Reco_mu_size]
Float_t         Reco_mu_dxyErr[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_QQ_mupl_dxy;   //!
TBranch        *b_Reco_QQ_mumi_dxy;   //!
TBranch        *b_Reco_QQ_mupl_dxyErr;   //!
TBranch        *b_Reco_QQ_mumi_dxyErr;   //!
TBranch        *b_Reco_mu_dxy;   //!
TBranch        *b_Reco_mu_dxyErr;   //!
Float_t         Reco_QQ_mupl_dz[maxBranchSize];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_dz[maxBranchSize];   //[Reco_QQ_size]
Float_t         Reco_QQ_mupl_dzErr[maxBranchSize];   //[Reco_QQ_size]
Float_t         Reco_QQ_mumi_dzErr[maxBranchSize];   //[Reco_QQ_size]
Float_t         Reco_mu_dz[maxBranchSize];   //[Reco_mu_size]
Float_t         Reco_mu_dzErr[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_QQ_mupl_dz;   //!
TBranch        *b_Reco_QQ_mumi_dz;   //!
TBranch        *b_Reco_QQ_mupl_dzErr;   //!
TBranch        *b_Reco_QQ_mumi_dzErr;   //!
TBranch        *b_Reco_mu_dz;   //!
TBranch        *b_Reco_mu_dzErr;   //!
Int_t           Reco_QQ_mupl_nTrkWMea[maxBranchSize];   //[Reco_QQ_size]
Int_t           Reco_QQ_mumi_nTrkWMea[maxBranchSize];   //[Reco_QQ_size]
Int_t           Reco_mu_nTrkWMea[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_QQ_mupl_nTrkWMea;   //!
TBranch        *b_Reco_QQ_mumi_nTrkWMea;   //!
TBranch        *b_Reco_mu_nTrkWMea;   //!
Bool_t          Reco_QQ_mupl_TMOneStaTight[maxBranchSize];   //[Reco_QQ_size]
Bool_t          Reco_QQ_mumi_TMOneStaTight[maxBranchSize];   //[Reco_QQ_size]
Bool_t          Reco_mu_TMOneStaTight[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_QQ_mupl_TMOneStaTight;   //!
TBranch        *b_Reco_QQ_mumi_TMOneStaTight;   //!
TBranch        *b_Reco_mu_TMOneStaTight;   //!
Int_t           Reco_QQ_mupl_nPixWMea[maxBranchSize];   //[Reco_QQ_size]
Int_t           Reco_QQ_mumi_nPixWMea[maxBranchSize];   //[Reco_QQ_size]
Int_t           Reco_mu_nPixWMea[maxBranchSize];   //[Reco_mu_size]
TBranch        *b_Reco_QQ_mupl_nPixWMea;   //!
TBranch        *b_Reco_QQ_mumi_nPixWMea;   //!
TBranch        *b_Reco_mu_nPixWMea;   //!
Int_t           Reco_QQ_sign[maxBranchSize];   //[Reco_QQ_size]
TBranch        *b_Reco_QQ_sign;   //!

Int_t           Reco_QQ_mupl_nPixValHits[maxBranchSize];   //[Reco_QQ_size]
TBranch        *b_Reco_QQ_mupl_nPixValHits;   //!
Int_t           Reco_QQ_mumi_nPixValHits[maxBranchSize];   //[Reco_QQ_size]
TBranch        *b_Reco_QQ_mumi_nPixValHits;   //!
Float_t         Reco_QQ_mupl_ptErr_global[maxBranchSize];   //[Reco_QQ_size]
TBranch        *b_Reco_QQ_mupl_ptErr_global;   //!
Float_t         Reco_QQ_mumi_ptErr_global[maxBranchSize];   //[Reco_QQ_size]
TBranch        *b_Reco_QQ_mumi_ptErr_global;   //!


Float_t         Reco_QQ_ctau[maxBranchSize];
Float_t         Reco_QQ_ctau3D[maxBranchSize];
TBranch        *b_Reco_QQ_ctau;
TBranch        *b_Reco_QQ_ctau3D;



/////////////////////////////////////////
////// Gen QQ 
/////////////////////////////////////////
Int_t           Gen_QQ_size;
Int_t           Gen_QQ_type[maxBranchSize];   //[Gen_QQ_size]
TClonesArray    *Gen_QQ_4mom;
TClonesArray    *Gen_QQ_mupl_4mom;
TClonesArray    *Gen_QQ_mumi_4mom;
TBranch        *b_Gen_QQ_size;   //!
TBranch        *b_Gen_QQ_type;   //!
TBranch        *b_Gen_QQ_4mom;   //!
TBranch        *b_Gen_QQ_mupl_4mom;   //!
TBranch        *b_Gen_QQ_mumi_4mom;   //!



//~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~//
//~*~*~*~*~*~*~*Conversion~*~*~*~*~*~*~*~*~//
//~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~//

const int MaxConvBranchSize = 1000;
Int_t          Reco_conv_size;                                                    
TClonesArray   *Reco_conv_4mom;                                                           
Int_t          Reco_conv_validVtx[MaxConvBranchSize];                         
Float_t        Reco_conv_chi2[MaxConvBranchSize];                           
Float_t        Reco_conv_chi2_probability[MaxConvBranchSize];         
Float_t        Reco_conv_vtx_xErr[MaxConvBranchSize];             
Float_t        Reco_conv_vtx_yErr[MaxConvBranchSize];                       
Float_t        Reco_conv_vtx_zErr[MaxConvBranchSize];           
Float_t        Reco_conv_ntracks[MaxConvBranchSize];                        
Float_t        Reco_conv_paircotthetasep[MaxConvBranchSize];            
Float_t        Reco_conv_zofprimvtxfromtrks[MaxConvBranchSize]; 
Float_t        Reco_conv_distofminapproach[MaxConvBranchSize];        
Float_t        Reco_conv_dphitrksatvtx[MaxConvBranchSize];        
Float_t        Reco_conv_dxy[MaxConvBranchSize];                            
Float_t        Reco_conv_dz[MaxConvBranchSize];                       
Float_t        Reco_conv_lxy[MaxConvBranchSize];                          
Float_t        Reco_conv_lz[MaxConvBranchSize];                     
Int_t          Reco_conv_nSharedHits[MaxConvBranchSize];        
Int_t          Reco_conv_tk1_charge[MaxConvBranchSize];             
Int_t          Reco_conv_tk2_charge[MaxConvBranchSize];             
Int_t          Reco_conv_tk1_nhit[MaxConvBranchSize];           
Float_t        Reco_conv_tk1_dz[MaxConvBranchSize];             
Float_t        Reco_conv_tk1_dzerr[MaxConvBranchSize];                    
Float_t        Reco_conv_tk1_pterr[MaxConvBranchSize];              
Float_t        Reco_conv_tk1_etaerr[MaxConvBranchSize];         
Float_t        Reco_conv_tk1_thetaerr[MaxConvBranchSize];           
Float_t        Reco_conv_tk1_phierr[MaxConvBranchSize];           
Float_t        Reco_conv_tk1_lambdaerr[MaxConvBranchSize];          
Float_t        Reco_conv_tk1_d0[MaxConvBranchSize];                 
Float_t        Reco_conv_tk1_pout[MaxConvBranchSize];               
Float_t        Reco_conv_tk1_pin[MaxConvBranchSize];            
Int_t          Reco_conv_tk2_nhit[MaxConvBranchSize];
Float_t        Reco_conv_tk2_dz[MaxConvBranchSize];               
Float_t        Reco_conv_tk2_dzerr[MaxConvBranchSize];        
Float_t        Reco_conv_tk2_pterr[MaxConvBranchSize];          
Float_t        Reco_conv_tk2_etaerr[MaxConvBranchSize];                       
Float_t        Reco_conv_tk2_thetaerr[MaxConvBranchSize];                     
Float_t        Reco_conv_tk2_phierr[MaxConvBranchSize];                   
Float_t        Reco_conv_tk2_lambdaerr[MaxConvBranchSize];          
Float_t        Reco_conv_tk2_d0[MaxConvBranchSize];               
Float_t        Reco_conv_tk2_pout[MaxConvBranchSize];           
Float_t        Reco_conv_tk2_pin[MaxConvBranchSize];                    
TClonesArray   *Reco_conv_vtx;                                          
vector<vector <unsigned short> > *Reco_conv_nHitsBeforeVtx;   
vector<vector <int> > *Reco_conv_quality;                         



TBranch   *b_Reco_conv_size;
TBranch   *b_Reco_conv_4mom;
TBranch   *b_Reco_conv_validVtx;
TBranch   *b_Reco_conv_chi2;
TBranch   *b_Reco_conv_chi2_probability;
TBranch   *b_Reco_conv_vtx_xErr;
TBranch   *b_Reco_conv_vtx_yErr;
TBranch   *b_Reco_conv_vtx_zErr;
TBranch   *b_Reco_conv_ntracks;
TBranch   *b_Reco_conv_paircotthetasep;
TBranch   *b_Reco_conv_zofprimvtxfromtrks;
TBranch   *b_Reco_conv_distofminapproach;
TBranch   *b_Reco_conv_dphitrksatvtx;
TBranch   *b_Reco_conv_dxy;
TBranch   *b_Reco_conv_dz;
TBranch   *b_Reco_conv_lxy;
TBranch   *b_Reco_conv_lz;
TBranch   *b_Reco_conv_nHitsBeforeVtx;
TBranch   *b_Reco_conv_quality;
TBranch   *b_Reco_conv_nSharedHits;
TBranch   *b_Reco_conv_tk1_charge;
TBranch   *b_Reco_conv_tk2_charge;
TBranch   *b_Reco_conv_tk1_nhit;
TBranch   *b_Reco_conv_tk1_dz;
TBranch   *b_Reco_conv_tk1_dzerr;
TBranch   *b_Reco_conv_tk1_pterr;
TBranch   *b_Reco_conv_tk1_etaerr;
TBranch   *b_Reco_conv_tk1_thetaerr;
TBranch   *b_Reco_conv_tk1_phierr;
TBranch   *b_Reco_conv_tk1_lambdaerr;
TBranch   *b_Reco_conv_tk1_d0;
TBranch   *b_Reco_conv_tk1_pout;
TBranch   *b_Reco_conv_tk1_pin;
TBranch   *b_Reco_conv_tk2_nhit;
TBranch   *b_Reco_conv_tk2_dz;
TBranch   *b_Reco_conv_tk2_dzerr;
TBranch   *b_Reco_conv_tk2_pterr;
TBranch   *b_Reco_conv_tk2_etaerr;
TBranch   *b_Reco_conv_tk2_thetaerr;
TBranch   *b_Reco_conv_tk2_phierr;
TBranch   *b_Reco_conv_tk2_lambdaerr;
TBranch   *b_Reco_conv_tk2_d0;
TBranch   *b_Reco_conv_tk2_pout;
TBranch   *b_Reco_conv_tk2_pin;
TBranch   *b_Reco_conv_vtx;


class SetTree
{
  public:
    SetTree(){};

    virtual ~SetTree();
    virtual void TreeSetting(TTree* tree, bool isMC, bool isConversion);
    Bool_t SoftMuIdCut(int irqq);
};

SetTree::~SetTree()
{
}


void SetTree::TreeSetting(TTree* tree, bool isMC, bool isConversion)
{
  Reco_QQ_4mom = 0;
  Reco_QQ_mupl_4mom = 0;
  Reco_QQ_mumi_4mom = 0;

  tree->SetBranchAddress("runNb", &runNb, &b_runNb);
  tree->SetBranchAddress("LS", &LS, &b_LS);
  tree->SetBranchAddress("eventNb", &eventNb, &b_eventNb);
  tree->SetBranchAddress("zVtx", &zVtx, &b_zVtx);
  tree->SetBranchAddress("Centrality", &Centrality, &b_Centrality);
  tree->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
  tree->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
  tree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
  tree->SetBranchAddress("Reco_QQ_mupl_4mom", &Reco_QQ_mupl_4mom, &b_Reco_QQ_mupl_4mom);
  tree->SetBranchAddress("Reco_QQ_mumi_4mom", &Reco_QQ_mumi_4mom, &b_Reco_QQ_mumi_4mom);
  tree->SetBranchAddress("Reco_QQ_trig", Reco_QQ_trig, &b_Reco_QQ_trig);
  tree->SetBranchAddress("Reco_QQ_VtxProb", Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
  tree->SetBranchAddress("Reco_QQ_mupl_highPurity", Reco_QQ_mupl_highPurity, &b_Reco_QQ_mupl_highPurity);
  tree->SetBranchAddress("Reco_QQ_mumi_highPurity", Reco_QQ_mumi_highPurity, &b_Reco_QQ_mumi_highPurity);
  tree->SetBranchAddress("Reco_QQ_mumi_ptErr_global", Reco_QQ_mumi_ptErr_global, &b_Reco_QQ_mumi_ptErr_global);
  tree->SetBranchAddress("Reco_QQ_mupl_ptErr_global", Reco_QQ_mupl_ptErr_global, &b_Reco_QQ_mupl_ptErr_global);
  tree->SetBranchAddress("Reco_QQ_mupl_nPixValHits", Reco_QQ_mupl_nPixValHits, &b_Reco_QQ_mupl_nPixValHits);
  tree->SetBranchAddress("Reco_QQ_mumi_nPixValHits", Reco_QQ_mumi_nPixValHits, &b_Reco_QQ_mumi_nPixValHits);
  tree->SetBranchAddress("Reco_QQ_mupl_nTrkHits", Reco_QQ_mupl_nTrkHits, &b_Reco_QQ_mupl_nTrkHits);
  tree->SetBranchAddress("Reco_QQ_mumi_nTrkHits", Reco_QQ_mumi_nTrkHits, &b_Reco_QQ_mumi_nTrkHits);
  tree->SetBranchAddress("Reco_mu_nTrkHits", Reco_mu_nTrkHits, &b_Reco_mu_nTrkHits);
  tree->SetBranchAddress("Reco_QQ_mupl_normChi2_global", Reco_QQ_mupl_normChi2_global, &b_Reco_QQ_mupl_normChi2_global);
  tree->SetBranchAddress("Reco_QQ_mumi_normChi2_global", Reco_QQ_mumi_normChi2_global, &b_Reco_QQ_mumi_normChi2_global);
  tree->SetBranchAddress("Reco_mu_normChi2_global", Reco_mu_normChi2_global, &b_Reco_mu_normChi2_global);
  tree->SetBranchAddress("Reco_QQ_mupl_nMuValHits", Reco_QQ_mupl_nMuValHits, &b_Reco_QQ_mupl_nMuValHits);
  tree->SetBranchAddress("Reco_QQ_mumi_nMuValHits", Reco_QQ_mumi_nMuValHits, &b_Reco_QQ_mumi_nMuValHits);
  tree->SetBranchAddress("Reco_mu_nMuValHits", Reco_mu_nMuValHits, &b_Reco_mu_nMuValHits);
  tree->SetBranchAddress("Reco_QQ_mupl_StationsMatched", Reco_QQ_mupl_StationsMatched, &b_Reco_QQ_mupl_StationsMatched);
  tree->SetBranchAddress("Reco_QQ_mumi_StationsMatched", Reco_QQ_mumi_StationsMatched, &b_Reco_QQ_mumi_StationsMatched);
  tree->SetBranchAddress("Reco_mu_StationsMatched", Reco_mu_StationsMatched, &b_Reco_mu_StationsMatched);
  tree->SetBranchAddress("Reco_QQ_mupl_dxy", Reco_QQ_mupl_dxy, &b_Reco_QQ_mupl_dxy);
  tree->SetBranchAddress("Reco_QQ_mumi_dxy", Reco_QQ_mumi_dxy, &b_Reco_QQ_mumi_dxy);
  tree->SetBranchAddress("Reco_QQ_mupl_dxyErr", Reco_QQ_mupl_dxyErr, &b_Reco_QQ_mupl_dxyErr);
  tree->SetBranchAddress("Reco_QQ_mumi_dxyErr", Reco_QQ_mumi_dxyErr, &b_Reco_QQ_mumi_dxyErr);
  tree->SetBranchAddress("Reco_mu_dxy", Reco_mu_dxy, &b_Reco_mu_dxy);
  tree->SetBranchAddress("Reco_mu_dxyErr", Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
  tree->SetBranchAddress("Reco_QQ_mupl_dz", Reco_QQ_mupl_dz, &b_Reco_QQ_mupl_dz);
  tree->SetBranchAddress("Reco_QQ_mumi_dz", Reco_QQ_mumi_dz, &b_Reco_QQ_mumi_dz);
  tree->SetBranchAddress("Reco_QQ_mupl_dzErr", Reco_QQ_mupl_dzErr, &b_Reco_QQ_mupl_dzErr);
  tree->SetBranchAddress("Reco_QQ_mumi_dzErr", Reco_QQ_mumi_dzErr, &b_Reco_QQ_mumi_dzErr);
  tree->SetBranchAddress("Reco_mu_dz", Reco_mu_dz, &b_Reco_mu_dz);
  tree->SetBranchAddress("Reco_mu_dzErr", Reco_mu_dzErr, &b_Reco_mu_dzErr);
  tree->SetBranchAddress("Reco_QQ_mupl_TMOneStaTight", Reco_QQ_mupl_TMOneStaTight, &b_Reco_QQ_mupl_TMOneStaTight);
  tree->SetBranchAddress("Reco_QQ_mumi_TMOneStaTight", Reco_QQ_mumi_TMOneStaTight, &b_Reco_QQ_mumi_TMOneStaTight);
  tree->SetBranchAddress("Reco_mu_TMOneStaTight", Reco_mu_TMOneStaTight, &b_Reco_mu_TMOneStaTight);
  tree->SetBranchAddress("Reco_QQ_mupl_nPixWMea", Reco_QQ_mupl_nPixWMea, &b_Reco_QQ_mupl_nPixWMea);
  tree->SetBranchAddress("Reco_QQ_mumi_nPixWMea", Reco_QQ_mumi_nPixWMea, &b_Reco_QQ_mumi_nPixWMea);
  tree->SetBranchAddress("Reco_mu_nPixWMea", Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
  tree->SetBranchAddress("Reco_QQ_sign", Reco_QQ_sign, &b_Reco_QQ_sign);
  tree->SetBranchAddress("Reco_QQ_mupl_nTrkWMea", Reco_QQ_mupl_nTrkWMea, &b_Reco_QQ_mupl_nTrkWMea);
  tree->SetBranchAddress("Reco_QQ_mumi_nTrkWMea", Reco_QQ_mumi_nTrkWMea, &b_Reco_QQ_mumi_nTrkWMea);
  tree->SetBranchAddress("Reco_mu_nTrkWMea", Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
  tree->SetBranchAddress("Reco_QQ_ctau",Reco_QQ_ctau,&b_Reco_QQ_ctau);
  tree->SetBranchAddress("Reco_QQ_ctau3D",Reco_QQ_ctau3D,&b_Reco_QQ_ctau3D);

  if (isMC) { 
    Gen_QQ_4mom = 0;
    Gen_QQ_mupl_4mom = 0;
    Gen_QQ_mumi_4mom = 0;
    tree->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
    tree->SetBranchAddress("Gen_QQ_type", Gen_QQ_type, &b_Gen_QQ_type);
    tree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
    tree->SetBranchAddress("Gen_QQ_mupl_4mom", &Gen_QQ_mupl_4mom, &b_Gen_QQ_mupl_4mom);
    tree->SetBranchAddress("Gen_QQ_mumi_4mom", &Gen_QQ_mumi_4mom, &b_Gen_QQ_mumi_4mom);
  }

  if(isConversion) {
    Reco_conv_4mom=0;
    Reco_conv_vtx=0;
    tree->SetBranchAddress("Reco_conv_size", &Reco_conv_size, &b_Reco_conv_size);
    tree->SetBranchAddress("Reco_conv_4mom", &Reco_conv_4mom, &b_Reco_conv_4mom);
    tree->SetBranchAddress("Reco_conv_validVtx",Reco_conv_validVtx ,&b_Reco_conv_validVtx);
    tree->SetBranchAddress("Reco_conv_chi2",Reco_conv_chi2 ,&b_Reco_conv_chi2);
    tree->SetBranchAddress("Reco_conv_chi2_probability",Reco_conv_chi2_probability ,&b_Reco_conv_chi2_probability);
    tree->SetBranchAddress("Reco_conv_vtx_xErr",Reco_conv_vtx_xErr ,&b_Reco_conv_vtx_xErr);
    tree->SetBranchAddress("Reco_conv_vtx_yErr",Reco_conv_vtx_yErr ,&b_Reco_conv_vtx_yErr);
    tree->SetBranchAddress("Reco_conv_vtx_zErr",Reco_conv_vtx_zErr ,&b_Reco_conv_vtx_zErr);
    tree->SetBranchAddress("Reco_conv_ntracks",Reco_conv_ntracks ,&b_Reco_conv_ntracks);
    tree->SetBranchAddress("Reco_conv_paircotthetasep",Reco_conv_paircotthetasep ,&b_Reco_conv_paircotthetasep);
    tree->SetBranchAddress("Reco_conv_zofprimvtxfromtrks",Reco_conv_zofprimvtxfromtrks,&b_Reco_conv_zofprimvtxfromtrks);
    tree->SetBranchAddress("Reco_conv_distofminapproach",Reco_conv_distofminapproach ,&b_Reco_conv_distofminapproach);
    tree->SetBranchAddress("Reco_conv_dphitrksatvtx",Reco_conv_dphitrksatvtx ,&b_Reco_conv_dphitrksatvtx);
    tree->SetBranchAddress("Reco_conv_dxy",Reco_conv_dxy ,&b_Reco_conv_distofminapproach);
    tree->SetBranchAddress("Reco_conv_dz",Reco_conv_dz ,&b_Reco_conv_dz);
    tree->SetBranchAddress("Reco_conv_lxy",Reco_conv_lxy ,&b_Reco_conv_lxy);
    tree->SetBranchAddress("Reco_conv_lz",Reco_conv_lz ,&b_Reco_conv_lz);
    tree->SetBranchAddress("Reco_conv_nSharedHits",Reco_conv_nSharedHits ,&b_Reco_conv_nSharedHits);
    tree->SetBranchAddress("Reco_conv_tk1_charge",Reco_conv_tk1_charge,&b_Reco_conv_tk1_charge);
    tree->SetBranchAddress("Reco_conv_tk1_nhit",Reco_conv_tk1_nhit ,&b_Reco_conv_tk1_nhit);
    tree->SetBranchAddress("Reco_conv_tk1_dz",Reco_conv_tk1_dz,&b_Reco_conv_tk1_dz);
    tree->SetBranchAddress("Reco_conv_tk1_dzerr",Reco_conv_tk1_dzerr ,&b_Reco_conv_tk1_dzerr);
    tree->SetBranchAddress("Reco_conv_tk1_pterr",Reco_conv_tk1_pterr ,&b_Reco_conv_tk1_pterr);
    tree->SetBranchAddress("Reco_conv_tk1_etaerr",Reco_conv_tk1_etaerr ,&b_Reco_conv_tk1_etaerr);
    tree->SetBranchAddress("Reco_conv_tk1_thetaerr",Reco_conv_tk1_thetaerr ,&b_Reco_conv_tk1_thetaerr);
    tree->SetBranchAddress("Reco_conv_tk1_phierr",Reco_conv_tk1_phierr ,&b_Reco_conv_tk1_phierr);
    tree->SetBranchAddress("Reco_conv_tk1_lambdaerr",Reco_conv_tk1_lambdaerr ,&b_Reco_conv_tk1_lambdaerr);
    tree->SetBranchAddress("Reco_conv_tk1_d0",Reco_conv_tk1_d0 ,&b_Reco_conv_tk1_d0);
    tree->SetBranchAddress("Reco_conv_tk1_pout",Reco_conv_tk1_pout ,&b_Reco_conv_tk1_pout);
    tree->SetBranchAddress("Reco_conv_tk1_pin",Reco_conv_tk1_pin ,&b_Reco_conv_tk1_pin);
    tree->SetBranchAddress("Reco_conv_tk2_charge",Reco_conv_tk2_charge,&b_Reco_conv_tk2_charge);
    tree->SetBranchAddress("Reco_conv_tk2_nhit",Reco_conv_tk2_nhit ,&b_Reco_conv_tk2_nhit);
    tree->SetBranchAddress("Reco_conv_tk2_dz",Reco_conv_tk2_dz ,&b_Reco_conv_tk2_dz);
    tree->SetBranchAddress("Reco_conv_tk2_dzerr",Reco_conv_tk2_dzerr ,&b_Reco_conv_tk2_dzerr);
    tree->SetBranchAddress("Reco_conv_tk2_pterr",Reco_conv_tk2_pterr ,&b_Reco_conv_tk2_pterr);
    tree->SetBranchAddress("Reco_conv_tk2_etaerr",Reco_conv_tk2_etaerr ,&b_Reco_conv_tk2_etaerr);
    tree->SetBranchAddress("Reco_conv_tk2_thetaerr",Reco_conv_tk2_thetaerr ,&b_Reco_conv_tk2_thetaerr);
    tree->SetBranchAddress("Reco_conv_tk2_phierr",Reco_conv_tk2_phierr ,&b_Reco_conv_tk2_phierr);
    tree->SetBranchAddress("Reco_conv_tk2_lambdaerr",Reco_conv_tk2_lambdaerr ,&b_Reco_conv_tk2_lambdaerr);
    tree->SetBranchAddress("Reco_conv_tk2_d0",Reco_conv_tk2_d0 ,&b_Reco_conv_tk2_d0);
    tree->SetBranchAddress("Reco_conv_tk2_pout",Reco_conv_tk2_pout ,&b_Reco_conv_tk2_pout);
    tree->SetBranchAddress("Reco_conv_tk2_pin",Reco_conv_tk2_pin ,&b_Reco_conv_tk2_pin);
    tree->SetBranchAddress("Reco_conv_vtx",&Reco_conv_vtx ,&b_Reco_conv_vtx);
    tree->SetBranchAddress("Reco_conv_nHitsBeforeVtx",&Reco_conv_nHitsBeforeVtx ,&b_Reco_conv_nHitsBeforeVtx);
    tree->SetBranchAddress("Reco_conv_quality",&Reco_conv_quality ,&b_Reco_conv_quality);
  }

};

Bool_t SetTree::SoftMuIdCut(int irqq)
{
  bool muplSoft = ( (Reco_QQ_mupl_TMOneStaTight[irqq]==true) &&
      (Reco_QQ_mupl_nTrkWMea[irqq] > 5) &&
      (Reco_QQ_mupl_nPixWMea[irqq] > 0) &&
      (Reco_QQ_mupl_dxy[irqq]<0.3) &&
      (Reco_QQ_mupl_dz[irqq]<20.)   			//			 &&  (Reco_QQ_mupl_highPurity[irqq]==true) 
      ) ; 

  bool mumiSoft = ( (Reco_QQ_mumi_TMOneStaTight[irqq]==true) &&
      (Reco_QQ_mumi_nTrkWMea[irqq] > 5) &&
      (Reco_QQ_mumi_nPixWMea[irqq] > 0) &&
      (Reco_QQ_mumi_dxy[irqq]<0.3) &&
      (Reco_QQ_mumi_dz[irqq]<20.)  // && (Reco_QQ_mumi_highPurity[irqq]==true)
      ) ; 

  return (muplSoft && mumiSoft);
}  



#endif 
