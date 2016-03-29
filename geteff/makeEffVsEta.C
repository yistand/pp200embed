//====================================================================================================
//
//	2016.03.24	Li YI
//	eta effiency for pp@200GeV P12id
//
//====================================================================================================



#include "TCanvas.h"
#include "TSystem.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TString.h"
#include "TMath.h"
#include "math.h"
#include <iostream>

using namespace std;

double myfun(double *x, double *par) {
  double eta = x[0];
  double y = 0;
  
  double turnpoint = 0.5;
  double maxeta = 5.3;
  double p0 = par[0];
  double p1 = par[1];
  double p2 = par[2];
  double p3 = par[3];
  double p4 = par[4];
  double p5 = par[5];
  if(eta<-turnpoint&&eta>=-1) {
    //y =  par[0]+par[1]*eta+par[2]*eta*eta;
    y = p0*exp(-pow(p1/(eta+maxeta),p2));
  }
  if(fabs(eta)<=turnpoint) {
    //y = par[0]-par[1]*0.8+par[2]*0.64;
    //double k = ((par[3]+par[4]*0.8+par[5]*0.8*0.8)-(par[0]+par[1]*0.8+par[2]*0.8*0.8))/1.6;
    //double b = ((par[3]+par[4]*0.8+par[5]*0.8*0.8)+(par[0]+par[1]*0.8+par[2]*0.8*0.8))/2.;
    //y = k*eta+b;
    double k = ((p3*exp(-pow(p4/(maxeta-fabs(turnpoint)),p5)))-(p0*exp(-pow(p1/(-fabs(turnpoint)+maxeta),p2))))/(fabs(turnpoint)*2.);
    double b = ((p3*exp(-pow(p4/(maxeta-fabs(turnpoint)),p5)))+(p0*exp(-pow(p1/(-fabs(turnpoint)+maxeta),p2))))/2.;
    y = k*eta+b;
  }
  if(eta>turnpoint&&eta<=1) {
    //y = (par[0]-par[1]*0.8+par[2]*0.64-par[3]*0.8-par[4]*0.64)+par[3]*eta+par[4]*eta*eta;
    //y = par[3]+par[4]*eta+par[5]*eta*eta;
    y = p3*exp(-pow(p4/(maxeta-eta),p5));
  }

 return y;

}

int ZdcArray(float zdc) {

  int array = 0;
  array = int(zdc/1000)-3;

  if(array>=10) array = 9;
  if(array<0) array = 0;

  return array; 

}


int GlobalArray(float globaltrk) {

  int array = 0;
  array = int(globaltrk/100);

  if(array>=10) array = 9;

  return array; 

}


void makeEffVsEta(const char* particle, int energy=200){	

  static const Double_t pi = TMath::Pi();

  gStyle->SetOptStat(111111);
  gStyle->SetPalette(1);
  //gStyle->SetFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetFrameBorderMode(0);
  //gStyle->SetPadTopMargin(0.05);
  //gStyle->SetPadRightMargin(0.05);
  //gStyle->SetPadBottomMargin(0.05);
  gStyle->SetPadLeftMargin(0.11);
  
  int Nbins = 100;
  float max_z = 30;
  float min_pt = 0.;
  float max_pt = 20;
  float min_fitpts = 20; // 25;
  float min_fitpts_nposs = .52;
  float max_eta = 1.0;
  float max_dca = 1.0;
  float PID = 0;
  char cbuff[10];
  char buffer[100];

  sprintf(cbuff, "%d",energy);
  TString NRG = TString(cbuff);
  
  cout << "particle = " << particle << endl;
  cout << "energy = " << NRG << endl;
  
  TString out_file, in_file;
  TString prefix; 
  TString postfix; 
  prefix = TString("/projects/rhig/ly247/mcembedding/pp200_P12id");
  postfix = TString(".");
  double mass;
  if(TString(particle) == TString("Piplus")){
    in_file = prefix + TString("/SinglePiplusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/") + TString("EtaDep_Piplus") + NRG + TString("GeV.root");
    mass = .13975;
    PID = 8;
  }
  else if(TString(particle) == TString("Piminus")){
    in_file = prefix + TString("/SinglePiminusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/") + TString("EtaDep_Piminus") + NRG + TString("GeV.root");
    mass = .13975;
    PID = 9;
  }
  else if(TString(particle) == TString("Kplus")){
    in_file = prefix + TString("/SingleKplusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/") + TString("EtaDep_Kplus") + NRG + TString("GeV.root");
    mass = .493677;
    PID = 11;
  }
  else if(TString(particle) == TString("Kminus")){
    in_file = prefix + TString("/SingleKminusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/") + TString("EtaDep_Kminus") + NRG + TString("GeV.root");
    mass = .493677;
    PID = 12;
  }
  else if(TString(particle) == TString("Proton")){
    in_file = prefix + TString("/SingleProtonNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/") + TString("EtaDep_Proton") + NRG + TString("GeV.root");
    mass = .93827;
    PID = 14;
  }
  else if(TString(particle) == TString("AntiProton")){
    in_file = prefix + TString("/SingleAntiProtonNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/") + TString("EtaDep_AntiProton") + NRG + TString("GeV.root");
    mass = .93827;
    PID = 15;
  }
  

//------------------------- Input ---------------------------
  cout << "opening input file: " << in_file << endl;
  TFile *f = new TFile(in_file, "READ");
  TNtuple *MatchedPairs = (TNtuple*) f->Get("MatchedPairs_NT");
  TNtuple *McTrack = (TNtuple*) f->Get("McTrack_NT");
  
  float Dedx, RefMult, RefMultCorrected, CentralityWeight, Centrality16, VertexX, VertexY, VertexZ, PtMc, PzMc, EtaMc, PhiMc, PtPr, PzPr, EtaPr, PhiPr, DcaGl, DcaZGl, DcaXYGl, Flag, FitPts, DedxPts, AllPts, NPossible, ParentGeantId, GeantId,mErrP, GlobalTr, ZdcRate, BbcRate;
  MatchedPairs->SetBranchAddress("Dedx", &Dedx);
  MatchedPairs->SetBranchAddress("RefMult", &RefMult);
  MatchedPairs->SetBranchAddress("RefMultCorrected", &RefMultCorrected);
  MatchedPairs->SetBranchAddress("CentralityWeight", &CentralityWeight);
  MatchedPairs->SetBranchAddress("Centrality16", &Centrality16);
  MatchedPairs->SetBranchAddress("VertexX", &VertexX);
  MatchedPairs->SetBranchAddress("VertexY", &VertexY);
  MatchedPairs->SetBranchAddress("VertexZ", &VertexZ);
  MatchedPairs->SetBranchAddress("PtMc", &PtMc);
  MatchedPairs->SetBranchAddress("PzMc", &PzMc);
  MatchedPairs->SetBranchAddress("EtaMc", &EtaMc);
  MatchedPairs->SetBranchAddress("PhiMc", &PhiMc);
  MatchedPairs->SetBranchAddress("PtPr", &PtPr);
  MatchedPairs->SetBranchAddress("mErrP", &mErrP);
  //MatchedPairs->SetBranchAddress("PzPr", &PzPr);
  MatchedPairs->SetBranchAddress("EtaPr", &EtaPr);
  MatchedPairs->SetBranchAddress("PhiPr", &PhiPr);
  MatchedPairs->SetBranchAddress("DcaGl", &DcaGl);
  MatchedPairs->SetBranchAddress("DcaZGl", &DcaZGl);
  MatchedPairs->SetBranchAddress("DcaXYGl", &DcaXYGl);
  MatchedPairs->SetBranchAddress("Flag", &Flag);
  MatchedPairs->SetBranchAddress("FitPts", &FitPts);
  MatchedPairs->SetBranchAddress("DedxPts", &DedxPts);
  MatchedPairs->SetBranchAddress("AllPts", &AllPts);
  MatchedPairs->SetBranchAddress("NPossible", &NPossible);
  MatchedPairs->SetBranchAddress("ParentGeantId", &ParentGeantId);
  MatchedPairs->SetBranchAddress("GeantId", &GeantId);
  MatchedPairs->SetBranchAddress("GlobalTr", &GlobalTr);
  MatchedPairs->SetBranchAddress("ZdcRate", &ZdcRate);
  MatchedPairs->SetBranchAddress("BbcRate", &BbcRate);
  
  float pRefMult, pRefMultCorrected, pCentralityWeight, pCentrality16, pVertexX, pVertexY, pVertexZ, pPtMc, pPzMc, pEtaMc, pPhiMc, pParentGeantId, pGeantId, pGlobalTr, pZdcRate, pBbcRate;
  McTrack->SetBranchAddress("RefMult", &pRefMult);
  McTrack->SetBranchAddress("RefMultCorrected", &pRefMultCorrected);
  McTrack->SetBranchAddress("CentralityWeight", &pCentralityWeight);
  McTrack->SetBranchAddress("Centrality16", &pCentrality16);
  McTrack->SetBranchAddress("VertexX", &pVertexX);
  McTrack->SetBranchAddress("VertexY", &pVertexY);
  McTrack->SetBranchAddress("VertexZ", &pVertexZ);
  McTrack->SetBranchAddress("PtMc", &pPtMc);
  McTrack->SetBranchAddress("PzMc", &pPzMc);
  McTrack->SetBranchAddress("EtaMc", &pEtaMc);
  McTrack->SetBranchAddress("PhiMc", &pPhiMc);
  McTrack->SetBranchAddress("ParentGeantId", &pParentGeantId);
  McTrack->SetBranchAddress("GeantId", &pGeantId);
  McTrack->SetBranchAddress("GlobalTr", &pGlobalTr);
  McTrack->SetBranchAddress("ZdcRate", &pZdcRate);
  McTrack->SetBranchAddress("BbcRate", &pBbcRate);
  
//------------------------- Output ---------------------------
  cout << "opening output file: " << out_file << endl;
  TFile *f_out = new TFile(out_file, "RECREATE");
  f_out->cd();
  
  TH1D* heta_mc;
  TH1D* heta_re;
  TH1D* heta;

  heta_mc = new TH1D("heta_mc", Form("MC Tracks vs eta"), Nbins, -max_eta, max_eta);
  heta_mc->Sumw2();
  heta_re = new TH1D("heta_re", Form("Matched Tracks vs eta"), Nbins, -max_eta,  max_eta);
  heta_re->Sumw2();

  cout << "histograms created" << endl;
  
//------------------------- Loop ---------------------------
  for(int i = 0; i<McTrack->GetEntries(); i++){
    McTrack->GetEntry(i);
 
    if( pParentGeantId!=0 ) continue;
   
    if( pGeantId!=PID ) continue;

    if( fabs(pVertexZ) > max_z ) continue;

    if( fabs(pEtaMc) > max_eta ) continue;
    if( fabs(pPtMc) < min_pt || fabs(pPtMc) > max_pt) continue;

    heta_mc->Fill(pEtaMc);

  }
  

  for(int i = 0; i <MatchedPairs->GetEntries(); i++){
    MatchedPairs->GetEntry(i);
    if( ParentGeantId!=0 ) continue;
    if( GeantId!=PID ) continue;
    if( Flag < 0 ) continue;
    if( fabs(VertexZ) > max_z ) continue;
    if( fabs(PtMc)< min_pt || fabs(PtMc)>max_pt) continue;
    if( FitPts < min_fitpts ) continue;
    if( FitPts/NPossible < min_fitpts_nposs ) continue;
    if(fabs(DcaGl) > max_dca ) continue;//change DcaGl-->GcaXYGl
    if( fabs(EtaMc) > max_eta ) continue;
    
    heta_re->Fill(EtaMc);
    
  }

  heta = (TH1D*)heta_re->Clone(Form("heta"));
  heta->SetTitle(Form("eff vs eta"));
  heta->Divide(heta_mc);

  f_out->cd();
  heta_mc->Write();
  heta_re->Write();
  heta->Write();
  f_out->Close();

  f->Close();

}



