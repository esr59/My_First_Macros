//Esha Rao
//June 15, 2017
//Rutgers University

//test macro to read hiForest and plot jet AND track variables with conditions for each


//Adapted from:

//Jennifer Coulter
// May 19th 2015
// Rutgers, jbc120@scarletmail.rutgers.edu
//
//
// test macro to read hiForest and plot jet variables
//
//May 19th 2015 --> added for loop rendition of histogram plotting
//  pushed to git
//
//
//
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <TH1F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TChain.h>
#include <TProfile.h>
#include <TStopwatch.h>
#include <TCut.h>
#include <cstdlib>
#include <cmath>
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"

using namespace std;

double deltaR(double eta1, double phi1, double eta2, double phi2){

  double delR = 0.0;
  double deleta = eta1 - eta2;
  double delphi = phi1 - phi2;

  if(delphi >= 3.1415)
    delphi = delphi - 3.1415;

  delR = sqrt(deleta*deleta + delphi*delphi);

  return delR;
  
}


void read_HiForest_jetsandtracks(Int_t radius = 3, char * algo = (char*)"PF"){
 
  bool printDebug = false;

  //Define file and tree
  TFile * fin = TFile::Open("/Users/rao/Root/Root_Files/HiForestAOD_ppMC_Py8_1.root");
  TFile * new_File = new TFile("/Users/rao/Root/Root_Files/test_pp_jetVariable_plots.root","RECREATE");
  TTree * jet = (TTree*)fin->Get(Form("ak%d%sJetAnalyzer/t",radius,algo));
  TTree * track = (TTree*)fin->Get(Form("ppTrack/trackTree"));
  jet->AddFriend(track);

  //Set branches of the tree for Track
  Float_t tpt[1000];
  Float_t teta[1000];
  Float_t tphi[1000];
  Int_t ntrk;
  
  track->SetBranchAddress("trkPt", &tpt);
  track->SetBranchAddress("trkEta", &teta);
  track->SetBranchAddress("trkPhi", &tphi);
  track->SetBranchAddress("nTrk", &ntrk);
  
  //Part 2: Drawing Using Loops

  //Definitions
  TH1F * hJet_pT = new TH1F("hJet_pT", "", 100,0,1000);
  TH1F * hJet_eta = new TH1F("hJet_eta", "", 60,-3,+3);
  TH1F * hJet_phi = new TH1F("hJet_phi", "", 60,-3.15,+3.15);
  TH2F * hJet_2d = new TH2F("hJet_2d","Jet pT vs Eta",30,-3,+3,100,0,500);
  TH1F * hTrk_pT = new TH1F("hTrk_pT", "", 50,0,500);
  TH1F * hTrk_eta = new TH1F("hTrk_eta", "", 60,-3,+3);
  TH1F * hTrk_phi = new TH1F("hTrk_phi", "", 60,-3.15,+3.15);
  TH2F * hTrk_2d = new TH2F("hTrk_2d","Track pT vs Track Eta",30,-3,+3,100,0,500);
   
  Float_t pt[1000];
  Float_t eta[1000];
  Float_t phi[1000];
  Int_t nref;
  
  //Set branches of the tree for Jet
  jet->SetBranchAddress("jtpt", &pt);
  jet->SetBranchAddress("jteta", &eta);
  jet->SetBranchAddress("jtphi", &phi);
  jet->SetBranchAddress("nref", &nref);

  Long64_t nentries =  jet->GetEntries();
  cout<< "Number of events "<<jet->GetEntries()<<endl;

  if(printDebug) nentries = 10;
  
  //Start entry loop
  for(Long64_t nentry = 0; nentry < nentries; ++nentry){
    
    jet->GetEvent(nentry);
    //if(nentry%1000 == 0) cout << nentry << "/" << nentries <<endl;
    
    //start the jet loop
    for(int jentry = 0; jentry < nref; ++jentry){

      if(fabs(eta[jentry])>2.0)
	continue;

      hJet_pT->Fill(pt[jentry]);
      if(printDebug) cout<<"pt = "<< pt[jentry] << endl;

      hJet_eta->Fill(eta[jentry]);
      if(printDebug) cout<<"eta = "<< eta[jentry] << endl;

      hJet_phi->Fill(phi[jentry]);
      if(printDebug) cout<<"phi = "<< phi[jentry] << endl;

      hJet_2d->Fill(eta[jentry],pt[jentry]);
      if(printDebug) cout<<"pt = "<< pt[jentry] << endl;
      if(printDebug) cout<<"eta  = "<< eta[jentry] << endl;

      //start the track loop
      for (int tentry = 0; tentry < ntrk; ++tentry){

	//! calculate delta R
	//double delR = sqrt(pow((tphi[tentry]-phi[jentry]),2)+pow((teta[tentry]-eta[jentry]),2));
	double delR = deltaR(teta[tentry], tphi[tentry], eta[jentry], phi[jentry]);
	
	//! check delR condition 
	if(delR <= 0.3) {

	  hTrk_pT->Fill(tpt[tentry]);
	  if(printDebug) cout<<"tpt = "<< pt[tentry] << endl;
	  
	  hTrk_eta->Fill(teta[tentry]);
	  if(printDebug) cout<<"teta = "<< eta[tentry] << endl;
	
	  hTrk_phi->Fill(tphi[tentry]);
	  if(printDebug) cout<<"tphi = "<< tphi[tentry] << endl;
	  
	  hTrk_2d->Fill(teta[tentry],tpt[tentry]);
	  if(printDebug) cout<<"tpt = "<< tpt[tentry] << endl;
	  if(printDebug) cout<<"teta  = "<< teta[tentry] << endl;

	}//! delR condition check
	
      }//end track loop
    }//end jet loop
   
  }//end entry loop

  // check histograms (Prints histogram information)
  hJet_2d->Print("base");
  cout<<"histogram mean = "<<hJet_pT->GetMean()<<endl;
  
  // plotting for jets
  TCanvas * c1 = new TCanvas("c1", "Simple Track Variables", 1200, 1200);
  c1->Divide(2,2);
  
  c1->cd(1)->SetLogy();
  hJet_pT->SetTitle("pp Jet pT");
  hJet_pT->SetXTitle("Jet pT (GeV/c)");
  hJet_pT->SetYTitle("Counts");
  hJet_pT->Draw();

  c1->cd(2);
  hJet_eta->SetTitle("pp Jet Eta");
  hJet_eta->SetXTitle("Jet Eta");
  hJet_eta->SetYTitle("Counts");
  hJet_eta->Draw();
  
  c1->cd(3);
  hJet_phi->SetTitle("pp Jet Phi");
  hJet_phi->SetXTitle("Jet Phi (Rad)");
  hJet_phi->SetYTitle("Counts");
  hJet_phi->Draw();

  c1->cd(4)->SetLogz();
  hJet_2d->SetTitle("pp Jet Eta vs. pT");
  hJet_2d->SetXTitle("Jet Eta");
  hJet_2d->SetYTitle("Jet pT (GeV/c)");
  hJet_2d->Draw("colz");

  //plotting for tracks
  TCanvas * c2 = new TCanvas("c2", "Simple Track Variables", 1200, 1200);
  c2->Divide(2,2);
  
  c2->cd(1)->SetLogy();
  hTrk_pT->SetTitle("pp Track pT");
  hTrk_pT->SetXTitle("Track pT (GeV/c)");
  hTrk_pT->SetYTitle("Counts");
  hTrk_pT->Draw();

  c2->cd(2);
  hTrk_eta->SetTitle("pp Track Eta");
  hTrk_eta->SetXTitle("Track Eta");
  hTrk_eta->SetYTitle("Counts");
  hTrk_eta->Draw();
  
  c2->cd(3);
  hTrk_phi->SetTitle("pp Track Phi");
  hTrk_phi->SetXTitle("Track Phi (Rad)");
  hTrk_phi->SetYTitle("Counts");
  hTrk_phi->Draw();

  c2->cd(4)->SetLogz();  hTrk_2d->SetTitle("pp Track Eta vs. pT");
  hTrk_2d->SetXTitle("Track Eta");
  hTrk_2d->SetYTitle("Track pT (GeV/c)");
  hTrk_2d->Draw("colz");

  //creates a pdf
  new_File->Write();
  c2->SaveAs("test_pp_jetvariables_with_fors.pdf","RECREATE");
 
}// end of macro
