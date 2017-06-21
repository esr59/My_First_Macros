//Esha Rao
//June 15, 2017
//Rutgers University

//test macro to read hiForest and plot histogram of z. 1/z and log(1/z) for specific tracks of jets

// testing edit on github

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

void zHiForest(Int_t radius = 3, char * algo = (char*)"PF"){
 
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
  TH1F * hZ = new TH1F("hZ", "",10,0,1);
  TH1F * hZinv = new TH1F("hZinv","",10,1,100);
  TH1F * hlog = new TH1F("hlog","",10,1,5);
   
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

  Long64_t nJets = 0;
  
  //Start entry loop
  for(Long64_t nentry = 0; nentry < nentries; ++nentry){
    
    jet->GetEvent(nentry);
    //if(nentry%1000 == 0) cout << nentry << "/" << nentries <<endl;
    
    //start the jet loop
    for(int jentry = 0; jentry < nref; ++jentry){

      if(fabs(eta[jentry])>2.0)
	continue;
      if(pt[jentry]<=50 || pt[jentry]>100)
	continue;

      nJets++;
    
      if(printDebug) cout<<"pt = "<< pt[jentry] << endl;

     
      if(printDebug) cout<<"eta = "<< eta[jentry] << endl;

     
      if(printDebug) cout<<"phi = "<< phi[jentry] << endl;

    
      if(printDebug) cout<<"pt = "<< pt[jentry] << endl;
      if(printDebug) cout<<"eta  = "<< eta[jentry] << endl;

      //start the track loop
      for (int tentry = 0; tentry < ntrk; ++tentry){

	//! calculate delta R
	//double delR = sqrt(pow((tphi[tentry]-phi[jentry]),2)+pow((teta[tentry]-eta[jentry]),2));
	double delR = deltaR(teta[tentry], tphi[tentry], eta[jentry], phi[jentry]);
	
	//! check delR condition 
	if(delR <= 0.3 && tpt[tentry]>1.) {
	 
     

	  double z= (tpt[tentry]*cos(delR)/pt[jentry]);
	  double zinv=(1/z);
	  double logzinv=(log(zinv));
	 
	  
	  if(printDebug) cout<<"tpt = "<< pt[tentry] << endl;
	  
	 
	  if(printDebug) cout<<"teta = "<< eta[tentry] << endl;
	
	  if(printDebug) cout<<"tphi = "<< tphi[tentry] << endl;
	  
       
	  if(printDebug) cout<<"tpt = "<< tpt[tentry] << endl;
	  if(printDebug) cout<<"teta  = "<< teta[tentry] << endl;

	  hZ->Fill(z);
	  hZinv->Fill(zinv);
	  hlog->Fill(logzinv);
	 

	}//! delR condition check
	
      }//end track loop
    }//end jet loop
   
  }//end entry loop

  // check histograms (Prints histogram information)
  hZ->Print("base");
  cout<<"histogram mean = "<<hZ->GetMean()<<endl;
  
  // plotting for jets
  TCanvas * c1 = new TCanvas("c1", "The z Histogram", 1200, 1200);
  c1->Divide(3,1);
  
  c1->cd(1)->SetLogy();
  hZ->SetTitle("Histogram of (#frac{p^{trk}_{T}}{ p^{jet}_{T}})");
  // hZ->SetXTitle("z (Track pT/Jet pT)");
  hZ->SetXTitle("#frac{p^{trk}_{T}}{p^{jet}_{T}}");
  // hZ->SetYTitle("Number of Tracks/Number of Bins/Number of Jets");
  hZ->SetYTitle("1/N_{jets} #frac{dN}{dp_{T}}");
  hZ->GetYaxis()->CenterTitle();
  hZ->Scale(1./(.1*nJets));
  hZ->Draw();

  c1->cd(2)->SetLogy();
  hZinv->SetTitle("Histogram of #frac{1}{Z}");
  // hZ->SetXTitle("z (Track pT/Jet pT)");
  hZinv->SetXTitle("#frac{p^{jet}_{T}}{p^{trk}_{T}}");
  // hZ->SetYTitle("Number of Tracks/Number of Bins/Number of Jets");
  hZinv->SetYTitle("1/N_{jets} #frac{dN}{dp_{T}}");
  hZinv->GetYaxis()->CenterTitle();
  hZinv->Scale(1./(.1*nJets));
  hZinv->Draw();

  c1->cd(3)->SetLogy();
  hlog->SetTitle("Histogram of log(#frac{1}{Z})");
  // hZ->SetXTitle("z (Track pT/Jet pT)");
  hlog->SetXTitle("log(#frac{p^{jet}_{T}}{p^{trk}_{T}})");
  // hZ->SetYTitle("Number of Tracks/Number of Bins/Number of Jets");
  hlog->SetYTitle("1/N_{jets} #frac{dN}{dp_{T}}");
  hlog->GetYaxis()->CenterTitle();
   hlog->Scale(1./(.1*nJets));
  hlog->Draw();

  gStyle->SetOptStat(0);

  //creates a pdf
  new_File->Write();
  c1->SaveAs("test_pp_z","RECREATE");
 
}// end of macro
