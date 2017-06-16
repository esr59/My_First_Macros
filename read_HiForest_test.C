//Esha Rao
//Rutgers University
//June 12, 2017

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
#include "TLatex.h"
#include "TMath.h"
#include "TLine.h"
#include "TAttLine.h"
#include "TROOT.h"
#include "THStack.h"

using namespace std;

void read_HiForest_test(Int_t radius = 3,
		        std::string algo = "PF"){

  gStyle->SetOptStat(0); //get rid of default legends
 
  bool printDebug = false;

  //Define file and tree
  TFile * finMC = TFile::Open("HiForestAOD_ppMC_Py8_1.root");
  TFile * finD = TFile::Open("HiForestAOD_ppData_Jet80_1.root");


  TFile * new_File = new TFile("test_pp_jetVariable_plots.root","RECREATE");

  TTree * jetMC = (TTree*)finMC->Get(Form("ak%d%sJetAnalyzer/t",radius, algo.c_str()));
  TTree * jetD = (TTree*)finD->Get(Form("ak%d%sJetAnalyzer/t",radius, algo.c_str()));

 
  //Part 1: Plotting Manually
  // updates from github

  /*
  //Define variables
  TH1F * hpT = new TH1F("hpT","",100,0,1000);
  TH1F * heta = new TH1F("heta","",60,-3,+3);
  TH1F * hphi = new TH1F("hphi","",60,-3.15,+3.15);
  TH2F * h2d = new TH2F("h2d","pT vs eta for pT>20",30,-2,+2,100,0,500);

  //Create canvas
  TCanvas * c1 = new TCanvas("c1","Simple Jet Variables 1",1200,1200);
  c1->Divide(2,2);

  //Plotting
  c1->cd(1);
  c1->cd(1)->SetLogy();
  jet->Draw("jtpt>>hpT");
  hpT->SetTitle("pp ak3PF jet pT");
  hpT->SetXTitle("Jet p_{T} (GeV/c)");
  hpT->SetYTitle("Counts");
  hpT->Draw();

  c1->cd(2);
  jet->Draw("jteta>>heta");
  heta->SetTitle("pp ak3PF jet eta");
  heta->SetXTitle("Jet Eta");
  heta->SetYTitle("Counts");
  heta->Draw();
  
  c1->cd(3);
  jet->Draw("jtphi>>hphi");
  hphi->SetTitle("pp ak3PF jet phi");
  hphi->SetXTitle("Jet Phi (rad)");
  hphi->SetYTitle("Counts");
  hphi->Draw();
  
  c1->cd(4);
  jet->Draw("jtpt:jteta>>h2d","jtpt>20","goff");
  h2d->SetTitle("pp ak3PF jet jteta v. jtpt");
  h2d->SetXTitle("jet Eta");
  h2d->SetYTitle("jet pT (GeV/c)");
  h2d->Draw("colz");

  c1->SaveAs("test_pp_jetvariables.pdf","RECREATE");
  */
  
  //Part 2: Drawing Using Loops

  //Definitions
  TH1F * hJet_pTMC = new TH1F("hJet_pTMC", "", 100,0,1000);
  TH1F * hJet_etaMC = new TH1F("hJet_etaMC", "", 60,-3,+3);
  TH1F * hJet_phiMC = new TH1F("hJet_phiMC", "", 60,-3.15,+3.15);
  TH2F * hJet_2dMC = new TH2F("hJet_2dMC","Jet pT vs Eta for pT>20",30,-3,+3,100,0,500);
   
  Float_t jtptMC[1000];
  Float_t jtetaMC[1000];
  Float_t jtphiMC[1000];
  Int_t nrefMC;
  
  //Set branches of the tree 
  jetMC->SetBranchAddress("jtpt", &jtptMC);
  jetMC->SetBranchAddress("jteta", &jtetaMC);
  jetMC->SetBranchAddress("jtphi", &jtphiMC);
  jetMC->SetBranchAddress("nref", &nrefMC);

  Long64_t nentriesMC =  jetMC->GetEntries();
  cout<< "Number of events "<<jetMC->GetEntries()<<endl;

  if(printDebug) nentriesMC = 10;
  
  //Start MC file entry loop
  for(Long64_t nentry = 0; nentry < nentriesMC; ++nentry){
    
    jetMC->GetEvent(nentry);

    for(int jentry = 0; jentry < nrefMC; ++jentry){
      
	if(printDebug)
	  cout << "jentry = " << jentry << endl;

	hJet_pTMC->Fill(jtptMC[jentry]);
	if(printDebug)
	  cout<<"jtptMC = "<< jtptMC[jentry] << endl;

	hJet_etaMC->Fill(jtetaMC[jentry]);
	if(printDebug)
	  cout<<"jtetaMC = "<< jtetaMC[jentry] << endl;

	hJet_phiMC->Fill(jtphiMC[jentry]);
	if(printDebug)
	  cout<<"jtphiMC = "<< jtphiMC[jentry] << endl;

	hJet_2dMC->Fill(jtetaMC[jentry],jtptMC[jentry]);

	// }//! jet pt and eta selection criteria 

    }//end MC jet loop
   
  }//end MC entry loop


  // Loops for the Data File

  TH1F * hJet_pTD = new TH1F("hJet_pTD", "", 100,0,1000);
  TH1F * hJet_etaD = new TH1F("hJet_etaD", "", 60,-3,+3);
  TH1F * hJet_phiD = new TH1F("hJet_phiD", "", 60,-3.15,+3.15);
  TH2F * hJet_2dD = new TH2F("hJet_2dD","Jet pT vs Eta for pT>20",30,-3,+3,100,0,500);
  
  Float_t jtptD[1000];
  Float_t jtetaD[1000];
  Float_t jtphiD[1000];
  Int_t nrefD;
  
  //Set branches of the tree 
  jetD->SetBranchAddress("jtpt", &jtptD);
  jetD->SetBranchAddress("jteta", &jtetaD);
  jetD->SetBranchAddress("jtphi", &jtphiD);
  jetD->SetBranchAddress("nref", &nrefD);

  Long64_t nentriesD =  jetD->GetEntries();
  cout<< "Number of events "<<jetD->GetEntries()<<endl;

  if(printDebug) nentriesD = 10;
  
  //Start data file entry loop
  for(Long64_t nentry = 0; nentry < nentriesD; ++nentry){
    
    jetD->GetEvent(nentry);
    // if(printDebug)
    // cout << nentry << "/" << nentries <<endl;

    /*

      nref = 4;
      jtpt[] = {136.24, 104.89, 18.04, 10.22};
      
     */
    
    //start data file jet loop
    for(int jentry = 0; jentry < nrefD; ++jentry){

      //if(jtpt[jentry] <= 100 || jteta[jentry] <= 1.0)
      //continue;
      
	// if(jtpt[jentry] > 100 && jteta[jentry] > 1.0) {
      
	if(printDebug)
	  cout << "jentry = " << jentry << endl;

	hJet_pTD->Fill(jtptD[jentry]);
	if(printDebug)
	  cout<<"jtptD = "<< jtptD[jentry] << endl;

	hJet_etaD->Fill(jtetaD[jentry]);
	if(printDebug)
	  cout<<"jtetaD = "<< jtetaD[jentry] << endl;

	hJet_phiD->Fill(jtphiD[jentry]);
	if(printDebug)
	  cout<<"jtphiD = "<< jtphiD[jentry] << endl;

	hJet_2dD->Fill(jtetaD[jentry],jtptD[jentry]);

	// }//! jet pt and eta selection criteria 

    }//end data jet loop
   
  }//end data entry loop

  // check histograms (Prints histogram information)
  hJet_2dMC->Print("base");
  cout<<"histogram mean = "<<hJet_pTMC->GetMean()<<endl;

  hJet_2dD->Print("base");
  cout<<"histogram mean = "<<hJet_pTD->GetMean()<<endl;
  
  // plotting
  TCanvas * c2 = new TCanvas("c2", "Simple Jet Variables 2", 1200,1200);
  c2->Divide(2,2);
  
  c2->cd(1)->SetLogy();
  hJet_pTMC->SetTitle("pp ak3PF Jet pT");
  hJet_pTMC->SetXTitle("Jet pT (GeV/c)");
  hJet_pTMC->SetYTitle("Counts");
  hJet_pTMC->SetLineColorAlpha(kBlue,0.9);
  hJet_pTMC->Draw();
  hJet_pTD->SetLineColorAlpha(kGreen, 0.9);
  hJet_pTD->Draw("same");
  TLegend * ltest1 = new TLegend(0.7,0.75,0.8,0.85);
  ltest1->AddEntry(hJet_pTMC,"MC","line");
  ltest1->AddEntry(hJet_pTD,"Data","line");
  ltest1->SetTextSize(0.04);
  ltest1->SetBorderSize(0);
  ltest1->SetFillStyle(0);
  ltest1->Draw();

  c2->cd(2);
  c2->cd(2)->SetLogy();
  hJet_etaMC->SetTitle("pp ak3PF Jet Eta");
  hJet_etaMC->SetXTitle("Jet Eta");
  hJet_etaMC->SetYTitle("Counts");
  hJet_etaMC->SetLineColorAlpha(kBlue,0.9);
  hJet_etaMC->Scale(1./nentriesMC);
  hJet_etaD->Scale(1./nentriesD);
  hJet_etaMC->Draw();
  hJet_etaD->SetLineColorAlpha(kGreen, 0.9);
  hJet_etaD->Draw("same");
  TLegend * ltest2 = new TLegend(0.7,0.75,0.8,0.85);
  ltest2->AddEntry(hJet_etaMC,"MC","line");
  ltest2->AddEntry(hJet_etaD,"Data","line");
  ltest2->SetTextSize(0.04);
  ltest2->SetBorderSize(0);
  ltest2->SetFillStyle(0);
  ltest2->Draw();
    
  c2->cd(3);
  c2->cd(3)->SetLogy();
  hJet_phiMC->SetTitle("pp ak3PF Jet Phi");
  hJet_phiMC->SetXTitle("Jet Phi (Rad)");
  hJet_phiMC->SetYTitle("Counts");
  hJet_phiMC->SetLineColorAlpha(kBlue,0.9);
  hJet_phiMC->Scale(1./nentriesMC);
  hJet_phiD->Scale(1./nentriesD);
  hJet_phiMC->Draw();
  hJet_phiD->SetLineColorAlpha(kGreen, 0.9);
  hJet_phiD->Draw("same");
  TLegend * ltest3 = new TLegend(0.7,0.75,0.8,0.85);
  ltest3->AddEntry(hJet_phiMC,"MC","line");
  ltest3->AddEntry(hJet_phiD,"Data","line");
  ltest3->SetTextSize(0.04);
  ltest3->SetBorderSize(0);
  ltest3->SetFillStyle(0);
  ltest3->Draw();

  c2->cd(4);
  c2->cd(4)->SetLogz();
  hJet_2dMC->SetTitle("pp ak3PF Jet Eta vs. pT for pT>20");
  hJet_2dMC->SetXTitle("Jet Eta");
  hJet_2dMC->SetYTitle("Jet pT (GeV/c)");
  hJet_phiMC->SetFillColor(kBlue);
  hJet_2dD->SetFillColor(kGreen);
  TLegend * ltest4 = new TLegend(0.7,0.75,0.8,0.85);
  ltest4->AddEntry(hJet_2dMC,"MC","");
  ltest4->AddEntry(hJet_2dD,"Data","");
  ltest4->SetTextSize(0.04);
  ltest4->SetBorderSize(0);
  ltest4->SetFillStyle(0);
  ltest4->Draw();
  THStack * theStack = new THStack("theStack","pp ak3PF Jet Eta vs. pT");
  theStack->Add(hJet_2dMC);
  theStack->Add(hJet_2dD);
  theStack->Draw();



  //creates a pdf
  new_File->Write();
  c2->SaveAs("test_pp_jetvariables_with_fors.pdf","RECREATE");
 
}// end of macro
