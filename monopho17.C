//05.11.2020, De-Lin Macive Xiong, Florida State University
//monophoton non-collision background estimate 2017 dataset


#define monopho17_cxx
#include "monopho17.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>
#include <ctime>
#include <chrono>
#include <fstream>

using namespace std;


Float_t MIPChi2_over_MIPNhitCone(Float_t MIPChi2, Float_t MIPNhitCone)
{
	Float_t phoMIPChi2 = MIPChi2;
	Float_t phoMIPNhitCone = MIPNhitCone;
	Float_t MIPChi2_over_MIPNhitCone = MIPChi2 / MIPNhitCone;

	return MIPChi2_over_MIPNhitCone;
}

Float_t newmMET(Float_t met, Float_t metphi, Float_t phi, Float_t et)
{
	Float_t phoEt;
	Float_t phoPhi;
	Float_t pfMET = met;
	Float_t pfMETPhi = metphi;
	Float_t ophoEt = et;
	Float_t ophoPhi = phi;

	Float_t pfMETX = met * cos(pfMETPhi);
	Float_t pfMETY = met * sin(pfMETPhi);
	Float_t ITEtX = phoEt * cos(phoPhi);
	Float_t ITEtY = phoEt * sin(phoPhi);

	pfMETX += ITEtX;
	pfMETY += ITEtY;

	Float_t OOTEtX = et * cos(ophoPhi);
	Float_t OOTEtY = et * sin(ophoPhi);

	pfMETX -= OOTEtX;
	pfMETY -= OOTEtY;

	Float_t newmMET = sqrt(pfMETX * pfMETX + pfMETY * pfMETY);

	return newmMET;
}

Float_t newuMET(Float_t met, Float_t metphi, Float_t phi, Float_t et)
{
	Float_t ophoEt = et;
	Float_t ophoPhi = phi;
	Float_t pfMET = met;
	Float_t pfMETPhi = metphi;

	Float_t pfMETX = met * cos(pfMETPhi);
	Float_t pfMETY = met * sin(pfMETPhi);
	Float_t OOTEtX = et * cos(ophoPhi);
	Float_t OOTEtY = et * sin(ophoPhi);

	pfMETX -= OOTEtX;
	pfMETY -= OOTEtY;

	Float_t newuMET = sqrt(pfMETX * pfMETX + pfMETY * pfMETY);

	return newuMET;
}



//calulate invariant mass Z

Float_t InvariMass(Float_t Et1, Float_t Et2, Float_t Phi1, Float_t Phi2, Float_t Eta1, Float_t Eta2)
{
	Float_t Theta1 = 2 * atan(exp(-1.*Eta1));
	Float_t Theta2 = 2 * atan(exp(-1.*Eta2));
	Float_t phoPhi1 = Phi1;
	Float_t phoPhi2 = Phi2;
	Float_t Etot1 = Et1 / sin(Theta1); //E_tot1
	Float_t Etot2 = Et2 / sin(Theta2); //E_tot2
	
	//reconstruct the vectors for x, y and z

	
	Float_t phoX1 = Etot1 * cos(Phi1) * sin(Theta1); 
	Float_t phoY1 = Etot1 * sin(Phi1) * sin(Theta1);
	Float_t phoZ1 = Etot1 * cos(Theta1);

	Float_t phoX2 = Etot2 * cos(Phi2) * sin(Theta2);
	Float_t phoY2 = Etot2 * sin(Phi2) * sin(Theta2);
	Float_t phoZ2 = Etot2 * cos(Theta2);

	Float_t EX1 = Et1*sin(Phi1);
	Float_t EY1 = Et1*cos(Phi1);
	Float_t EZ1 = Etot1*cos(Theta1);

	Float_t EX2 = Et2*sin(Phi2);
	Float_t EY2 = Et2*cos(Phi2);
	Float_t EZ2 = Etot2*cos(Theta2);

	Float_t E1 = sqrt(Etot1*Etot1);
	Float_t E2 = sqrt(Etot2*Etot2);

	//Float_t InvariMassSq = 2 * E1*E2 - phoMag1 * phoMag1 - phoMag2 * phoMag2 - 2 * DotProd12;

	//1	
	Float_t InvariMassSq = (Etot1 + Etot2)*(Etot1 + Etot2) - (EX1 + EX2)*(EX1 + EX2) - (EY1 + EY2)*(EY1 + EY2) - (EZ1 + EZ2)*(EZ1 + EZ2);

	Float_t InvMass = sqrt(InvariMassSq);

	
	
	return InvMass;
	
}


//calulate transverse mass W (we cannot calcuate the invariant mass of W, since we don't know the mass of neutrino)

Float_t TransMassW(Float_t Et1, Float_t Phi1, Float_t Eta1)
{
	Float_t Theta1 = 2 * atan(exp(-1.*Eta1));
	Float_t phoPhi1 = Phi1;
	Float_t Etot1 = Et1 / sin(Theta1); //E_tot1

	//reconstruct the vectors for x, y and z

	Float_t phoX1 = Etot1 * cos(Phi1) * sin(Theta1);
	Float_t phoY1 = Etot1 * sin(Phi1) * sin(Theta1);
	Float_t phoZ1 = Etot1 * cos(Theta1);

	Float_t EX1 = Et1 * sin(Phi1);
	Float_t EY1 = Et1 * cos(Phi1);
	Float_t EZ1 = Etot1 * cos(Theta1);

	Float_t TransMassSqW = (Etot1)*(Etot1) - (EZ1)*(EZ1);

	Float_t TranMassW = sqrt(TransMassSqW);



	return TranMassW;

}


//for checking matching

Float_t DeltaR(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2)
{

	Float_t ophoPhi = phi1;
	Float_t phoPhi = phi2;
	Float_t dphi = fabs(phoPhi - ophoPhi);
	Float_t tdphi = dphi;
    //phi wrap
	if (dphi > TMath::Pi()) tdphi = TMath::Pi()*2.0 - dphi;
	dphi = tdphi;

	Float_t ophoEta = eta1;
	Float_t phoEta = eta2;
	Float_t deta = fabs(phoEta - ophoEta);

	Float_t deltaR = sqrt(deta*deta + dphi * dphi);
	return deltaR;
}


void monopho17::Loop()
{
//   In a ROOT session, you can do:
//      root> .L monopho17.C
//      root> monopho17 t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
    TFile *monophoseedtime = new TFile("2017Dataset_eta-wing_spikeincand.root", "RECREATE");
	

    //----------------------------------detid---------------------------------
    //checking info
    TH1F *ophoseedieta = new TH1F("ophoseedieta", "ophoseedieta", 200, 100.0, 100.0);
	ophoseedieta->GetXaxis()->SetTitle("ieta");
	ophoseedieta->GetYaxis()->SetTitle("Entries");

    TH1F *ophoseediphi = new TH1F("ophoseediphi", "ophoseediphi", 200, 100.0, 100.0);
	ophoseediphi->GetXaxis()->SetTitle("iphi");
	ophoseediphi->GetYaxis()->SetTitle("Entries");

    
	TH1F *phoseedieta = new TH1F("phoseedieta", "phoseedieta", 2000, 100.0, 100.0);
	phoseedieta->GetXaxis()->SetTitle("ieta");
	phoseedieta->GetYaxis()->SetTitle("Entries");

    TH1F *phoseediphi = new TH1F("phoseediphi", "phoseediphi", 2000, 100.0, 100.0);
	phoseediphi->GetXaxis()->SetTitle("iphi");
	phoseediphi->GetYaxis()->SetTitle("Entries");


    //----------------------------------candidate template----------------------------------
    TH1F *candidate_events_case1 = new TH1F("candidate_events_case1", "candidate_events_case1", 200, -25, 25);
	candidate_events_case1->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_case1->GetYaxis()->SetTitle("Entries");	  

	TH1F *candidate_events_case2 = new TH1F("candidate_events_case2", "candidate_events_case2", 200, -25, 25);
	candidate_events_case2->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_case2->GetYaxis()->SetTitle("Entries");

	TH1F *candidate_events_case3 = new TH1F("candidate_events_case3", "candidate_events_case3", 200, -25, 25);
	candidate_events_case3->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_case3->GetYaxis()->SetTitle("Entries");

	TH1F *candidate_events_case4 = new TH1F("candidate_events_case4", "candidate_events_case4", 200, -25, 25);
	candidate_events_case4->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_case4->GetYaxis()->SetTitle("Entries");

	TH1F *candidate_events_case5 = new TH1F("candidate_events_case5", "candidate_events_case5", 200, -25, 25);
	candidate_events_case5->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_case5->GetYaxis()->SetTitle("Entries");

	TH1F *candidate_events = new TH1F("candidate_events", "candidate_events", 200, -25, 25);
	candidate_events->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events->GetYaxis()->SetTitle("Entries");

    TH1F *candidate_events_3ns = new TH1F("candidate_events_3ns", "candidate_events_3ns", 200, -25, 25);
	candidate_events_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_3ns->GetYaxis()->SetTitle("Entries");



    //candidate etawing calculation
    //if there's only one
    TH1F *candidate_eta_wing_rightonly1 = new TH1F("candidate_eta_wing_rightonly1", "candidate_eta_wing_rightonly1", 2000, 0, 1);
	candidate_eta_wing_rightonly1->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing_rightonly1->GetYaxis()->SetTitle("Entries");

    TH1F *candidate_eta_wing_rightonly2 = new TH1F("candidate_eta_wing_rightonly2", "candidate_eta_wing_rightonly2", 2000, 0, 1);
	candidate_eta_wing_rightonly2->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing_rightonly2->GetYaxis()->SetTitle("Entries");

    TH1F *candidate_eta_wing_rightonly3 = new TH1F("candidate_eta_wing_rightonly3", "candidate_eta_wing_rightonly3", 2000, 0, 1);
	candidate_eta_wing_rightonly3->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing_rightonly3->GetYaxis()->SetTitle("Entries");

    TH1F *candidate_eta_wing_rightonly4 = new TH1F("candidate_eta_wing_rightonly4", "candidate_eta_wing_rightonly4", 2000, 0, 1);
	candidate_eta_wing_rightonly4->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing_rightonly4->GetYaxis()->SetTitle("Entries");

    TH1F *candidate_eta_wing_rightonly5 = new TH1F("candidate_eta_wing_rightonly5", "candidate_eta_wing_rightonly5", 2000, 0, 1);
	candidate_eta_wing_rightonly5->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing_rightonly5->GetYaxis()->SetTitle("Entries");


    TH1F *candidate_eta_wing_leftonly1 = new TH1F("candidate_eta_wing_leftonly1", "candidate_eta_wing_leftonly1", 2000, 0, 1);
	candidate_eta_wing_leftonly1->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing_leftonly1->GetYaxis()->SetTitle("Entries");

    TH1F *candidate_eta_wing_leftonly2 = new TH1F("candidate_eta_wing_leftonly2", "candidate_eta_wing_leftonly2", 2000, 0, 1);
	candidate_eta_wing_leftonly2->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing_leftonly2->GetYaxis()->SetTitle("Entries");

    TH1F *candidate_eta_wing_leftonly3 = new TH1F("candidate_eta_wing_leftonly3", "candidate_eta_wing_leftonly3", 2000, 0, 1);
	candidate_eta_wing_leftonly3->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing_leftonly3->GetYaxis()->SetTitle("Entries");

    TH1F *candidate_eta_wing_leftonly4 = new TH1F("candidate_eta_wing_leftonly4", "candidate_eta_wing_leftonly4", 2000, 0, 1);
	candidate_eta_wing_leftonly4->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing_leftonly4->GetYaxis()->SetTitle("Entries");

    TH1F *candidate_eta_wing_leftonly5 = new TH1F("candidate_eta_wing_leftonly5", "candidate_eta_wing_leftonly5", 2000, 0, 1);
	candidate_eta_wing_leftonly5->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing_leftonly5->GetYaxis()->SetTitle("Entries");

    //if there are both
    TH1F *candidate_eta_wing_both1 = new TH1F("candidate_eta_wing_both1", "candidate_eta_wing_both1", 2000, 0, 1);
	candidate_eta_wing_both1->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing_both1->GetYaxis()->SetTitle("Entries");

    TH1F *candidate_eta_wing_both2 = new TH1F("candidate_eta_wing_both2", "candidate_eta_wing_both2", 2000, 0, 1);
	candidate_eta_wing_both2->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing_both2->GetYaxis()->SetTitle("Entries");

    TH1F *candidate_eta_wing_both3 = new TH1F("candidate_eta_wing_both3", "candidate_eta_wing_both3", 2000, 0, 1);
	candidate_eta_wing_both3->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing_both3->GetYaxis()->SetTitle("Entries");

    TH1F *candidate_eta_wing_both4 = new TH1F("candidate_eta_wing_both4", "candidate_eta_wing_both4", 2000, 0, 1);
	candidate_eta_wing_both4->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing_both4->GetYaxis()->SetTitle("Entries");

    TH1F *candidate_eta_wing_both5 = new TH1F("candidate_eta_wing_both5", "candidate_eta_wing_both5", 2000, 0, 1);
	candidate_eta_wing_both5->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing_both5->GetYaxis()->SetTitle("Entries");


    //combine all of the etawing
    TH1F *candidate_eta_wing1 = new TH1F("candidate_eta_wing1", "candidate_eta_wing1", 2000, 0, 1);
	candidate_eta_wing1->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing1->GetYaxis()->SetTitle("Entries");

    TH1F *candidate_eta_wing2 = new TH1F("candidate_eta_wing2", "candidate_eta_wing2", 2000, 0, 1);
	candidate_eta_wing2->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing2->GetYaxis()->SetTitle("Entries");

    TH1F *candidate_eta_wing3 = new TH1F("candidate_eta_wing3", "candidate_eta_wing3", 2000, 0, 1);
	candidate_eta_wing3->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing3->GetYaxis()->SetTitle("Entries");

    TH1F *candidate_eta_wing4 = new TH1F("candidate_eta_wing4", "candidate_eta_wing4", 2000, 0, 1);
	candidate_eta_wing4->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing4->GetYaxis()->SetTitle("Entries");

    TH1F *candidate_eta_wing5 = new TH1F("candidate_eta_wing5", "candidate_eta_wing5", 2000, 0, 1);
	candidate_eta_wing5->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing5->GetYaxis()->SetTitle("Entries");

    TH1F *candidate_eta_wing = new TH1F("candidate_eta_wing", "candidate_eta_wing", 2000, 0, 1);
	candidate_eta_wing->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing->GetYaxis()->SetTitle("Entries");

    TH1F *candidate_eta_wing_prompt = new TH1F("candidate_eta_wing_prompt", "candidate_eta_wing_prompt", 2000, 0, 1);
	candidate_eta_wing_prompt->GetXaxis()->SetTitle("E_wing");
	candidate_eta_wing_prompt->GetYaxis()->SetTitle("Entries");



    //----------------------------------spikeincandidate template----------------------------------
    TH1F *spikeincandidate_events_case1 = new TH1F("spikeincandidate_events_case1", "spikeincandidate_events_case1", 200, -25, 25);
	spikeincandidate_events_case1->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_case1->GetYaxis()->SetTitle("Entries");	  

	TH1F *spikeincandidate_events_case2 = new TH1F("spikeincandidate_events_case2", "spikeincandidate_events_case2", 200, -25, 25);
	spikeincandidate_events_case2->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_case2->GetYaxis()->SetTitle("Entries");

	TH1F *spikeincandidate_events_case3 = new TH1F("spikeincandidate_events_case3", "spikeincandidate_events_case3", 200, -25, 25);
	spikeincandidate_events_case3->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_case3->GetYaxis()->SetTitle("Entries");

	TH1F *spikeincandidate_events_case4 = new TH1F("spikeincandidate_events_case4", "spikeincandidate_events_case4", 200, -25, 25);
	spikeincandidate_events_case4->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_case4->GetYaxis()->SetTitle("Entries");

	TH1F *spikeincandidate_events_case5 = new TH1F("spikeincandidate_events_case5", "spikeincandidate_events_case5", 200, -25, 25);
	spikeincandidate_events_case5->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_case5->GetYaxis()->SetTitle("Entries");

	TH1F *spikeincandidate_events = new TH1F("spikeincandidate_events", "spikeincandidate_events", 200, -25, 25);
	spikeincandidate_events->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events->GetYaxis()->SetTitle("Entries");

    TH1F *spikeincandidate_events_3ns = new TH1F("spikeincandidate_events_3ns", "spikeincandidate_events_3ns", 200, -25, 25);
	spikeincandidate_events_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_3ns->GetYaxis()->SetTitle("Entries");


    
    //spikeincandidate etawing calculation
    //if there's only one
    TH1F *spikeincandidate_eta_wing_rightonly1 = new TH1F("spikeincandidate_eta_wing_rightonly1", "spikeincandidate_eta_wing_rightonly1", 2000, 0, 1);
	spikeincandidate_eta_wing_rightonly1->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing_rightonly1->GetYaxis()->SetTitle("Entries");

    TH1F *spikeincandidate_eta_wing_rightonly2 = new TH1F("spikeincandidate_eta_wing_rightonly2", "spikeincandidate_eta_wing_rightonly2", 2000, 0, 1);
	spikeincandidate_eta_wing_rightonly2->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing_rightonly2->GetYaxis()->SetTitle("Entries");

    TH1F *spikeincandidate_eta_wing_rightonly3 = new TH1F("spikeincandidate_eta_wing_rightonly3", "spikeincandidate_eta_wing_rightonly3", 2000, 0, 1);
	spikeincandidate_eta_wing_rightonly3->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing_rightonly3->GetYaxis()->SetTitle("Entries");

    TH1F *spikeincandidate_eta_wing_rightonly4 = new TH1F("spikeincandidate_eta_wing_rightonly4", "spikeincandidate_eta_wing_rightonly4", 2000, 0, 1);
	spikeincandidate_eta_wing_rightonly4->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing_rightonly4->GetYaxis()->SetTitle("Entries");

    TH1F *spikeincandidate_eta_wing_rightonly5 = new TH1F("spikeincandidate_eta_wing_rightonly5", "spikeincandidate_eta_wing_rightonly5", 2000, 0, 1);
	spikeincandidate_eta_wing_rightonly5->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing_rightonly5->GetYaxis()->SetTitle("Entries");


    TH1F *spikeincandidate_eta_wing_leftonly1 = new TH1F("spikeincandidate_eta_wing_leftonly1", "spikeincandidate_eta_wing_leftonly1", 2000, 0, 1);
	spikeincandidate_eta_wing_leftonly1->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing_leftonly1->GetYaxis()->SetTitle("Entries");

    TH1F *spikeincandidate_eta_wing_leftonly2 = new TH1F("spikeincandidate_eta_wing_leftonly2", "spikeincandidate_eta_wing_leftonly2", 2000, 0, 1);
	spikeincandidate_eta_wing_leftonly2->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing_leftonly2->GetYaxis()->SetTitle("Entries");

    TH1F *spikeincandidate_eta_wing_leftonly3 = new TH1F("spikeincandidate_eta_wing_leftonly3", "spikeincandidate_eta_wing_leftonly3", 2000, 0, 1);
	spikeincandidate_eta_wing_leftonly3->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing_leftonly3->GetYaxis()->SetTitle("Entries");

    TH1F *spikeincandidate_eta_wing_leftonly4 = new TH1F("spikeincandidate_eta_wing_leftonly4", "spikeincandidate_eta_wing_leftonly4", 2000, 0, 1);
	spikeincandidate_eta_wing_leftonly4->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing_leftonly4->GetYaxis()->SetTitle("Entries");

    TH1F *spikeincandidate_eta_wing_leftonly5 = new TH1F("spikeincandidate_eta_wing_leftonly5", "spikeincandidate_eta_wing_leftonly5", 2000, 0, 1);
	spikeincandidate_eta_wing_leftonly5->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing_leftonly5->GetYaxis()->SetTitle("Entries");

    //if there are both
    TH1F *spikeincandidate_eta_wing_both1 = new TH1F("spikeincandidate_eta_wing_both1", "spikeincandidate_eta_wing_both1", 2000, 0, 1);
	spikeincandidate_eta_wing_both1->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing_both1->GetYaxis()->SetTitle("Entries");

    TH1F *spikeincandidate_eta_wing_both2 = new TH1F("spikeincandidate_eta_wing_both2", "spikeincandidate_eta_wing_both2", 2000, 0, 1);
	spikeincandidate_eta_wing_both2->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing_both2->GetYaxis()->SetTitle("Entries");

    TH1F *spikeincandidate_eta_wing_both3 = new TH1F("spikeincandidate_eta_wing_both3", "spikeincandidate_eta_wing_both3", 2000, 0, 1);
	spikeincandidate_eta_wing_both3->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing_both3->GetYaxis()->SetTitle("Entries");

    TH1F *spikeincandidate_eta_wing_both4 = new TH1F("spikeincandidate_eta_wing_both4", "spikeincandidate_eta_wing_both4", 2000, 0, 1);
	spikeincandidate_eta_wing_both4->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing_both4->GetYaxis()->SetTitle("Entries");

    TH1F *spikeincandidate_eta_wing_both5 = new TH1F("spikeincandidate_eta_wing_both5", "spikeincandidate_eta_wing_both5", 2000, 0, 1);
	spikeincandidate_eta_wing_both5->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing_both5->GetYaxis()->SetTitle("Entries");


    //combine all of the etawing
    TH1F *spikeincandidate_eta_wing1 = new TH1F("spikeincandidate_eta_wing1", "spikeincandidate_eta_wing1", 2000, 0, 1);
	spikeincandidate_eta_wing1->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing1->GetYaxis()->SetTitle("Entries");

    TH1F *spikeincandidate_eta_wing2 = new TH1F("spikeincandidate_eta_wing2", "spikeincandidate_eta_wing2", 2000, 0, 1);
	spikeincandidate_eta_wing2->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing2->GetYaxis()->SetTitle("Entries");

    TH1F *spikeincandidate_eta_wing3 = new TH1F("spikeincandidate_eta_wing3", "spikeincandidate_eta_wing3", 2000, 0, 1);
	spikeincandidate_eta_wing3->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing3->GetYaxis()->SetTitle("Entries");

    TH1F *spikeincandidate_eta_wing4 = new TH1F("spikeincandidate_eta_wing4", "spikeincandidate_eta_wing4", 2000, 0, 1);
	spikeincandidate_eta_wing4->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing4->GetYaxis()->SetTitle("Entries");

    TH1F *spikeincandidate_eta_wing5 = new TH1F("spikeincandidate_eta_wing5", "spikeincandidate_eta_wing5", 2000, 0, 1);
	spikeincandidate_eta_wing5->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing5->GetYaxis()->SetTitle("Entries");

    TH1F *spikeincandidate_eta_wing = new TH1F("spikeincandidate_eta_wing", "spikeincandidate_eta_wing", 2000, 0, 1);
	spikeincandidate_eta_wing->GetXaxis()->SetTitle("E_wing");
	spikeincandidate_eta_wing->GetYaxis()->SetTitle("Entries");



    //----------------------------------promptZ template----------------------------------

    TH1F *promptZ_events_case1 = new TH1F("promptZ_events_case1", "promptZ_events_case1", 200, -25, 25);
	promptZ_events_case1->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_case1->GetYaxis()->SetTitle("Entries");	  

	TH1F *promptZ_events_case2 = new TH1F("promptZ_events_case2", "promptZ_events_case2", 200, -25, 25);
	promptZ_events_case2->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_case2->GetYaxis()->SetTitle("Entries");

	TH1F *promptZ_events_case3 = new TH1F("promptZ_events_case3", "promptZ_events_case3", 200, -25, 25);
	promptZ_events_case3->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_case3->GetYaxis()->SetTitle("Entries");

	TH1F *promptZ_events_case4 = new TH1F("promptZ_events_case4", "promptZ_events_case4", 200, -25, 25);
	promptZ_events_case4->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_case4->GetYaxis()->SetTitle("Entries");

	TH1F *promptZ_events_case5 = new TH1F("promptZ_events_case5", "promptZ_events_case5", 200, -25, 25);
	promptZ_events_case5->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_case5->GetYaxis()->SetTitle("Entries");

	TH1F *promptZ_events = new TH1F("promptZ_events", "promptZ_events", 200, -25, 25);
	promptZ_events->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_events_3ns = new TH1F("promptZ_events_3ns", "promptZ_events_3ns", 200, -25, 25);
	promptZ_events_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_3ns->GetYaxis()->SetTitle("Entries");



    TH1F *InvMass_Z = new TH1F("InvMass_Z", "InvMass_Z", 2000, 25, 25);
	InvMass_Z->GetXaxis()->SetTitle("Invaraint mass of Z (GeV)");
	InvMass_Z->GetYaxis()->SetTitle("Entries");

    //promptZ etawing calculation
    //if there's only one
    TH1F *promptZ_eta_wing_rightonly1 = new TH1F("promptZ_eta_wing_rightonly1", "promptZ_eta_wing_rightonly1", 2000, 0, 1);
	promptZ_eta_wing_rightonly1->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing_rightonly1->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_eta_wing_rightonly2 = new TH1F("promptZ_eta_wing_rightonly2", "promptZ_eta_wing_rightonly2", 2000, 0, 1);
	promptZ_eta_wing_rightonly2->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing_rightonly2->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_eta_wing_rightonly3 = new TH1F("promptZ_eta_wing_rightonly3", "promptZ_eta_wing_rightonly3", 2000, 0, 1);
	promptZ_eta_wing_rightonly3->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing_rightonly3->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_eta_wing_rightonly4 = new TH1F("promptZ_eta_wing_rightonly4", "promptZ_eta_wing_rightonly4", 2000, 0, 1);
	promptZ_eta_wing_rightonly4->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing_rightonly4->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_eta_wing_rightonly5 = new TH1F("promptZ_eta_wing_rightonly5", "promptZ_eta_wing_rightonly5", 2000, 0, 1);
	promptZ_eta_wing_rightonly5->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing_rightonly5->GetYaxis()->SetTitle("Entries");


    TH1F *promptZ_eta_wing_leftonly1 = new TH1F("promptZ_eta_wing_leftonly1", "promptZ_eta_wing_leftonly1", 2000, 0, 1);
	promptZ_eta_wing_leftonly1->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing_leftonly1->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_eta_wing_leftonly2 = new TH1F("promptZ_eta_wing_leftonly2", "promptZ_eta_wing_leftonly2", 2000, 0, 1);
	promptZ_eta_wing_leftonly2->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing_leftonly2->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_eta_wing_leftonly3 = new TH1F("promptZ_eta_wing_leftonly3", "promptZ_eta_wing_leftonly3", 2000, 0, 1);
	promptZ_eta_wing_leftonly3->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing_leftonly3->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_eta_wing_leftonly4 = new TH1F("promptZ_eta_wing_leftonly4", "promptZ_eta_wing_leftonly4", 2000, 0, 1);
	promptZ_eta_wing_leftonly4->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing_leftonly4->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_eta_wing_leftonly5 = new TH1F("promptZ_eta_wing_leftonly5", "promptZ_eta_wing_leftonly5", 2000, 0, 1);
	promptZ_eta_wing_leftonly5->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing_leftonly5->GetYaxis()->SetTitle("Entries");

    //if there are both
    TH1F *promptZ_eta_wing_both1 = new TH1F("promptZ_eta_wing_both1", "promptZ_eta_wing_both1", 2000, 0, 1);
	promptZ_eta_wing_both1->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing_both1->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_eta_wing_both2 = new TH1F("promptZ_eta_wing_both2", "promptZ_eta_wing_both2", 2000, 0, 1);
	promptZ_eta_wing_both2->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing_both2->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_eta_wing_both3 = new TH1F("promptZ_eta_wing_both3", "promptZ_eta_wing_both3", 2000, 0, 1);
	promptZ_eta_wing_both3->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing_both3->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_eta_wing_both4 = new TH1F("promptZ_eta_wing_both4", "promptZ_eta_wing_both4", 2000, 0, 1);
	promptZ_eta_wing_both4->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing_both4->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_eta_wing_both5 = new TH1F("promptZ_eta_wing_both5", "promptZ_eta_wing_both5", 2000, 0, 1);
	promptZ_eta_wing_both5->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing_both5->GetYaxis()->SetTitle("Entries");


    //combine all of the etawing
    TH1F *promptZ_eta_wing1 = new TH1F("promptZ_eta_wing1", "promptZ_eta_wing1", 2000, 0, 1);
	promptZ_eta_wing1->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing1->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_eta_wing2 = new TH1F("promptZ_eta_wing2", "promptZ_eta_wing2", 2000, 0, 1);
	promptZ_eta_wing2->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing2->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_eta_wing3 = new TH1F("promptZ_eta_wing3", "promptZ_eta_wing3", 2000, 0, 1);
	promptZ_eta_wing3->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing3->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_eta_wing4 = new TH1F("promptZ_eta_wing4", "promptZ_eta_wing4", 2000, 0, 1);
	promptZ_eta_wing4->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing4->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_eta_wing5 = new TH1F("promptZ_eta_wing5", "promptZ_eta_wing5", 2000, 0, 1);
	promptZ_eta_wing5->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing5->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_eta_wing = new TH1F("promptZ_eta_wing", "promptZ_eta_wing", 2000, 0, 1);
	promptZ_eta_wing->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing->GetYaxis()->SetTitle("Entries");

    TH1F *promptZ_eta_wing_prompt = new TH1F("promptZ_eta_wing_prompt", "promptZ_eta_wing_prompt", 2000, 0, 1);
	promptZ_eta_wing_prompt->GetXaxis()->SetTitle("E_wing");
	promptZ_eta_wing_prompt->GetYaxis()->SetTitle("Entries");
 


    //----------------------------------spike template----------------------------------
    TH1F *spike_events_case1 = new TH1F("spike_events_case1", "spike_events_case1", 200, -25, 25);
	spike_events_case1->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_case1->GetYaxis()->SetTitle("Entries");	  

	TH1F *spike_events_case2 = new TH1F("spike_events_case2", "spike_events_case2", 200, -25, 25);
	spike_events_case2->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_case2->GetYaxis()->SetTitle("Entries");

	TH1F *spike_events_case3 = new TH1F("spike_events_case3", "spike_events_case3", 200, -25, 25);
	spike_events_case3->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_case3->GetYaxis()->SetTitle("Entries");

	TH1F *spike_events_case4 = new TH1F("spike_events_case4", "spike_events_case4", 200, -25, 25);
	spike_events_case4->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_case4->GetYaxis()->SetTitle("Entries");

	TH1F *spike_events_case5 = new TH1F("spike_events_case5", "spike_events_case5", 200, -25, 25);
	spike_events_case5->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_case5->GetYaxis()->SetTitle("Entries");

	TH1F *spike_events = new TH1F("spike_events", "spike_events", 200, -25, 25);
	spike_events->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events->GetYaxis()->SetTitle("Entries");

    TH1F *spike_events_3ns = new TH1F("spike_events_3ns", "spike_events_3ns", 200, -25, 25);
	spike_events_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_3ns->GetYaxis()->SetTitle("Entries");



    //eta vs phi for spike
    TH2F *spike_events_etaphi_case1 = new TH2F("spike_events_etaphi_case1", "spike_events_etaphi_case1", 200, 20, 20, 200, 20, 20);
	spike_events_etaphi_case1->GetXaxis()->SetTitle("eta");
	spike_events_etaphi_case1->GetYaxis()->SetTitle("phi");

    TH2F *spike_events_etaphi_case2 = new TH2F("spike_events_etaphi_case2", "spike_events_etaphi_case2", 200, 20, 20, 200, 20, 20);
	spike_events_etaphi_case2->GetXaxis()->SetTitle("eta");
	spike_events_etaphi_case2->GetYaxis()->SetTitle("phi");

    TH2F *spike_events_etaphi_case3 = new TH2F("spike_events_etaphi_case3", "spike_events_etaphi_case3", 200, 20, 20, 200, 20, 20);
	spike_events_etaphi_case3->GetXaxis()->SetTitle("eta");
	spike_events_etaphi_case3->GetYaxis()->SetTitle("phi");

    TH2F *spike_events_etaphi_case4 = new TH2F("spike_events_etaphi_case4", "spike_events_etaphi_case4", 200, 20, 20, 200, 20, 20);
	spike_events_etaphi_case4->GetXaxis()->SetTitle("eta");
	spike_events_etaphi_case4->GetYaxis()->SetTitle("phi");

    TH2F *spike_events_etaphi_case5 = new TH2F("spike_events_etaphi_case5", "spike_events_etaphi_case5", 200, 20, 20, 200, 20, 20);
	spike_events_etaphi_case5->GetXaxis()->SetTitle("eta");
	spike_events_etaphi_case5->GetYaxis()->SetTitle("phi");

    TH2F *spike_events_etaphi = new TH2F("spike_events_etaphi", "spike_events_etaphi", 200, 20, 20, 200, 20, 20);
	spike_events_etaphi->GetXaxis()->SetTitle("eta");
	spike_events_etaphi->GetYaxis()->SetTitle("phi");

    //spike etawing calculation
    //if there's only one
    TH1F *spike_eta_wing_rightonly1 = new TH1F("spike_eta_wing_rightonly1", "spike_eta_wing_rightonly1", 2000, 0, 1);
	spike_eta_wing_rightonly1->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_rightonly1->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_rightonly2 = new TH1F("spike_eta_wing_rightonly2", "spike_eta_wing_rightonly2", 2000, 0, 1);
	spike_eta_wing_rightonly2->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_rightonly2->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_rightonly3 = new TH1F("spike_eta_wing_rightonly3", "spike_eta_wing_rightonly3", 2000, 0, 1);
	spike_eta_wing_rightonly3->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_rightonly3->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_rightonly4 = new TH1F("spike_eta_wing_rightonly4", "spike_eta_wing_rightonly4", 2000, 0, 1);
	spike_eta_wing_rightonly4->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_rightonly4->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_rightonly5 = new TH1F("spike_eta_wing_rightonly5", "spike_eta_wing_rightonly5", 2000, 0, 1);
	spike_eta_wing_rightonly5->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_rightonly5->GetYaxis()->SetTitle("Entries");


    TH1F *spike_eta_wing_leftonly1 = new TH1F("spike_eta_wing_leftonly1", "spike_eta_wing_leftonly1", 2000, 0, 1);
	spike_eta_wing_leftonly1->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_leftonly1->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_leftonly2 = new TH1F("spike_eta_wing_leftonly2", "spike_eta_wing_leftonly2", 2000, 0, 1);
	spike_eta_wing_leftonly2->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_leftonly2->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_leftonly3 = new TH1F("spike_eta_wing_leftonly3", "spike_eta_wing_leftonly3", 2000, 0, 1);
	spike_eta_wing_leftonly3->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_leftonly3->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_leftonly4 = new TH1F("spike_eta_wing_leftonly4", "spike_eta_wing_leftonly4", 2000, 0, 1);
	spike_eta_wing_leftonly4->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_leftonly4->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_leftonly5 = new TH1F("spike_eta_wing_leftonly5", "spike_eta_wing_leftonly5", 2000, 0, 1);
	spike_eta_wing_leftonly5->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_leftonly5->GetYaxis()->SetTitle("Entries");

    //if there are both
    TH1F *spike_eta_wing_both1 = new TH1F("spike_eta_wing_both1", "spike_eta_wing_both1", 2000, 0, 1);
	spike_eta_wing_both1->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_both1->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_both2 = new TH1F("spike_eta_wing_both2", "spike_eta_wing_both2", 2000, 0, 1);
	spike_eta_wing_both2->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_both2->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_both3 = new TH1F("spike_eta_wing_both3", "spike_eta_wing_both3", 2000, 0, 1);
	spike_eta_wing_both3->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_both3->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_both4 = new TH1F("spike_eta_wing_both4", "spike_eta_wing_both4", 2000, 0, 1);
	spike_eta_wing_both4->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_both4->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_both5 = new TH1F("spike_eta_wing_both5", "spike_eta_wing_both5", 2000, 0, 1);
	spike_eta_wing_both5->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_both5->GetYaxis()->SetTitle("Entries");


    //combine all of the etawing
    TH1F *spike_eta_wing1 = new TH1F("spike_eta_wing1", "spike_eta_wing1", 2000, 0, 1);
	spike_eta_wing1->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing1->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing2 = new TH1F("spike_eta_wing2", "spike_eta_wing2", 2000, 0, 1);
	spike_eta_wing2->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing2->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing3 = new TH1F("spike_eta_wing3", "spike_eta_wing3", 2000, 0, 1);
	spike_eta_wing3->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing3->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing4 = new TH1F("spike_eta_wing4", "spike_eta_wing4", 2000, 0, 1);
	spike_eta_wing4->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing4->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing5 = new TH1F("spike_eta_wing5", "spike_eta_wing5", 2000, 0, 1);
	spike_eta_wing5->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing5->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing = new TH1F("spike_eta_wing", "spike_eta_wing", 2000, 0, 1);
	spike_eta_wing->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing->GetYaxis()->SetTitle("Entries");




    //with time < -12.5 ns cut
    TH1F *spike_events_time_case1 = new TH1F("spike_events_time_case1", "spike_events_time_case1", 200, -25, 25);
	spike_events_time_case1->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_time_case1->GetYaxis()->SetTitle("Entries");	  

	TH1F *spike_events_time_case2 = new TH1F("spike_events_time_case2", "spike_events_time_case2", 200, -25, 25);
	spike_events_time_case2->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_time_case2->GetYaxis()->SetTitle("Entries");

	TH1F *spike_events_time_case3 = new TH1F("spike_events_time_case3", "spike_events_time_case3", 200, -25, 25);
	spike_events_time_case3->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_time_case3->GetYaxis()->SetTitle("Entries");

	TH1F *spike_events_time_case4 = new TH1F("spike_events_time_case4", "spike_events_time_case4", 200, -25, 25);
	spike_events_time_case4->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_time_case4->GetYaxis()->SetTitle("Entries");

	TH1F *spike_events_time_case5 = new TH1F("spike_events_time_case5", "spike_events_time_case5", 200, -25, 25);
	spike_events_time_case5->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_time_case5->GetYaxis()->SetTitle("Entries");

	TH1F *spike_events_time = new TH1F("spike_events_time", "spike_events_time", 200, -25, 25);
	spike_events_time->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_time->GetYaxis()->SetTitle("Entries");

    TH1F *spike_events_time_3ns = new TH1F("spike_events_time_3ns", "spike_events_time_3ns", 200, -25, 25);
	spike_events_time_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_time_3ns->GetYaxis()->SetTitle("Entries");


    //eta vs phi for spike
    TH2F *spike_events_time_etaphi_case1 = new TH2F("spike_events_time_etaphi_case1", "spike_events_time_etaphi_case1", 200, 20, 20, 200, 20, 20);
	spike_events_time_etaphi_case1->GetXaxis()->SetTitle("eta");
	spike_events_time_etaphi_case1->GetYaxis()->SetTitle("phi");

    TH2F *spike_events_time_etaphi_case2 = new TH2F("spike_events_time_etaphi_case2", "spike_events_time_etaphi_case2", 200, 20, 20, 200, 20, 20);
	spike_events_time_etaphi_case2->GetXaxis()->SetTitle("eta");
	spike_events_time_etaphi_case2->GetYaxis()->SetTitle("phi");

    TH2F *spike_events_time_etaphi_case3 = new TH2F("spike_events_time_etaphi_case3", "spike_events_time_etaphi_case3", 200, 20, 20, 200, 20, 20);
	spike_events_time_etaphi_case3->GetXaxis()->SetTitle("eta");
	spike_events_time_etaphi_case3->GetYaxis()->SetTitle("phi");

    TH2F *spike_events_time_etaphi_case4 = new TH2F("spike_events_time_etaphi_case4", "spike_events_time_etaphi_case4", 200, 20, 20, 200, 20, 20);
	spike_events_time_etaphi_case4->GetXaxis()->SetTitle("eta");
	spike_events_time_etaphi_case4->GetYaxis()->SetTitle("phi");

    TH2F *spike_events_time_etaphi_case5 = new TH2F("spike_events_time_etaphi_case5", "spike_events_time_etaphi_case5", 200, 20, 20, 200, 20, 20);
	spike_events_time_etaphi_case5->GetXaxis()->SetTitle("eta");
	spike_events_time_etaphi_case5->GetYaxis()->SetTitle("phi");

    TH2F *spike_events_time_etaphi = new TH2F("spike_events_time_etaphi", "spike_events_time_etaphi", 200, 20, 20, 200, 20, 20);
	spike_events_time_etaphi->GetXaxis()->SetTitle("eta");
	spike_events_time_etaphi->GetYaxis()->SetTitle("phi");


    //spike etawing calculation
    //if there's only one
    TH1F *spike_eta_wing_time_rightonly1 = new TH1F("spike_eta_wing_time_rightonly1", "spike_eta_wing_time_rightonly1", 2000, 0, 1);
	spike_eta_wing_time_rightonly1->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time_rightonly1->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_time_rightonly2 = new TH1F("spike_eta_wing_time_rightonly2", "spike_eta_wing_time_rightonly2", 2000, 0, 1);
	spike_eta_wing_time_rightonly2->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time_rightonly2->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_time_rightonly3 = new TH1F("spike_eta_wing_time_rightonly3", "spike_eta_wing_time_rightonly3", 2000, 0, 1);
	spike_eta_wing_time_rightonly3->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time_rightonly3->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_time_rightonly4 = new TH1F("spike_eta_wing_time_rightonly4", "spike_eta_wing_time_rightonly4", 2000, 0, 1);
	spike_eta_wing_time_rightonly4->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time_rightonly4->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_time_rightonly5 = new TH1F("spike_eta_wing_time_rightonly5", "spike_eta_wing_time_rightonly5", 2000, 0, 1);
	spike_eta_wing_time_rightonly5->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time_rightonly5->GetYaxis()->SetTitle("Entries");


    TH1F *spike_eta_wing_time_leftonly1 = new TH1F("spike_eta_wing_time_leftonly1", "spike_eta_wing_time_leftonly1", 2000, 0, 1);
	spike_eta_wing_time_leftonly1->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time_leftonly1->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_time_leftonly2 = new TH1F("spike_eta_wing_time_leftonly2", "spike_eta_wing_time_leftonly2", 2000, 0, 1);
	spike_eta_wing_time_leftonly2->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time_leftonly2->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_time_leftonly3 = new TH1F("spike_eta_wing_time_leftonly3", "spike_eta_wing_time_leftonly3", 2000, 0, 1);
	spike_eta_wing_time_leftonly3->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time_leftonly3->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_time_leftonly4 = new TH1F("spike_eta_wing_time_leftonly4", "spike_eta_wing_time_leftonly4", 2000, 0, 1);
	spike_eta_wing_time_leftonly4->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time_leftonly4->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_time_leftonly5 = new TH1F("spike_eta_wing_time_leftonly5", "spike_eta_wing_time_leftonly5", 2000, 0, 1);
	spike_eta_wing_time_leftonly5->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time_leftonly5->GetYaxis()->SetTitle("Entries");

    //if there are both
    TH1F *spike_eta_wing_time_both1 = new TH1F("spike_eta_wing_time_both1", "spike_eta_wing_time_both1", 2000, 0, 1);
	spike_eta_wing_time_both1->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time_both1->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_time_both2 = new TH1F("spike_eta_wing_time_both2", "spike_eta_wing_time_both2", 2000, 0, 1);
	spike_eta_wing_time_both2->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time_both2->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_time_both3 = new TH1F("spike_eta_wing_time_both3", "spike_eta_wing_time_both3", 2000, 0, 1);
	spike_eta_wing_time_both3->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time_both3->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_time_both4 = new TH1F("spike_eta_wing_time_both4", "spike_eta_wing_time_both4", 2000, 0, 1);
	spike_eta_wing_time_both4->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time_both4->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_time_both5 = new TH1F("spike_eta_wing_time_both5", "spike_eta_wing_time_both5", 2000, 0, 1);
	spike_eta_wing_time_both5->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time_both5->GetYaxis()->SetTitle("Entries");


    //combine all of the etawing
    TH1F *spike_eta_wing_time1 = new TH1F("spike_eta_wing_time1", "spike_eta_wing_time1", 2000, 0, 1);
	spike_eta_wing_time1->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time1->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_time2 = new TH1F("spike_eta_wing_time2", "spike_eta_wing_time2", 2000, 0, 1);
	spike_eta_wing_time2->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time2->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_time3 = new TH1F("spike_eta_wing_time3", "spike_eta_wing_time3", 2000, 0, 1);
	spike_eta_wing_time3->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time3->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_time4 = new TH1F("spike_eta_wing_time4", "spike_eta_wing_time4", 2000, 0, 1);
	spike_eta_wing_time4->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time4->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_time5 = new TH1F("spike_eta_wing_time5", "spike_eta_wing_time5", 2000, 0, 1);
	spike_eta_wing_time5->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time5->GetYaxis()->SetTitle("Entries");

    TH1F *spike_eta_wing_time = new TH1F("spike_eta_wing_time", "spike_eta_wing_time", 2000, 0, 1);
	spike_eta_wing_time->GetXaxis()->SetTitle("E_wing");
	spike_eta_wing_time->GetYaxis()->SetTitle("Entries");



    //----------------------------------candidate events with Ewing > 0.01----------------------------------
    TH1F *candidate_events_ewing1_case1 = new TH1F("candidate_events_ewing1_case1", "candidate_events_ewing1_case1", 200, -25, 25);
	candidate_events_ewing1_case1->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing1_case1->GetYaxis()->SetTitle("Entries");	 

    TH1F *candidate_events_ewing1_case2 = new TH1F("candidate_events_ewing1_case2", "candidate_events_ewing1_case2", 200, -25, 25);
	candidate_events_ewing1_case2->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing1_case2->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing1_case3 = new TH1F("candidate_events_ewing1_case3", "candidate_events_ewing1_case3", 200, -25, 25);
	candidate_events_ewing1_case3->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing1_case3->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing1_case4 = new TH1F("candidate_events_ewing1_case4", "candidate_events_ewing1_case4", 200, -25, 25);
	candidate_events_ewing1_case4->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing1_case4->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing1_case5 = new TH1F("candidate_events_ewing1_case5", "candidate_events_ewing1_case5", 200, -25, 25);
	candidate_events_ewing1_case5->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing1_case5->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing1 = new TH1F("candidate_events_ewing1", "candidate_events_ewing1", 200, -25, 25);
	candidate_events_ewing1->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing1->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing1_3ns = new TH1F("candidate_events_ewing1_3ns", "candidate_events_ewing1_3ns", 200, -25, 25);
	candidate_events_ewing1_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing1_3ns->GetYaxis()->SetTitle("Entries");	

    //----------------------------------candidate events with Ewing > 0.03----------------------------------
    TH1F *candidate_events_ewing3_case1 = new TH1F("candidate_events_ewing3_case1", "candidate_events_ewing3_case1", 200, -25, 25);
	candidate_events_ewing3_case1->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing3_case1->GetYaxis()->SetTitle("Entries");	 

    TH1F *candidate_events_ewing3_case2 = new TH1F("candidate_events_ewing3_case2", "candidate_events_ewing3_case2", 200, -25, 25);
	candidate_events_ewing3_case2->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing3_case2->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing3_case3 = new TH1F("candidate_events_ewing3_case3", "candidate_events_ewing3_case3", 200, -25, 25);
	candidate_events_ewing3_case3->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing3_case3->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing3_case4 = new TH1F("candidate_events_ewing3_case4", "candidate_events_ewing3_case4", 200, -25, 25);
	candidate_events_ewing3_case4->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing3_case4->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing3_case5 = new TH1F("candidate_events_ewing3_case5", "candidate_events_ewing3_case5", 200, -25, 25);
	candidate_events_ewing3_case5->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing3_case5->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing3 = new TH1F("candidate_events_ewing3", "candidate_events_ewing3", 200, -25, 25);
	candidate_events_ewing3->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing3->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing3_3ns = new TH1F("candidate_events_ewing3_3ns", "candidate_events_ewing3_3ns", 200, -25, 25);
	candidate_events_ewing3_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing3_3ns->GetYaxis()->SetTitle("Entries");	



    //----------------------------------candidate events with Ewing > 0.02----------------------------------

    TH1F *candidate_events_ewing2 = new TH1F("candidate_events_ewing2", "candidate_events_ewing2", 200, -25, 25);
	candidate_events_ewing2->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing2->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing2_3ns = new TH1F("candidate_events_ewing2_3ns", "candidate_events_ewing2_3ns", 200, -25, 25);
	candidate_events_ewing2_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing2_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.025----------------------------------

    TH1F *candidate_events_ewing25 = new TH1F("candidate_events_ewing25", "candidate_events_ewing25", 200, -25, 25);
	candidate_events_ewing25->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing25->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing25_3ns = new TH1F("candidate_events_ewing25_3ns", "candidate_events_ewing25_3ns", 200, -25, 25);
	candidate_events_ewing25_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing25_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.031----------------------------------

    TH1F *candidate_events_ewing31 = new TH1F("candidate_events_ewing31", "candidate_events_ewing31", 200, -25, 25);
	candidate_events_ewing31->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing31->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing31_3ns = new TH1F("candidate_events_ewing31_3ns", "candidate_events_ewing31_3ns", 200, -25, 25);
	candidate_events_ewing31_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing31_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.032----------------------------------

    TH1F *candidate_events_ewing32 = new TH1F("candidate_events_ewing32", "candidate_events_ewing32", 200, -25, 25);
	candidate_events_ewing32->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing32->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing32_3ns = new TH1F("candidate_events_ewing32_3ns", "candidate_events_ewing32_3ns", 200, -25, 25);
	candidate_events_ewing32_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing32_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.033----------------------------------

    TH1F *candidate_events_ewing33 = new TH1F("candidate_events_ewing33", "candidate_events_ewing33", 200, -25, 25);
	candidate_events_ewing33->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing33->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing33_3ns = new TH1F("candidate_events_ewing33_3ns", "candidate_events_ewing33_3ns", 200, -25, 25);
	candidate_events_ewing33_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing33_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.034----------------------------------

    TH1F *candidate_events_ewing34 = new TH1F("candidate_events_ewing34", "candidate_events_ewing34", 200, -25, 25);
	candidate_events_ewing34->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing34->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing34_3ns = new TH1F("candidate_events_ewing34_3ns", "candidate_events_ewing34_3ns", 200, -25, 25);
	candidate_events_ewing34_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing34_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.035----------------------------------

    TH1F *candidate_events_ewing35 = new TH1F("candidate_events_ewing35", "candidate_events_ewing35", 200, -25, 25);
	candidate_events_ewing35->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing35->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing35_3ns = new TH1F("candidate_events_ewing35_3ns", "candidate_events_ewing35_3ns", 200, -25, 25);
	candidate_events_ewing35_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing35_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.036----------------------------------

    TH1F *candidate_events_ewing36 = new TH1F("candidate_events_ewing36", "candidate_events_ewing36", 200, -25, 25);
	candidate_events_ewing36->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing36->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing36_3ns = new TH1F("candidate_events_ewing36_3ns", "candidate_events_ewing36_3ns", 200, -25, 25);
	candidate_events_ewing36_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing36_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.037----------------------------------

    TH1F *candidate_events_ewing37 = new TH1F("candidate_events_ewing37", "candidate_events_ewing37", 200, -25, 25);
	candidate_events_ewing37->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing37->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing37_3ns = new TH1F("candidate_events_ewing37_3ns", "candidate_events_ewing37_3ns", 200, -25, 25);
	candidate_events_ewing37_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing37_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.038----------------------------------

    TH1F *candidate_events_ewing38 = new TH1F("candidate_events_ewing38", "candidate_events_ewing38", 200, -25, 25);
	candidate_events_ewing38->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing38->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing38_3ns = new TH1F("candidate_events_ewing38_3ns", "candidate_events_ewing38_3ns", 200, -25, 25);
	candidate_events_ewing38_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing38_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.039----------------------------------

    TH1F *candidate_events_ewing39 = new TH1F("candidate_events_ewing39", "candidate_events_ewing39", 200, -25, 25);
	candidate_events_ewing39->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing39->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing39_3ns = new TH1F("candidate_events_ewing39_3ns", "candidate_events_ewing39_3ns", 200, -25, 25);
	candidate_events_ewing39_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing39_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.04----------------------------------

    TH1F *candidate_events_ewing4 = new TH1F("candidate_events_ewing4", "candidate_events_ewing4", 200, -25, 25);
	candidate_events_ewing4->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing4->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing4_3ns = new TH1F("candidate_events_ewing4_3ns", "candidate_events_ewing4_3ns", 200, -25, 25);
	candidate_events_ewing4_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing4_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.041----------------------------------

    TH1F *candidate_events_ewing41 = new TH1F("candidate_events_ewing41", "candidate_events_ewing41", 200, -25, 25);
	candidate_events_ewing41->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing41->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing41_3ns = new TH1F("candidate_events_ewing41_3ns", "candidate_events_ewing41_3ns", 200, -25, 25);
	candidate_events_ewing41_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing41_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.042----------------------------------

    TH1F *candidate_events_ewing42 = new TH1F("candidate_events_ewing42", "candidate_events_ewing42", 200, -25, 25);
	candidate_events_ewing42->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing42->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing42_3ns = new TH1F("candidate_events_ewing42_3ns", "candidate_events_ewing42_3ns", 200, -25, 25);
	candidate_events_ewing42_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing42_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.043----------------------------------

    TH1F *candidate_events_ewing43 = new TH1F("candidate_events_ewing43", "candidate_events_ewing43", 200, -25, 25);
	candidate_events_ewing43->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing43->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing43_3ns = new TH1F("candidate_events_ewing43_3ns", "candidate_events_ewing43_3ns", 200, -25, 25);
	candidate_events_ewing43_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing43_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.044----------------------------------

    TH1F *candidate_events_ewing44 = new TH1F("candidate_events_ewing44", "candidate_events_ewing44", 200, -25, 25);
	candidate_events_ewing44->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing44->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing44_3ns = new TH1F("candidate_events_ewing44_3ns", "candidate_events_ewing44_3ns", 200, -25, 25);
	candidate_events_ewing44_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing44_3ns->GetYaxis()->SetTitle("Entries");
    

    //----------------------------------candidate events with Ewing > 0.045----------------------------------

    TH1F *candidate_events_ewing45 = new TH1F("candidate_events_ewing45", "candidate_events_ewing45", 200, -25, 25);
	candidate_events_ewing45->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing45->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing45_3ns = new TH1F("candidate_events_ewing45_3ns", "candidate_events_ewing45_3ns", 200, -25, 25);
	candidate_events_ewing45_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing45_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.046----------------------------------

    TH1F *candidate_events_ewing46 = new TH1F("candidate_events_ewing46", "candidate_events_ewing46", 200, -25, 25);
	candidate_events_ewing46->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing46->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing46_3ns = new TH1F("candidate_events_ewing46_3ns", "candidate_events_ewing46_3ns", 200, -25, 25);
	candidate_events_ewing46_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing46_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.047----------------------------------

    TH1F *candidate_events_ewing47 = new TH1F("candidate_events_ewing47", "candidate_events_ewing47", 200, -25, 25);
	candidate_events_ewing47->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing47->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing47_3ns = new TH1F("candidate_events_ewing47_3ns", "candidate_events_ewing47_3ns", 200, -25, 25);
	candidate_events_ewing47_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing47_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.048----------------------------------

    TH1F *candidate_events_ewing48 = new TH1F("candidate_events_ewing48", "candidate_events_ewing48", 200, -25, 25);
	candidate_events_ewing48->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing48->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing48_3ns = new TH1F("candidate_events_ewing48_3ns", "candidate_events_ewing48_3ns", 200, -25, 25);
	candidate_events_ewing48_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing48_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.049----------------------------------

    TH1F *candidate_events_ewing49 = new TH1F("candidate_events_ewing49", "candidate_events_ewing49", 200, -25, 25);
	candidate_events_ewing49->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing49->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing49_3ns = new TH1F("candidate_events_ewing49_3ns", "candidate_events_ewing49_3ns", 200, -25, 25);
	candidate_events_ewing49_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing49_3ns->GetYaxis()->SetTitle("Entries");
    



    //----------------------------------candidate events with Ewing > 0.05----------------------------------

    TH1F *candidate_events_ewing5 = new TH1F("candidate_events_ewing5", "candidate_events_ewing5", 200, -25, 25);
	candidate_events_ewing5->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing5->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing5_3ns = new TH1F("candidate_events_ewing5_3ns", "candidate_events_ewing5_3ns", 200, -25, 25);
	candidate_events_ewing5_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing5_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.07----------------------------------

    TH1F *candidate_events_ewing7 = new TH1F("candidate_events_ewing7", "candidate_events_ewing7", 200, -25, 25);
	candidate_events_ewing7->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing7->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing7_3ns = new TH1F("candidate_events_ewing7_3ns", "candidate_events_ewing7_3ns", 200, -25, 25);
	candidate_events_ewing7_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing7_3ns->GetYaxis()->SetTitle("Entries");	

    //----------------------------------candidate events with Ewing > 0.09----------------------------------

    TH1F *candidate_events_ewing9 = new TH1F("candidate_events_ewing9", "candidate_events_ewing9", 200, -25, 25);
	candidate_events_ewing9->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing9->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing9_3ns = new TH1F("candidate_events_ewing9_3ns", "candidate_events_ewing9_3ns", 200, -25, 25);
	candidate_events_ewing9_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing9_3ns->GetYaxis()->SetTitle("Entries");	

    //----------------------------------candidate events with Ewing > 0.1----------------------------------

    TH1F *candidate_events_ewing10 = new TH1F("candidate_events_ewing10", "candidate_events_ewing10", 200, -25, 25);
	candidate_events_ewing10->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing10->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing10_3ns = new TH1F("candidate_events_ewing10_3ns", "candidate_events_ewing10_3ns", 200, -25, 25);
	candidate_events_ewing10_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing10_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.12----------------------------------

    TH1F *candidate_events_ewing12 = new TH1F("candidate_events_ewing12", "candidate_events_ewing12", 200, -25, 25);
	candidate_events_ewing12->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing12->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing12_3ns = new TH1F("candidate_events_ewing12_3ns", "candidate_events_ewing12_3ns", 200, -25, 25);
	candidate_events_ewing12_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing12_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.14----------------------------------

    TH1F *candidate_events_ewing14 = new TH1F("candidate_events_ewing14", "candidate_events_ewing14", 200, -25, 25);
	candidate_events_ewing14->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing14->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing14_3ns = new TH1F("candidate_events_ewing14_3ns", "candidate_events_ewing14_3ns", 200, -25, 25);
	candidate_events_ewing14_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing14_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------candidate events with Ewing > 0.15----------------------------------

    TH1F *candidate_events_ewing15 = new TH1F("candidate_events_ewing15", "candidate_events_ewing15", 200, -25, 25);
	candidate_events_ewing15->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing15->GetYaxis()->SetTitle("Entries");	

    TH1F *candidate_events_ewing15_3ns = new TH1F("candidate_events_ewing15_3ns", "candidate_events_ewing15_3ns", 200, -25, 25);
	candidate_events_ewing15_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	candidate_events_ewing15_3ns->GetYaxis()->SetTitle("Entries");






    //----------------------------------spikeincandidate events with Ewing > 0.01----------------------------------
    TH1F *spikeincandidate_events_ewing1_case1 = new TH1F("spikeincandidate_events_ewing1_case1", "spikeincandidate_events_ewing1_case1", 200, -25, 25);
	spikeincandidate_events_ewing1_case1->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing1_case1->GetYaxis()->SetTitle("Entries");	 

    TH1F *spikeincandidate_events_ewing1_case2 = new TH1F("spikeincandidate_events_ewing1_case2", "spikeincandidate_events_ewing1_case2", 200, -25, 25);
	spikeincandidate_events_ewing1_case2->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing1_case2->GetYaxis()->SetTitle("Entries");	

    TH1F *spikeincandidate_events_ewing1_case3 = new TH1F("spikeincandidate_events_ewing1_case3", "spikeincandidate_events_ewing1_case3", 200, -25, 25);
	spikeincandidate_events_ewing1_case3->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing1_case3->GetYaxis()->SetTitle("Entries");	

    TH1F *spikeincandidate_events_ewing1_case4 = new TH1F("spikeincandidate_events_ewing1_case4", "spikeincandidate_events_ewing1_case4", 200, -25, 25);
	spikeincandidate_events_ewing1_case4->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing1_case4->GetYaxis()->SetTitle("Entries");	

    TH1F *spikeincandidate_events_ewing1_case5 = new TH1F("spikeincandidate_events_ewing1_case5", "spikeincandidate_events_ewing1_case5", 200, -25, 25);
	spikeincandidate_events_ewing1_case5->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing1_case5->GetYaxis()->SetTitle("Entries");	

    TH1F *spikeincandidate_events_ewing1 = new TH1F("spikeincandidate_events_ewing1", "spikeincandidate_events_ewing1", 200, -25, 25);
	spikeincandidate_events_ewing1->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing1->GetYaxis()->SetTitle("Entries");	


    //----------------------------------spikeincandidate events with Ewing > 0.03----------------------------------
    TH1F *spikeincandidate_events_ewing3_case1 = new TH1F("spikeincandidate_events_ewing3_case1", "spikeincandidate_events_ewing3_case1", 200, -25, 25);
	spikeincandidate_events_ewing3_case1->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing3_case1->GetYaxis()->SetTitle("Entries");	 

    TH1F *spikeincandidate_events_ewing3_case2 = new TH1F("spikeincandidate_events_ewing3_case2", "spikeincandidate_events_ewing3_case2", 200, -25, 25);
	spikeincandidate_events_ewing3_case2->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing3_case2->GetYaxis()->SetTitle("Entries");	

    TH1F *spikeincandidate_events_ewing3_case3 = new TH1F("spikeincandidate_events_ewing3_case3", "spikeincandidate_events_ewing3_case3", 200, -25, 25);
	spikeincandidate_events_ewing3_case3->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing3_case3->GetYaxis()->SetTitle("Entries");	

    TH1F *spikeincandidate_events_ewing3_case4 = new TH1F("spikeincandidate_events_ewing3_case4", "spikeincandidate_events_ewing3_case4", 200, -25, 25);
	spikeincandidate_events_ewing3_case4->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing3_case4->GetYaxis()->SetTitle("Entries");	

    TH1F *spikeincandidate_events_ewing3_case5 = new TH1F("spikeincandidate_events_ewing3_case5", "spikeincandidate_events_ewing3_case5", 200, -25, 25);
	spikeincandidate_events_ewing3_case5->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing3_case5->GetYaxis()->SetTitle("Entries");	

    TH1F *spikeincandidate_events_ewing3 = new TH1F("spikeincandidate_events_ewing3", "spikeincandidate_events_ewing3", 200, -25, 25);
	spikeincandidate_events_ewing3->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing3->GetYaxis()->SetTitle("Entries");	


    //----------------------------------spikeincandidate events with Ewing > 0.02----------------------------------

    TH1F *spikeincandidate_events_ewing2 = new TH1F("spikeincandidate_events_ewing2", "spikeincandidate_events_ewing2", 200, -25, 25);
	spikeincandidate_events_ewing2->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing2->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.025----------------------------------

    TH1F *spikeincandidate_events_ewing25 = new TH1F("spikeincandidate_events_ewing25", "spikeincandidate_events_ewing25", 200, -25, 25);
	spikeincandidate_events_ewing25->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing25->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.031----------------------------------

    TH1F *spikeincandidate_events_ewing31 = new TH1F("spikeincandidate_events_ewing31", "spikeincandidate_events_ewing31", 200, -25, 25);
	spikeincandidate_events_ewing31->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing31->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.032----------------------------------

    TH1F *spikeincandidate_events_ewing32 = new TH1F("spikeincandidate_events_ewing32", "spikeincandidate_events_ewing32", 200, -25, 25);
	spikeincandidate_events_ewing32->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing32->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.033----------------------------------

    TH1F *spikeincandidate_events_ewing33 = new TH1F("spikeincandidate_events_ewing33", "spikeincandidate_events_ewing33", 200, -25, 25);
	spikeincandidate_events_ewing33->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing33->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.034----------------------------------

    TH1F *spikeincandidate_events_ewing34 = new TH1F("spikeincandidate_events_ewing34", "spikeincandidate_events_ewing34", 200, -25, 25);
	spikeincandidate_events_ewing34->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing34->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.035----------------------------------

    TH1F *spikeincandidate_events_ewing35 = new TH1F("spikeincandidate_events_ewing35", "spikeincandidate_events_ewing35", 200, -25, 25);
	spikeincandidate_events_ewing35->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing35->GetYaxis()->SetTitle("Entries");	


    //----------------------------------spikeincandidate events with Ewing > 0.036----------------------------------

    TH1F *spikeincandidate_events_ewing36 = new TH1F("spikeincandidate_events_ewing36", "spikeincandidate_events_ewing36", 200, -25, 25);
	spikeincandidate_events_ewing36->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing36->GetYaxis()->SetTitle("Entries");	


    //----------------------------------spikeincandidate events with Ewing > 0.037----------------------------------

    TH1F *spikeincandidate_events_ewing37 = new TH1F("spikeincandidate_events_ewing37", "spikeincandidate_events_ewing37", 200, -25, 25);
	spikeincandidate_events_ewing37->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing37->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.038----------------------------------

    TH1F *spikeincandidate_events_ewing38 = new TH1F("spikeincandidate_events_ewing38", "spikeincandidate_events_ewing38", 200, -25, 25);
	spikeincandidate_events_ewing38->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing38->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.039----------------------------------

    TH1F *spikeincandidate_events_ewing39 = new TH1F("spikeincandidate_events_ewing39", "spikeincandidate_events_ewing39", 200, -25, 25);
	spikeincandidate_events_ewing39->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing39->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.04----------------------------------

    TH1F *spikeincandidate_events_ewing4 = new TH1F("spikeincandidate_events_ewing4", "spikeincandidate_events_ewing4", 200, -25, 25);
	spikeincandidate_events_ewing4->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing4->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.041----------------------------------

    TH1F *spikeincandidate_events_ewing41 = new TH1F("spikeincandidate_events_ewing41", "spikeincandidate_events_ewing41", 200, -25, 25);
	spikeincandidate_events_ewing41->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing41->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.042----------------------------------

    TH1F *spikeincandidate_events_ewing42 = new TH1F("spikeincandidate_events_ewing42", "spikeincandidate_events_ewing42", 200, -25, 25);
	spikeincandidate_events_ewing42->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing42->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.043----------------------------------

    TH1F *spikeincandidate_events_ewing43 = new TH1F("spikeincandidate_events_ewing43", "spikeincandidate_events_ewing43", 200, -25, 25);
	spikeincandidate_events_ewing43->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing43->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.044----------------------------------

    TH1F *spikeincandidate_events_ewing44 = new TH1F("spikeincandidate_events_ewing44", "spikeincandidate_events_ewing44", 200, -25, 25);
	spikeincandidate_events_ewing44->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing44->GetYaxis()->SetTitle("Entries");	
    

    //----------------------------------spikeincandidate events with Ewing > 0.045----------------------------------

    TH1F *spikeincandidate_events_ewing45 = new TH1F("spikeincandidate_events_ewing45", "spikeincandidate_events_ewing45", 200, -25, 25);
	spikeincandidate_events_ewing45->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing45->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.046----------------------------------

    TH1F *spikeincandidate_events_ewing46 = new TH1F("spikeincandidate_events_ewing46", "spikeincandidate_events_ewing46", 200, -25, 25);
	spikeincandidate_events_ewing46->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing46->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.047----------------------------------

    TH1F *spikeincandidate_events_ewing47 = new TH1F("spikeincandidate_events_ewing47", "spikeincandidate_events_ewing47", 200, -25, 25);
	spikeincandidate_events_ewing47->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing47->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.048----------------------------------

    TH1F *spikeincandidate_events_ewing48 = new TH1F("spikeincandidate_events_ewing48", "spikeincandidate_events_ewing48", 200, -25, 25);
	spikeincandidate_events_ewing48->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing48->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.049----------------------------------

    TH1F *spikeincandidate_events_ewing49 = new TH1F("spikeincandidate_events_ewing49", "spikeincandidate_events_ewing49", 200, -25, 25);
	spikeincandidate_events_ewing49->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing49->GetYaxis()->SetTitle("Entries");	



    //----------------------------------spikeincandidate events with Ewing > 0.05----------------------------------

    TH1F *spikeincandidate_events_ewing5 = new TH1F("spikeincandidate_events_ewing5", "spikeincandidate_events_ewing5", 200, -25, 25);
	spikeincandidate_events_ewing5->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing5->GetYaxis()->SetTitle("Entries");	


    //----------------------------------spikeincandidate events with Ewing > 0.07----------------------------------

    TH1F *spikeincandidate_events_ewing7 = new TH1F("spikeincandidate_events_ewing7", "spikeincandidate_events_ewing7", 200, -25, 25);
	spikeincandidate_events_ewing7->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing7->GetYaxis()->SetTitle("Entries");	


    //----------------------------------spikeincandidate events with Ewing > 0.09----------------------------------

    TH1F *spikeincandidate_events_ewing9 = new TH1F("spikeincandidate_events_ewing9", "spikeincandidate_events_ewing9", 200, -25, 25);
	spikeincandidate_events_ewing9->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing9->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.1----------------------------------

    TH1F *spikeincandidate_events_ewing10 = new TH1F("spikeincandidate_events_ewing10", "spikeincandidate_events_ewing10", 200, -25, 25);
	spikeincandidate_events_ewing10->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing10->GetYaxis()->SetTitle("Entries");	


    //----------------------------------spikeincandidate events with Ewing > 0.12----------------------------------

    TH1F *spikeincandidate_events_ewing12 = new TH1F("spikeincandidate_events_ewing12", "spikeincandidate_events_ewing12", 200, -25, 25);
	spikeincandidate_events_ewing12->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing12->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.14----------------------------------

    TH1F *spikeincandidate_events_ewing14 = new TH1F("spikeincandidate_events_ewing14", "spikeincandidate_events_ewing14", 200, -25, 25);
	spikeincandidate_events_ewing14->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing14->GetYaxis()->SetTitle("Entries");	

    //----------------------------------spikeincandidate events with Ewing > 0.15----------------------------------

    TH1F *spikeincandidate_events_ewing15 = new TH1F("spikeincandidate_events_ewing15", "spikeincandidate_events_ewing15", 200, -25, 25);
	spikeincandidate_events_ewing15->GetXaxis()->SetTitle("SeedTime (ns)");
	spikeincandidate_events_ewing15->GetYaxis()->SetTitle("Entries");	




    //----------------------------------promptZ events with Ewing > 0.01----------------------------------
    TH1F *promptZ_events_ewing1_case1 = new TH1F("promptZ_events_ewing1_case1", "promptZ_events_ewing1_case1", 200, -25, 25);
	promptZ_events_ewing1_case1->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing1_case1->GetYaxis()->SetTitle("Entries");	 

    TH1F *promptZ_events_ewing1_case2 = new TH1F("promptZ_events_ewing1_case2", "promptZ_events_ewing1_case2", 200, -25, 25);
	promptZ_events_ewing1_case2->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing1_case2->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing1_case3 = new TH1F("promptZ_events_ewing1_case3", "promptZ_events_ewing1_case3", 200, -25, 25);
	promptZ_events_ewing1_case3->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing1_case3->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing1_case4 = new TH1F("promptZ_events_ewing1_case4", "promptZ_events_ewing1_case4", 200, -25, 25);
	promptZ_events_ewing1_case4->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing1_case4->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing1_case5 = new TH1F("promptZ_events_ewing1_case5", "promptZ_events_ewing1_case5", 200, -25, 25);
	promptZ_events_ewing1_case5->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing1_case5->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing1 = new TH1F("promptZ_events_ewing1", "promptZ_events_ewing1", 200, -25, 25);
	promptZ_events_ewing1->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing1->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing1_3ns = new TH1F("promptZ_events_ewing1_3ns", "promptZ_events_ewing1_3ns", 200, -25, 25);
	promptZ_events_ewing1_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing1_3ns->GetYaxis()->SetTitle("Entries");	


    //----------------------------------promptZ events with Ewing > 0.03----------------------------------
    TH1F *promptZ_events_ewing3_case1 = new TH1F("promptZ_events_ewing3_case1", "promptZ_events_ewing3_case1", 200, -25, 25);
	promptZ_events_ewing3_case1->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing3_case1->GetYaxis()->SetTitle("Entries");	 

    TH1F *promptZ_events_ewing3_case2 = new TH1F("promptZ_events_ewing3_case2", "promptZ_events_ewing3_case2", 200, -25, 25);
	promptZ_events_ewing3_case2->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing3_case2->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing3_case3 = new TH1F("promptZ_events_ewing3_case3", "promptZ_events_ewing3_case3", 200, -25, 25);
	promptZ_events_ewing3_case3->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing3_case3->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing3_case4 = new TH1F("promptZ_events_ewing3_case4", "promptZ_events_ewing3_case4", 200, -25, 25);
	promptZ_events_ewing3_case4->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing3_case4->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing3_case5 = new TH1F("promptZ_events_ewing3_case5", "promptZ_events_ewing3_case5", 200, -25, 25);
	promptZ_events_ewing3_case5->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing3_case5->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing3 = new TH1F("promptZ_events_ewing3", "promptZ_events_ewing3", 200, -25, 25);
	promptZ_events_ewing3->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing3->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing3_3ns = new TH1F("promptZ_events_ewing3_3ns", "promptZ_events_ewing3_3ns", 200, -25, 25);
	promptZ_events_ewing3_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing3_3ns->GetYaxis()->SetTitle("Entries");	
    


    //----------------------------------promptZ events with Ewing > 0.02----------------------------------

    TH1F *promptZ_events_ewing2 = new TH1F("promptZ_events_ewing2", "promptZ_events_ewing2", 200, -25, 25);
	promptZ_events_ewing2->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing2->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing2_3ns = new TH1F("promptZ_events_ewing2_3ns", "promptZ_events_ewing2_3ns", 200, -25, 25);
	promptZ_events_ewing2_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing2_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.025----------------------------------

    TH1F *promptZ_events_ewing25 = new TH1F("promptZ_events_ewing25", "promptZ_events_ewing25", 200, -25, 25);
	promptZ_events_ewing25->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing25->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing25_3ns = new TH1F("promptZ_events_ewing25_3ns", "promptZ_events_ewing25_3ns", 200, -25, 25);
	promptZ_events_ewing25_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing25_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.031----------------------------------

    TH1F *promptZ_events_ewing31 = new TH1F("promptZ_events_ewing31", "promptZ_events_ewing31", 200, -25, 25);
	promptZ_events_ewing31->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing31->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing31_3ns = new TH1F("promptZ_events_ewing31_3ns", "promptZ_events_ewing31_3ns", 200, -25, 25);
	promptZ_events_ewing31_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing31_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.032----------------------------------

    TH1F *promptZ_events_ewing32 = new TH1F("promptZ_events_ewing32", "promptZ_events_ewing32", 200, -25, 25);
	promptZ_events_ewing32->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing32->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing32_3ns = new TH1F("promptZ_events_ewing32_3ns", "promptZ_events_ewing32_3ns", 200, -25, 25);
	promptZ_events_ewing32_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing32_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.033----------------------------------

    TH1F *promptZ_events_ewing33 = new TH1F("promptZ_events_ewing33", "promptZ_events_ewing33", 200, -25, 25);
	promptZ_events_ewing33->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing33->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing33_3ns = new TH1F("promptZ_events_ewing33_3ns", "promptZ_events_ewing33_3ns", 200, -25, 25);
	promptZ_events_ewing33_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing33_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.034----------------------------------

    TH1F *promptZ_events_ewing34 = new TH1F("promptZ_events_ewing34", "promptZ_events_ewing34", 200, -25, 25);
	promptZ_events_ewing34->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing34->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing34_3ns = new TH1F("promptZ_events_ewing34_3ns", "promptZ_events_ewing34_3ns", 200, -25, 25);
	promptZ_events_ewing34_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing34_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.035----------------------------------

    TH1F *promptZ_events_ewing35 = new TH1F("promptZ_events_ewing35", "promptZ_events_ewing35", 200, -25, 25);
	promptZ_events_ewing35->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing35->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing35_3ns = new TH1F("promptZ_events_ewing35_3ns", "promptZ_events_ewing35_3ns", 200, -25, 25);
	promptZ_events_ewing35_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing35_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.036----------------------------------

    TH1F *promptZ_events_ewing36 = new TH1F("promptZ_events_ewing36", "promptZ_events_ewing36", 200, -25, 25);
	promptZ_events_ewing36->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing36->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing36_3ns = new TH1F("promptZ_events_ewing36_3ns", "promptZ_events_ewing36_3ns", 200, -25, 25);
	promptZ_events_ewing36_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing36_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.037----------------------------------

    TH1F *promptZ_events_ewing37 = new TH1F("promptZ_events_ewing37", "promptZ_events_ewing37", 200, -25, 25);
	promptZ_events_ewing37->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing37->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing37_3ns = new TH1F("promptZ_events_ewing37_3ns", "promptZ_events_ewing37_3ns", 200, -25, 25);
	promptZ_events_ewing37_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing37_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.038----------------------------------

    TH1F *promptZ_events_ewing38 = new TH1F("promptZ_events_ewing38", "promptZ_events_ewing38", 200, -25, 25);
	promptZ_events_ewing38->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing38->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing38_3ns = new TH1F("promptZ_events_ewing38_3ns", "promptZ_events_ewing38_3ns", 200, -25, 25);
	promptZ_events_ewing38_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing38_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.039----------------------------------

    TH1F *promptZ_events_ewing39 = new TH1F("promptZ_events_ewing39", "promptZ_events_ewing39", 200, -25, 25);
	promptZ_events_ewing39->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing39->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing39_3ns = new TH1F("promptZ_events_ewing39_3ns", "promptZ_events_ewing39_3ns", 200, -25, 25);
	promptZ_events_ewing39_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing39_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.04----------------------------------

    TH1F *promptZ_events_ewing4 = new TH1F("promptZ_events_ewing4", "promptZ_events_ewing4", 200, -25, 25);
	promptZ_events_ewing4->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing4->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing4_3ns = new TH1F("promptZ_events_ewing4_3ns", "promptZ_events_ewing4_3ns", 200, -25, 25);
	promptZ_events_ewing4_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing4_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.041----------------------------------

    TH1F *promptZ_events_ewing41 = new TH1F("promptZ_events_ewing41", "promptZ_events_ewing41", 200, -25, 25);
	promptZ_events_ewing41->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing41->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing41_3ns = new TH1F("promptZ_events_ewing41_3ns", "promptZ_events_ewing41_3ns", 200, -25, 25);
	promptZ_events_ewing41_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing41_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.042----------------------------------

    TH1F *promptZ_events_ewing42 = new TH1F("promptZ_events_ewing42", "promptZ_events_ewing42", 200, -25, 25);
	promptZ_events_ewing42->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing42->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing42_3ns = new TH1F("promptZ_events_ewing42_3ns", "promptZ_events_ewing42_3ns", 200, -25, 25);
	promptZ_events_ewing42_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing42_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.043----------------------------------

    TH1F *promptZ_events_ewing43 = new TH1F("promptZ_events_ewing43", "promptZ_events_ewing43", 200, -25, 25);
	promptZ_events_ewing43->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing43->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing43_3ns = new TH1F("promptZ_events_ewing43_3ns", "promptZ_events_ewing43_3ns", 200, -25, 25);
	promptZ_events_ewing43_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing43_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.044----------------------------------

    TH1F *promptZ_events_ewing44 = new TH1F("promptZ_events_ewing44", "promptZ_events_ewing44", 200, -25, 25);
	promptZ_events_ewing44->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing44->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing44_3ns = new TH1F("promptZ_events_ewing44_3ns", "promptZ_events_ewing44_3ns", 200, -25, 25);
	promptZ_events_ewing44_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing44_3ns->GetYaxis()->SetTitle("Entries");
    

    //----------------------------------promptZ events with Ewing > 0.045----------------------------------

    TH1F *promptZ_events_ewing45 = new TH1F("promptZ_events_ewing45", "promptZ_events_ewing45", 200, -25, 25);
	promptZ_events_ewing45->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing45->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing45_3ns = new TH1F("promptZ_events_ewing45_3ns", "promptZ_events_ewing45_3ns", 200, -25, 25);
	promptZ_events_ewing45_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing45_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.046----------------------------------

    TH1F *promptZ_events_ewing46 = new TH1F("promptZ_events_ewing46", "promptZ_events_ewing46", 200, -25, 25);
	promptZ_events_ewing46->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing46->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing46_3ns = new TH1F("promptZ_events_ewing46_3ns", "promptZ_events_ewing46_3ns", 200, -25, 25);
	promptZ_events_ewing46_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing46_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.047----------------------------------

    TH1F *promptZ_events_ewing47 = new TH1F("promptZ_events_ewing47", "promptZ_events_ewing47", 200, -25, 25);
	promptZ_events_ewing47->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing47->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing47_3ns = new TH1F("promptZ_events_ewing47_3ns", "promptZ_events_ewing47_3ns", 200, -25, 25);
	promptZ_events_ewing47_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing47_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.048----------------------------------

    TH1F *promptZ_events_ewing48 = new TH1F("promptZ_events_ewing48", "promptZ_events_ewing48", 200, -25, 25);
	promptZ_events_ewing48->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing48->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing48_3ns = new TH1F("promptZ_events_ewing48_3ns", "promptZ_events_ewing48_3ns", 200, -25, 25);
	promptZ_events_ewing48_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing48_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.049----------------------------------

    TH1F *promptZ_events_ewing49 = new TH1F("promptZ_events_ewing49", "promptZ_events_ewing49", 200, -25, 25);
	promptZ_events_ewing49->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing49->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing49_3ns = new TH1F("promptZ_events_ewing49_3ns", "promptZ_events_ewing49_3ns", 200, -25, 25);
	promptZ_events_ewing49_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing49_3ns->GetYaxis()->SetTitle("Entries");



    //----------------------------------promptZ events with Ewing > 0.05----------------------------------

    TH1F *promptZ_events_ewing5 = new TH1F("promptZ_events_ewing5", "promptZ_events_ewing5", 200, -25, 25);
	promptZ_events_ewing5->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing5->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing5_3ns = new TH1F("promptZ_events_ewing5_3ns", "promptZ_events_ewing5_3ns", 200, -25, 25);
	promptZ_events_ewing5_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing5_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.07----------------------------------

    TH1F *promptZ_events_ewing7 = new TH1F("promptZ_events_ewing7", "promptZ_events_ewing7", 200, -25, 25);
	promptZ_events_ewing7->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing7->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing7_3ns = new TH1F("promptZ_events_ewing7_3ns", "promptZ_events_ewing7_3ns", 200, -25, 25);
	promptZ_events_ewing7_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing7_3ns->GetYaxis()->SetTitle("Entries");	



    //----------------------------------promptZ events with Ewing > 0.09----------------------------------

    TH1F *promptZ_events_ewing9 = new TH1F("promptZ_events_ewing9", "promptZ_events_ewing9", 200, -25, 25);
	promptZ_events_ewing9->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing9->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing9_3ns = new TH1F("promptZ_events_ewing9_3ns", "promptZ_events_ewing9_3ns", 200, -25, 25);
	promptZ_events_ewing9_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing9_3ns->GetYaxis()->SetTitle("Entries");	

    //----------------------------------promptZ events with Ewing > 0.1----------------------------------

    TH1F *promptZ_events_ewing10 = new TH1F("promptZ_events_ewing10", "promptZ_events_ewing10", 200, -25, 25);
	promptZ_events_ewing10->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing10->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing10_3ns = new TH1F("promptZ_events_ewing10_3ns", "promptZ_events_ewing10_3ns", 200, -25, 25);
	promptZ_events_ewing10_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing10_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.12----------------------------------

    TH1F *promptZ_events_ewing12 = new TH1F("promptZ_events_ewing12", "promptZ_events_ewing12", 200, -25, 25);
	promptZ_events_ewing12->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing12->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing12_3ns = new TH1F("promptZ_events_ewing12_3ns", "promptZ_events_ewing12_3ns", 200, -25, 25);
	promptZ_events_ewing12_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing12_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.14----------------------------------

    TH1F *promptZ_events_ewing14 = new TH1F("promptZ_events_ewing14", "promptZ_events_ewing14", 200, -25, 25);
	promptZ_events_ewing14->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing14->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing14_3ns = new TH1F("promptZ_events_ewing14_3ns", "promptZ_events_ewing14_3ns", 200, -25, 25);
	promptZ_events_ewing14_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing14_3ns->GetYaxis()->SetTitle("Entries");

    //----------------------------------promptZ events with Ewing > 0.15----------------------------------

    TH1F *promptZ_events_ewing15 = new TH1F("promptZ_events_ewing15", "promptZ_events_ewing15", 200, -25, 25);
	promptZ_events_ewing15->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing15->GetYaxis()->SetTitle("Entries");	

    TH1F *promptZ_events_ewing15_3ns = new TH1F("promptZ_events_ewing15_3ns", "promptZ_events_ewing15_3ns", 200, -25, 25);
	promptZ_events_ewing15_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	promptZ_events_ewing15_3ns->GetYaxis()->SetTitle("Entries");




    //----------------------------------spike events with Ewing > 0.01----------------------------------
    TH1F *spike_events_ewing1_case1 = new TH1F("spike_events_ewing1_case1", "spike_events_ewing1_case1", 200, -25, 25);
	spike_events_ewing1_case1->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing1_case1->GetYaxis()->SetTitle("Entries");	 

    TH1F *spike_events_ewing1_case2 = new TH1F("spike_events_ewing1_case2", "spike_events_ewing1_case2", 200, -25, 25);
	spike_events_ewing1_case2->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing1_case2->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing1_case3 = new TH1F("spike_events_ewing1_case3", "spike_events_ewing1_case3", 200, -25, 25);
	spike_events_ewing1_case3->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing1_case3->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing1_case4 = new TH1F("spike_events_ewing1_case4", "spike_events_ewing1_case4", 200, -25, 25);
	spike_events_ewing1_case4->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing1_case4->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing1_case5 = new TH1F("spike_events_ewing1_case5", "spike_events_ewing1_case5", 200, -25, 25);
	spike_events_ewing1_case5->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing1_case5->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing1 = new TH1F("spike_events_ewing1", "spike_events_ewing1", 200, -25, 25);
	spike_events_ewing1->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing1->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing1_3ns = new TH1F("spike_events_ewing1_3ns", "spike_events_ewing1_3ns", 200, -25, 25);
	spike_events_ewing1_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing1_3ns->GetYaxis()->SetTitle("Entries");	


    //----------------------------------spike events with Ewing > 0.03----------------------------------
    TH1F *spike_events_ewing3_case1 = new TH1F("spike_events_ewing3_case1", "spike_events_ewing3_case1", 200, -25, 25);
	spike_events_ewing3_case1->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing3_case1->GetYaxis()->SetTitle("Entries");	 

    TH1F *spike_events_ewing3_case2 = new TH1F("spike_events_ewing3_case2", "spike_events_ewing3_case2", 200, -25, 25);
	spike_events_ewing3_case2->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing3_case2->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing3_case3 = new TH1F("spike_events_ewing3_case3", "spike_events_ewing3_case3", 200, -25, 25);
	spike_events_ewing3_case3->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing3_case3->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing3_case4 = new TH1F("spike_events_ewing3_case4", "spike_events_ewing3_case4", 200, -25, 25);
	spike_events_ewing3_case4->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing3_case4->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing3_case5 = new TH1F("spike_events_ewing3_case5", "spike_events_ewing3_case5", 200, -25, 25);
	spike_events_ewing3_case5->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing3_case5->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing3 = new TH1F("spike_events_ewing3", "spike_events_ewing3", 200, -25, 25);
	spike_events_ewing3->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing3->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing3_3ns = new TH1F("spike_events_ewing3_3ns", "spike_events_ewing3_3ns", 200, -25, 25);
	spike_events_ewing3_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing3_3ns->GetYaxis()->SetTitle("Entries");	


    //----------------------------------spike events with Ewing > 0.01 time----------------------------------
    TH1F *spike_events_ewing1_time_case1 = new TH1F("spike_events_ewing1_time_case1", "spike_events_ewing1_time_case1", 200, -25, 25);
	spike_events_ewing1_time_case1->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing1_time_case1->GetYaxis()->SetTitle("Entries");	 

    TH1F *spike_events_ewing1_time_case2 = new TH1F("spike_events_ewing1_time_case2", "spike_events_ewing1_time_case2", 200, -25, 25);
	spike_events_ewing1_time_case2->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing1_time_case2->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing1_time_case3 = new TH1F("spike_events_ewing1_time_case3", "spike_events_ewing1_time_case3", 200, -25, 25);
	spike_events_ewing1_time_case3->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing1_time_case3->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing1_time_case4 = new TH1F("spike_events_ewing1_time_case4", "spike_events_ewing1_time_case4", 200, -25, 25);
	spike_events_ewing1_time_case4->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing1_time_case4->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing1_time_case5 = new TH1F("spike_events_ewing1_time_case5", "spike_events_ewing1_time_case5", 200, -25, 25);
	spike_events_ewing1_time_case5->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing1_time_case5->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing1_time = new TH1F("spike_events_ewing1_time", "spike_events_ewing1_time", 200, -25, 25);
	spike_events_ewing1_time->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing1_time->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing1_time_3ns = new TH1F("spike_events_ewing1_time_3ns", "spike_events_ewing1_time_3ns", 200, -25, 25);
	spike_events_ewing1_time_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing1_time_3ns->GetYaxis()->SetTitle("Entries");	


    //----------------------------------spike events with Ewing _time> 0.03----------------------------------
    TH1F *spike_events_ewing3_time_case1 = new TH1F("spike_events_ewing3_time_case1", "spike_events_ewing3_time_case1", 200, -25, 25);
	spike_events_ewing3_time_case1->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing3_time_case1->GetYaxis()->SetTitle("Entries");	 

    TH1F *spike_events_ewing3_time_case2 = new TH1F("spike_events_ewing3_time_case2", "spike_events_ewing3_time_case2", 200, -25, 25);
	spike_events_ewing3_time_case2->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing3_time_case2->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing3_time_case3 = new TH1F("spike_events_ewing3_time_case3", "spike_events_ewing3_time_case3", 200, -25, 25);
	spike_events_ewing3_time_case3->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing3_time_case3->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing3_time_case4 = new TH1F("spike_events_ewing3_time_case4", "spike_events_ewing3_time_case4", 200, -25, 25);
	spike_events_ewing3_time_case4->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing3_time_case4->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing3_time_case5 = new TH1F("spike_events_ewing3_time_case5", "spike_events_ewing3_time_case5", 200, -25, 25);
	spike_events_ewing3_time_case5->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing3_time_case5->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing3_time = new TH1F("spike_events_ewing3_time", "spike_events_ewing3_time", 200, -25, 25);
	spike_events_ewing3_time->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing3_time->GetYaxis()->SetTitle("Entries");	

    TH1F *spike_events_ewing3_time_3ns = new TH1F("spike_events_ewing3_time_3ns", "spike_events_ewing3_time_3ns", 200, -25, 25);
	spike_events_ewing3_time_3ns->GetXaxis()->SetTitle("SeedTime (ns)");
	spike_events_ewing3_time_3ns->GetYaxis()->SetTitle("Entries");	



    //----------------------------------BeamHalo template----------------------------------
	TH1F *BeamHalo_template_case1 = new TH1F("BeamHalo_template_case1", "BeamHalo_template_case1", 200, -25, 25);
	BeamHalo_template_case1->GetXaxis()->SetTitle("SeedTime (ns)");
	BeamHalo_template_case1->GetYaxis()->SetTitle("Entries");

	TH1F *BeamHalo_template_case2 = new TH1F("BeamHalo_template_case2", "BeamHalo_template_case2", 200, -25, 25);
	BeamHalo_template_case2->GetXaxis()->SetTitle("SeedTime (ns)");
	BeamHalo_template_case2->GetYaxis()->SetTitle("Entries");

	TH1F *BeamHalo_template_case3 = new TH1F("BeamHalo_template_case3", "BeamHalo_template_case3", 200, -25, 25);
	BeamHalo_template_case3->GetXaxis()->SetTitle("SeedTime (ns)");
	BeamHalo_template_case3->GetYaxis()->SetTitle("Entries");

	TH1F *BeamHalo_template_case4 = new TH1F("BeamHalo_template_case4", "BeamHalo_template_case4", 200, -25, 25);
	BeamHalo_template_case4->GetXaxis()->SetTitle("SeedTime (ns)");
	BeamHalo_template_case4->GetYaxis()->SetTitle("Entries");

	TH1F *BeamHalo_template_case5 = new TH1F("BeamHalo_template_case5", "BeamHalo_template_case5", 200, -25, 25);
	BeamHalo_template_case5->GetXaxis()->SetTitle("SeedTime (ns)");
	BeamHalo_template_case5->GetYaxis()->SetTitle("Entries");

	TH1F *BeamHalo_template = new TH1F("BeamHalo_template", "BeamHalo_template", 200, -25, 25);
	BeamHalo_template->GetXaxis()->SetTitle("SeedTime (ns)");
	BeamHalo_template->GetYaxis()->SetTitle("Entries");

    TH1F *BeamHalo_template_sieieb176 = new TH1F("BeamHalo_template_sieieb176", "BeamHalo template sieie > 0.0176", 200, -25, 25);
	BeamHalo_template_sieieb176->GetXaxis()->SetTitle("SeedTime (ns)");
	BeamHalo_template_sieieb176->GetYaxis()->SetTitle("Entries");

    TH1F *BeamHalo_template_sieieb206 = new TH1F("BeamHalo_template_sieieb206", "BeamHalo template sieie > 0.0206", 200, -25, 25);
	BeamHalo_template_sieieb206->GetXaxis()->SetTitle("SeedTime (ns)");
	BeamHalo_template_sieieb206->GetYaxis()->SetTitle("Entries");

    TH1F *BeamHalo_template_sieies176 = new TH1F("BeamHalo_template_sieies176", "BeamHalo template sieie < 0.0176", 200, -25, 25);
	BeamHalo_template_sieies176->GetXaxis()->SetTitle("SeedTime (ns)");
	BeamHalo_template_sieies176->GetYaxis()->SetTitle("Entries");

    TH1F *BeamHalo_template_sieies206 = new TH1F("BeamHalo_template_sieies206", "BeamHalo template sieie < 0.0206", 200, -25, 25);
	BeamHalo_template_sieies206->GetXaxis()->SetTitle("SeedTime (ns)");
	BeamHalo_template_sieies206->GetYaxis()->SetTitle("Entries");

    TH1F *sieie_BeamHalo = new TH1F("sieie_BeamHalo", "sieie_BeamHalo", 200, 25, 25);
    sieie_BeamHalo->GetXaxis()->SetTitle("Sieie");
    sieie_BeamHalo->GetYaxis()->SetTitle("Entries");


    
    cout << "Start running analysis on monopho17.C...." << endl;

	auto timenow_start = chrono::system_clock::to_time_t(chrono::system_clock::now());

	cout << "start time: " << ctime(&timenow_start) << " (CT)" << endl;


	ofstream file1, file2, file3, file4, file5, file6, file7, file8, file9, file10, file11, file12, file13, file14, file15, file16;

    file1.open("Analysis_log1.txt", ios::out);
    file1 << "Analysis started at " << ctime(&timenow_start) << " (CT)" << endl;
	file1 << "Analysis in process....." << endl;

    
    file2.open("Candidate_etawing_result.txt", ios::out);

    file3.open("promptZ_etawing_result.txt", ios::out);

    file4.open("spike_etawing_result.txt", ios::out);
    file5.open("spike_time_etawing_result.txt", ios::out);
    file6.open("pick_spike_timing_events.txt", ios::out);

    file7.open("err_etawing_cal.txt", ios::out);

    file8.open("candidate_ietamid.txt", ios::out);
    file9.open("promptZ_ietamid.txt", ios::out);
    file10.open("spike_ietamid.txt", ios::out);
    file11.open("extra_spike_events.txt", ios::out);
    file11 << "observed there are 6 events (case3, 4, 5 each one has 2) but not appeared in etawing dist" << endl;

    file12.open("event_display_info.txt", ios::out);

    file13.open("pick_candidate_failed_3.txt", ios::out);
    file14.open("pick_candidate_failed_1.txt", ios::out);
    file15.open("pick_spike_passed_3.txt", ios::out);
    file16.open("pick_spike_passed_1.txt", ios::out);


    
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        if ((jentry % 10000) == 0)
                // to print the number of processed entries
                std::cout << "Processed: " << jentry << std::endl;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        vector<Int_t> itpho;
		vector<Int_t> itootmatched;
		vector<Int_t> ootpho;
		vector<Int_t> ootitmatched;

        vector<Int_t> iiZpho;
		vector<Int_t> itootZmatched;
		vector<Int_t> ooZpho;
		vector<Int_t> ootitZmatched;

		//HLTPho == 1024 <== bitPho = 10 for HLT_Photon200_v (trigger info in ggNtuplizer_globalEvent.cc)
	
        /*
        The “eta-wing” variable is defined as the ratio W_eta = E1/E0. 
        where E0 is the energy of the seed crystal of the cluster and 
        E1 is the higher of the energies of two crystals adjacent to the seed in eta
        */

		//in-time collection
		for (int ipho = 0; ipho < nPho; ipho++)
		{
            if ((*phoEt)[ipho] > 230.0 && fabs((*phoEta)[ipho]) < 1.442 && (*phoR9)[ipho] > 0.8 && 
		(*phoHoverE)[ipho] < 0.02197 && ((HLTPho >> 10 & 1) == 1))//&& HLTPho < 2048
            {
                file12 << run << " " << event << " " << lumis << endl;
				itpho.push_back(ipho);
				itootmatched.push_back(-1);
                phoseedieta->Fill((*phoSeedIEta)[ipho]);
                phoseediphi->Fill((*phoSeedIPhi)[ipho]);
				//itootseparate.push_back(-1);
                
                //prompt Z in-time collection
                float PhoSep = 999;
				for (int iipho = ipho + 1; iipho < nPho; iipho++)
				{
					Float_t deltaR = DeltaR((*phoEta)[ipho], (*phoPhi)[ipho], 
								(*phoEta)[iipho], (*phoPhi)[iipho]);
					if (deltaR < PhoSep)
					{
						PhoSep = deltaR;
					}
                
                    //collect the second photon
                    if ((*phoEt)[iipho] > 10.0 && fabs((*phoEta)[iipho]) < 1.442 && 
			(*phoR9)[iipho] > 0.8 && (*phoHoverE)[iipho] < 0.02197 && PhoSep > 0.2)
                    {
                        iiZpho.push_back(iipho);
                        itootZmatched.push_back(-1);
                    }
                }
			}
		}


		//out-of-time collection
		for (int opho = 0; opho < onPho; opho++)
		{
            if ((*ophoEt)[opho] > 230.0 && fabs((*ophoEta)[opho]) < 1.442 && (*ophoR9)[opho] > 0.8 && 
		(*ophoHoverE)[opho] < 0.02197 && ((HLTPho >> 10 & 1) == 1))//&& HLTPho < 2048
            {
				ootpho.push_back(opho);
				ootitmatched.push_back(-1);

                ophoseedieta->Fill((*ophoSeedIEta)[opho]);
                ophoseediphi->Fill((*ophoSeedIPhi)[opho]);
				//ootitseparate.push_back(-1);

                //prompt Z out-of-time collection
                float PhoSep = 999;
				for (int oopho = opho + 1; oopho < nPho; oopho++)
				{
					Float_t deltaR = DeltaR((*ophoEta)[opho], (*ophoPhi)[opho], 
								(*ophoEta)[oopho], (*ophoPhi)[oopho]);
					if (deltaR < PhoSep)
					{
						PhoSep = deltaR;
					}

                    //collect the second photon
                    if ((*ophoEt)[oopho] > 10.0 && fabs((*ophoEta)[oopho]) < 1.442 && 
			(*ophoR9)[oopho] > 0.8 && (*ophoHoverE)[oopho] < 0.02197)
                    {
                        ooZpho.push_back(oopho);
                        ootitZmatched.push_back(-1);
                    }
                }
			}
		}

		
		//for getting the minimum deltaR
		//for (unsigned i = 0; i < itpho.size(); i++)
		for (unsigned i = 0; i < itpho.size(); ++i)
		{
			float PhoSep = 999;
			Int_t itootmatchedidx = -1; //defualt that they are not matched
			//Int_t itootseparateidx = 1;

			//for (unsigned j = 0; j < ootpho.size(); j++)
			for (unsigned j = 0; j < ootpho.size(); ++j)
			{
				Float_t deltaR = DeltaR((*phoEta)[i], (*phoPhi)[i], (*ophoEta)[j], (*ophoPhi)[j]);
				if (deltaR < PhoSep)
				{
					PhoSep = deltaR;
				}
				
				if (PhoSep < 0.2)
				{
					itootmatchedidx = j;
					itootmatched[i] = itootmatchedidx;
					ootitmatched[itootmatchedidx] = i;
				}
			}
		}

        //second photon matched for Z prompt
        for (unsigned i = 0; i < iiZpho.size(); ++i)
		{
			float PhoSep = 999;
			Int_t itootZmatchedidx = -1; //defualt that they are not matched
			//Int_t itootseparateidx = 1;

			//for (unsigned j = 0; j < ootpho.size(); j++)
			for (unsigned j = 0; j < ooZpho.size(); ++j)
			{
				Float_t deltaR = DeltaR((*phoEta)[i], (*phoPhi)[i], (*ophoEta)[j], (*ophoPhi)[j]);
				if (deltaR < PhoSep)
				{
					PhoSep = deltaR;
				}
				
				if (PhoSep < 0.2)
				{
					itootZmatchedidx = j;
					itootZmatched[i] = itootZmatchedidx;
					ootitZmatched[itootZmatchedidx] = i;
				}
			}
		}

        Int_t mitIDXSEL = -1; //set index for photon, initial value could be something non-physical, negative something.
		Int_t mootIDXSEL = -1;
        Int_t mootZIDXSEL = -1;
        Int_t mitZIDXSEL = -1;
        

		Bool_t misOOT = kFALSE; //if it's OOT photon, select the OOT, can be true or false but false makes more sense.
		Int_t mIDXSEL = -1;
        Int_t mZIDXSEL = -1;
    

		/*
		basic photon ID selection: (medium)
		(*phoEt)[iPho] > 230.0
		fabs((*phoEta)[iPho]) < 1.442
		(*phoR9)[iPho] > 0.8
		(*phoHoverE) < 0.02197


		Candidate events:
		MET cut > 210 GeV
		has no pixel seed -> photon
		sieie and sipip > 0.01 for spike cleaning
		sieie < 0.01015 for discriminating jet
		MIP TotE < 4.9 GeV for Beam Halo reducing

		Beam Halo Template:
		basic photon ID selection
		MET cut > 210 GeV
		has no pixel seed -> photon
		sieie and sipip > 0.01
		MIP TotE > 4.9 GeV 

		Spike Template:
		basic photon ID selection
		MET cut > 210 GeV
		has no pixel seed -> photon
		sieie or sipip < 0.01
		MIP TotE < 4.9 GeV for Beam Halo reducing

		Prompt W Template:
		basic photon ID selection
		MET cut > 210 GeV
		has pixel seed -> electron
		sieie and sipip > 0.01
		PhoSep > 0.2
		Require 2nd photon ET > 10 GeV
		Invariant mass ~ 80 GeV (70~90 window)
		*/

        Int_t CellIEtaEB_right = 0;
        Float_t AllCellsE_EB_right = 0;
        Int_t CellIEtaEB_left = 0;
        Float_t AllCellsE_EB_left = 0;
        Float_t maxCellsE = 0;

        Int_t rightcellsIDX = -1;
        Int_t leftcellsIDX = -1;

        Int_t rightcellsIDX_only = -1;
        Int_t leftcellsIDX_only = -1;

        vector<Int_t> cellsEB_right;
        vector<Int_t> cellsEB_left;
        



        //---------------------------case 1, only it, no oot---------------------------
		if (itpho.size() != 0 && ootpho.size() == 0 && itootmatched.size() != 0)
		{
			mIDXSEL = itpho[0];

            //candidate event
			if (pfMET > 210 &&
				(*phohasPixelSeed)[mIDXSEL] == 0 &&
				(*phoMIPTotEnergy)[mIDXSEL] < 4.9 &&
				(*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 &&
				(*phoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001 &&
				(*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.01015 &&
				metFilters == 0)
			{
                candidate_events_case1->Fill((*phoSeedTime)[mIDXSEL]);
				candidate_events->Fill((*phoSeedTime)[mIDXSEL]);

                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                {
                    candidate_events_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                }

                if ((*phoSeedTime)[mIDXSEL] < -12.5)
                {
                    spikeincandidate_events_case1->Fill((*phoSeedTime)[mIDXSEL]);
				    spikeincandidate_events->Fill((*phoSeedTime)[mIDXSEL]);
                }

                //cand1 searching adjencent cells
                for (int icell = 0; icell < nAllCellsEB; icell++)
                {                    
                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mIDXSEL] + 1)
                    {
                        cellsEB_right.push_back(icell); 

                        //if there's any err
                        if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mIDXSEL]) == 86 || (*phoSeedIEta)[mIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                        {
                            file7 << "------------------cand1_err_right---------------------" << "\n"
                            << "phoSeediphi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mIDXSEL] << "\n"
                            << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                        }
                    }

                    //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && (*phoSeedIEta)[mIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                    {
                        cellsEB_right.push_back(icell); 

                        file8 << "-------------Cand1_ietamid_right-------------" << "\n"
                        << "phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;
                    }

                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mIDXSEL] - 1)
                    {                    
                        cellsEB_left.push_back(icell);

                        //if there's any err
                        if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mIDXSEL]) == 86 || (*phoSeedIEta)[mIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                        {
                            file7 << "------------------cand1_err_left---------------------" << "\n"
                            << "phoSeediphi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mIDXSEL] << "\n"
                            << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                        }
                    }            

                    //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && (*phoSeedIEta)[mIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                    {
                        cellsEB_left.push_back(icell); 

                        file8 << "-------------Cand1_ietamid_left-------------" << "\n"
                        << "phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;
                    }        
                }

                //Now I have to take if there's only one and both
                //If there's only right
                if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                {
                    rightcellsIDX_only = cellsEB_right[0];
                    maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                    file2 << "-------------Cand1_rightonly-------------" << "\n"
                    << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                    << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                    << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                    << run << " " << lumis << " " << event << "\n"
                    << endl;

                    candidate_eta_wing_rightonly1->Fill(Ewing);
                    candidate_eta_wing1->Fill(Ewing);
                    candidate_eta_wing->Fill(Ewing);

                    if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                    {
                        candidate_eta_wing_prompt->Fill(Ewing);
                    }

                    if (Ewing > 0.03)
                    {
                        candidate_events_ewing3_case1->Fill((*phoSeedTime)[mIDXSEL]);
                        candidate_events_ewing3->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing3_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.01)
                    {
                        candidate_events_ewing1_case1->Fill((*phoSeedTime)[mIDXSEL]);
                        candidate_events_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing1_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }


                    //---------------------------more ewing cuts-------------------------------



                    if (Ewing > 0.02)
                    {
                        candidate_events_ewing2->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing2_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.025)
                    {
                        candidate_events_ewing25->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing25_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.031)
                    {
                        candidate_events_ewing31->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing31_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.032)
                    {
                        candidate_events_ewing32->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing32_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.033)
                    {
                        candidate_events_ewing33->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing33_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.034)
                    {
                        candidate_events_ewing34->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing34_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.035)
                    {
                        candidate_events_ewing35->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing35_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.036)
                    {
                        candidate_events_ewing36->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing36_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.037)
                    {
                        candidate_events_ewing37->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing37_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.038)
                    {
                        candidate_events_ewing38->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing38_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.039)
                    {
                        candidate_events_ewing39->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing39_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }
                    

                    if (Ewing > 0.040)
                    {
                        candidate_events_ewing4->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing4_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.041)
                    {
                        candidate_events_ewing41->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing41_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.042)
                    {
                        candidate_events_ewing42->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing42_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.043)
                    {
                        candidate_events_ewing43->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing43_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.044)
                    {
                        candidate_events_ewing44->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing44_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.045)
                    {
                        candidate_events_ewing45->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing45_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.046)
                    {
                        candidate_events_ewing46->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing46_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.047)
                    {
                        candidate_events_ewing47->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing47_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.048)
                    {
                        candidate_events_ewing48->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing48_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.049)
                    {
                        candidate_events_ewing49->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing49_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }



                    if (Ewing > 0.05)
                    {
                        candidate_events_ewing5->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing5_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.07)
                    {
                        candidate_events_ewing7->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing7_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.09)
                    {
                        candidate_events_ewing9->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing9_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.10)
                    {
                        candidate_events_ewing10->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing10_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.12)
                    {
                        candidate_events_ewing12->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing12_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.14)
                    {
                        candidate_events_ewing14->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing14_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.15)
                    {
                        candidate_events_ewing15->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing15_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }




                    if ((*phoSeedTime)[mIDXSEL] < -12.5)
                    {
                        spikeincandidate_eta_wing_rightonly1->Fill(Ewing);
                        spikeincandidate_eta_wing1->Fill(Ewing);
                        spikeincandidate_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spikeincandidate_events_ewing3_case1->Fill((*phoSeedTime)[mIDXSEL]);
                            spikeincandidate_events_ewing3->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.01)
                        {
                            spikeincandidate_events_ewing1_case1->Fill((*phoSeedTime)[mIDXSEL]);
                            spikeincandidate_events_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        //---------------------------more ewing cuts-------------------------------

                        if (Ewing > 0.02)
                        {
                            spikeincandidate_events_ewing2->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.025)
                        {
                            spikeincandidate_events_ewing25->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.031)
                        {
                            spikeincandidate_events_ewing31->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.032)
                        {
                            spikeincandidate_events_ewing32->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.033)
                        {
                            spikeincandidate_events_ewing33->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.034)
                        {
                            spikeincandidate_events_ewing34->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.035)
                        {
                            spikeincandidate_events_ewing35->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.036)
                        {
                            spikeincandidate_events_ewing36->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.037)
                        {
                            spikeincandidate_events_ewing37->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.038)
                        {
                            spikeincandidate_events_ewing38->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.039)
                        {
                            spikeincandidate_events_ewing39->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                        

                        if (Ewing > 0.040)
                        {
                            spikeincandidate_events_ewing4->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.041)
                        {
                            spikeincandidate_events_ewing41->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.042)
                        {
                            spikeincandidate_events_ewing42->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.043)
                        {
                            spikeincandidate_events_ewing43->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.044)
                        {
                            spikeincandidate_events_ewing44->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.045)
                        {
                            spikeincandidate_events_ewing45->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.046)
                        {
                            spikeincandidate_events_ewing46->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.047)
                        {
                            spikeincandidate_events_ewing47->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.048)
                        {
                            spikeincandidate_events_ewing48->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.049)
                        {
                            spikeincandidate_events_ewing49->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.05)
                        {
                            spikeincandidate_events_ewing5->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.07)
                        {
                            spikeincandidate_events_ewing7->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.09)
                        {
                            spikeincandidate_events_ewing9->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.10)
                        {
                            spikeincandidate_events_ewing10->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.12)
                        {
                            spikeincandidate_events_ewing12->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.14)
                        {
                            spikeincandidate_events_ewing14->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.15)
                        {
                            spikeincandidate_events_ewing15->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    //to pick the can failed the cut
                    if (Ewing < 0.03)
                    {
                        file13 << "--------cand case 1--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }

                    if (Ewing < 0.01)
                    {
                        file14 << "--------cand case 1--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }

                }

                //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                    file2 << "-------------Cand1_leftonly-------------" << "\n"
                    << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                    << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                    << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                    << run << " " << lumis << " " << event << "\n"
                    << endl;

                    candidate_eta_wing_leftonly1->Fill(Ewing);
                    candidate_eta_wing1->Fill(Ewing);
                    candidate_eta_wing->Fill(Ewing);

                    if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                    {
                        candidate_eta_wing_prompt->Fill(Ewing);
                    }

                    if (Ewing > 0.03)
                    {
                        candidate_events_ewing3_case1->Fill((*phoSeedTime)[mIDXSEL]);
                        candidate_events_ewing3->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing3_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.01)
                    {
                        candidate_events_ewing1_case1->Fill((*phoSeedTime)[mIDXSEL]);
                        candidate_events_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing1_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    //---------------------------more ewing cuts-------------------------------

                    if (Ewing > 0.02)
                    {
                        candidate_events_ewing2->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing2_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.025)
                    {
                        candidate_events_ewing25->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing25_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.031)
                    {
                        candidate_events_ewing31->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing31_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.032)
                    {
                        candidate_events_ewing32->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing32_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.033)
                    {
                        candidate_events_ewing33->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing33_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.034)
                    {
                        candidate_events_ewing34->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing34_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.035)
                    {
                        candidate_events_ewing35->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing35_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.036)
                    {
                        candidate_events_ewing36->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing36_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.037)
                    {
                        candidate_events_ewing37->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing37_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.038)
                    {
                        candidate_events_ewing38->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing38_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.039)
                    {
                        candidate_events_ewing39->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing39_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }
                    

                    if (Ewing > 0.040)
                    {
                        candidate_events_ewing4->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing4_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.041)
                    {
                        candidate_events_ewing41->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing41_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.042)
                    {
                        candidate_events_ewing42->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing42_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.043)
                    {
                        candidate_events_ewing43->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing43_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.044)
                    {
                        candidate_events_ewing44->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing44_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.045)
                    {
                        candidate_events_ewing45->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing45_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.046)
                    {
                        candidate_events_ewing46->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing46_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.047)
                    {
                        candidate_events_ewing47->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing47_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.048)
                    {
                        candidate_events_ewing48->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing48_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.049)
                    {
                        candidate_events_ewing49->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing49_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.05)
                    {
                        candidate_events_ewing5->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing5_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.07)
                    {
                        candidate_events_ewing7->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing7_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.09)
                    {
                        candidate_events_ewing9->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing9_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.10)
                    {
                        candidate_events_ewing10->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing10_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.12)
                    {
                        candidate_events_ewing12->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing12_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.14)
                    {
                        candidate_events_ewing14->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing14_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.15)
                    {
                        candidate_events_ewing15->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing15_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if ((*phoSeedTime)[mIDXSEL] < -12.5)
                    {
                        spikeincandidate_eta_wing_leftonly1->Fill(Ewing);
                        spikeincandidate_eta_wing1->Fill(Ewing);
                        spikeincandidate_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spikeincandidate_events_ewing3_case1->Fill((*phoSeedTime)[mIDXSEL]);
                            spikeincandidate_events_ewing3->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.01)
                        {
                            spikeincandidate_events_ewing1_case1->Fill((*phoSeedTime)[mIDXSEL]);
                            spikeincandidate_events_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        //---------------------------more ewing cuts-------------------------------

                        if (Ewing > 0.02)
                        {
                            spikeincandidate_events_ewing2->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.025)
                        {
                            spikeincandidate_events_ewing25->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.031)
                        {
                            spikeincandidate_events_ewing31->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.032)
                        {
                            spikeincandidate_events_ewing32->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.033)
                        {
                            spikeincandidate_events_ewing33->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.034)
                        {
                            spikeincandidate_events_ewing34->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.035)
                        {
                            spikeincandidate_events_ewing35->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.036)
                        {
                            spikeincandidate_events_ewing36->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.037)
                        {
                            spikeincandidate_events_ewing37->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.038)
                        {
                            spikeincandidate_events_ewing38->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.039)
                        {
                            spikeincandidate_events_ewing39->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                        

                        if (Ewing > 0.040)
                        {
                            spikeincandidate_events_ewing4->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.041)
                        {
                            spikeincandidate_events_ewing41->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.042)
                        {
                            spikeincandidate_events_ewing42->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.043)
                        {
                            spikeincandidate_events_ewing43->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.044)
                        {
                            spikeincandidate_events_ewing44->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.045)
                        {
                            spikeincandidate_events_ewing45->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.046)
                        {
                            spikeincandidate_events_ewing46->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.047)
                        {
                            spikeincandidate_events_ewing47->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.048)
                        {
                            spikeincandidate_events_ewing48->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.049)
                        {
                            spikeincandidate_events_ewing49->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.05)
                        {
                            spikeincandidate_events_ewing5->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.07)
                        {
                            spikeincandidate_events_ewing7->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.09)
                        {
                            spikeincandidate_events_ewing9->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.10)
                        {
                            spikeincandidate_events_ewing10->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.12)
                        {
                            spikeincandidate_events_ewing12->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.14)
                        {
                            spikeincandidate_events_ewing14->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.15)
                        {
                            spikeincandidate_events_ewing15->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing < 0.03)
                    {
                        file13 << "--------cand case 1--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }

                    if (Ewing < 0.01)
                    {
                        file14 << "--------cand case 1--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }
                }

                //if there are both
                if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                {
                    rightcellsIDX = cellsEB_right[0];
                    leftcellsIDX = cellsEB_left[0];

                    if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[rightcellsIDX];
                    }
                    else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[leftcellsIDX];
                    }

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                    file2 << "-------------Cand1_both-------------" << "\n"
                    << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                    << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                    << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                    << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                    << run << " " << lumis << " " << event << "\n"
                    << endl;

                    candidate_eta_wing_both1->Fill(Ewing);
                    candidate_eta_wing1->Fill(Ewing);
                    candidate_eta_wing->Fill(Ewing);

                    if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                    {
                        candidate_eta_wing_prompt->Fill(Ewing);
                    }

                    if (Ewing > 0.03)
                    {
                        candidate_events_ewing3_case1->Fill((*phoSeedTime)[mIDXSEL]);
                        candidate_events_ewing3->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing3_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.01)
                    {
                        candidate_events_ewing1_case1->Fill((*phoSeedTime)[mIDXSEL]);
                        candidate_events_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing1_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    //---------------------------more ewing cuts-------------------------------

                    if (Ewing > 0.02)
                    {
                        candidate_events_ewing2->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing2_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.025)
                    {
                        candidate_events_ewing25->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing25_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.031)
                    {
                        candidate_events_ewing31->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing31_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.032)
                    {
                        candidate_events_ewing32->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing32_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.033)
                    {
                        candidate_events_ewing33->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing33_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.034)
                    {
                        candidate_events_ewing34->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing34_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.035)
                    {
                        candidate_events_ewing35->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing35_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.036)
                    {
                        candidate_events_ewing36->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing36_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.037)
                    {
                        candidate_events_ewing37->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing37_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.038)
                    {
                        candidate_events_ewing38->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing38_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.039)
                    {
                        candidate_events_ewing39->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing39_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }
                    

                    if (Ewing > 0.040)
                    {
                        candidate_events_ewing4->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing4_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.041)
                    {
                        candidate_events_ewing41->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing41_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.042)
                    {
                        candidate_events_ewing42->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing42_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.043)
                    {
                        candidate_events_ewing43->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing43_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.044)
                    {
                        candidate_events_ewing44->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing44_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.045)
                    {
                        candidate_events_ewing45->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing45_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.046)
                    {
                        candidate_events_ewing46->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing46_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.047)
                    {
                        candidate_events_ewing47->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing47_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.048)
                    {
                        candidate_events_ewing48->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing48_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.049)
                    {
                        candidate_events_ewing49->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing49_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.05)
                    {
                        candidate_events_ewing5->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing5_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.07)
                    {
                        candidate_events_ewing7->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing7_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.09)
                    {
                        candidate_events_ewing9->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing9_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.10)
                    {
                        candidate_events_ewing10->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing10_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.12)
                    {
                        candidate_events_ewing12->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing12_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.14)
                    {
                        candidate_events_ewing14->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing14_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.15)
                    {
                        candidate_events_ewing15->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing15_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }


                    if ((*phoSeedTime)[mIDXSEL] < -12.5)
                    {
                        spikeincandidate_eta_wing_both1->Fill(Ewing);
                        spikeincandidate_eta_wing1->Fill(Ewing);
                        spikeincandidate_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spikeincandidate_events_ewing3_case1->Fill((*phoSeedTime)[mIDXSEL]);
                            spikeincandidate_events_ewing3->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.01)
                        {
                            spikeincandidate_events_ewing1_case1->Fill((*phoSeedTime)[mIDXSEL]);
                            spikeincandidate_events_ewing1->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        //---------------------------more ewing cuts-------------------------------

                        if (Ewing > 0.02)
                        {
                            spikeincandidate_events_ewing2->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.025)
                        {
                            spikeincandidate_events_ewing25->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.031)
                        {
                            spikeincandidate_events_ewing31->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.032)
                        {
                            spikeincandidate_events_ewing32->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.033)
                        {
                            spikeincandidate_events_ewing33->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.034)
                        {
                            spikeincandidate_events_ewing34->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.035)
                        {
                            spikeincandidate_events_ewing35->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.036)
                        {
                            spikeincandidate_events_ewing36->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.037)
                        {
                            spikeincandidate_events_ewing37->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.038)
                        {
                            spikeincandidate_events_ewing38->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.039)
                        {
                            spikeincandidate_events_ewing39->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                        

                        if (Ewing > 0.040)
                        {
                            spikeincandidate_events_ewing4->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.041)
                        {
                            spikeincandidate_events_ewing41->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.042)
                        {
                            spikeincandidate_events_ewing42->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.043)
                        {
                            spikeincandidate_events_ewing43->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.044)
                        {
                            spikeincandidate_events_ewing44->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.045)
                        {
                            spikeincandidate_events_ewing45->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.046)
                        {
                            spikeincandidate_events_ewing46->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.047)
                        {
                            spikeincandidate_events_ewing47->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.048)
                        {
                            spikeincandidate_events_ewing48->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.049)
                        {
                            spikeincandidate_events_ewing49->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.05)
                        {
                            spikeincandidate_events_ewing5->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.07)
                        {
                            spikeincandidate_events_ewing7->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.09)
                        {
                            spikeincandidate_events_ewing9->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.10)
                        {
                            spikeincandidate_events_ewing10->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.12)
                        {
                            spikeincandidate_events_ewing12->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.14)
                        {
                            spikeincandidate_events_ewing14->Fill((*phoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.15)
                        {
                            spikeincandidate_events_ewing15->Fill((*phoSeedTime)[mIDXSEL]);
                        }
                    }


                    if (Ewing < 0.03)
                    {
                        file13 << "--------cand case 1--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }

                    if (Ewing < 0.01)
                    {
                        file14 << "--------cand case 1--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }
                }
			}


            //promptZ selection
            if (iiZpho.size() != 0 && ooZpho.size() == 0 && itootZmatched.size() != 0)
            {
                mZIDXSEL = iiZpho[0];
                //Prompt Z template:
                if ((*phohasPixelSeed)[mIDXSEL] == 1 &&
                    (*phoMIPTotEnergy)[mIDXSEL] < 4.9 &&
                    (*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 &&
                    (*phoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001 &&
                    (*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.01015 &&
                    metFilters == 0 &&
                    (*phohasPixelSeed)[mZIDXSEL] == 1 &&
                    (*phoMIPTotEnergy)[mZIDXSEL] < 4.9 &&
                    (*phoSigmaIEtaIEtaFull5x5)[mZIDXSEL] > 0.001 &&
                    (*phoSigmaIPhiIPhiFull5x5)[mZIDXSEL] > 0.001 &&
                    (*phoSigmaIEtaIEtaFull5x5)[mZIDXSEL] < 0.01015)
                {
                    Float_t InvM = InvariMass((*phoEt)[mIDXSEL], (*phoEt)[mZIDXSEL], (*phoPhi)[mIDXSEL], (*phoPhi)[mZIDXSEL], (*phoEta)[mIDXSEL], (*phoEta)[mZIDXSEL]);
                    InvMass_Z->Fill(InvM);

                    if (InvM > 85 && InvM < 100)
                    {
                        promptZ_events_case1->Fill((*phoSeedTime)[mIDXSEL]);
                        promptZ_events->Fill((*phoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            promptZ_events_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                        }
  
                        //promptZ1 searching adjencent cells
                        for (int icell = 0; icell < nAllCellsEB; icell++)
                        {                        
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mIDXSEL] + 1)
                            {
                                cellsEB_right.push_back(icell); 

                                //if there's any err
                                if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mIDXSEL]) == 86 || (*phoSeedIEta)[mIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                                {
                                    file7 << "------------------promptZ_err_right---------------------" << "\n"
                                    << "phoSeediphi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mIDXSEL] << "\n"
                                    << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                                }
                            }

                            //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && (*phoSeedIEta)[mIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                            {
                                cellsEB_right.push_back(icell); 

                                file9 << "-------------promptZ_ietamid_right-------------" << "\n"
                                << "phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                                << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                                << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                                << run << " " << lumis << " " << event << "\n"
                                << endl;
                            }

                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mIDXSEL] - 1)
                            {                    
                                cellsEB_left.push_back(icell);

                                //if there's any err
                                if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mIDXSEL]) == 86 || (*phoSeedIEta)[mIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                                {
                                    file7 << "------------------promptZ_err_left---------------------" << "\n"
                                    << "phoSeediphi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mIDXSEL] << "\n"
                                    << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                                }
                            }            

                            //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && (*phoSeedIEta)[mIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                            {
                                cellsEB_left.push_back(icell); 

                                file9 << "-------------promptZ_ietamid_left-------------" << "\n"
                                << "phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                                << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                                << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                                << run << " " << lumis << " " << event << "\n"
                                << endl;
                            }  
                        }

                        //Now I have to take if there's only one and both
                        //If there's only right
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                        {
                            rightcellsIDX_only = cellsEB_right[0];
                            maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                            file3 << "-------------promptZ1_rightonly-------------" << "\n"
                            << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;

                            promptZ_eta_wing_rightonly1->Fill(Ewing);
                            promptZ_eta_wing1->Fill(Ewing);
                            promptZ_eta_wing->Fill(Ewing);

                            if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                promptZ_eta_wing_prompt->Fill(Ewing);
                            }

                            if (Ewing > 0.03)
                            {
                                promptZ_events_ewing3_case1->Fill((*phoSeedTime)[mIDXSEL]);
                                promptZ_events_ewing3->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing3_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.01)
                            {
                                promptZ_events_ewing1_case1->Fill((*phoSeedTime)[mIDXSEL]);
                                promptZ_events_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing1_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            //---------------------------more ewing cuts-------------------------------

                            if (Ewing > 0.02)
                            {
                                promptZ_events_ewing2->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing2_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.025)
                            {
                                promptZ_events_ewing25->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing25_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.031)
                            {
                                promptZ_events_ewing31->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing31_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.032)
                            {
                                promptZ_events_ewing32->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing32_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.033)
                            {
                                promptZ_events_ewing33->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing33_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.034)
                            {
                                promptZ_events_ewing34->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing34_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.035)
                            {
                                promptZ_events_ewing35->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing35_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.036)
                            {
                                promptZ_events_ewing36->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing36_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.037)
                            {
                                promptZ_events_ewing37->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing37_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.038)
                            {
                                promptZ_events_ewing38->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing38_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.039)
                            {
                                promptZ_events_ewing39->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing39_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }
                            

                            if (Ewing > 0.040)
                            {
                                promptZ_events_ewing4->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing4_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.041)
                            {
                                promptZ_events_ewing41->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing41_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.042)
                            {
                                promptZ_events_ewing42->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing42_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.043)
                            {
                                promptZ_events_ewing43->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing43_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.044)
                            {
                                promptZ_events_ewing44->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing44_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.045)
                            {
                                promptZ_events_ewing45->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing45_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.046)
                            {
                                promptZ_events_ewing46->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing46_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.047)
                            {
                                promptZ_events_ewing47->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing47_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.048)
                            {
                                promptZ_events_ewing48->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing48_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.049)
                            {
                                promptZ_events_ewing49->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing49_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.05)
                            {
                                promptZ_events_ewing5->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing5_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.07)
                            {
                                promptZ_events_ewing7->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing7_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.09)
                            {
                                promptZ_events_ewing9->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing9_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.10)
                            {
                                promptZ_events_ewing10->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing10_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.12)
                            {
                                promptZ_events_ewing12->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing12_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.14)
                            {
                                promptZ_events_ewing14->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing14_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.15)
                            {
                                promptZ_events_ewing15->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing15_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }
                        }

                        //if there is only left
                        if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                        {
                            leftcellsIDX_only = cellsEB_left[0];
                            maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                            file3 << "-------------promptZ1_leftonly-------------" << "\n"
                            << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;

                            promptZ_eta_wing_leftonly1->Fill(Ewing);
                            promptZ_eta_wing1->Fill(Ewing);
                            promptZ_eta_wing->Fill(Ewing);

                            if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                promptZ_eta_wing_prompt->Fill(Ewing);
                            }

                            if (Ewing > 0.03)
                            {
                                promptZ_events_ewing3_case1->Fill((*phoSeedTime)[mIDXSEL]);
                                promptZ_events_ewing3->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing3_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.01)
                            {
                                promptZ_events_ewing1_case1->Fill((*phoSeedTime)[mIDXSEL]);
                                promptZ_events_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing1_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            //---------------------------more ewing cuts-------------------------------

                            if (Ewing > 0.02)
                            {
                                promptZ_events_ewing2->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing2_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.025)
                            {
                                promptZ_events_ewing25->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing25_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.031)
                            {
                                promptZ_events_ewing31->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing31_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.032)
                            {
                                promptZ_events_ewing32->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing32_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.033)
                            {
                                promptZ_events_ewing33->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing33_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.034)
                            {
                                promptZ_events_ewing34->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing34_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.035)
                            {
                                promptZ_events_ewing35->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing35_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.036)
                            {
                                promptZ_events_ewing36->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing36_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.037)
                            {
                                promptZ_events_ewing37->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing37_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.038)
                            {
                                promptZ_events_ewing38->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing38_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.039)
                            {
                                promptZ_events_ewing39->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing39_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }
                            

                            if (Ewing > 0.040)
                            {
                                promptZ_events_ewing4->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing4_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.041)
                            {
                                promptZ_events_ewing41->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing41_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.042)
                            {
                                promptZ_events_ewing42->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing42_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.043)
                            {
                                promptZ_events_ewing43->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing43_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.044)
                            {
                                promptZ_events_ewing44->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing44_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.045)
                            {
                                promptZ_events_ewing45->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing45_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.046)
                            {
                                promptZ_events_ewing46->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing46_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.047)
                            {
                                promptZ_events_ewing47->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing47_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.048)
                            {
                                promptZ_events_ewing48->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing48_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.049)
                            {
                                promptZ_events_ewing49->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing49_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.05)
                            {
                                promptZ_events_ewing5->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing5_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.07)
                            {
                                promptZ_events_ewing7->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing7_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.09)
                            {
                                promptZ_events_ewing9->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing9_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.10)
                            {
                                promptZ_events_ewing10->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing10_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.12)
                            {
                                promptZ_events_ewing12->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing12_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.14)
                            {
                                promptZ_events_ewing14->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing14_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.15)
                            {
                                promptZ_events_ewing15->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing15_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }
                        }

                        //if there are both
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                        {
                            rightcellsIDX = cellsEB_right[0];
                            leftcellsIDX = cellsEB_left[0];

                            if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[rightcellsIDX];
                            }
                            else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[leftcellsIDX];
                            }

                            Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                            file3 << "-------------promptZ1_both-------------" << "\n"
                            << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;

                            promptZ_eta_wing_both1->Fill(Ewing);
                            promptZ_eta_wing1->Fill(Ewing);
                            promptZ_eta_wing->Fill(Ewing);

                            if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                promptZ_eta_wing_prompt->Fill(Ewing);
                            }

                            if (Ewing > 0.03)
                            {
                                promptZ_events_ewing3_case1->Fill((*phoSeedTime)[mIDXSEL]);
                                promptZ_events_ewing3->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing3_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.01)
                            {
                                promptZ_events_ewing1_case1->Fill((*phoSeedTime)[mIDXSEL]);
                                promptZ_events_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing1_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            //---------------------------more ewing cuts-------------------------------

                            if (Ewing > 0.02)
                            {
                                promptZ_events_ewing2->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing2_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.025)
                            {
                                promptZ_events_ewing25->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing25_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.031)
                            {
                                promptZ_events_ewing31->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing31_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.032)
                            {
                                promptZ_events_ewing32->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing32_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.033)
                            {
                                promptZ_events_ewing33->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing33_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.034)
                            {
                                promptZ_events_ewing34->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing34_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.035)
                            {
                                promptZ_events_ewing35->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing35_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.036)
                            {
                                promptZ_events_ewing36->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing36_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.037)
                            {
                                promptZ_events_ewing37->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing37_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.038)
                            {
                                promptZ_events_ewing38->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing38_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.039)
                            {
                                promptZ_events_ewing39->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing39_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }
                            

                            if (Ewing > 0.040)
                            {
                                promptZ_events_ewing4->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing4_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.041)
                            {
                                promptZ_events_ewing41->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing41_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.042)
                            {
                                promptZ_events_ewing42->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing42_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.043)
                            {
                                promptZ_events_ewing43->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing43_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.044)
                            {
                                promptZ_events_ewing44->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing44_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.045)
                            {
                                promptZ_events_ewing45->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing45_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.046)
                            {
                                promptZ_events_ewing46->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing46_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.047)
                            {
                                promptZ_events_ewing47->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing47_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.048)
                            {
                                promptZ_events_ewing48->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing48_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.049)
                            {
                                promptZ_events_ewing49->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing49_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.05)
                            {
                                promptZ_events_ewing5->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing5_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.07)
                            {
                                promptZ_events_ewing7->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing7_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.09)
                            {
                                promptZ_events_ewing9->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing9_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.10)
                            {
                                promptZ_events_ewing10->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing10_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.12)
                            {
                                promptZ_events_ewing12->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing12_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.14)
                            {
                                promptZ_events_ewing14->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing14_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.15)
                            {
                                promptZ_events_ewing15->Fill((*phoSeedTime)[mIDXSEL]);

                                if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing15_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                                }
                            }
                        } 
                    }
                }
            }

            //BeamHalo template
			if (pfMET > 210 &&
				(*phohasPixelSeed)[mIDXSEL] == 0 &&
				(*phoMIPTotEnergy)[mIDXSEL] > 4.9 &&
				(*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 &&
				(*phoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001 &&
				metFilters == 8)
			{
				BeamHalo_template_case1->Fill((*phoSeedTime)[mIDXSEL]);
				BeamHalo_template->Fill((*phoSeedTime)[mIDXSEL]);

                sieie_BeamHalo->Fill((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL]);


                if ((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.0176)
                {
                    BeamHalo_template_sieieb176->Fill((*phoSeedTime)[mIDXSEL]);
                }
                
                if ((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.0206)
                {
                    BeamHalo_template_sieieb206->Fill((*phoSeedTime)[mIDXSEL]);
                }

                if ((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.0176)
                {
                    BeamHalo_template_sieies176->Fill((*phoSeedTime)[mIDXSEL]);
                }

                if ((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.0206)
                {
                    BeamHalo_template_sieies206->Fill((*phoSeedTime)[mIDXSEL]);
                }
            }

            //Spikes template
            if (pfMET > 210 &&
                (*phohasPixelSeed)[mIDXSEL] == 0 &&
                (*phoMIPTotEnergy)[mIDXSEL] < 4.9 &&
                ((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.001 ||
                (*phoSigmaIPhiIPhiFull5x5)[mIDXSEL] < 0.001) &&
                metFilters == 0)
            {
                //original showershap cut
                if ((*phoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.01015)
                {
                    spike_events_case1->Fill((*phoSeedTime)[mIDXSEL]);
                    spike_events->Fill((*phoSeedTime)[mIDXSEL]);
                    
                    spike_events_etaphi_case1->Fill((*phoEta)[mIDXSEL], (*phoPhi)[mIDXSEL]);
                    spike_events_etaphi->Fill((*phoEta)[mIDXSEL], (*phoPhi)[mIDXSEL]);

                    if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                    {
                        spike_events_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                    }

                    //cand1 searching adjencent cells
                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {                    
                        if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mIDXSEL] + 1)
                        {
                            cellsEB_right.push_back(icell); 

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mIDXSEL]) == 86 || (*phoSeedIEta)[mIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike1_err_right---------------------" << "\n"
                                << "phoSeediphi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }

                        //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && (*phoSeedIEta)[mIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                        {
                            cellsEB_right.push_back(icell); 

                            file9 << "-------------spike1_ietamid_right-------------" << "\n"
                            << "phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }

                        if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mIDXSEL] - 1)
                        {                    
                            cellsEB_left.push_back(icell);

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mIDXSEL]) == 86 || (*phoSeedIEta)[mIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike1_err_left---------------------" << "\n"
                                << "phoSeediphi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }            

                        //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && (*phoSeedIEta)[mIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                        {
                            cellsEB_left.push_back(icell); 

                            file9 << "-------------spike1_ietamid_left-------------" << "\n"
                            << "phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }  
                    
                    }

                    //Now I have to take if there's only one and both
                    //If there's only right
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                    {
                        rightcellsIDX_only = cellsEB_right[0];
                        maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                        file4 << "-------------spike1_rightonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_rightonly1->Fill(Ewing);
                        spike_eta_wing1->Fill(Ewing);
                        spike_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_case1->Fill((*phoSeedTime)[mIDXSEL]);
                            spike_events_ewing3->Fill((*phoSeedTime)[mIDXSEL]);

                            if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                spike_events_ewing3_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                            }
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_case1->Fill((*phoSeedTime)[mIDXSEL]);
                            spike_events_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                            if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                spike_events_ewing1_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                            }
                        }
                    }

                    //if there is only left
                    if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                    {
                        leftcellsIDX_only = cellsEB_left[0];
                        maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                        file4 << "-------------spike1_leftonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_leftonly1->Fill(Ewing);
                        spike_eta_wing1->Fill(Ewing);
                        spike_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_case1->Fill((*phoSeedTime)[mIDXSEL]);
                            spike_events_ewing3->Fill((*phoSeedTime)[mIDXSEL]);

                            if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                spike_events_ewing3_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                            }
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_case1->Fill((*phoSeedTime)[mIDXSEL]);
                            spike_events_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                            if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                spike_events_ewing1_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                            }
                        }
                    }

                    //if there are both
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                    {
                        rightcellsIDX = cellsEB_right[0];
                        leftcellsIDX = cellsEB_left[0];

                        if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[rightcellsIDX];
                        }
                        else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[leftcellsIDX];
                        }

                        Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                        file4 << "-------------spike1_both-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_both1->Fill(Ewing);
                        spike_eta_wing1->Fill(Ewing);
                        spike_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_case1->Fill((*phoSeedTime)[mIDXSEL]);
                            spike_events_ewing3->Fill((*phoSeedTime)[mIDXSEL]);

                            if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                spike_events_ewing3_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                            }
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_case1->Fill((*phoSeedTime)[mIDXSEL]);
                            spike_events_ewing1->Fill((*phoSeedTime)[mIDXSEL]);

                            if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                spike_events_ewing1_3ns->Fill((*phoSeedTime)[mIDXSEL]);
                            }
                        }
                    }
                }

                //time cut no showershape cut
                if ((*phoSeedTime)[mIDXSEL] < -12.5)
                {
                    spike_events_time_case1->Fill((*phoSeedTime)[mIDXSEL]);
                    spike_events_time->Fill((*phoSeedTime)[mIDXSEL]);
                    
                    spike_events_time_etaphi_case1->Fill((*phoEta)[mIDXSEL], (*phoPhi)[mIDXSEL]);
                    spike_events_time_etaphi->Fill((*phoEta)[mIDXSEL], (*phoPhi)[mIDXSEL]);

                    file6 << run << " " << lumis << " " << event << endl;

                    //cand1 searching adjencent cells
                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {                    
                        if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mIDXSEL] + 1)
                        {
                            cellsEB_right.push_back(icell); 

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mIDXSEL]) == 86 || (*phoSeedIEta)[mIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike1_time_err_right---------------------" << "\n"
                                << "phoSeediphi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }

                        //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && (*phoSeedIEta)[mIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                        {
                            cellsEB_right.push_back(icell); 

                            file10 << "-------------spike1_time_ietamid_right-------------" << "\n"
                            << "phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }

                        if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mIDXSEL] - 1)
                        {                    
                            cellsEB_left.push_back(icell);

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mIDXSEL]) == 86 || (*phoSeedIEta)[mIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike1_time_err_left---------------------" << "\n"
                                << "phoSeediphi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }            

                        //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && (*phoSeedIEta)[mIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                        {
                            cellsEB_left.push_back(icell); 

                            file10 << "-------------spike1_time_ietamid_left-------------" << "\n"
                            << "phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }  
                    
                    }

                    //Now I have to take if there's only one and both
                    //If there's only right
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                    {
                        rightcellsIDX_only = cellsEB_right[0];
                        maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                        file5 << "-------------spike1_time_rightonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_time_rightonly1->Fill(Ewing);
                        spike_eta_wing_time1->Fill(Ewing);
                        spike_eta_wing_time->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_time_case1->Fill((*phoSeedTime)[mIDXSEL]);
                            spike_events_ewing3_time->Fill((*phoSeedTime)[mIDXSEL]);

                            file15 << "--------spike case 1--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_time_case1->Fill((*phoSeedTime)[mIDXSEL]);
                            spike_events_ewing1_time->Fill((*phoSeedTime)[mIDXSEL]);

                            file16 << "--------spike case 1--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }
                    }

                    //if there is only left
                    if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                    {
                        leftcellsIDX_only = cellsEB_left[0];
                        maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                        file5 << "-------------spike1_time_leftonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_time_leftonly1->Fill(Ewing);
                        spike_eta_wing_time1->Fill(Ewing);
                        spike_eta_wing_time->Fill(Ewing);
                        
                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_time_case1->Fill((*phoSeedTime)[mIDXSEL]);
                            spike_events_ewing3_time->Fill((*phoSeedTime)[mIDXSEL]);

                            file15 << "--------spike case 1--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_time_case1->Fill((*phoSeedTime)[mIDXSEL]);
                            spike_events_ewing1_time->Fill((*phoSeedTime)[mIDXSEL]);

                            file16 << "--------spike case 1--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }
                    }

                    //if there are both
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                    {
                        rightcellsIDX = cellsEB_right[0];
                        leftcellsIDX = cellsEB_left[0];

                        if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[rightcellsIDX];
                        }
                        else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[leftcellsIDX];
                        }

                        Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mIDXSEL]);

                        file5 << "-------------spike1_time_both-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_time_both1->Fill(Ewing);
                        spike_eta_wing_time1->Fill(Ewing);
                        spike_eta_wing_time->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_time_case1->Fill((*phoSeedTime)[mIDXSEL]);
                            spike_events_ewing3_time->Fill((*phoSeedTime)[mIDXSEL]);

                            file15 << "--------spike case 1--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_time_case1->Fill((*phoSeedTime)[mIDXSEL]);
                            spike_events_ewing1_time->Fill((*phoSeedTime)[mIDXSEL]);

                            file16 << "--------spike case 1--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }
                    }
                }
            }
        }

        
        //---------------------------case 2, only oot, no it---------------------------
		if (itpho.size() == 0 && ootpho.size() != 0 && ootitmatched.size() != 0)
		{
			mIDXSEL = ootpho[0];
			misOOT = kTRUE;
			Float_t NewMET = newuMET(pfMET, pfMETPhi, (*ophoPhi)[mIDXSEL], (*ophoEt)[mIDXSEL]);

			//Candidate Events
			if (NewMET > 210 &&
				(*ophohasPixelSeed)[mIDXSEL] == 0 &&
				(*ophoMIPTotEnergy)[mIDXSEL] < 4.9 &&
				(*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 &&
				(*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001 &&
				(*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.01015 &&
				metFilters == 0)
			{
                candidate_events_case2->Fill((*ophoSeedTime)[mIDXSEL]);
				candidate_events->Fill((*ophoSeedTime)[mIDXSEL]);

                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                {
                    candidate_events_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                }

                if ((*ophoSeedTime)[mIDXSEL] < -12.5)
                {
                    spikeincandidate_events_case2->Fill((*ophoSeedTime)[mIDXSEL]);
				    spikeincandidate_events->Fill((*ophoSeedTime)[mIDXSEL]);
                }

                

                //cand1 searching adjencent cells
                for (int icell = 0; icell < nAllCellsEB; icell++)
                {
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mIDXSEL] + 1)
                    {
                        cellsEB_right.push_back(icell); 

                        //if there's any err
                        if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*ophoSeedIEta)[mIDXSEL]) == 86 || (*ophoSeedIEta)[mIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                        {
                            file7 << "------------------cand2_err_right---------------------" << "\n"
                            << "phoSeediphi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedieta = " << (*ophoSeedIEta)[mIDXSEL] << "\n"
                            << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                        }
                    }

                    //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && (*ophoSeedIEta)[mIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                    {
                        cellsEB_right.push_back(icell); 

                        file8 << "-------------Cand2_ietamid_right-------------" << "\n"
                        << "phoseedTime = " << (*ophoSeedTime)[mIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mIDXSEL] << "\n" 
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;
                    }

                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mIDXSEL] - 1)
                    {                    
                        cellsEB_left.push_back(icell);

                        //if there's any err
                        if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*ophoSeedIEta)[mIDXSEL]) == 86 || (*ophoSeedIEta)[mIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                        {
                            file7 << "------------------cand2_err_left---------------------" << "\n"
                            << "phoSeediphi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedieta = " << (*ophoSeedIEta)[mIDXSEL] << "\n"
                            << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                        }
                    }            

                    //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && (*ophoSeedIEta)[mIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                    {
                        cellsEB_left.push_back(icell); 

                        file8 << "-------------Cand2_ietamid_left-------------" << "\n"
                        << "phoseedTime = " << (*ophoSeedTime)[mIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;
                    }  
                
                }

                //Now I have to take if there's only one and both
                //If there's only right
                if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                {
                    rightcellsIDX_only = cellsEB_right[0];
                    maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                    file2 << "-------------Cand2_rightonly-------------" << "\n"
                    << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mIDXSEL] << "\n"
                    << "phoSeedIPhi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mIDXSEL] << "\n" 
                    << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                    << run << " " << lumis << " " << event << "\n"
                    << endl;

                    candidate_eta_wing_rightonly2->Fill(Ewing);
                    candidate_eta_wing2->Fill(Ewing);
                    candidate_eta_wing->Fill(Ewing);

                    if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                    {
                        candidate_eta_wing_prompt->Fill(Ewing);
                    }

                    if (Ewing > 0.03)
                    {
                        candidate_events_ewing3_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                        candidate_events_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing3_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.01)
                    {
                        candidate_events_ewing1_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                        candidate_events_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing1_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    //---------------------------more ewing cuts-------------------------------

                    if (Ewing > 0.02)
                    {
                        candidate_events_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing2_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.025)
                    {
                        candidate_events_ewing25->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing25_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.031)
                    {
                        candidate_events_ewing31->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing31_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.032)
                    {
                        candidate_events_ewing32->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing32_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.033)
                    {
                        candidate_events_ewing33->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing33_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.034)
                    {
                        candidate_events_ewing34->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing34_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.035)
                    {
                        candidate_events_ewing35->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing35_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.036)
                    {
                        candidate_events_ewing36->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing36_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.037)
                    {
                        candidate_events_ewing37->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing37_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.038)
                    {
                        candidate_events_ewing38->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing38_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.039)
                    {
                        candidate_events_ewing39->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing39_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }
                    

                    if (Ewing > 0.040)
                    {
                        candidate_events_ewing4->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing4_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.041)
                    {
                        candidate_events_ewing41->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing41_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.042)
                    {
                        candidate_events_ewing42->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing42_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.043)
                    {
                        candidate_events_ewing43->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing43_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.044)
                    {
                        candidate_events_ewing44->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing44_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.045)
                    {
                        candidate_events_ewing45->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing45_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.046)
                    {
                        candidate_events_ewing46->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing46_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.047)
                    {
                        candidate_events_ewing47->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing47_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.048)
                    {
                        candidate_events_ewing48->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing48_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.049)
                    {
                        candidate_events_ewing49->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing49_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }



                    if (Ewing > 0.05)
                    {
                        candidate_events_ewing5->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing5_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.07)
                    {
                        candidate_events_ewing7->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing7_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.09)
                    {
                        candidate_events_ewing9->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing9_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.10)
                    {
                        candidate_events_ewing10->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing10_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.12)
                    {
                        candidate_events_ewing12->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing12_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.14)
                    {
                        candidate_events_ewing14->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing14_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.15)
                    {
                        candidate_events_ewing15->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing15_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if ((*ophoSeedTime)[mIDXSEL] < -12.5)
                    {
                        spikeincandidate_eta_wing_rightonly2->Fill(Ewing);
                        spikeincandidate_eta_wing2->Fill(Ewing);
                        spikeincandidate_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spikeincandidate_events_ewing3_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                            spikeincandidate_events_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.01)
                        {
                            spikeincandidate_events_ewing1_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                            spikeincandidate_events_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        //---------------------------more ewing cuts-------------------------------

                        if (Ewing > 0.02)
                        {
                            spikeincandidate_events_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.025)
                        {
                            spikeincandidate_events_ewing25->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.031)
                        {
                            spikeincandidate_events_ewing31->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.032)
                        {
                            spikeincandidate_events_ewing32->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.033)
                        {
                            spikeincandidate_events_ewing33->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.034)
                        {
                            spikeincandidate_events_ewing34->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.035)
                        {
                            spikeincandidate_events_ewing35->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.036)
                        {
                            spikeincandidate_events_ewing36->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.037)
                        {
                            spikeincandidate_events_ewing37->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.038)
                        {
                            spikeincandidate_events_ewing38->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.039)
                        {
                            spikeincandidate_events_ewing39->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                        

                        if (Ewing > 0.040)
                        {
                            spikeincandidate_events_ewing4->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.041)
                        {
                            spikeincandidate_events_ewing41->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.042)
                        {
                            spikeincandidate_events_ewing42->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.043)
                        {
                            spikeincandidate_events_ewing43->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.044)
                        {
                            spikeincandidate_events_ewing44->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.045)
                        {
                            spikeincandidate_events_ewing45->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.046)
                        {
                            spikeincandidate_events_ewing46->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.047)
                        {
                            spikeincandidate_events_ewing47->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.048)
                        {
                            spikeincandidate_events_ewing48->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.049)
                        {
                            spikeincandidate_events_ewing49->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.05)
                        {
                            spikeincandidate_events_ewing5->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.07)
                        {
                            spikeincandidate_events_ewing7->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.09)
                        {
                            spikeincandidate_events_ewing9->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.10)
                        {
                            spikeincandidate_events_ewing10->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.12)
                        {
                            spikeincandidate_events_ewing12->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.14)
                        {
                            spikeincandidate_events_ewing14->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.15)
                        {
                            spikeincandidate_events_ewing15->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing < 0.03)
                    {
                        file13 << "--------cand case 2--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }

                    if (Ewing < 0.01)
                    {
                        file14 << "--------cand case 2--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }
                }

                //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                    file2 << "-------------Cand2_leftonly-------------" << "\n"
                    << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mIDXSEL] << "\n"
                    << "phoSeedIPhi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mIDXSEL] << "\n" 
                    << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                    << run << " " << lumis << " " << event << "\n"
                    << endl;

                    candidate_eta_wing_leftonly2->Fill(Ewing);
                    candidate_eta_wing2->Fill(Ewing);
                    candidate_eta_wing->Fill(Ewing);

                    if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                    {
                        candidate_eta_wing_prompt->Fill(Ewing);
                    }

                    if (Ewing > 0.03)
                    {
                        candidate_events_ewing3_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                        candidate_events_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing3_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.01)
                    {
                        candidate_events_ewing1_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                        candidate_events_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing1_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    //---------------------------more ewing cuts-------------------------------


                    if (Ewing > 0.02)
                    {
                        candidate_events_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing2_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.025)
                    {
                        candidate_events_ewing25->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing25_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.031)
                    {
                        candidate_events_ewing31->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing31_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.032)
                    {
                        candidate_events_ewing32->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing32_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.033)
                    {
                        candidate_events_ewing33->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing33_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.034)
                    {
                        candidate_events_ewing34->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing34_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.035)
                    {
                        candidate_events_ewing35->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing35_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.036)
                    {
                        candidate_events_ewing36->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing36_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.037)
                    {
                        candidate_events_ewing37->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing37_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.038)
                    {
                        candidate_events_ewing38->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing38_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.039)
                    {
                        candidate_events_ewing39->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing39_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }
                    

                    if (Ewing > 0.040)
                    {
                        candidate_events_ewing4->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing4_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.041)
                    {
                        candidate_events_ewing41->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing41_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.042)
                    {
                        candidate_events_ewing42->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing42_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.043)
                    {
                        candidate_events_ewing43->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing43_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.044)
                    {
                        candidate_events_ewing44->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing44_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.045)
                    {
                        candidate_events_ewing45->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing45_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.046)
                    {
                        candidate_events_ewing46->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing46_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.047)
                    {
                        candidate_events_ewing47->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing47_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.048)
                    {
                        candidate_events_ewing48->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing48_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.049)
                    {
                        candidate_events_ewing49->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing49_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.05)
                    {
                        candidate_events_ewing5->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing5_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.07)
                    {
                        candidate_events_ewing7->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing7_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.09)
                    {
                        candidate_events_ewing9->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing9_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.10)
                    {
                        candidate_events_ewing10->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing10_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.12)
                    {
                        candidate_events_ewing12->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing12_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.14)
                    {
                        candidate_events_ewing14->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing14_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.15)
                    {
                        candidate_events_ewing15->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing15_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if ((*ophoSeedTime)[mIDXSEL] < -12.5)
                    {
                        spikeincandidate_eta_wing_leftonly2->Fill(Ewing);
                        spikeincandidate_eta_wing2->Fill(Ewing);
                        spikeincandidate_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spikeincandidate_events_ewing3_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                            spikeincandidate_events_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.01)
                        {
                            spikeincandidate_events_ewing1_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                            spikeincandidate_events_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        //---------------------------more ewing cuts-------------------------------

                        if (Ewing > 0.02)
                        {
                            spikeincandidate_events_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.025)
                        {
                            spikeincandidate_events_ewing25->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.031)
                        {
                            spikeincandidate_events_ewing31->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.032)
                        {
                            spikeincandidate_events_ewing32->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.033)
                        {
                            spikeincandidate_events_ewing33->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.034)
                        {
                            spikeincandidate_events_ewing34->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.035)
                        {
                            spikeincandidate_events_ewing35->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.036)
                        {
                            spikeincandidate_events_ewing36->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.037)
                        {
                            spikeincandidate_events_ewing37->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.038)
                        {
                            spikeincandidate_events_ewing38->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.039)
                        {
                            spikeincandidate_events_ewing39->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                        

                        if (Ewing > 0.040)
                        {
                            spikeincandidate_events_ewing4->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.041)
                        {
                            spikeincandidate_events_ewing41->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.042)
                        {
                            spikeincandidate_events_ewing42->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.043)
                        {
                            spikeincandidate_events_ewing43->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.044)
                        {
                            spikeincandidate_events_ewing44->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.045)
                        {
                            spikeincandidate_events_ewing45->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.046)
                        {
                            spikeincandidate_events_ewing46->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.047)
                        {
                            spikeincandidate_events_ewing47->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.048)
                        {
                            spikeincandidate_events_ewing48->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.049)
                        {
                            spikeincandidate_events_ewing49->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.05)
                        {
                            spikeincandidate_events_ewing5->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.07)
                        {
                            spikeincandidate_events_ewing7->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.09)
                        {
                            spikeincandidate_events_ewing9->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.10)
                        {
                            spikeincandidate_events_ewing10->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.12)
                        {
                            spikeincandidate_events_ewing12->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.14)
                        {
                            spikeincandidate_events_ewing14->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.15)
                        {
                            spikeincandidate_events_ewing15->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing < 0.03)
                    {
                        file13 << "--------cand case 2--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }

                    if (Ewing < 0.01)
                    {
                        file14 << "--------cand case 2--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }
                }

                //if there are both
                if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                {
                    rightcellsIDX = cellsEB_right[0];
                    leftcellsIDX = cellsEB_left[0];

                    if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[rightcellsIDX];
                    }
                    else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[leftcellsIDX];
                    }

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                    file2 << "-------------Cand2_both-------------" << "\n"
                    << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mIDXSEL] << "\n"
                    << "phoSeedIPhi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mIDXSEL] << "\n" 
                    << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                    << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                    << run << " " << lumis << " " << event << "\n"
                    << endl;

                    candidate_eta_wing_both2->Fill(Ewing);
                    candidate_eta_wing2->Fill(Ewing);
                    candidate_eta_wing->Fill(Ewing);

                    if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                    {
                        candidate_eta_wing_prompt->Fill(Ewing);
                    }

                    if (Ewing > 0.03)
                    {
                        candidate_events_ewing3_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                        candidate_events_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing3_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.01)
                    {
                        candidate_events_ewing1_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                        candidate_events_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing1_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    //---------------------------more ewing cuts-------------------------------

                    if (Ewing > 0.02)
                    {
                        candidate_events_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing2_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.025)
                    {
                        candidate_events_ewing25->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing25_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.031)
                    {
                        candidate_events_ewing31->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing31_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.032)
                    {
                        candidate_events_ewing32->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing32_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.033)
                    {
                        candidate_events_ewing33->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing33_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.034)
                    {
                        candidate_events_ewing34->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing34_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.035)
                    {
                        candidate_events_ewing35->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing35_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.036)
                    {
                        candidate_events_ewing36->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing36_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.037)
                    {
                        candidate_events_ewing37->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing37_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.038)
                    {
                        candidate_events_ewing38->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing38_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.039)
                    {
                        candidate_events_ewing39->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing39_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }
                    

                    if (Ewing > 0.040)
                    {
                        candidate_events_ewing4->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing4_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.041)
                    {
                        candidate_events_ewing41->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing41_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.042)
                    {
                        candidate_events_ewing42->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing42_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.043)
                    {
                        candidate_events_ewing43->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing43_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.044)
                    {
                        candidate_events_ewing44->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing44_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.045)
                    {
                        candidate_events_ewing45->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing45_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.046)
                    {
                        candidate_events_ewing46->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing46_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.047)
                    {
                        candidate_events_ewing47->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing47_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.048)
                    {
                        candidate_events_ewing48->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing48_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.049)
                    {
                        candidate_events_ewing49->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing49_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.05)
                    {
                        candidate_events_ewing5->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing5_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.07)
                    {
                        candidate_events_ewing7->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing7_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.09)
                    {
                        candidate_events_ewing9->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing9_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.10)
                    {
                        candidate_events_ewing10->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing10_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.12)
                    {
                        candidate_events_ewing12->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing12_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.14)
                    {
                        candidate_events_ewing14->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing14_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing > 0.15)
                    {
                        candidate_events_ewing15->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing15_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if ((*ophoSeedTime)[mIDXSEL] < -12.5)
                    {
                        spikeincandidate_eta_wing_both2->Fill(Ewing);
                        spikeincandidate_eta_wing2->Fill(Ewing);
                        spikeincandidate_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spikeincandidate_events_ewing3_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                            spikeincandidate_events_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.01)
                        {
                            spikeincandidate_events_ewing1_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                            spikeincandidate_events_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        //---------------------------more ewing cuts-------------------------------

                        if (Ewing > 0.02)
                        {
                            spikeincandidate_events_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.025)
                        {
                            spikeincandidate_events_ewing25->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.031)
                        {
                            spikeincandidate_events_ewing31->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.032)
                        {
                            spikeincandidate_events_ewing32->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.033)
                        {
                            spikeincandidate_events_ewing33->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.034)
                        {
                            spikeincandidate_events_ewing34->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.035)
                        {
                            spikeincandidate_events_ewing35->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.036)
                        {
                            spikeincandidate_events_ewing36->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.037)
                        {
                            spikeincandidate_events_ewing37->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.038)
                        {
                            spikeincandidate_events_ewing38->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.039)
                        {
                            spikeincandidate_events_ewing39->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                        

                        if (Ewing > 0.040)
                        {
                            spikeincandidate_events_ewing4->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.041)
                        {
                            spikeincandidate_events_ewing41->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.042)
                        {
                            spikeincandidate_events_ewing42->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.043)
                        {
                            spikeincandidate_events_ewing43->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.044)
                        {
                            spikeincandidate_events_ewing44->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.045)
                        {
                            spikeincandidate_events_ewing45->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.046)
                        {
                            spikeincandidate_events_ewing46->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.047)
                        {
                            spikeincandidate_events_ewing47->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.048)
                        {
                            spikeincandidate_events_ewing48->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.049)
                        {
                            spikeincandidate_events_ewing49->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.05)
                        {
                            spikeincandidate_events_ewing5->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.07)
                        {
                            spikeincandidate_events_ewing7->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.09)
                        {
                            spikeincandidate_events_ewing9->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.10)
                        {
                            spikeincandidate_events_ewing10->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.12)
                        {
                            spikeincandidate_events_ewing12->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.14)
                        {
                            spikeincandidate_events_ewing14->Fill((*ophoSeedTime)[mIDXSEL]);
                        }

                        if (Ewing > 0.15)
                        {
                            spikeincandidate_events_ewing15->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
                    }

                    if (Ewing < 0.03)
                    {
                        file13 << "--------cand case 2--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }

                    if (Ewing < 0.01)
                    {
                        file14 << "--------cand case 2--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }
                }
            }

            if (iiZpho.size() == 0 && ooZpho.size() != 0 && ootitZmatched.size() != 0)
            {
                mZIDXSEL = ooZpho[0];
                //Prompt Z template:
                if ((*ophohasPixelSeed)[mIDXSEL] == 1 &&
                    (*ophoMIPTotEnergy)[mIDXSEL] < 4.9 &&
                    (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 &&
                    (*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001 &&
                    (*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.01015 &&
                    metFilters == 0 &&
                    (*ophohasPixelSeed)[mZIDXSEL] == 1 &&
                    (*ophoMIPTotEnergy)[mZIDXSEL] < 4.9 &&
                    (*ophoSigmaIEtaIEtaFull5x5)[mZIDXSEL] > 0.001 &&
                    (*ophoSigmaIPhiIPhiFull5x5)[mZIDXSEL] > 0.001 &&
                    (*ophoSigmaIEtaIEtaFull5x5)[mZIDXSEL] < 0.01015)
                {
                    Float_t InvM = InvariMass((*ophoEt)[mIDXSEL], (*ophoEt)[mZIDXSEL], (*ophoPhi)[mIDXSEL], (*ophoPhi)[mZIDXSEL], (*ophoEta)[mIDXSEL], (*ophoEta)[mZIDXSEL]);
                    InvMass_Z->Fill(InvM);

                    if (InvM > 85 && InvM < 100)
                    {
                        promptZ_events_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                        promptZ_events->Fill((*ophoSeedTime)[mIDXSEL]);

                        if (fabs((*phoSeedTime)[mIDXSEL]) < 3.0)
                        {
                            promptZ_events_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                        }
    
                        //promptZ1 searching adjencent cells
                        for (int icell = 0; icell < nAllCellsEB; icell++)
                        {                        
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mIDXSEL] + 1)
                            {
                                cellsEB_right.push_back(icell); 

                                //if there's any err
                                if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mIDXSEL]) == 86 || (*phoSeedIEta)[mIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                                {
                                    file7 << "------------------promptZ2_err_right---------------------" << "\n"
                                    << "phoSeediphi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mIDXSEL] << "\n"
                                    << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                                }
                            }

                            //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && (*phoSeedIEta)[mIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                            {
                                cellsEB_right.push_back(icell); 

                                file9 << "-------------promptZ2_ietamid_right-------------" << "\n"
                                << "phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                                << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                                << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                                << run << " " << lumis << " " << event << "\n"
                                << endl;
                            }

                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mIDXSEL] - 1)
                            {                    
                                cellsEB_left.push_back(icell);

                                //if there's any err
                                if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mIDXSEL]) == 86 || (*phoSeedIEta)[mIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                                {
                                    file7 << "------------------promptZ2_err_left---------------------" << "\n"
                                    << "phoSeediphi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mIDXSEL] << "\n"
                                    << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                                }
                            }            

                            //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mIDXSEL] && (*phoSeedIEta)[mIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                            {
                                cellsEB_left.push_back(icell); 

                                file9 << "-------------promptZ2_ietamid_left-------------" << "\n"
                                << "phoseedTime = " << (*phoSeedTime)[mIDXSEL] << "\n"
                                << "phoSeedIPhi = " << (*phoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mIDXSEL] << "\n" 
                                << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                                << run << " " << lumis << " " << event << "\n"
                                << endl;
                            }  
                        }

                        //Now I have to take if there's only one and both
                        //If there's only right
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                        {
                            rightcellsIDX_only = cellsEB_right[0];
                            maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                            file3 << "-------------promptZ2_rightonly-------------" << "\n"
                            << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mIDXSEL] << "\n" 
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;

                            promptZ_eta_wing_rightonly2->Fill(Ewing);
                            promptZ_eta_wing2->Fill(Ewing);
                            promptZ_eta_wing->Fill(Ewing);

                            if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                promptZ_eta_wing_prompt->Fill(Ewing);
                            }

                            if (Ewing > 0.03)
                            {
                                promptZ_events_ewing3_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                                promptZ_events_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing3_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.01)
                            {
                                promptZ_events_ewing1_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                                promptZ_events_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing1_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            //---------------------------more ewing cuts-------------------------------

                            if (Ewing > 0.02)
                            {
                                promptZ_events_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing2_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.025)
                            {
                                promptZ_events_ewing25->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing25_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.031)
                            {
                                promptZ_events_ewing31->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing31_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.032)
                            {
                                promptZ_events_ewing32->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing32_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.033)
                            {
                                promptZ_events_ewing33->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing33_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.034)
                            {
                                promptZ_events_ewing34->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing34_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.035)
                            {
                                promptZ_events_ewing35->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing35_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.036)
                            {
                                promptZ_events_ewing36->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing36_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.037)
                            {
                                promptZ_events_ewing37->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing37_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.038)
                            {
                                promptZ_events_ewing38->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing38_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.039)
                            {
                                promptZ_events_ewing39->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing39_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }
                            

                            if (Ewing > 0.040)
                            {
                                promptZ_events_ewing4->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing4_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.041)
                            {
                                promptZ_events_ewing41->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing41_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.042)
                            {
                                promptZ_events_ewing42->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing42_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.043)
                            {
                                promptZ_events_ewing43->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing43_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.044)
                            {
                                promptZ_events_ewing44->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing44_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.045)
                            {
                                promptZ_events_ewing45->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing45_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.046)
                            {
                                promptZ_events_ewing46->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing46_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.047)
                            {
                                promptZ_events_ewing47->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing47_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.048)
                            {
                                promptZ_events_ewing48->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing48_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.049)
                            {
                                promptZ_events_ewing49->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing49_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.05)
                            {
                                promptZ_events_ewing5->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing5_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.07)
                            {
                                promptZ_events_ewing7->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing7_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.09)
                            {
                                promptZ_events_ewing9->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing9_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.10)
                            {
                                promptZ_events_ewing10->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing10_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.12)
                            {
                                promptZ_events_ewing12->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing12_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.14)
                            {
                                promptZ_events_ewing14->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing14_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.15)
                            {
                                promptZ_events_ewing15->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing15_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }
                        }

                        //if there is only left
                        if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                        {
                            leftcellsIDX_only = cellsEB_left[0];
                            maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                            file3 << "-------------promptZ2_leftonly-------------" << "\n"
                            << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;

                            promptZ_eta_wing_leftonly2->Fill(Ewing);
                            promptZ_eta_wing2->Fill(Ewing);
                            promptZ_eta_wing->Fill(Ewing);

                            if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                promptZ_eta_wing_prompt->Fill(Ewing);
                            }

                            if (Ewing > 0.03)
                            {
                                promptZ_events_ewing3_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                                promptZ_events_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing3_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.01)
                            {
                                promptZ_events_ewing1_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                                promptZ_events_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing1_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            //---------------------------more ewing cuts-------------------------------

                            if (Ewing > 0.02)
                            {
                                promptZ_events_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing2_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.025)
                            {
                                promptZ_events_ewing25->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing25_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.031)
                            {
                                promptZ_events_ewing31->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing31_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.032)
                            {
                                promptZ_events_ewing32->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing32_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.033)
                            {
                                promptZ_events_ewing33->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing33_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.034)
                            {
                                promptZ_events_ewing34->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing34_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.035)
                            {
                                promptZ_events_ewing35->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing35_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.036)
                            {
                                promptZ_events_ewing36->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing36_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.037)
                            {
                                promptZ_events_ewing37->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing37_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.038)
                            {
                                promptZ_events_ewing38->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing38_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.039)
                            {
                                promptZ_events_ewing39->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing39_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }
                            

                            if (Ewing > 0.040)
                            {
                                promptZ_events_ewing4->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing4_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.041)
                            {
                                promptZ_events_ewing41->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing41_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.042)
                            {
                                promptZ_events_ewing42->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing42_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.043)
                            {
                                promptZ_events_ewing43->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing43_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.044)
                            {
                                promptZ_events_ewing44->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing44_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.045)
                            {
                                promptZ_events_ewing45->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing45_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.046)
                            {
                                promptZ_events_ewing46->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing46_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.047)
                            {
                                promptZ_events_ewing47->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing47_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.048)
                            {
                                promptZ_events_ewing48->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing48_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.049)
                            {
                                promptZ_events_ewing49->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing49_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.05)
                            {
                                promptZ_events_ewing5->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing5_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.07)
                            {
                                promptZ_events_ewing7->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing7_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.09)
                            {
                                promptZ_events_ewing9->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing9_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.10)
                            {
                                promptZ_events_ewing10->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing10_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.12)
                            {
                                promptZ_events_ewing12->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing12_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.14)
                            {
                                promptZ_events_ewing14->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing14_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.15)
                            {
                                promptZ_events_ewing15->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing15_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }
                        }

                        //if there are both
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                        {
                            rightcellsIDX = cellsEB_right[0];
                            leftcellsIDX = cellsEB_left[0];

                            if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[rightcellsIDX];
                            }
                            else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[leftcellsIDX];
                            }

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                            file3 << "-------------promptZ2_both-------------" << "\n"
                            << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;

                            promptZ_eta_wing_both2->Fill(Ewing);
                            promptZ_eta_wing2->Fill(Ewing);
                            promptZ_eta_wing->Fill(Ewing);

                            if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                promptZ_eta_wing_prompt->Fill(Ewing);
                            }

                            if (Ewing > 0.03)
                            {
                                promptZ_events_ewing3_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                                promptZ_events_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing3_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.01)
                            {
                                promptZ_events_ewing1_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                                promptZ_events_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing1_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            //---------------------------more ewing cuts-------------------------------

                            if (Ewing > 0.02)
                            {
                                promptZ_events_ewing2->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing2_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.025)
                            {
                                promptZ_events_ewing25->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing25_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.031)
                            {
                                promptZ_events_ewing31->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing31_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.032)
                            {
                                promptZ_events_ewing32->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing32_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.033)
                            {
                                promptZ_events_ewing33->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing33_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.034)
                            {
                                promptZ_events_ewing34->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing34_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.035)
                            {
                                promptZ_events_ewing35->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing35_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.036)
                            {
                                promptZ_events_ewing36->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing36_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.037)
                            {
                                promptZ_events_ewing37->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing37_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.038)
                            {
                                promptZ_events_ewing38->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing38_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.039)
                            {
                                promptZ_events_ewing39->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing39_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }
                            

                            if (Ewing > 0.040)
                            {
                                promptZ_events_ewing4->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing4_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.041)
                            {
                                promptZ_events_ewing41->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing41_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.042)
                            {
                                promptZ_events_ewing42->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing42_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.043)
                            {
                                promptZ_events_ewing43->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing43_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.044)
                            {
                                promptZ_events_ewing44->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing44_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.045)
                            {
                                promptZ_events_ewing45->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing45_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.046)
                            {
                                promptZ_events_ewing46->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing46_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.047)
                            {
                                promptZ_events_ewing47->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing47_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.048)
                            {
                                promptZ_events_ewing48->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing48_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.049)
                            {
                                promptZ_events_ewing49->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing49_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.05)
                            {
                                promptZ_events_ewing5->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing5_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.07)
                            {
                                promptZ_events_ewing7->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing7_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.09)
                            {
                                promptZ_events_ewing9->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing9_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.10)
                            {
                                promptZ_events_ewing10->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing10_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.12)
                            {
                                promptZ_events_ewing12->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing12_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.14)
                            {
                                promptZ_events_ewing14->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing14_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }

                            if (Ewing > 0.15)
                            {
                                promptZ_events_ewing15->Fill((*ophoSeedTime)[mIDXSEL]);

                                if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing15_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                                }
                            }
                        } 
                    }
                }
            }

            //BeamHalo template
			if (NewMET > 210 &&
				(*ophohasPixelSeed)[mIDXSEL] == 0 &&
				(*ophoMIPTotEnergy)[mIDXSEL] > 4.9 &&
				(*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.001 &&
				(*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] > 0.001 &&
				metFilters == 8)
			{
				BeamHalo_template_case2->Fill((*ophoSeedTime)[mIDXSEL]);
				BeamHalo_template->Fill((*ophoSeedTime)[mIDXSEL]);

                sieie_BeamHalo->Fill((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL]);

                if ((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.0176)
                {
                    BeamHalo_template_sieieb176->Fill((*ophoSeedTime)[mIDXSEL]);
                }
                
                if ((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] > 0.0206)
                {
                    BeamHalo_template_sieieb206->Fill((*ophoSeedTime)[mIDXSEL]);
                }

                if ((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.0176)
                {
                    BeamHalo_template_sieies176->Fill((*ophoSeedTime)[mIDXSEL]);
                }

                if ((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.0206)
                {
                    BeamHalo_template_sieies206->Fill((*ophoSeedTime)[mIDXSEL]);
                }
            }

            //Spikes template
            if (NewMET > 210 &&
                (*ophohasPixelSeed)[mIDXSEL] == 0 &&
                (*ophoMIPTotEnergy)[mIDXSEL] < 4.9 &&
                ((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.001 ||
                (*ophoSigmaIPhiIPhiFull5x5)[mIDXSEL] < 0.001) &&
                metFilters == 0)
            {
                //original showershap cut
                if ((*ophoSigmaIEtaIEtaFull5x5)[mIDXSEL] < 0.01015)
                {  
                    spike_events_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                    spike_events->Fill((*ophoSeedTime)[mIDXSEL]);
                    
                    spike_events_etaphi_case2->Fill((*ophoEta)[mIDXSEL], (*ophoPhi)[mIDXSEL]);
                    spike_events_etaphi->Fill((*ophoEta)[mIDXSEL], (*ophoPhi)[mIDXSEL]);

                    if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                    {
                        spike_events_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                    }

                    //cand1 searching adjencent cells
                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {                    
                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mIDXSEL] + 1)
                        {
                            cellsEB_right.push_back(icell); 

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*ophoSeedIEta)[mIDXSEL]) == 86 || (*ophoSeedIEta)[mIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike2_err_right---------------------" << "\n"
                                << "phoSeediphi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedieta = " << (*ophoSeedIEta)[mIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }

                        //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && (*ophoSeedIEta)[mIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                        {
                            cellsEB_right.push_back(icell); 

                            file9 << "-------------spike2_ietamid_right-------------" << "\n"
                            << "phoseedTime = " << (*ophoSeedTime)[mIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mIDXSEL] << "\n" 
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }

                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mIDXSEL] - 1)
                        {                    
                            cellsEB_left.push_back(icell);

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*ophoSeedIEta)[mIDXSEL]) == 86 || (*ophoSeedIEta)[mIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike2_err_left---------------------" << "\n"
                                << "phoSeediphi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedieta = " << (*ophoSeedIEta)[mIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }            

                        //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && (*ophoSeedIEta)[mIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                        {
                            cellsEB_left.push_back(icell); 

                            file9 << "-------------spike2_ietamid_left-------------" << "\n"
                            << "phoseedTime = " << (*ophoSeedTime)[mIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }  
                    
                    }

                    //Now I have to take if there's only one and both
                    //If there's only right
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                    {
                        rightcellsIDX_only = cellsEB_right[0];
                        maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                        file4 << "-------------spike2_rightonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mIDXSEL] << "\n" 
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_rightonly2->Fill(Ewing);
                        spike_eta_wing2->Fill(Ewing);
                        spike_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                            spike_events_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);

                            if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                spike_events_ewing3_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                            }
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                            spike_events_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                            if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                spike_events_ewing1_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                            }
                        }
                    }

                    //if there is only left
                    if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                    {
                        leftcellsIDX_only = cellsEB_left[0];
                        maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                        file4 << "-------------spike2_leftonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_leftonly2->Fill(Ewing);
                        spike_eta_wing2->Fill(Ewing);
                        spike_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                            spike_events_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);

                            if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                spike_events_ewing3_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                            }
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                            spike_events_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                            if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                spike_events_ewing1_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                            }
                        }
                    }

                    //if there are both
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                    {
                        rightcellsIDX = cellsEB_right[0];
                        leftcellsIDX = cellsEB_left[0];

                        if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[rightcellsIDX];
                        }
                        else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[leftcellsIDX];
                        }

                        Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                        file4 << "-------------spike2_both-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_both2->Fill(Ewing);
                        spike_eta_wing2->Fill(Ewing);
                        spike_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                            spike_events_ewing3->Fill((*ophoSeedTime)[mIDXSEL]);

                            if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                spike_events_ewing3_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                            }
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                            spike_events_ewing1->Fill((*ophoSeedTime)[mIDXSEL]);

                            if (fabs((*ophoSeedTime)[mIDXSEL]) < 3.0)
                            {
                                spike_events_ewing1_3ns->Fill((*ophoSeedTime)[mIDXSEL]);
                            }
                        }
                    }
                }

                //time cut no showershape cut
                if ((*ophoSeedTime)[mIDXSEL] < -12.5)
                {
                    spike_events_time_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                    spike_events_time->Fill((*ophoSeedTime)[mIDXSEL]);
                    
                    spike_events_time_etaphi_case2->Fill((*ophoEta)[mIDXSEL], (*ophoPhi)[mIDXSEL]);
                    spike_events_time_etaphi->Fill((*ophoEta)[mIDXSEL], (*ophoPhi)[mIDXSEL]);

                    file6 << run << " " << lumis << " " << event << endl;

                    //cand1 searching adjencent cells
                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {                    
                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mIDXSEL] + 1)
                        {
                            cellsEB_right.push_back(icell); 

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*ophoSeedIEta)[mIDXSEL]) == 86 || (*ophoSeedIEta)[mIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike2_time_err_right---------------------" << "\n"
                                << "phoSeediphi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedieta = " << (*ophoSeedIEta)[mIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }

                        //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && (*ophoSeedIEta)[mIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                        {
                            cellsEB_right.push_back(icell); 

                            file10 << "-------------spike2_time_ietamid_right-------------" << "\n"
                            << "phoseedTime = " << (*ophoSeedTime)[mIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mIDXSEL] << "\n" 
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }

                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mIDXSEL] - 1)
                        {                    
                            cellsEB_left.push_back(icell);

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*ophoSeedIEta)[mIDXSEL]) == 86 || (*ophoSeedIEta)[mIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike2_time_err_left---------------------" << "\n"
                                << "phoSeediphi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedieta = " << (*ophoSeedIEta)[mIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }            

                        //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mIDXSEL] && (*ophoSeedIEta)[mIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                        {
                            cellsEB_left.push_back(icell); 

                            file10 << "-------------spike2_time_ietamid_left-------------" << "\n"
                            << "phoseedTime = " << (*ophoSeedTime)[mIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }
                    }

                    //Now I have to take if there's only one and both
                    //If there's only right
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                    {
                        rightcellsIDX_only = cellsEB_right[0];
                        maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                        file5 << "-------------spike2_time_rightonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mIDXSEL] << "\n" 
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_time_rightonly2->Fill(Ewing);
                        spike_eta_wing_time2->Fill(Ewing);
                        spike_eta_wing_time->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_time_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                            spike_events_ewing3_time->Fill((*ophoSeedTime)[mIDXSEL]);

                            file15 << "--------spike case 2--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_time_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                            spike_events_ewing1_time->Fill((*ophoSeedTime)[mIDXSEL]);

                            file16 << "--------spike case 2--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }
                    }

                    //if there is only left
                    if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                    {
                        leftcellsIDX_only = cellsEB_left[0];
                        maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                        file5 << "-------------spike2_time_leftonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_time_leftonly2->Fill(Ewing);
                        spike_eta_wing_time2->Fill(Ewing);
                        spike_eta_wing_time->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_time_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                            spike_events_ewing3_time->Fill((*ophoSeedTime)[mIDXSEL]);

                            file15 << "--------spike case 2--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_time_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                            spike_events_ewing1_time->Fill((*ophoSeedTime)[mIDXSEL]);

                            file16 << "--------spike case 2--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }
                    }

                    //if there are both
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                    {
                        rightcellsIDX = cellsEB_right[0];
                        leftcellsIDX = cellsEB_left[0];

                        if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[rightcellsIDX];
                        }
                        else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[leftcellsIDX];
                        }

                        Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mIDXSEL]);

                        file5 << "-------------spike2_time_both-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_time_both2->Fill(Ewing);
                        spike_eta_wing_time2->Fill(Ewing);
                        spike_eta_wing_time->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_time_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                            spike_events_ewing3_time->Fill((*ophoSeedTime)[mIDXSEL]);

                            file15 << "--------spike case 2--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_time_case2->Fill((*ophoSeedTime)[mIDXSEL]);
                            spike_events_ewing1_time->Fill((*ophoSeedTime)[mIDXSEL]);

                            file16 << "--------spike case 2--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }
                    }
                }
            }
        }

        //cases when there are both
		if (itpho.size() != 0 && ootpho.size() != 0)
		{
			//case 3, only it pass trigger
			mootIDXSEL = ootitmatched[0];
			mitIDXSEL = itootmatched[0];

			//Candidate Events
			if (pfMET > 210 &&
				(*phohasPixelSeed)[mitIDXSEL] == 0 &&
				(*phoMIPTotEnergy)[mitIDXSEL] < 4.9 &&
				(*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.001 &&
				(*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] > 0.001 &&
				(*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] < 0.01015 &&
				metFilters == 0)
			{
                candidate_events_case3->Fill((*phoSeedTime)[mitIDXSEL]);
				candidate_events->Fill((*phoSeedTime)[mitIDXSEL]);

                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                {
                    candidate_events_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                }

                if ((*phoSeedTime)[mitIDXSEL] < -12.5)
                {
                    spikeincandidate_events_case3->Fill((*phoSeedTime)[mitIDXSEL]);
				    spikeincandidate_events->Fill((*phoSeedTime)[mitIDXSEL]);
                }

                //cand1 searching adjencent cells
                for (int icell = 0; icell < nAllCellsEB; icell++)
                {                    
                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mitIDXSEL] + 1)
                    {
                        cellsEB_right.push_back(icell); 

                        //if there's any err
                        if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mitIDXSEL]) == 86 || (*phoSeedIEta)[mitIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                        {
                            file7 << "------------------cand3_err_right---------------------" << "\n"
                            << "phoSeediphi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mitIDXSEL] << "\n"
                            << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                        }
                    }

                    //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && (*phoSeedIEta)[mitIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                    {
                        cellsEB_right.push_back(icell); 

                        file8 << "-------------Cand3_ietamid_right-------------" << "\n"
                        << "phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;
                    }

                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mitIDXSEL] - 1)
                    {                    
                        cellsEB_left.push_back(icell);

                        //if there's any err
                        if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mitIDXSEL]) == 86 || (*phoSeedIEta)[mitIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                        {
                            file7 << "------------------cand3_err_left---------------------" << "\n"
                            << "phoSeediphi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mitIDXSEL] << "\n"
                            << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                        }
                    }            

                    //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && (*phoSeedIEta)[mitIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                    {
                        cellsEB_left.push_back(icell); 

                        file8 << "-------------Cand3_ietamid_left-------------" << "\n"
                        << "phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;
                    }        
                }

                //Now I have to take if there's only one and both
                //If there's only right
                if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                {
                    rightcellsIDX_only = cellsEB_right[0];
                    maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);

                    file2 << "-------------Cand3_rightonly-------------" << "\n"
                    << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                    << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                    << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                    << run << " " << lumis << " " << event << "\n"
                    << endl;

                    candidate_eta_wing_rightonly3->Fill(Ewing);
                    candidate_eta_wing3->Fill(Ewing);
                    candidate_eta_wing->Fill(Ewing);

                    if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                    {
                        candidate_eta_wing_prompt->Fill(Ewing);
                    }

                    if (Ewing > 0.03)
                    {
                        candidate_events_ewing3_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                        candidate_events_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing3_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.01)
                    {
                        candidate_events_ewing1_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                        candidate_events_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing1_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    //---------------------------more ewing cuts-------------------------------


                    if (Ewing > 0.02)
                    {
                        candidate_events_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing2_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.025)
                    {
                        candidate_events_ewing25->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing25_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.031)
                    {
                        candidate_events_ewing31->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing31_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.032)
                    {
                        candidate_events_ewing32->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing32_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.033)
                    {
                        candidate_events_ewing33->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing33_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.034)
                    {
                        candidate_events_ewing34->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing34_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.035)
                    {
                        candidate_events_ewing35->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing35_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.036)
                    {
                        candidate_events_ewing36->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing36_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.037)
                    {
                        candidate_events_ewing37->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing37_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.038)
                    {
                        candidate_events_ewing38->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing38_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.039)
                    {
                        candidate_events_ewing39->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing39_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }
                    

                    if (Ewing > 0.040)
                    {
                        candidate_events_ewing4->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing4_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.041)
                    {
                        candidate_events_ewing41->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing41_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.042)
                    {
                        candidate_events_ewing42->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing42_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.043)
                    {
                        candidate_events_ewing43->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing43_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.044)
                    {
                        candidate_events_ewing44->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing44_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.045)
                    {
                        candidate_events_ewing45->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing45_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.046)
                    {
                        candidate_events_ewing46->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing46_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.047)
                    {
                        candidate_events_ewing47->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing47_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.048)
                    {
                        candidate_events_ewing48->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing48_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.049)
                    {
                        candidate_events_ewing49->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing49_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }



                    if (Ewing > 0.05)
                    {
                        candidate_events_ewing5->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing5_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.07)
                    {
                        candidate_events_ewing7->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing7_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.09)
                    {
                        candidate_events_ewing9->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing9_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.10)
                    {
                        candidate_events_ewing10->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing10_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.12)
                    {
                        candidate_events_ewing12->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing12_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.14)
                    {
                        candidate_events_ewing14->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing14_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.15)
                    {
                        candidate_events_ewing15->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing15_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if ((*phoSeedTime)[mitIDXSEL] < -12.5)
                    {
                        spikeincandidate_eta_wing_rightonly3->Fill(Ewing);
                        spikeincandidate_eta_wing3->Fill(Ewing);
                        spikeincandidate_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spikeincandidate_events_ewing3_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                            spikeincandidate_events_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.01)
                        {
                            spikeincandidate_events_ewing1_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                            spikeincandidate_events_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        //---------------------------more ewing cuts-------------------------------

                        if (Ewing > 0.02)
                        {
                            spikeincandidate_events_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.025)
                        {
                            spikeincandidate_events_ewing25->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.031)
                        {
                            spikeincandidate_events_ewing31->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.032)
                        {
                            spikeincandidate_events_ewing32->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.033)
                        {
                            spikeincandidate_events_ewing33->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.034)
                        {
                            spikeincandidate_events_ewing34->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.035)
                        {
                            spikeincandidate_events_ewing35->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.036)
                        {
                            spikeincandidate_events_ewing36->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.037)
                        {
                            spikeincandidate_events_ewing37->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.038)
                        {
                            spikeincandidate_events_ewing38->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.039)
                        {
                            spikeincandidate_events_ewing39->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                        

                        if (Ewing > 0.040)
                        {
                            spikeincandidate_events_ewing4->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.041)
                        {
                            spikeincandidate_events_ewing41->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.042)
                        {
                            spikeincandidate_events_ewing42->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.043)
                        {
                            spikeincandidate_events_ewing43->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.044)
                        {
                            spikeincandidate_events_ewing44->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.045)
                        {
                            spikeincandidate_events_ewing45->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.046)
                        {
                            spikeincandidate_events_ewing46->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.047)
                        {
                            spikeincandidate_events_ewing47->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.048)
                        {
                            spikeincandidate_events_ewing48->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.049)
                        {
                            spikeincandidate_events_ewing49->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.05)
                        {
                            spikeincandidate_events_ewing5->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.07)
                        {
                            spikeincandidate_events_ewing7->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.09)
                        {
                            spikeincandidate_events_ewing9->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.10)
                        {
                            spikeincandidate_events_ewing10->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.12)
                        {
                            spikeincandidate_events_ewing12->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.14)
                        {
                            spikeincandidate_events_ewing14->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.15)
                        {
                            spikeincandidate_events_ewing15->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing < 0.03)
                    {
                        file13 << "--------cand case 3--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }

                    if (Ewing < 0.01)
                    {
                        file14 << "--------cand case 3--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }
                }

                //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);

                    file2 << "-------------Cand3_leftonly-------------" << "\n"
                    << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                    << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                    << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                    << run << " " << lumis << " " << event << "\n"
                    << endl;

                    candidate_eta_wing_leftonly3->Fill(Ewing);
                    candidate_eta_wing3->Fill(Ewing);
                    candidate_eta_wing->Fill(Ewing);

                    if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                    {
                        candidate_eta_wing_prompt->Fill(Ewing);
                    }

                    if (Ewing > 0.03)
                    {
                        candidate_events_ewing3_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                        candidate_events_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing3_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.01)
                    {
                        candidate_events_ewing1_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                        candidate_events_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing1_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    //---------------------------more ewing cuts-------------------------------

                    if (Ewing > 0.02)
                    {
                        candidate_events_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing2_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.025)
                    {
                        candidate_events_ewing25->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing25_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.031)
                    {
                        candidate_events_ewing31->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing31_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.032)
                    {
                        candidate_events_ewing32->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing32_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.033)
                    {
                        candidate_events_ewing33->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing33_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.034)
                    {
                        candidate_events_ewing34->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing34_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.035)
                    {
                        candidate_events_ewing35->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing35_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.036)
                    {
                        candidate_events_ewing36->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing36_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.037)
                    {
                        candidate_events_ewing37->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing37_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.038)
                    {
                        candidate_events_ewing38->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing38_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.039)
                    {
                        candidate_events_ewing39->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing39_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }
                    

                    if (Ewing > 0.040)
                    {
                        candidate_events_ewing4->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing4_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.041)
                    {
                        candidate_events_ewing41->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing41_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.042)
                    {
                        candidate_events_ewing42->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing42_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.043)
                    {
                        candidate_events_ewing43->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing43_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.044)
                    {
                        candidate_events_ewing44->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing44_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.045)
                    {
                        candidate_events_ewing45->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing45_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.046)
                    {
                        candidate_events_ewing46->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing46_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.047)
                    {
                        candidate_events_ewing47->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing47_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.048)
                    {
                        candidate_events_ewing48->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing48_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.049)
                    {
                        candidate_events_ewing49->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing49_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.05)
                    {
                        candidate_events_ewing5->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing5_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.07)
                    {
                        candidate_events_ewing7->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing7_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.09)
                    {
                        candidate_events_ewing9->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing9_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.10)
                    {
                        candidate_events_ewing10->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing10_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.12)
                    {
                        candidate_events_ewing12->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing12_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.14)
                    {
                        candidate_events_ewing14->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing14_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.15)
                    {
                        candidate_events_ewing15->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing15_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if ((*phoSeedTime)[mitIDXSEL] < -12.5)
                    {
                        spikeincandidate_eta_wing_leftonly3->Fill(Ewing);
                        spikeincandidate_eta_wing3->Fill(Ewing);
                        spikeincandidate_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spikeincandidate_events_ewing3_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                            spikeincandidate_events_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.01)
                        {
                            spikeincandidate_events_ewing1_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                            spikeincandidate_events_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        //---------------------------more ewing cuts-------------------------------

                        if (Ewing > 0.02)
                        {
                            spikeincandidate_events_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.025)
                        {
                            spikeincandidate_events_ewing25->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.031)
                        {
                            spikeincandidate_events_ewing31->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.032)
                        {
                            spikeincandidate_events_ewing32->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.033)
                        {
                            spikeincandidate_events_ewing33->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.034)
                        {
                            spikeincandidate_events_ewing34->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.035)
                        {
                            spikeincandidate_events_ewing35->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.036)
                        {
                            spikeincandidate_events_ewing36->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.037)
                        {
                            spikeincandidate_events_ewing37->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.038)
                        {
                            spikeincandidate_events_ewing38->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.039)
                        {
                            spikeincandidate_events_ewing39->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                        

                        if (Ewing > 0.040)
                        {
                            spikeincandidate_events_ewing4->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.041)
                        {
                            spikeincandidate_events_ewing41->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.042)
                        {
                            spikeincandidate_events_ewing42->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.043)
                        {
                            spikeincandidate_events_ewing43->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.044)
                        {
                            spikeincandidate_events_ewing44->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.045)
                        {
                            spikeincandidate_events_ewing45->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.046)
                        {
                            spikeincandidate_events_ewing46->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.047)
                        {
                            spikeincandidate_events_ewing47->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.048)
                        {
                            spikeincandidate_events_ewing48->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.049)
                        {
                            spikeincandidate_events_ewing49->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.05)
                        {
                            spikeincandidate_events_ewing5->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.07)
                        {
                            spikeincandidate_events_ewing7->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.09)
                        {
                            spikeincandidate_events_ewing9->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.10)
                        {
                            spikeincandidate_events_ewing10->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.12)
                        {
                            spikeincandidate_events_ewing12->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.14)
                        {
                            spikeincandidate_events_ewing14->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.15)
                        {
                            spikeincandidate_events_ewing15->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing < 0.03)
                    {
                        file13 << "--------cand case 3--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }

                    if (Ewing < 0.01)
                    {
                        file14 << "--------cand case 3--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }
                }

                //if there are both
                if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                {
                    rightcellsIDX = cellsEB_right[0];
                    leftcellsIDX = cellsEB_left[0];

                    if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[rightcellsIDX];
                    }
                    else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[leftcellsIDX];
                    }

                    Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);

                    file2 << "-------------Cand3_both-------------" << "\n"
                    << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                    << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                    << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                    << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                    << run << " " << lumis << " " << event << "\n"
                    << endl;

                    candidate_eta_wing_both3->Fill(Ewing);
                    candidate_eta_wing3->Fill(Ewing);
                    candidate_eta_wing->Fill(Ewing);

                    if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                    {
                        candidate_eta_wing_prompt->Fill(Ewing);
                    }

                    if (Ewing > 0.03)
                    {
                        candidate_events_ewing3_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                        candidate_events_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing3_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.01)
                    {
                        candidate_events_ewing1_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                        candidate_events_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing1_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    //---------------------------more ewing cuts-------------------------------

                    if (Ewing > 0.02)
                    {
                        candidate_events_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing2_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.025)
                    {
                        candidate_events_ewing25->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing25_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.031)
                    {
                        candidate_events_ewing31->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing31_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.032)
                    {
                        candidate_events_ewing32->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing32_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.033)
                    {
                        candidate_events_ewing33->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing33_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.034)
                    {
                        candidate_events_ewing34->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing34_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.035)
                    {
                        candidate_events_ewing35->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing35_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.036)
                    {
                        candidate_events_ewing36->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing36_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.037)
                    {
                        candidate_events_ewing37->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing37_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.038)
                    {
                        candidate_events_ewing38->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing38_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.039)
                    {
                        candidate_events_ewing39->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing39_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }
                    

                    if (Ewing > 0.040)
                    {
                        candidate_events_ewing4->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing4_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.041)
                    {
                        candidate_events_ewing41->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing41_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.042)
                    {
                        candidate_events_ewing42->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing42_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.043)
                    {
                        candidate_events_ewing43->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing43_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.044)
                    {
                        candidate_events_ewing44->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing44_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.045)
                    {
                        candidate_events_ewing45->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing45_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.046)
                    {
                        candidate_events_ewing46->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing46_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.047)
                    {
                        candidate_events_ewing47->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing47_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.048)
                    {
                        candidate_events_ewing48->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing48_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.049)
                    {
                        candidate_events_ewing49->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing49_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.05)
                    {
                        candidate_events_ewing5->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing5_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.07)
                    {
                        candidate_events_ewing7->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing7_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.09)
                    {
                        candidate_events_ewing9->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing9_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.10)
                    {
                        candidate_events_ewing10->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing10_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.12)
                    {
                        candidate_events_ewing12->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing12_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.14)
                    {
                        candidate_events_ewing14->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing14_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing > 0.15)
                    {
                        candidate_events_ewing15->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing15_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if ((*phoSeedTime)[mitIDXSEL] < -12.5)
                    {
                        spikeincandidate_eta_wing_both3->Fill(Ewing);
                        spikeincandidate_eta_wing3->Fill(Ewing);
                        spikeincandidate_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spikeincandidate_events_ewing3_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                            spikeincandidate_events_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.01)
                        {
                            spikeincandidate_events_ewing1_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                            spikeincandidate_events_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        //---------------------------more ewing cuts-------------------------------

                        if (Ewing > 0.02)
                        {
                            spikeincandidate_events_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.025)
                        {
                            spikeincandidate_events_ewing25->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.031)
                        {
                            spikeincandidate_events_ewing31->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.032)
                        {
                            spikeincandidate_events_ewing32->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.033)
                        {
                            spikeincandidate_events_ewing33->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.034)
                        {
                            spikeincandidate_events_ewing34->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.035)
                        {
                            spikeincandidate_events_ewing35->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.036)
                        {
                            spikeincandidate_events_ewing36->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.037)
                        {
                            spikeincandidate_events_ewing37->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.038)
                        {
                            spikeincandidate_events_ewing38->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.039)
                        {
                            spikeincandidate_events_ewing39->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                        

                        if (Ewing > 0.040)
                        {
                            spikeincandidate_events_ewing4->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.041)
                        {
                            spikeincandidate_events_ewing41->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.042)
                        {
                            spikeincandidate_events_ewing42->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.043)
                        {
                            spikeincandidate_events_ewing43->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.044)
                        {
                            spikeincandidate_events_ewing44->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.045)
                        {
                            spikeincandidate_events_ewing45->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.046)
                        {
                            spikeincandidate_events_ewing46->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.047)
                        {
                            spikeincandidate_events_ewing47->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.048)
                        {
                            spikeincandidate_events_ewing48->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.049)
                        {
                            spikeincandidate_events_ewing49->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.05)
                        {
                            spikeincandidate_events_ewing5->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.07)
                        {
                            spikeincandidate_events_ewing7->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.09)
                        {
                            spikeincandidate_events_ewing9->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.10)
                        {
                            spikeincandidate_events_ewing10->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.12)
                        {
                            spikeincandidate_events_ewing12->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.14)
                        {
                            spikeincandidate_events_ewing14->Fill((*phoSeedTime)[mitIDXSEL]);
                        }

                        if (Ewing > 0.15)
                        {
                            spikeincandidate_events_ewing15->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
                    }

                    if (Ewing < 0.03)
                    {
                        file13 << "--------cand case 3--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }

                    if (Ewing < 0.01)
                    {
                        file14 << "--------cand case 3--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }
                }
            }


            //-----------------promptZ cases 3, 4 ,5-----------------
            if (iiZpho.size() != 0 && ooZpho.size() != 0)
            {
                mootZIDXSEL = ootitZmatched[0];
		    	mitZIDXSEL = itootZmatched[0];
                //Prompt Z template case 3
                if ((*phohasPixelSeed)[mitIDXSEL] == 1 &&
                    (*phoMIPTotEnergy)[mitIDXSEL] < 4.9 &&
                    (*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.001 &&
                    (*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] > 0.001 &&
                    (*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] < 0.01015 &&
                    metFilters == 0 &&
                    (*phohasPixelSeed)[mitZIDXSEL] == 1 &&
                    (*phoMIPTotEnergy)[mitZIDXSEL] < 4.9 &&
                    (*phoSigmaIEtaIEtaFull5x5)[mitZIDXSEL] > 0.001 &&
                    (*phoSigmaIPhiIPhiFull5x5)[mitZIDXSEL] > 0.001 &&
                    (*phoSigmaIEtaIEtaFull5x5)[mitZIDXSEL] < 0.01015)
                {
                    Float_t InvM = InvariMass((*phoEt)[mitIDXSEL], (*phoEt)[mitZIDXSEL], (*phoPhi)[mitIDXSEL], (*phoPhi)[mitZIDXSEL], (*phoEta)[mitIDXSEL], (*phoEta)[mitZIDXSEL]);
                    InvMass_Z->Fill(InvM);

                    if (InvM > 85 && InvM < 100)
                    {
                        promptZ_events_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                        promptZ_events->Fill((*phoSeedTime)[mitIDXSEL]);

                        if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                        {
                            promptZ_events_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                        }
  
                        //promptZ3 searching adjencent cells
                        for (int icell = 0; icell < nAllCellsEB; icell++)
                        {                        
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mitIDXSEL] + 1)
                            {
                                cellsEB_right.push_back(icell); 

                                //if there's any err
                                if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mitIDXSEL]) == 86 || (*phoSeedIEta)[mitIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                                {
                                    file7 << "------------------promptZ3_err_right---------------------" << "\n"
                                    << "phoSeediphi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mitIDXSEL] << "\n"
                                    << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                                }
                            }

                            //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && (*phoSeedIEta)[mitIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                            {
                                cellsEB_right.push_back(icell); 

                                file9 << "-------------promptZ3_ietamid_right-------------" << "\n"
                                << "phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                                << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                                << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                                << run << " " << lumis << " " << event << "\n"
                                << endl;
                            }

                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mitIDXSEL] - 1)
                            {                    
                                cellsEB_left.push_back(icell);

                                //if there's any err
                                if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mitIDXSEL]) == 86 || (*phoSeedIEta)[mitIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                                {
                                    file7 << "------------------promptZ3_err_left---------------------" << "\n"
                                    << "phoSeediphi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mitIDXSEL] << "\n"
                                    << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                                }
                            }            

                            //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && (*phoSeedIEta)[mitIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                            {
                                cellsEB_left.push_back(icell); 

                                file9 << "-------------promptZ3_ietamid_left-------------" << "\n"
                                << "phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                                << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                                << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                                << run << " " << lumis << " " << event << "\n"
                                << endl;
                            }  
                        }

                        //Now I have to take if there's only one and both
                        //If there's only right
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                        {
                            rightcellsIDX_only = cellsEB_right[0];
                            maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);

                            file3 << "-------------promptZ3_rightonly-------------" << "\n"
                            << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;

                            promptZ_eta_wing_rightonly3->Fill(Ewing);
                            promptZ_eta_wing3->Fill(Ewing);
                            promptZ_eta_wing->Fill(Ewing);

                            if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                            {
                                promptZ_eta_wing_prompt->Fill(Ewing);
                            }

                            if (Ewing > 0.03)
                            {
                                promptZ_events_ewing3_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                                promptZ_events_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing3_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.01)
                            {
                                promptZ_events_ewing1_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                                promptZ_events_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing1_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            //---------------------------more ewing cuts-------------------------------

                            if (Ewing > 0.02)
                            {
                                promptZ_events_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing2_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.025)
                            {
                                promptZ_events_ewing25->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing25_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.031)
                            {
                                promptZ_events_ewing31->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing31_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.032)
                            {
                                promptZ_events_ewing32->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing32_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.033)
                            {
                                promptZ_events_ewing33->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing33_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.034)
                            {
                                promptZ_events_ewing34->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing34_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.035)
                            {
                                promptZ_events_ewing35->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing35_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.036)
                            {
                                promptZ_events_ewing36->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing36_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.037)
                            {
                                promptZ_events_ewing37->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing37_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.038)
                            {
                                promptZ_events_ewing38->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing38_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.039)
                            {
                                promptZ_events_ewing39->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing39_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }
                            

                            if (Ewing > 0.040)
                            {
                                promptZ_events_ewing4->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing4_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.041)
                            {
                                promptZ_events_ewing41->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing41_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.042)
                            {
                                promptZ_events_ewing42->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing42_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.043)
                            {
                                promptZ_events_ewing43->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing43_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.044)
                            {
                                promptZ_events_ewing44->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing44_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.045)
                            {
                                promptZ_events_ewing45->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing45_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.046)
                            {
                                promptZ_events_ewing46->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing46_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.047)
                            {
                                promptZ_events_ewing47->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing47_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.048)
                            {
                                promptZ_events_ewing48->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing48_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.049)
                            {
                                promptZ_events_ewing49->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing49_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.05)
                            {
                                promptZ_events_ewing5->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing5_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.07)
                            {
                                promptZ_events_ewing7->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing7_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.09)
                            {
                                promptZ_events_ewing9->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing9_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.10)
                            {
                                promptZ_events_ewing10->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing10_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.12)
                            {
                                promptZ_events_ewing12->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing12_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.14)
                            {
                                promptZ_events_ewing14->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing14_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.15)
                            {
                                promptZ_events_ewing15->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing15_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }
                        }

                        //if there is only left
                        if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                        {
                            leftcellsIDX_only = cellsEB_left[0];
                            maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);

                            file3 << "-------------promptZ3_leftonly-------------" << "\n"
                            << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;

                            promptZ_eta_wing_leftonly3->Fill(Ewing);
                            promptZ_eta_wing3->Fill(Ewing);
                            promptZ_eta_wing->Fill(Ewing);

                            if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                            {
                                promptZ_eta_wing_prompt->Fill(Ewing);
                            }

                            if (Ewing > 0.03)
                            {
                                promptZ_events_ewing3_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                                promptZ_events_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing3_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.01)
                            {
                                promptZ_events_ewing1_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                                promptZ_events_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing1_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            //---------------------------more ewing cuts-------------------------------

                            if (Ewing > 0.02)
                            {
                                promptZ_events_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing2_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.025)
                            {
                                promptZ_events_ewing25->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing25_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.031)
                            {
                                promptZ_events_ewing31->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing31_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.032)
                            {
                                promptZ_events_ewing32->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing32_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.033)
                            {
                                promptZ_events_ewing33->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing33_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.034)
                            {
                                promptZ_events_ewing34->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing34_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.035)
                            {
                                promptZ_events_ewing35->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing35_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.036)
                            {
                                promptZ_events_ewing36->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing36_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.037)
                            {
                                promptZ_events_ewing37->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing37_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.038)
                            {
                                promptZ_events_ewing38->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing38_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.039)
                            {
                                promptZ_events_ewing39->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing39_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }
                            

                            if (Ewing > 0.040)
                            {
                                promptZ_events_ewing4->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing4_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.041)
                            {
                                promptZ_events_ewing41->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing41_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.042)
                            {
                                promptZ_events_ewing42->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing42_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.043)
                            {
                                promptZ_events_ewing43->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing43_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.044)
                            {
                                promptZ_events_ewing44->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing44_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.045)
                            {
                                promptZ_events_ewing45->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing45_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.046)
                            {
                                promptZ_events_ewing46->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing46_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.047)
                            {
                                promptZ_events_ewing47->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing47_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.048)
                            {
                                promptZ_events_ewing48->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing48_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.049)
                            {
                                promptZ_events_ewing49->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing49_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.05)
                            {
                                promptZ_events_ewing5->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing5_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.07)
                            {
                                promptZ_events_ewing7->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing7_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.09)
                            {
                                promptZ_events_ewing9->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing9_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.10)
                            {
                                promptZ_events_ewing10->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing10_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.12)
                            {
                                promptZ_events_ewing12->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing12_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.14)
                            {
                                promptZ_events_ewing14->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing14_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.15)
                            {
                                promptZ_events_ewing15->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing15_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }
                        }

                        //if there are both
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                        {
                            rightcellsIDX = cellsEB_right[0];
                            leftcellsIDX = cellsEB_left[0];

                            if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[rightcellsIDX];
                            }
                            else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[leftcellsIDX];
                            }

                            Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);

                            file3 << "-------------promptZ3_both-------------" << "\n"
                            << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;

                            promptZ_eta_wing_both3->Fill(Ewing);
                            promptZ_eta_wing3->Fill(Ewing);
                            promptZ_eta_wing->Fill(Ewing);

                            if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                            {
                                promptZ_eta_wing_prompt->Fill(Ewing);
                            }

                            if (Ewing > 0.03)
                            {
                                promptZ_events_ewing3_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                                promptZ_events_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing3_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.01)
                            {
                                promptZ_events_ewing1_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                                promptZ_events_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing1_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            //---------------------------more ewing cuts-------------------------------

                            if (Ewing > 0.02)
                            {
                                promptZ_events_ewing2->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing2_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.025)
                            {
                                promptZ_events_ewing25->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing25_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.031)
                            {
                                promptZ_events_ewing31->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing31_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.032)
                            {
                                promptZ_events_ewing32->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing32_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.033)
                            {
                                promptZ_events_ewing33->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing33_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.034)
                            {
                                promptZ_events_ewing34->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing34_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.035)
                            {
                                promptZ_events_ewing35->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing35_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.036)
                            {
                                promptZ_events_ewing36->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing36_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.037)
                            {
                                promptZ_events_ewing37->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing37_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.038)
                            {
                                promptZ_events_ewing38->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing38_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.039)
                            {
                                promptZ_events_ewing39->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing39_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }
                            

                            if (Ewing > 0.040)
                            {
                                promptZ_events_ewing4->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing4_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.041)
                            {
                                promptZ_events_ewing41->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing41_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.042)
                            {
                                promptZ_events_ewing42->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing42_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.043)
                            {
                                promptZ_events_ewing43->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing43_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.044)
                            {
                                promptZ_events_ewing44->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing44_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.045)
                            {
                                promptZ_events_ewing45->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing45_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.046)
                            {
                                promptZ_events_ewing46->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing46_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.047)
                            {
                                promptZ_events_ewing47->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing47_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.048)
                            {
                                promptZ_events_ewing48->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing48_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.049)
                            {
                                promptZ_events_ewing49->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing49_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.05)
                            {
                                promptZ_events_ewing5->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing5_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.07)
                            {
                                promptZ_events_ewing7->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing7_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.09)
                            {
                                promptZ_events_ewing9->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing9_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.10)
                            {
                                promptZ_events_ewing10->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing10_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.12)
                            {
                                promptZ_events_ewing12->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing12_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.14)
                            {
                                promptZ_events_ewing14->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing14_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }

                            if (Ewing > 0.15)
                            {
                                promptZ_events_ewing15->Fill((*phoSeedTime)[mitIDXSEL]);

                                if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing15_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                                }
                            }
                        } 
                    }
                }

                //Prompt Z template case 4
                if ((*ophohasPixelSeed)[mootIDXSEL] == 1 &&
                    (*ophoMIPTotEnergy)[mootIDXSEL] < 4.9 &&
                    (*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 &&
                    (*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001 &&
                    (*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] < 0.01015 &&
                    metFilters == 0 &&
                    (*ophohasPixelSeed)[mootZIDXSEL] == 1 &&
                    (*ophoMIPTotEnergy)[mootZIDXSEL] < 4.9 &&
                    (*ophoSigmaIEtaIEtaFull5x5)[mootZIDXSEL] > 0.001 &&
                    (*ophoSigmaIPhiIPhiFull5x5)[mootZIDXSEL] > 0.001 &&
                    (*ophoSigmaIEtaIEtaFull5x5)[mootZIDXSEL] < 0.01015)
                {
                    Float_t InvM = InvariMass((*ophoEt)[mootIDXSEL], (*ophoEt)[mootZIDXSEL], (*ophoPhi)[mootIDXSEL], (*ophoPhi)[mootZIDXSEL], (*ophoEta)[mootIDXSEL], (*ophoEta)[mootZIDXSEL]);
                    InvMass_Z->Fill(InvM);

                    if (InvM > 85 && InvM < 100)
                    {
                        promptZ_events_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                        promptZ_events->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            promptZ_events_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
  
                        //promptZ1 searching adjencent cells
                        for (int icell = 0; icell < nAllCellsEB; icell++)
                        {                        
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mootIDXSEL] + 1)
                            {
                                cellsEB_right.push_back(icell); 

                                //if there's any err
                                if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mootIDXSEL]) == 86 || (*phoSeedIEta)[mootIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                                {
                                    file7 << "------------------promptZ4_err_right---------------------" << "\n"
                                    << "phoSeediphi = " << (*phoSeedIPhi)[mootIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mootIDXSEL] << "\n"
                                    << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                                }
                            }

                            //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mootIDXSEL] && (*phoSeedIEta)[mootIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                            {
                                cellsEB_right.push_back(icell); 

                                file9 << "-------------promptZ4_ietamid_right-------------" << "\n"
                                << "phoseedTime = " << (*phoSeedTime)[mootIDXSEL] << "\n"
                                << "phoSeedIPhi = " << (*phoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mootIDXSEL] << "\n" 
                                << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                                << run << " " << lumis << " " << event << "\n"
                                << endl;
                            }

                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mootIDXSEL] - 1)
                            {                    
                                cellsEB_left.push_back(icell);

                                //if there's any err
                                if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mootIDXSEL]) == 86 || (*phoSeedIEta)[mootIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                                {
                                    file7 << "------------------promptZ4_err_left---------------------" << "\n"
                                    << "phoSeediphi = " << (*phoSeedIPhi)[mootIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mootIDXSEL] << "\n"
                                    << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                                }
                            }            

                            //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mootIDXSEL] && (*phoSeedIEta)[mootIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                            {
                                cellsEB_left.push_back(icell); 

                                file9 << "-------------promptZ4_ietamid_left-------------" << "\n"
                                << "phoseedTime = " << (*phoSeedTime)[mootIDXSEL] << "\n"
                                << "phoSeedIPhi = " << (*phoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mootIDXSEL] << "\n" 
                                << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                                << run << " " << lumis << " " << event << "\n"
                                << endl;
                            }  
                        }

                        //Now I have to take if there's only one and both
                        //If there's only right
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                        {
                            rightcellsIDX_only = cellsEB_right[0];
                            maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                            file3 << "-------------promptZ4_rightonly-------------" << "\n"
                            << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;

                            promptZ_eta_wing_rightonly4->Fill(Ewing);
                            promptZ_eta_wing4->Fill(Ewing);
                            promptZ_eta_wing->Fill(Ewing);

                            if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                            {
                                promptZ_eta_wing_prompt->Fill(Ewing);
                            }

                            if (Ewing > 0.03)
                            {
                                promptZ_events_ewing3_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                                promptZ_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing3_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.01)
                            {
                                promptZ_events_ewing1_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                                promptZ_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing1_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            //---------------------------more ewing cuts-------------------------------

                            if (Ewing > 0.02)
                            {
                                promptZ_events_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing2_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.025)
                            {
                                promptZ_events_ewing25->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing25_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.031)
                            {
                                promptZ_events_ewing31->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing31_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.032)
                            {
                                promptZ_events_ewing32->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing32_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.033)
                            {
                                promptZ_events_ewing33->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing33_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.034)
                            {
                                promptZ_events_ewing34->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing34_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.035)
                            {
                                promptZ_events_ewing35->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing35_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.036)
                            {
                                promptZ_events_ewing36->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing36_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.037)
                            {
                                promptZ_events_ewing37->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing37_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.038)
                            {
                                promptZ_events_ewing38->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing38_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.039)
                            {
                                promptZ_events_ewing39->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing39_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }
                            

                            if (Ewing > 0.040)
                            {
                                promptZ_events_ewing4->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing4_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.041)
                            {
                                promptZ_events_ewing41->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing41_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.042)
                            {
                                promptZ_events_ewing42->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing42_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.043)
                            {
                                promptZ_events_ewing43->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing43_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.044)
                            {
                                promptZ_events_ewing44->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing44_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.045)
                            {
                                promptZ_events_ewing45->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing45_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.046)
                            {
                                promptZ_events_ewing46->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing46_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.047)
                            {
                                promptZ_events_ewing47->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing47_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.048)
                            {
                                promptZ_events_ewing48->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing48_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.049)
                            {
                                promptZ_events_ewing49->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing49_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.05)
                            {
                                promptZ_events_ewing5->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing5_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.07)
                            {
                                promptZ_events_ewing7->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing7_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.09)
                            {
                                promptZ_events_ewing9->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing9_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.10)
                            {
                                promptZ_events_ewing10->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing10_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.12)
                            {
                                promptZ_events_ewing12->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing12_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.14)
                            {
                                promptZ_events_ewing14->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing14_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.15)
                            {
                                promptZ_events_ewing15->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing15_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }
                        }

                        //if there is only left
                        if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                        {
                            leftcellsIDX_only = cellsEB_left[0];
                            maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                            file3 << "-------------promptZ4_leftonly-------------" << "\n"
                            << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;

                            promptZ_eta_wing_leftonly4->Fill(Ewing);
                            promptZ_eta_wing4->Fill(Ewing);
                            promptZ_eta_wing->Fill(Ewing);

                            if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                            {
                                promptZ_eta_wing_prompt->Fill(Ewing);
                            }

                            if (Ewing > 0.03)
                            {
                                promptZ_events_ewing3_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                                promptZ_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing3_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.01)
                            {
                                promptZ_events_ewing1_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                                promptZ_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing1_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            //---------------------------more ewing cuts-------------------------------

                            if (Ewing > 0.02)
                            {
                                promptZ_events_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing2_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.025)
                            {
                                promptZ_events_ewing25->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing25_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.031)
                            {
                                promptZ_events_ewing31->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing31_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.032)
                            {
                                promptZ_events_ewing32->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing32_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.033)
                            {
                                promptZ_events_ewing33->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing33_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.034)
                            {
                                promptZ_events_ewing34->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing34_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.035)
                            {
                                promptZ_events_ewing35->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing35_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.036)
                            {
                                promptZ_events_ewing36->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing36_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.037)
                            {
                                promptZ_events_ewing37->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing37_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.038)
                            {
                                promptZ_events_ewing38->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing38_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.039)
                            {
                                promptZ_events_ewing39->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing39_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }
                            

                            if (Ewing > 0.040)
                            {
                                promptZ_events_ewing4->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing4_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.041)
                            {
                                promptZ_events_ewing41->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing41_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.042)
                            {
                                promptZ_events_ewing42->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing42_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.043)
                            {
                                promptZ_events_ewing43->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing43_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.044)
                            {
                                promptZ_events_ewing44->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing44_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.045)
                            {
                                promptZ_events_ewing45->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing45_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.046)
                            {
                                promptZ_events_ewing46->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing46_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.047)
                            {
                                promptZ_events_ewing47->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing47_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.048)
                            {
                                promptZ_events_ewing48->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing48_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.049)
                            {
                                promptZ_events_ewing49->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing49_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.05)
                            {
                                promptZ_events_ewing5->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing5_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.07)
                            {
                                promptZ_events_ewing7->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing7_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.09)
                            {
                                promptZ_events_ewing9->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing9_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.10)
                            {
                                promptZ_events_ewing10->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing10_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.12)
                            {
                                promptZ_events_ewing12->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing12_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.14)
                            {
                                promptZ_events_ewing14->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing14_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.15)
                            {
                                promptZ_events_ewing15->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing15_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }
                        }

                        //if there are both
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                        {
                            rightcellsIDX = cellsEB_right[0];
                            leftcellsIDX = cellsEB_left[0];

                            if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[rightcellsIDX];
                            }
                            else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[leftcellsIDX];
                            }

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                            file3 << "-------------promptZ4_both-------------" << "\n"
                            << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;

                            promptZ_eta_wing_both4->Fill(Ewing);
                            promptZ_eta_wing4->Fill(Ewing);
                            promptZ_eta_wing->Fill(Ewing);

                            if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                            {
                                promptZ_eta_wing_prompt->Fill(Ewing);
                            }

                            if (Ewing > 0.03)
                            {
                                promptZ_events_ewing3_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                                promptZ_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing3_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.01)
                            {
                                promptZ_events_ewing1_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                                promptZ_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing1_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            //---------------------------more ewing cuts-------------------------------

                            if (Ewing > 0.02)
                            {
                                promptZ_events_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing2_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.025)
                            {
                                promptZ_events_ewing25->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing25_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.031)
                            {
                                promptZ_events_ewing31->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing31_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.032)
                            {
                                promptZ_events_ewing32->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing32_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.033)
                            {
                                promptZ_events_ewing33->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing33_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.034)
                            {
                                promptZ_events_ewing34->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing34_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.035)
                            {
                                promptZ_events_ewing35->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing35_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.036)
                            {
                                promptZ_events_ewing36->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing36_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.037)
                            {
                                promptZ_events_ewing37->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing37_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.038)
                            {
                                promptZ_events_ewing38->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing38_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.039)
                            {
                                promptZ_events_ewing39->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing39_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }
                            

                            if (Ewing > 0.040)
                            {
                                promptZ_events_ewing4->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing4_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.041)
                            {
                                promptZ_events_ewing41->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing41_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.042)
                            {
                                promptZ_events_ewing42->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing42_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.043)
                            {
                                promptZ_events_ewing43->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing43_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.044)
                            {
                                promptZ_events_ewing44->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing44_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.045)
                            {
                                promptZ_events_ewing45->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing45_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.046)
                            {
                                promptZ_events_ewing46->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing46_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.047)
                            {
                                promptZ_events_ewing47->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing47_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.048)
                            {
                                promptZ_events_ewing48->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing48_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.049)
                            {
                                promptZ_events_ewing49->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing49_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.05)
                            {
                                promptZ_events_ewing5->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing5_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.07)
                            {
                                promptZ_events_ewing7->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing7_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.09)
                            {
                                promptZ_events_ewing9->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing9_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.10)
                            {
                                promptZ_events_ewing10->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing10_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.12)
                            {
                                promptZ_events_ewing12->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing12_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.14)
                            {
                                promptZ_events_ewing14->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing14_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.15)
                            {
                                promptZ_events_ewing15->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing15_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }
                        } 
                    }
                }

                //Prompt Z template case 5
                if ((*phohasPixelSeed)[mitIDXSEL] == 1 &&
                    (*phoMIPTotEnergy)[mitIDXSEL] < 4.9 &&
                    (*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.001 &&
                    (*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] > 0.001 &&
                    (*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] < 0.01015 &&
                    metFilters == 0 &&
                    (*phohasPixelSeed)[mitZIDXSEL] == 1 &&
                    (*phoMIPTotEnergy)[mitZIDXSEL] < 4.9 &&
                    (*phoSigmaIEtaIEtaFull5x5)[mitZIDXSEL] > 0.001 &&
                    (*phoSigmaIPhiIPhiFull5x5)[mitZIDXSEL] > 0.001 &&
                    (*phoSigmaIEtaIEtaFull5x5)[mitZIDXSEL] < 0.01015 &&

                    (*ophohasPixelSeed)[mootIDXSEL] == 1 &&
                    (*ophoMIPTotEnergy)[mootIDXSEL] < 4.9 &&
                    (*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 &&
                    (*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001 &&
                    (*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] < 0.01015 &&

                    (*ophohasPixelSeed)[mootZIDXSEL] == 1 &&
                    (*ophoMIPTotEnergy)[mootZIDXSEL] < 4.9 &&
                    (*ophoSigmaIEtaIEtaFull5x5)[mootZIDXSEL] > 0.001 &&
                    (*ophoSigmaIPhiIPhiFull5x5)[mootZIDXSEL] > 0.001 &&
                    (*ophoSigmaIEtaIEtaFull5x5)[mootZIDXSEL] < 0.01015)
                {
                    Float_t InvM = InvariMass((*ophoEt)[mootIDXSEL], (*ophoEt)[mootZIDXSEL], (*ophoPhi)[mootIDXSEL], (*ophoPhi)[mootZIDXSEL], (*ophoEta)[mootIDXSEL], (*ophoEta)[mootZIDXSEL]);
                    InvMass_Z->Fill(InvM);

                    if (InvM > 85 && InvM < 100)
                    {
                        promptZ_events_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                        promptZ_events->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            promptZ_events_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
  
                        //promptZ1 searching adjencent cells
                        for (int icell = 0; icell < nAllCellsEB; icell++)
                        {                        
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mootIDXSEL] + 1)
                            {
                                cellsEB_right.push_back(icell); 

                                //if there's any err
                                if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mootIDXSEL]) == 86 || (*phoSeedIEta)[mootIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                                {
                                    file7 << "------------------promptZ5_err_right---------------------" << "\n"
                                    << "phoSeediphi = " << (*phoSeedIPhi)[mootIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mootIDXSEL] << "\n"
                                    << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                                }
                            }

                            //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mootIDXSEL] && (*phoSeedIEta)[mootIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                            {
                                cellsEB_right.push_back(icell); 

                                file9 << "-------------promptZ5_ietamid_right-------------" << "\n"
                                << "phoseedTime = " << (*phoSeedTime)[mootIDXSEL] << "\n"
                                << "phoSeedIPhi = " << (*phoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mootIDXSEL] << "\n" 
                                << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                                << run << " " << lumis << " " << event << "\n"
                                << endl;
                            }

                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mootIDXSEL] - 1)
                            {                    
                                cellsEB_left.push_back(icell);

                                //if there's any err
                                if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mootIDXSEL]) == 86 || (*phoSeedIEta)[mootIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                                {
                                    file7 << "------------------promptZ5_err_left---------------------" << "\n"
                                    << "phoSeediphi = " << (*phoSeedIPhi)[mootIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mootIDXSEL] << "\n"
                                    << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                                }
                            }            

                            //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                            if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mootIDXSEL] && (*phoSeedIEta)[mootIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                            {
                                cellsEB_left.push_back(icell); 

                                file9 << "-------------promptZ5_ietamid_left-------------" << "\n"
                                << "phoseedTime = " << (*phoSeedTime)[mootIDXSEL] << "\n"
                                << "phoSeedIPhi = " << (*phoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mootIDXSEL] << "\n" 
                                << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                                << run << " " << lumis << " " << event << "\n"
                                << endl;
                            }  
                        }

                        //Now I have to take if there's only one and both
                        //If there's only right
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                        {
                            rightcellsIDX_only = cellsEB_right[0];
                            maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                            file3 << "-------------promptZ5_rightonly-------------" << "\n"
                            << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;

                            promptZ_eta_wing_rightonly5->Fill(Ewing);
                            promptZ_eta_wing5->Fill(Ewing);
                            promptZ_eta_wing->Fill(Ewing);

                            if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                            {
                                promptZ_eta_wing_prompt->Fill(Ewing);
                            }

                            if (Ewing > 0.03)
                            {
                                promptZ_events_ewing3_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                                promptZ_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing3_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.01)
                            {
                                promptZ_events_ewing1_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                                promptZ_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing1_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            //---------------------------more ewing cuts-------------------------------

                            if (Ewing > 0.02)
                            {
                                promptZ_events_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing2_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.025)
                            {
                                promptZ_events_ewing25->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing25_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.031)
                            {
                                promptZ_events_ewing31->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing31_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.032)
                            {
                                promptZ_events_ewing32->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing32_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.033)
                            {
                                promptZ_events_ewing33->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing33_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.034)
                            {
                                promptZ_events_ewing34->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing34_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.035)
                            {
                                promptZ_events_ewing35->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing35_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.036)
                            {
                                promptZ_events_ewing36->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing36_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.037)
                            {
                                promptZ_events_ewing37->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing37_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.038)
                            {
                                promptZ_events_ewing38->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing38_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.039)
                            {
                                promptZ_events_ewing39->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing39_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }
                            

                            if (Ewing > 0.040)
                            {
                                promptZ_events_ewing4->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing4_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.041)
                            {
                                promptZ_events_ewing41->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing41_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.042)
                            {
                                promptZ_events_ewing42->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing42_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.043)
                            {
                                promptZ_events_ewing43->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing43_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.044)
                            {
                                promptZ_events_ewing44->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing44_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.045)
                            {
                                promptZ_events_ewing45->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing45_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.046)
                            {
                                promptZ_events_ewing46->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing46_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.047)
                            {
                                promptZ_events_ewing47->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing47_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.048)
                            {
                                promptZ_events_ewing48->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing48_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.049)
                            {
                                promptZ_events_ewing49->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing49_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.05)
                            {
                                promptZ_events_ewing5->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing5_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.07)
                            {
                                promptZ_events_ewing7->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing7_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.09)
                            {
                                promptZ_events_ewing9->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing9_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.10)
                            {
                                promptZ_events_ewing10->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing10_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.12)
                            {
                                promptZ_events_ewing12->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing12_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.14)
                            {
                                promptZ_events_ewing14->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing14_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.15)
                            {
                                promptZ_events_ewing15->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing15_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }
                        }

                        //if there is only left
                        if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                        {
                            leftcellsIDX_only = cellsEB_left[0];
                            maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                            file3 << "-------------promptZ5_leftonly-------------" << "\n"
                            << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;

                            promptZ_eta_wing_leftonly5->Fill(Ewing);
                            promptZ_eta_wing5->Fill(Ewing);
                            promptZ_eta_wing->Fill(Ewing);

                            if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                            {
                                promptZ_eta_wing_prompt->Fill(Ewing);
                            }

                            if (Ewing > 0.03)
                            {
                                promptZ_events_ewing3_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                                promptZ_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing3_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.01)
                            {
                                promptZ_events_ewing1_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                                promptZ_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing1_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            //---------------------------more ewing cuts-------------------------------

                            if (Ewing > 0.02)
                            {
                                promptZ_events_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing2_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.025)
                            {
                                promptZ_events_ewing25->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing25_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.031)
                            {
                                promptZ_events_ewing31->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing31_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.032)
                            {
                                promptZ_events_ewing32->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing32_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.033)
                            {
                                promptZ_events_ewing33->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing33_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.034)
                            {
                                promptZ_events_ewing34->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing34_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.035)
                            {
                                promptZ_events_ewing35->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing35_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.036)
                            {
                                promptZ_events_ewing36->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing36_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.037)
                            {
                                promptZ_events_ewing37->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing37_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.038)
                            {
                                promptZ_events_ewing38->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing38_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.039)
                            {
                                promptZ_events_ewing39->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing39_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }
                            

                            if (Ewing > 0.040)
                            {
                                promptZ_events_ewing4->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing4_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.041)
                            {
                                promptZ_events_ewing41->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing41_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.042)
                            {
                                promptZ_events_ewing42->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing42_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.043)
                            {
                                promptZ_events_ewing43->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing43_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.044)
                            {
                                promptZ_events_ewing44->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing44_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.045)
                            {
                                promptZ_events_ewing45->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing45_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.046)
                            {
                                promptZ_events_ewing46->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing46_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.047)
                            {
                                promptZ_events_ewing47->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing47_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.048)
                            {
                                promptZ_events_ewing48->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing48_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.049)
                            {
                                promptZ_events_ewing49->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing49_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.05)
                            {
                                promptZ_events_ewing5->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing5_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.07)
                            {
                                promptZ_events_ewing7->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing7_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.09)
                            {
                                promptZ_events_ewing9->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing9_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.10)
                            {
                                promptZ_events_ewing10->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing10_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.12)
                            {
                                promptZ_events_ewing12->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing12_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.14)
                            {
                                promptZ_events_ewing14->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing14_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.15)
                            {
                                promptZ_events_ewing15->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing15_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }
                        }

                        //if there are both
                        if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                        {
                            rightcellsIDX = cellsEB_right[0];
                            leftcellsIDX = cellsEB_left[0];

                            if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[rightcellsIDX];
                            }
                            else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                            {
                                maxCellsE = AllCellsE_EB[leftcellsIDX];
                            }

                            Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                            file3 << "-------------promptZ5_both-------------" << "\n"
                            << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;

                            promptZ_eta_wing_both5->Fill(Ewing);
                            promptZ_eta_wing5->Fill(Ewing);
                            promptZ_eta_wing->Fill(Ewing);

                            if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                            {
                                promptZ_eta_wing_prompt->Fill(Ewing);
                            }

                            if (Ewing > 0.03)
                            {
                                promptZ_events_ewing3_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                                promptZ_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing3_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.01)
                            {
                                promptZ_events_ewing1_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                                promptZ_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing1_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            //---------------------------more ewing cuts-------------------------------

                            if (Ewing > 0.02)
                            {
                                promptZ_events_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing2_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.025)
                            {
                                promptZ_events_ewing25->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing25_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.031)
                            {
                                promptZ_events_ewing31->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing31_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.032)
                            {
                                promptZ_events_ewing32->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing32_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.033)
                            {
                                promptZ_events_ewing33->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing33_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.034)
                            {
                                promptZ_events_ewing34->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing34_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.035)
                            {
                                promptZ_events_ewing35->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing35_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.036)
                            {
                                promptZ_events_ewing36->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing36_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.037)
                            {
                                promptZ_events_ewing37->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing37_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.038)
                            {
                                promptZ_events_ewing38->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing38_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.039)
                            {
                                promptZ_events_ewing39->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing39_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }
                            

                            if (Ewing > 0.040)
                            {
                                promptZ_events_ewing4->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing4_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.041)
                            {
                                promptZ_events_ewing41->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing41_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.042)
                            {
                                promptZ_events_ewing42->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing42_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.043)
                            {
                                promptZ_events_ewing43->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing43_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.044)
                            {
                                promptZ_events_ewing44->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing44_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.045)
                            {
                                promptZ_events_ewing45->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing45_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.046)
                            {
                                promptZ_events_ewing46->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing46_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.047)
                            {
                                promptZ_events_ewing47->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing47_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.048)
                            {
                                promptZ_events_ewing48->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing48_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.049)
                            {
                                promptZ_events_ewing49->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing49_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.05)
                            {
                                promptZ_events_ewing5->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing5_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.07)
                            {
                                promptZ_events_ewing7->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing7_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.09)
                            {
                                promptZ_events_ewing9->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing9_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.10)
                            {
                                promptZ_events_ewing10->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing10_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.12)
                            {
                                promptZ_events_ewing12->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing12_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.14)
                            {
                                promptZ_events_ewing14->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing14_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }

                            if (Ewing > 0.15)
                            {
                                promptZ_events_ewing15->Fill((*ophoSeedTime)[mootIDXSEL]);

                                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                                {
                                    promptZ_events_ewing15_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                                }
                            }
                        } 
                    }
                }
            }

            //BeamHalo template
			if (pfMET > 210 &&
				(*phohasPixelSeed)[mitIDXSEL] == 0 &&
				(*phoMIPTotEnergy)[mitIDXSEL] > 4.9 &&
				(*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.001 &&
				(*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] > 0.001 &&
				metFilters == 8)
			{
				BeamHalo_template_case3->Fill((*phoSeedTime)[mitIDXSEL]);
				BeamHalo_template->Fill((*phoSeedTime)[mitIDXSEL]);
                
                sieie_BeamHalo->Fill((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL]);


                if ((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.0176)
                {
                    BeamHalo_template_sieieb176->Fill((*phoSeedTime)[mitIDXSEL]);
                }
                
                if ((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.0206)
                {
                    BeamHalo_template_sieieb206->Fill((*phoSeedTime)[mitIDXSEL]);
                }

                if ((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] < 0.0176)
                {
                    BeamHalo_template_sieies176->Fill((*phoSeedTime)[mitIDXSEL]);
                }

                if ((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] < 0.0206)
                {
                    BeamHalo_template_sieies206->Fill((*phoSeedTime)[mitIDXSEL]);
                }
            }

            //Spikes template
            if (pfMET > 210 &&
                (*phohasPixelSeed)[mitIDXSEL] == 0 &&
                (*phoMIPTotEnergy)[mitIDXSEL] < 4.9 &&
                ((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] < 0.001 ||
                (*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] < 0.001) &&
                metFilters == 0)
            {
                //original showershap cut
                if ((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] < 0.01015)
                {  
                    spike_events_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                    //there are some seed without energy and timing info
                    //spike_events->Fill((*phoSeedTime)[mitIDXSEL]);
                    
                    spike_events_etaphi_case3->Fill((*phoEta)[mitIDXSEL], (*phoPhi)[mitIDXSEL]);
                    spike_events_etaphi->Fill((*phoEta)[mitIDXSEL], (*phoPhi)[mitIDXSEL]);

                    if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                    {
                        //spike_events_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                    }

                    file11 << "-------------spike3_extra-------------" << "\n"
                    << "phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n"
                    << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << "\n" 
                    << run << " " << lumis << " " << event << "\n"
                    << endl;

                    //cand1 searching adjencent cells
                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {                    
                        if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mitIDXSEL] + 1)
                        {
                            cellsEB_right.push_back(icell); 

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mitIDXSEL]) == 86 || (*phoSeedIEta)[mitIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike3_err_right---------------------" << "\n"
                                << "phoSeediphi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mitIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }

                        //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && (*phoSeedIEta)[mitIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                        {
                            cellsEB_right.push_back(icell); 

                            file9 << "-------------spike3_ietamid_right-------------" << "\n"
                            << "phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }

                        if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mitIDXSEL] - 1)
                        {                    
                            cellsEB_left.push_back(icell);

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mitIDXSEL]) == 86 || (*phoSeedIEta)[mitIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike3_err_left---------------------" << "\n"
                                << "phoSeediphi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mitIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }            

                        //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && (*phoSeedIEta)[mitIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                        {
                            cellsEB_left.push_back(icell); 

                            file9 << "-------------spike3_ietamid_left-------------" << "\n"
                            << "phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }  
                    
                    }

                    //Now I have to take if there's only one and both
                    //If there's only right
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                    {
                        rightcellsIDX_only = cellsEB_right[0];
                        maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);

                        file4 << "-------------spike3_rightonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_rightonly3->Fill(Ewing);
                        spike_eta_wing3->Fill(Ewing);
                        spike_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                            spike_events_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);

                            if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                            {
                                spike_events_ewing3_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                            }
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                            spike_events_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                            if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                            {
                                spike_events_ewing1_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                            }
                        }
                    }

                    //if there is only left
                    if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                    {
                        leftcellsIDX_only = cellsEB_left[0];
                        maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);

                        file4 << "-------------spike3_leftonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_leftonly3->Fill(Ewing);
                        spike_eta_wing3->Fill(Ewing);
                        spike_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                            spike_events_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);

                            if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                            {
                                spike_events_ewing3_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                            }
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                            spike_events_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                            if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                            {
                                spike_events_ewing1_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                            }
                        }
                    }

                    //if there are both
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                    {
                        rightcellsIDX = cellsEB_right[0];
                        leftcellsIDX = cellsEB_left[0];

                        if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[rightcellsIDX];
                        }
                        else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[leftcellsIDX];
                        }

                        Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);

                        file4 << "-------------spike3_both-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_both3->Fill(Ewing);
                        spike_eta_wing3->Fill(Ewing);
                        spike_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                            spike_events_ewing3->Fill((*phoSeedTime)[mitIDXSEL]);

                            if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                            {
                                spike_events_ewing3_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                            }
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                            spike_events_ewing1->Fill((*phoSeedTime)[mitIDXSEL]);

                            if (fabs((*phoSeedTime)[mitIDXSEL]) < 3.0)
                            {
                                spike_events_ewing1_3ns->Fill((*phoSeedTime)[mitIDXSEL]);
                            }
                        }
                    }
                }

                //time cut no showershape cut
                if ((*phoSeedTime)[mitIDXSEL] < -12.5)
                {
                    spike_events_time_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                    spike_events_time->Fill((*phoSeedTime)[mitIDXSEL]);
                    
                    spike_events_time_etaphi_case3->Fill((*phoEta)[mitIDXSEL], (*phoPhi)[mitIDXSEL]);
                    spike_events_time_etaphi->Fill((*phoEta)[mitIDXSEL], (*phoPhi)[mitIDXSEL]);

                    file6 << run << " " << lumis << " " << event << endl;

                    //cand1 searching adjencent cells
                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {                    
                        if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mitIDXSEL] + 1)
                        {
                            cellsEB_right.push_back(icell); 

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mitIDXSEL]) == 86 || (*phoSeedIEta)[mitIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike3_time_err_right---------------------" << "\n"
                                << "phoSeediphi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mitIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }

                        //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && (*phoSeedIEta)[mitIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                        {
                            cellsEB_right.push_back(icell); 

                            file10 << "-------------spike3_time_ietamid_right-------------" << "\n"
                            << "phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }

                        if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && AllCellsIEtaEB[icell] == (*phoSeedIEta)[mitIDXSEL] - 1)
                        {                    
                            cellsEB_left.push_back(icell);

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*phoSeedIEta)[mitIDXSEL]) == 86 || (*phoSeedIEta)[mitIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike3_time_err_left---------------------" << "\n"
                                << "phoSeediphi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedieta = " << (*phoSeedIEta)[mitIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }            

                        //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*phoSeedIPhi)[mitIDXSEL] && (*phoSeedIEta)[mitIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                        {
                            cellsEB_left.push_back(icell); 

                            file10 << "-------------spike3_time_ietamid_left-------------" << "\n"
                            << "phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }  
                    
                    }

                    //Now I have to take if there's only one and both
                    //If there's only right
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                    {
                        rightcellsIDX_only = cellsEB_right[0];
                        maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);

                        file5 << "-------------spike3_time_rightonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_time_rightonly3->Fill(Ewing);
                        spike_eta_wing_time3->Fill(Ewing);
                        spike_eta_wing_time->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_time_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                            spike_events_ewing3_time->Fill((*phoSeedTime)[mitIDXSEL]);

                            file15 << "--------spike case 3--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_time_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                            spike_events_ewing1_time->Fill((*phoSeedTime)[mitIDXSEL]);

                            file16 << "--------spike case 3--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }
                    }

                    //if there is only left
                    if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                    {
                        leftcellsIDX_only = cellsEB_left[0];
                        maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);

                        file5 << "-------------spike3_time_leftonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_time_leftonly3->Fill(Ewing);
                        spike_eta_wing_time3->Fill(Ewing);
                        spike_eta_wing_time->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_time_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                            spike_events_ewing3_time->Fill((*phoSeedTime)[mitIDXSEL]);

                            file15 << "--------spike case 3--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_time_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                            spike_events_ewing1_time->Fill((*phoSeedTime)[mitIDXSEL]);

                            file16 << "--------spike case 3--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }
                    }

                    //if there are both
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                    {
                        rightcellsIDX = cellsEB_right[0];
                        leftcellsIDX = cellsEB_left[0];

                        if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[rightcellsIDX];
                        }
                        else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[leftcellsIDX];
                        }

                        Float_t Ewing = maxCellsE/((*phoSeedEnergy)[mitIDXSEL]);

                        file5 << "-------------spike3_time_both-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*phoSeedTime)[mitIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*phoSeedIPhi)[mitIDXSEL] << " phoSeedIEta = " << (*phoSeedIEta)[mitIDXSEL] << " phoSeedEnergy = " << (*phoSeedEnergy)[mitIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_time_both3->Fill(Ewing);
                        spike_eta_wing_time3->Fill(Ewing);
                        spike_eta_wing_time->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_time_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                            spike_events_ewing3_time->Fill((*phoSeedTime)[mitIDXSEL]);

                            file15 << "--------spike case 3--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_time_case3->Fill((*phoSeedTime)[mitIDXSEL]);
                            spike_events_ewing1_time->Fill((*phoSeedTime)[mitIDXSEL]);

                            file16 << "--------spike case 3--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }
                    }
                }
            }


            //Case 4, only OOT
            Float_t NewMET = newuMET(pfMET, pfMETPhi, (*ophoPhi)[mootIDXSEL], (*ophoEt)[mootIDXSEL]);

            //Candidate Events
            if (NewMET > 210 &&
                (*ophohasPixelSeed)[mootIDXSEL] == 0 &&
                (*ophoMIPTotEnergy)[mootIDXSEL] < 4.9 &&
                (*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 &&
                (*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001 &&
                (*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] < 0.01015 &&
                metFilters == 0)
            {
                candidate_events_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
				candidate_events->Fill((*ophoSeedTime)[mootIDXSEL]);

                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                {
                    candidate_events_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                }

                if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                {
                    spikeincandidate_events_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
				    spikeincandidate_events->Fill((*ophoSeedTime)[mootIDXSEL]);
                }

                //cand1 searching adjencent cells
                for (int icell = 0; icell < nAllCellsEB; icell++)
                {
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] + 1)
                    {
                        cellsEB_right.push_back(icell); 

                        //if there's any err
                        if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*ophoSeedIEta)[mootIDXSEL]) == 86 || (*ophoSeedIEta)[mootIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                        {
                            file7 << "------------------cand4_err_right---------------------" << "\n"
                            << "phoSeediphi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedieta = " << (*ophoSeedIEta)[mootIDXSEL] << "\n"
                            << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                        }
                    }

                    //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                    {
                        cellsEB_right.push_back(icell); 

                        file8 << "-------------Cand4_ietamid_right-------------" << "\n"
                        << "phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;
                    }

                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] - 1)
                    {                    
                        cellsEB_left.push_back(icell);

                        //if there's any err
                        if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*ophoSeedIEta)[mootIDXSEL]) == 86 || (*ophoSeedIEta)[mootIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                        {
                            file7 << "------------------cand4_err_left---------------------" << "\n"
                            << "phoSeediphi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedieta = " << (*ophoSeedIEta)[mootIDXSEL] << "\n"
                            << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                        }
                    }            

                    //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                    {
                        cellsEB_left.push_back(icell); 

                        file8 << "-------------Cand4_ietamid_left-------------" << "\n"
                        << "phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;
                    }  
                
                }

                //Now I have to take if there's only one and both
                //If there's only right
                if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                {
                    rightcellsIDX_only = cellsEB_right[0];
                    maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                    file2 << "-------------Cand4_rightonly-------------" << "\n"
                    << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                    << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                    << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                    << run << " " << lumis << " " << event << "\n"
                    << endl;

                    candidate_eta_wing_rightonly4->Fill(Ewing);
                    candidate_eta_wing4->Fill(Ewing);
                    candidate_eta_wing->Fill(Ewing);

                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                    {
                        candidate_eta_wing_prompt->Fill(Ewing);
                    }

                    if (Ewing > 0.03)
                    {
                        candidate_events_ewing3_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                        candidate_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing3_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.01)
                    {
                        candidate_events_ewing1_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                        candidate_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing1_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    //---------------------------more ewing cuts-------------------------------

                    if (Ewing > 0.02)
                    {
                        candidate_events_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing2_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.025)
                    {
                        candidate_events_ewing25->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing25_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.031)
                    {
                        candidate_events_ewing31->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing31_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.032)
                    {
                        candidate_events_ewing32->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing32_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.033)
                    {
                        candidate_events_ewing33->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing33_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.034)
                    {
                        candidate_events_ewing34->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing34_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.035)
                    {
                        candidate_events_ewing35->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing35_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.036)
                    {
                        candidate_events_ewing36->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing36_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.037)
                    {
                        candidate_events_ewing37->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing37_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.038)
                    {
                        candidate_events_ewing38->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing38_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.039)
                    {
                        candidate_events_ewing39->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing39_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }
                    

                    if (Ewing > 0.040)
                    {
                        candidate_events_ewing4->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing4_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.041)
                    {
                        candidate_events_ewing41->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing41_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.042)
                    {
                        candidate_events_ewing42->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing42_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.043)
                    {
                        candidate_events_ewing43->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing43_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.044)
                    {
                        candidate_events_ewing44->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing44_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.045)
                    {
                        candidate_events_ewing45->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing45_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.046)
                    {
                        candidate_events_ewing46->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing46_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.047)
                    {
                        candidate_events_ewing47->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing47_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.048)
                    {
                        candidate_events_ewing48->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing48_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.049)
                    {
                        candidate_events_ewing49->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing49_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.05)
                    {
                        candidate_events_ewing5->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing5_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.07)
                    {
                        candidate_events_ewing7->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing7_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.09)
                    {
                        candidate_events_ewing9->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing9_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.10)
                    {
                        candidate_events_ewing10->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing10_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.12)
                    {
                        candidate_events_ewing12->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing12_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.14)
                    {
                        candidate_events_ewing14->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing14_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.15)
                    {
                        candidate_events_ewing15->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing15_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                    {
                        spikeincandidate_eta_wing_rightonly4->Fill(Ewing);
                        spikeincandidate_eta_wing4->Fill(Ewing);
                        spikeincandidate_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spikeincandidate_events_ewing3_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spikeincandidate_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.01)
                        {
                            spikeincandidate_events_ewing1_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spikeincandidate_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        //---------------------------more ewing cuts-------------------------------

                        if (Ewing > 0.02)
                        {
                            spikeincandidate_events_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.025)
                        {
                            spikeincandidate_events_ewing25->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.031)
                        {
                            spikeincandidate_events_ewing31->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.032)
                        {
                            spikeincandidate_events_ewing32->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.033)
                        {
                            spikeincandidate_events_ewing33->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.034)
                        {
                            spikeincandidate_events_ewing34->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.035)
                        {
                            spikeincandidate_events_ewing35->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.036)
                        {
                            spikeincandidate_events_ewing36->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.037)
                        {
                            spikeincandidate_events_ewing37->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.038)
                        {
                            spikeincandidate_events_ewing38->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.039)
                        {
                            spikeincandidate_events_ewing39->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                        

                        if (Ewing > 0.040)
                        {
                            spikeincandidate_events_ewing4->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.041)
                        {
                            spikeincandidate_events_ewing41->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.042)
                        {
                            spikeincandidate_events_ewing42->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.043)
                        {
                            spikeincandidate_events_ewing43->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.044)
                        {
                            spikeincandidate_events_ewing44->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.045)
                        {
                            spikeincandidate_events_ewing45->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.046)
                        {
                            spikeincandidate_events_ewing46->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.047)
                        {
                            spikeincandidate_events_ewing47->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.048)
                        {
                            spikeincandidate_events_ewing48->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.049)
                        {
                            spikeincandidate_events_ewing49->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.05)
                        {
                            spikeincandidate_events_ewing5->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.07)
                        {
                            spikeincandidate_events_ewing7->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.09)
                        {
                            spikeincandidate_events_ewing9->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.10)
                        {
                            spikeincandidate_events_ewing10->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.12)
                        {
                            spikeincandidate_events_ewing12->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.14)
                        {
                            spikeincandidate_events_ewing14->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.15)
                        {
                            spikeincandidate_events_ewing15->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing < 0.03)
                    {
                        file13 << "--------cand case 4--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }

                    if (Ewing < 0.01)
                    {
                        file14 << "--------cand case 4--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }
                }

                //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                    file2 << "-------------Cand4_leftonly-------------" << "\n"
                    << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                    << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                    << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                    << run << " " << lumis << " " << event << "\n"
                    << endl;

                    candidate_eta_wing_leftonly4->Fill(Ewing);
                    candidate_eta_wing4->Fill(Ewing);
                    candidate_eta_wing->Fill(Ewing);

                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                    {
                        candidate_eta_wing_prompt->Fill(Ewing);
                    }

                    if (Ewing > 0.03)
                    {
                        candidate_events_ewing3_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                        candidate_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing3_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.01)
                    {
                        candidate_events_ewing1_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                        candidate_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing1_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }
                    //---------------------------more ewing cuts-------------------------------

                    if (Ewing > 0.02)
                    {
                        candidate_events_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing2_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.025)
                    {
                        candidate_events_ewing25->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing25_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.031)
                    {
                        candidate_events_ewing31->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing31_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.032)
                    {
                        candidate_events_ewing32->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing32_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.033)
                    {
                        candidate_events_ewing33->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing33_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.034)
                    {
                        candidate_events_ewing34->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing34_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.035)
                    {
                        candidate_events_ewing35->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing35_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.036)
                    {
                        candidate_events_ewing36->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing36_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.037)
                    {
                        candidate_events_ewing37->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing37_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.038)
                    {
                        candidate_events_ewing38->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing38_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.039)
                    {
                        candidate_events_ewing39->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing39_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }
                    

                    if (Ewing > 0.040)
                    {
                        candidate_events_ewing4->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing4_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.041)
                    {
                        candidate_events_ewing41->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing41_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.042)
                    {
                        candidate_events_ewing42->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing42_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.043)
                    {
                        candidate_events_ewing43->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing43_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.044)
                    {
                        candidate_events_ewing44->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing44_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.045)
                    {
                        candidate_events_ewing45->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing45_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.046)
                    {
                        candidate_events_ewing46->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing46_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.047)
                    {
                        candidate_events_ewing47->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing47_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.048)
                    {
                        candidate_events_ewing48->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing48_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.049)
                    {
                        candidate_events_ewing49->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing49_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }


                    if (Ewing > 0.05)
                    {
                        candidate_events_ewing5->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing5_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.07)
                    {
                        candidate_events_ewing7->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing7_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.09)
                    {
                        candidate_events_ewing9->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing9_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.10)
                    {
                        candidate_events_ewing10->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing10_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.12)
                    {
                        candidate_events_ewing12->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing12_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.14)
                    {
                        candidate_events_ewing14->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing14_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.15)
                    {
                        candidate_events_ewing15->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing15_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }



                    if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                    {
                        spikeincandidate_eta_wing_leftonly4->Fill(Ewing);
                        spikeincandidate_eta_wing4->Fill(Ewing);
                        spikeincandidate_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spikeincandidate_events_ewing3_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spikeincandidate_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.01)
                        {
                            spikeincandidate_events_ewing1_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spikeincandidate_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        //---------------------------more ewing cuts-------------------------------

                        if (Ewing > 0.02)
                        {
                            spikeincandidate_events_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.025)
                        {
                            spikeincandidate_events_ewing25->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.031)
                        {
                            spikeincandidate_events_ewing31->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.032)
                        {
                            spikeincandidate_events_ewing32->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.033)
                        {
                            spikeincandidate_events_ewing33->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.034)
                        {
                            spikeincandidate_events_ewing34->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.035)
                        {
                            spikeincandidate_events_ewing35->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.036)
                        {
                            spikeincandidate_events_ewing36->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.037)
                        {
                            spikeincandidate_events_ewing37->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.038)
                        {
                            spikeincandidate_events_ewing38->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.039)
                        {
                            spikeincandidate_events_ewing39->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                        

                        if (Ewing > 0.040)
                        {
                            spikeincandidate_events_ewing4->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.041)
                        {
                            spikeincandidate_events_ewing41->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.042)
                        {
                            spikeincandidate_events_ewing42->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.043)
                        {
                            spikeincandidate_events_ewing43->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.044)
                        {
                            spikeincandidate_events_ewing44->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.045)
                        {
                            spikeincandidate_events_ewing45->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.046)
                        {
                            spikeincandidate_events_ewing46->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.047)
                        {
                            spikeincandidate_events_ewing47->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.048)
                        {
                            spikeincandidate_events_ewing48->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.049)
                        {
                            spikeincandidate_events_ewing49->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.05)
                        {
                            spikeincandidate_events_ewing5->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.07)
                        {
                            spikeincandidate_events_ewing7->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.09)
                        {
                            spikeincandidate_events_ewing9->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.10)
                        {
                            spikeincandidate_events_ewing10->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.12)
                        {
                            spikeincandidate_events_ewing12->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.14)
                        {
                            spikeincandidate_events_ewing14->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.15)
                        {
                            spikeincandidate_events_ewing15->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing < 0.03)
                    {
                        file13 << "--------cand case 4--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }

                    if (Ewing < 0.01)
                    {
                        file14 << "--------cand case 4--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }
                }

                //if there are both
                if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                {
                    rightcellsIDX = cellsEB_right[0];
                    leftcellsIDX = cellsEB_left[0];

                    if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[rightcellsIDX];
                    }
                    else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[leftcellsIDX];
                    }

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                    file2 << "-------------Cand4_both-------------" << "\n"
                    << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                    << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                    << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                    << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                    << run << " " << lumis << " " << event << "\n"
                    << endl;

                    candidate_eta_wing_both4->Fill(Ewing);
                    candidate_eta_wing4->Fill(Ewing);
                    candidate_eta_wing->Fill(Ewing);

                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                    {
                        candidate_eta_wing_prompt->Fill(Ewing);
                    }

                    if (Ewing > 0.03)
                    {
                        candidate_events_ewing3_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                        candidate_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing3_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.01)
                    {
                        candidate_events_ewing1_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                        candidate_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing1_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    //---------------------------more ewing cuts-------------------------------

                    if (Ewing > 0.02)
                    {
                        candidate_events_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing2_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.025)
                    {
                        candidate_events_ewing25->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing25_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.031)
                    {
                        candidate_events_ewing31->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing31_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.032)
                    {
                        candidate_events_ewing32->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing32_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.033)
                    {
                        candidate_events_ewing33->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing33_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.034)
                    {
                        candidate_events_ewing34->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing34_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.035)
                    {
                        candidate_events_ewing35->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing35_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.036)
                    {
                        candidate_events_ewing36->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing36_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.037)
                    {
                        candidate_events_ewing37->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing37_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.038)
                    {
                        candidate_events_ewing38->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing38_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.039)
                    {
                        candidate_events_ewing39->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing39_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }
                    

                    if (Ewing > 0.040)
                    {
                        candidate_events_ewing4->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing4_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.041)
                    {
                        candidate_events_ewing41->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing41_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.042)
                    {
                        candidate_events_ewing42->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing42_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.043)
                    {
                        candidate_events_ewing43->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing43_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.044)
                    {
                        candidate_events_ewing44->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing44_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.045)
                    {
                        candidate_events_ewing45->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing45_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.046)
                    {
                        candidate_events_ewing46->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing46_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.047)
                    {
                        candidate_events_ewing47->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing47_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.048)
                    {
                        candidate_events_ewing48->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing48_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.049)
                    {
                        candidate_events_ewing49->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing49_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }


                    if (Ewing > 0.05)
                    {
                        candidate_events_ewing5->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing5_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.07)
                    {
                        candidate_events_ewing7->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing7_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.09)
                    {
                        candidate_events_ewing9->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing9_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.10)
                    {
                        candidate_events_ewing10->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing10_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.12)
                    {
                        candidate_events_ewing12->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing12_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.14)
                    {
                        candidate_events_ewing14->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing14_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.15)
                    {
                        candidate_events_ewing15->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing15_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }



                    if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                    {
                        spikeincandidate_eta_wing_both4->Fill(Ewing);
                        spikeincandidate_eta_wing4->Fill(Ewing);
                        spikeincandidate_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spikeincandidate_events_ewing3_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spikeincandidate_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.01)
                        {
                            spikeincandidate_events_ewing1_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spikeincandidate_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        //---------------------------more ewing cuts-------------------------------

                        if (Ewing > 0.02)
                        {
                            spikeincandidate_events_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.025)
                        {
                            spikeincandidate_events_ewing25->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.031)
                        {
                            spikeincandidate_events_ewing31->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.032)
                        {
                            spikeincandidate_events_ewing32->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.033)
                        {
                            spikeincandidate_events_ewing33->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.034)
                        {
                            spikeincandidate_events_ewing34->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.035)
                        {
                            spikeincandidate_events_ewing35->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.036)
                        {
                            spikeincandidate_events_ewing36->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.037)
                        {
                            spikeincandidate_events_ewing37->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.038)
                        {
                            spikeincandidate_events_ewing38->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.039)
                        {
                            spikeincandidate_events_ewing39->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                        

                        if (Ewing > 0.040)
                        {
                            spikeincandidate_events_ewing4->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.041)
                        {
                            spikeincandidate_events_ewing41->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.042)
                        {
                            spikeincandidate_events_ewing42->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.043)
                        {
                            spikeincandidate_events_ewing43->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.044)
                        {
                            spikeincandidate_events_ewing44->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.045)
                        {
                            spikeincandidate_events_ewing45->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.046)
                        {
                            spikeincandidate_events_ewing46->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.047)
                        {
                            spikeincandidate_events_ewing47->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.048)
                        {
                            spikeincandidate_events_ewing48->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.049)
                        {
                            spikeincandidate_events_ewing49->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.05)
                        {
                            spikeincandidate_events_ewing5->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.07)
                        {
                            spikeincandidate_events_ewing7->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.09)
                        {
                            spikeincandidate_events_ewing9->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.10)
                        {
                            spikeincandidate_events_ewing10->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.12)
                        {
                            spikeincandidate_events_ewing12->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.14)
                        {
                            spikeincandidate_events_ewing14->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.15)
                        {
                            spikeincandidate_events_ewing15->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing < 0.03)
                    {
                        file13 << "--------cand case 4--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }

                    if (Ewing < 0.01)
                    {
                        file14 << "--------cand case 4--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }
                }
            }

            //BeamHalo template
			if (NewMET > 210 &&
				(*ophohasPixelSeed)[mootIDXSEL] == 0 &&
				(*ophoMIPTotEnergy)[mootIDXSEL] > 4.9 &&
				(*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 &&
				(*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001 &&
				metFilters == 8)
			{
				BeamHalo_template_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
				BeamHalo_template->Fill((*ophoSeedTime)[mootIDXSEL]);

				sieie_BeamHalo->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL]);


                if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.0176)
                {
                    BeamHalo_template_sieieb176->Fill((*ophoSeedTime)[mootIDXSEL]);
                }
                
                if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.0206)
                {
                    BeamHalo_template_sieieb206->Fill((*ophoSeedTime)[mootIDXSEL]);
                }

                if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] < 0.0176)
                {
                    BeamHalo_template_sieies176->Fill((*ophoSeedTime)[mootIDXSEL]);
                }

                if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] < 0.0206)
                {
                    BeamHalo_template_sieies206->Fill((*ophoSeedTime)[mootIDXSEL]);
                }
			}

            //Spikes template
            if (NewMET > 210 &&
                (*ophohasPixelSeed)[mootIDXSEL] == 0 &&
                (*ophoMIPTotEnergy)[mootIDXSEL] < 4.9 &&
                ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] < 0.001 ||
                (*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] < 0.001) &&
                metFilters == 0)
            {
                //original showershap cut
                if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] < 0.01015)
                {  
                    spike_events_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                    //there are some seed without energy and timing info
                    //spike_events->Fill((*ophoSeedTime)[mootIDXSEL]);
                    
                    spike_events_etaphi_case4->Fill((*ophoEta)[mootIDXSEL], (*ophoPhi)[mootIDXSEL]);
                    spike_events_etaphi->Fill((*ophoEta)[mootIDXSEL], (*ophoPhi)[mootIDXSEL]);

                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                    {
                        //spike_events_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                    }

                    file11 << "-------------spike4_extra-------------" << "\n"
                    << "phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n"
                    << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << "\n" 
                    << run << " " << lumis << " " << event << "\n"
                    << endl;

                    //cand1 searching adjencent cells
                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {                    
                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] + 1)
                        {
                            cellsEB_right.push_back(icell); 

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*ophoSeedIEta)[mootIDXSEL]) == 86 || (*ophoSeedIEta)[mootIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike4_err_right---------------------" << "\n"
                                << "phoSeediphi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedieta = " << (*ophoSeedIEta)[mootIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }

                        //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                        {
                            cellsEB_right.push_back(icell); 

                            file9 << "-------------spike4_ietamid_right-------------" << "\n"
                            << "phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }

                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] - 1)
                        {                    
                            cellsEB_left.push_back(icell);

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*ophoSeedIEta)[mootIDXSEL]) == 86 || (*ophoSeedIEta)[mootIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike4_err_left---------------------" << "\n"
                                << "phoSeediphi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedieta = " << (*ophoSeedIEta)[mootIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }            

                        //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                        {
                            cellsEB_left.push_back(icell); 

                            file9 << "-------------spike4_ietamid_left-------------" << "\n"
                            << "phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }  
                    
                    }

                    //Now I have to take if there's only one and both
                    //If there's only right
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                    {
                        rightcellsIDX_only = cellsEB_right[0];
                        maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                        file4 << "-------------spike4_rightonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_rightonly4->Fill(Ewing);
                        spike_eta_wing4->Fill(Ewing);
                        spike_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);

                            if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                            {
                                spike_events_ewing3_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                            }
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                            if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                            {
                                spike_events_ewing1_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                            }
                        }
                    }

                    //if there is only left
                    if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                    {
                        leftcellsIDX_only = cellsEB_left[0];
                        maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                        file4 << "-------------spike4_leftonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_leftonly4->Fill(Ewing);
                        spike_eta_wing4->Fill(Ewing);
                        spike_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);

                            if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                            {
                                spike_events_ewing3_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                            }
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                            if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                            {
                                spike_events_ewing1_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                            }
                        }
                    }

                    //if there are both
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                    {
                        rightcellsIDX = cellsEB_right[0];
                        leftcellsIDX = cellsEB_left[0];

                        if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[rightcellsIDX];
                        }
                        else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[leftcellsIDX];
                        }

                        Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                        file4 << "-------------spike4_both-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_both4->Fill(Ewing);
                        spike_eta_wing4->Fill(Ewing);
                        spike_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);

                            if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                            {
                                spike_events_ewing3_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                            }
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                            if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                            {
                                spike_events_ewing1_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                            }
                        }
                    }
                }

                //time cut no showershape cut
                if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                {
                    spike_events_time_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                    spike_events_time->Fill((*ophoSeedTime)[mootIDXSEL]);
                    
                    spike_events_time_etaphi_case4->Fill((*ophoEta)[mootIDXSEL], (*ophoPhi)[mootIDXSEL]);
                    spike_events_time_etaphi->Fill((*ophoEta)[mootIDXSEL], (*ophoPhi)[mootIDXSEL]);

                    file6 << run << " " << lumis << " " << event << endl;

                    //cand1 searching adjencent cells
                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {                    
                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] + 1)
                        {
                            cellsEB_right.push_back(icell); 

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*ophoSeedIEta)[mootIDXSEL]) == 86 || (*ophoSeedIEta)[mootIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike4_time_err_right---------------------" << "\n"
                                << "phoSeediphi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedieta = " << (*ophoSeedIEta)[mootIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }

                        //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                        {
                            cellsEB_right.push_back(icell); 

                            file10 << "-------------spike4_time_ietamid_right-------------" << "\n"
                            << "phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }

                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] - 1)
                        {                    
                            cellsEB_left.push_back(icell);

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*ophoSeedIEta)[mootIDXSEL]) == 86 || (*ophoSeedIEta)[mootIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike4_time_err_left---------------------" << "\n"
                                << "phoSeediphi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedieta = " << (*ophoSeedIEta)[mootIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }            

                        //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                        {
                            cellsEB_left.push_back(icell); 

                            file10 << "-------------spike4_time_ietamid_left-------------" << "\n"
                            << "phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }
                    }

                    //Now I have to take if there's only one and both
                    //If there's only right
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                    {
                        rightcellsIDX_only = cellsEB_right[0];
                        maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                        file5 << "-------------spike4_time_rightonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_time_rightonly4->Fill(Ewing);
                        spike_eta_wing_time4->Fill(Ewing);
                        spike_eta_wing_time->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_time_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing3_time->Fill((*ophoSeedTime)[mootIDXSEL]);

                            file15 << "--------spike case 4--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_time_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing1_time->Fill((*ophoSeedTime)[mootIDXSEL]);

                            file16 << "--------spike case 4--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }
                    }

                    //if there is only left
                    if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                    {
                        leftcellsIDX_only = cellsEB_left[0];
                        maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                        file5 << "-------------spike4_time_leftonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_time_leftonly4->Fill(Ewing);
                        spike_eta_wing_time4->Fill(Ewing);
                        spike_eta_wing_time->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_time_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing3_time->Fill((*ophoSeedTime)[mootIDXSEL]);

                            file15 << "--------spike case 4--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_time_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing1_time->Fill((*ophoSeedTime)[mootIDXSEL]);

                            file16 << "--------spike case 4--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }
                    }

                    //if there are both
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                    {
                        rightcellsIDX = cellsEB_right[0];
                        leftcellsIDX = cellsEB_left[0];

                        if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[rightcellsIDX];
                        }
                        else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[leftcellsIDX];
                        }

                        Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                        file5 << "-------------spike4_time_both-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_time_both4->Fill(Ewing);
                        spike_eta_wing_time4->Fill(Ewing);
                        spike_eta_wing_time->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_time_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing3_time->Fill((*ophoSeedTime)[mootIDXSEL]);

                            file15 << "--------spike case 4--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_time_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing1_time->Fill((*ophoSeedTime)[mootIDXSEL]);

                            file16 << "--------spike case 4--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }
                    }
                }
            }

            //Case 5, it + oot but take oot
            //Candidate Events
            if (pfMET > 210 &&
                (*phohasPixelSeed)[mitIDXSEL] == 0 &&
                (*phoMIPTotEnergy)[mitIDXSEL] < 4.9 &&
                (*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.001 &&
                (*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] > 0.001 &&
                (*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] < 0.01015 &&
                metFilters == 0 &&
                NewMET > 210 &&
                (*ophohasPixelSeed)[mootIDXSEL] == 0 &&
                (*ophoMIPTotEnergy)[mootIDXSEL] < 4.9 &&
                (*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 &&
                (*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001 &&
                (*ophoSigmaIEtaIEtaFull5x5)[mitIDXSEL] < 0.01015)
            {
                candidate_events_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
				candidate_events->Fill((*ophoSeedTime)[mootIDXSEL]);

                if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                {
                    candidate_events_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                }

                if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                {
                    spikeincandidate_events_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
				    spikeincandidate_events->Fill((*ophoSeedTime)[mootIDXSEL]);
                }

                //cand1 searching adjencent cells
                for (int icell = 0; icell < nAllCellsEB; icell++)
                {
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] + 1)
                    {
                        cellsEB_right.push_back(icell); 

                        //if there's any err
                        if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*ophoSeedIEta)[mootIDXSEL]) == 86 || (*ophoSeedIEta)[mootIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                        {
                            file7 << "------------------cand5_err_right---------------------" << "\n"
                            << "phoSeediphi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedieta = " << (*ophoSeedIEta)[mootIDXSEL] << "\n"
                            << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                        }
                    }

                    //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                    {
                        cellsEB_right.push_back(icell); 

                        file8 << "-------------Cand5_ietamid_right-------------" << "\n"
                        << "phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;
                    }

                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] - 1)
                    {                    
                        cellsEB_left.push_back(icell);

                        //if there's any err
                        if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*ophoSeedIEta)[mootIDXSEL]) == 86 || (*ophoSeedIEta)[mootIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                        {
                            file7 << "------------------cand5_err_left---------------------" << "\n"
                            << "phoSeediphi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedieta = " << (*ophoSeedIEta)[mootIDXSEL] << "\n"
                            << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                        }
                    }            

                    //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                    if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                    {
                        cellsEB_left.push_back(icell); 

                        file8 << "-------------Cand5_ietamid_left-------------" << "\n"
                        << "phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;
                    }  
                
                }

                //Now I have to take if there's only one and both
                //If there's only right
                if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                {
                    rightcellsIDX_only = cellsEB_right[0];
                    maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                    file2 << "-------------Cand5_rightonly-------------" << "\n"
                    << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                    << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                    << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                    << run << " " << lumis << " " << event << "\n"
                    << endl;

                    candidate_eta_wing_rightonly5->Fill(Ewing);
                    candidate_eta_wing5->Fill(Ewing);
                    candidate_eta_wing->Fill(Ewing);

                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                    {
                        candidate_eta_wing_prompt->Fill(Ewing);
                    }

                    if (Ewing > 0.03)
                    {
                        candidate_events_ewing3_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                        candidate_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing3_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.01)
                    {
                        candidate_events_ewing1_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                        candidate_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing1_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    //---------------------------more ewing cuts-------------------------------

                    if (Ewing > 0.02)
                    {
                        candidate_events_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing2_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.025)
                    {
                        candidate_events_ewing25->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing25_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.031)
                    {
                        candidate_events_ewing31->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing31_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.032)
                    {
                        candidate_events_ewing32->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing32_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.033)
                    {
                        candidate_events_ewing33->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing33_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.034)
                    {
                        candidate_events_ewing34->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing34_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.035)
                    {
                        candidate_events_ewing35->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing35_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.036)
                    {
                        candidate_events_ewing36->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing36_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.037)
                    {
                        candidate_events_ewing37->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing37_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.038)
                    {
                        candidate_events_ewing38->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing38_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.039)
                    {
                        candidate_events_ewing39->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing39_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }
                    

                    if (Ewing > 0.040)
                    {
                        candidate_events_ewing4->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing4_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.041)
                    {
                        candidate_events_ewing41->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing41_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.042)
                    {
                        candidate_events_ewing42->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing42_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.043)
                    {
                        candidate_events_ewing43->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing43_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.044)
                    {
                        candidate_events_ewing44->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing44_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.045)
                    {
                        candidate_events_ewing45->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing45_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.046)
                    {
                        candidate_events_ewing46->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing46_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.047)
                    {
                        candidate_events_ewing47->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing47_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.048)
                    {
                        candidate_events_ewing48->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing48_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.049)
                    {
                        candidate_events_ewing49->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing49_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }


                    if (Ewing > 0.05)
                    {
                        candidate_events_ewing5->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing5_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.07)
                    {
                        candidate_events_ewing7->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing7_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.09)
                    {
                        candidate_events_ewing9->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing9_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.10)
                    {
                        candidate_events_ewing10->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing10_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.12)
                    {
                        candidate_events_ewing12->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing12_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.14)
                    {
                        candidate_events_ewing14->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing14_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.15)
                    {
                        candidate_events_ewing15->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing15_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }


                    
                    if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                    {
                        spikeincandidate_eta_wing_rightonly5->Fill(Ewing);
                        spikeincandidate_eta_wing5->Fill(Ewing);
                        spikeincandidate_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spikeincandidate_events_ewing3_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spikeincandidate_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.01)
                        {
                            spikeincandidate_events_ewing1_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spikeincandidate_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        //---------------------------more ewing cuts-------------------------------

                        if (Ewing > 0.02)
                        {
                            spikeincandidate_events_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.025)
                        {
                            spikeincandidate_events_ewing25->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.031)
                        {
                            spikeincandidate_events_ewing31->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.032)
                        {
                            spikeincandidate_events_ewing32->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.033)
                        {
                            spikeincandidate_events_ewing33->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.034)
                        {
                            spikeincandidate_events_ewing34->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.035)
                        {
                            spikeincandidate_events_ewing35->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.036)
                        {
                            spikeincandidate_events_ewing36->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.037)
                        {
                            spikeincandidate_events_ewing37->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.038)
                        {
                            spikeincandidate_events_ewing38->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.039)
                        {
                            spikeincandidate_events_ewing39->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                        

                        if (Ewing > 0.040)
                        {
                            spikeincandidate_events_ewing4->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.041)
                        {
                            spikeincandidate_events_ewing41->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.042)
                        {
                            spikeincandidate_events_ewing42->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.043)
                        {
                            spikeincandidate_events_ewing43->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.044)
                        {
                            spikeincandidate_events_ewing44->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.045)
                        {
                            spikeincandidate_events_ewing45->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.046)
                        {
                            spikeincandidate_events_ewing46->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.047)
                        {
                            spikeincandidate_events_ewing47->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.048)
                        {
                            spikeincandidate_events_ewing48->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.049)
                        {
                            spikeincandidate_events_ewing49->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }



                        if (Ewing > 0.05)
                        {
                            spikeincandidate_events_ewing5->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.07)
                        {
                            spikeincandidate_events_ewing7->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.09)
                        {
                            spikeincandidate_events_ewing9->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.10)
                        {
                            spikeincandidate_events_ewing10->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.12)
                        {
                            spikeincandidate_events_ewing12->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.14)
                        {
                            spikeincandidate_events_ewing14->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.15)
                        {
                            spikeincandidate_events_ewing15->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing < 0.03)
                    {
                        file13 << "--------cand case 5--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }

                    if (Ewing < 0.01)
                    {
                        file14 << "--------cand case 5--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }
                }

                //if there is only left
                if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                {
                    leftcellsIDX_only = cellsEB_left[0];
                    maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                    file2 << "-------------Cand5_leftonly-------------" << "\n"
                    << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                    << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                    << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                    << run << " " << lumis << " " << event << "\n"
                    << endl;

                    candidate_eta_wing_leftonly5->Fill(Ewing);
                    candidate_eta_wing5->Fill(Ewing);
                    candidate_eta_wing->Fill(Ewing);

                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                    {
                        candidate_eta_wing_prompt->Fill(Ewing);
                    }

                    if (Ewing > 0.03)
                    {
                        candidate_events_ewing3_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                        candidate_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing3_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.01)
                    {
                        candidate_events_ewing1_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                        candidate_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing1_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    //---------------------------more ewing cuts-------------------------------

                    if (Ewing > 0.02)
                    {
                        candidate_events_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing2_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.025)
                    {
                        candidate_events_ewing25->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing25_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.031)
                    {
                        candidate_events_ewing31->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing31_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.032)
                    {
                        candidate_events_ewing32->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing32_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.033)
                    {
                        candidate_events_ewing33->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing33_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.034)
                    {
                        candidate_events_ewing34->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing34_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.035)
                    {
                        candidate_events_ewing35->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing35_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.036)
                    {
                        candidate_events_ewing36->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing36_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.037)
                    {
                        candidate_events_ewing37->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing37_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.038)
                    {
                        candidate_events_ewing38->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing38_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.039)
                    {
                        candidate_events_ewing39->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing39_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }
                    

                    if (Ewing > 0.040)
                    {
                        candidate_events_ewing4->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing4_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.041)
                    {
                        candidate_events_ewing41->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing41_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.042)
                    {
                        candidate_events_ewing42->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing42_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.043)
                    {
                        candidate_events_ewing43->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing43_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.044)
                    {
                        candidate_events_ewing44->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing44_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.045)
                    {
                        candidate_events_ewing45->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing45_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.046)
                    {
                        candidate_events_ewing46->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing46_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.047)
                    {
                        candidate_events_ewing47->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing47_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.048)
                    {
                        candidate_events_ewing48->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing48_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.049)
                    {
                        candidate_events_ewing49->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing49_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }


                    if (Ewing > 0.05)
                    {
                        candidate_events_ewing5->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing5_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.07)
                    {
                        candidate_events_ewing7->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing7_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.09)
                    {
                        candidate_events_ewing9->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing9_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.10)
                    {
                        candidate_events_ewing10->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing10_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.12)
                    {
                        candidate_events_ewing12->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing12_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.14)
                    {
                        candidate_events_ewing14->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing14_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.15)
                    {
                        candidate_events_ewing15->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing15_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }



                    if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                    {
                        spikeincandidate_eta_wing_leftonly5->Fill(Ewing);
                        spikeincandidate_eta_wing5->Fill(Ewing);
                        spikeincandidate_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spikeincandidate_events_ewing3_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spikeincandidate_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.01)
                        {
                            spikeincandidate_events_ewing1_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spikeincandidate_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        //---------------------------more ewing cuts-------------------------------

                        if (Ewing > 0.02)
                        {
                            spikeincandidate_events_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.025)
                        {
                            spikeincandidate_events_ewing25->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.031)
                        {
                            spikeincandidate_events_ewing31->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.032)
                        {
                            spikeincandidate_events_ewing32->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.033)
                        {
                            spikeincandidate_events_ewing33->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.034)
                        {
                            spikeincandidate_events_ewing34->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.035)
                        {
                            spikeincandidate_events_ewing35->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.036)
                        {
                            spikeincandidate_events_ewing36->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.037)
                        {
                            spikeincandidate_events_ewing37->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.038)
                        {
                            spikeincandidate_events_ewing38->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.039)
                        {
                            spikeincandidate_events_ewing39->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                        

                        if (Ewing > 0.040)
                        {
                            spikeincandidate_events_ewing4->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.041)
                        {
                            spikeincandidate_events_ewing41->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.042)
                        {
                            spikeincandidate_events_ewing42->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.043)
                        {
                            spikeincandidate_events_ewing43->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.044)
                        {
                            spikeincandidate_events_ewing44->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.045)
                        {
                            spikeincandidate_events_ewing45->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.046)
                        {
                            spikeincandidate_events_ewing46->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.047)
                        {
                            spikeincandidate_events_ewing47->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.048)
                        {
                            spikeincandidate_events_ewing48->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.049)
                        {
                            spikeincandidate_events_ewing49->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.05)
                        {
                            spikeincandidate_events_ewing5->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.07)
                        {
                            spikeincandidate_events_ewing7->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.09)
                        {
                            spikeincandidate_events_ewing9->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.10)
                        {
                            spikeincandidate_events_ewing10->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.12)
                        {
                            spikeincandidate_events_ewing12->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.14)
                        {
                            spikeincandidate_events_ewing14->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.15)
                        {
                            spikeincandidate_events_ewing15->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing < 0.03)
                    {
                        file13 << "--------cand case 5--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }

                    if (Ewing < 0.01)
                    {
                        file14 << "--------cand case 5--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }
                }

                //if there are both
                if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                {
                    rightcellsIDX = cellsEB_right[0];
                    leftcellsIDX = cellsEB_left[0];

                    if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[rightcellsIDX];
                    }
                    else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                    {
                        maxCellsE = AllCellsE_EB[leftcellsIDX];
                    }

                    Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                    file2 << "-------------Cand5_both-------------" << "\n"
                    << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                    << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                    << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                    << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                    << run << " " << lumis << " " << event << "\n"
                    << endl;

                    candidate_eta_wing_both5->Fill(Ewing);
                    candidate_eta_wing5->Fill(Ewing);
                    candidate_eta_wing->Fill(Ewing);

                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                    {
                        candidate_eta_wing_prompt->Fill(Ewing);
                    }

                    if (Ewing > 0.03)
                    {
                        candidate_events_ewing3_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                        candidate_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing3_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.01)
                    {
                        candidate_events_ewing1_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                        candidate_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing1_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    //---------------------------more ewing cuts-------------------------------

                    if (Ewing > 0.02)
                    {
                        candidate_events_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing2_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.025)
                    {
                        candidate_events_ewing25->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing25_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.031)
                    {
                        candidate_events_ewing31->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing31_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.032)
                    {
                        candidate_events_ewing32->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing32_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.033)
                    {
                        candidate_events_ewing33->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing33_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.034)
                    {
                        candidate_events_ewing34->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing34_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.035)
                    {
                        candidate_events_ewing35->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing35_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.036)
                    {
                        candidate_events_ewing36->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing36_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.037)
                    {
                        candidate_events_ewing37->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing37_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.038)
                    {
                        candidate_events_ewing38->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing38_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.039)
                    {
                        candidate_events_ewing39->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing39_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }
                    

                    if (Ewing > 0.040)
                    {
                        candidate_events_ewing4->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing4_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.041)
                    {
                        candidate_events_ewing41->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing41_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.042)
                    {
                        candidate_events_ewing42->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing42_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.043)
                    {
                        candidate_events_ewing43->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing43_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.044)
                    {
                        candidate_events_ewing44->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing44_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.045)
                    {
                        candidate_events_ewing45->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing45_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.046)
                    {
                        candidate_events_ewing46->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing46_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.047)
                    {
                        candidate_events_ewing47->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing47_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.048)
                    {
                        candidate_events_ewing48->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing48_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.049)
                    {
                        candidate_events_ewing49->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing49_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }


                    if (Ewing > 0.05)
                    {
                        candidate_events_ewing5->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing5_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.07)
                    {
                        candidate_events_ewing7->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing7_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.09)
                    {
                        candidate_events_ewing9->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing9_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.10)
                    {
                        candidate_events_ewing10->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing10_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.12)
                    {
                        candidate_events_ewing12->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing12_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.14)
                    {
                        candidate_events_ewing14->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing14_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing > 0.15)
                    {
                        candidate_events_ewing15->Fill((*ophoSeedTime)[mootIDXSEL]);

                        if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                        {
                            candidate_events_ewing15_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }


                    

                    if ((*ophoSeedTime)[mootIDXSEL] < -12.5)
                    {
                        spikeincandidate_eta_wing_both5->Fill(Ewing);
                        spikeincandidate_eta_wing5->Fill(Ewing);
                        spikeincandidate_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spikeincandidate_events_ewing3_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spikeincandidate_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.01)
                        {
                            spikeincandidate_events_ewing1_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spikeincandidate_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        //---------------------------more ewing cuts-------------------------------

                        if (Ewing > 0.02)
                        {
                            spikeincandidate_events_ewing2->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.025)
                        {
                            spikeincandidate_events_ewing25->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.031)
                        {
                            spikeincandidate_events_ewing31->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.032)
                        {
                            spikeincandidate_events_ewing32->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.033)
                        {
                            spikeincandidate_events_ewing33->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.034)
                        {
                            spikeincandidate_events_ewing34->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.035)
                        {
                            spikeincandidate_events_ewing35->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.036)
                        {
                            spikeincandidate_events_ewing36->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.037)
                        {
                            spikeincandidate_events_ewing37->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.038)
                        {
                            spikeincandidate_events_ewing38->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.039)
                        {
                            spikeincandidate_events_ewing39->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                        

                        if (Ewing > 0.040)
                        {
                            spikeincandidate_events_ewing4->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.041)
                        {
                            spikeincandidate_events_ewing41->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.042)
                        {
                            spikeincandidate_events_ewing42->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.043)
                        {
                            spikeincandidate_events_ewing43->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.044)
                        {
                            spikeincandidate_events_ewing44->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.045)
                        {
                            spikeincandidate_events_ewing45->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.046)
                        {
                            spikeincandidate_events_ewing46->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.047)
                        {
                            spikeincandidate_events_ewing47->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.048)
                        {
                            spikeincandidate_events_ewing48->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.049)
                        {
                            spikeincandidate_events_ewing49->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.05)
                        {
                            spikeincandidate_events_ewing5->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.07)
                        {
                            spikeincandidate_events_ewing7->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.09)
                        {
                            spikeincandidate_events_ewing9->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.10)
                        {
                            spikeincandidate_events_ewing10->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.12)
                        {
                            spikeincandidate_events_ewing12->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.14)
                        {
                            spikeincandidate_events_ewing14->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }

                        if (Ewing > 0.15)
                        {
                            spikeincandidate_events_ewing15->Fill((*ophoSeedTime)[mootIDXSEL]);
                        }
                    }

                    if (Ewing < 0.03)
                    {
                        file13 << "--------cand case 5--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }

                    if (Ewing < 0.01)
                    {
                        file14 << "--------cand case 5--------" << "\n"
                        << "Ewing = " << Ewing << "\n"
                        << run << lumis << event << endl;
                    }
                }
            }

            //BeamHalo template
			if (pfMET > 210 &&
				(*phohasPixelSeed)[mitIDXSEL] == 0 &&
				(*phoMIPTotEnergy)[mitIDXSEL] > 4.9 &&
				(*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] > 0.001 &&
				(*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] > 0.001 &&
				metFilters == 8 &&
				NewMET > 210 &&
				(*ophohasPixelSeed)[mootIDXSEL] == 0 &&
				(*ophoMIPTotEnergy)[mootIDXSEL] > 4.9 &&
				(*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.001 &&
				(*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] > 0.001)
			{
				BeamHalo_template_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
				BeamHalo_template->Fill((*ophoSeedTime)[mootIDXSEL]);

				sieie_BeamHalo->Fill((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL]);


                if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.0176)
                {
                    BeamHalo_template_sieieb176->Fill((*ophoSeedTime)[mootIDXSEL]);
                }
                
                if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] > 0.0206)
                {
                    BeamHalo_template_sieieb206->Fill((*ophoSeedTime)[mootIDXSEL]);
                }

                if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] < 0.0176)
                {
                    BeamHalo_template_sieies176->Fill((*ophoSeedTime)[mootIDXSEL]);
                }

                if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] < 0.0206)
                {
                    BeamHalo_template_sieies206->Fill((*ophoSeedTime)[mootIDXSEL]);
                }
			}

            //Spikes template
			if (pfMET > 210 &&
				(*phohasPixelSeed)[mitIDXSEL] == 0 &&
				(*phoMIPTotEnergy)[mitIDXSEL] < 4.9 &&
				((*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] < 0.001 ||
				(*phoSigmaIPhiIPhiFull5x5)[mitIDXSEL] < 0.001) &&
				metFilters == 0 &&
				NewMET > 210 &&
				(*ophohasPixelSeed)[mootIDXSEL] == 0 &&
				(*ophoMIPTotEnergy)[mootIDXSEL] < 4.9 &&
				((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] < 0.001 ||
				(*ophoSigmaIPhiIPhiFull5x5)[mootIDXSEL] < 0.001))
			{
                //original showershap cut
                if ((*ophoSigmaIEtaIEtaFull5x5)[mootIDXSEL] < 0.01015 && (*phoSigmaIEtaIEtaFull5x5)[mitIDXSEL] < 0.01015)
                {
                    spike_events_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                    //there are some seed without energy and timing info
                    //spike_events->Fill((*ophoSeedTime)[mootIDXSEL]);
                    
                    spike_events_etaphi_case5->Fill((*ophoEta)[mootIDXSEL], (*ophoPhi)[mootIDXSEL]);
                    spike_events_etaphi->Fill((*ophoEta)[mootIDXSEL], (*ophoPhi)[mootIDXSEL]);

                    if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                    {
                        //spike_events_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                    }

                    file11 << "-------------spike4_extra-------------" << "\n"
                    << "phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n"
                    << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << "\n" 
                    << run << " " << lumis << " " << event << "\n"
                    << endl;

                    //cand1 searching adjencent cells
                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {                    
                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] + 1)
                        {
                            cellsEB_right.push_back(icell); 

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*ophoSeedIEta)[mootIDXSEL]) == 86 || (*ophoSeedIEta)[mootIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike5_err_right---------------------" << "\n"
                                << "phoSeediphi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedieta = " << (*ophoSeedIEta)[mootIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }

                        //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                        {
                            cellsEB_right.push_back(icell); 

                            file9 << "-------------spike5_ietamid_right-------------" << "\n"
                            << "phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }

                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] - 1)
                        {                    
                            cellsEB_left.push_back(icell);

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*ophoSeedIEta)[mootIDXSEL]) == 86 || (*ophoSeedIEta)[mootIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike5_err_left---------------------" << "\n"
                                << "phoSeediphi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedieta = " << (*ophoSeedIEta)[mootIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }            

                        //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                        {
                            cellsEB_left.push_back(icell); 

                            file9 << "-------------spike5_ietamid_left-------------" << "\n"
                            << "phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }  
                    
                    }

                    //Now I have to take if there's only one and both
                    //If there's only right
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                    {
                        rightcellsIDX_only = cellsEB_right[0];
                        maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                        file4 << "-------------spike5_rightonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_rightonly5->Fill(Ewing);
                        spike_eta_wing5->Fill(Ewing);
                        spike_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);

                            if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                            {
                                spike_events_ewing3_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                            }
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                            if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                            {
                                spike_events_ewing1_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                            }
                        }
                    }

                    //if there is only left
                    if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                    {
                        leftcellsIDX_only = cellsEB_left[0];
                        maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                        file4 << "-------------spike5_leftonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_leftonly5->Fill(Ewing);
                        spike_eta_wing5->Fill(Ewing);
                        spike_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);

                            if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                            {
                                spike_events_ewing3_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                            }
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                            if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                            {
                                spike_events_ewing1_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                            }
                        }
                    }

                    //if there are both
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                    {
                        rightcellsIDX = cellsEB_right[0];
                        leftcellsIDX = cellsEB_left[0];

                        if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[rightcellsIDX];
                        }
                        else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[leftcellsIDX];
                        }

                        Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                        file4 << "-------------spike5_both-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_both5->Fill(Ewing);
                        spike_eta_wing5->Fill(Ewing);
                        spike_eta_wing->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing3->Fill((*ophoSeedTime)[mootIDXSEL]);

                            if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                            {
                                spike_events_ewing3_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                            }
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_case4->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing1->Fill((*ophoSeedTime)[mootIDXSEL]);

                            if (fabs((*ophoSeedTime)[mootIDXSEL]) < 3.0)
                            {
                                spike_events_ewing1_3ns->Fill((*ophoSeedTime)[mootIDXSEL]);
                            }
                        }
                    }
                }

                //time cut no showershape cut
                if ((*ophoSeedTime)[mootIDXSEL] < -12.5 && (*phoSeedTime)[mitIDXSEL] < -12.5)
                {
                    spike_events_time_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                    spike_events_time->Fill((*ophoSeedTime)[mootIDXSEL]);
                    
                    spike_events_time_etaphi_case5->Fill((*ophoEta)[mootIDXSEL], (*ophoPhi)[mootIDXSEL]);
                    spike_events_time_etaphi->Fill((*ophoEta)[mootIDXSEL], (*ophoPhi)[mootIDXSEL]);

                    file6 << run << " " << lumis << " " << event << endl;

                    //cand1 searching adjencent cells
                    for (int icell = 0; icell < nAllCellsEB; icell++)
                    {                    
                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] + 1)
                        {
                            cellsEB_right.push_back(icell); 

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*ophoSeedIEta)[mootIDXSEL]) == 86 || (*ophoSeedIEta)[mootIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike5_time_err_right---------------------" << "\n"
                                << "phoSeediphi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedieta = " << (*ophoSeedIEta)[mootIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }

                        //since there's no 0, right of the -1 is 1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == -1 && AllCellsIEtaEB[icell] == 1)
                        {
                            cellsEB_right.push_back(icell); 

                            file10 << "-------------spike5_time_ietamid_right-------------" << "\n"
                            << "phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                            << "CellIPhiEB_right = " << AllCellsIPhiEB[icell] << " CellIEtaEB_right = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_right = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }

                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && AllCellsIEtaEB[icell] == (*ophoSeedIEta)[mootIDXSEL] - 1)
                        {                    
                            cellsEB_left.push_back(icell);

                            //if there's any err
                            if (fabs(AllCellsIEtaEB[icell]) == 86 || fabs((*ophoSeedIEta)[mootIDXSEL]) == 86 || (*ophoSeedIEta)[mootIDXSEL] == 0 || AllCellsIEtaEB[icell] == 0)
                            {
                                file7 << "------------------spike5_time_err_left---------------------" << "\n"
                                << "phoSeediphi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedieta = " << (*ophoSeedIEta)[mootIDXSEL] << "\n"
                                << "cellsiphi = " << AllCellsIPhiEB[icell] << " cellsieta = " << AllCellsIEtaEB[icell] << endl;
                            }
                        }            

                        //since there's no 0, left of the 1 is -1 => -2 -1 1 2
                        if (AllCellsIPhiEB[icell] == (*ophoSeedIPhi)[mootIDXSEL] && (*ophoSeedIEta)[mootIDXSEL] == 1 && AllCellsIEtaEB[icell] == -1)
                        {
                            cellsEB_left.push_back(icell); 

                            file10 << "-------------spike5_time_ietamid_left-------------" << "\n"
                            << "phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                            << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                            << "CellIPhiEB_left = " << AllCellsIPhiEB[icell] << " CellIEtaEB_left = " << AllCellsIEtaEB[icell] << " AllCellsE_EB_left = " << AllCellsE_EB[icell] << "\n"
                            << run << " " << lumis << " " << event << "\n"
                            << endl;
                        }
                    }

                    //Now I have to take if there's only one and both
                    //If there's only right
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() == 0)
                    {
                        rightcellsIDX_only = cellsEB_right[0];
                        maxCellsE = AllCellsE_EB[rightcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                        file5 << "-------------spike5_time_rightonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX_only] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX_only] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_time_rightonly5->Fill(Ewing);
                        spike_eta_wing_time5->Fill(Ewing);
                        spike_eta_wing_time->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_time_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing3_time->Fill((*ophoSeedTime)[mootIDXSEL]);

                            file15 << "--------spike case 5--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_time_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing1_time->Fill((*ophoSeedTime)[mootIDXSEL]);

                            file16 << "--------spike case 5--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }
                    }

                    //if there is only left
                    if (cellsEB_right.size() == 0 && cellsEB_left.size() != 0)
                    {
                        leftcellsIDX_only = cellsEB_left[0];
                        maxCellsE = AllCellsE_EB[leftcellsIDX_only];

                        Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                        file5 << "-------------spike5_time_leftonly-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX_only] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX_only] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX_only] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_time_leftonly5->Fill(Ewing);
                        spike_eta_wing_time5->Fill(Ewing);
                        spike_eta_wing_time->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_time_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing3_time->Fill((*ophoSeedTime)[mootIDXSEL]);

                            file15 << "--------spike case 5--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_time_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing1_time->Fill((*ophoSeedTime)[mootIDXSEL]);

                            file16 << "--------spike case 5--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }
                    }

                    //if there are both
                    if (cellsEB_right.size() != 0 && cellsEB_left.size() != 0)
                    {
                        rightcellsIDX = cellsEB_right[0];
                        leftcellsIDX = cellsEB_left[0];

                        if (AllCellsE_EB[rightcellsIDX] > AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[rightcellsIDX];
                        }
                        else if (AllCellsE_EB[rightcellsIDX] < AllCellsE_EB[leftcellsIDX])
                        {
                            maxCellsE = AllCellsE_EB[leftcellsIDX];
                        }

                        Float_t Ewing = maxCellsE/((*ophoSeedEnergy)[mootIDXSEL]);

                        file5 << "-------------spike5_time_both-------------" << "\n"
                        << "Ewing = " << Ewing << " maxCellsE = " << maxCellsE << " phoseedTime = " << (*ophoSeedTime)[mootIDXSEL] << "\n"
                        << "phoSeedIPhi = " << (*ophoSeedIPhi)[mootIDXSEL] << " phoSeedIEta = " << (*ophoSeedIEta)[mootIDXSEL] << " phoSeedEnergy = " << (*ophoSeedEnergy)[mootIDXSEL] << "\n" 
                        << "CellIPhiEB_left = " << AllCellsIPhiEB[leftcellsIDX] << " CellIEtaEB_left = " << AllCellsIEtaEB[leftcellsIDX] << " AllCellsE_EB_left = " << AllCellsE_EB[leftcellsIDX] << "\n"
                        << "CellIPhiEB_right = " << AllCellsIPhiEB[rightcellsIDX] << " CellIEtaEB_right = " << AllCellsIEtaEB[rightcellsIDX] << " AllCellsE_EB_right = " << AllCellsE_EB[rightcellsIDX] << "\n"
                        << run << " " << lumis << " " << event << "\n"
                        << endl;

                        spike_eta_wing_time_both5->Fill(Ewing);
                        spike_eta_wing_time5->Fill(Ewing);
                        spike_eta_wing_time->Fill(Ewing);

                        if (Ewing > 0.03)
                        {
                            spike_events_ewing3_time_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing3_time->Fill((*ophoSeedTime)[mootIDXSEL]);

                            file15 << "--------spike case 5--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }

                        if (Ewing > 0.01)
                        {
                            spike_events_ewing1_time_case5->Fill((*ophoSeedTime)[mootIDXSEL]);
                            spike_events_ewing1_time->Fill((*ophoSeedTime)[mootIDXSEL]);

                            file16 << "--------spike case 5--------" << "\n"
                            << "Ewing = " << Ewing << "\n"
                            << run << lumis << event << endl;
                        }
                    }
                }
            }
        }
    }

    phoseedieta->SetLineColor(1);
	phoseedieta->Write();

    phoseediphi->SetLineColor(1);
	phoseediphi->Write();

    ophoseedieta->SetLineColor(1);
	ophoseedieta->Write();

    ophoseediphi->SetLineColor(1);
	ophoseediphi->Write();


    //-------------------candidate--------------------
    candidate_events_case1->SetLineColor(1);
	candidate_events_case1->Write();

    candidate_events_case2->SetLineColor(1);
	candidate_events_case2->Write();

    candidate_events_case3->SetLineColor(1);
	candidate_events_case3->Write();

    candidate_events_case4->SetLineColor(1);
	candidate_events_case4->Write();

    candidate_events_case5->SetLineColor(1);
	candidate_events_case5->Write();

    candidate_events->SetLineColor(1);
	candidate_events->Write();

    candidate_events_3ns->SetLineColor(1);
	candidate_events_3ns->Write();



    candidate_eta_wing_rightonly1->SetLineColor(1);
	candidate_eta_wing_rightonly1->Write();

    candidate_eta_wing_rightonly2->SetLineColor(1);
	candidate_eta_wing_rightonly2->Write();

    candidate_eta_wing_rightonly3->SetLineColor(1);
	candidate_eta_wing_rightonly3->Write();

    candidate_eta_wing_rightonly4->SetLineColor(1);
	candidate_eta_wing_rightonly4->Write();

    candidate_eta_wing_rightonly5->SetLineColor(1);
	candidate_eta_wing_rightonly5->Write();


    candidate_eta_wing_leftonly1->SetLineColor(1);
	candidate_eta_wing_leftonly1->Write();

    candidate_eta_wing_leftonly2->SetLineColor(1);
	candidate_eta_wing_leftonly2->Write();

    candidate_eta_wing_leftonly3->SetLineColor(1);
	candidate_eta_wing_leftonly3->Write();

    candidate_eta_wing_leftonly4->SetLineColor(1);
	candidate_eta_wing_leftonly4->Write();

    candidate_eta_wing_leftonly5->SetLineColor(1);
	candidate_eta_wing_leftonly5->Write();


    candidate_eta_wing_both1->SetLineColor(1);
	candidate_eta_wing_both1->Write();

    candidate_eta_wing_both2->SetLineColor(1);
	candidate_eta_wing_both2->Write();

    candidate_eta_wing_both3->SetLineColor(1);
	candidate_eta_wing_both3->Write();

    candidate_eta_wing_both4->SetLineColor(1);
	candidate_eta_wing_both4->Write();

    candidate_eta_wing_both5->SetLineColor(1);
	candidate_eta_wing_both5->Write();


    candidate_eta_wing1->SetLineColor(1);
	candidate_eta_wing1->Write();

    candidate_eta_wing2->SetLineColor(1);
	candidate_eta_wing2->Write();

    candidate_eta_wing3->SetLineColor(1);
	candidate_eta_wing3->Write();

    candidate_eta_wing4->SetLineColor(1);
	candidate_eta_wing4->Write();

    candidate_eta_wing5->SetLineColor(1);
	candidate_eta_wing5->Write();

    candidate_eta_wing->SetLineColor(1);
	candidate_eta_wing->Write();

    candidate_eta_wing_prompt->SetLineColor(1);
	candidate_eta_wing_prompt->Write();



    //-------------------spikeincandidate--------------------
    spikeincandidate_events_case1->SetLineColor(1);
	spikeincandidate_events_case1->Write();

    spikeincandidate_events_case2->SetLineColor(1);
	spikeincandidate_events_case2->Write();

    spikeincandidate_events_case3->SetLineColor(1);
	spikeincandidate_events_case3->Write();

    spikeincandidate_events_case4->SetLineColor(1);
	spikeincandidate_events_case4->Write();

    spikeincandidate_events_case5->SetLineColor(1);
	spikeincandidate_events_case5->Write();

    spikeincandidate_events->SetLineColor(1);
	spikeincandidate_events->Write();

    spikeincandidate_events_3ns->SetLineColor(1);
	spikeincandidate_events_3ns->Write();



    spikeincandidate_eta_wing_rightonly1->SetLineColor(1);
	spikeincandidate_eta_wing_rightonly1->Write();

    spikeincandidate_eta_wing_rightonly2->SetLineColor(1);
	spikeincandidate_eta_wing_rightonly2->Write();

    spikeincandidate_eta_wing_rightonly3->SetLineColor(1);
	spikeincandidate_eta_wing_rightonly3->Write();

    spikeincandidate_eta_wing_rightonly4->SetLineColor(1);
	spikeincandidate_eta_wing_rightonly4->Write();

    spikeincandidate_eta_wing_rightonly5->SetLineColor(1);
	spikeincandidate_eta_wing_rightonly5->Write();


    spikeincandidate_eta_wing_leftonly1->SetLineColor(1);
	spikeincandidate_eta_wing_leftonly1->Write();

    spikeincandidate_eta_wing_leftonly2->SetLineColor(1);
	spikeincandidate_eta_wing_leftonly2->Write();

    spikeincandidate_eta_wing_leftonly3->SetLineColor(1);
	spikeincandidate_eta_wing_leftonly3->Write();

    spikeincandidate_eta_wing_leftonly4->SetLineColor(1);
	spikeincandidate_eta_wing_leftonly4->Write();

    spikeincandidate_eta_wing_leftonly5->SetLineColor(1);
	spikeincandidate_eta_wing_leftonly5->Write();


    spikeincandidate_eta_wing_both1->SetLineColor(1);
	spikeincandidate_eta_wing_both1->Write();

    spikeincandidate_eta_wing_both2->SetLineColor(1);
	spikeincandidate_eta_wing_both2->Write();

    spikeincandidate_eta_wing_both3->SetLineColor(1);
	spikeincandidate_eta_wing_both3->Write();

    spikeincandidate_eta_wing_both4->SetLineColor(1);
	spikeincandidate_eta_wing_both4->Write();

    spikeincandidate_eta_wing_both5->SetLineColor(1);
	spikeincandidate_eta_wing_both5->Write();


    spikeincandidate_eta_wing1->SetLineColor(1);
	spikeincandidate_eta_wing1->Write();

    spikeincandidate_eta_wing2->SetLineColor(1);
	spikeincandidate_eta_wing2->Write();

    spikeincandidate_eta_wing3->SetLineColor(1);
	spikeincandidate_eta_wing3->Write();

    spikeincandidate_eta_wing4->SetLineColor(1);
	spikeincandidate_eta_wing4->Write();

    spikeincandidate_eta_wing5->SetLineColor(1);
	spikeincandidate_eta_wing5->Write();

    spikeincandidate_eta_wing->SetLineColor(1);
	spikeincandidate_eta_wing->Write();




    //-------------------promptZ--------------------
    promptZ_events_case1->SetLineColor(1);
	promptZ_events_case1->Write();

    promptZ_events_case2->SetLineColor(1);
	promptZ_events_case2->Write();

    promptZ_events_case3->SetLineColor(1);
	promptZ_events_case3->Write();

    promptZ_events_case4->SetLineColor(1);
	promptZ_events_case4->Write();

    promptZ_events_case5->SetLineColor(1);
	promptZ_events_case5->Write();

    promptZ_events->SetLineColor(1);
	promptZ_events->Write();

    promptZ_events_3ns->SetLineColor(1);
	promptZ_events_3ns->Write();


    InvMass_Z->SetLineColor(1);
	InvMass_Z->Write();


    promptZ_eta_wing_rightonly1->SetLineColor(1);
	promptZ_eta_wing_rightonly1->Write();

    promptZ_eta_wing_rightonly2->SetLineColor(1);
	promptZ_eta_wing_rightonly2->Write();

    promptZ_eta_wing_rightonly3->SetLineColor(1);
	promptZ_eta_wing_rightonly3->Write();

    promptZ_eta_wing_rightonly4->SetLineColor(1);
	promptZ_eta_wing_rightonly4->Write();

    promptZ_eta_wing_rightonly5->SetLineColor(1);
	promptZ_eta_wing_rightonly5->Write();


    promptZ_eta_wing_leftonly1->SetLineColor(1);
	promptZ_eta_wing_leftonly1->Write();

    promptZ_eta_wing_leftonly2->SetLineColor(1);
	promptZ_eta_wing_leftonly2->Write();

    promptZ_eta_wing_leftonly3->SetLineColor(1);
	promptZ_eta_wing_leftonly3->Write();

    promptZ_eta_wing_leftonly4->SetLineColor(1);
	promptZ_eta_wing_leftonly4->Write();

    promptZ_eta_wing_leftonly5->SetLineColor(1);
	promptZ_eta_wing_leftonly5->Write();


    promptZ_eta_wing_both1->SetLineColor(1);
	promptZ_eta_wing_both1->Write();

    promptZ_eta_wing_both2->SetLineColor(1);
	promptZ_eta_wing_both2->Write();

    promptZ_eta_wing_both3->SetLineColor(1);
	promptZ_eta_wing_both3->Write();

    promptZ_eta_wing_both4->SetLineColor(1);
	promptZ_eta_wing_both4->Write();

    promptZ_eta_wing_both5->SetLineColor(1);
	promptZ_eta_wing_both5->Write();


    promptZ_eta_wing1->SetLineColor(1);
	promptZ_eta_wing1->Write();

    promptZ_eta_wing2->SetLineColor(1);
	promptZ_eta_wing2->Write();

    promptZ_eta_wing3->SetLineColor(1);
	promptZ_eta_wing3->Write();

    promptZ_eta_wing4->SetLineColor(1);
	promptZ_eta_wing4->Write();

    promptZ_eta_wing5->SetLineColor(1);
	promptZ_eta_wing5->Write();

    promptZ_eta_wing->SetLineColor(1);
	promptZ_eta_wing->Write();

    promptZ_eta_wing_prompt->SetLineColor(1);
	promptZ_eta_wing_prompt->Write();




    //-------------------spike--------------------
    spike_events_case1->SetLineColor(1);
	spike_events_case1->Write();

    spike_events_case2->SetLineColor(1);
	spike_events_case2->Write();

    spike_events_case3->SetLineColor(1);
	spike_events_case3->Write();

    spike_events_case4->SetLineColor(1);
	spike_events_case4->Write();

    spike_events_case5->SetLineColor(1);
	spike_events_case5->Write();

    spike_events->SetLineColor(1);
	spike_events->Write();

    spike_events_3ns->SetLineColor(1);
	spike_events_3ns->Write();



    spike_eta_wing_rightonly1->SetLineColor(1);
	spike_eta_wing_rightonly1->Write();

    spike_eta_wing_rightonly2->SetLineColor(1);
	spike_eta_wing_rightonly2->Write();

    spike_eta_wing_rightonly3->SetLineColor(1);
	spike_eta_wing_rightonly3->Write();

    spike_eta_wing_rightonly4->SetLineColor(1);
	spike_eta_wing_rightonly4->Write();

    spike_eta_wing_rightonly5->SetLineColor(1);
	spike_eta_wing_rightonly5->Write();


    spike_eta_wing_leftonly1->SetLineColor(1);
	spike_eta_wing_leftonly1->Write();

    spike_eta_wing_leftonly2->SetLineColor(1);
	spike_eta_wing_leftonly2->Write();

    spike_eta_wing_leftonly3->SetLineColor(1);
	spike_eta_wing_leftonly3->Write();

    spike_eta_wing_leftonly4->SetLineColor(1);
	spike_eta_wing_leftonly4->Write();

    spike_eta_wing_leftonly5->SetLineColor(1);
	spike_eta_wing_leftonly5->Write();


    spike_eta_wing_both1->SetLineColor(1);
	spike_eta_wing_both1->Write();

    spike_eta_wing_both2->SetLineColor(1);
	spike_eta_wing_both2->Write();

    spike_eta_wing_both3->SetLineColor(1);
	spike_eta_wing_both3->Write();

    spike_eta_wing_both4->SetLineColor(1);
	spike_eta_wing_both4->Write();

    spike_eta_wing_both5->SetLineColor(1);
	spike_eta_wing_both5->Write();


    spike_eta_wing1->SetLineColor(1);
	spike_eta_wing1->Write();

    spike_eta_wing2->SetLineColor(1);
	spike_eta_wing2->Write();

    spike_eta_wing3->SetLineColor(1);
	spike_eta_wing3->Write();

    spike_eta_wing4->SetLineColor(1);
	spike_eta_wing4->Write();

    spike_eta_wing5->SetLineColor(1);
	spike_eta_wing5->Write();

    spike_eta_wing->SetLineColor(1);
	spike_eta_wing->Write();



    spike_events_etaphi_case1->SetOption("COLZ");
	spike_events_etaphi_case1->Write();

    spike_events_etaphi_case2->SetOption("COLZ");
	spike_events_etaphi_case2->Write();

    spike_events_etaphi_case3->SetOption("COLZ");
	spike_events_etaphi_case3->Write();

    spike_events_etaphi_case4->SetOption("COLZ");
	spike_events_etaphi_case4->Write();

    spike_events_etaphi_case5->SetOption("COLZ");
	spike_events_etaphi_case5->Write();

    spike_events_etaphi->SetOption("COLZ");
	spike_events_etaphi->Write();



    //-------------------spike_time--------------------
    spike_events_time_case1->SetLineColor(1);
	spike_events_time_case1->Write();

    spike_events_time_case2->SetLineColor(1);
	spike_events_time_case2->Write();

    spike_events_time_case3->SetLineColor(1);
	spike_events_time_case3->Write();

    spike_events_time_case4->SetLineColor(1);
	spike_events_time_case4->Write();

    spike_events_time_case5->SetLineColor(1);
	spike_events_time_case5->Write();

    spike_events_time->SetLineColor(1);
	spike_events_time->Write();

    spike_events_time_3ns->SetLineColor(1);
	spike_events_time_3ns->Write();




    spike_eta_wing_time_rightonly1->SetLineColor(1);
	spike_eta_wing_time_rightonly1->Write();

    spike_eta_wing_time_rightonly2->SetLineColor(1);
	spike_eta_wing_time_rightonly2->Write();

    spike_eta_wing_time_rightonly3->SetLineColor(1);
	spike_eta_wing_time_rightonly3->Write();

    spike_eta_wing_time_rightonly4->SetLineColor(1);
	spike_eta_wing_time_rightonly4->Write();

    spike_eta_wing_time_rightonly5->SetLineColor(1);
	spike_eta_wing_time_rightonly5->Write();


    spike_eta_wing_time_leftonly1->SetLineColor(1);
	spike_eta_wing_time_leftonly1->Write();

    spike_eta_wing_time_leftonly2->SetLineColor(1);
	spike_eta_wing_time_leftonly2->Write();

    spike_eta_wing_time_leftonly3->SetLineColor(1);
	spike_eta_wing_time_leftonly3->Write();

    spike_eta_wing_time_leftonly4->SetLineColor(1);
	spike_eta_wing_time_leftonly4->Write();

    spike_eta_wing_time_leftonly5->SetLineColor(1);
	spike_eta_wing_time_leftonly5->Write();


    spike_eta_wing_time_both1->SetLineColor(1);
	spike_eta_wing_time_both1->Write();

    spike_eta_wing_time_both2->SetLineColor(1);
	spike_eta_wing_time_both2->Write();

    spike_eta_wing_time_both3->SetLineColor(1);
	spike_eta_wing_time_both3->Write();

    spike_eta_wing_time_both4->SetLineColor(1);
	spike_eta_wing_time_both4->Write();

    spike_eta_wing_time_both5->SetLineColor(1);
	spike_eta_wing_time_both5->Write();


    spike_eta_wing_time1->SetLineColor(1);
	spike_eta_wing_time1->Write();

    spike_eta_wing_time2->SetLineColor(1);
	spike_eta_wing_time2->Write();

    spike_eta_wing_time3->SetLineColor(1);
	spike_eta_wing_time3->Write();

    spike_eta_wing_time4->SetLineColor(1);
	spike_eta_wing_time4->Write();

    spike_eta_wing_time5->SetLineColor(1);
	spike_eta_wing_time5->Write();

    spike_eta_wing_time->SetLineColor(1);
	spike_eta_wing_time->Write();



    spike_events_time_etaphi_case1->SetOption("COLZ");
	spike_events_time_etaphi_case1->Write();

    spike_events_time_etaphi_case2->SetOption("COLZ");
	spike_events_time_etaphi_case2->Write();

    spike_events_time_etaphi_case3->SetOption("COLZ");
	spike_events_time_etaphi_case3->Write();

    spike_events_time_etaphi_case4->SetOption("COLZ");
	spike_events_time_etaphi_case4->Write();

    spike_events_time_etaphi_case5->SetOption("COLZ");
	spike_events_time_etaphi_case5->Write();

    spike_events_time_etaphi->SetOption("COLZ");
	spike_events_time_etaphi->Write();



    //-------------------candidate ewing cut > 0.03--------------------
    candidate_events_ewing3_case1->SetLineColor(1);
	candidate_events_ewing3_case1->Write();

    candidate_events_ewing3_case2->SetLineColor(1);
	candidate_events_ewing3_case2->Write();

    candidate_events_ewing3_case3->SetLineColor(1);
	candidate_events_ewing3_case3->Write();

    candidate_events_ewing3_case4->SetLineColor(1);
	candidate_events_ewing3_case4->Write();

    candidate_events_ewing3_case5->SetLineColor(1);
	candidate_events_ewing3_case5->Write();

    candidate_events_ewing3->SetLineColor(1);
	candidate_events_ewing3->Write();

    candidate_events_ewing3_3ns->SetLineColor(1);
	candidate_events_ewing3_3ns->Write();

    //-------------------candidate ewing cut > 0.01--------------------
    candidate_events_ewing1_case1->SetLineColor(1);
	candidate_events_ewing1_case1->Write();

    candidate_events_ewing1_case2->SetLineColor(1);
	candidate_events_ewing1_case2->Write();

    candidate_events_ewing1_case3->SetLineColor(1);
	candidate_events_ewing1_case3->Write();

    candidate_events_ewing1_case4->SetLineColor(1);
	candidate_events_ewing1_case4->Write();

    candidate_events_ewing1_case5->SetLineColor(1);
	candidate_events_ewing1_case5->Write();

    candidate_events_ewing1->SetLineColor(1);
	candidate_events_ewing1->Write();

    candidate_events_ewing1_3ns->SetLineColor(1);
	candidate_events_ewing1_3ns->Write();


    //-------------------candidate ewing cut > 0.02--------------------
    candidate_events_ewing2->SetLineColor(1);
	candidate_events_ewing2->Write();

    candidate_events_ewing2_3ns->SetLineColor(1);
	candidate_events_ewing2_3ns->Write();

    //-------------------candidate ewing cut > 0.025--------------------
    candidate_events_ewing25->SetLineColor(1);
	candidate_events_ewing25->Write();

    candidate_events_ewing25_3ns->SetLineColor(1);
	candidate_events_ewing25_3ns->Write();

    //-------------------candidate ewing cut > 0.031--------------------
    candidate_events_ewing31->SetLineColor(1);
	candidate_events_ewing31->Write();

    candidate_events_ewing31_3ns->SetLineColor(1);
	candidate_events_ewing31_3ns->Write();

    //-------------------candidate ewing cut > 0.032--------------------
    candidate_events_ewing32->SetLineColor(1);
	candidate_events_ewing32->Write();

    candidate_events_ewing32_3ns->SetLineColor(1);
	candidate_events_ewing32_3ns->Write();

    //-------------------candidate ewing cut > 0.033--------------------
    candidate_events_ewing33->SetLineColor(1);
	candidate_events_ewing33->Write();

    candidate_events_ewing33_3ns->SetLineColor(1);
	candidate_events_ewing33_3ns->Write();

    //-------------------candidate ewing cut > 0.034--------------------
    candidate_events_ewing34->SetLineColor(1);
	candidate_events_ewing34->Write();

    candidate_events_ewing34_3ns->SetLineColor(1);
	candidate_events_ewing34_3ns->Write();

    //-------------------candidate ewing cut > 0.035--------------------
    candidate_events_ewing35->SetLineColor(1);
	candidate_events_ewing35->Write();

    candidate_events_ewing35_3ns->SetLineColor(1);
	candidate_events_ewing35_3ns->Write();

    //-------------------candidate ewing cut > 0.036--------------------
    candidate_events_ewing36->SetLineColor(1);
	candidate_events_ewing36->Write();

    candidate_events_ewing36_3ns->SetLineColor(1);
	candidate_events_ewing36_3ns->Write();

    //-------------------candidate ewing cut > 0.037--------------------
    candidate_events_ewing37->SetLineColor(1);
	candidate_events_ewing37->Write();

    candidate_events_ewing37_3ns->SetLineColor(1);
	candidate_events_ewing37_3ns->Write();

    //-------------------candidate ewing cut > 0.038--------------------
    candidate_events_ewing38->SetLineColor(1);
	candidate_events_ewing38->Write();

    candidate_events_ewing38_3ns->SetLineColor(1);
	candidate_events_ewing38_3ns->Write();

    //-------------------candidate ewing cut > 0.039--------------------
    candidate_events_ewing39->SetLineColor(1);
	candidate_events_ewing39->Write();

    candidate_events_ewing39_3ns->SetLineColor(1);
	candidate_events_ewing39_3ns->Write();

    //-------------------candidate ewing cut > 0.040--------------------
    candidate_events_ewing4->SetLineColor(1);
	candidate_events_ewing4->Write();

    candidate_events_ewing4_3ns->SetLineColor(1);
	candidate_events_ewing4_3ns->Write();

    //-------------------candidate ewing cut > 0.041--------------------
    candidate_events_ewing41->SetLineColor(1);
	candidate_events_ewing41->Write();

    candidate_events_ewing41_3ns->SetLineColor(1);
	candidate_events_ewing41_3ns->Write();

    //-------------------candidate ewing cut > 0.042--------------------
    candidate_events_ewing42->SetLineColor(1);
	candidate_events_ewing42->Write();

    candidate_events_ewing42_3ns->SetLineColor(1);
	candidate_events_ewing42_3ns->Write();

    //-------------------candidate ewing cut > 0.043--------------------
    candidate_events_ewing43->SetLineColor(1);
	candidate_events_ewing43->Write();

    candidate_events_ewing43_3ns->SetLineColor(1);
	candidate_events_ewing43_3ns->Write();

    //-------------------candidate ewing cut > 0.044--------------------
    candidate_events_ewing44->SetLineColor(1);
	candidate_events_ewing44->Write();

    candidate_events_ewing44_3ns->SetLineColor(1);
	candidate_events_ewing44_3ns->Write();

    //-------------------candidate ewing cut > 0.045--------------------
    candidate_events_ewing45->SetLineColor(1);
	candidate_events_ewing45->Write();

    candidate_events_ewing45_3ns->SetLineColor(1);
	candidate_events_ewing45_3ns->Write();

    //-------------------candidate ewing cut > 0.046--------------------
    candidate_events_ewing46->SetLineColor(1);
	candidate_events_ewing46->Write();

    candidate_events_ewing46_3ns->SetLineColor(1);
	candidate_events_ewing46_3ns->Write();

    //-------------------candidate ewing cut > 0.047--------------------
    candidate_events_ewing47->SetLineColor(1);
	candidate_events_ewing47->Write();

    candidate_events_ewing47_3ns->SetLineColor(1);
	candidate_events_ewing47_3ns->Write();

    //-------------------candidate ewing cut > 0.048--------------------
    candidate_events_ewing48->SetLineColor(1);
	candidate_events_ewing48->Write();

    candidate_events_ewing48_3ns->SetLineColor(1);
	candidate_events_ewing48_3ns->Write();

    //-------------------candidate ewing cut > 0.049--------------------
    candidate_events_ewing49->SetLineColor(1);
	candidate_events_ewing49->Write();

    candidate_events_ewing49_3ns->SetLineColor(1);
	candidate_events_ewing49_3ns->Write();



    //-------------------candidate ewing cut > 0.05--------------------
    candidate_events_ewing5->SetLineColor(1);
	candidate_events_ewing5->Write();

    candidate_events_ewing5_3ns->SetLineColor(1);
	candidate_events_ewing5_3ns->Write();

    //-------------------candidate ewing cut > 0.07--------------------
    candidate_events_ewing7->SetLineColor(1);
	candidate_events_ewing7->Write();

    candidate_events_ewing7_3ns->SetLineColor(1);
	candidate_events_ewing7_3ns->Write();

    //-------------------candidate ewing cut > 0.09--------------------
    candidate_events_ewing9->SetLineColor(1);
	candidate_events_ewing9->Write();

    candidate_events_ewing9_3ns->SetLineColor(1);
	candidate_events_ewing9_3ns->Write();

    //-------------------candidate ewing cut > 0.1--------------------
    candidate_events_ewing10->SetLineColor(1);
	candidate_events_ewing10->Write();

    candidate_events_ewing10_3ns->SetLineColor(1);
	candidate_events_ewing10_3ns->Write();

    //-------------------candidate ewing cut > 0.12--------------------
    candidate_events_ewing12->SetLineColor(1);
	candidate_events_ewing12->Write();

    candidate_events_ewing12_3ns->SetLineColor(1);
	candidate_events_ewing12_3ns->Write();

    //-------------------candidate ewing cut > 0.14--------------------
    candidate_events_ewing14->SetLineColor(1);
	candidate_events_ewing14->Write();

    candidate_events_ewing14_3ns->SetLineColor(1);
	candidate_events_ewing14_3ns->Write();

    //-------------------candidate ewing cut > 0.15--------------------
    candidate_events_ewing15->SetLineColor(1);
	candidate_events_ewing15->Write();

    candidate_events_ewing15_3ns->SetLineColor(1);
	candidate_events_ewing15_3ns->Write();



    //-------------------spikeincandidate ewing cut > 0.03--------------------
    spikeincandidate_events_ewing3_case1->SetLineColor(1);
	spikeincandidate_events_ewing3_case1->Write();

    spikeincandidate_events_ewing3_case2->SetLineColor(1);
	spikeincandidate_events_ewing3_case2->Write();

    spikeincandidate_events_ewing3_case3->SetLineColor(1);
	spikeincandidate_events_ewing3_case3->Write();

    spikeincandidate_events_ewing3_case4->SetLineColor(1);
	spikeincandidate_events_ewing3_case4->Write();

    spikeincandidate_events_ewing3_case5->SetLineColor(1);
	spikeincandidate_events_ewing3_case5->Write();

    spikeincandidate_events_ewing3->SetLineColor(1);
	spikeincandidate_events_ewing3->Write();

    //-------------------spikeincandidate ewing cut > 0.01--------------------
    spikeincandidate_events_ewing1_case1->SetLineColor(1);
	spikeincandidate_events_ewing1_case1->Write();

    spikeincandidate_events_ewing1_case2->SetLineColor(1);
	spikeincandidate_events_ewing1_case2->Write();

    spikeincandidate_events_ewing1_case3->SetLineColor(1);
	spikeincandidate_events_ewing1_case3->Write();

    spikeincandidate_events_ewing1_case4->SetLineColor(1);
	spikeincandidate_events_ewing1_case4->Write();

    spikeincandidate_events_ewing1_case5->SetLineColor(1);
	spikeincandidate_events_ewing1_case5->Write();

    spikeincandidate_events_ewing1->SetLineColor(1);
	spikeincandidate_events_ewing1->Write();


    
    //-------------------spikeincandidate ewing cut > 0.02--------------------
    spikeincandidate_events_ewing2->SetLineColor(1);
	spikeincandidate_events_ewing2->Write();

    //-------------------spikeincandidate ewing cut > 0.025--------------------
    spikeincandidate_events_ewing25->SetLineColor(1);
	spikeincandidate_events_ewing25->Write();

    //-------------------spikeincandidate ewing cut > 0.031--------------------
    spikeincandidate_events_ewing31->SetLineColor(1);
	spikeincandidate_events_ewing31->Write();

    //-------------------spikeincandidate ewing cut > 0.032--------------------
    spikeincandidate_events_ewing32->SetLineColor(1);
	spikeincandidate_events_ewing32->Write();

    //-------------------spikeincandidate ewing cut > 0.033--------------------
    spikeincandidate_events_ewing33->SetLineColor(1);
	spikeincandidate_events_ewing33->Write();

    //-------------------spikeincandidate ewing cut > 0.034--------------------
    spikeincandidate_events_ewing34->SetLineColor(1);
	spikeincandidate_events_ewing34->Write();

    //-------------------spikeincandidate ewing cut > 0.035--------------------
    spikeincandidate_events_ewing35->SetLineColor(1);
	spikeincandidate_events_ewing35->Write();

    //-------------------spikeincandidate ewing cut > 0.036--------------------
    spikeincandidate_events_ewing36->SetLineColor(1);
	spikeincandidate_events_ewing36->Write();

    //-------------------spikeincandidate ewing cut > 0.037--------------------
    spikeincandidate_events_ewing37->SetLineColor(1);
	spikeincandidate_events_ewing37->Write();

    //-------------------spikeincandidate ewing cut > 0.038--------------------
    spikeincandidate_events_ewing38->SetLineColor(1);
	spikeincandidate_events_ewing38->Write();

    //-------------------spikeincandidate ewing cut > 0.039--------------------
    spikeincandidate_events_ewing39->SetLineColor(1);
	spikeincandidate_events_ewing39->Write();

    //-------------------spikeincandidate ewing cut > 0.040--------------------
    spikeincandidate_events_ewing4->SetLineColor(1);
	spikeincandidate_events_ewing4->Write();

    //-------------------spikeincandidate ewing cut > 0.041--------------------
    spikeincandidate_events_ewing41->SetLineColor(1);
	spikeincandidate_events_ewing41->Write();

    //-------------------spikeincandidate ewing cut > 0.042--------------------
    spikeincandidate_events_ewing42->SetLineColor(1);
	spikeincandidate_events_ewing42->Write();

    //-------------------spikeincandidate ewing cut > 0.043--------------------
    spikeincandidate_events_ewing43->SetLineColor(1);
	spikeincandidate_events_ewing43->Write();

    //-------------------spikeincandidate ewing cut > 0.044--------------------
    spikeincandidate_events_ewing44->SetLineColor(1);
	spikeincandidate_events_ewing44->Write();

    //-------------------spikeincandidate ewing cut > 0.045--------------------
    spikeincandidate_events_ewing45->SetLineColor(1);
	spikeincandidate_events_ewing45->Write();

    //-------------------spikeincandidate ewing cut > 0.046--------------------
    spikeincandidate_events_ewing46->SetLineColor(1);
	spikeincandidate_events_ewing46->Write();

    //-------------------spikeincandidate ewing cut > 0.047--------------------
    spikeincandidate_events_ewing47->SetLineColor(1);
	spikeincandidate_events_ewing47->Write();

    //-------------------spikeincandidate ewing cut > 0.048--------------------
    spikeincandidate_events_ewing48->SetLineColor(1);
	spikeincandidate_events_ewing48->Write();

    //-------------------spikeincandidate ewing cut > 0.049--------------------
    spikeincandidate_events_ewing49->SetLineColor(1);
	spikeincandidate_events_ewing49->Write();

    //-------------------spikeincandidate ewing cut > 0.05--------------------
    spikeincandidate_events_ewing5->SetLineColor(1);
	spikeincandidate_events_ewing5->Write();

    //-------------------spikeincandidate ewing cut > 0.07--------------------
    spikeincandidate_events_ewing7->SetLineColor(1);
	spikeincandidate_events_ewing7->Write();

    //-------------------spikeincandidate ewing cut > 0.09--------------------
    spikeincandidate_events_ewing9->SetLineColor(1);
	spikeincandidate_events_ewing9->Write();

    //-------------------spikeincandidate ewing cut > 0.10--------------------
    spikeincandidate_events_ewing10->SetLineColor(1);
	spikeincandidate_events_ewing10->Write();

    //-------------------spikeincandidate ewing cut > 0.12--------------------
    spikeincandidate_events_ewing12->SetLineColor(1);
	spikeincandidate_events_ewing12->Write();

    //-------------------spikeincandidate ewing cut > 0.14--------------------
    spikeincandidate_events_ewing14->SetLineColor(1);
	spikeincandidate_events_ewing14->Write();

    //-------------------spikeincandidate ewing cut > 0.15--------------------
    spikeincandidate_events_ewing15->SetLineColor(1);
	spikeincandidate_events_ewing15->Write();


    //-------------------promptZ ewing cut > 0.03--------------------
    promptZ_events_ewing3_case1->SetLineColor(1);
	promptZ_events_ewing3_case1->Write();

    promptZ_events_ewing3_case2->SetLineColor(1);
	promptZ_events_ewing3_case2->Write();

    promptZ_events_ewing3_case3->SetLineColor(1);
	promptZ_events_ewing3_case3->Write();

    promptZ_events_ewing3_case4->SetLineColor(1);
	promptZ_events_ewing3_case4->Write();

    promptZ_events_ewing3_case5->SetLineColor(1);
	promptZ_events_ewing3_case5->Write();

    promptZ_events_ewing3->SetLineColor(1);
	promptZ_events_ewing3->Write();

    promptZ_events_ewing3_3ns->SetLineColor(1);
	promptZ_events_ewing3_3ns->Write();

    //-------------------promptZ ewing cut > 0.01--------------------
    promptZ_events_ewing1_case1->SetLineColor(1);
	promptZ_events_ewing1_case1->Write();

    promptZ_events_ewing1_case2->SetLineColor(1);
	promptZ_events_ewing1_case2->Write();

    promptZ_events_ewing1_case3->SetLineColor(1);
	promptZ_events_ewing1_case3->Write();

    promptZ_events_ewing1_case4->SetLineColor(1);
	promptZ_events_ewing1_case4->Write();

    promptZ_events_ewing1_case5->SetLineColor(1);
	promptZ_events_ewing1_case5->Write();

    promptZ_events_ewing1->SetLineColor(1);
	promptZ_events_ewing1->Write();

    promptZ_events_ewing1_3ns->SetLineColor(1);
	promptZ_events_ewing1_3ns->Write();


    
    //-------------------promptZ ewing cut > 0.02--------------------
    promptZ_events_ewing2->SetLineColor(1);
	promptZ_events_ewing2->Write();

    promptZ_events_ewing2_3ns->SetLineColor(1);
	promptZ_events_ewing2_3ns->Write();

    //-------------------promptZ ewing cut > 0.025--------------------
    promptZ_events_ewing25->SetLineColor(1);
	promptZ_events_ewing25->Write();

    promptZ_events_ewing25_3ns->SetLineColor(1);
	promptZ_events_ewing25_3ns->Write();

    //-------------------promptZ ewing cut > 0.031--------------------
    promptZ_events_ewing31->SetLineColor(1);
	promptZ_events_ewing31->Write();

    promptZ_events_ewing31_3ns->SetLineColor(1);
	promptZ_events_ewing31_3ns->Write();

    //-------------------promptZ ewing cut > 0.032--------------------
    promptZ_events_ewing32->SetLineColor(1);
	promptZ_events_ewing32->Write();

    promptZ_events_ewing32_3ns->SetLineColor(1);
	promptZ_events_ewing32_3ns->Write();

    //-------------------promptZ ewing cut > 0.033--------------------
    promptZ_events_ewing33->SetLineColor(1);
	promptZ_events_ewing33->Write();

    promptZ_events_ewing33_3ns->SetLineColor(1);
	promptZ_events_ewing33_3ns->Write();

    //-------------------promptZ ewing cut > 0.034--------------------
    promptZ_events_ewing34->SetLineColor(1);
	promptZ_events_ewing34->Write();

    promptZ_events_ewing34_3ns->SetLineColor(1);
	promptZ_events_ewing34_3ns->Write();

    //-------------------promptZ ewing cut > 0.035--------------------
    promptZ_events_ewing35->SetLineColor(1);
	promptZ_events_ewing35->Write();

    promptZ_events_ewing35_3ns->SetLineColor(1);
	promptZ_events_ewing35_3ns->Write();

    //-------------------promptZ ewing cut > 0.036--------------------
    promptZ_events_ewing36->SetLineColor(1);
	promptZ_events_ewing36->Write();

    promptZ_events_ewing36_3ns->SetLineColor(1);
	promptZ_events_ewing36_3ns->Write();

    //-------------------promptZ ewing cut > 0.037--------------------
    promptZ_events_ewing37->SetLineColor(1);
	promptZ_events_ewing37->Write();

    promptZ_events_ewing37_3ns->SetLineColor(1);
	promptZ_events_ewing37_3ns->Write();

    //-------------------promptZ ewing cut > 0.038--------------------
    promptZ_events_ewing38->SetLineColor(1);
	promptZ_events_ewing38->Write();

    promptZ_events_ewing38_3ns->SetLineColor(1);
	promptZ_events_ewing38_3ns->Write();

    //-------------------promptZ ewing cut > 0.039--------------------
    promptZ_events_ewing39->SetLineColor(1);
	promptZ_events_ewing39->Write();

    promptZ_events_ewing39_3ns->SetLineColor(1);
	promptZ_events_ewing39_3ns->Write();

    //-------------------promptZ ewing cut > 0.040--------------------
    promptZ_events_ewing4->SetLineColor(1);
	promptZ_events_ewing4->Write();

    promptZ_events_ewing4_3ns->SetLineColor(1);
	promptZ_events_ewing4_3ns->Write();

    //-------------------promptZ ewing cut > 0.041--------------------
    promptZ_events_ewing41->SetLineColor(1);
	promptZ_events_ewing41->Write();

    promptZ_events_ewing41_3ns->SetLineColor(1);
	promptZ_events_ewing41_3ns->Write();

    //-------------------promptZ ewing cut > 0.042--------------------
    promptZ_events_ewing42->SetLineColor(1);
	promptZ_events_ewing42->Write();

    promptZ_events_ewing42_3ns->SetLineColor(1);
	promptZ_events_ewing42_3ns->Write();

    //-------------------promptZ ewing cut > 0.043--------------------
    promptZ_events_ewing43->SetLineColor(1);
	promptZ_events_ewing43->Write();

    promptZ_events_ewing43_3ns->SetLineColor(1);
	promptZ_events_ewing43_3ns->Write();

    //-------------------promptZ ewing cut > 0.044--------------------
    promptZ_events_ewing44->SetLineColor(1);
	promptZ_events_ewing44->Write();

    promptZ_events_ewing44_3ns->SetLineColor(1);
	promptZ_events_ewing44_3ns->Write();

    //-------------------promptZ ewing cut > 0.045--------------------
    promptZ_events_ewing45->SetLineColor(1);
	promptZ_events_ewing45->Write();

    promptZ_events_ewing45_3ns->SetLineColor(1);
	promptZ_events_ewing45_3ns->Write();

    //-------------------promptZ ewing cut > 0.046--------------------
    promptZ_events_ewing46->SetLineColor(1);
	promptZ_events_ewing46->Write();

    promptZ_events_ewing46_3ns->SetLineColor(1);
	promptZ_events_ewing46_3ns->Write();

    //-------------------promptZ ewing cut > 0.047--------------------
    promptZ_events_ewing47->SetLineColor(1);
	promptZ_events_ewing47->Write();

    promptZ_events_ewing47_3ns->SetLineColor(1);
	promptZ_events_ewing47_3ns->Write();

    //-------------------promptZ ewing cut > 0.048--------------------
    promptZ_events_ewing48->SetLineColor(1);
	promptZ_events_ewing48->Write();

    promptZ_events_ewing48_3ns->SetLineColor(1);
	promptZ_events_ewing48_3ns->Write();

    //-------------------promptZ ewing cut > 0.049--------------------
    promptZ_events_ewing49->SetLineColor(1);
	promptZ_events_ewing49->Write();

    promptZ_events_ewing49_3ns->SetLineColor(1);
	promptZ_events_ewing49_3ns->Write();


    //-------------------promptZ ewing cut > 0.05--------------------
    promptZ_events_ewing5->SetLineColor(1);
	promptZ_events_ewing5->Write();

    promptZ_events_ewing5_3ns->SetLineColor(1);
	promptZ_events_ewing5_3ns->Write();

    //-------------------promptZ ewing cut > 0.07--------------------
    promptZ_events_ewing7->SetLineColor(1);
	promptZ_events_ewing7->Write();

    promptZ_events_ewing7_3ns->SetLineColor(1);
	promptZ_events_ewing7_3ns->Write();

    //-------------------promptZ ewing cut > 0.09--------------------
    promptZ_events_ewing9->SetLineColor(1);
	promptZ_events_ewing9->Write();

    promptZ_events_ewing9_3ns->SetLineColor(1);
	promptZ_events_ewing9_3ns->Write();

    //-------------------promptZ ewing cut > 0.1--------------------
    promptZ_events_ewing10->SetLineColor(1);
	promptZ_events_ewing10->Write();

    promptZ_events_ewing10_3ns->SetLineColor(1);
	promptZ_events_ewing10_3ns->Write();

    //-------------------promptZ ewing cut > 0.12--------------------
    promptZ_events_ewing12->SetLineColor(1);
	promptZ_events_ewing12->Write();

    promptZ_events_ewing12_3ns->SetLineColor(1);
	promptZ_events_ewing12_3ns->Write();

    //-------------------promptZ ewing cut > 0.14--------------------
    promptZ_events_ewing14->SetLineColor(1);
	promptZ_events_ewing14->Write();

    promptZ_events_ewing14_3ns->SetLineColor(1);
	promptZ_events_ewing14_3ns->Write();

    //-------------------promptZ ewing cut > 0.15--------------------
    promptZ_events_ewing15->SetLineColor(1);
	promptZ_events_ewing15->Write();

    promptZ_events_ewing15_3ns->SetLineColor(1);
	promptZ_events_ewing15_3ns->Write();





    //-------------------spike ewing cut > 0.03--------------------
    spike_events_ewing3_case1->SetLineColor(1);
	spike_events_ewing3_case1->Write();

    spike_events_ewing3_case2->SetLineColor(1);
	spike_events_ewing3_case2->Write();

    spike_events_ewing3_case3->SetLineColor(1);
	spike_events_ewing3_case3->Write();

    spike_events_ewing3_case4->SetLineColor(1);
	spike_events_ewing3_case4->Write();

    spike_events_ewing3_case5->SetLineColor(1);
	spike_events_ewing3_case5->Write();

    spike_events_ewing3->SetLineColor(1);
	spike_events_ewing3->Write();

    spike_events_ewing3_3ns->SetLineColor(1);
	spike_events_ewing3_3ns->Write();

    //-------------------spike ewing cut > 0.01--------------------
    spike_events_ewing1_case1->SetLineColor(1);
	spike_events_ewing1_case1->Write();

    spike_events_ewing1_case2->SetLineColor(1);
	spike_events_ewing1_case2->Write();

    spike_events_ewing1_case3->SetLineColor(1);
	spike_events_ewing1_case3->Write();

    spike_events_ewing1_case4->SetLineColor(1);
	spike_events_ewing1_case4->Write();

    spike_events_ewing1_case5->SetLineColor(1);
	spike_events_ewing1_case5->Write();

    spike_events_ewing1->SetLineColor(1);
	spike_events_ewing1->Write();

    spike_events_ewing1_3ns->SetLineColor(1);
	spike_events_ewing1_3ns->Write();


    //-------------------spike ewing cut > 0.03--------------------
    spike_events_ewing3_time_case1->SetLineColor(1);
	spike_events_ewing3_time_case1->Write();

    spike_events_ewing3_time_case2->SetLineColor(1);
	spike_events_ewing3_time_case2->Write();

    spike_events_ewing3_time_case3->SetLineColor(1);
	spike_events_ewing3_time_case3->Write();

    spike_events_ewing3_time_case4->SetLineColor(1);
	spike_events_ewing3_time_case4->Write();

    spike_events_ewing3_time_case5->SetLineColor(1);
	spike_events_ewing3_time_case5->Write();

    spike_events_ewing3_time->SetLineColor(1);
	spike_events_ewing3_time->Write();

    spike_events_ewing3_time_3ns->SetLineColor(1);
	spike_events_ewing3_time_3ns->Write();

    //-------------------spike ewing _timecut > 0.01--------------------
    spike_events_ewing1_time_case1->SetLineColor(1);
	spike_events_ewing1_time_case1->Write();

    spike_events_ewing1_time_case2->SetLineColor(1);
	spike_events_ewing1_time_case2->Write();

    spike_events_ewing1_time_case3->SetLineColor(1);
	spike_events_ewing1_time_case3->Write();

    spike_events_ewing1_time_case4->SetLineColor(1);
	spike_events_ewing1_time_case4->Write();

    spike_events_ewing1_time_case5->SetLineColor(1);
	spike_events_ewing1_time_case5->Write();

    spike_events_ewing1_time->SetLineColor(1);
	spike_events_ewing1_time->Write();

    spike_events_ewing1_time_3ns->SetLineColor(1);
	spike_events_ewing1_time_3ns->Write();


    //-------------------BeamHalo tempalte--------------------
    BeamHalo_template_case1->SetLineColor(1);
	BeamHalo_template_case1->Write();

	BeamHalo_template_case2->SetLineColor(1);
	BeamHalo_template_case2->Write();

	BeamHalo_template_case3->SetLineColor(1);
	BeamHalo_template_case3->Write();

	BeamHalo_template_case4->SetLineColor(1);
	BeamHalo_template_case4->Write();

	BeamHalo_template_case5->SetLineColor(1);
	BeamHalo_template_case5->Write();	

	BeamHalo_template->SetLineColor(1);
	BeamHalo_template->Write();

    BeamHalo_template_sieieb176->SetLineColor(1);
	BeamHalo_template_sieieb176->Write();

    BeamHalo_template_sieieb206->SetLineColor(1);
	BeamHalo_template_sieieb206->Write();

    BeamHalo_template_sieies176->SetLineColor(1);
	BeamHalo_template_sieies176->Write();

    BeamHalo_template_sieies206->SetLineColor(1);
	BeamHalo_template_sieies206->Write();

    sieie_BeamHalo->SetLineColor(1);
    sieie_BeamHalo->Write();

    


    file1 << "Done filling histograms" << endl;



	file1 << "Analysis completed" << endl;


    cout << "Done filling histograms" << endl;



	cout << "Analysis completed" << endl;
	auto timenow_end = chrono::system_clock::to_time_t(chrono::system_clock::now());

	cout << "start time: " << ctime(&timenow_start) << " (CT)" << endl;

	cout << "End time: " << ctime(&timenow_end) << " (CT)" << endl;

    cout << "Output file: 2017Dataset_eta-wing_spikeincand.root" << endl;

    file1 << "Analysis completed!" << endl;
	file1 << "Completed time: " << ctime(&timenow_end) << endl;

    file1.close();

	file2.close();
	file3.close();
	file4.close();

	file5.close();
	file6.close();
    file7.close();

	file8.close();
	file9.close();
    file10.close();
    file11.close();

    file12.close();
    file13.close();
    file14.close();
    file15.close();
    file16.close();


}
                
            
//04.24.2020
//etawing calculation remake

//04.26.2020
//added plotting with eta-wing cut (0.03 and 0.01)

//04.27.2020
//fixed spike's issue. dont fill 3,4,5 into all events (zero energy and zero time)

//04.28.2020
//the efficiency of spike should come from the timing cut

//05.11.2020
//spot the difference between cmsRun and crab jobs, re-run this.
//adding BH template for timing fit


//05.19.2020
//cmslpc account back
//not applying HLTPho < 2048, just HLTPho >> 10 & 1 == 1
//correct the second oot photon pt (it was 230) = 10.0

//06.03.2020
//adding event picking information for 
//1. candidate who failed the cut
//2. spike who passed the 0.01 cut

//06.04.2020
//adding another tempalte - spikeincandidate
//which this template is a subset of candidate events
//with exactly the same candidate cut and additionally has t < -12.5 ns


//06.08.2020
//adding template eta-wing distribution in prompt window
//adding more template eta-wing cuts: 0.01, 0.03, 0.05, 0.07, 0.09, 0.1, 0.12, 0.14, 0.15
//removed the spikeincandidate subset fill in prompt window


//06.09.2020
//increased the second photon pt > 10 to pt > 15.
//start from 0.05, Ns goes negative... Look into 0.03 - 0.05 region and see what we get.
//also add 0.02



