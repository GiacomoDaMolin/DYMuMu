#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

// include user defined histograms and auxiliary macros
#include "Histodef.cpp"
#include "Auxiliary.cpp"
#include "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/roccor/RoccoR.cc"

using namespace std;

#define MAX_ARRAY_SIZE 128

void DataAnalysis(string inputFile, string ofile)
{

    TFile *fin = TFile::Open(inputFile.c_str());
    TTree *tin = static_cast<TTree *>(fin->Get("Events"));

    // Set all branches to 0
    tin->SetBranchStatus("*", 0);

    Float_t Muon_pt[MAX_ARRAY_SIZE], Jet_pt[MAX_ARRAY_SIZE],Muon_eta[MAX_ARRAY_SIZE], Jet_eta[MAX_ARRAY_SIZE], Muon_phi[MAX_ARRAY_SIZE], Jet_phi[MAX_ARRAY_SIZE],Muon_mass[MAX_ARRAY_SIZE], Jet_mass[MAX_ARRAY_SIZE];
    UInt_t nMuon;
    tin->SetBranchStatus("Muon_pt", 1);
    tin->SetBranchAddress("Muon_pt", &Muon_pt);
    tin->SetBranchStatus("Jet_pt", 1);
    tin->SetBranchAddress("Jet_pt", &Jet_pt);
    tin->SetBranchStatus("nMuon", 1);
    tin->SetBranchAddress("nMuon", &nMuon);
    tin->SetBranchStatus("Muon_eta", 1);
    tin->SetBranchAddress("Muon_eta", &Muon_eta);
    tin->SetBranchStatus("Jet_eta", 1);
    tin->SetBranchAddress("Jet_eta", &Jet_eta);
    tin->SetBranchStatus("Muon_phi", 1);
    tin->SetBranchAddress("Muon_phi", &Muon_phi);
    tin->SetBranchStatus("Jet_phi", 1);
    tin->SetBranchAddress("Jet_phi", &Jet_phi);
    tin->SetBranchStatus("Muon_mass", 1);
    tin->SetBranchAddress("Muon_mass", &Muon_mass);
    tin->SetBranchStatus("Jet_mass", 1);
    tin->SetBranchAddress("Jet_mass", &Jet_mass);

    // collect the trigger information
    Bool_t HLT_IsoMu24;
    tin->SetBranchStatus("HLT_IsoMu24", 1);
    tin->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24);

    Int_t Muon_charge[MAX_ARRAY_SIZE];
    Bool_t Muon_tightId[MAX_ARRAY_SIZE];
    Float_t Muon_pfRelIso04_all[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Muon_tightId", 1);
    tin->SetBranchStatus("Muon_charge", 1);
    tin->SetBranchStatus("Muon_pfRelIso04_all", 1);
    tin->SetBranchAddress("Muon_tightId", &Muon_tightId);
    tin->SetBranchAddress("Muon_charge", &Muon_charge);
    tin->SetBranchAddress("Muon_pfRelIso04_all", &Muon_pfRelIso04_all);

    Float_t Jet_btagDeepFlavB[MAX_ARRAY_SIZE];
    UInt_t nJet;
    Int_t Jet_jetId[MAX_ARRAY_SIZE], Jet_puId[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Jet_btagDeepFlavB", 1);
    tin->SetBranchStatus("nJet", 1);
    tin->SetBranchStatus("Jet_jetId", 1);
    tin->SetBranchStatus("Jet_puId", 1);
    tin->SetBranchAddress("nJet", &nJet);
    tin->SetBranchAddress("Jet_btagDeepFlavB", &Jet_btagDeepFlavB);
    tin->SetBranchAddress("Jet_jetId", &Jet_jetId);
    tin->SetBranchAddress("Jet_puId", &Jet_puId);

    int non_matching_muon = 0;
    int n_dropped = 0;
    int trigger_dropped = 0;
    const auto nEv = tin->GetEntries();
    TLorentzVector *Muon1_p4 = new TLorentzVector();
    TLorentzVector *Muon2_p4 = new TLorentzVector();

       // allow pt, inv mass, and eta to be stored in a Branch
    Float_t invMass, muon1_eta, muon1_pt, muon2_pt, muon2_eta, dphi;
    
    TFile *fout =new TFile(ofile.c_str(),"RECREATE");
    TTree *tout = new TTree("tout","tout");
    tout->Branch("dphi", &dphi);
    tout->Branch("invMass", &invMass);
    tout->Branch("muon2_eta", &muon2_eta);
    tout->Branch("muon2_pt", &muon2_pt);
    tout->Branch("muon1_eta", &muon1_eta);
    tout->Branch("muon1_pt", &muon1_pt);

    RoccoR rc;
    rc.init("/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/roccor/RoccoR2018UL.txt");
    
    #pragma omp parallel for
    for (UInt_t i = 0; i < nEv; i++)
    {
        tin->GetEntry(i);
        if (i % 100000 == 0)
            std::cout << "Processing entry " << i << " of " << nEv << std::endl;
        if (!(HLT_IsoMu24))
        {
            trigger_dropped++;
            continue;
        };
	bool gotmuplus=false,gotmuminus=false;
        for (UInt_t j = 0; j < nMuon; j++){
            if (((Muon_pt[j]>27.||( (gotmuplus||gotmuminus) && Muon_pt[j]>25.)) && abs(Muon_eta[j])<2.4 && Muon_tightId[j] && Muon_pfRelIso04_all[j] < 0.15)){
                if (!gotmuplus && Muon_charge[j]==1){
			double scmDT=rc.kScaleDT(Muon_charge[j],Muon_pt[j],Muon_eta[j],Muon_phi[j]);
			Muon1_p4->SetPtEtaPhiM(Muon_pt[j]*scmDT,Muon_eta[j],Muon_phi[j],Muon_mass[j]);
			gotmuplus=true;
			}
		if (!gotmuminus && Muon_charge[j]==-1){
			double scmDT=rc.kScaleDT(Muon_charge[j],Muon_pt[j],Muon_eta[j],Muon_phi[j]);
			Muon2_p4->SetPtEtaPhiM(Muon_pt[j]*scmDT,Muon_eta[j],Muon_phi[j],Muon_mass[j]);
			gotmuminus=true;
			}
            }
        }
 
	if(!(gotmuplus && gotmuminus)) {n_dropped++; continue;}
	if(Muon1_p4->DeltaR(*Muon2_p4)<0.4) {n_dropped++; continue;}

	int njet=0,nbjet=0;
        for (size_t j = 0; j < nJet; j++){
	  bool passesPUID=(Jet_puId[j]>=4);
          if((abs(Jet_eta[j]) < 2.4) && Jet_pt[j]>25 && (Jet_jetId[j]==2 || Jet_jetId[j]==6) && (Jet_pt[j]>50 || passesPUID))  {
		njet++;
           
	  }
        }
	
       
        dphi=Muon1_p4->DeltaPhi(*Muon2_p4);
       
        muon1_pt = Muon1_p4->Pt();
        muon1_eta = Muon1_p4->Eta();
	muon2_pt = Muon2_p4->Pt();
        muon2_eta = Muon2_p4->Eta();

        h_Muon1_pt->Fill(muon1_pt);
        h_Muon1_eta->Fill(muon1_eta);
	h_Muon2_pt->Fill(muon2_pt);
        h_Muon2_eta->Fill(muon2_eta);
        h_acopla_mumu->Fill(M_PI-dphi);
        h_NJets->Fill(njet);


        invMass = (*(Muon1_p4) + *(Muon2_p4)).M();
        h_Muon_Muon_invariant_mass->Fill(invMass);
	tout->Fill();
    }
    std::cout << "Total number of events: " << nEv << std::endl;
    int NumbEv=nEv;
    std::cout << "trigger dropped = " << trigger_dropped << endl;
    std::cout << "selections dropped = " << n_dropped << endl; 
    std::cout << "Fraction of events discarded by trigger = " << (trigger_dropped * 1. / NumbEv) << endl;
    int Rem_trigger=NumbEv-trigger_dropped; 
    std::cout << "Fraction of events removed by selections = " << (n_dropped * 1. / Rem_trigger) << endl;
    std::cout << "Final number of events "<< Rem_trigger - n_dropped<<endl;
    // Write the histograms to the file
    h_Muon1_eta->Write();
    h_Muon1_pt->Write();
    h_Muon2_eta->Write();
    h_Muon2_pt->Write();

    h_Muon_Muon_invariant_mass->Write();
    h_acopla_mumu->Write();
    h_NJets->Write();


    fout->Write();
    fout->Close();
}

int main(int argc, char **argv)
{
    string inputFile = argv[1];
    string outputFile = argv[2];

    h_Muon1_pt->Sumw2();
    h_Muon1_eta->Sumw2();
    h_Muon2_pt->Sumw2();
    h_Muon2_eta->Sumw2();
    h_Muon_Muon_invariant_mass->Sumw2();    
    h_leading_lepton_pt->Sumw2();
    h_NJets->Sumw2();
    h_acopla_mumu->Sumw2();
   
    DataAnalysis(inputFile, outputFile);
}
