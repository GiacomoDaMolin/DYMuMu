#ifndef Histodef_cpp
#define Histodef_cpp
#include "TH1.h"

    // Define the histogram objects
    TH1D* h_Muon1_pt = new TH1D("h_Muon1_pt","h_Muon_positive_pt",100,0,200);
    TH1D* h_Muon1_eta = new TH1D("h_Muon1_eta","h_Muon_positive_eta",100,-5,5);
    TH1D* h_Muon2_pt = new TH1D("h_Muon2_pt","h_Muon_negative_pt",100,0,200);
    TH1D* h_Muon2_eta = new TH1D("h_Muon2_eta","h_Muon_negative_eta",100,-5,5);

    TH1D* h_Muon_Muon_invariant_mass = new TH1D("Muon_Muon_invariant_mass", "Muon_Muon_invariant_mass", 50, 12, 212);
    
    TH1D* h_leading_lepton_pt = new TH1D("leading_lepton_pt", "leading_lepton_pt", 45, 20, 200);
 
    TH1D *h_NJets = new TH1D("N_Jets","N_Jets",12,0,12);

    TH1D * h_acopla_mumu =new TH1D("h_Acopla_emu","h_Acopla_emu",50,0, 2*M_PI);


#endif // Histodef_cpp
