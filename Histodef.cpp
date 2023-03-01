#ifndef Histodef_cpp
#define Histodef_cpp
#include "TH1.h"

    // Define the histogram objects
    TH1D* h_Muon1_pt = new TH1D("h_Muon1_pt","h_Muon_positive_pt",40,0,200);
    TH1D* h_Muon1_eta = new TH1D("h_Muon1_eta","h_Muon_positive_eta",40,-4,4);
    TH1D* h_Muon2_pt = new TH1D("h_Muon2_pt","h_Muon_negative_pt",40,0,200);
    TH1D* h_Muon2_eta = new TH1D("h_Muon2_eta","h_Muon_negative_eta",40,-4,4);

    TH1D* h_Muon_Muon_invariant_mass = new TH1D("Muon_Muon_invariant_mass", "Muon_Muon_invariant_mass", 50, 12, 212);
    
    TH1D* h_leading_lepton_pt = new TH1D("leading_lepton_pt", "leading_lepton_pt", 36, 20, 200);
 
    TH1D *h_NJets = new TH1D("N_Jets","N_Jets",12,0,12);

    TH1D * h_acopla_mumu =new TH1D("h_Acopla_emu","h_Acopla_emu",30,0, 2*M_PI);
    TH1D* h_mu_3dsig = new TH1D("h_mu_3dsig","h_mu_3dsig",30,0,30);
    TH1D* h_mu_3d = new TH1D("h_mu_3d","h_mu_3d",30,0,0.06);
    TH1D * h_mu_dxy =new TH1D("h_mu_dxy","h_mu_dxy",30,0, 0.03);


#endif // Histodef_cpp
