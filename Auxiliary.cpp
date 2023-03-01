#ifndef Auxiliary_cpp
#define Auxiliary_cpp
#include <iostream>
// inputs are: size of arrays GenID and GenPar
// the pointer at the start of the array of GenParticles_PDGID
// the pointer at the start of the array of GenParticles_Parent
// initialID: last index of the particle (ex: Electron_genPartIdx[j]=initialID)

// should be called from func as isFromW(MAX_ARRAY_SIZE,GenPart_pdgId,GenPart_genPartIdxMother,Electron_genPartIdx[j])
//arguments: (array of PdgID of Gen Particles (Int_t), array of GenParent of Gen Particles (Int_t, contains the index to the parent of the selected GenParticle
// InitialID  is the index of the starting muon in the GenParticles array)
bool isFromW(int size, Int_t *GenId, Int_t *GenParent, int initialID)
{
	if (initialID < 0)
	{
		return false;
	}
	// retrieve first PDG ID number
	int startPdg = GenId[initialID];
	int newID = initialID, newPdg = startPdg;
	// look for the parent; if the parent is of same PDGID of starting particle, iterate until parent is different particle
	while (newPdg == startPdg)
	{
		if (newID > size)
		{
			std::cout << "WARNING: index " << newID << " exceeding max size " << size << std::endl;
		}
		newID = GenParent[newID];
		newPdg = GenId[newID];
		if (abs(newPdg) == 24)
			return true;
	}
	return false;
}

bool isFromTau(int size, Int_t *GenId, Int_t *GenParent, int initialID)
{
	if (initialID < 0)
	{
		return false;
	}
	// retrieve first PDG ID number
	int startPdg = GenId[initialID];
	int newID = initialID, newPdg = startPdg;
	// look for the parent; if the parent is of same PDGID of starting particle, iterate until parent is different particle
	while (newPdg == startPdg)
	{
		if (newID > size)
		{
			std::cout << "WARNING: index " << newID << " exceeding max size " << size << std::endl;
		}
		newID = GenParent[newID];
		newPdg = GenId[newID];
		//std::cout<< "PDG corrente "<< newPdg<<std::endl;
		if (abs(newPdg) == 15)
			{return true;}
	}
	return false;
}

double InvertPhi(double phi){
 double invphi=phi+M_PI;
 if (invphi>M_PI){invphi=invphi-2*M_PI;}
 return invphi;
}


void HistWrite(){
h_Muon_Muon_invariant_mass->SetBinContent(h_Muon_Muon_invariant_mass->GetNbinsX(), h_Muon_Muon_invariant_mass->GetBinContent(h_Muon_Muon_invariant_mass->GetNbinsX()) + h_Muon_Muon_invariant_mass->GetBinContent(h_Muon_Muon_invariant_mass->GetNbinsX() + 1));
h_Muon1_pt->SetBinContent(h_Muon1_pt->GetNbinsX(), h_Muon1_pt->GetBinContent(h_Muon1_pt->GetNbinsX()) + h_Muon1_pt->GetBinContent(h_Muon1_pt->GetNbinsX() + 1));
h_Muon2_pt->SetBinContent(h_Muon2_pt->GetNbinsX(), h_Muon2_pt->GetBinContent(h_Muon2_pt->GetNbinsX()) + h_Muon2_pt->GetBinContent(h_Muon1_pt->GetNbinsX() + 1));
    h_Muon1_eta->Write();
    h_Muon1_pt->Write();
    h_Muon2_eta->Write();
    h_Muon2_pt->Write();

    h_Muon_Muon_invariant_mass->Write();
    h_acopla_mumu->Write();
    h_NJets->Write();
    h_mu_3dsig->Write();
    h_mu_3d->Write();
    h_mu_dxy->Write();
}

void HistPrep(){
h_Muon1_pt->Sumw2();
    h_Muon1_eta->Sumw2();
    h_Muon2_pt->Sumw2();
    h_Muon2_eta->Sumw2();
    h_Muon_Muon_invariant_mass->Sumw2();    
    h_leading_lepton_pt->Sumw2();
    h_NJets->Sumw2();
    h_acopla_mumu->Sumw2();
    h_mu_3dsig->Sumw2();
    h_mu_3d->Sumw2();
    h_mu_dxy->Sumw2();
   
}
#endif // Auxiliary_cpp
