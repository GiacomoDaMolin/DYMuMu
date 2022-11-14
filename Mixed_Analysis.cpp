#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

// include user defined histograms and auxiliary macros
#include "Auxiliary.cpp"
#include "Histodef.cpp"
#include "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/roccor/RoccoR.cc"

// correctionlib
#include "correction.h"
using namespace std;
using correction::CorrectionSet;

#define MAX_ARRAY_SIZE 128
#define GEN_MAX_ARRAY_SIZE 1024

// function to calculate the weight for each event
// the weight is calculated as the product of luminosity and cross section of the process times the genWeight,
// LATER TO BE divided by the number of generated events OF ALL FILES OF THE DATASET(S)
double getWeight(double luminosity, double crossSection, Float_t genWeight, double SumWeights)
{
    return (luminosity * crossSection * genWeight); // / SumWeights;
}

void Mixed_Analysis(string inputFile, string ofile, double crossSection = -1, double IntLuminosity = 59.827879506, bool Signal = false)
{

    if (crossSection < 0. || IntLuminosity < 0.)
    {
        std::cout << "WARNING: crossection " << crossSection << " and Integrated luminosity " << IntLuminosity << endl;
    }

cout<<"Call completed!"<<endl;

    TFile *fin = TFile::Open(inputFile.c_str());
    TTree *trun = static_cast<TTree *>(fin->Get("Runs"));
    Long64_t genEventCount;
    Double_t genEventSumw;
    trun->SetBranchStatus("*", 0);
    trun->SetBranchStatus("genEventSumw", 1);
    trun->SetBranchStatus("genEventCount", 1);
    trun->SetBranchAddress("genEventSumw", &genEventSumw);
    trun->SetBranchAddress("genEventCount", &genEventCount);


    trun->GetEntry(0);

    TTree *tin = static_cast<TTree *>(fin->Get("Events"));

    // Set all branches to 0
    tin->SetBranchStatus("*", 0);
    // get the pt
    Float_t Muon_pt[MAX_ARRAY_SIZE],  Jet_pt[MAX_ARRAY_SIZE],Muon_mass[MAX_ARRAY_SIZE], Jet_mass[MAX_ARRAY_SIZE];
    Float_t Muon_eta[MAX_ARRAY_SIZE],  Jet_eta[MAX_ARRAY_SIZE], Muon_phi[MAX_ARRAY_SIZE], Jet_phi[MAX_ARRAY_SIZE];
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

    // get gen quantities

    Int_t Jet_genJetIdx[MAX_ARRAY_SIZE];
    Float_t GenPart_pt[GEN_MAX_ARRAY_SIZE];

    tin->SetBranchStatus("Jet_genJetIdx",1);
    tin->SetBranchAddress("Jet_genJetIdx",&Jet_genJetIdx);
    tin->SetBranchStatus("GenPart_pt",1);
    tin->SetBranchAddress("GenPart_pt",&GenPart_pt);

    // collect the trigger information
    Bool_t HLT_IsoMu24;
    tin->SetBranchStatus("HLT_IsoMu24", 1);
    tin->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24);

    // collect the triggger Ids
    Int_t Muon_charge[MAX_ARRAY_SIZE],Muon_nTrackerLayers[MAX_ARRAY_SIZE], Muon_genPartIdx[MAX_ARRAY_SIZE];
    Bool_t  Muon_triggerIdLoose[MAX_ARRAY_SIZE], Muon_tightId[MAX_ARRAY_SIZE];
    Float_t Muon_pfRelIso04_all[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Muon_tightId", 1);
    tin->SetBranchStatus("Muon_charge", 1);
    tin->SetBranchStatus("Muon_triggerIdLoose", 1);
    tin->SetBranchStatus("Muon_pfRelIso04_all", 1);
    tin->SetBranchStatus("Muon_nTrackerLayers", 1);
    tin->SetBranchStatus("Muon_genPartIdx", 1);
    tin->SetBranchAddress("Muon_tightId", &Muon_tightId);
    tin->SetBranchAddress("Muon_charge", &Muon_charge);
    tin->SetBranchAddress("Muon_triggerIdLoose", &Muon_triggerIdLoose);
    tin->SetBranchAddress("Muon_pfRelIso04_all", &Muon_pfRelIso04_all);
    tin->SetBranchAddress("Muon_nTrackerLayers", &Muon_nTrackerLayers);
    tin->SetBranchAddress("Muon_genPartIdx", &Muon_genPartIdx);

    // Jet tagging and ID, FlavB is the recomended one, DeepB was used by Anup
    Float_t Jet_btagDeepFlavB[MAX_ARRAY_SIZE];
    UInt_t nJet;
    Int_t Jet_jetId[MAX_ARRAY_SIZE], Jet_puId[MAX_ARRAY_SIZE],Jet_hadronFlavour[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Jet_btagDeepFlavB", 1);
    tin->SetBranchStatus("nJet", 1);
    tin->SetBranchStatus("Jet_jetId", 1);
    tin->SetBranchStatus("Jet_puId", 1);
    tin->SetBranchStatus("Jet_hadronFlavour", 1);
    tin->SetBranchAddress("nJet", &nJet);
    tin->SetBranchAddress("Jet_btagDeepFlavB", &Jet_btagDeepFlavB);   
    tin->SetBranchAddress("Jet_jetId", &Jet_jetId);
    tin->SetBranchAddress("Jet_puId", &Jet_puId);
    tin->SetBranchAddress("Jet_hadronFlavour", &Jet_hadronFlavour);

    // pu stuff
    Float_t N_pu_vertices;
    tin->SetBranchStatus("Pileup_nTrueInt", 1);
    tin->SetBranchAddress("Pileup_nTrueInt", &N_pu_vertices);

    // gen weight
    Float_t genWeight;
    tin->SetBranchStatus("genWeight", 1);
    tin->SetBranchAddress("genWeight", &genWeight);

    int non_matching_muon = 0, non_matching_electron = 0;
    int n_dropped = 0;
    int trigger_dropped = 0;
    UInt_t nEv = tin->GetEntries();
    unsigned int n_events = nEv;
    TLorentzVector *Muon1_p4 = new TLorentzVector();
    TLorentzVector *Muon2_p4 = new TLorentzVector();
    float Weight;

    // open correctionfiles
    
    string muon_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/muon_Z.json.gz";
    string jets_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/jet_jmar.json";
    string b_tag_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/btagging.json.gz";
    string pileup_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/puWeights.json.gz";

    
    auto muon_c_set = CorrectionSet::from_file(muon_json);
    auto jet_c_set = CorrectionSet::from_file(jets_json);
    auto btag_c_set = CorrectionSet::from_file(b_tag_json);
    auto pu_c_set = CorrectionSet::from_file(pileup_json);

    auto muon_trigger = muon_c_set->at("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight");
    auto muon_id = muon_c_set->at("NUM_TightID_DEN_genTracks");
    auto muon_iso = muon_c_set->at("NUM_TightRelIso_DEN_TightIDandIPCut");
    auto jet_pu = jet_c_set->at("PUJetID_eff");
    auto b_tag = btag_c_set->at("deepJet_mujets");
    auto b_mistag= btag_c_set->at("deepJet_incl"); //only for light jets
    auto pu_correction = pu_c_set->at("Collisions18_UltraLegacy_goldenJSON");


    TFile *fb_eff = new TFile("/afs/cern.ch/user/g/gdamolin/public/Beff_puLoose.root");
    TH2D * l_eff= static_cast<TH2D *>(fb_eff->Get("l_jets_tagged")); 
    TH2D * c_eff= static_cast<TH2D *>(fb_eff->Get("c_jets_tagged")); 
    TH2D * b_eff= static_cast<TH2D *>(fb_eff->Get("b_jets_tagged")); 
   

    RoccoR rc;
    rc.init("/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/roccor/RoccoR2018UL.txt");
      
    // save the histograms in a new File

    TFile *fout = new TFile(ofile.c_str(), "RECREATE");
    
    // create a new tree for the output
    TTree *tout = new TTree("tout", "tout");
    TTree *trun_out = new TTree("Run_out", "Run_out");
    Float_t invMass, muon1_eta, muon1_pt, muon2_pt, muon2_eta, dphi;

    tout->Branch("dphi", &dphi);
    tout->Branch("invMass", &invMass);
    tout->Branch("muon2_eta", &muon2_eta);
    tout->Branch("muon2_pt", &muon2_pt);
    tout->Branch("muon1_eta", &muon1_eta);
    tout->Branch("muon1_pt", &muon1_pt);
    tout->Branch("Weight", &Weight);

    trun_out->Branch("genEventSumw", &genEventSumw);
    trun_out->Branch("IntLumi", &IntLuminosity);
    trun_out->Branch("xs", &crossSection);
    trun_out->Branch("nEv", &n_events);

    trun_out->Fill(); // we already called trun->GetEntry(0);

    TRandom3* RndGen= new TRandom3();

    #pragma omp parallel for
    for (UInt_t i = 0; i <nEv; i++)
    {
        tin->GetEntry(i);
        if (i % 100000 == 0)
            std::cout << "Processing entry " << i << " of " << nEv << endl;
        // apply triggers

        if (!(HLT_IsoMu24)){
            trigger_dropped++;
            continue;
        };

        bool gotmuplus=false,gotmuminus=false;
        for (UInt_t j = 0; j < nMuon; j++){
            if ((Muon_pt[j]>27.|| ((gotmuplus||gotmuminus) && Muon_pt[j]>25.)) && abs(Muon_eta[j])<2.4 && Muon_tightId[j] && Muon_pfRelIso04_all[j] < 0.15){
                if (!gotmuplus && Muon_charge[j]==1){
			int NMCparticle=Muon_genPartIdx[j];
			double scmDT;
			if(NMCparticle>=0) {
				scmDT=rc.kSpreadMC(Muon_charge[j],Muon_pt[j],Muon_eta[j],Muon_phi[j],GenPart_pt[NMCparticle]);
				}
			else {
				scmDT=rc.kSmearMC(Muon_charge[j],Muon_pt[j],Muon_eta[j],Muon_phi[j],Muon_nTrackerLayers[j],RndGen->Rndm());
				}
			Muon1_p4->SetPtEtaPhiM(Muon_pt[j]*scmDT,Muon_eta[j],Muon_phi[j],Muon_mass[j]);
			gotmuplus=true;
			}
		if (!gotmuminus && Muon_charge[j]==-1){
			int NMCparticle=Muon_genPartIdx[j];
			double scmDT;
			if(NMCparticle>=0) {
				scmDT=rc.kSpreadMC(Muon_charge[j],Muon_pt[j],Muon_eta[j],Muon_phi[j],GenPart_pt[NMCparticle]);
				}
			else {
				scmDT=rc.kSmearMC(Muon_charge[j],Muon_pt[j],Muon_eta[j],Muon_phi[j],Muon_nTrackerLayers[j],RndGen->Rndm()); //TODO: Rndm may not be safe for stat application
				}
			Muon2_p4->SetPtEtaPhiM(Muon_pt[j]*scmDT,Muon_eta[j],Muon_phi[j],Muon_mass[j]);
			gotmuminus=true;
			}
            }
        }
 
	if(!(gotmuplus && gotmuminus)) {n_dropped++; continue;}

        Weight = getWeight(IntLuminosity, crossSection, genWeight, genEventSumw);
        Weight *= pu_correction->evaluate({N_pu_vertices, "nominal"}); 

        if(Muon1_p4->Pt()>Muon2_p4->Pt()) {Weight *= muon_trigger->evaluate({"2018_UL", abs(Muon1_p4->Eta()), Muon1_p4->Pt(), "sf"});}
	else {Weight *= muon_trigger->evaluate({"2018_UL", abs(Muon2_p4->Eta()), Muon2_p4->Pt(), "sf"});}
        Weight *= muon_id->evaluate({"2018_UL", abs(Muon1_p4->Eta()),Muon1_p4->Pt(), "sf"});
	Weight *= muon_id->evaluate({"2018_UL", abs(Muon2_p4->Eta()),Muon2_p4->Pt(), "sf"}); 
        Weight *= muon_iso->evaluate({"2018_UL", abs(Muon1_p4->Eta()), Muon1_p4->Pt(), "sf"}); 
	Weight *= muon_iso->evaluate({"2018_UL", abs(Muon2_p4->Eta()), Muon2_p4->Pt(), "sf"}); 
        
        Float_t jet_btag_deepFlav_wp = 0.2783;
        bool one_Bjet = false;
	//vectors for applying b-tag corrections
	vector<int> njet_in_collection, flavor;
	vector<bool> tagged; 
	double t_weight=1.;
	int Njets=0;
        for (size_t j = 0; j < nJet; j++)
        {
          if((abs(Jet_eta[j]) < 2.4) && Jet_pt[j]>25 && (Jet_jetId[j]==2 || Jet_jetId[j]==6)){
            //correction for pileupID
            int MC_pu = Jet_genJetIdx[j];
            float tempSF,tempEff;
            //if is pileUpjet
            if (MC_pu<0 ) {
		tempSF=1.;
            	tempEff= 0;
            	}
            //if is truly a jet
            else { if (Jet_pt[j]<=50){
            	     tempSF= jet_pu->evaluate({Jet_eta[j],Jet_pt[j],"nom", "L"});
            	     tempEff= jet_pu->evaluate({Jet_eta[j],Jet_pt[j],"MCEff", "L"});
		    }
            	}
            bool passesPUID=(Jet_puId[j]==4 || Jet_puId[j]==6 ||Jet_puId[j]==7);
            if(!(Jet_pt[j]>50 || passesPUID))	{t_weight*=(1-tempSF*tempEff)/(1-tempEff); }
            if((Jet_pt[j]>50 || passesPUID)) { 
             if(Jet_pt[j]<50) t_weight*=tempSF; //else you are in pT>50 case: apply no sf
             
             njet_in_collection.push_back(j);
             flavor.push_back(abs(Jet_hadronFlavour[j]));
             tagged.push_back((Jet_btagDeepFlavB[j] > jet_btag_deepFlav_wp));
	     Njets++;
            
             if (Jet_btagDeepFlavB[j] > jet_btag_deepFlav_wp)  {one_Bjet = true;}
             }//passes PUID
            }//passes kin
        }
        Weight*=t_weight; 
            
	for(int jj=0;jj<flavor.size();jj++){
		int convflav=flavor[jj];
		if (flavor[jj]<4) convflav==0;
		if (!(convflav==0 || convflav==4 || convflav==5)) {cout<<"Something weird in the flavor of jet"<<endl;}
		if(tagged[jj]){
			if (convflav!=0) Weight *= b_tag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});
			else  Weight *= b_mistag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});
			continue;}
		if(!tagged[jj]) {
			double Eff=1.;
			double SF=1;
			if (convflav!=0) SF=b_tag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});
			else SF=b_mistag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});
			
			//Get Eff
			if(convflav==0) {
				int bin =l_eff->FindBin(Jet_pt[njet_in_collection[jj]],abs(Jet_eta[njet_in_collection[jj]]));
				Eff=l_eff->GetBinContent(bin);
				}
			if(convflav==4) {
				int bin =c_eff->FindBin(Jet_pt[njet_in_collection[jj]],abs(Jet_eta[njet_in_collection[jj]]));
				Eff=c_eff->GetBinContent(bin);
				}
			if(convflav==5) {
				int bin =b_eff->FindBin(Jet_pt[njet_in_collection[jj]],abs(Jet_eta[njet_in_collection[jj]]));
				Eff=b_eff->GetBinContent(bin);
				}
			Weight*=(1-SF*Eff)/(1-Eff);
			}
		
		}
        if (one_Bjet){ n_dropped++; continue;}

        dphi=Muon1_p4->DeltaPhi(*Muon2_p4);

        h_NJets->Fill(Njets,Weight);
       
        muon1_pt = Muon1_p4->Pt();
        muon1_eta = Muon1_p4->Eta();
	muon2_pt = Muon2_p4->Pt();
        muon2_eta = Muon2_p4->Eta();

        h_Muon1_pt->Fill(muon1_pt,Weight);
        h_Muon1_eta->Fill(muon1_eta,Weight);
	h_Muon2_pt->Fill(muon2_pt,Weight);
        h_Muon2_eta->Fill(muon2_eta,Weight);
        h_acopla_mumu->Fill(M_PI-dphi,Weight);

        invMass = (*(Muon1_p4) + *(Muon2_p4)).M();
        h_Muon_Muon_invariant_mass->Fill(invMass,Weight);
	tout->Fill();
    }


    std::cout << "NeV = " << nEv << endl;
    std::cout << "trigger dropped = " << trigger_dropped << endl;
    std::cout << "selections dropped = " << n_dropped << endl; //remember the cross trigger in Data

    std::cout << "Fraction of events discarded by trigger = " << (trigger_dropped * 1. / nEv) << endl;
    int Rem_trigger=nEv-trigger_dropped; //remember the cross trigger in Data
    std::cout << "Fraction of events removed by selections = " << (n_dropped * 1. / Rem_trigger) << endl;
    std::cout << "Final number of events "<< Rem_trigger - n_dropped<<endl;

    tout->Write();
    trun_out->Write();
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
    double crossSection = atof(argv[3]);
    double IntLuminosity = atof(argv[4]);
    string boolstr = argv[5];
    bool Signal = (boolstr == "true");

    h_Muon1_pt->Sumw2();
    h_Muon1_eta->Sumw2();
    h_Muon2_pt->Sumw2();
    h_Muon2_eta->Sumw2();
    h_Muon_Muon_invariant_mass->Sumw2();    
    h_leading_lepton_pt->Sumw2();
    h_NJets->Sumw2();
    h_acopla_mumu->Sumw2();

    Mixed_Analysis(inputFile, outputFile, crossSection, IntLuminosity, Signal);

    return 0;
}
