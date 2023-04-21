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
#include "Auxiliary.cc"
#include "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/roccor/RoccoR.cc"

// correctionlib
#include "correction.h"
using namespace std;
using correction::CorrectionSet;

#define MAX_ARRAY_SIZE 128
#define GEN_MAX_ARRAY_SIZE 1024


double getWeight(double luminosity, double crossSection, Float_t genWeight, double SumWeights){
    return (luminosity * crossSection * genWeight);
}

void Mixed_Analysis(string inputFile, string ofile, double crossSection = -1, double IntLuminosity = 59.827879506, bool Data= false, bool systematics=false, string processname="A"){
    if (crossSection < 0. || IntLuminosity < 0.){
        std::cout << "WARNING: crossection " << crossSection << " and Integrated luminosity " << IntLuminosity << endl;
    }

cout<<"Call completed!"<<endl;
 if (Data) {systematics=false;}

    TFile *fin = TFile::Open(inputFile.c_str());
    TTree *trun;
    Long64_t genEventCount;
    Double_t genEventSumw;
    if(!Data){
	    trun = static_cast<TTree *>(fin->Get("Runs"));
	    trun->SetBranchStatus("*", 0);
	    trun->SetBranchStatus("genEventSumw", 1);
	    trun->SetBranchStatus("genEventCount", 1);
	    trun->SetBranchAddress("genEventSumw", &genEventSumw);
	    trun->SetBranchAddress("genEventCount", &genEventCount);

	    trun->GetEntry(0);
    }

    TTree *tin = static_cast<TTree *>(fin->Get("Events"));

    tin->SetBranchStatus("*", 0);
    Float_t Muon_pt[MAX_ARRAY_SIZE],  Jet_pt[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Muon_pt", 1);
    tin->SetBranchAddress("Muon_pt", &Muon_pt);
    tin->SetBranchStatus("Jet_pt", 1);
    tin->SetBranchAddress("Jet_pt", &Jet_pt);
    UInt_t nMuon;
    tin->SetBranchStatus("nMuon", 1);
    tin->SetBranchAddress("nMuon", &nMuon);
    Float_t Muon_eta[MAX_ARRAY_SIZE],  Jet_eta[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Muon_eta", 1);
    tin->SetBranchAddress("Muon_eta", &Muon_eta);
    tin->SetBranchStatus("Jet_eta", 1);
    tin->SetBranchAddress("Jet_eta", &Jet_eta);
    Float_t Muon_phi[MAX_ARRAY_SIZE], Jet_phi[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Muon_phi", 1);
    tin->SetBranchAddress("Muon_phi", &Muon_phi);
    tin->SetBranchStatus("Jet_phi", 1);
    tin->SetBranchAddress("Jet_phi", &Jet_phi);
    Float_t Muon_mass[MAX_ARRAY_SIZE], Jet_mass[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Muon_mass", 1);
    tin->SetBranchAddress("Muon_mass", &Muon_mass);
    tin->SetBranchStatus("Jet_mass", 1);
    tin->SetBranchAddress("Jet_mass", &Jet_mass);

    // get gen quantities
    Int_t Muon_genPartIdx[MAX_ARRAY_SIZE], Jet_hadronFlavour[MAX_ARRAY_SIZE];
    Int_t GenPart_pdgId[GEN_MAX_ARRAY_SIZE], GenPart_genPartIdxMother[GEN_MAX_ARRAY_SIZE], GenPart_statusFlags[GEN_MAX_ARRAY_SIZE], Jet_genJetIdx[MAX_ARRAY_SIZE];
    UChar_t Muon_genPartFlav[MAX_ARRAY_SIZE];
    UInt_t nGenPart;
    Float_t GenPart_pt[GEN_MAX_ARRAY_SIZE];
    Float_t genWeight, N_pu_vertices, L1PreFiringWeight_Nom;
    if(!Data){
	    tin->SetBranchStatus("Muon_genPartIdx", 1);
	    tin->SetBranchStatus("Muon_genPartFlav", 1);
	    tin->SetBranchStatus("GenPart_pdgId", 1);
	    tin->SetBranchStatus("GenPart_genPartIdxMother", 1);
	    tin->SetBranchStatus("nGenPart", 1);
	    tin->SetBranchStatus("Jet_genJetIdx",1);
	    tin->SetBranchStatus("GenPart_pt",1);
	    tin->SetBranchStatus("GenPart_statusFlags",1);
	    tin->SetBranchStatus("Jet_hadronFlavour", 1);
	    tin->SetBranchAddress("nGenPart", &nGenPart);
	    tin->SetBranchAddress("Muon_genPartIdx", &Muon_genPartIdx);
	    tin->SetBranchAddress("Muon_genPartFlav", &Muon_genPartFlav);
	    tin->SetBranchAddress("GenPart_pdgId", &GenPart_pdgId);
	    tin->SetBranchAddress("GenPart_genPartIdxMother", &GenPart_genPartIdxMother);
	    tin->SetBranchAddress("Jet_genJetIdx",&Jet_genJetIdx);
	    tin->SetBranchAddress("GenPart_pt",&GenPart_pt);
	    tin->SetBranchAddress("GenPart_statusFlags",&GenPart_statusFlags);
	    tin->SetBranchAddress("Jet_hadronFlavour", &Jet_hadronFlavour);
	    tin->SetBranchStatus("Pileup_nTrueInt", 1);
	    tin->SetBranchAddress("Pileup_nTrueInt", &N_pu_vertices);
	    tin->SetBranchStatus("genWeight", 1);
	    tin->SetBranchAddress("genWeight", &genWeight);
	    tin->SetBranchStatus("L1PreFiringWeight_Nom", 1);
	    tin->SetBranchAddress("L1PreFiringWeight_Nom", &L1PreFiringWeight_Nom);
    }

    // collect the trigger information
    Bool_t HLT_IsoMu24;
    tin->SetBranchStatus("HLT_IsoMu24", 1);
    tin->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24);

    // collect the triggger Ids
    Int_t Muon_charge[MAX_ARRAY_SIZE], Muon_nTrackerLayers[MAX_ARRAY_SIZE];
    Bool_t Muon_triggerIdLoose[MAX_ARRAY_SIZE], Muon_tightId[MAX_ARRAY_SIZE];
    Float_t Muon_pfRelIso04_all[MAX_ARRAY_SIZE], Muon_dxy[MAX_ARRAY_SIZE], Muon_dz[MAX_ARRAY_SIZE], Muon_ip3d[MAX_ARRAY_SIZE], Muon_sip3d[MAX_ARRAY_SIZE];

    tin->SetBranchStatus("Muon_tightId", 1);
    tin->SetBranchStatus("Muon_charge", 1);
    tin->SetBranchStatus("Muon_triggerIdLoose", 1);
    tin->SetBranchStatus("Muon_pfRelIso04_all", 1);
    tin->SetBranchStatus("Muon_nTrackerLayers", 1);
    tin->SetBranchAddress("Muon_tightId", &Muon_tightId);
    tin->SetBranchAddress("Muon_charge", &Muon_charge);
    tin->SetBranchAddress("Muon_triggerIdLoose", &Muon_triggerIdLoose);
    tin->SetBranchAddress("Muon_pfRelIso04_all", &Muon_pfRelIso04_all);
    tin->SetBranchAddress("Muon_nTrackerLayers", &Muon_nTrackerLayers);
    tin->SetBranchStatus("Muon_dxy", 1);
    tin->SetBranchStatus("Muon_dz", 1);
    tin->SetBranchStatus("Muon_ip3d", 1);
    tin->SetBranchStatus("Muon_sip3d", 1);
    tin->SetBranchAddress("Muon_dxy", &Muon_dxy);
    tin->SetBranchAddress("Muon_dz", &Muon_dz);
    tin->SetBranchAddress("Muon_ip3d", &Muon_ip3d);
    tin->SetBranchAddress("Muon_sip3d", &Muon_sip3d);

    // Jet tagging and ID, FlavB is the recomended one, DeepB was used by Anup
    Float_t Jet_btagDeepFlavB[MAX_ARRAY_SIZE], Jet_btagDeepB[MAX_ARRAY_SIZE];
    UInt_t nJet;
    Int_t Jet_jetId[MAX_ARRAY_SIZE], Jet_puId[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Jet_btagDeepB", 1);
    tin->SetBranchStatus("Jet_btagDeepFlavB", 1);
    tin->SetBranchStatus("nJet", 1);
    tin->SetBranchStatus("Jet_jetId", 1);
    tin->SetBranchStatus("Jet_puId", 1);
    tin->SetBranchAddress("nJet", &nJet);
    tin->SetBranchAddress("Jet_btagDeepFlavB", &Jet_btagDeepFlavB);
    tin->SetBranchAddress("Jet_btagDeepB", &Jet_btagDeepB);   
    tin->SetBranchAddress("Jet_jetId", &Jet_jetId);
    tin->SetBranchAddress("Jet_puId", &Jet_puId);
   

    int non_matching_muon = 0, non_matching_electron = 0, evenottrigMatch=0,n_dropped = 0, trigger_dropped = 0;
    UInt_t nEv = tin->GetEntries();
    unsigned int n_events = nEv;
    TLorentzVector *Muon1_p4 = new TLorentzVector();
    TLorentzVector *Muon2_p4 = new TLorentzVector();

    // allow pt, inv mass, and eta to be stored in a Branch
    Float_t invMass, muon1_eta, muon1_pt, muon2_pt, muon2_eta, Acoplanarity;

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
   
    TRandom3 * RndGen=new TRandom3();
    RoccoR rc;
    rc.init("/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/roccor/RoccoR2018UL.txt");
     
    TFile *fout = new TFile(ofile.c_str(), "RECREATE");
    
    // create a new tree for the output
    TTree *tout = new TTree("tout", "tout");
    TTree *trun_out;
    
    // set the branches for the output tree
    tout->Branch("Acoplanarity", &Acoplanarity);
    tout->Branch("invMass", &invMass);
    tout->Branch("muon2_eta", &muon2_eta);
    tout->Branch("muon2_pt", &muon2_pt);
    tout->Branch("muon1_eta", &muon1_eta);
    tout->Branch("muon1_pt", &muon1_pt);

    int Nloose = 0, Nmedium = 0, Ntight = 0, JetsNotB=0, Njets=0;

    tout->Branch("Njets", &Njets);
    tout->Branch("Nloose", &Nloose);
    tout->Branch("Nmedium", &Nmedium);
    tout->Branch("Ntight", &Ntight);
    tout->Branch("JetNotB", &JetsNotB);

    if(!Data){
	    trun_out = new TTree("Run_out", "Run_out");
	    trun_out->Branch("genEventSumw", &genEventSumw);
	    trun_out->Branch("IntLumi", &IntLuminosity);
	    trun_out->Branch("xs", &crossSection);
	    trun_out->Branch("nEv", &n_events);
	    trun_out->Fill();
    }
    vector<string> observables2,observables3;
    if(processname=="TT2LL") {
	processname="TT2ll";
    	observables2= {"TT2lT_Muon1_pt","TT2lT_Acoplanarity","TT2lT_Invariant_Mass","TT2lT_Njets"}; 
    	observables3= {"TT2TT_Muon1_pt","TT2TT_Acoplanarity","TT2TT_Invariant_Mass","TT2TT_Njets"}; 
    }
    if(processname=="DY") {
	processname="DY2ee";
    	observables2= {"DY2mm_Muon1_pt","DY2mm_Acoplanarity","DY2mm_Invariant_Mass","DY2mm_Njets"}; 
    	observables3= {"DY2TT_Muon1_pt","DY2TT_Acoplanarity","DY2TT_Invariant_Mass","DY2TT_Njets"}; 
    }  

    vector<string> observables= {processname+"_Muon1_pt",processname+"_Acoplanarity",processname+"_Invariant_Mass",processname+"_Njets"}; 
    vector<string> systs= {"Muon_Id","Muon_Iso","Muon_trigger"};
    vector<string> shift={"","Up","Down"};
  

    TH1D* temp=NULL;
    vector<TH1D*>  Histos; //defaulted to LL if TT2LL && ee if DY
    vector<TH1D*>  HistosTL,HistosTT; //correspond to mumu && TauTau in DY
   if(systematics){ 
	for(int i=0;i<observables.size()*(systs.size()*2+1);i++) {Histos.push_back(temp); 
		if(observables2.size()>0){HistosTL.push_back(temp); HistosTT.push_back(temp);}
	}
     }
    else for(int i=0;i<observables.size();i++) {Histos.push_back(temp); if(observables2.size()>0){HistosTL.push_back(temp); HistosTT.push_back(temp);}}

    //get histos nominal values
    temp= new TH1D(observables[0].c_str(),observables[0].c_str(),40,25,205);
    Histos[0] = (TH1D*)temp->Clone();
    temp = new TH1D(observables[1].c_str(),observables[1].c_str(),40,0, 2*M_PI);
    Histos[1] = (TH1D*)temp->Clone();
    temp= new TH1D(observables[2].c_str(),observables[2].c_str(),40, 12, 412);
    Histos[2] = (TH1D*)temp->Clone();
    temp = new TH1D(observables[3].c_str(),observables[3].c_str(),12,0,12);
    Histos[3] = (TH1D*)temp->Clone();
    if(HistosTT.size()>0){
	temp= new TH1D(observables2[0].c_str(),observables2[0].c_str(),40,25,205);
	HistosTL[0] = (TH1D*)temp->Clone();
	temp = new TH1D(observables2[1].c_str(),observables2[1].c_str(),40,0, 2*M_PI);
	HistosTL[1] = (TH1D*)temp->Clone();
	temp = new TH1D(observables2[2].c_str(),observables2[2].c_str(),40, 12, 412);
	HistosTL[2]= (TH1D*)temp->Clone(); 
	temp= new TH1D(observables2[3].c_str(),observables2[3].c_str(),12,0,12);
	HistosTL[3] = (TH1D*)temp->Clone();

	temp= new TH1D(observables3[0].c_str(),observables3[0].c_str(),40,25,205);
	HistosTT[0] = (TH1D*)temp->Clone();
	temp = new TH1D(observables3[1].c_str(),observables3[1].c_str(),40,0, 2*M_PI);
	HistosTT[1] = (TH1D*)temp->Clone();
	temp = new TH1D(observables3[2].c_str(),observables3[2].c_str(),40, 12, 412);
	HistosTT[2]= (TH1D*)temp->Clone(); 
	temp= new TH1D(observables3[3].c_str(),observables3[3].c_str(),12,0,12);
	HistosTT[3] = (TH1D*)temp->Clone();
    }

    int auxindex=observables.size()-1;
    int idx2=auxindex,idx3=auxindex;
    //here loop on systematics, clone the histo and do your things
    if(systematics){
    	for(int i=0; i<systs.size(); i++){
	    	 for(int j=0;j<observables.size();j++){
		    	 Histos[++auxindex]=cloneDims1d(Histos[j],(systs[i]+"Up").c_str());	
	    	}
		for(int j=0;j<observables.size();j++){
		    	 Histos[++auxindex]=cloneDims1d(Histos[j],(systs[i]+"Down").c_str());	
	    	}
    	}
	if(HistosTT.size()>0){
		for(int i=0; i<systs.size(); i++){
			for(int j=0;j<observables2.size();j++){
				    	 HistosTL[++idx2]=cloneDims1d(HistosTL[j],(systs[i]+"Up").c_str());	
			    	}
			for(int j=0;j<observables2.size();j++){
				    	 HistosTL[++idx2]=cloneDims1d(HistosTL[j],(systs[i]+"Down").c_str());	
			    	}
			for(int j=0;j<observables3.size();j++){
				    	 HistosTT[++idx3]=cloneDims1d(HistosTT[j],(systs[i]+"Up").c_str());	
			    	}
			for(int j=0;j<observables3.size();j++){
				    	 HistosTT[++idx3]=cloneDims1d(HistosTT[j],(systs[i]+"Down").c_str());	
			    	}
		}
	}
   }


   if(systematics && auxindex!=observables.size()*(systs.size()*2+1)-1) cout<<"Something may go terribly wrong here!"<< auxindex<<" vs "<<observables.size()*(systs.size()*2+1)-1 <<endl;
	int cat=0; //category of final state index
    #pragma omp parallel for
    for (UInt_t i = 0; i <nEv; i++){
	cat=0;
        tin->GetEntry(i);
        if (i % 100000 == 0) std::cout << "Processing entry " << i << " of " << nEv << endl;
        
        if (!(HLT_IsoMu24)){trigger_dropped++; continue;}
        
        int Master_index=1; //the first entry is the nominal one!
	vector<double> VecWeights(systs.size()*2+1,1.);
	
        bool gotmuplus=false,gotmuminus=false;
	int mu1idx=-1, mu2idx=-1;
        for (UInt_t j = 0; j < nMuon; j++){
            if (abs(Muon_eta[j])<2.4 && Muon_tightId[j] && Muon_pfRelIso04_all[j] < 0.15){
		if(!Data){
			int NMCparticle=Muon_genPartIdx[j];
			double scmDT;
			if(NMCparticle>=0) {scmDT=rc.kSpreadMC(Muon_charge[j],Muon_pt[j],Muon_eta[j],Muon_phi[j],GenPart_pt[NMCparticle]);}
			else {scmDT=rc.kSmearMC(Muon_charge[j],Muon_pt[j],Muon_eta[j],Muon_phi[j],Muon_nTrackerLayers[j],RndGen->Rndm());}
			Muon_pt[j]*=scmDT;
			}
		else if(Data){
			double scmDT=rc.kScaleDT(Muon_charge[j],Muon_pt[j],Muon_eta[j],Muon_phi[j]);
			Muon_pt[j]*=scmDT;
			}
		if(Muon_pt[j]>27.|| ((gotmuplus||gotmuminus) && Muon_pt[j]>25.)){
		        if (!gotmuplus && Muon_charge[j]==1){
				Muon1_p4->SetPtEtaPhiM(Muon_pt[j],Muon_eta[j],Muon_phi[j],Muon_mass[j]);
				if(gotmuminus) {if(Muon1_p4->DeltaR(*Muon2_p4)<0.4) {continue;} }
				gotmuplus=true;
				mu1idx=j;
				}
			if (!gotmuminus && Muon_charge[j]==-1){
				Muon2_p4->SetPtEtaPhiM(Muon_pt[j],Muon_eta[j],Muon_phi[j],Muon_mass[j]);
				if(gotmuplus) {if(Muon1_p4->DeltaR(*Muon2_p4)<0.4) {continue;} }
				gotmuminus=true;
				mu2idx=j;
				}
		}
            }
        }
        
        if(!(gotmuplus && gotmuminus)) {n_dropped++; continue;}

        Float_t jet_btag_deepFlav_wp = 0.2783;
        bool one_Bjet = false;
        int id_m_jet = -1;
        Nloose = 0, Nmedium = 0, Ntight = 0, JetsNotB=0, Njets=0;
	//vectors for applying b-tag corrections
	vector<int> njet_in_collection; vector<int> flavor; vector<bool> tagged;
	double t_weight=1.;
        for (size_t j = 0; j < nJet; j++) {
            if((abs(Jet_eta[j]) < 2.4) && Jet_pt[j]>30 && (Jet_jetId[j]==2 || Jet_jetId[j]==6)){
	    TLorentzVector *Tjet_p4 = new TLorentzVector();
	    Tjet_p4->SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
	    if((Tjet_p4->DeltaR(*Muon1_p4)<0.4) || (Tjet_p4->DeltaR(*Muon2_p4)<0.4)) {delete Tjet_p4; continue;}
	    else {delete Tjet_p4;}
	
            //correction for pileupID
            float tempSF=1.,tempEff=0.;
            if(!Data){
		    int MC_pu = Jet_genJetIdx[j];
		    if (MC_pu<0 ) {tempSF=1.; tempEff= 0;}
		    //if is truly a jet
		    else { if (Jet_pt[j]<=50){
		    	     tempSF= jet_pu->evaluate({Jet_eta[j],Jet_pt[j],"nom", "L"});
		    	     tempEff= jet_pu->evaluate({Jet_eta[j],Jet_pt[j],"MCEff", "L"});
			    }
		    	}
            }
            bool passesPUID=(Jet_puId[j]>=4);
		
            if(!(Jet_pt[j]>50 || passesPUID ))	{t_weight*=(1-tempSF*tempEff)/(1-tempEff); }
            if((Jet_pt[j]>50 || passesPUID)) { 
             if(Jet_pt[j]<=50) t_weight*=tempSF; //else you are in pT>50 case: apply no sf
              //correction for b-tag
              njet_in_collection.push_back(j);
              flavor.push_back(abs(Jet_hadronFlavour[j]));
              tagged.push_back((Jet_btagDeepFlavB[j] > jet_btag_deepFlav_wp));
	      Njets++;
        
	      if (Jet_btagDeepFlavB[j] < 0.0490) JetsNotB++;
	      if (Jet_btagDeepFlavB[j] > 0.0490) Nloose++;
              if (Jet_btagDeepFlavB[j] > 0.2783) Nmedium++;
              if (Jet_btagDeepFlavB[j] > 0.71) Ntight++;
            
              if (Jet_btagDeepFlavB[j] > jet_btag_deepFlav_wp){
                 if (!one_Bjet){
                      	one_Bjet = true;
                      	id_m_jet = j;
                      }
               }
            }//end if(jetpt>50 !!puid==7)
          }//end kinematic if
        }//end for

        if (one_Bjet){ n_dropped++;  continue;}

	
	if (processname=="TT2ll") {cat=getCategoryTT(GenPart_pdgId,GenPart_genPartIdxMother,nGenPart);}

	if (processname=="DY2ee") {cat=getCategoryDY(GenPart_pdgId,GenPart_genPartIdxMother,nGenPart);}
         
        if(!Data){ 
		 
		Weight = getWeight(IntLuminosity, crossSection, genWeight, genEventSumw);
		Weight *=  getTopPtWeight(GenPart_pdgId,GenPart_statusFlags,GenPart_pt,nGenPart);
		Weight *= pu_correction->evaluate({N_pu_vertices, "nominal"});
		Weight*= L1PreFiringWeight_Nom;
		
		for(auto &k: VecWeights) k*=Weight;
				
		for(int j=0;j<VecWeights.size();j++){
			if(systematics){
				if(j==Master_index) {VecWeights[j]*= muon_id->evaluate({"2018_UL", abs(Muon1_p4->Eta()), Muon1_p4->Pt(), "systup"}); 
							VecWeights[j]*= muon_id->evaluate({"2018_UL", abs(Muon2_p4->Eta()), Muon2_p4->Pt(), "systup"}); continue;}
				if(j==Master_index+1) {VecWeights[j]*= muon_id->evaluate({"2018_UL", abs(Muon1_p4->Eta()), Muon1_p4->Pt(), "systdown"});
							VecWeights[j]*= muon_id->evaluate({"2018_UL", abs(Muon2_p4->Eta()), Muon2_p4->Pt(), "systdown"}); continue;}
				else {VecWeights[j] *= muon_id->evaluate({"2018_UL", abs(Muon1_p4->Eta()), Muon1_p4->Pt(), "sf"});
					VecWeights[j]*= muon_id->evaluate({"2018_UL", abs(Muon2_p4->Eta()), Muon2_p4->Pt(), "sf"});}
			}
			else {VecWeights[j] *= muon_id->evaluate({"2018_UL", abs(Muon1_p4->Eta()), Muon1_p4->Pt(), "sf"});
				VecWeights[j]*= muon_id->evaluate({"2018_UL", abs(Muon2_p4->Eta()), Muon2_p4->Pt(), "sf"});}
		} 
		Master_index+=2;
		for(int j=0;j<VecWeights.size();j++){
			if(systematics){
				if(j==Master_index) {VecWeights[j]*= muon_iso->evaluate({"2018_UL", abs(Muon1_p4->Eta()), Muon1_p4->Pt(), "systup"}); 
							VecWeights[j]*= muon_iso->evaluate({"2018_UL", abs(Muon2_p4->Eta()), Muon2_p4->Pt(), "systup"}); continue;}
				if(j==Master_index+1) {VecWeights[j]*= muon_iso->evaluate({"2018_UL", abs(Muon1_p4->Eta()), Muon1_p4->Pt(), "systdown"}); 
							VecWeights[j]*= muon_iso->evaluate({"2018_UL", abs(Muon2_p4->Eta()), Muon2_p4->Pt(), "systdown"}); continue;}
				else {VecWeights[j]*= muon_iso->evaluate({"2018_UL", abs(Muon1_p4->Eta()), Muon1_p4->Pt(), "sf"});
					VecWeights[j]*= muon_iso->evaluate({"2018_UL", abs(Muon2_p4->Eta()), Muon2_p4->Pt(), "sf"});}
			}
			else {VecWeights[j]*= muon_iso->evaluate({"2018_UL", abs(Muon1_p4->Eta()), Muon1_p4->Pt(), "sf"});
				VecWeights[j]*= muon_iso->evaluate({"2018_UL", abs(Muon2_p4->Eta()), Muon2_p4->Pt(), "sf"});}
		} 
		Master_index+=2;
		float temp_pt,temp_eta;
		if(Muon1_p4->Pt()>Muon2_p4->Pt()){temp_pt=Muon1_p4->Pt(); temp_eta=Muon1_p4->Eta();}
		else {temp_pt=Muon2_p4->Pt(); temp_eta=Muon2_p4->Eta();}
		for(int j=0;j<VecWeights.size();j++){
			if(systematics){
				if(j==Master_index) {VecWeights[j]*= muon_trigger->evaluate({"2018_UL", abs(temp_eta), temp_pt, "systup"}); continue;}
				if(j==Master_index+1) {VecWeights[j]*= muon_trigger->evaluate({"2018_UL", abs(temp_eta), temp_pt, "systdown"}); continue;}
				else VecWeights[j]*= muon_trigger->evaluate({"2018_UL", abs(temp_eta),temp_pt, "sf"});
			}
			else VecWeights[j]*= muon_trigger->evaluate({"2018_UL", abs(temp_eta), temp_pt, "sf"});
		} 
		Master_index+=2;
		
		for(auto &k: VecWeights) k*=t_weight; 
		double Weight2=1.;
		for(int jj=0;jj<flavor.size();jj++){
			int convflav=flavor[jj];
			if (flavor[jj]<4) convflav==0;
			if (!(convflav==0 || convflav==4 || convflav==5)) {cout<<"Something weird in the flavor of jet"<<endl;}
			if(tagged[jj]){
				if (convflav!=0) {Weight2 *= b_tag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});}
				else  {Weight2 *= b_mistag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});}
				continue;}

			if(!tagged[jj]) {
				double Eff=1.;
				double SF=1.;
				if (convflav!=0) SF=b_tag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});
				else SF=b_mistag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});
				
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
				Weight2*=(1-SF*Eff)/(1-Eff);
				}
			
			}
		for(auto &k: VecWeights) k*=Weight2;
	}
	
	muon1_pt = Muon1_p4->Pt();
        muon1_eta = Muon1_p4->Eta();
        muon2_pt = Muon2_p4->Pt();
        muon2_eta = Muon2_p4->Eta();
        invMass = (*(Muon1_p4) + *(Muon2_p4)).M();
	Acoplanarity=M_PI-(Muon1_p4->DeltaPhi(*Muon2_p4));
	
	for(int k=0;k<VecWeights.size();k++){
		if(cat==0 || cat==2){
			Histos[k*observables.size()] ->Fill(muon1_pt,VecWeights[k]);
	    		Histos[k*observables.size()+1] ->Fill(Acoplanarity,VecWeights[k]);
	    		Histos[k*observables.size()+2] ->Fill(invMass,VecWeights[k]);
	    		Histos[k*observables.size()+3] ->Fill(Njets,VecWeights[k]);
			}
		if(HistosTT.size()>0){
		if(cat==3){
			HistosTL[k*observables.size()] ->Fill(muon1_pt,VecWeights[k]);
	    		HistosTL[k*observables.size()+1] ->Fill(Acoplanarity,VecWeights[k]);
	    		HistosTL[k*observables.size()+2] ->Fill(invMass,VecWeights[k]);
	    		HistosTL[k*observables.size()+3] ->Fill(Njets,VecWeights[k]);
			}
		if(cat==4){
			HistosTT[k*observables.size()] ->Fill(muon1_pt,VecWeights[k]);
	    		HistosTT[k*observables.size()+1] ->Fill(Acoplanarity,VecWeights[k]);
	    		HistosTT[k*observables.size()+2] ->Fill(invMass,VecWeights[k]);
	    		HistosTT[k*observables.size()+3] ->Fill(Njets,VecWeights[k]);
			}
		 }
		if(!systematics){break;}
		}

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
    if(!Data) trun_out->Write();
    for(auto &k: Histos) HistWrite(k);
    if( HistosTT.size()>0){
	 for(auto &k: HistosTL) HistWrite(k);
	 for(auto &k: HistosTT) HistWrite(k);
    }

    fout->Close();
}

int main(int argc, char **argv){
    string inputFile = argv[1];
    string outputFile = argv[2];
    double crossSection = atof(argv[3]);
    double IntLuminosity = atof(argv[4]);
    string boolstr = argv[5];
    bool Data = (boolstr == "true" || boolstr == "True");
    string boolstr2 = argv[6];
    bool systematics = (boolstr2 == "true" || boolstr2 == "True");
    string processname = argv[7];

cout<<"Process "<<processname<<endl;

    Mixed_Analysis(inputFile, outputFile, crossSection, IntLuminosity, Data, systematics,processname);

    return 0;
}
