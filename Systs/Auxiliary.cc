#ifndef Auxiliary_cpp
#define Auxiliary_cpp
#include <iostream>
#include<vector>
#include<string>
#include<TH1D.h>
#include<string>
using std::vector;
using std::string;
using std::cout; using std::endl;
using json = nlohmann::json;

bool lumicheck(json& goldenjson,UInt_t run,UInt_t luminosityBlock){
 string run_s = std::to_string(run);
 //string lumi_s = std::to_string(luminosityBlock);
 if (goldenjson[run_s]==nullptr) return true; //run not in golden json

 for(int i=0;i<goldenjson[run_s].size();i++){	
	if(luminosityBlock>=goldenjson[run_s][i][0] && luminosityBlock<=goldenjson[run_s][i][1])  return false; //if run found, look if lumisec is in good interval
 }
 
 return true; //if lumiblock is not found, it's not in the json --> discard the event

}


float getTopPtWeight(Int_t * pdgId,Int_t *statusFlags,Float_t * pt, Int_t Ngen) {
    float wgt=1.0;
    vector<Int_t> topquarks_lastcopy;
    for(Int_t i=0;i<Ngen;i++){
	if(abs(pdgId[i])==6 && ((statusFlags[i]/8192)%2))
		topquarks_lastcopy.push_back(i);
	}
    if(topquarks_lastcopy.size()!=2) return wgt;

    //NNLO / NLO SF parameterization
    for(auto i : topquarks_lastcopy) {
      wgt *= 0.103 * exp(-0.0118 * pt[i]) - 0.000134 * pt[i] + 0.973;
    }
    return sqrt(wgt);
  }

TH1D* cloneDims1d(TH1D* hist, string newname){
  if (hist == NULL) return NULL;
  TH1D* cloneHist = new TH1D(Form("%s_%s", hist->GetName(),newname.c_str()), 
                             Form("%s_%s", hist->GetName(),newname.c_str()), 
                             hist->GetNbinsX(),
                             hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
                             
  cloneHist->Sumw2();
  return cloneHist;
}

void HistWrite(TH1D *a){
 a->SetBinContent(a->GetNbinsX(), a->GetBinContent(a->GetNbinsX()) + a->GetBinContent(a->GetNbinsX() + 1));
 a->Write();
}

double InvertPhi(double phi){
 double invphi=phi+M_PI;
 if (invphi>M_PI){invphi=invphi-2*M_PI;}
 return invphi;
}

bool isFromtop(int size, Int_t *GenId, Int_t *GenParent, int initialID){
	if (initialID < 0){return false;}
	// retrieve first PDG ID number
	int startPdg = GenId[initialID];
	int newID = initialID, newPdg = startPdg;
	// look for the parent; if the parent is of same PDGID of starting particle, iterate until parent is different particle
	while (newPdg == startPdg){
		if (newID > size){std::cout << "WARNING: index " << newID << " exceeding max size " << size << std::endl;}
		newID = GenParent[newID];
		newPdg = GenId[newID];
		if (abs(newPdg) == 6)
			return true;
	}
	return false;
}

int getCategoryTT(Int_t * pdgId,Int_t *IdxMother, Int_t Ngen) {
   int cat=0;
   bool foundWp=false,foundWm=false;
   bool TauWp=false,TauWm=false;
   //first 
   for(int i=Ngen-1;i>-1;i--){ //start from last
   	if(pdgId[i]==24 && isFromtop(Ngen, pdgId, IdxMother, i) && !foundWp){
   	foundWp=true;
   	for(int j=i;j<Ngen;j++){
   		if(abs(pdgId[j])==15 && IdxMother[j]==i) {TauWp=true; break;}
   		if((abs(pdgId[j])==11 || abs(pdgId[j])==13) && IdxMother[j]==i) {TauWp=false; break;}
   		if(j==Ngen-1) cout<<"Daughter not found!"<<endl;
		}
   	}
   	if(pdgId[i]==-24 && isFromtop(Ngen, pdgId, IdxMother, i) && !foundWm){
   	foundWm=true;
   	for(int j=i;j<Ngen;j++){
   		if(abs(pdgId[j])==15 && IdxMother[j]==i) {TauWm=true; break;}
   		if((abs(pdgId[j])==11 || abs(pdgId[j])==13) && IdxMother[j]==i) {TauWm=false; break;}
   		if(j==Ngen-1) cout<<"Daughter not found!"<<endl;
		}
   	}
   if(foundWp&&foundWm) break;
   }
  if(TauWp && TauWm) cat=4; //-->TauTau
  else if(!TauWp && !TauWm) cat=2; //->LL
  else if((TauWp && !TauWm)||(!TauWp && TauWm)) cat=3; //->LTau

 return cat;
}

int getCategoryDY(Int_t * pdgId,Int_t *IdxMother, Int_t Ngen) {
   int cat=0;
   bool foundZ=false;
   bool recover=false;
   //fZ is always index==2, daughter of index 0
   for(int i=0;i<Ngen;i++){
   	if((pdgId[i]==22 || pdgId[i]==23) && (IdxMother[i]==0 || IdxMother[i]==1)){
   	foundZ=true;
   	int idx=i,startPdg=pdgId[i];
   	for(int j=i;j<Ngen;j++){
   		if(IdxMother[j]==idx){
   			if(abs(pdgId[j])==startPdg) {idx=j; continue;}
   		if(abs(pdgId[j])==11 && IdxMother[j]==idx) { cat=2; break;} //->ee
   		if(abs(pdgId[j])==13 && IdxMother[j]==idx) { cat=3; break;} //->mumu
   		if(abs(pdgId[j])==15 && IdxMother[j]==idx) { cat=4; break;} //->tautau
   		}
   	}
   if(foundZ) break;
   }
 }  
 if (cat==0) {cout<<"Z classification not found "; recover=true;}
 if(recover) { cout<<" ..trying to recover.. ";
	for(int i=0;i<Ngen;i++){
		if(abs(pdgId[i])==11 && IdxMother[i]<2) { cat=2; break;} //->ee
   		if(abs(pdgId[i])==13 && IdxMother[i]<2) { cat=3; break;} //->mumu
   		if(abs(pdgId[i])==15 && IdxMother[i]<2) { cat=4; break;} //->tautau
		}	
 }
 if(recover) {
 	if(cat==0) cout<<" FAILED!!"<<endl;
 	else cout<<" SUCCESS!! Recovered a "<<cat<<endl;
 }	
 return cat;
}
#endif // Auxiliary_cpp
