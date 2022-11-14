#include<iostream>
#include<TH1.h>
#include<TFile.h>
#include<array>
#include<THStack.h>
#include<vector>
#include<stdio.h>

#define NFILES 6
#define NMCFILES 5

using namespace std;

void Distributions(){
 TFile *FileData = TFile::Open("Data.root");

 TFile *FileDB = TFile::Open("Diboson.root");
 TFile *FileDY = TFile::Open("DYM.root");
 TFile *FileDYm = TFile::Open("DYm.root");
 TFile *FileTT = TFile::Open("TTbkg.root");
 TFile *FileW = TFile::Open("W.root");

 array<TFile*, NFILES> Files={FileData,FileDB,FileDY,FileTT,FileW,FileDYm};

 array<TTree*,NFILES> Trees;
 //std::vector<std::vector<TH1D>> FILEHIST; //a vector (one entry per histo) of vectors (one entry per file)
 std::vector<TH1*> vec;
 TH1D* a;
 std::vector<std::string> names;

      /*a=(TH1D*)key->ReadObj();*/ 

 int numberhistos=0;
   for(auto k : *FileData->GetListOfKeys()) {
      TKey *key = static_cast<TKey*>(k);
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1")) continue;
      numberhistos++;
      std::string temp=(std::string)key->GetName(); 
      names.push_back(temp);
   }
 std::vector<double> Integral;
 cout<<"N histos "<<numberhistos<<endl;
for( int k=0;k<numberhistos;k++){
 cout<<" Histo: "<< names[k] <<endl;
 char nameA[30];
 for(int i=0;i<NFILES;i++){
        cout<<" File: "<< Files[i]->GetName() <<endl;
  	TTree * tout= static_cast<TTree *>(Files[i]->Get("tout"));
	if (k==0) {
		if(i==0){Integral.push_back(tout->GetEntries());}
		else{TH1F *htemp= new TH1F("htemp","htemp",1000,tout->GetMinimum("CorrWeight"),tout->GetMaximum("CorrWeight"));
cout<<"debug"<<endl;
		     tout->Draw("CorrWeight>>htemp","","goff");
 		     double contr=htemp->GetMean()*htemp->GetEntries();
cout<<"debug"<<endl;
		     Integral.push_back(contr);
		     }
		} 
  	string dset;
   	switch(i){
       		case 0: dset="Data";break;
		case 1: dset="DiBoson"; break;
        	case 2: dset="DYM"; break;
		case 3: dset="TTbkg"; break;
		case 4: dset="Wjets"; break;
		case 5: dset="DYm"; break;
  	 }
	cout<<"prestring"<<endl;
  	sprintf(nameA,"%s_%s",names[k].c_str(),dset.c_str());
        cout<<"poststring"<<endl;
   	a = static_cast<TH1D *>(Files[i]->Get(names[k].c_str()));
	a->SetNameTitle(nameA,nameA);
cout<<"Done"<<endl;
	vec.push_back(a);
 	
 }/*
 TCanvas *d1 = new TCanvas(names[k].c_str());
        THStack *hs1 = new THStack("hs1", "");
        vec[1]->SetFillColor(kGreen);
        vec[2]->SetFillColor(kBlack);
        vec[3]->SetFillColor(kRed);
        vec[4]->SetFillColor(kOrange);
	vec[5]->SetFillColor(kBlue);
        vec[1]->SetMarkerColor(kGreen);
        vec[2]->SetMarkerColor(kBlack);
        vec[3]->SetMarkerColor(kRed);
        vec[4]->SetMarkerColor(kOrange);
        hs1->Add(vec[1]); hs1->Add(vec[2]);
        hs1->Add(vec[3]); hs1->Add(vec[4]); hs1->Add(vec[5]);
        hs1->Draw("HIST");
        hs1->GetXaxis()->SetTitle(names[k].c_str());
        vec[0]->SetMarkerColor(kRed);
        vec[0]->SetMarkerSize(2);
        vec[0]->SetMarkerStyle(8);
        vec[0]->Draw("P SAME");
	TRatioPlot* a1 = new TRatioPlot(hs1,vec[0]);
	a1->Draw();
        d1->Update();
        gPad->BuildLegend(0.7,0.7,0.9,0.9);
        gPad->Update();
	d1->SaveAs((names[k]+".pdf").c_str());
        */
	vec.clear();
    
}
 

//use trees->GetEntries()
 cout<<" Watch out: DY correction factors in use, still yield imprecise because no sepration DYm and DYM"<<endl;
 cout<<"Green: Dibosons         "<< Integral[1]<<endl;
 cout<<"Black: DY               "<< Integral[2]*6.077/5.34<<endl;
 cout<<"Orange: ttbkg           "<< Integral[3]<<endl;
 cout<<"Yellow: Wjets           "<< Integral[4]<<endl;
 cout<<"DYm "<<Integral[5]*18610./15880<<endl;

 cout<<"MCsum "<< Integral[1]+Integral[2]*6.077/5.34+Integral[3]+Integral[4]+Integral[5]*18610./15880 <<endl;
 cout<<"Data "<< Integral[0]<<endl;
        


return;
}
