/*
This code produces histogras for the opposite sign dilepton channel.
The tree uploaded can be changed and corrsponfs to different 'Events' files.
The pre-selection criteria ans can changed.
To run upload to root, create an Events object (e.g. type 'Events t'), the loop over the object (type t.Loop() ).
*/
#define Events_cxx
#include "Events.h"   //Change the 'Events' file here which encodes the tree infomation obtained via make class function (i.e. Events->MakeClass()' ) 
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>
void Events::Loop()
{
TH1F *myHisto  = new TH1F("myHisto","Title of Histogram",6,2,8);  //title and size of object 
float NN[101]={};
Long64_t n = 101; 
float x[101] = {};
if (fChain == 0) return;
Long64_t nentries = fChain->GetEntriesFast();
Long64_t nbytes=0,nb=0;
		

float Btag[960869] = {}; //Set up b-tagging array 


Int_t counter = 0;
Long64_t N = 0;
//loop over all events 
for (Long64_t jentry=0; jentry<nentries;jentry++) { 
	Long64_t ientry = LoadTree(jentry);
	//define and initialise parameters for each event 
	Long64_t elec =0;
	Long64_t muon = 0;
	Long64_t u = 0;
        float H_t=0;
	Long64_t p = 0;
        float pt_misselecx = 0;
	float pt_misselecy=0;
	float pt_missmuonx =0;
	float pt_missmuony =0;
	float pt_missx =0;
	float pt_missy =0;
	float pt_missjetx =0;
	float pt_missjety = 0; 
	float pt_miss = 0;
	float x_L = 0.67;   //Discriinators of each working point
	float x_M = 0.93;
	float x_T = 0.988;
	float x_I =0;
	float Lep_pt = 0;
	float First_Discrim = 0;
	float Second_Discrim = 0;
	float Third_Discrim =0;
	Long64_t Loose = 0;
	Long64_t Medium = 0;
	Long64_t Tight = 0;
	Long64_t Inter = 0;
	float Elec_charge[10] = {};
	float Mu_charge[10] ={};
	float Elec_pt[10] = {};
	float Mu_pt[10] = {};
	float Elec_mass[10] = {};
	float Mu_mass[10] = {};	
	float Invariant_mass = 0;
	Long64_t Di_muon = -10;
	Long64_t Di_elec = -10;
	Long64_t second_elec = -10;
	Long64_t second_muon = -10;
	Long64_t Electron_number =0;	
        Long64_t Muon_number =0;
	float Discrim[40] ={};
	Long64_t repeat =0; 


		for (Long64_t Jet =0;Jet<nJet;++Jet) { //loop over jets
			if (Jet_nElectrons[Jet] > 0) {                 //Identification of isolated leptons
				float I_rel = 0;                       //Electron conditions 
                                float pT = 0;   
				for (Long64_t i=0;i<Jet_nElectrons[Jet];++i) {  //loop over each lepton 
						for (Long64_t NumPart = 0;NumPart<nGenPart;++NumPart){
                					float deltaphi = Jet_phi[Jet]-GenPart_phi[NumPart]; //Calculate delta_R to match each particle to jet to jet 
               						if (abs(deltaphi)>= M_PI) deltaphi = 2*M_PI - abs(deltaphi);
                					float deltaeta = Jet_eta[Jet]-GenPart_eta[NumPart];
                					float deltaR = sqrt(deltaphi * deltaphi + deltaeta * deltaeta);
                					if (deltaR<0.4){ //if paticle in that jet, does it correpsond to this electron
                						deltaphi = Electron_phi[Electron_number]-GenPart_phi[NumPart];
								if (abs(deltaphi)>= M_PI) deltaphi = 2*M_PI - abs(deltaphi);
                                        	                deltaeta = Electron_eta[Electron_number]-GenPart_eta[NumPart];
								deltaR = sqrt(deltaphi * deltaphi + deltaeta * deltaeta);
								if (deltaR<0.4){
									pT += GenPart_pt[NumPart]; //if particle and lepton from same vertex take particles transverse momentum 
								}
                					}
       						 }
					I_rel = (pT)/Electron_pt[Electron_number]; //Calculate I_rel
					if (I_rel < 0.15 && Electron_pt[Electron_number] >20){ //Requirments for Electron 
                                        	if (Di_elec ==-10) Di_elec = Electron_number; //records the placments in the array of the first electron of intreset 
                                        	else second_elec = Electron_number;  //records the placments in the array of the second electron of intreset
						//Records infomation of cuurent electron of interest 
                                        	Lep_pt += Electron_pt[Electron_number];
                                        	Elec_charge[Electron_number] =Electron_charge[Electron_number];
                                        	Elec_pt[Electron_number] = Electron_pt[Electron_number];
                                        	Elec_mass[Electron_number] = Electron_mass[Electron_number];
                                        	elec=elec+1;
                                        	pt_misselecx = Electron_pt[Electron_number]*cos(Electron_phi[Electron_number]);
                                        	pt_misselecy = Electron_pt[Electron_number]*sin(Electron_phi[Electron_number]);
					}
					Electron_number += 1; 
				}
			}
		}

		for (Long64_t Jet =0;Jet<nJet;++Jet) {
                        if (Jet_nMuons[Jet] > 0) {                   //Muon conditons same process as for electron but with different conditions 
                                float I_rel = 0;
                                float pT = 0;
                                for (Long64_t i=0;i<Jet_nMuons[Jet];++i) {
                                                for (Long64_t NumPart = 0;NumPart<nGenPart;++NumPart){
                                                        float deltaphi = Jet_phi[Jet]-GenPart_phi[NumPart];
                                                        if (abs(deltaphi)>= M_PI) deltaphi = 2*M_PI - abs(deltaphi);
                                                        float deltaeta = Jet_eta[Jet]-GenPart_eta[NumPart];
                                                        float deltaR = sqrt(deltaphi * deltaphi + deltaeta * deltaeta);
                                                        if (deltaR<0.4){
                                                                deltaphi = Muon_phi[Muon_number]-GenPart_phi[NumPart];
                                                                if (abs(deltaphi)>= M_PI) deltaphi = 2*M_PI - abs(deltaphi);
                                                                deltaeta = Muon_eta[Muon_number]-GenPart_eta[NumPart];
                                                                deltaR = sqrt(deltaphi * deltaphi + deltaeta * deltaeta);
                                                                if (deltaR<0.4){
                                                                        pT += GenPart_pt[NumPart];
                                                                }
                                                        }
                                                 }
                                        I_rel = (pT)/Muon_pt[Muon_number];
                                        if (I_rel < 0.15 && Muon_pt[Muon_number] >20){
                                                if (Di_muon ==-10) Di_muon = Muon_number;
                                                else second_muon = Muon_number;
                                                Lep_pt += Muon_pt[Muon_number];
                                                Mu_charge[Muon_number] =Muon_charge[Muon_number];
                                                Mu_pt[Muon_number] = Muon_pt[Muon_number];
                                                Mu_mass[Muon_number] = Muon_mass[Muon_number];
                                                muon=muon+1;
                                                pt_missmuonx = Muon_pt[Muon_number]*cos(Muon_phi[Muon_number]);
                                                pt_missmuony = Muon_pt[Muon_number]*sin(Muon_phi[Muon_number]);
                                        }
                                        Muon_number += 1;
                                }
                        }
                }


		Long64_t Pos_Lead_Elec = 0;
		Long64_t Neg_Lead_Elec = 0;
		Long64_t Pos_Sub_Elec = 0;                          //Identifying leading and subleading leptons
		Long64_t Neg_Sub_Elec = 0;  			    //Differentiating between, generation and subleading or leading  
		Long64_t Pos_Lead_Muon = 0;
                Long64_t Neg_Lead_Muon = 0;
                Long64_t Pos_Sub_Muon = 0;
                Long64_t Neg_Sub_Muon = 0;
		for (Long64_t k=0;k<10;++k) {
			if (Elec_charge[k] ==1 && Elec_pt[k] >25) Pos_Lead_Elec +=1;
			if (Elec_charge[k] ==1 && Elec_pt[k] > 20 && Elec_pt[k] <=25) Pos_Sub_Elec+=1;
			if (Elec_charge[k] ==-1 && Elec_pt[k] >25) Neg_Lead_Elec +=1;
                        if (Elec_charge[k] ==-1 && Elec_pt[k] > 20 && Elec_pt[k] <=25) Neg_Sub_Elec+=1;
			if (Mu_charge[k] ==1 && Mu_pt[k] >25) Pos_Lead_Muon +=1;
                        if (Mu_charge[k] ==1 && Mu_pt[k] > 20 && Mu_pt[k] <=25) Pos_Sub_Muon+=1;
                        if (Mu_charge[k] ==-1 && Mu_pt[k] >25) Neg_Lead_Muon +=1;
                        if (Mu_charge[k] ==-1 && Mu_pt[k] > 20 && Mu_pt[k] <=25) Neg_Sub_Muon+=1;

					//Calculating invarient masses 
		}     
		if (Di_elec != -10 && second_elec != -10) Invariant_mass = sqrt(2*Electron_pt[Di_elec]*Electron_pt[second_elec]*(cosh(Electron_eta[Di_elec]-Electron_eta[second_elec]) - 
			cos(Electron_phi[Di_elec] - Electron_phi[second_elec])));
		if (Di_elec != -10 && Di_muon != -10) Invariant_mass = sqrt(2*Electron_pt[Di_elec]*Muon_pt[Di_muon]*(cosh(Electron_eta[Di_elec]-Muon_eta[Di_muon]) - cos(Electron_phi[Di_elec]
                        -Muon_phi[Di_muon])));
		if (Di_muon != -10 && second_muon != -10) Invariant_mass = sqrt(2*Muon_pt[Di_muon]*Muon_pt[second_muon]*(cosh(Muon_eta[Di_muon]-Muon_eta[second_muon]) - cos(Muon_phi[Di_muon]
                        -Muon_phi[second_muon])));	
	
		//Different allowed combinations of dileptons 

		if ((Pos_Lead_Elec==1 && Neg_Lead_Muon==1) or (Neg_Lead_Elec==1 && Pos_Lead_Muon==1) or (Pos_Lead_Elec==1 && Neg_Sub_Muon==1) or 
			(Neg_Lead_Elec==1 && Pos_Sub_Muon==1) or (Pos_Sub_Elec==1 && Neg_Lead_Muon ==1) or (Neg_Sub_Elec==1 && Pos_Lead_Muon==1) or (Pos_Lead_Elec==1 && Neg_Lead_Elec==1) 
			or (Pos_Lead_Elec==1 && Neg_Sub_Elec==1) or (Pos_Sub_Elec==1 && Neg_Lead_Elec==1) or (Pos_Lead_Muon==1 && Neg_Lead_Muon==1)   
                        or (Pos_Lead_Muon==1 && Neg_Sub_Muon==1) or (Pos_Sub_Muon==1 && Neg_Lead_Muon==1)){
			if (Pos_Lead_Elec + Neg_Lead_Elec + Pos_Sub_Elec + Neg_Sub_Elec+ Pos_Lead_Muon + Neg_Lead_Muon + Pos_Sub_Muon + Neg_Sub_Muon == 2) { 
					//Calculatign missing transverse momentum 
					for (Long64_t l =0;l<nJet;++l) { 
						pt_missjetx += Jet_pt[l]*cos(Jet_phi[l]);
                                		pt_missjety += Jet_pt[l]*sin(Jet_phi[l]);
					}
				pt_missx = pow(pt_missjetx +pt_misselecx + pt_missmuonx,2);
				pt_missy = pow(pt_missjety +pt_misselecy + pt_missmuony,2);
				pt_miss = sqrt(pt_missx+pt_missy);
				if (nJet>=4){ //requiered number of jets 
					for (Long64_t NumJet= 0;NumJet<nJet;++NumJet) {
						H_t += Jet_pt[NumJet];   //Total Jet Transverse Momentum 
						counter = counter + 1;  //Placement holder for cross-event arrays 
						Btag[counter] = Jet_btagCSVV2[NumJet];  //Change b-tagger here 
						//Btag2[counter] = Jet_btagCMVA[NumJet];
						if (abs(Jet_eta[NumJet])<2.5 && Jet_pt[NumJet] > 30) {
							p=p+1;   //p is the total  number of jets 
							Discrim[NumJet] = Btag[counter];
						}					
					}
					//Distinguishes number of Loose,Medium and tight b-jets
                                	for(Long64_t i = 0;i < 40; ++i)    {			
                                        	if (Discrim[i] > x_L) Loose+=1;
						if (Discrim[i] > x_M) Medium+=1;
						if (Discrim[i] > x_T) Tight+=1;
						if (Discrim[i] > x_I) Inter+=1;
                                	}

					//Find top 3 discriminant values and set them to the last three values in the array
					for(Long64_t i = 0;i < 40; ++i)    {
						if(Discrim[39] < Discrim[i]) Discrim[39] = Discrim[i];
    					}
					repeat = 0;
					for(Long64_t i = 0;i < 40; ++i)    {
						if (Discrim[i] == Discrim[39]) repeat +=1;
						if (repeat ==3) Discrim[38] = Discrim[39];
                                		if(Discrim[38] < Discrim[i] && Discrim[i]<Discrim[39]) Discrim[38] = Discrim[i];
                        		}    
					repeat = 0; 
					if (Discrim[39] == Discrim[38]) repeat = -2;
					for(Long64_t i = 0;i < 40; ++i)    {
						if (Discrim[i] == Discrim[38]) repeat +=1;
                                        	if (repeat ==3) Discrim[37] = Discrim[38];
                                		if(Discrim[37] < Discrim[i] && Discrim[i]<Discrim[38]) Discrim[37] = Discrim[i];
                        		}     
		    	
					First_Discrim = Discrim[39];
					Second_Discrim = Discrim[38];
					Third_Discrim = Discrim[37];

				}
				else {

				for (Long64_t NumJet= 0;NumJet<nJet;++NumJet) {
				counter=counter+1;
				Btag[counter]=-10;
				}
			}
		}
   //Pre-selection criteria 			
if( H_t>500 && p>=4 && pt_miss>50 && Tight>= 0 && Medium>=0 && Loose >= 0) {
N=N+1;
//Prints small breakdown of each events 
cout<<Loose<<"    "<<Medium<<"    "<<Tight<<"     "<<First_Discrim<<"     "<<Second_Discrim<<"     "<<Third_Discrim<<"   "<<nJet<<endl;
myHisto->Fill(Tight);
}
//cout<<"counter: "<<counter<<"    "<<"u: "<<u<<"    "<<"Event Number: "<<"   "<<jentry<<endl;
}
if (ientry < 0) break;
nb=fChain->GetEntry(jentry);nbytes+=nb;
}

//NN[j]=N;


TCanvas *c1 = new TCanvas("c1");
c1->SetLogy();
myHisto->Draw();

//Saves graph into directory, overwrites if there is a graph with the same name 
TFile myGraphFile("Name_of_file.root", "RECREATE");
myGraphFile.cd();
myHisto->Write();
myGraphFile.Close();

}

