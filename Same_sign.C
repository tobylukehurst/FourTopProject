#define Events_cxx
#include "Events.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>
void Events::Loop()
{
TH1F *myHisto  = new TH1F("myHisto","4 top Number of Jets, Single muon channel, bjets>=2, H_t>500, p>=7 with pt>30 and |eta|<2.5, pt_miss>50",66,0.67,1);
float NN[101]={};
Long64_t n = 101; 
float x[101] = {};
if (fChain == 0) return;
Long64_t nentries = fChain->GetEntriesFast();
Long64_t nbytes=0,nb=0;
		
float y_DeepFlavB[101] = {};
float x_erDeepFlavB[101]={};
float y_erDeepFlavB[101]={};
float DeepFlavB[960869] = {};
float Btag2[960869] = {};

Int_t counter = 0;
Long64_t N = 0;
for (Long64_t jentry=0; jentry<nentries;jentry++) {
	Long64_t ientry = LoadTree(jentry);
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
	float x_L = 0.67;
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
	float muoneta[8] = {};
	float muonphi[8] = {};
	float eleceta[8] = {};
	float elecphi[8] = {};
	float Invariant_mass = 0;
	Long64_t Di_muon = -10;
	Long64_t Di_elec = -10;
	Long64_t second_elec = -10;
	Long64_t second_muon = -10;
	Long64_t Electron_number =0;	
        Long64_t Muon_number =0;
	float Discrim[40] ={};
	Long64_t repeat =0; 


		for (Long64_t Jet =0;Jet<nJet;++Jet) {
			if (Jet_nElectrons[Jet] > 0) {
				float I_rel = 0;
                                float pT = 0;   
				for (Long64_t i=0;i<Jet_nElectrons[Jet];++i) {
						for (Long64_t NumPart = 0;NumPart<nGenPart;++NumPart){
                					float deltaphi = Jet_phi[Jet]-GenPart_phi[NumPart];
               						if (abs(deltaphi)>= M_PI) deltaphi = 2*M_PI - abs(deltaphi);
                					float deltaeta = Jet_eta[Jet]-GenPart_eta[NumPart];
                					float deltaR = sqrt(deltaphi * deltaphi + deltaeta * deltaeta);
                					if (deltaR<0.4){
                						deltaphi = Electron_phi[Electron_number]-GenPart_phi[NumPart];
								if (abs(deltaphi)>= M_PI) deltaphi = 2*M_PI - abs(deltaphi);
                                        	                deltaeta = Electron_eta[Electron_number]-GenPart_eta[NumPart];
								deltaR = sqrt(deltaphi * deltaphi + deltaeta * deltaeta);
								if (deltaR<0.4){
									pT += GenPart_pt[NumPart];
								}
                					}
       						 }
					I_rel = (pT)/Electron_pt[Electron_number];
					if (I_rel < 0.15 && Electron_pt[Electron_number] >7){
                                        	if (Di_elec ==-10) Di_elec = Electron_number;
                                        	else second_elec = Electron_number;
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
                        if (Jet_nMuons[Jet] > 0) {
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
                                        if (I_rel < 0.15 && Muon_pt[Muon_number] >5){
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
		Long64_t Pos_Sub_Elec = 0;
		Long64_t Neg_Sub_Elec = 0;
		Long64_t Pos_Lead_Muon = 0;
		Long64_t Neg_Lead_Muon = 0;
		Long64_t Pos_Sub_Muon = 0;
		Long64_t Neg_Sub_Muon = 0;
		Long64_t pe = 0; Long64_t ne =0; Long64_t pm =0; Long64_t nm =0;
		float lpm_pt=0; float lpm_eta=0; float lpm_phi=0; float lpm_pt1=0; float lpm_eta1=0; float lpm_phi1=0;
		float spm_pt=0; float spm_eta=0; float spm_phi=0; float spm_pt1=0; float spm_eta1=0; float spm_phi1=0;
		float lnm_pt=0; float lnm_eta=0; float lnm_phi=0; float lnm_pt1=0; float lnm_eta1=0; float lnm_phi1=0;
		float snm_pt=0; float snm_eta=0; float snm_phi=0; float snm_pt1=0; float snm_eta1=0; float snm_phi1=0;
		float lpe_pt=0; float lpe_eta=0; float lpe_phi=0; float lpe_pt1=0; float lpe_eta1=0; float lpe_phi1=0;
		float lne_pt=0; float lne_eta=0; float lne_phi=0; float lne_pt1=0; float lne_eta1=0; float lne_phi1=0;
		float spe_pt=0; float spe_eta=0; float spe_phi=0; float spe_pt1=0; float spe_eta1=0; float spe_phi1=0;
		float sne_pt=0; float sne_eta=0; float sne_phi=0; float sne_pt1=0; float sne_eta1=0; float sne_phi1=0;
		float pe_pt=0; float pe_eta=0; float pe_phi=0; float ne_pt=0; float ne_eta=0; float ne_phi=0;
		float pm_pt=0; float pm_eta=0; float pm_phi=0; float nm_pt=0; float nm_eta=0; float nm_phi=0;
		for(Long64_t k=0;k<8;++k){
			if(Mu_charge[k]==1 && Mu_pt[k]>25){
				Pos_Lead_Muon +=1;
				if(Pos_Lead_Muon ==1){
					lpm_pt = Mu_pt[k];
					lpm_eta = muoneta[k];
					lpm_phi =muonphi[k];
				}
				if(Pos_Lead_Muon==2){
					lpm_pt1 = Mu_pt[k];
					lpm_eta1 = muoneta[k];
					lpm_phi1 = muonphi[k];
				}
			}
			if(Mu_charge[k]==1 && Mu_pt[k]>20 && Mu_pt[k]<=25){
				Pos_Sub_Muon +=1;
				if (Pos_Sub_Muon ==1){
					spm_pt = Mu_pt[k];
					spm_eta = muoneta[k];
					spm_phi =muoneta[k];
				}
				if(Pos_Sub_Muon ==2){
					spm_pt1 = Mu_pt[k];
					spm_eta1 = muoneta[k];
					spm_phi1 =muoneta[k];
				}
			}
			if(Mu_charge[k]==-1 && Mu_pt[k]>25){
				Neg_Lead_Muon +=1;
				if(Neg_Lead_Muon ==1){
					lnm_pt = Mu_pt[k];
					lnm_eta = muoneta[k];
					lnm_phi = muonphi[k];
				}
				if(Neg_Lead_Muon ==2){
					lnm_pt1 = Mu_pt[k];
					lnm_eta1 = muoneta[k];
					lnm_phi1 = muonphi[k];
				}
			}
			if(Mu_charge[k]==-1 && Mu_pt[k]>20 && Mu_pt[k]<=25){
				Neg_Sub_Muon +=1;
				if(Neg_Sub_Muon ==1){
					snm_pt = Mu_pt[k];
					snm_eta= muoneta[k];
					snm_phi= muonphi[k];
				}
				if(Neg_Sub_Muon ==2){
					snm_pt1 = Mu_pt[k];
					snm_eta1= muoneta[k];
					snm_phi1= muonphi[k];
				}
			}
			if(Elec_charge[k]==1 && Elec_pt[k]>25){
				Pos_Lead_Elec +=1;
				if(Pos_Lead_Elec ==1){
					lpe_pt = Elec_pt[k];
					lpe_eta = eleceta[k];
					lpe_phi =elecphi[k];
				}
				if(Pos_Lead_Elec ==2){
					lpe_pt1 = Elec_pt[k];
					lpe_eta1 = eleceta[k];
					lpe_phi1 =elecphi[k];
				}
			}
			if(Elec_charge[k]==1 && Elec_pt[k]>20 && Elec_pt[k]<=25){
				Pos_Sub_Elec +=1;
				if( Pos_Sub_Elec ==1){
					spe_pt = Elec_pt[k];
					spe_eta = eleceta[k];
					spe_phi = elecphi[k];
				}
				if(Pos_Sub_Elec ==2){
					spe_pt1 =Elec_pt[k];
					spe_eta1=eleceta[k];
					spe_phi1=elecphi[k];
				}
			}
			if(Elec_charge[k]==-1 && Elec_pt[k]>=25){
				Neg_Lead_Elec +=1;
				if(Neg_Lead_Elec ==1){
					lne_pt = Elec_pt[k];
					lne_eta = eleceta[k];
					lne_phi = elecphi[k];
				}
				if(Neg_Lead_Elec ==2){
					lne_pt1 = Elec_pt[k];
					lne_eta1 = eleceta[k];
					lne_phi1 = elecphi[k];
				}
			}
			if(Elec_charge[k]==-1 && Elec_pt[k]>20 && Elec_pt[k]<=25){
				Neg_Sub_Elec +=1;
				if(Neg_Sub_Elec ==1){
					sne_pt=Elec_pt[k];
					sne_eta = eleceta[k];
					sne_phi = elecphi[k];
				}
				if(Neg_Sub_Elec ==2){
					sne_pt1 = Elec_pt[k];
					sne_eta1 = eleceta[k];
					sne_phi1 = elecphi[k];
				}
			}
			if( (Elec_charge[k]==1 && Elec_pt[k]<20 && Elec_pt[k]>pe_pt) or ( Elec_charge[k]==1 && Pos_Lead_Elec>2 && Elec_pt[k] > 25) or (Elec_charge[k] ==1 && 
				Pos_Sub_Elec>2 && Elec_pt[k]>20 && Elec_pt[k]<=25)){
				pe +=1;
				if ( pe_pt < Elec_pt[k]) {
					pe_pt = Elec_pt[k];
					pe_eta =eleceta[k];
					pe_phi = elecphi[k];
				}
			}
				if ((Elec_charge[k]==-1 && Elec_pt[k]<20 && Elec_pt[k] >ne_pt) or ( Elec_charge[k]==-1 && Neg_Lead_Elec>2 && Elec_pt[k]>25) or ( Elec_charge[k]==-1 && 
				Neg_Sub_Elec>2 && Elec_pt[k]>20 && Elec_pt[k]<=25)){
				ne +=1;
				if( ne_pt < Elec_pt[k]){
					ne_pt = Elec_pt[k];
					ne_eta = eleceta[k];
					ne_phi = elecphi[k];
				}
			}
			if ( (Mu_charge[k]==1 && Mu_pt[k]<20 && Mu_pt[k]>pm_pt) or (Mu_charge[k]==1 && Pos_Lead_Muon>2 && Mu_pt[k]>25) or ( Mu_charge[k]==1 && Pos_Sub_Muon>2 
				&& Mu_pt[k]>20 && Mu_pt[k]<=25)){
				pm += 1;
				if(pm_pt<Mu_pt[k]){
					pm_pt = Mu_pt[k];
					pm_eta = muoneta[k];
					pm_phi = muonphi[k];
				}
			}
			if ((Mu_charge[k]==-1 && Mu_pt[k]<20 && Mu_pt[k]>pm_pt) or (Mu_charge[k]==-1 && Neg_Lead_Muon>2 && Mu_pt[k]>25) or ( Mu_charge[k]==-1 && 
				Neg_Sub_Muon>2 && Mu_pt[k]>20 && Mu_pt[k]<=25)){
				nm +=1;
				if ( nm_pt < Mu_pt[k]){
					nm_pt = Mu_pt[k];
					nm_eta = muoneta[k];
					nm_phi = muonphi[k]; 
				}
			}
		}

			




		
		if (Di_elec != -10 && second_elec != -10) Invariant_mass = sqrt(2*Electron_pt[Di_elec]*Electron_pt[second_elec]*(cosh(Electron_eta[Di_elec]-Electron_eta[second_elec]) - 
			cos(Electron_phi[Di_elec] - Electron_phi[second_elec])));
		if (Di_elec != -10 && Di_muon != -10) Invariant_mass = sqrt(2*Electron_pt[Di_elec]*Muon_pt[Di_muon]*(cosh(Electron_eta[Di_elec]-Muon_eta[Di_muon]) - cos(Electron_phi[Di_elec]
                        -Muon_phi[Di_muon])));
		if (Di_muon != -10 && second_muon != -10) Invariant_mass = sqrt(2*Muon_pt[Di_muon]*Muon_pt[second_muon]*(cosh(Muon_eta[Di_muon]-Muon_eta[second_muon])-cos(Muon_phi[Di_muon] 
			-Muon_phi[second_muon])));



		float M_lpe_spe = sqrt(2*lpe_pt*spe_pt*(cosh(lpe_eta-spe_eta)-cos(lpe_phi-spe_phi)));
		float M_lne_lne = sqrt(2*lne_pt*lne_pt1*(cosh(lne_eta-lne_eta1)-cos(lne_phi-lne_phi1)));
		float M_lpe_lpe = sqrt(2*lpe_pt*lpe_pt1*(cosh(lpe_eta-lpe_eta1)-cos(lpe_phi-lpe_phi1)));
		float M_lne_sne = sqrt(2*lne_pt*sne_pt*(cosh(lne_eta-sne_eta)-cos(lne_phi-sne_phi)));
		float M_lpm_pm = sqrt(2*lpm_pt*pm_pt*(cosh(lpm_eta-pm_eta)-cos(lpm_phi-pm_phi)));
		float M_lpe_pe = sqrt(2*lpe_pt*pe_pt*(cosh(lpe_eta-pe_eta)-cos(lpe_phi-pe_phi)));
		float M_spe_pe = sqrt(2*spe_pt*pe_pt*(cosh(spe_eta-pe_eta)-cos(spe_phi-pe_phi)));
		float M_spm_pm = sqrt(2*spm_pt*pm_pt*(cosh(spm_eta-pm_eta)-cos(spm_phi-pm_phi)));
		float M_lnm_nm = sqrt(2*lnm_pt*nm_pt*(cosh(lnm_eta-nm_eta)-cos(lnm_phi-nm_phi)));
		float M_lne_ne = sqrt(2*lne_pt*ne_pt*(cosh(lne_eta-ne_eta)-cos(lne_phi-ne_phi)));
		float M_snm_nm = sqrt(2*snm_pt*nm_pt*(cosh(snm_eta-nm_eta)-cos(snm_phi-nm_phi)));
		float M_sne_ne = sqrt(2*sne_pt*ne_pt*(cosh(sne_eta-ne_eta)-cos(sne_phi-ne_phi)));
////////////////
		Long64_t Q = Neg_Lead_Elec+ Pos_Lead_Elec + Pos_Lead_Muon+Neg_Lead_Muon+Neg_Sub_Elec+Pos_Sub_Elec+Pos_Sub_Muon+Neg_Sub_Muon;
			if ((Neg_Lead_Elec==2 && M_lne_lne >12 && Q==2) or
			(Pos_Lead_Elec==2 && M_lpe_lpe>12 && Q==2) or 
			(Pos_Lead_Muon==2 && Q==2) or
			(Neg_Lead_Muon==2 && Q==2) or
			(Pos_Lead_Muon==1 && Pos_Lead_Elec==1 && Q==2 && ( (M_lpm_pm>12 && (M_lpm_pm<76 or M_lpm_pm >106)) or pm ==0 or (M_lpe_pe>12 && (M_lpe_pe<76 or M_lpe_pe >106))
			or pe ==0)) or 
			(Pos_Lead_Muon==1 && Pos_Sub_Elec==1 && Q==2 && ( (M_lpm_pm>12 && (M_lpm_pm<76 or M_lpm_pm >106)) or pm==0 or (M_spe_pe>12 && (M_spe_pe<76 or M_spe_pe >106))
			or pe ==0)) or
			(Pos_Lead_Muon==1 && Pos_Sub_Muon==1 && Q==2 && ( (M_lpm_pm>12 && (M_lpm_pm<76 or M_lpm_pm >106)) or pm==0 or (M_spm_pm>12 && (M_spm_pm<76 or M_spm_pm >106))
			or pm ==0)) or
			(Neg_Lead_Muon==1 && Neg_Lead_Elec==1 && Q==2 && ( (M_lnm_nm>12 && (M_lnm_nm<76 or M_lnm_nm >106)) or nm==0 or (M_lne_ne>12 && (M_lne_ne<76 or M_lne_ne >106))
			or ne ==0)) or
			(Neg_Lead_Muon==1 && Neg_Sub_Elec==1 && Q==2 && ( (M_lnm_nm>12 && (M_lnm_nm<76 or M_lnm_nm >106)) or nm==0 or (M_sne_ne>12 && (M_sne_ne<76 or M_sne_ne >106))
			or ne ==0)) or
			(Neg_Lead_Muon==1 && Neg_Sub_Muon==1 && Q==2 && ( (M_lnm_nm>12 && (M_lnm_nm<76 or M_lnm_nm >106)) or nm==0 or (M_snm_nm>12 && (M_snm_nm<76 or M_snm_nm >106))
			or nm ==0)) or
			(Pos_Lead_Elec==1 && Pos_Sub_Muon==1 && Q==2 && ( (M_lpe_pe>12 && (M_lpe_pe<76 or M_lpe_pe>106)) or pe==0 or (M_spm_pm>12 && (M_spm_pm<76 or M_spm_pm >106))
			or pm ==0)) or
			(Pos_Lead_Elec==1 && Pos_Sub_Elec==1 && M_lpe_spe>12 && Q==2) or
			(Neg_Lead_Elec==1 && Neg_Sub_Muon==1 && Q==2 && ( (M_lne_ne>12 && (M_lne_ne<76 or M_lne_ne >106)) or ne==0 or (M_snm_nm>12 && (M_snm_nm<76 or M_snm_nm >106))
			or nm ==0)) or (Neg_Lead_Elec==1 && Neg_Sub_Elec==1 && M_lne_sne>12 && Q==2)){


					for (Long64_t l =0;l<nJet;++l) { 
						pt_missjetx += Jet_pt[l]*cos(Jet_phi[l]);
                                		pt_missjety += Jet_pt[l]*sin(Jet_phi[l]);
					}
				pt_missx = pow(pt_missjetx +pt_misselecx + pt_missmuonx,2);
				pt_missy = pow(pt_missjety +pt_misselecy + pt_missmuony,2);
				pt_miss = sqrt(pt_missx+pt_missy);
				if (nJet>=4){
					for (Long64_t NumJet= 0;NumJet<nJet;++NumJet) {
						H_t += Jet_pt[NumJet];                                   //////////////////Total Jet Transverse Momentum 
						counter = counter + 1;
						DeepFlavB[counter] = Jet_btagCSVV2[NumJet];
						if (abs(Jet_eta[NumJet])<2.5 && Jet_pt[NumJet] > 30) {
							p=p+1;
							Discrim[NumJet] = DeepFlavB[counter];
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
						//cout<<Discrim[i]<<"    ";
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
				DeepFlavB[counter]=-10;
				}

	
		}
			
if( H_t>300 && p>=2 && pt_miss>50 && Tight>= 0 && Medium>=2 && Loose >= 3) {
N=N+1;
cout<<Loose<<"    "<<Medium<<"    "<<Tight<<"     "<<First_Discrim<<"     "<<Second_Discrim<<"     "<<Third_Discrim<<"   "<<nJet<<endl;
myHisto->Fill(Third_Discrim);
}
}
if (ientry < 0) break;
nb=fChain->GetEntry(jentry);nbytes+=nb;
}

TCanvas *c1 = new TCanvas("c1");
c1->SetLogy();
myHisto->Draw();


TFile myGraphFile("4_SAME_MM+L", "RECREATE");
myGraphFile.cd();
myHisto->Write();
myGraphFile.Close();

}

