#define Events_cxx
#include "Events.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <cmath>


////TEST TEST TEST TEST

float Error_plus (float TotalJets,float B)
{
Int_t NewB =0;
float n =5;
float m=6;
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);
for (Long64_t inc=1;inc<1000;++inc){
	if ((B+inc)>=TotalJets) break;
        float prob = (B+inc)/TotalJets; 
	std::binomial_distribution<int> distribution (TotalJets,prob);
        float D =0;
	float NumberofB[1000]={}; 
   	for (int i=0; i<1000; ++i){
        	NumberofB[i]= distribution(generator);
         	//cout<<D<<"   ";
		if (NumberofB[i]>B) D=D+1;}
        if (D/1000 >= (n/m)){
                NewB=B+inc;
		break;}          
  
      
        }
return (NewB - B);
}


float Error_minus (float TotalJets,float B)
{
Int_t newB =0;
float n =1;
float m=6;
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);
for (Long64_t inc=1;inc<1000;++inc){
        if ((B-inc)<=0) break;
        float prob = (B-inc)/TotalJets;
        std::binomial_distribution<int> distribution (TotalJets,prob);
        float D =0;
        float NumberofB[1000]={};
        for (int i=0; i<1000; ++i){
                NumberofB[i]= distribution(generator);
               // cout<<D<<"   ";
                if (NumberofB[i]>B) D=D+1;}
        if (D/1000 <= (n/m)){
                newB=B-inc;
                break;}


        }
return newB;
}




void Events::Loop()
{
if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
TCanvas *c1 = new TCanvas("c1","title",200,10,500,300);
   Long64_t nbytes = 0, nb = 0;
Long64_t n = 101; 
float x_CMVA[101] = {};
float x_CSVV2[101] = {};
float x_DeepB[101] = {};
float x_DeepC[101] = {};
float x_DeepFlavB[101] = {};
float y_CMVA[101] = {};
float y_CSVV2[101] = {};
float y_DeepB[101] = {};
float y_DeepC[101] = {};
float y_DeepFlavB[101] = {};
float CMVA[8363] = {};
float CSVV2[8363] = {};
float DeepB[8363] = {};
float DeepC[8363] = {};
float DeepFlavB[8368] = {};
float x_er_CMVA[101] = {};
float y_er_CMVA[101] = {};
float x_er_CSVV2[101] = {};
float y_er_CSVV2[101] = {};
float x_er_DeepB[101] = {};
float y_er_DeepB[101] = {};
float x_er_DeepC[101] = {};
float y_er_DeepC[101] = {};
float x_er_DeepFlavB[101] = {};
float y_er_DeepFlavB[101] = {};
float Discrim[101] = {};
float Id[8363] ={};
for (Long64_t j=0;j<n;++j){
      x_CMVA[j] = 1 - j*0.02;
      x_DeepB[j] = 1 - j*0.01;
      x_CSVV2[j] = 1 - j*0.01;
      Discrim[j] = x_CSVV2[j];
      x_DeepC[j] = 1 - j*0.01;
      x_DeepFlavB[j] = 1 - j*0.01;
}
   Int_t counter = 0;
   Int_t N = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) { 
   Long64_t ientry = LoadTree(jentry);
	for ( Long64_t NumJet= 0;NumJet<nJet;++NumJet) {
	int Btag = 0;
	int Ctag = 0;
	int Ltag = 0;
        counter = counter + 1;
	CMVA[counter] = Jet_btagCMVA[NumJet];
	CSVV2[counter] = Jet_btagCSVV2[NumJet];
        DeepB[counter] = Jet_btagDeepB[NumJet];
	DeepC[counter] = Jet_btagDeepC[NumJet];
	DeepFlavB[counter] = Jet_btagDeepFlavB[NumJet];
	for (Long64_t NumPart = 0;NumPart<nGenPart;++NumPart){
	if ((abs(GenPart_pdgId[NumPart])>4999 && abs(GenPart_pdgId[NumPart]) <6000) or (abs(GenPart_pdgId[NumPart])>499 && 
		abs(GenPart_pdgId[NumPart])<600) or (abs(GenPart_pdgId[NumPart])>10499 && abs(GenPart_pdgId[NumPart])<10600) or (abs(GenPart_pdgId[NumPart])==5) or
		 (abs(GenPart_pdgId[NumPart])>49 && abs(GenPart_pdgId[NumPart])<60) ){
		float deltaphi = Jet_phi[NumJet]-GenPart_phi[NumPart];
		if (abs(deltaphi)>= M_PI) deltaphi = 2*M_PI - abs(deltaphi);
		float deltaeta = Jet_eta[NumJet]-GenPart_eta[NumPart];
		float deltaR = sqrt(deltaphi * deltaphi + deltaeta * deltaeta); 
        	if (deltaR<0.4){
        	Btag = Btag+1;
		}
        }
	else if ((abs(GenPart_pdgId[NumPart])>3999 && abs(GenPart_pdgId[NumPart]) <5000) or (abs(GenPart_pdgId[NumPart])>399 && 
		abs(GenPart_pdgId[NumPart])<500) or (abs(GenPart_pdgId[NumPart])>10399 && abs(GenPart_pdgId[NumPart])<10500) or (abs(GenPart_pdgId[NumPart])==4) or 
		(abs(GenPart_pdgId[NumPart])>39 && abs(GenPart_pdgId[NumPart])<50) ){
		float deltaphi = Jet_phi[NumJet]-GenPart_phi[NumPart];
		if (abs(deltaphi)>= M_PI) deltaphi = 2*M_PI - abs(deltaphi);
                float deltaeta = Jet_eta[NumJet]-GenPart_eta[NumPart];
                float deltaR = sqrt(deltaphi * deltaphi + deltaeta * deltaeta);
                if (deltaR<0.4){
                Ctag = Ctag+1;}
	}
	
	else
                Ltag = Ltag +1;

}
if (Btag>0)
Id[counter]=5;
else if (Btag==0 && Ctag>0)
Id[counter]=4; 
else if (Btag==0 && Ctag==0 && Ltag>0)
Id[counter]=3;
else 
Id[counter]=0;
}
//MyHist->Draw(); 
      // cout << "Event Number" << jentry << endl;
      N =N+ nJet;
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);  nbytes += nb;
      // if (Cut(ientry) < 0) continue;  
  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

for(Int_t xval=0;xval<n;++xval){
        float counter = 0;
        float Btag_Num = 0;
	float Ctag_Num = 0;
	float Ltag_Num= 0;
        for (Int_t yval=0;yval<N;++yval){
           if (DeepB[yval]>x_DeepB[xval]){
            counter=counter+1;
	    if (Id[yval] == 5) Btag_Num = Btag_Num +1;
            else if (Id[yval] == 4) Ctag_Num = Ctag_Num +1;
            else if (Id[yval] == 3) Ltag_Num = Ltag_Num +1;
           }
	   else{ counter=counter;
                 Btag_Num = Btag_Num;
                 Ctag_Num = Ctag_Num;
                 Ltag_Num = Ltag_Num;
                }
            
	}




	//cout<<"BJets: "<<Btag_Num<<"     CJets: "<<Ctag_Num<<"     LJets: "<<Ltag_Num<<"      Discrim: "<<x_DeepB[xval]<<"    Total: "<<counter<<endl;
        y_DeepB[xval] = (Ctag_Num+Ltag_Num);
        x_DeepB[xval]=Btag_Num;
        x_er_DeepB[xval] =Error_plus(counter, Btag_Num);
        y_er_DeepB[xval] =Error_plus(counter, (Ctag_Num + Ltag_Num));
 



}       

for(Int_t xval=0;xval<n;++xval){
        float counter = 0;
        float Btag_Num = 0;
        float Ctag_Num = 0;
        float Ltag_Num= 0;
        for (Int_t yval=0;yval<N;++yval){
           if (CMVA[yval]>x_CMVA[xval]){
            counter=counter+1;
            if (Id[yval] == 5) Btag_Num = Btag_Num +1;
            else if (Id[yval] == 4) Ctag_Num = Ctag_Num +1;
            else if (Id[yval] == 3) Ltag_Num = Ltag_Num +1;
           }  
           else{ counter=counter;
                 Btag_Num = Btag_Num;
                 Ctag_Num = Ctag_Num; 
                 Ltag_Num = Ltag_Num;
                }

        }
	//cout<<"BJets: "<<Btag_Num<<"     CJets: "<<Ctag_Num<<"     LJets: "<<Ltag_Num<<"      Discrim: "<<x_CMVA[xval]<<"    Total: "<<counter<<endl;
        y_CMVA[xval] = (Ctag_Num+Ltag_Num);
        x_CMVA[xval]=Btag_Num;
        x_er_CMVA[xval] = Error_plus(counter, Btag_Num);
        y_er_CMVA[xval] = Error_plus(counter, (Ctag_Num + Ltag_Num));   

  // cout<<"B Jets: "<<Btag_Num<<"   "<<"C Jets: "<<Ctag_Num<<"      "<<"L Jets: "<<Ltag_Num<<"      "<<"Total: "<<counter<<"        CMVA"<<endl;

}


for(Int_t xval=0;xval<n;++xval){
        float counter = 0;
        float Btag_Num = 0;
        float Ctag_Num = 0;
        float Ltag_Num= 0;
        for (Int_t yval=0;yval<N;++yval){
           if (CSVV2[yval]>x_CSVV2[xval]){
            counter=counter+1;
            if (Id[yval] == 5) Btag_Num = Btag_Num +1;
            else if (Id[yval] == 4) Ctag_Num = Ctag_Num +1;
            else if (Id[yval] == 3) Ltag_Num = Ltag_Num +1;
           }  
           else{ counter=counter;
                 Btag_Num = Btag_Num;
                 Ctag_Num = Ctag_Num; 
                 Ltag_Num = Ltag_Num;////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                }

        }
	cout<<"BJets: "<<Btag_Num<<"     CJets: "<<Ctag_Num<<"     LJets: "<<Ltag_Num<<"      Discrim: "<<x_CSVV2[xval]<<"    Total: "<<counter<<endl;
        y_CSVV2[xval] = (Ctag_Num+Ltag_Num);
        x_CSVV2[xval]=Btag_Num;
        x_er_CSVV2[xval] = Error_plus(counter, Btag_Num);
        y_er_CSVV2[xval] = Error_plus(counter, (Ctag_Num + Ltag_Num));     

    //    cout<<"B Jets: "<<Btag_Num<<"	"<<"C Jets: "<<Ctag_Num<<"	"<<"L Jets: "<<Ltag_Num<<"	"<<"Total: "<<counter<<"	CSVV2"<<endl;
}


for(Int_t xval=0;xval<n;++xval){
        float counter = 0;
        float Btag_Num = 0;
        float Ctag_Num = 0;
        float Ltag_Num= 0; 
        for (Int_t yval=0;yval<N;++yval){
           if (DeepC[yval]>x_DeepC[xval]){
            counter=counter+1;
            if (Id[yval] == 5) Btag_Num = Btag_Num +1;
            else if (Id[yval] == 4) Ctag_Num = Ctag_Num +1;
            else if (Id[yval] == 3) Ltag_Num = Ltag_Num +1;
           }
           else{ counter=counter;
                 Btag_Num = Btag_Num;
                 Ctag_Num = Ctag_Num;
                 Ltag_Num = Ltag_Num;
                }
                 
        }        
        //cout<<"btag num: "<<Btag_Num<<" "<<"ctag num: "<<Ctag_Num<<"    "<<"ltag num: "<<Ltag_Num<<endl;
        y_DeepC[xval] = (Ctag_Num+Ltag_Num);
        x_DeepC[xval]=Btag_Num;
        x_er_DeepC[xval] = Error_plus(counter, Btag_Num);
        y_er_DeepC[xval] = Error_plus(counter, (Ctag_Num + Ltag_Num));     

  // cout<<"B Jets: "<<Btag_Num<<"   "<<"C Jets: "<<Ctag_Num<<"      "<<"L Jets: "<<Ltag_Num<<"      "<<"Total: "<<counter<<"        DeepC"<<endl;

}

for(Int_t xval=0;xval<n;++xval){
        float counter = 0;
        float Btag_Num = 0;
        float Ctag_Num = 0;
        float Ltag_Num= 0; 
        for (Int_t yval=0;yval<N;++yval){
           if (DeepFlavB[yval]>x_DeepFlavB[xval]){
            counter=counter+1;
            if (Id[yval] == 5) Btag_Num = Btag_Num +1;
            else if (Id[yval] == 4) Ctag_Num = Ctag_Num +1;
            else if (Id[yval] == 3) Ltag_Num = Ltag_Num +1;
           }
           else{ counter=counter;
                 Btag_Num = Btag_Num;
                 Ctag_Num = Ctag_Num;
                 Ltag_Num = Ltag_Num;
                }
                 
        }        
        //cout<<"btag num: "<<Btag_Num<<" "<<"ctag num: "<<Ctag_Num<<"    "<<"ltag num: "<<Ltag_Num<<endl;
        y_DeepFlavB[xval] = (Ctag_Num+Ltag_Num);
        x_DeepFlavB[xval]=Btag_Num;
        x_er_DeepFlavB[xval] = Error_plus(counter, Btag_Num);
        y_er_DeepFlavB[xval] = Error_plus(counter, (Ctag_Num + Ltag_Num));     

//cout<<"B Jets: "<<Btag_Num<<"   "<<"C Jets: "<<Ctag_Num<<"      "<<"L Jets: "<<Ltag_Num<<"      "<<"Total: "<<counter<<"        DeepFlavB"<<endl;

}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






for (Int_t i=0;i<n;++i){
    x_er_CMVA[i] =x_er_CMVA[i]/x_CMVA[n-1];
    y_er_CMVA[i] =y_er_CMVA[i]/y_CMVA[n-1] ;
    float Norm_x_CMVA= x_CMVA[n-1];
    float Norm_y_CMVA= y_CMVA[n-1];
    x_CMVA[i]=x_CMVA[i]/x_CMVA[n-1];
    y_CMVA[i]=y_CMVA[i]/y_CMVA[n-1];
    //x_er_CMVA[i] = sqrt((x_CMVA[i]*(1-x_CMVA[i]))/Norm_x_CMVA);                             
    //y_er_CMVA[i] = sqrt((y_CMVA[i]*(1-y_CMVA[i]))/Norm_y_CMVA);

     //if (y_CSVV2[i] <= 485){
       // cout<<"Miss Identification rate: "<< y_CSVV2[i]<<"      "<<"Discrim Value: "<<Discrim[i]<<endl;
   // }
//  cout<<"X Error: "<<x_er_CMVA[i]<<"          "<<"Y Error: "<<y_er_CMVA[i]<<endl;}
    x_er_CSVV2[i] =x_er_CSVV2[i]/x_CSVV2[n-1];
    y_er_CSVV2[i] =y_er_CSVV2[i]/y_CSVV2[n-1];
    float Norm_x_CSVV2= x_CSVV2[n-1];
    float Norm_y_CSVV2= y_CSVV2[n-1];
    x_CSVV2[i]=x_CSVV2[i]/x_CSVV2[n-1]; 
    y_CSVV2[i]=y_CSVV2[i]/y_CSVV2[n-1];
    //x_er_CSVV2[i] = sqrt((x_CSVV2[i]*(1-x_CSVV2[i]))/Norm_x_CSVV2);
    //y_er_CSVV2[i] = sqrt((y_CSVV2[i]*(1-y_CSVV2[i]))/Norm_y_CSVV2);    

    x_er_DeepB[i] =x_er_DeepB[i]/x_DeepB[n-1];
    y_er_DeepB[i] =y_er_DeepB[i]/y_DeepB[n-1];
    float Norm_x_DeepB=DeepB[i] = x_DeepB[n-1];
    float Norm_y_DeepB=DeepB[i] = y_DeepB[n-1];
    x_DeepB[i]=x_DeepB[i]/x_DeepB[n-1];
    y_DeepB[i]=y_DeepB[i]/y_DeepB[n-1];
    //x_er_DeepB[i] = sqrt((x_DeepB[i]*(1-x_DeepB[i]))/Norm_x_DeepB);
    //y_er_DeepB[i] = sqrt((y_DeepB[i]*(1-y_DeepB[i]))/Norm_y_DeepB);
//    cout<<"x Error: "<<x_er_DeepB[i]<<"    y Error: "<<y_er_DeepB[i]<<"x Efficiency"<<x_DeepB[i]<<"  y Efficiency"<< y_DeepB[i]<< "    Norm x"<<Norm_x<<"    Norm y"<<Norm_y <<endl;

    x_er_DeepC[i] =x_er_DeepC[i]/x_DeepC[n-1];
    y_er_DeepC[i] =y_er_DeepC[i]/y_DeepC[n-1];
    float Norm_x_DeepC= x_DeepC[n-1];
    float Norm_y_DeepC= y_DeepC[n-1];
    x_DeepC[i]=x_DeepC[i]/x_DeepC[n-1];
    y_DeepC[i]=y_DeepC[i]/y_DeepC[n-1];
    //x_er_DeepC[i] = sqrt((x_DeepC[i]*(1-x_DeepC[i]))/Norm_x_DeepC);
    //y_er_DeepC[i] = sqrt((y_DeepC[i]*(1-y_DeepC[i]))/Norm_y_DeepC);

//    cout<<"x: "<<x_DeepC[i]<<"      y: "<<y_DeepC[i]<<endl;
   
    x_er_DeepFlavB[i] =x_er_DeepFlavB[i]/x_DeepFlavB[n-1];
    y_er_DeepFlavB[i] =y_er_DeepFlavB[i]/y_DeepFlavB[n-1];
    float Norm_x_DeepFlavB= x_DeepFlavB[n-1];
    float Norm_y_DeepFlavB= y_DeepFlavB[n-1];
    x_DeepFlavB[i]=x_DeepFlavB[i]/x_DeepFlavB[n-1];
    y_DeepFlavB[i]=y_DeepFlavB[i]/y_DeepFlavB[n-1];
    //x_er_DeepFlavB[i] = sqrt((x_DeepFlavB[i]*(1-x_DeepFlavB[i]))/Norm_x_DeepFlavB);
    //y_er_DeepFlavB[i] = sqrt((y_DeepFlavB[i]*(1-y_DeepFlavB[i]))/Norm_y_DeepFlavB);

}

for (Long64_t k =0;k<101;++k) {
	float t = y_CMVA[k];
//	cout<<t<<endl;
}

c1->DrawFrame(0,0.001,1,1,"; B_tag Efficiency; Fake rate");
TGraph *gr= new TGraphErrors(n,x_CMVA,y_CMVA,x_er_CMVA,y_er_CMVA);
gr->SetLineColorAlpha(kBlue, 0.6);

TGraph *gr1= new TGraphErrors(n,x_CSVV2,y_CSVV2,x_er_CSVV2,y_er_CSVV2);
gr1->SetLineColorAlpha(kGreen, 0.6);

TGraph *gr2= new TGraphErrors(n,x_DeepB,y_DeepB,x_er_DeepB,y_er_DeepB);
gr2->SetLineColorAlpha(kBlack, 0.6);

TGraph *gr3= new TGraphErrors(n,x_DeepC,y_DeepC,x_er_DeepC,y_er_DeepC);
gr3->SetLineColorAlpha(kMagenta , 0.6);

TGraph *gr4= new TGraphErrors(n,x_DeepFlavB,y_DeepFlavB,x_er_DeepFlavB,y_er_DeepFlavB);
gr4->SetLineColorAlpha(kRed, 0.6);
gPad->SetLogy();
gPad->SetGrid();
//gr->SetLineWidth(2);
//gr1->SetLineWidth(2);
//gr2->SetLineWidth(2);
//gr3->SetLineWidth(2);
//gr4->SetLineWidth(2);

//gr2->SetTitle("Title; x axis; y axis");


TLegend *legend = new TLegend(0.1,0.65,0.25,0.95);
legend->AddEntry(gr,"CMVA","l");
legend->AddEntry(gr1,"CSVV2","l");
legend->AddEntry(gr2,"DeepB","l");
legend->AddEntry(gr3,"DeepC","l");
legend->AddEntry(gr4,"DeepFlavB","l");

gr->Draw("PC;E1");
gr1->Draw("PC;E1");
gr2->Draw("PC;E1");
gr3->Draw("PC;E1");
gr4->Draw("PC;E1");
legend->Draw();


}
