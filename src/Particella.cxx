#include <Riostream.h>
#include "TObject.h"
#include "TMath.h"
#include "Particella.h"
#include "TFile.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TRandom3.h"





//________________________________________________________________________
Particella::Particella()
 {
   // default constructor
 }

//___________________________________________________________________________
Particella::~Particella()	 {
  // destructor
}

//___________________________________________________________________________
//Definizione funzioni
//funzione che serve a generare la molteplicità delle particelle
int Particella::molteplicita(){
	TFile F("kinem.root");
  	TH1F *dismult = (TH1F*)F.Get("hmul");
  	dismult->SetDirectory(0);
  	dismult->SetMinimum(0);
	int mult=(int)(dismult->GetRandom());
	F.Close();
	delete dismult;
  return mult;
	}

//funzione che genera gli angoli di emissione della particellla
double Particella::theta(){
	TFile F("kinem.root"); 
	TH1F *disheta = (TH1F*)F.Get("heta");
  	disheta->SetDirectory(0);
  	disheta->SetMinimum(0);
	/*
	TH1F *hc = (TH1F*)disheta->DrawClone();
  	cout<<"min i: "<<hc->GetXaxis()->GetXmin()<<endl;
  	cout<<"max i: "<<hc->GetXaxis()->GetXmax()<<endl;
  	hc->GetXaxis()->SetLimits(-1, 1);
  	cout<<"min f: "<<hc->GetXaxis()->GetXmin()<<endl;
  	cout<<"max f: "<<hc->GetXaxis()->GetXmax()<<endl;
	*/
	TH1F *hc = (TH1F*)disheta->DrawClone();
	hc->GetXaxis()->SetLimits(-1, 1);
	double hetaval=hc->GetRandom();
	double thetaval=2.*TMath::ATan(TMath::Exp(-hetaval));
	//cout<<"angolo heta: "<<thetaval<<endl;
	F.Close();
	delete disheta;
  return thetaval;
	}

//funzione che genera angolo phi
double Particella::phi(){
	double phival=2.*TMath::Pi()*(gRandom->Rndm());
	//cout<<"angolo phi: "<<phival<<endl;
  return phival;
	}

