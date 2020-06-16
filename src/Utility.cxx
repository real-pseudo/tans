#include "Utility.h"

#define DEBUG 0

//___________________________________________________________________________
//Genera la molteplicità delle particelle da una distribuzione data
int getMultiplicity(){
	TFile F("kinem.root");
  	TH1F *dismult = (TH1F*) F.Get("hmul");
  	dismult->SetDirectory(0);
  	dismult->SetMinimum(0);
	int mult = (int) dismult->GetRandom();
	F.Close();
	delete dismult;
	return mult;
}

//___________________________________________________________________________
//Smearing delle coordinate sui rivelatori
/* 	i valori deltaz e deltaphi vanno a modificare i punti di intersezione con una 
	gaussiana avente media il punto trovato e sigma la delta fornita, in questo modo 
    si in considerazione la risposta reale del rivelatore.*/

void smearing(hit& interaction, const Cilindro& c){
	//Dati-variabili utile a valutare lo smearing(descrive la risposta del rivelatore)
	double delta_z = 0.00012;//120micrometri in cm
	double deltaphidet = DELTAS/(c.getRadius());//0.003 è deltaS

	double z_smear = 0;
	do {
        z_smear=gRandom->Gaus(interaction.getZ(),delta_z);
        }	
	while(abs(z_smear)>=c.getLenght()/2.);
  					
    double phi_smear=gRandom->Gaus(interaction.getPhi(),deltaphidet);

	interaction.cartesian(c,phi_smear,z_smear);
}

//Si estrae z e phi da una distribuzione uniforme,si passa in cartesiane e si riempe HIT come al solito
void add_noise(hit& interact, const Cilindro& detector, double count_hit, ClonesArray& clonearray){
	for(int k=count_hit; k<count_hit+NOISE; k++){
		double z_noise=detector.getLenght()*(gRandom->Rndm()) - detector.getLenght()/2; //metto - L/2 perchè noi ragioniamo nell'accettanza
		double phi_noise=2*TMath::Pi()*(gRandom->Rndm());
		interact.cartesian(detector,phi_noise,z_noise);
		new(clonearray.array[k]) hit(interact);
		}		
	}

bool more_peaks(TH1D* histo, const int nbin, const int peak){

	for(int v=0;v<nbin+1;v++){
		if(v!=peak && histo->GetBinContent(v) == histo->GetBinContent(peak)&& abs(v-peak)>1){
			return true;
		}
	}
	return false;
}



 
