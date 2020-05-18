#include "Utility.h"
#include <Riostream.h>
#include "TFile.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TRandom3.h"
#include "hit.h"
#include "Cilindro.h"
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

	return interaction.cartesian(c,phi_smear,z_smear);
}
