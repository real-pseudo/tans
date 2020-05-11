#include "Utility.h"
#include <Riostream.h>
#include "TFile.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TRandom3.h"
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