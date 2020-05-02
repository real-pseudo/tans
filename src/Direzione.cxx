#include <Riostream.h>
#include "TObject.h"
#include "TMath.h"
#include "Direzione.h"

ClassImp(Direzione)

//________________________________________________________________________
Direzione::Direzione():TObject(),
 fLabel(0.),
 fPhi(0.),
 fTheta(0.){
   // default constructor
 }


//___________________________________________________________________________
Direzione::Direzione(int Label,double Phi, double Theta):TObject(),
 fLabel(Label),
 fPhi(Phi),
 fTheta(Theta){
	//standard constructor 
}	     

//___________________________________________________________________________
Direzione::~Direzione()	 {
  // destructor
}

