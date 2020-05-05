#include <Riostream.h>
#include "TObject.h"
#include "TMath.h"
#include "Direzione.h"

ClassImp(Direzione)

//________________________________________________________________________
Direzione::Direzione() : TObject(),
 fLabel(0.),
 fPhi(0.),
 fTheta(0.) {
   // default constructor
 }


//___________________________________________________________________________
Direzione::Direzione(int Label,double Phi, double Theta):TObject(),
 fLabel(Label),
 fPhi(Phi),
 fTheta(Theta){
	//standard constructor 
}	     


Direzione::Direzione(const Direzione& source) : TObject(), 
  fLabel(source.fLabel), 
  fPhi(source.fPhi), 
  fTheta(source.fTheta) {

} 

//___________________________________________________________________________
Direzione::~Direzione()	 {
  // destructor
}

