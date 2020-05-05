#include <Riostream.h>
#include "hit.h"
#include "TString.h"
#include "TMath.h"
#include "Cilindro.h"
#include "Direzione.h"
ClassImp(hit)
//_________________________________________________________________
hit::hit(): TObject(),
	fX(0.),
 	fY(0.),
 	fZ(0.){
   // default constructor
 }

//_________________________________________________________________
hit::hit(double X, double Y, double Z): TObject(),
 fX(X),
 fY(Y),
 fZ(Z){
	//standard constructor 
}

//_________________________________________________________________
hit::hit(const hit &source): TObject(source),
  fX(source.fX),
  fY(source.fY),
  fZ(source.fZ)
{
  
}

//_________________________________________________________________
hit::~hit() 
{
}

//_________________________________________________________________
void hit::intersezione(double x, double y, double z, Cilindro *c, Direzione *d) {
	double 	c1=sin(d->GetTheta())*cos(d->GetPhi()),
			c2=sin(d->GetTheta())*sin(d->GetPhi()),
			c3=cos(d->GetTheta());	
	double add1=x*c1+y*c2, coeff=c1*c1+c2*c2;
	double radDelta=sqrt(add1*add1-coeff*(x*x+y*y-(c->Getraggio())*(c->Getraggio())));
	double t1=(-add1+radDelta)/coeff;
	double t2=(-add1-radDelta)/coeff;

	double t;
	if (t1>0) 
		t=t1; 
	else 	
		t=t2;

	fX=x+c1*t;
	fY=y+c2*t;
	fZ=z+c3*t;			
}
void hit::PrintStatus() const {
		cout<<fZ<<endl;
	
}
int hit::condizione(double lunghezza) const {
	return (abs(fZ)<=lunghezza/2); //se la condizione è vera restituisce true 
}
//_________________________________________________________________

