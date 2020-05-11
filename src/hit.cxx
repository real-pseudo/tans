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

void hit::intersezione(const Vertex& vrx, const Cilindro& c, const Direzione& d) {
	intersezione(vrx.X, vrx.Y, vrx.Z, &c, &d); 			
}

//_________________________________________________________________
void hit::intersezione(double x, double y, double z, const Cilindro *c, const Direzione *d) {
	double 	c1=sin(d->getTheta())*cos(d->getPhi()),
			c2=sin(d->getTheta())*sin(d->getPhi()),
			c3=cos(d->getTheta());	
	double add=x*c1+y*c2, coeff=c1*c1+c2*c2;
	double radDelta=sqrt(add*add-coeff*(x*x+y*y-(c->getRadius())*(c->getRadius())));
	double t1=(-add+radDelta)/coeff;
	double t2=(-add-radDelta)/coeff;

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
bool hit::accettanza(double lunghezza) const {
	return (abs(fZ)<=lunghezza/2); //se l'urto avviene entro la lunghezza del rivelatore restituisce True

}
//_________________________________________________________________

