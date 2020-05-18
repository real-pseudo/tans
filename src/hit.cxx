#include <Riostream.h>
#include "hit.h"
#include "TString.h"
#include "TMath.h"


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

void hit::intersezione(const Vertex& vrx, const Cilindro& c, const Particella& particle) {
	//intersezione(vrx.X, vrx.Y, vrx.Z, &c, &direction); 			
	double x = vrx.X, y = vrx.Y, z =vrx.Z;
	double 	c1=sin(particle.getTheta())*cos(particle.getPhi()),
			c2=sin(particle.getTheta())*sin(particle.getPhi()),
			c3=cos(particle.getTheta());	
	double add=x*c1+y*c2, coeff=c1*c1+c2*c2;
	double radDelta=sqrt(add*add-coeff*(x*x+y*y-(c.getRadius())*(c.getRadius())));
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



bool hit::accettanza(const Cilindro& c) const {
	return (abs(fZ)<=c.getLenght()/2); //se l'urto avviene entro la lunghezza del rivelatore restituisce True
}

void hit::cartesian(const Cilindro& c, double phi, double z){
	fX=c.getRadius()*cos(phi);
	fY=c.getRadius()*sin(phi);
	fZ=z;
}

void change_vertex(Vertex& point, const hit& intersection){
	point.X = intersection.getX();
	point.Y = intersection.getY();
	point.Z = intersection.getZ();
}



