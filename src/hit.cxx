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
/*hit::hit(const hit &source): TObject(source),
  fstrato(source.fstrato),
  fraggio(source.fraggio),
  fspessore(source.fspessore)
{
  
}*/

//_________________________________________________________________
hit::~hit() 
{
}

//_________________________________________________________________
void hit::intersezione(double x, double y, double z, Cilindro *c, Direzione *d) {
					double c1,c2,c3;
					c1=sin(d->GetTheta())*cos(d->GetPhi());
					c2=sin(d->GetTheta())*sin(d->GetPhi());
					c3=cos(d->GetTheta());
					double t1,t2,delta;
					delta=(x*c1+y*c2)*(x*c1+y*c2)-(c1*c1+c2*c2)*(x*x+y*y-(c->Getraggio())*(c->Getraggio()));
					t1=(-(x*c1+y*c2)+sqrt(delta))/(c1*c1+c2*c2);
					t2=(-(x*c1+y*c2)-sqrt(delta))/(c1*c1+c2*c2);
					
					if(t1>0){
							fX=x+c1*t1;
							fY=y+c2*t1;
							fZ=z+c3*t1;
					}
					else{
							fX=x+c1*t2;
							fY=y+c2*t2;
							fZ=z+c3*t2;
					}				
}
void hit::PrintStatus(){
		cout<<fZ<<endl;
	
}
int hit::condizione(double lunghezza){
		int flag=0;
		if(abs(fZ)<=lunghezza/2.){
			
			return flag;
		}
		else{
			flag=1;
			return flag;
		}
	
	
	
}
//_________________________________________________________________

