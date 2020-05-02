#ifndef HIT_H
#define HIT_H

#include "TObject.h"
#include "Cilindro.h"
#include "Direzione.h"
class hit : public TObject
{

public:

hit();
hit(double X, double Y, double Z);

virtual ~hit();
void intersezione(double X,double Y,double Z,Cilindro *c,Direzione *d);
double GetX() const {return fX;} 
double GetY() const {return fY;}
double GetZ() const {return fZ;}
//double GetPhi() const {return fPhi;}
int condizione(double);
void PrintStatus();

private:
    double fX;
    double fY;
    double fZ;
ClassDef(hit,1)
};


#endif 
