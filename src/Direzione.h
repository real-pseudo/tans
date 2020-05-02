#ifndef DIREZIONE_H
#define DIREZIONE_H

#include "TObject.h"

class Direzione : public TObject
{

public:

Direzione();
Direzione(int Label, double Phi, double Theta);

virtual ~Direzione();

 int GetLabel() const {return fLabel;} 
 double GetPhi() const {return fPhi;}
 double GetTheta() const {return fTheta;}


private:

int fLabel;
double fPhi;
double fTheta;



ClassDef(Direzione,1)
};


#endif 


