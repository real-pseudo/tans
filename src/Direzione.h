#ifndef DIREZIONE_H
#define DIREZIONE_H

#include "TObject.h"

class Direzione : public TObject {
    int fLabel;
    double fPhi;
    double fTheta;
public:

    Direzione();
    Direzione(int Label, double Phi, double Theta);
    Direzione(const Direzione& source); 

    virtual ~Direzione();

    int GetLabel() const {return fLabel;} 
    double GetPhi() const {return fPhi;}
    double GetTheta() const {return fTheta;}

ClassDef(Direzione,1)
};


#endif 


