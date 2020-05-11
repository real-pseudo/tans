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

    int getLabel() const {return fLabel;} 
    double getPhi() const {return fPhi;}
    double getTheta() const {return fTheta;}
    void rotation(const Particella& particle, const Cilindro& layer);

ClassDef(Direzione,1)
};


#endif 


