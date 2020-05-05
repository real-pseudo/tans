#ifndef HIT_H
#define HIT_H

#include "TObject.h"
#include "Cilindro.h"
#include "Direzione.h"


class hit : public TObject {
    double fX;
    double fY;
    double fZ;
public:

    hit();
    hit(double X, double Y, double Z);
    hit(const hit& source); 

    virtual ~hit();
    void intersezione(double X, double Y, double Z, Cilindro *c, Direzione *d);
    double GetX() const {return fX;} 
    double GetY() const {return fY;}
    double GetZ() const {return fZ;}
    //double GetPhi() const {return fPhi;}
    int condizione(double) const;
    void PrintStatus() const;

ClassDef(hit,1)
};


#endif 
