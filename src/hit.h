#ifndef HIT_H
#define HIT_H

#include "TObject.h"
#include "Cilindro.h"
#include "Direzione.h"

// Definizione di una struct (come nuovo tipo grazie a typedef di nome VTX) il mio vertice
typedef struct {
    Double_t X,Y,Z;
    Int_t mult;
} Vertex;


class hit : public TObject {
    double fX;
    double fY;
    double fZ;

    void intersezione(double X, double Y, double Z, const Cilindro *c, const Direzione *d);
public:

    hit();
    hit(double X, double Y, double Z);
    hit(const hit& source); 

    virtual ~hit();
    
    void intersezione(const Vertex& vrx, const Cilindro& c, const Direzione& d); 
    double getX() const {return fX;} 
    double getY() const {return fY;}
    double getZ() const {return fZ;}

    double get_Phi() const {return atan2(fY,fX);}
    bool accettanza(double) const;
    void cylindrical( const Cilindro c,double phi, double z);

    void PrintStatus() const;

    Direzione* rotate() const;

ClassDef(hit,1)
};


#endif 
