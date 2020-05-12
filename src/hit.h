#ifndef HIT_H
#define HIT_H

#include "TObject.h"
#include "Cilindro.h"
#include "Particella.h"


// Definizione di una struct per il vertice
typedef struct {
    Double_t X,Y,Z;
    Int_t mult;
    }
    Vertex;


class hit : public TObject {
    double fX;
    double fY;
    double fZ;


public:

    hit();
    hit(double X, double Y, double Z);
    hit(const hit& source); 

    virtual ~hit();
    
    void intersezione(const Vertex& vrx, const Cilindro& c, const Particella& particle); 
    double getX() const {return fX;} 
    double getY() const {return fY;}
    double getZ() const {return fZ;}

    double getPhi() const {return atan2(fY,fX);}
    bool accettanza(double) const;
    void cartesian( const Cilindro& c,double phi, double z);

    void PrintStatus() const;

   

ClassDef(hit,1)
};


#endif 
