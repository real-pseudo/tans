#ifndef HIT_H
#define HIT_H

#include "TObject.h"
#include "geometry.h"
#include "particle.h"

class hit;

// Definizione di una struct per il vertice
typedef struct {
    Double_t X,Y,Z;
    Int_t mult;
    }
    Vertex;

void change_vertex(Vertex& point, const hit& intersection);
void reconstruction_vtx(Vertex& point, const hit& hit_1, const hit& hit_2, const Cilindro& c, double delta_R);

class hit : public TObject {
    double fX;
    double fY;
    double fZ;
    double fPhi;


public:

    hit();
    hit(double X, double Y, double Z);
    hit(const hit& source); 

    virtual ~hit();
    
    void traject_intersection(const Vertex& vrx, const Cilindro& c, const Particella& particle);
    double getX() const {return fX;} 
    double getY() const {return fY;}
    double getZ() const {return fZ;}
    double getPhi() const;
    //double getPhi() const{return atan2(fY,fX);}
    bool acceptance(const Cilindro& c) const;
    void cartesian( const Cilindro& c,double phi, double z);

    void PrintStatus() const;

   

ClassDef(hit,1)
};


#endif 
