#ifndef CILINDRO_H
#define CILINDRO_H


#include "TObject.h"

class Cilindro : public TObject {
    int fstrato;  //strato
    double fraggio;  //raggio
    double fspessore;   //spessore
    double flunghezza;
    
public: 
 
  Cilindro();//COSTRUTTORE
  Cilindro(int strato, double raggio, double spessore, double lunghezza);
  Cilindro(const Cilindro& source); //COSTRUTTORE COPIA

  virtual ~Cilindro();//DISTRUTTORE

  //MyClass& operator=(const MyClass& source);

  //definizione funzioni della classe cilindro
  
  void setLayer(int);
  void setRadius(double);
  void setThickness(double);
  void setLunghezza(double);

  int getLayer() const;
  double getRadius() const;
  double getThickness() const;
  double getLenght() const;

  
  void printStatus() const; 
  //double mult_scat() const;
          
  ClassDef(Cilindro,1) 
};


#endif
