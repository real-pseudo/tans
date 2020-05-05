#ifndef CILINDRO_H
#define CILINDRO_H


#include "TObject.h"

class Cilindro : public TObject {
    int fstrato;  //strato
    double fraggio;  //raggio
    double fspessore;   //spessore
    
public: 
 
  Cilindro();//COSTRUTTORE
  Cilindro(int strato, double raggio, double spessore);
  Cilindro(const Cilindro& source); //COSTRUTTORE COPIA

  virtual ~Cilindro();//DISTRUTTORE

  //MyClass& operator=(const MyClass& source);

  //definizione funzioni della classe cilindro
  
  void Setstrato(int);
  void Setraggio(double);
  void Setspessore(double);

  int Getstrato() const;
  double Getraggio() const;
  double Getspessore() const;
  
  void PrintStatus() const; 
  
          
  ClassDef(Cilindro,1) 
};


#endif
