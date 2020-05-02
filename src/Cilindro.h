#ifndef CILINDRO_H
#define CILINDRO_H


#include "TObject.h"
///attenzione quando ho // può vuol dire altro oltre il commento


class Cilindro : public TObject {
////////////////////////////////////////////////////////////////////////
// Class just for example
////////////////////////////////////////////////////////////////////////
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


 //
 private:
 //
  int fstrato;  //strato
  double fraggio;  //raggio
  double fspessore;   //spessore  
             


  ClassDef(Cilindro,1) 
};


#endif
