#include <Riostream.h>
#include "Cilindro.h"
#include "TString.h"

ClassImp(Cilindro)

Cilindro::Cilindro():TObject(),
fstrato(0),
fraggio(0.),
fspessore(0.),
flunghezza(0.)
{
	
}

Cilindro::Cilindro(int strato, double raggio, double spessore, double lunghezza): TObject(),
fstrato(strato), //inizializzazione variabili come int a=5 a(5)
fraggio(raggio),
fspessore(spessore),
flunghezza(lunghezza) {
}

//_________________________________________________________________
Cilindro::Cilindro(const Cilindro &source) : TObject(source),
  fstrato(source.fstrato),
  fraggio(source.fraggio),
  fspessore(source.fspessore),
  flunghezza(source.flunghezza) {
}

//_________________________________________________________________
Cilindro::~Cilindro() {
}

//__SET________________________________________________________________
void Cilindro::setLayer(int strato) {
 	fstrato = strato;
}

void Cilindro::setRadius(double raggio) {
    fraggio = raggio;
}

void Cilindro::setThickness(double spessore) {
    fspessore = spessore;
}

void Cilindro::setLunghezza(double lunghezza) {
    flunghezza = lunghezza;
}

//_________________________________________________________________
int Cilindro::getLayer() const {
    return fstrato;
}

//_________________________________________________________________
double Cilindro::getRadius() const {
    return fraggio;
}

//_________________________________________________________________
double Cilindro::getThickness() const {
    return fspessore;
}

double Cilindro::getLenght() const {
    return flunghezza;
}


//_________________________________________________________________
void Cilindro::printStatus() const {
  // printout
  cout<<"\n Funzione PrintStatus di Cilindro\n";
  cout<<"======================================================\n";
  
  switch(fstrato){

    case 0:
          cout<<"Beam pipe di Berilio di raggio: "<<fraggio<<" cm e di spessore: "<<fspessore<<" mm ."<<endl;
          break;

    case 1:
          cout<<"Layer 1 di raggio: "<<fraggio<< " cm e di spessore: "<<fspessore<<" mm ."<<endl;
          break;

    case 2:
          cout<<"Layer 2 di raggio: "<<fraggio<<" cm e di spessore: "<<fspessore<<" mm ."<<endl;
          break;

    default:
      cout << "hai inserito una scelta non valida" << endl ;
      break;			


  }
}
  // inutile???
/*
double Cilindro::mult_scat(double p) const {
    double theta0 = 0, X_0 = 0;
    switch (fstrato) {
        case 0:
            X_0 = 65.19; //lunghezza di radiazione Be
            break;

        case 1:
            X_0 = 21.82; //lunghezza di radiazione Si          
            break;
    }

    if(X_0) //Calcola theta_0 solo per i primi due strati
        theta0 = ((13.6/p) * (TMath::Sqrt(fspessore/X_0)) * (1 + 0.038*TMath::Log(fspessore/X_0)));

    return theta0; 

} */

//_________________________________________________________________
  /*void Cilindro::SetArray(int *vt, int sz){
    // initialization of the array
    if(fSize>0)delete []fInfo; // clean memory
    fInfo = new int[sz];
    fSize=sz;
    for(int i=0;i<fSize;i++)fInfo[i]=vt[i];
}*/



