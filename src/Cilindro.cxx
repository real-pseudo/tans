#include <Riostream.h>
#include "Cilindro.h"
#include "TString.h"

ClassImp(Cilindro)
  ///////////////////////////////////////////////////////
  // 
  // This is a trivial example
  // 
  //
  ////////////////////////////////////////////////////////
//_________________________________________________________________
  Cilindro::Cilindro():TObject(),
fstrato(0),
fraggio(0.),
fspessore(0.)
{
	
}

//_________________________________________________________________
  Cilindro::Cilindro(int strato, double raggio, double spessore): TObject(),
fstrato(strato), //inizializzazione variabili come int a=5 a(5)
fraggio(raggio),
fspessore(spessore)
{
	
}

//_________________________________________________________________
Cilindro::Cilindro(const Cilindro &source): TObject(source),
  fstrato(source.fstrato),
  fraggio(source.fraggio),
  fspessore(source.fspessore)
{
  
}

//_________________________________________________________________
Cilindro::~Cilindro() 
{
}

//__SET________________________________________________________________
void Cilindro::Setstrato(int strato)
	{
 		fstrato = strato;
	}
void Cilindro::Setraggio(double raggio)
	{
		fraggio = raggio;
	}
void Cilindro::Setspessore(double spessore)
	{
		fspessore = spessore;
	}


//_________________________________________________________________
int Cilindro::Getstrato() const
	{
		return fstrato;
	}

//_________________________________________________________________
double Cilindro::Getraggio() const {
  
  return fraggio;
}

//_________________________________________________________________
double Cilindro::Getspessore() const {
  
  return fspessore;
}

//_________________________________________________________________
void Cilindro::PrintStatus() const {
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
          printf("hai inserito una scelta non valida \n");
          break;			


}
		


  }

//_________________________________________________________________
  /*void Cilindro::SetArray(int *vt, int sz){
    // initialization of the array
    if(fSize>0)delete []fInfo; // clean memory
    fInfo = new int[sz];
    fSize=sz;
    for(int i=0;i<fSize;i++)fInfo[i]=vt[i];
}*/



