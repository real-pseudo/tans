#ifndef PARTICELLA_H
#define PARTICELLA_H

#include "TObject.h"
#include "TFile.h"

class Particella : public TObject {
	int pLabel;
	double phi, theta;

	
	void setTheta(TFile *F);
	void setPhi();

public:
	Particella();
	Particella(int label,TFile *P);
	Particella(const Particella& source);

	virtual ~Particella();

	void scattering();
	double getTheta() const {return theta;}
	double getPhi() const {return phi;}
	int getLabel() const {return pLabel;}



   

ClassDef(Particella,1)
};


#endif 


