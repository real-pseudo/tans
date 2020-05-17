#ifndef PARTICELLA_H
#define PARTICELLA_H

#include "TObject.h"

class Particella : public TObject {
	int pLabel;
	double phi, theta;

	
	void setTheta();
	void setPhi();

public:
	Particella();
	Particella(int label);
	Particella(const Particella& source);
	void scattering();
	double getTheta() const {return theta;}
	double getPhi() const {return phi;}
	int getLabel() const {return pLabel;}

	virtual ~Particella() {}

   

ClassDef(Particella,1)
};


#endif 


