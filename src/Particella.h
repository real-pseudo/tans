#ifndef PARTICELLA_H
#define PARTICELLA_H

#include "TObject.h"

class Particella : public TObject
{

public:

Particella();


virtual ~Particella();
	int molteplicita();
    	double theta();
	double phi();


private:

ClassDef(Particella,1)
};


#endif 


