#ifndef UTILITY_H
#define UTILITY_H
#define DELTAS 0.003
#define NOISE 10

#include <Riostream.h>
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TRandom3.h"
#include "TMath.h"
#include "hit.h"
#include "geometry.h"

class ClonesArray {
public:
	TClonesArray* ptr = nullptr; 
	TClonesArray& array;

	ClonesArray(const std::string& classname, unsigned size) : //sto facendo come al solito con il costruttore
		ptr(new TClonesArray(classname.c_str(), size)), //creo un puntatore ad un oggetto TClonesArray
		array(*ptr) {//faccio un alias di ptr
	}

	~ClonesArray() {
		delete ptr; //per deallocare la memoria allocata con new in ptr
	}

	void clear() {
		array.Clear(); //non dovrei far clear su ptr? (è uguale)
	}
}; 

int getMultiplicity();
void smearing(hit& interaction, const Cilindro& c);
void add_noise(hit& interact, const Cilindro& detector, double count_hit, ClonesArray& clonearray);
bool more_peaks(TH1D* histo, const int nbin, const int peak);


#endif