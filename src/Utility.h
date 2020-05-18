#ifndef UTILITY_H
#define UTILITY_H
#define DELTAS 0.003

#include "TObject.h"
#include "hit.h"
#include "Cilindro.h"

int getMultiplicity();
void smearing(hit& interaction, const Cilindro& c);


#endif