#include <Riostream.h>
#include "TObject.h"
#include "TRandom3.h"
#include "TMath.h"
#include "Direzione.h"

ClassImp(Direzione)

//________________________________________________________________________
Direzione::Direzione() : TObject(),
 fLabel(0.),
 fPhi(0.),
 fTheta(0.) {
 }


//___________________________________________________________________________
Direzione::Direzione(int Label,double Phi, double Theta): TObject(),
 fLabel(Label),
 fPhi(Phi),
 fTheta(Theta){
}	     


Direzione::Direzione(const Direzione& source) : TObject(), 
  fLabel(source.fLabel), 
  fPhi(source.fPhi), 
  fTheta(source.fTheta) {

} 

//___________________________________________________________________________
Direzione::~Direzione()	 {
}

void Direzione::rotation(const Particella& particle, const Cilindro& angle){
    double cd[3] = {0, 0, 0}; //componenti nuova direzione 

    double sphi = TMath::Sin(fPhi), cphi = TMath::Cos(fPhi);
    double stheta = TMath::Sin(fTheta), ctheta = TMath::Cos(fTheta);

    double phi_0 = 2.*TMath::Pi()*(gRandom->Rndm());
    double sphi_0 = TMath::Sin(phi_0), cphi_0 = TMath::Cos(phi_0);
    double stheta_0 = 1, ctheta_0 = 1;
   // double stheta_0 = TMath::Sin(angle.mult_scat()), ctheta_0 = TMath::Cos(angle.mult_scat());
    
    //matrice di rotazione
    double m_rot[3][3] = {
        {-sphi, -cphi * ctheta, stheta * cphi},
        {cphi, -sphi * ctheta, stheta * sphi},
        {0., stheta, ctheta}
    };
    double cdp[3] = {stheta_0 * cphi_0, stheta_0 * sphi_0, ctheta_0};
    
    for(int i=0; i<3; ++i){
        for(int j=0; j<3; ++j){
            cd[i] += m_rot[i][j]*cdp[j];
        }
    }

}

