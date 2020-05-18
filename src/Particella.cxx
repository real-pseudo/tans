#include "Particella.h"
#include <Riostream.h>
#include "TObject.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TRandom3.h"

ClassImp(Particella)

 Particella::Particella() : TObject(),
	pLabel(0.),
 	phi(0.),
 	theta(0.){
   // default constructor
 }

Particella::Particella(int label) : TObject(),
	pLabel(label){
	setTheta();
	setPhi();
}


//Costruttore di copia
Particella::Particella(const Particella& source) : TObject(), 
  pLabel(source.pLabel),
  phi(source.phi),
  theta(source.theta) {
} 

//Genera gli angoli di emissione theta da distribuzione di pseudorapidità
void Particella::setTheta(){
	TFile F("kinem.root"); 
	TH1F *disheta = (TH1F*)F.Get("heta");
  	disheta->SetDirectory(0);
  	disheta->SetMinimum(0);
	
	TH1F *hc = (TH1F*) disheta->Clone(); //check prima avevo messo DrawClone, ma dovrebbe fungere uguale
	hc->GetXaxis()->SetLimits(-1, 1);
	double hetaval = hc->GetRandom();
	theta = 2.*TMath::ATan(TMath::Exp(-hetaval));
	//cout<<"angolo heta: "<<thetaval<<endl;

	F.Close();
	//hc->Close();
	delete disheta;
	//return theta;
}

//Genera angolo phi da una distribuzione uniforme
void Particella::setPhi(){
	phi=2.*TMath::Pi()*(gRandom->Rndm());
	//cout<<"angolo phi: "<<phival<<endl;
	//return phi;
}
/**Genera la nuova direzione in seguito al multiple scattering */
void Particella::scattering(){
    double cd[3] = {0, 0, 0}; //componenti nuova direzione 

    double sphi = TMath::Sin(phi), cphi = TMath::Cos(phi);
    double stheta = TMath::Sin(theta), ctheta = TMath::Cos(theta);

	//angolo di scattering con distribuzione uniforme
    double phi_0 = 2.*TMath::Pi()*(gRandom->Rndm()); 
    double sphi_0 = TMath::Sin(phi_0), cphi_0 = TMath::Cos(phi_0);

	//angolo di scattering con rms 1 mrad
    double theta_0 = gRandom->Gaus(0.,0.001);
    double stheta_0 = TMath::Sin(theta_0), ctheta_0 = TMath::Cos(phi_0);;
    //cout<<"theta0 :  "<<theta_0<<endl;
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
	double r = TMath::Sqrt(cd[0]*cd[0]+cd[1]*cd[1]+cd[2]*cd[2]); //modulo nuova direzione
	//versori nuova direzione
    double ux = cd[0]/r;
    double uy = cd[1]/r;
    double uz = cd[2]/r;

	//nuove direzioni in seguito allo scattering
    theta = TMath::ACos(uz);
    phi = TMath::ATan2(uy,ux);

}

