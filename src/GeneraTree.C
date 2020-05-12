#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "Direzione.h"
#include "TClonesArray.h"
#include "Riostream.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TAxis.h"
#include "Cilindro.h"
#include "Particella.h"
#include "hit.h"
#include "Utility.h"

#define ARRAY_SIZE 100 
#define NUMBER_OF_EVENTS 100


class ClonesArray {
public:
	TClonesArray* ptr = nullptr; 
	TClonesArray& array;

	ClonesArray(const std::string& classname, unsigned size) :
		ptr(new TClonesArray(classname.c_str(), size)), 
		array(*ptr) {
	}

	~ClonesArray() {
		delete ptr;
	}

	void clear() {
		array.Clear(); 
	}
}; 


void GeneraTree() {
	bool multScattering = true;
	static Vertex point;
	Vertex vtx_hit;

	//Dati di partenza 
	double lunghezza = 27;
	hit hits1(0, 0, 0), 
		hits2(0, 0, 0),
		hitBP(0, 0, 0);
	Cilindro beampipe(0, 3, 0.08), 
			 det1(1, 4, 0.02), 
			 det2(2, 7, 0.02); 

	double smear_z=0.00012;//120micrometri in cm
	double smear_phi=0.003;
	double Z_smear1,Z_smear2,phi_smear1,phi_smear2;



  	// Apertura del file di output
  	TFile hfile("htree.root","RECREATE");
  	// dichiarazione del TTree
	std::string n = "4"; //numero branches del TTree
	if(multScattering) 
		n = "6"; 
	
  	TTree *tree = new TTree("T", ("TTree con " + n + " branches").c_str()); //puntatore ad un oggetto ttree

	ClonesArray particles("Particella", ARRAY_SIZE); 
	ClonesArray hit_det1("hit", ARRAY_SIZE);
	ClonesArray hit_det2("hit", ARRAY_SIZE);

	ClonesArray *scatter_det1 = nullptr, *scatter_bp = nullptr;
	
	if(multScattering){
		scatter_det1 = new ClonesArray("Direzione", ARRAY_SIZE);
		scatter_bp = new ClonesArray("Direzione", ARRAY_SIZE);
	}

	 
  	// Dichiarazione dei branch del TTree
	tree->Branch("VertMult", &point.X, "X/D:Y:Z:mult/I");
	tree->Branch("Particella", &particles.ptr);
	tree->Branch("hit1", &hit_det1.ptr);
	tree->Branch("hit2", &hit_det2.ptr);
	if(multScattering){
		tree->Branch("scattering_bp", &scatter_bp->ptr);
		tree->Branch("scattering_det1", &scatter_det1->ptr);
	}

	for(Int_t i = 0; i < NUMBER_OF_EVENTS; i++) { // loop sugli eventi

		point.mult = getMultiplicity();
		point.X = gRandom->Gaus(0, 0.01);
		point.Y = gRandom->Gaus(0, 0.01);
		point.Z = gRandom->Gaus(0, 5.3);

		

		cout << "Evento #" << i << "----" << 
				"molteplicità = " << point.mult << endl;
		int count_hit1 = 0, count_hit2 = 0;

		//loop su ogni evento
		for(int j = 0; j < point.mult; j++) {
				vtx_hit = point;
				Particella& particle = *(new(particles.array[j]) Particella(j));  
			cout << "\nDirezione # " << particle.getLabel() 
				 << ": phi = " << particle.getPhi()
				 << ": theta = " << particle.getTheta() << endl; 
			
			if(multScattering){
				hitBP.intersezione(vtx_hit, beampipe, particle);
				particle.scattering();// deviazione particella
				vtx_hit.X = hitBP.getX(); 
				vtx_hit.Y = hitBP.getY();
				vtx_hit.Z = hitBP.getZ();
				//MANCA ACCETTANZA E RIEMPIMENTO ISTOGRAMMI
				}
			cout <<"Z beam pipe:"<<vtx_hit.Z <<endl;
			cout << ": phi = " << particle.getPhi()
				 << ": theta = " << particle.getTheta() << endl; 

			hits1.intersezione(vtx_hit, det1, particle);
			if(multScattering){
				particle.scattering();
				vtx_hit.X = hits1.getX();
				vtx_hit.Y = hits1.getY();
				vtx_hit.Z = hits1.getZ();
				//MANCA ACCETTANZA E RIEMPIMENTO ISTOGRAMMI
				}
			cout <<"Z1:"<<endl;
			hits1.PrintStatus();
			cout << ": phi = " << particle.getPhi()
				 << ": theta = " << particle.getTheta() << endl;
				
			hits2.intersezione(vtx_hit, det2, particle); 
			cout <<"Z2:" <<endl;
			hits2.PrintStatus();
			

			if(hits1.accettanza(lunghezza)){
				Z_smear1=gRandom->Gaus(hits1.getZ(),smear_z);
				while(abs(Z_smear1)>=lunghezza/2.){
					Z_smear1=gRandom->Gaus(hits1.getZ(),smear_z);
				}
				phi_smear1=gRandom->Gaus(hits1.getPhi(),smear_phi);

				hits1.cartesian(det1,phi_smear1,Z_smear1);

				new(hit_det1.array[count_hit1]) hit(hits1) ; //costruttore di copia
				count_hit1++;
				
			}
			cout << "Z smearing 1: "<<endl;
			hits1.PrintStatus();
			
			
			
			if(hits2.accettanza(lunghezza)){
				Z_smear2=gRandom->Gaus(hits2.getZ(),smear_z);
				while(abs(Z_smear2)>=lunghezza/2.){
					Z_smear2=gRandom->Gaus(hits2.getZ(),smear_z);
				}
				phi_smear2=gRandom->Gaus(hits2.getPhi(),smear_phi);

				hits2.cartesian(det2, phi_smear2,Z_smear2);

				new(hit_det2.array[count_hit2])hit(hits2);
				count_hit2++;	
			}
			cout << "Z smearing 2: "<<endl;
			hits2.PrintStatus();
		}

		tree->Fill();

		particles.clear(); 
		hit_det1.clear(); 
		hit_det2.clear(); 
	}

	// Save all objects in this file
	hfile.Write();

	// Close the file. 
	hfile.Close();

	if(multScattering) {
		delete scatter_bp;
		delete scatter_det1;
	} 
}











