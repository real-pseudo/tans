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

	//Dati di partenza 
	double lunghezza = 27;
	hit hits1(0, 0, 0), 
		hits2(0,0, 0);
	Cilindro beampipe(0, 3, 0.08), 
			 det1(1, 4, 0.02), 
			 det2(2, 7, 0.02); 
	Particella p;
	
	double smear_z=0.00012;//120micrometri in cm
	double smear_phi=0.003;
	double Z_smear1,Z_smear2,phi_smear1,phi_smear2;
//	Particella *p = new Particella();


  	// Apertura del file di output
  	TFile hfile("htree.root","RECREATE");
  	// dichiarazione del TTree
	std::string n = "4"; //numero branches del TTree
	if(multScattering) 
		n = "6"; 
	
  	TTree *tree = new TTree("T", ("TTree con " + n + " branches").c_str()); //puntatore ad un oggetto ttree

	ClonesArray directions("Direzione", ARRAY_SIZE); 
	ClonesArray hit_det1("hit", ARRAY_SIZE);
	ClonesArray hit_det2("hit", ARRAY_SIZE);

	ClonesArray *scatter_det1 = nullptr, *scatter_bp = nullptr;
	
	if(multScattering){
		scatter_det1 = new ClonesArray("Direzione", ARRAY_SIZE);
		scatter_bp = new ClonesArray("Direzione", ARRAY_SIZE);
	}

	
	
	/*
  	TClonesArray *ptrdirez = new TClonesArray("Direzione", ARRAY_SIZE); //vettore di classi
  	TClonesArray& direz = *ptrdirez; //direz è un puntatore ad oggetti della classe Tclones array è una copia di ptrdirez
 	TClonesArray *ptrHIT1 = new TClonesArray("hit", ARRAY_SIZE);
  	TClonesArray& HIT1 = *ptrHIT1;
  	TClonesArray *ptrHIT2 = new TClonesArray("hit", ARRAY_SIZE);
  	TClonesArray& HIT2 = *ptrHIT2;
	*/


	
     
	 
  	// Dichiarazione dei branch del TTree
	tree->Branch("VertMult", &point.X, "X/D:Y:Z:mult/I");
	tree->Branch("Direzione", &directions.ptr);
	tree->Branch("hit1", &hit_det1.ptr);
	tree->Branch("hit2", &hit_det2.ptr);
	if(multScattering){
		tree->Branch("scattering_bp", &scatter_bp->ptr);
		tree->Branch("scattering_det1", &scatter_det1->ptr);
	}

	for(Int_t i = 0; i < NUMBER_OF_EVENTS; i++) { // loop sugli eventi
	
    // Genero una molteplicita di particelle e un vertice vtx point con coordinate.
    // Int_t numpart=0; //la condizione serve perchè se viene zero non è una moltepicità di particelle
    // while(numpart<=0){
	//    numpart=(Int_t)(dismult->GetRandom());	//molteplicita di particelle generate a partire dalla distribuzione fornita
    //   cout<<i<<"__ "<<numpart<<endl;
    // }

		point.mult = p.molteplicita();
		point.X = gRandom->Gaus(0, 0.01);
		point.Y = gRandom->Gaus(0, 0.01);
		point.Z = gRandom->Gaus(0, 5.3);

		cout << "Evento #" << i << "----" << 
				"molteplicità = " << point.mult << endl;
		int count_hit1 = 0, count_hit2 = 0;

		//loop su ogni evento
		for(int j = 0; j < point.mult; j++) {
			Direzione& tst = *(new(directions.array[j]) Direzione(j, p.phi(), p.theta()));  
			//Direzione& tst = *(new Direzione(j, p.phi(), p.theta()));
			//Direzione *tst = (Direzione*) ptrdirez->At(j); 
			cout << "Direzione # " << tst.getLabel() 
				 << ": phi = " << tst.getPhi()
				 << ", theta = " << tst.getTheta() << endl; 
			
				 
			hits1.intersezione(point, det1, tst); 
			hits2.intersezione(point, det2, tst); 
			//hits1.intersezione(point.X, point.Y, point.Z, &riv1, tst);
			//hits2.intersezione(point.X, point.Y, point.Z, &riv2, tst);

			hits1.PrintStatus();
		
			if(hits1.accettanza(lunghezza)){
				//cout<<"ok"<<endl;
			//	new(HIT1[a]) hit(hits1.GetX(), hits1.GetY(), hits1.GetZ());
				Z_smear1=gRandom->Gaus(hits1.getZ(),smear_z);
				while(abs(Z_smear1)>=lunghezza/2.){
					Z_smear1=gRandom->Gaus(hits1.getZ(),smear_z);
					}
				phi_smear1=gRandom->Gaus(hits1.get_Phi(),smear_phi);

				hits1.cylindrical(det1,phi_smear1,Z_smear1);

				new(hit_det1.array[count_hit1]) hit(hits1) ; //costruttore di copia
				count_hit1++;

				
			}
			if(hits2.accettanza(lunghezza)){
				//cout<<"ok!!!!!!!!"<<endl;
				
			//	new(HIT2[b]) hit(hits2.GetX(), hits2.GetY(), hits2.GetZ());
				Z_smear2=gRandom->Gaus(hits2.getZ(),smear_z);
				while(abs(Z_smear2)>=lunghezza/2.){
						Z_smear2=gRandom->Gaus(hits2.getZ(),smear_z);
						}
				phi_smear2=gRandom->Gaus(hits2.get_Phi(),smear_phi);

							//	new(HIT2[b]) hit(hits2.GetX(), hits2.GetY(), hits2.getZ());
							//new(HIT2[b]) hit(hits2);
				hits2.cylindrical(det2, phi_smear2,Z_smear2);

				new(hit_det2.array[count_hit2])hit(hits2);
				count_hit2++;	
			}
		}

		tree->Fill();

		directions.clear(); 
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











