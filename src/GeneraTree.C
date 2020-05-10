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

//prova commento
void GeneraTree() {
	static Vertex point;

	//Dati di partenza 
	double lunghezza = 27;
	hit hits1(0, 0, 0), 
		hits2(0,0, 0);
	Cilindro beampipe(0, 3, 0.08), 
			 riv1(1, 4, 0.02), 
			 riv2(2, 7, 0.02); 
	Particella p;
//	Particella *p = new Particella();

  	// Apertura del file di output
  	TFile hfile("htree.root","RECREATE");
  	// dichiarazione del TTree
  	TTree *tree = new TTree("T","TTree con 4 branches"); //puntatore ad un oggetto ttree
    // se invertissi l'ordine dovrei scrivere
    // tree->SetDirectory(&hfile);

  	TClonesArray *ptrdirez = new TClonesArray("Direzione", ARRAY_SIZE); //vettore di classi
  	TClonesArray& direz = *ptrdirez; //direz è un puntatore ad oggetti della classe Tclones array è una copia di ptrdirez
 	TClonesArray *ptrHIT1 = new TClonesArray("hit", ARRAY_SIZE);
  	TClonesArray& HIT1 = *ptrHIT1;
  	TClonesArray *ptrHIT2 = new TClonesArray("hit", ARRAY_SIZE);
  	TClonesArray& HIT2 = *ptrHIT2;
     
	 
  	// Dichiarazione dei branch del TTree
	tree->Branch("VertMult", &point.X, "X/D:Y:Z:mult/I");
	tree->Branch("Direzione", &ptrdirez);
	tree->Branch("hit1", &ptrHIT1);
	tree->Branch("hit2", &ptrHIT2);
	
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
		int a = 0, b = 0;

		//loop su ogni evento
		for(int j = 0; j < point.mult; j++) {
			Direzione& tst = *(new(direz[j]) Direzione(j, p.phi(), p.theta())); 
			//Direzione *tst = (Direzione*) ptrdirez->At(j); 
			cout << "Direzione # " << tst.getLabel() 
				 << ": phi = " << tst.getPhi()
				 << ", theta = " << tst.getTheta() << endl; 
				 
			//per ogni particella devo controllare la condizione su z
			hits1.intersezione(point, riv1, tst); 
			hits2.intersezione(point, riv2, tst); 
			//hits1.intersezione(point.X, point.Y, point.Z, &riv1, tst);
			//hits2.intersezione(point.X, point.Y, point.Z, &riv2, tst);

			hits1.PrintStatus();
			
			if(hits1.condizione(lunghezza)){
				//cout<<"ok"<<endl;
			//	new(HIT1[a]) hit(hits1.GetX(), hits1.GetY(), hits1.GetZ());

				new(HIT1[a]) hit(hits1); //costruttore di copia 
				a++;
				
			}
			if(hits2.condizione(lunghezza)){
				//cout<<"ok!!!!!!!!"<<endl;
				
			//	new(HIT2[b]) hit(hits2.GetX(), hits2.GetY(), hits2.GetZ());
				new(HIT2[b]) hit(hits2); 
				b++;
				
			}
		}
		
	

		tree->Fill();
		ptrdirez->Clear();
		ptrHIT1->Clear();
        ptrHIT2->Clear();

  }

  // Save all objects in this file
  hfile.Write();

  // Close the file. 
  hfile.Close();


}











