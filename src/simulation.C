#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "Riostream.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TH1F.h"
#include "TAxis.h"
#include "geometry.h"
#include "particle.h"
#include "hit.h"
#include "Utility.h"
#include <time.h>

#define ARRAY_SIZE 100
#define NUMBER_OF_EVENTS 100


void simulation() {
	clock_t start=clock();
	bool multScattering = true;
	bool noise = true;
	static Vertex point; //vertice che verrà generato
	
	
	//numero di conteggi per ciascun cilindro
	hit		hits1(0, 0, 0), 
			hits2(0, 0, 0),
			hitBP(0, 0, 0);
	//Dati di partenza: dimensioni dei rivelatori
	Cilindro	beampipe(0, 3, 0.08, 27), 
				det1(1, 4, 0.02, 27), 
				det2(2, 7, 0.02, 27); 


		// Apertura del file di output
	TFile hfile("htree.root","RECREATE");

		// dichiarazione del TTree
	std::string n = "4"; //numero branches del TTree
	if(multScattering){
		n = "6"; 
	}

	TTree *tree = new TTree("Tree", ("TTree con " + n + " branches").c_str()); //puntatore ad un oggetto ttree

	ClonesArray particles("Particella", ARRAY_SIZE); 
	ClonesArray hit_det1("hit", ARRAY_SIZE);
	ClonesArray hit_det2("hit", ARRAY_SIZE);
	//Utili solo in caso di Multiple Scattering

	ClonesArray scatter_det1("Particella", ARRAY_SIZE);
	ClonesArray scatter_bp("Particella", ARRAY_SIZE);

	// Dichiarazione dei branch del TTree ____________________________________________________________________
	tree->Branch("VertMult", &point.X, "X/D:Y:Z:mult/I");
	tree->Branch("DirParticella", &particles.ptr);//devo fare .ptr perchè devo accedere al ptr nella classe
	tree->Branch("hit1", &hit_det1.ptr);
	tree->Branch("hit2", &hit_det2.ptr);  

	if(multScattering){ 
		tree->Branch("scattering_bp", &scatter_bp.ptr);
		tree->Branch("scattering_det1", &scatter_det1.ptr);
	}
	

	//INIZIO SIMULAZIONE______________________________________________________________________________________
	for(int i = 0; i < NUMBER_OF_EVENTS; i++) { // loop sugli eventi
		/*Generazione di un vertice casuale:
		dispersione su z di alcuni centimetri, in x e y dell’ordine del decimo di millimetro */
		point.X = gRandom->Gaus(0, 0.01);
		point.Y = gRandom->Gaus(0, 0.01);
		point.Z = gRandom->Gaus(0, 5.3);

		point.mult = getMultiplicity(); //molteplicità dell'evento, generata dalla distribuzione data

		if (i % 1000 == 0)
			cout << "Simulo evento #" << i << " con molteplicità = " << point.mult << endl;

		/*Conteggio delle interazioni con i rivelatori e con la beampipe*/
		int count_hit1 = 0, count_hit2 = 0, count_bp = 0, count_det1 = 0;

		/*Loop per ogni evento*/
		for(int j = 0; j < point.mult; j++) {
			Vertex vtx_hit = point;
			//Metto nell'array particles, la particella j-esima e mi creo una copia di essa in direction
			Particella direction(*new(particles.array[j]) Particella(j));


			/*cout  << "\nDirezione_Particella #" << direction.getLabel() <<
			"(event " << i << ") " <<
			": phi = " << direction.getPhi() <<
			"; theta = " << direction.getTheta() << endl;
		
*/
			
			if(multScattering){
				/*Considero l'interazione con la beampipe solo in caso di Multiple Scattering*/
				hitBP.traject_intersection(point, beampipe, direction);

				/*Se il punto di intersezione su BP soddisfa accettanza non si riempe nulla: 
				non ci sono HIT su nessuno dei due rivelatori*/
				if(hitBP.acceptance(beampipe)){
					direction.scattering();//cambio angoli in seguito a MS
					change_vertex(vtx_hit, hitBP); //prende l'intersezione come nuovo vertice
					//crea particella che è copia di direction e la salva nell'array
					new(scatter_bp.array[count_bp]) Particella(direction);
					#if DEBUG
					cout << "newtheta:" << direction.getTheta() << "\nnewPhi:" <<direction.getPhi() <<endl;
					#endif

					count_bp++;

					/*Interazione con il primo rivelatore*/
					hits1.traject_intersection(vtx_hit, det1, direction);
					
					/*Se il punto di intersezione sul primo rivelatore soddisfa l'accettanza si riempe hit_det1.
					Se non è dentro al primo non può essere nel secondo: man mano che si allontana dal
					vertice è sempre più esterno come punto su asse z*/

					if(hits1.acceptance(det1)){
		//				cout <<"Z1:"<<endl;
		//				hits1.PrintStatus();
						direction.scattering(); //modifica la direzione della particella in seguito al MS
						//cout << "2)newtheta:" << direction.getTheta() << "\nnewPhi:" <<direction.getPhi() <<endl;
						change_vertex(vtx_hit, hits1); //prende l'intersezione come nuovo vertice
						new(scatter_det1.array[count_det1]) Particella(direction);
						count_det1++;

						//Smearing sul primo rivelatore
						smearing(hits1, det1);
						new(hit_det1.array[count_hit1]) hit(hits1) ;
						count_hit1++;
						#if DEBUG
						cout << "newtheta:" << direction.getTheta() << "\nnewPhi:" <<direction.getPhi() <<endl;
						#endif

						hits2.traject_intersection(vtx_hit, det2, direction);
						if(hits2.acceptance(det2)){
		//					cout <<"Z2:" <<endl;
		//					hits2.PrintStatus();
							

							smearing(hits2, det2);
							new(hit_det2.array[count_hit2])hit(hits2);
							count_hit2++;
						}

					}
				}
			}

		/* MS è FALSE si passa direttamente qui dove si effettuare prima l'intersezione 
		sul primo layer, si verifica accettanza,si valuta smearing e si riempe hit_det1.
		Se si è passato il primo rivelatore si passerà poi al secondo e si controlla
		accettanza, smearing e si riempe hit_det2*/

			else {

				hits1.traject_intersection(vtx_hit, det1, direction);
				if(hits1.acceptance(det1)){
		//			cout <<"Z1:"<<endl;
		//			hits1.PrintStatus();
					smearing(hits1, det1);
					new(hit_det1.array[count_hit1]) hit(hits1);
					count_hit1++;

					/*Interazione con rivelatore 2*/
					hits2.traject_intersection(vtx_hit, det2, direction);
					if(hits2.acceptance(det2)){
		//				cout <<"Z2:" <<endl;
		//				hits2.PrintStatus();

						smearing(hits2, det2);
						new(hit_det2.array[count_hit2])hit(hits2);
						count_hit2++;
					}
				}
			}
			/*
			cout << "Z smearing 1: "<<endl;
			hits1.PrintStatus();
			
			cout << "Z smearing 2: "<<endl;
			hits2.PrintStatus();*/
		}

		/*RUMORE: punti aggiuntivi che si mettono per ogni evento simulato oltre ai conteggi veri. 
		Sono conteggi spuri dovuti all'elettronica.Se non dovessi soddisfare le condizioni di accettanza
		del primo detector(senza MS) o della BP(con MS) si salverebbero solo i conteggi di noise che sono
		messi dentro al ciclo degli eventi.*/
		
			if(noise) {
				add_noise(hits1, det1, count_hit1, hit_det1); 
				add_noise(hits2, det2, count_hit2, hit_det2);
			}

//cout<< "evento: "<< i<<endl;
		tree->Fill();

		particles.clear(); 
		hit_det1.clear(); 
		hit_det2.clear();
		if(multScattering) {
			scatter_bp.clear();
			scatter_det1.clear();
		}
	}



	// Save all objects in this file
	hfile.Write();

	// Close the file. 
	hfile.Close();

	clock_t end=clock();
	cout<<"Simulation time: "<<((double)(end-start))<<endl;

}










