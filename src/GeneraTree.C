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
#include "Cilindro.h"
#include "Particella.h"
#include "hit.h"
#include "Utility.h"

#define ARRAY_SIZE 100 
#define NUMBER_OF_EVENTS 100
#define NOISE 10



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
		array.Clear(); //non dovrei far clear su ptr?
	}
}; 


void GeneraTree() {
	bool multScattering = true;
	bool noise = true;
	static Vertex point; //vertice
	
	
	//numero di hit per ciascun cilindro
	hit hits1(0, 0, 0), 
		hits2(0, 0, 0),
		hitBP(0, 0, 0);
	//Dati di partenza: dimensioni dei rivelatori
	double lunghezza = 27;
	Cilindro beampipe(0, 3, 0.08, 27), 
			 det1(1, 4, 0.02, 27), 
			 det2(2, 7, 0.02, 27); 

	//Dati-variabili utile a valutare lo smearing(descrive la risposta del rivelatore)
	double delta_z = 0.00012;//120micrometri in cm
	double deltaphidet1 = 0.003/(det1.getRadius()),deltaphidet2 = 0.003/(det2.getRadius());//0.003 è deltaS
	double z_smear1,z_smear2,phi_smear1,phi_smear2;

	//Variabili utili ad aggiungere il noise (se necessario)
	double z_noise1,z_noise2,phi_noise1,phi_noise2;

  	// Apertura del file di output
  	TFile hfile("htree.root","RECREATE");

  	// dichiarazione del TTree
	std::string n = "4"; //numero branches del TTree
	if(multScattering){
		n = "6"; 
	}

	TTree *tree = new TTree("T", ("TTree con " + n + " branches").c_str()); //puntatore ad un oggetto ttree

	ClonesArray particles("Particella", ARRAY_SIZE); 
	ClonesArray hit_det1("hit", ARRAY_SIZE);
	ClonesArray hit_det2("hit", ARRAY_SIZE);

	//ClonesArray *scatter_det1=nullptr , *scatter_bp = nullptr;//perchè li definisco puntatori?
  //Perchè prima venivano creati i ClonesArray dentro l'if
	
	ClonesArray scatter_det1("Particella", ARRAY_SIZE);
	ClonesArray scatter_bp("Particella", ARRAY_SIZE);

	/*if(multScattering){
		scatter_det1 = new ClonesArray("Particella", ARRAY_SIZE);
		scatter_bp = new ClonesArray("Particella", ARRAY_SIZE);
	}
*/
	 
  // Dichiarazione dei branch del TTree
	tree->Branch("VertMult", &point.X, "X/D:Y:Z:mult/I");
	tree->Branch("Particella", &particles.ptr);//devo fare .ptr perchè devo accedere al ptr nella classe
	tree->Branch("hit1", &hit_det1.ptr);
	tree->Branch("hit2", &hit_det2.ptr);

	if(multScattering){
    tree->Branch("scattering_bp", &scatter_bp.ptr);
    tree->Branch("scattering_det1", &scatter_det1.ptr);
		}


	for(int i = 0; i < NUMBER_OF_EVENTS; i++) { // loop sugli eventi
    //Generazione di un vertice casuale
		point.X = gRandom->Gaus(0, 0.01);
		point.Y = gRandom->Gaus(0, 0.01);
		point.Z = gRandom->Gaus(0, 5.3);

    point.mult = getMultiplicity(); //molteplicità dell'evento, generata dalla distribuzione data


		cout << "Evento #" << i << "----" << 
				"molteplicità = " << point.mult << endl;

		int count_hit1 = 0, count_hit2 = 0, count_bp = 0, count_det1 = 0;

		//loop su ogni evento
		for(int j = 0; j < point.mult; j++) {
		  Vertex vtx_hit = point;
      /*Metto nell'array particles, la particella j-esima e mi creo una copia di essa in direction*/
      Particella direction(*new(particles.array[j]) Particella(j)); 
//			Particella& direction = *(new(particles.array[j]) Particella(j));
			cout  << "\nDirezione # " << direction.getLabel()
			      << ": phi = " << direction.getPhi()
				    << ": theta = " << direction.getTheta() << endl;
			
			if(multScattering){

				hitBP.intersezione(point, beampipe, direction);

				/*Si verifica che il punto di intersezione su BP soddisfi accettanza,
        si cambia direzione e si sostituiscono i valori di vtx con quelli dell'intersezione:
				se non è soddisfatta la condizione, non si riempe nulla: non ci sono HIT 
        nè sul primo layer, nè sul secondo*/

				if(hitBP.accettanza(beampipe)){
          
					direction.scattering();//cambio angoli in seguito a MS
          cout << "newtheta:" << direction.getTheta() << "\nnewPhi:" <<direction.getPhi() <<endl;
          change_vertex(vtx_hit, hitBP); //cambio coordinate cartesiane del vertice
          new(scatter_bp.array[count_bp]) Particella(direction);

					//Particella& changedirect_bp = *(new(scatter_bp.array[count_bp]) Particella(direction));

					/*cout << "\nDirezione # " << changedirect_bp.getLabel() << ": phi = " <<
           changedirect_bp.getPhi() << ": theta = " << changedirect_bp.getTheta() << endl*/
					count_bp++;

					hits1.intersezione(vtx_hit, det1, direction);
				/*Una volta determinato il punto di intersezione sul primo layer, si controlla
          che l'accettanza sia soddisfatta e solo in quel caso si riempe hit_det1
          se non sono dentro al primo non posso essere dentro al secondo, infatti man mano
          che mi allontano dal vertice sono sempre più esterno come punto su asse z*/

					if(hits1.accettanza(det1)){
						cout <<"Z1:"<<endl;
						hits1.PrintStatus();

						direction.scattering();
            change_vertex(vtx_hit, hits1);
            new(scatter_det1.array[count_det1]) Particella(direction);

						//Particella& changedirect_1 = *(new(scatter_det1.array[count_det1]) Particella(changedirect_bp));
						count_det1++;

				//SMEARING:
				/* i valori deltaz e deltaphi vanno a modificare i punti di intersezione con una 
        gaussiana avente media il punto trovato e sigma la delta fornita, in questo modo 
        si in considerazione la risposta reale del rivelatore.*/
		smearing(hits1, det1);

						/*do {
              z_smear1=gRandom->Gaus(hits1.getZ(),delta_z);
            }	while(abs(z_smear1)>=lunghezza/2.);
  					
            phi_smear1=gRandom->Gaus(hits1.getPhi(),deltaphidet1);

						hits1.cartesian(det1,phi_smear1,z_smear1);*/
						new(hit_det1.array[count_hit1]) hit(hits1) ;
						count_hit1++;

						hits2.intersezione(vtx_hit, det2, direction);
            /*Una volta determinato il punto di intersezione sul secondo layer, si controlla
            che l'interazione sia dentro il det2 e in quel caso si riempe hit_det2*/
						if(hits2.accettanza(det2)){
							cout <<"Z2:" <<endl;
							hits2.PrintStatus();

							do{
                			z_smear2=gRandom->Gaus(hits2.getZ(),delta_z);} //SMEARING come prima per il punto trovato
							while(abs(z_smear2)>=det2.getLenght()/2.);
                            
							phi_smear2=gRandom->Gaus(hits2.getPhi(),deltaphidet2);

							hits2.cartesian(det2, phi_smear2,z_smear2);
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

				hits1.intersezione(vtx_hit, det1, direction);

				if(hits1.accettanza(det1)){

					cout <<"Z1:"<<endl;
					hits1.PrintStatus();
				//SMEARING
					z_smear1=gRandom->Gaus(hits1.getZ(),delta_z);
					while(abs(z_smear1)>=det1.getLenght()/2.){
						z_smear1=gRandom->Gaus(hits1.getZ(),delta_z); //se il punto sul rivelatore era dentro e lo z trovato dopo lo smearing non lo è, riestraggo il punto di z_smear
					}
					phi_smear1=gRandom->Gaus(hits1.getPhi(),deltaphidet1);
				//salvataggio punti
					hits1.cartesian(det1,phi_smear1,z_smear1);
					new(hit_det1.array[count_hit1]) hit(hits1);
					count_hit1++;

					hits2.intersezione(vtx_hit, det2, direction);

					if(hits2.accettanza(det2)){

						cout <<"Z2:" <<endl;
						hits2.PrintStatus();

						z_smear2=gRandom->Gaus(hits2.getZ(),delta_z);
						while(abs(z_smear2)>=det2.getLenght()/2.){
							z_smear2=gRandom->Gaus(hits2.getZ(),delta_z);
						}
						phi_smear2=gRandom->Gaus(hits2.getPhi(),deltaphidet2);

						hits2.cartesian(det2, phi_smear2,z_smear2);
						new(hit_det2.array[count_hit2])hit(hits2);
						count_hit2++;
					}
				}
			}

			cout << "Z smearing 1: "<<endl;
			hits1.PrintStatus();
			
			cout << "Z smearing 2: "<<endl;
			hits2.PrintStatus();
		}

		/*RUMORE: punti aggiuntivi che si mettono per ogni evento simulato oltre ai conteggi veri. 
    Sono conteggi spuri dovuti all'elettronica.Se non dovessi soddisfare le condizioni di accettanza
    del primo detector(senza MS) o della BP(con MS) si salverebbero solo i conteggi di noise che sono
    messi dentro al ciclo degli eventi.*/

		//Si estrae z e phi da una distribuzione uniforme,si passa in cartesiane e si riempe HIT come al solito
		  if(noise){
			  for(int k=count_hit1;k<count_hit1+NOISE;k++){
				  z_noise1=lunghezza*(gRandom->Rndm()) - lunghezza/2; //metto - L/2 perchè noi ragioniamo nell'accettanza
				  phi_noise1=2*TMath::Pi()*(gRandom->Rndm());

				  hits1.cartesian(det1,phi_noise1,z_noise1);
				  new(hit_det1.array[k]) hit(hits1);
		 		}

			  for(int i=count_hit2;i<count_hit2+NOISE;i++){
				  z_noise2=lunghezza*(gRandom->Rndm()) - lunghezza/2;
				  phi_noise2=2*(TMath::Pi())*(gRandom->Rndm());

				  hits2.cartesian(det2,phi_noise2,z_noise2);
				  new(hit_det2.array[i]) hit(hits2);
		  		}
		  }


		tree->Fill();

		particles.clear(); 
		hit_det1.clear(); 
		hit_det2.clear();
		if(multScattering) {
				//delete scatter_bp;
				//delete scatter_det1;
			scatter_bp.clear();
			scatter_det1.clear();
			}
	}



	// Save all objects in this file
	hfile.Write();

	// Close the file. 
	hfile.Close();


}










