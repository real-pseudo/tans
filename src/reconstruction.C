#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "Riostream.h"
#include "TRandom3.h"
#include "geometry.h"
#include "particle.h"
#include "hit.h"
#include "Utility.h"
#include "TCanvas.h"
#include "TH1D.h"
#include <algorithm>
#include <vector>
#include "TGraph.h"

#define DEBUG 1

#define ARRAY_SIZE 100
#define DELTAPHI 0.006 //6 mrad differenza di phi oltre la quale i due hit non possono appartenere allo stesso vertice

void reconstruction() {
  //Dati di partenza: dimensioni dei rivelatori
  Cilindro beampipe(0, 3, 0.08, 27),
          det1(1, 4, 0.02, 27),
          det2(2, 7, 0.02, 27);
  //Dichiarazione NTupla e file in cui salvare i vertici
  TFile fout("vtxreco.root","recreate");
//  TNtuple *nt = new TNtuple("nt","vertices","zsim:zrec:diff");
	TNtuple nt_rec("nt_rec","vertices","zsim:zrec:diff:mult");
	TNtuple *nt_sim = new TNtuple("z_sim","vertices","zsimtot:multtot");

  //Dichiarazione TClonesArray
  static Vertex point;
  ClonesArray particles("Particella", ARRAY_SIZE);
  ClonesArray hit_det1("hit", ARRAY_SIZE);
  ClonesArray hit_det2("hit", ARRAY_SIZE);

  int c=0;
	int success=0;

  //Apertura file di input
  TFile hfile("htree.root");
  //Lettura TTree  e branch
  TTree *tree = (TTree*) hfile.Get("Tree");

  /*Definizione dei branch e degli indirizzi per la lettura dei dati su ttree
  tree->SetBranchAddress("VertMult",&point.X);
  tree->SetBranchAddress("DirParticella",&particles.ptr);
  tree->SetBranchAddress("hit1",&hit_det1.ptr);
  tree->SetBranchAddress("hit2",&hit_det2.ptr);*/

  TBranch *b1=tree->GetBranch("VertMult");
  TBranch *b2=tree->GetBranch("DirParticella");
  TBranch *b3=tree->GetBranch("hit1");
  TBranch *b4=tree->GetBranch("hit2");
  // Definizione degli indirizzi per la lettura dei dati su ttree
  b1->SetAddress(&point.X);
  b2->SetAddress(&particles.ptr);
  b3->SetAddress(&hit_det1.ptr);
  b4->SetAddress(&hit_det2.ptr);

//varibili utili alla ricostruzione
  int evsim=tree->GetEntries();
  int entries1 = 0, entries2 = 0;
  double deltaR=(det2.getRadius()- det1.getRadius());
  double deltaphi;
  TH1D *z = new TH1D("z","istogramma ",201,-0.2,0.2);//istogramma zrec-zsim
  z->SetDirectory(0);
 

  for(int i=0 ;i < evsim; i++) {
    tree->GetEvent(i);
    entries1=hit_det1.ptr->GetEntries();
    entries2=hit_det2.ptr->GetEntries();
    nt_sim->Fill(point.Z,point.mult);
		
    cout<<"evento "<<i<<endl;
		#if DEBUG
    //cout<<"Entries-hit1: "<<entries1<<endl;
    //cout<<"Entries-hit2: "<<entries2<<endl;
    cout<<"vtx: "<<point.Z<<endl;
		#endif

    //definisco un vettore di vertici in cui raccogliere i vertici ricostruiti
    Vertex rec_vtx[entries1*entries2];
    int count = 0;

    //ciclo per verificare quali punti rimandano ad un vertice,partendo dagli hit sul più esterno
    //prendo contenuto della cella j-esima e dell'i esima e controllo la differenza in phi per capire se possono appartenere alla stessa traiettoria
    for(int j=0;j<entries2;j++) {
      hit *hit2_event = (hit*) hit_det2.ptr->At(j);
      for(int k=0;k<entries1;k++){
        hit *hit1_event = (hit*) hit_det1.ptr->At(k);
        deltaphi = abs((hit2_event->getPhi() - hit1_event->getPhi()));

        //determino il vertice come intersezione della retta passante per i due hit con l asse del fascio(z)
        if(deltaphi<DELTAPHI){
          reconstruction_vtx(rec_vtx[count],*hit1_event,*hit2_event,det1,deltaR);
					#if DEBUG
          cout<<"vtx ricostruito: "<<rec_vtx[count].Z<<endl;
					#endif
          count++;
        }
      }
    }
    //MANCA IL CALCOLO DEL MASSIMO-tracklets DEI VARI VALORI DI Zrec, PER RIEMPIRE L'NTUPLA HO USATO L'ULTIMO VALORE DI Z RICOSTRUITO
    //PRIMA DI PASSARE ALL'EVENTO SUCCESSIVO
    double width = 0.001;
    double bin_extreme = 121;
    double nbin=(2*bin_extreme)/width;
    //TH1D *trackZ = new TH1D("VertrecZ","tracklets",nbin+1,-bin_extreme-(width/2.),bin_extreme+(width/2.));
    TH1D *trackZ = new TH1D("VertrecZ","tracklets",nbin+1,-bin_extreme-(width/2.),bin_extreme+(width/2.));
    trackZ->SetDirectory(0);
    int bin_peak = 0;
    double most_prob_Z = 0;

		//Riempimento istogramma dei tracklets
    for(int i=0;i<count;i++){	
      trackZ->Fill(rec_vtx[i].Z);
      #if DEBUG
      TAxis *xaxis = trackZ->GetXaxis();
      Int_t binx = xaxis->FindBin(rec_vtx[i].Z); 
      cout<<"vtx ricostruito: "<<rec_vtx[i].Z << "pos:"<<binx<< "  centro del bin:" << trackZ->GetBinCenter(binx) <<endl;
      #endif 
		}
		int nbins=0;

		//trackZ->Draw();
		//new TCanvas();

    //Trovo il massimo
    bin_peak = trackZ->GetMaximumBin();
		//Controllo se c'è più di un massimo
		bool twopeaks = more_peaks(trackZ, nbin, bin_peak);
		//Se c'è più di un max fai un rebin massimo 2 volte 
		for(int rebin=1; twopeaks && rebin < 3; rebin++){
			trackZ->Rebin(3);
			//trackZ->Rebin(2);
			nbins=trackZ->GetNbinsX();
			cout << "new nbins: " << nbins<< endl;
    	bin_peak = trackZ->GetMaximumBin();
    	//Verifico nuovamente se il picco è unico
    	twopeaks = more_peaks(trackZ, nbin, bin_peak);
			if(twopeaks==false){
				++success;
			}

		}


    //Calcolo la media intorno al bin con il picco
    if(twopeaks == false) {
    //cout<<"calcolomedia"<<endl;
    	double sum=0;
    	double ncounts=0;
    	for(int i=bin_peak-2; i<bin_peak+3; i++) {
    		sum += trackZ->GetBinContent(i)*trackZ->GetBinCenter(i);
    		ncounts += trackZ->GetBinContent(i);
    	}

    	most_prob_Z = sum/ncounts;

      //if (i % 1000 == 0)
        cout << "vtx originale:"<< point.Z << " -- vtx ricostruito: "<< most_prob_Z <<endl;

        nt_rec.Fill(point.Z, most_prob_Z,most_prob_Z-point.Z,point.mult);
        z->Fill(most_prob_Z-point.Z);
    	}

      else{
    	  c++;
        //cout << "vtx originale:"<< point.Z << " -- vtx non ricostruito??: "<< most_prob_Z <<endl;
      }


   delete trackZ;

  }
  //z->SetMarkerStyle(21);
  //z->SetMarkerSize(.4);
  //z->Draw();

  cout<<"picchi non ric: "<<c<<endl;
	cout<< "rebin utili" << success<< endl;
  fout.Write();
  fout.Close();

/*//grafica
  //leggo ntuple
  TFile filin("vtxreco.root");
  float sim,rec,diff_z,mult_rec,mult_sim,sim_z;
  TNtuple *pippo = (TNtuple*)filin.Get("nt");
      pippo->SetBranchAddress("zsim",&sim);
      pippo->SetBranchAddress("zrec",&rec);
      pippo->SetBranchAddress("diff",&diff_z);
      pippo->SetBranchAddress("mult",&mult_rec);
  TNtuple *nt1 = (TNtuple*)filin.Get("z_sim");
      nt1->SetBranchAddress("zsimtot",&sim_z);
      nt1->SetBranchAddress("multtot",&mult_sim);

  double sigma_z=5.3; //caratteristico del vtx generato
  int nsim=nt1->GetEntries();
  int nrec=pippo->GetEntries();
  int n=12;
  double efficiency[12]={0.} , eff[14]={0};
  double mult_int[12]={0.};
  int multiplicity[12]={4,6,10,15,20,25,30,35,40,45,50,55};//range di molteplicità
  for(int range=0;range<12;range++){
	  double count_sim = 0. , count_rec = 0. , count_sigmatot = 0., count_sigma = 0;
//se il valore di molteplicità ricada nei range stabiliti allora si aumenta il contatore degli eventi simulati
	  for(int i=0;i<nsim;i++){
		  nt1->GetEvent(i);
		  if((mult_sim>multiplicity[range]) && ((mult_sim<=multiplicity[range+1]))){
			  count_sim++;
			  if(abs(sim_z)<sigma_z){
				  count_sigmatot++;
			  }
		  }
	  }
	  //se il valore di molteplicità ricada nei range stabiliti allora si aumenta il contatore degli eventi ricostruiti
	  for(int j=0;j<nrec;j++){
		  pippo->GetEvent(j);
		  if((mult_rec>multiplicity[range]) && ((mult_rec<=multiplicity[range+1]))){
		  	  count_rec++;
		  	if(abs(sim)<sigma_z){
		  					  count_sigma++;
		  				  }

		  }
	  }

	  if(count_sim!=0 && count_rec!=0){
		  cout<<"ev sim= "<<count_sim<<"ev_rec= "<<count_rec<<endl;
		  efficiency[range]=count_rec/count_sim;

	  }

	  if(count_sigmatot!=0 && count_sigma!=0){
		  eff[range]=count_sigma/count_sigmatot;

	  }
	  mult_int[range]=(multiplicity[range]+multiplicity[range+1])/2.;

  }
  new TCanvas("Graph1","Efficienza Vs Molteplicita'",200,10,800,500);
  TGraph *gr1=new TGraph(n-1,mult_int,efficiency);
      gr1->SetTitle("Effficienza vs molteplcita'");
      gr1->GetXaxis()->SetTitle("Molteplicita'");
      gr1->GetYaxis()->SetTitle("Efficienza [cm]");
      gr1->SetMarkerColor(kBlack);
      gr1->SetMarkerStyle(21);

  //    new TCanvas("Graph1","efficienza Vs Molteplicita'",200,10,800,500);

      gr1->Draw("APL");

   new TCanvas("Graph2","efficienza Vs Molteplicita'(Z < 3sigma)",200,10,800,500);
   TGraph *gr2=new TGraph(n-1,mult_int,eff);
   gr2->SetTitle("Efficienza vs molteplcita'(Z < 3sigma)");
   gr2->GetXaxis()->SetTitle("Molteplicita'");
   gr2->GetYaxis()->SetTitle("Efficienza [cm]");
   gr2->SetMarkerStyle(21);
   gr2->Draw("APL");
*/




}
