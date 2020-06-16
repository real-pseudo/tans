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
  TH1D *z = new TH1D("z","istogramma ",91,-0.18,0.18);//istogramma zrec-zsim
  z->SetDirectory(0);


  for(int i=0 ;i < evsim; i++) {
	  std::vector<double> recz;
    tree->GetEvent(i);
    entries1=hit_det1.ptr->GetEntries();
    entries2=hit_det2.ptr->GetEntries();
    nt_sim->Fill(point.Z,point.mult);
    if(i%1000==0)
      cout<<"Ricostruisco evento "<<i<<endl;
		#if DEBUG
    cout<<"Entries-hit1: "<<entries1<<endl;
    cout<<"Entries-hit2: "<<entries2<<endl;
    cout<<"vtx: "<<point.Z<< "---"<<point.mult <<endl;
		#endif

    double width = 0.1;
    double bin_extreme = 36.4;
    double nbin=(2*bin_extreme)/width;

    TH1D *trackZ = new TH1D("VertrecZ","tracklets",nbin+1,-bin_extreme-(width/2.),bin_extreme+(width/2.));
    trackZ->SetDirectory(0);
    //definisco un vettore di vertici in cui raccogliere i vertici ricostruiti
    Vertex rec_vtx[entries1*entries2];
    int count = 0;

    //ciclo per verificare quali punti rimandano ad un vertice,partendo dagli hit sul più esterno
    //prendo contenuto della cella j-esima e dell'i esima e controllo la differenza in phi per capire se possono appartenere alla stessa traiettoria
    for(int j=0;j<entries2;j++) {
      hit *hit2_event = (hit*) hit_det2.ptr->At(j);
      for(int k=0;k<entries1;k++){
        hit *hit1_event = (hit*) hit_det1.ptr->At(k);
        //cout<<"phi2: "<<hit2_event->getPhi()<<"|||||||||"<<"phi1: "<<hit1_event->getPhi()<<endl;
        deltaphi = abs((hit2_event->getPhi() - hit1_event->getPhi()));
        
        //determino il vertice come intersezione della retta passante per i due hit con l asse del fascio(z)
        if(deltaphi<=DELTAPHI){
        	//cout<<"phi2: "<<hit2_event->getPhi()<<"|||||||||"<<"phi1: "<<hit1_event->getPhi()<<endl;
        	//cout<<deltaphi<<endl;
          reconstruction_vtx(rec_vtx[count],*hit1_event,*hit2_event,det1,deltaR);
          trackZ->Fill(rec_vtx[count].Z);
          recz.push_back(rec_vtx[count].Z);

          count++;
        }
      }
    }
    
    //RICERCA DEL PICCO
   
    int bin_peak = 0;
    double most_prob_Z = 0;
    double average =0;

		//Riempimento istogramma dei tracklets
    #if DEBUG
    for(int i=0;i<count;i++){
      TAxis *xaxis = trackZ->GetXaxis();
      Int_t binx = xaxis->FindBin(rec_vtx[i].Z); 
      cout<<"vtx ricostruito: "<<rec_vtx[i].Z << "pos:"<<binx<< "  centro del bin:" << trackZ->GetBinCenter(binx) <<endl;
      
		}
    #endif
    std::sort(recz.begin(),recz.end());
    

    //Trovo il massimo
    bin_peak = trackZ->GetMaximumBin();
		//Controllo se c'è più di un massimo
		bool twopeaks = more_peaks(trackZ, nbin, bin_peak);
    if(recz.size()>1 && twopeaks==false && (recz.size()>2 || abs(recz[0]-recz[1])<(2*width))){

      double bin_center = trackZ->GetBinCenter(bin_peak);
      //cout << "centro del bin"<< bin_center <<endl;
      //media usando il vector  
      
      double left_value = bin_center - width,
            right_value = bin_center + width; 

      #if DEBUG
      cout << "il mio vettore ha " << recz.size() << " elementi: "; 
      
      for (auto v: recz) {
        cout << v << " ";
      }
      cout << endl; 
  
      cout << "cerco intervallo " << left_value << " " << right_value << endl; 
      #endif
      auto it_left = std::lower_bound(recz.begin(), recz.end(), left_value);
      auto it_right = std::upper_bound(recz.begin(), recz.end(), right_value);


      if (it_right == recz.end())
        it_right = recz.end() - 1; 
      else 
        --it_right; 
        
      int n_elements = it_right - it_left + 1; 
      average = std::accumulate(it_left, it_right + 1, 0.0f) / n_elements; 
      #if DEBUG
      cout << "inizia in " << *it_left << " e finisce in " << *it_right << endl;
      cout << "nell'intervallo ci sono " << n_elements << "elementi" << endl; 
      cout << "valore medio nell'intervallo: " << average << endl; 
      #endif
    
    //#if DEBUG 
    if(i%1000==0)
      cout << "vtx originale:"<< point.Z << " -- vtx ricostruito: "<< average << endl;
    //#endif
    
    //if(point.mult > 1)
    z->Fill(average-point.Z);
    nt_rec.Fill(point.Z, average,average-point.Z,point.mult);
  }

  else{
    c++;
    //cout << "vtx originale:"<< point.Z << " -- vtx non ricostruito??: "<< most_prob_Z <<endl;
  }


  delete trackZ;

  }

 
  
  //z->SetMarkerStyle(21);
  //z->SetMarkerSize(.4);
  z->Draw();

  cout<<"picchi non ric: "<<c<<endl;
  fout.Write();
  fout.Close();






}
