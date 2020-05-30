#include <Riostream.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "Riostream.h"
#include "TRandom3.h"
#include "Cilindro.h"
#include "Particella.h"
#include "hit.h"
#include "Utility.h"
#include "TCanvas.h"
#include "TH1D.h"
#include <algorithm>
#include <vector>

#define ARRAY_SIZE 100
#define DELTAPHI 0.006 //6 mrad differenza di phi oltre la quale i due hit non possono appartenere allo stesso vertice

void LeggiTree() {
  //Dati di partenza: dimensioni dei rivelatori
  Cilindro beampipe(0, 3, 0.08, 27),
          det1(1, 4, 0.02, 27),
          det2(2, 7, 0.02, 27);
  //Dichiarazione NTupla e file in cui salvare i vertici
  TFile fout("vtxreco.root","recreate");
  TNtuple *nt = new TNtuple("nt","vertices","zsim:zrec:diff");

  //Dichiarazione TClonesArray
  static Vertex point;
  ClonesArray particles("Particella", ARRAY_SIZE);
  ClonesArray hit_det1("hit", ARRAY_SIZE);
  ClonesArray hit_det2("hit", ARRAY_SIZE);


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
  //TH1D *z = new TH1D("z","istogramma ",15,-10.5,10.5);
  //z->SetDirectory(0);
  /*****
  double width = 0.1;
  double bin_extreme = 200;//det2.getRadius()*2*100;
  double nbin=(2*bin_extreme)/width;
//  TH1D *trackZ = new TH1D("VertrecZ","tracklets",nbin+5,-bin_extreme-2-(width/2.),bin_extreme+2+(width/2.));
  TH1D *trackZ = new TH1D("VertrecZ","tracklets",nbin+5,-bin_extreme-(width/2.),bin_extreme+(width/2.)); ***/

  cout << "porcoddio" << endl; 

  for(int i=0 ;i < evsim; i++) {
    tree->GetEvent(i);
    entries1=hit_det1.ptr->GetEntries();
    entries2=hit_det2.ptr->GetEntries();
    //cout<<"evento "<<i<<endl;
    //cout<<"Entries-hit1: "<<entries1<<endl;
    //cout<<"Entries-hit2: "<<entries2<<endl;
    //cout<<"vtx: "<<point.Z<<endl;

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
          //cout<<"vtx ricostruito: "<<rec_vtx[count].Z<<endl;
          count++;
        }
      }
    }
    //MANCA IL CALCOLO DEL MASSIMO-tracklets DEI VARI VALORI DI Zrec, PER RIEMPIRE L'NTUPLA HO USATO L'ULTIMO VALORE DI Z RICOSTRUITO
    //PRIMA DI PASSARE ALL'EVENTO SUCCESSIVO
    double width = 0.1;
    double bin_extreme = 200;//det2.getRadius()*2*100;
    double nbin=(2*bin_extreme)/width;
    TH1D *trackZ = new TH1D("VertrecZ","tracklets",nbin+5,-bin_extreme-2-(width/2.),bin_extreme+2+(width/2.));
    int bin_peak = 0, nc = 0;
    double most_prob_Z = 0;
    bool twopeaks = false; //variabile di controllo per la presenza di più picchi

    
    for(int i=0;i<count;i++){	
      trackZ->Fill(rec_vtx[i].Z);
      TAxis *xaxis = trackZ->GetXaxis();
      Int_t binx = xaxis->FindBin(rec_vtx[i].Z);
      //cout<<"vtx ricostruito: "<<rec_vtx[i].Z << "pos:"<<binx<<endl;
    }
    //Trovo il massimo
    bin_peak = trackZ->GetMaximumBin();
    
    //Verifico che il massimo sia unico
    for(int v=0;v<nbin+1;v++){
      if(trackZ->GetBinContent(v) == trackZ->GetBinContent(bin_peak) && v!=bin_peak && abs(v-bin_peak) >= 10) {
        twopeaks = true; 
        break;
      }
    }
    
    //Ricerco i limiti per il rebin (da fare solo se twopeaks=true)
    if(twopeaks){
      trackZ->Rebin(5);
      twopeaks = false;
      bin_peak = trackZ->GetMaximumBin();
      //Verifico nuovamente se il picco è unico
      for(int v=0;v<nbin+1;v++){
        if(trackZ->GetBinContent(v) == trackZ->GetBinContent(bin_peak) && v!=bin_peak && abs(v-bin_peak) >= 10) 
          twopeaks = true; 
          break;
      }    
    
    //Calcolo la media intorno al bin con il picco
    if(twopeaks == false) {
      double sum=0;
      double ncounts=0;
      for(int i=bin_peak-3; i<bin_peak+4; i++) {
        sum += trackZ->GetBinContent(i)*trackZ->GetBinCenter(i);
        ncounts += trackZ->GetBinContent(i);
      }
      most_prob_Z = sum/ncounts;

      if (i % 1000 == 0)
        cout << "vtx originale:"<< point.Z << " -- vtx ricostruito: "<< most_prob_Z <<endl;
    
      //nt->Fill(point.Z, most_prob_Z,  (point.Z - most_prob_Z));
    //QUESTO NON VA BENE PERCHÈ STO SCEGLIENDO
    if(abs(point.Z - most_prob_Z)<1)
        nt->Fill(point.Z, most_prob_Z,  (point.Z - most_prob_Z)*100);
      else{
        cout << "vtx originale:"<< point.Z << " -- vtx ricostruito: "<< most_prob_Z <<endl;
      }
    }
    
    //nt->Fill(point.Z, most_prob_Z, point.Z -most_prob_Z);

    delete trackZ;
    
  }

  fout.Write();
  fout.Close();
}
