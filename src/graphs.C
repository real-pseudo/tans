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

void graphs() {
  
  //grafica
  //leggo ntuple
  TFile filin("vtxreco.root");
  float rec_sim,rec_z,rec_diff,mult,mult_tot,sim_z;
  TNtuple *simulat = (TNtuple*)filin.Get("z_sim");
        simulat->SetBranchAddress("zsimtot",&sim_z);
        simulat->SetBranchAddress("multtot",&mult_tot);
  TNtuple *reconstr = (TNtuple*)filin.Get("nt_rec");
      reconstr->SetBranchAddress("zsim",&rec_sim);
      reconstr->SetBranchAddress("zrec",&rec_z);
      reconstr->SetBranchAddress("diff",&rec_diff);
      reconstr->SetBranchAddress("mult",&mult);


  double sigma_z = 5.3; //caratteristico del vtx generato
  int nsim=simulat->GetEntries();
  int nrec=reconstr->GetEntries();
  double efficiency[52]={0.};
  double mult_eff[52]={0.};
  //double multiplicity[n]={0.};//{2,4,6,8,10,15,20,25,30,35,40,45,50,55};//range di molteplicità
  TH1D *eff = new TH1D("eff","Efficienza",55,0.5,55.5);
  eff->SetDirectory(0);

  for(int range=1;range<=52;range++){
	  double count_sim = 0. , count_rec = 0. , count_sigmatot = 0., count_sigma = 0.;
//se il valore di molteplicità corrisponde si riempe il contatore degli eventi simulati
	  for(int i=0;i<nsim;i++){
		  simulat->GetEvent(i);
		  if(mult_tot==range){
			  count_sim++;
			  /*if(abs(sim_z) < (sigma_z)){
			  		count_sigmatot++;
		  	  }*/
		  }
	  }
	  //se il valore di molteplicità corrisponde si aumenta il contatore degli eventi ricostruiti
	  for(int j=0;j<nrec;j++){
		  reconstr->GetEvent(j);
		  if(mult==range){
		  	 count_rec++;
		  	/*if(abs(rec_sim)< (sigma_z)){
		  		count_sigma++;
		  	  	 }*/
		  	 }
	  	  }
	  //SERVE PER EFFICIENZA dipendente da SIGMA_Z
	  //cout<<"count_sim "<<count_sigmatot<<"count_rec:"<<count_sigma<<endl;
	  /*if(count_sigmatot!=0 && count_sigma!=0){
	  mult_int[range-1]=range;
	  efficiency[range-1]=count_sigma/count_sigmatot;
	  eff->Fill(mult_eff[range-1],efficiency[range-1]);
	  }*/
	  if(count_sim!=0 && count_rec!=0){
	  	  mult_eff[range-1]=range;
	  	  cout<<mult_eff[range-1]<<endl;
	  	  efficiency[range-1]=count_rec/count_sim;
	  	  eff->Fill(mult_eff[range-1],efficiency[range-1]);
	  }
  }
  new TCanvas();
  eff->Draw("hist");
  new TCanvas("Graph1","Efficienza Vs Molteplicita'",200,10,800,500);
  TGraph *gr1=new TGraph(52,mult_eff,efficiency);
  gr1->SetTitle("Efficienza vs molteplicita'");
  gr1->GetXaxis()->SetTitle("Molteplicita'");
  gr1->GetYaxis()->SetTitle("Efficienza");
  gr1->Draw("AP*");

  /* new TCanvas("Graph2","efficienza Vs Molteplicita'(Z < sigma)",200,10,800,500);
   TGraph *gr2=new TGraph(n-1,mult_int,eff);
   gr2->SetTitle("Efficienza vs molteplcita'(Z < sigma)");
   gr2->GetXaxis()->SetTitle("Molteplicita'");
   gr2->GetYaxis()->SetTitle("Efficienza ");
   gr2->SetMarkerStyle(21);
   gr2->Draw("APL");*/
}
