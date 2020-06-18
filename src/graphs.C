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

void graphs2() {
  
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
  double efficiency[52]={0.} , eff_multz[52]={0.};
  double mult_eff[52]={0.} , mult_int[52]={0.}  ;
  //double multiplicity[n]={0.};//{2,4,6,8,10,15,20,25,30,35,40,45,50,55};//range di molteplicità
  TH1D *eff = new TH1D("eff","Efficienza",55,0.5,55.5);
  eff->SetDirectory(0);
  TH1D *ris = new TH1D("ris","risol",201,-0.5,0.5);
    ris->SetDirectory(0);
    TH1D *resmult = new TH1D("resmult","risolmul",201,-0.3,0.3);
        resmult->SetDirectory(0);
  double zrange[12] = {-17.0,-13.0,-9.0,-5.0,-3.0,-1.0,1.0,3.0,5.0,9.0,13.0,17.0};//"Resolution vs Simulated Z"
  double zsimul[11] = {0.};
  double resrz[11] = {0.} , resrm[52] = {0.} ,resrm_error[52] = {0.} ;
   double resrz_error[11] = {0.} , eff_z[11] = {0.};
   for(int k=0;k<11;k++){
	   double countsim_effz = 0. , countrec_effz = 0. ;
	   for(int m=0;m<nrec;m++){
		   reconstr->GetEvent(m);

		   if((rec_sim>=zrange[k])&&(rec_sim<zrange[k+1])){
			   ris->Fill(rec_z-rec_sim);
			   countrec_effz++;
		   }
	   }
	   for(int l=0;l<nsim;l++){
		   simulat->GetEvent(l);
		   if((sim_z >= zrange[k]) && (sim_z < zrange[k+1])){
			   countsim_effz++;
	   	  	}
	   }
	resrz[k] =10*(ris->GetRMS());
	resrz_error[k] =10*ris->GetRMSError();
	zsimul[k] = (zrange[k] + zrange[k+1])/2.;
	if(countsim_effz != 0){
		eff_z[k]=countrec_effz/countsim_effz;
	}
ris->Reset();
}


  for(int range=1;range<=52;range++){
	  double count_sim = 0. , count_rec = 0. , count_sigmatot = 0., count_sigma = 0.;
//se il valore di molteplicità corrisponde si riempe il contatore degli eventi simulati
	  for(int i=0;i<nsim;i++){
		  simulat->GetEvent(i);
		  if(mult_tot==range){
			  count_sim++;
			  if(abs(sim_z) < (sigma_z)){
			  	count_sigmatot++;
		  	  }
		  }
	  }
	  //se il valore di molteplicità corrisponde si aumenta il contatore degli eventi ricostruiti
	  for(int j=0;j<nrec;j++){
		  reconstr->GetEvent(j);
		  if(mult==range){
		  	 count_rec++;
		  	resmult->Fill(rec_sim-rec_z);
		  	if(abs(rec_sim)< (sigma_z)){
		  		count_sigma++;
		  	  	 }
		  	 }
	  	  }
	  //SERVE PER EFFICIENZA dipendente da SIGMA_Z
	  //cout<<"count_sim "<<count_sigmatot<<"count_rec:"<<count_sigma<<endl;
	  if(count_sigmatot!=0 ){
	  mult_int[range-1]=range;
	  eff_multz[range-1]=count_sigma/count_sigmatot;

	  }
	  if(count_sim!=0){
	  	  mult_eff[range-1]=range;
	  	  //cout<< "molteplicita" << mult_eff[range-1]<<endl;
	  	  efficiency[range-1]= ((double) count_rec) / count_sim;
        /*cout << "conteggi registrati: " << count_rec << "---- conteggi simulati: "<<
                  count_sim << "efficienza" <<efficiency[range-1]<<endl;*/
	  	  eff->Fill(mult_eff[range-1],efficiency[range-1]);
	  }
	  resrm[range-1] = 10*resmult->GetRMS();//resolution in mm
	  		resrm_error[range-1] = 10*resmult->GetRMSError();//resolution error in mm
	  		resmult->Reset();

  }
  new TCanvas();
  eff->Draw("hist");
  new TCanvas("Graph1","Efficienza Vs Molteplicita'",200,10,800,500);
  TGraph *gr1=new TGraph(51,mult_eff,efficiency);
  gr1->SetTitle("Efficienza vs molteplicita'");
  gr1->GetXaxis()->SetTitle("Molteplicita'");
  gr1->GetYaxis()->SetTitle("Efficienza");
  gr1->Draw("AP*");
  /*new TCanvas("Graph4","res Vs Molteplicita'",200,10,800,500);
  TGraph *gr1=new TGraph(51,mult_eff,efficiency);
  gr1->SetTitle("Efficienza vs molteplicita'");
  gr1->GetXaxis()->SetTitle("Molteplicita'");
  gr1->GetYaxis()->SetTitle("Efficienza");
  gr1->Draw("AP*");*/

  new TCanvas("Graph","Resolution Vs Z",200,10,800,500);
  TGraphErrors *grrz = new TGraphErrors(11,zsimul,resrz,0,resrz_error);//"Resolution vs Simulated Z" TGraphErrors
  grrz->SetTitle("Resolution vs Simulated Z");
  grrz->GetXaxis()->SetTitle("Simulated Z [cm]");
  grrz->GetYaxis()->SetTitle("Resolution ");
  grrz->Draw("AP*");

  new TCanvas("Graph3","Efficiency Vs Z",200,10,800,500);
  TGraph *gr3=new TGraph(11,zsimul,eff_z);
  gr3->SetTitle("Efficiency vs Simulated Z");
  gr3->GetXaxis()->SetTitle("Simulated Z [cm]");
  gr3->GetYaxis()->SetTitle("Efficiency");
  gr3->Draw("AP*");


   new TCanvas("Graph2","efficienza Vs Molteplicita'(Z < sigma)",200,10,800,500);
   TGraph *gr2=new TGraph(51,mult_int,eff_multz);
   gr2->SetTitle("Efficienza vs molteplcita'(Z < sigma)");
   gr2->GetXaxis()->SetTitle("Molteplicita'");
   gr2->GetYaxis()->SetTitle("Efficienza ");
   gr2->SetMarkerStyle(21);
   gr2->Draw("APL");
}
