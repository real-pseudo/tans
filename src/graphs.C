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

#define M 52
#define Z 12
void graphs() {

  //lettura ntuple
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

  //double multiplicity[n]={0.};//{2,4,6,8,10,15,20,25,30,35,40,45,50,55};//range di molteplicità

  	  	  //Valutazione della risoluzione e dell'efficienza in funzione dello Z simulato
  TH1D *resolution = new TH1D("ris","risol",201,-0.5,0.5);
  resolution->SetDirectory(0);
  double zrange[Z] = {-17.0,-13.0,-9.0,-5.0,-3.0,-1.0,1.0,3.0,5.0,9.0,13.0,17.0};
  double resrz[Z-1] = {0.} ;
  double resrz_error[Z-1] = {0.} ;
  double zsimul[Z-1] = {0.} , eff_z[Z-1] = {0.} , effz_err[Z-1]={0.};
  for(int k=0;k<Z-1;k++){
	   double countsim_effz = 0. , countrec_effz = 0. ;
	   for(int m=0;m<nrec;m++){
		   reconstr->GetEvent(m);
		   if((rec_sim>=zrange[k])&&(rec_sim<zrange[k+1])){
			   resolution->Fill(rec_z-rec_sim);
			   countrec_effz++;
		   }
	   }
	   for(int l=0;l<nsim;l++){
		   simulat->GetEvent(l);
		   if((sim_z >= zrange[k]) && (sim_z < zrange[k+1])){
			   countsim_effz++;
	   	  	}
	   }
	   resrz[k] =10*(resolution->GetRMS());
	   resrz_error[k] =10*resolution->GetRMSError();
	   resolution->Reset();
	   zsimul[k] = (zrange[k] + zrange[k+1])/2.;
	   if(countsim_effz != 0){
		   eff_z[k]=((double)countrec_effz)/countsim_effz;
		   effz_err[k]=sqrt(eff_z[k]*(1-eff_z[k])/countsim_effz);
	   }
}

  //Valutazione della risoluzione e dell'efficienza in funzione della molteplicità
  TH1D *eff = new TH1D("eff","Efficienza",55,0.5,55.5);
  eff->SetDirectory(0);
  double efficiency[M]={0.} , eff_mz[M]={0.} , resrm[M] = {0.} , resrm_error[M] = {0.};
  double mult_eff[M]={0.} , mult_effz[M]={0.} , mult_res[M]={0.} , multeff_err[M]={0.} , multeff_errz[M]={0.} ;
  TH1D *resmult = new TH1D("resmult","risolmul",201,-0.5,0.5);
  resmult->SetDirectory(0);

  for(int range=1;range<=M;range++){
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
		  	 if(abs(rec_sim)< (sigma_z)){  //EFFICIENZA dipendente da SIGMA_Z
		  		count_sigma++;
		  	  	 }
		  	 }
	  	  }
	  //EFFICINZA VS MOLTEPLICITÀ
	  if(count_sim!=0){
	 	  mult_eff[range-1]=range;
	 	  //cout<< "molteplicita" << mult_eff[range-1]<<endl;
	 	   efficiency[range-1]= ((double) count_rec) / count_sim;
	 	   multeff_err[range-1]=sqrt(efficiency[range-1]*(1-efficiency[range-1])/count_sim);
	 	  	  eff->Fill(mult_eff[range-1],efficiency[range-1]);
	 	  }
	  //EFFICINZA VS MOLTEPLICITÀ (Z<SIGMA)
	  if(count_sigmatot!=0 ){
		  mult_effz[range-1]=range;
		  eff_mz[range-1]=count_sigma/count_sigmatot;
		  multeff_errz[range-1]=sqrt(eff_mz[range-1]*(1-eff_mz[range-1])/count_sigmatot);
	  }

	  mult_res[range-1]=range;
	  resrm[range-1] = 100*resmult->GetRMS();
	  resrm_error[range-1] = 100*resmult->GetRMSError();
	  resmult->Reset();

  }
  new TCanvas();
  eff->Draw("hist");

  new TCanvas("Graph1","Efficiency vs Multiplicity",200,10,800,500);
  TGraphErrors *gr1=new TGraphErrors(51,mult_eff,efficiency,0,multeff_err);
  gr1->SetTitle("Efficiency Vs Multiplicity");
  gr1->GetXaxis()->SetTitle("Multiplicity");
  gr1->GetYaxis()->SetTitle("Efficiency");
  gr1->Draw("AP*");


  new TCanvas("Graph2","Resolution vs Z",200,10,800,500);
  TGraphErrors *grrz = new TGraphErrors(11,zsimul,resrz,0,resrz_error);//"Resolution vs Simulated Z" TGraphErrors
  grrz->SetTitle("Resolution vs Simulated Z");
  grrz->GetXaxis()->SetTitle("Simulated Z [cm]");
  grrz->GetYaxis()->SetTitle("Resolution ");
  grrz->Draw("AP*");

  new TCanvas("Graph3","Efficiency vs Z",200,10,800,500);
  TGraphErrors *gr3=new TGraphErrors(11,zsimul,eff_z,0,effz_err);
  gr3->SetTitle("Efficiency vs Simulated Z");
  gr3->GetXaxis()->SetTitle("Simulated Z [cm]");
  gr3->GetYaxis()->SetTitle("Efficiency");
  gr3->Draw("AP*");


  new TCanvas("Graph4","Resolution Vs Multiplicity",200,10,800,500);
  TGraphErrors *grrm = new TGraphErrors(51,mult_res,resrm,0,resrm_error);
  grrm->SetTitle("Resolution vs Multiplicity");
  grrm->GetXaxis()->SetTitle("Multiplicity");
  grrm->GetYaxis()->SetTitle("Resolution ");
  grrm->Draw("AP*");

   new TCanvas("Graph5","Efficiency Vs Multiplicity(Z < sigma)",200,10,800,500);
   TGraphErrors *gr5=new TGraphErrors(51,mult_effz,eff_mz,0,multeff_errz);
   gr5->SetTitle("Efficiency Vs Multiplicity(Z < sigma)");
   gr5->GetXaxis()->SetTitle("Multiplicity");
   gr5->GetYaxis()->SetTitle("Efficiency");
   gr5->Draw("AP*");
}
