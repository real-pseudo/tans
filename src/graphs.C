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
#define Nsigma 1
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
  TH1D *resolz_3sigma = new TH1D("ris","risol",201,-0.5,0.5);
    resolz_3sigma->SetDirectory(0);
  double zrange[Z] = {-17.0,-13.0,-9.0,-5.0,-3.0,-1.0,1.0,3.0,5.0,9.0,13.0,17.0};
  double resrz[Z-1] = {0.} , resrz_3sigma[Z-1] = {0.} ;
  double resrz_error[Z-1] = {0.} , resrz_error_3sigma[Z-1] = {0.} ;
  double zsimul[Z-1] = {0.} , eff_z[Z-1] = {0.} , effz_err[Z-1]={0.} , eff_z_3sigma[Z-1] = {0.} , effz_err_3sigma[Z-1]={0.};
  for(int k=0;k<Z-1;k++){
	   double countsim_effz = 0. , countrec_effz = 0. , countrec_effz_3sigma = 0. ,countsim_effz_3sigma = 0. ;
	   for(int m=0;m<nrec;m++){
		   reconstr->GetEvent(m);
		   if((rec_sim>=zrange[k])&&(rec_sim<zrange[k+1])){
			   resolution->Fill(rec_sim-rec_z);
			   countrec_effz++;

			   if(abs(rec_sim)< (3.*sigma_z)){
				   resolz_3sigma->Fill(rec_sim-rec_z);
				   countrec_effz_3sigma++;
			   }
		   }
	   }
	   for(int l=0;l<nsim;l++){
		   simulat->GetEvent(l);
		   if((sim_z >= zrange[k]) && (sim_z < zrange[k+1])){
			   countsim_effz++;
			   if(abs(sim_z) < (3.*sigma_z)){
				   countsim_effz_3sigma++;
			   }
	   	  	}
	   }
	   zsimul[k] = (zrange[k] + zrange[k+1])/2.;
	   //EFFICIENZA VS Z
	   if(countsim_effz != 0){
	   		   eff_z[k]=((double)countrec_effz)/countsim_effz;
	   		   effz_err[k]=sqrt(eff_z[k]*(1-eff_z[k])/countsim_effz);
	   	   }
	   //EFFICIENZA VS Z (Z<3SIGMA)
	   	   if(countsim_effz_3sigma != 0){
	   	   		   eff_z_3sigma[k]=((double)countrec_effz_3sigma)/countsim_effz_3sigma;
	   	   		   effz_err_3sigma[k]=sqrt(eff_z[k]*(1-eff_z[k])/countsim_effz_3sigma);
	   	   	   }
	   //RISOLUZIONE VS Z
	   resrz[k] =10*(resolution->GetRMS());
	   resrz_error[k] =10*resolution->GetRMSError();
	   resolution->Reset();
	   //RISOLUZIONE VS Z (Z<3SIGMA)
	   resrz_3sigma[k] =10*(resolz_3sigma->GetRMS());
	   resrz_error_3sigma[k] =10*resolz_3sigma->GetRMSError();
	   resolz_3sigma->Reset();


}

  //Valutazione della risoluzione e dell'efficienza in funzione della molteplicità
  TH1D *eff = new TH1D("eff","Efficienza",55,0.5,55.5);
  eff->SetDirectory(0);
  double efficiency[M]={0.} , eff_mz[M]={0.} , resrm[M] = {0.} , resrm_error[M] = {0.} , resrm_3sigma[M] = {0.} , resrm_error_3sigma[M] = {0.};
  double mult_eff[M]={0.} , mult_effz[M]={0.} , mult_res[M]={0.} , multeff_err[M]={0.} , multeff_errz[M]={0.} ;
  TH1D *resmult_sigma = new TH1D("resmult_sigma","risolmul",201,-0.5,0.5);
  resmult_sigma->SetDirectory(0);
  TH1D *resmult_3sigma = new TH1D("resmult_3sigma","risolmul3sigma",201,-0.5,0.5);
   resmult_3sigma->SetDirectory(0);

  for(int range=1;range<=M;range++){
	  double count_sim_sigma = 0. , count_rec_sigma = 0. , count_sim_3sigma = 0., count_rec_3sigma = 0.;
//se il valore di molteplicità corrisponde si riempe il contatore degli eventi simulati
	  for(int i=0;i<nsim;i++){
		  simulat->GetEvent(i);
		  if(mult_tot==range){
			  if(abs(sim_z) < (sigma_z)){
				  count_sim_sigma++;
			  }
			  if(abs(sim_z) < (3.*sigma_z)){
			  	count_sim_3sigma++;
		  	  }
		  }
	  }
	  //se il valore di molteplicità corrisponde si aumenta il contatore degli eventi ricostruiti
	  for(int j=0;j<nrec;j++){
		  reconstr->GetEvent(j);
		  if(mult==range){
			if(abs(rec_sim) < (sigma_z)){
				resmult_sigma->Fill(rec_sim-rec_z);
				count_rec_sigma++;
			 }
			if(abs(rec_sim)< (3.*sigma_z)){
		  		resmult_3sigma->Fill(rec_sim-rec_z);
		  		 count_rec_3sigma++;
		  	  	 }
		  	 }
	  	  }
	  //EFFICINZA VS MOLTEPLICITÀ(Z<SIGMA)
	  if(count_sim_sigma!=0){
	 	  mult_eff[range-1]=range;
	 	  //cout<< "molteplicita" << mult_eff[range-1]<<endl;
	 	   efficiency[range-1]= ((double) count_rec_sigma) / count_sim_sigma;
	 	   multeff_err[range-1]=sqrt(efficiency[range-1]*(1-efficiency[range-1])/count_sim_sigma);
	 	   eff->Fill(mult_eff[range-1],efficiency[range-1]);
	 	  }
	  //EFFICINZA VS MOLTEPLICITÀ (Z<3SIGMA)
	  if(count_sim_3sigma!=0 ){
		  mult_effz[range-1]=range;
		  eff_mz[range-1]=count_rec_3sigma/count_sim_3sigma;
		  multeff_errz[range-1]=sqrt(eff_mz[range-1]*(1-eff_mz[range-1])/count_sim_3sigma);
	  }
	  mult_res[range-1]=range;
	  //RISOLUZIONE VS MOLTEPLICITÀ (Z<SIGMA)
	  resrm[range-1] = 100*resmult_sigma->GetRMS();
	  resrm_error[range-1] = 100*resmult_sigma->GetRMSError();
	  resmult_sigma->Reset();

	  //RISOLUZIONE VS MOLTEPLICITÀ (Z<3SIGMA)
	  resrm_3sigma[range-1] = 100*resmult_3sigma->GetRMS();
	  resrm_error_3sigma[range-1] = 100*resmult_3sigma->GetRMSError();
	  resmult_3sigma->Reset();

  }
  //new TCanvas();
  //eff->Draw("hist");

  new TCanvas("Graph1","Efficiency vs Multiplicity (Z < sigma)",200,10,800,500);
  TGraphErrors *gr1=new TGraphErrors(51,mult_eff,efficiency,0,multeff_err);
  gr1->SetTitle("Efficiency Vs Multiplicity (Z < sigma)");
  gr1->GetXaxis()->SetTitle("Multiplicity");
  gr1->GetYaxis()->SetTitle("Efficiency");
  gr1->Draw("AP*");


  new TCanvas("Graph2","Resolution vs Z",200,10,800,500);
  TGraphErrors *grrz = new TGraphErrors(11,zsimul,resrz,0,resrz_error);
  grrz->SetTitle("Resolution vs Simulated Z");
  grrz->GetXaxis()->SetTitle("Simulated Z [cm]");
  grrz->GetYaxis()->SetTitle("Resolution ");
  grrz->Draw("AP*");

  new TCanvas("Graph3","Resolution vs Z (Z < 3sigma)",200,10,800,500);
   TGraphErrors *grrz_3sigma = new TGraphErrors(11,zsimul,resrz_3sigma,0,resrz_error_3sigma);
   grrz_3sigma->SetTitle("Resolution vs Simulated Z (Z < 3sigma)");
   grrz_3sigma->GetXaxis()->SetTitle("Simulated Z [cm]");
   grrz_3sigma->GetYaxis()->SetTitle("Resolution ");
   grrz_3sigma->Draw("AP*");

  new TCanvas("Graph4","Efficiency vs Z",200,10,800,500);
  TGraphErrors *gr4=new TGraphErrors(11,zsimul,eff_z,0,effz_err);
  gr4->SetTitle("Efficiency vs Simulated Z");
  gr4->GetXaxis()->SetTitle("Simulated Z [cm]");
  gr4->GetYaxis()->SetTitle("Efficiency");
  gr4->Draw("AP*");

  new TCanvas("Graph5","Efficiency vs Z (Z < 3sigma)",200,10,800,500);
    TGraphErrors *gr5=new TGraphErrors(11,zsimul,eff_z_3sigma,0,effz_err_3sigma);
    gr5->SetTitle("Efficiency vs Simulated Z (Z < 3sigma)");
    gr5->GetXaxis()->SetTitle("Simulated Z [cm]");
    gr5->GetYaxis()->SetTitle("Efficiency");
    gr5->Draw("AP*");


  new TCanvas("Graph6","Resolution Vs Multiplicity (Z < sigma)",200,10,800,500);
  TGraphErrors *grrmsigma = new TGraphErrors(51,mult_res,resrm,0,resrm_error);
  grrmsigma->SetTitle("Resolution vs Multiplicity(Z < sigma)");
  grrmsigma->GetXaxis()->SetTitle("Multiplicity");
  grrmsigma->GetYaxis()->SetTitle("Resolution ");
  grrmsigma->Draw("AP*");

  new TCanvas("Graph7","Resolution Vs Multiplicity (Z < 3sigma)",200,10,800,500);
  TGraphErrors *grrm3sigma = new TGraphErrors(51,mult_res,resrm_3sigma,0,resrm_error_3sigma);
  grrm3sigma->SetTitle("Resolution vs Multiplicity(Z < 3sigma)");
  grrm3sigma->GetXaxis()->SetTitle("Multiplicity");
  grrm3sigma->GetYaxis()->SetTitle("Resolution ");
  grrm3sigma->Draw("AP*");

   new TCanvas("Graph8","Efficiency Vs Multiplicity (Z < 3sigma)",200,10,800,500);
   TGraphErrors *gr8=new TGraphErrors(51,mult_effz,eff_mz,0,multeff_errz);
   gr8->SetTitle("Efficiency Vs Multiplicity (Z < 3sigma)");
   gr8->GetXaxis()->SetTitle("Multiplicity");
   gr8->GetYaxis()->SetTitle("Efficiency");
   gr8->Draw("AP*");
}
