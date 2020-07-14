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

#define M 25
#define Z 14
#define ZS 10
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
  double zrange[Z] = {-17.0,-15.0,-13.0,-9.0,-5.0,-3.0,-1.0,1.0,3.0,5.0,9.0,13.0,15.0,17.0};
  double zrange_3sigma[ZS] = {-17.0,-13.0,-9.0,-5.0,-3.0,3.0,5.0,9.0,13.0,17.0};
  double resrz[Z-1] = {0.} , resrz_3sigma[ZS-1] = {0.} ;
  double resrz_error[Z-1] = {0.} , resrz_error_3sigma[ZS-1] = {0.} ;
  double zsimul[Z-1] = {0.} , zsimul_3sigma[ZS-1] = {0.} , eff_z[Z-1] = {0.} , effz_err[Z-1]={0.} , eff_z_3sigma[ZS-1] = {0.} , effz_err_3sigma[ZS-1]={0.};
  for(int k=0;k<Z-1;k++){
	   double countsim_effz = 0. , countrec_effz = 0. ;
	   for(int m=0;m<nrec;m++){
		   reconstr->GetEvent(m);
		   if((rec_sim>=zrange[k])&&(rec_sim<zrange[k+1])){
			   resolution->Fill(rec_sim-rec_z);
			   countrec_effz++;

			   /*if(abs(rec_sim)<= (3*sigma_z)){
				   resolz_3sigma->Fill(rec_sim-rec_z);
				   countrec_effz_3sigma++;
			   }*/
		   }
	   }
	   for(int l=0;l<nsim;l++){
		   simulat->GetEvent(l);
		   if((sim_z >= zrange[k]) && (sim_z < zrange[k+1])){
			   countsim_effz++;
			   /*if(abs(sim_z) <= (3*sigma_z)){
				   countsim_effz_3sigma++;
			   }*/
	   	  	}
	   }
	   zsimul[k] = (zrange[k] + zrange[k+1])/2.;
	   //zsimul_3sigma[k] = (zrange_3sigma[k] + zrange_3sigma[k+1])/2.;

	   //EFFICIENZA VS Z
	   if(countsim_effz != 0){
	   		   eff_z[k]=((double)countrec_effz)/countsim_effz;
	   		   effz_err[k]=sqrt(eff_z[k]*(1-eff_z[k])/countsim_effz);
	   	   }
	   //RISOLUZIONE VS Z
	   	   resrz[k] =1000*(resolution->GetRMS());
	   	   resrz_error[k] =1000*resolution->GetRMSError();
	   	   resolution->Reset();
  }
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________	 

//EFFICIENZA VS Z (Z<3SIGMA)

	   for(int k=0;k<ZS-1;k++){
	   	   double countrec_effz_3sigma = 0. ,countsim_effz_3sigma = 0. ;
	   	   for(int m=0;m<nrec;m++){
	   		   reconstr->GetEvent(m);
	   		   if((rec_sim>=zrange_3sigma[k])&&(rec_sim<zrange_3sigma[k+1])){
	   			   if(abs(rec_sim)<= (3*sigma_z)){
	   				   resolz_3sigma->Fill(rec_sim-rec_z);
	   				   countrec_effz_3sigma++;
	   			   }
	   		   }
	   	   }
	   	   for(int l=0;l<nsim;l++){
	   		   simulat->GetEvent(l);
	   		   if((sim_z >= zrange_3sigma[k]) && (sim_z < zrange_3sigma[k+1])){
	   			   if(abs(sim_z) <= (3*sigma_z)){
	   				   countsim_effz_3sigma++;
	   			   }
	   	   	  	}
	   	   }

	   	   zsimul_3sigma[k] = (zrange_3sigma[k] + zrange_3sigma[k+1])/2.;

	   	   if(countsim_effz_3sigma != 0){
	   	   		   eff_z_3sigma[k]=((double)countrec_effz_3sigma)/countsim_effz_3sigma;
	   	   		   effz_err_3sigma[k]=sqrt(eff_z[k]*(1-eff_z[k])/countsim_effz_3sigma);
	   	   	   }

	   //RISOLUZIONE VS Z (Z<3SIGMA)
	   resrz_3sigma[k] =1000*(resolz_3sigma->GetRMS());
	   resrz_error_3sigma[k] =1000*resolz_3sigma->GetRMSError();
	   resolz_3sigma->Reset();


}

  //Valutazione della risoluzione e dell'efficienza in funzione della molteplicità
  TH1D *eff = new TH1D("eff","Efficienza",55,0.5,55.5);
  eff->SetDirectory(0);
  double efficiency[M-1]={0.} , eff_mz[M-1]={0.} , resrm[M-1] = {0.} , resrm_error[M-1] = {0.} , resrm_3sigma[M-1] = {0.} , resrm_error_3sigma[M-1] = {0.};
  double mult_eff[M-1]={0.} , mult_effz[M-1]={0.} , mult_res[M-1]={0.} , multeff_err[M-1]={0.} , multeff_errz[M-1]={0.};
  TH1D *resmult_sigma = new TH1D("resmult_sigma","risolmul",201,-0.5,0.5);
  resmult_sigma->SetDirectory(0);
  TH1D *resmult_3sigma = new TH1D("resmult_3sigma","risolmul3sigma",201,-0.5,0.5);
   resmult_3sigma->SetDirectory(0);
double moltep[M]={2,3,5,6,7,9,10,12,15,17,18,20,22,25,28,30,32,35,38,40,42,45,48,50,52};

  for(int range=0;range<M-1;range++){
	  double count_sim_sigma = 0. , count_rec_sigma = 0. , count_sim_3sigma = 0., count_rec_3sigma = 0.;
//se il valore di molteplicità corrisponde si riempe il contatore degli eventi simulati
	  for(int i=0;i<nsim;i++){
		  simulat->GetEvent(i);
		  if((mult_tot>=moltep[range])&&(mult_tot<moltep[range+1])){
			  if(abs(sim_z) <= (sigma_z)){
				  count_sim_sigma++;
			  }
			  if(abs(sim_z) <= (3*sigma_z)){
			  	count_sim_3sigma++;
		  	  }
		  }
	  }
	  //se il valore di molteplicità corrisponde si aumenta il contatore degli eventi ricostruiti
	  for(int j=0;j<nrec;j++){
		  reconstr->GetEvent(j);
		  if((mult >= moltep[range])&&(mult < moltep[range+1])){
			if(abs(rec_sim) <= (sigma_z)){
				resmult_sigma->Fill(rec_sim-rec_z);
				count_rec_sigma++;
			 }
			if(abs(rec_sim)<= (3*sigma_z)){
		  		resmult_3sigma->Fill(rec_sim-rec_z);
		  		 count_rec_3sigma++;
		  	  	 }
		  	 }
	  	  }

	 mult_eff[range]=(moltep[range]+moltep[range+1])/2.;
	  //EFFICIENZA VS MOLTEPLICITÀ(Z<SIGMA)
	
	  if(count_sim_sigma!=0){
	 	  //cout<< "molteplicita" << mult_eff[range]<<endl;
	 	   efficiency[range]= ((double) count_rec_sigma) / count_sim_sigma;
	 	   multeff_err[range]=sqrt(efficiency[range]*(1-efficiency[range])/count_sim_sigma);
	 	   eff->Fill(mult_eff[range],efficiency[range]);
	 	  }
	
	mult_effz[range]=(moltep[range]+moltep[range+1])/2.;
	  //EFFICIENZA VS MOLTEPLICITÀ (Z<3SIGMA)
	  if(count_sim_3sigma!=0 ){
		  
		  eff_mz[range]=count_rec_3sigma/count_sim_3sigma;
		  multeff_errz[range]=sqrt(eff_mz[range]*(1-eff_mz[range])/count_sim_3sigma);
	  }
	 /* mult_res[range-1]=range;
cout<<mult_res[range-1]<<endl;*/

	  //RISOLUZIONE VS MOLTEPLICITÀ (Z<SIGMA)
//if(resmult_sigma->GetRMS()!=0){
		mult_res[range]=(moltep[range]+moltep[range+1])/2.;
	  resrm[range] = 10000*resmult_sigma->GetRMS();

	  resrm_error[range] = 10000*resmult_sigma->GetRMSError();
//}
	  resmult_sigma->Reset();

	  //RISOLUZIONE VS MOLTEPLICITÀ (Z<3SIGMA)
//if(resmult_3sigma->GetRMS()!=0){
	  resrm_3sigma[range] = 10000*resmult_3sigma->GetRMS();
	  resrm_error_3sigma[range] = 10000*resmult_3sigma->GetRMSError();
//}
	  resmult_3sigma->Reset();

  }
  //new TCanvas();
  //eff->Draw("hist");

  new TCanvas("Graph1","Efficiency vs Multiplicity (Z < sigma)",200,10,800,500);
  TGraphErrors *gr1=new TGraphErrors(M-1,mult_eff,efficiency,0,multeff_err);
  gr1->SetTitle("Efficiency Vs Multiplicity (Z < #sigma)");
  gr1->GetXaxis()->SetTitle("Multiplicity");
  gr1->GetYaxis()->SetTitle("Efficiency");
  gr1->SetMarkerColor(2);
  gr1->SetMarkerStyle(20);
  gr1->Draw("ALP");


  new TCanvas("Graph2","Resolution vs Z",200,10,800,500);
  TGraphErrors *grrz = new TGraphErrors(Z-1,zsimul,resrz,0,resrz_error);
  grrz->SetTitle("Resolution vs Simulated Z");
  grrz->GetXaxis()->SetTitle("Simulated Z [cm]");
  grrz->GetYaxis()->SetTitle("Resolution (:1000) [cm]");
  //grrz->SetMarkerColor(4);
  //grrz->Draw("A*L");
  grrz->SetMarkerColor(2);
  grrz->SetMarkerStyle(20);
  grrz->Draw("ALP");

  new TCanvas("Graph3","Resolution vs Z (Z < 3 sigma)",200,10,800,500);
   TGraphErrors *grrz_3sigma = new TGraphErrors(ZS-1,zsimul_3sigma,resrz_3sigma,0,resrz_error_3sigma);
   grrz_3sigma->SetTitle("Resolution vs Simulated Z (Z < 3#sigma)");
   grrz_3sigma->GetXaxis()->SetTitle("Simulated Z [cm]");
   grrz_3sigma->GetYaxis()->SetTitle("Resolution (:1000) [cm]");
   //grrz_3sigma->SetMarkerColor(4);
   //grrz_3sigma->Draw("A*L");
  grrz_3sigma->SetMarkerColor(2);
  grrz_3sigma->SetMarkerStyle(20);
  grrz_3sigma->Draw("ALP");

  new TCanvas("Graph4","Efficiency vs Z",200,10,800,500);
  TGraphErrors *gr4=new TGraphErrors(Z-1,zsimul,eff_z,0,effz_err);
  gr4->SetTitle("Efficiency vs Simulated Z");
  gr4->GetXaxis()->SetTitle("Simulated Z [cm]");
  gr4->GetYaxis()->SetTitle("Efficiency");
  gr4->SetMarkerColor(2);
	gr4->SetMarkerStyle(20);
  gr4->Draw("ALP");

  new TCanvas("Graph5","Efficiency vs Z (Z < 3sigma)",200,10,800,500);
    TGraphErrors *gr5=new TGraphErrors(ZS-1,zsimul_3sigma,eff_z_3sigma,0,effz_err_3sigma);
    gr5->SetTitle("Efficiency vs Simulated Z (Z < 3 #sigma)");
    gr5->GetXaxis()->SetTitle("Simulated Z [cm]");
    gr5->GetYaxis()->SetTitle("Efficiency");
    gr5->SetMarkerColor(2);
		gr5->SetMarkerStyle(20);
    gr5->Draw("ALP");


  new TCanvas("Graph6","Resolution Vs Multiplicity (Z < sigma)",200,10,800,500);
  TGraphErrors *grrmsigma = new TGraphErrors(M-1,mult_res,resrm,0,resrm_error);
  grrmsigma->SetTitle("Resolution vs Multiplicity(Z < #sigma)");
  grrmsigma->GetXaxis()->SetTitle("Multiplicity");
  grrmsigma->GetYaxis()->SetTitle("Resolution [#mum] ");
  grrmsigma->SetMarkerColor(2);
	grrmsigma->SetMarkerStyle(20);
  grrmsigma->Draw("ALP");

  new TCanvas("Graph7","Resolution Vs Multiplicity (Z < 3sigma)",200,10,800,500);
  TGraphErrors *grrm3sigma = new TGraphErrors(M-1,mult_res,resrm_3sigma,0,resrm_error_3sigma);
  grrm3sigma->SetTitle("Resolution vs Multiplicity(Z < 3#sigma)");
  grrm3sigma->GetXaxis()->SetTitle("Multiplicity");
  grrm3sigma->GetYaxis()->SetTitle("Resolution [#mum]");
	grrm3sigma->SetMarkerColor(2);
	grrm3sigma->SetMarkerStyle(20);
  grrm3sigma->Draw("ALP");

   new TCanvas("Graph8","Efficiency Vs Multiplicity (Z < 3sigma)",200,10,800,500);
   TGraphErrors *gr8=new TGraphErrors(M-1,mult_effz,eff_mz,0,multeff_errz);
   gr8->SetTitle("Efficiency Vs Multiplicity (Z < 3#sigma)");
   gr8->GetXaxis()->SetTitle("Multiplicity");
   gr8->GetYaxis()->SetTitle("Efficiency");
   gr8->SetMarkerColor(2);
	 gr8->SetMarkerStyle(20);
   gr8->Draw("ALP"); 
}
