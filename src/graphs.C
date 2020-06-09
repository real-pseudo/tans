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
  const int n=14;
  double efficiency[14]={0.} , eff[14]={0};
  double mult_int[14]={0.};
  int multiplicity[14]={2,4,6,8,10,15,20,25,30,35,40,45,50,55};//range di molteplicità
  for(int range=0;range<n;range++){
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
}