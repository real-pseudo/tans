#include <TROOT.h>
#include <Riostream.h>
#include <vector>
#include <cstdlib>
#include <time.h>
void loader() {   
	clock_t start=clock();
    std::vector<std::string> files {
        "geometry.cxx", 
        "Utility.cxx",
        "hit.cxx", 
        "particle.cxx"
        
    };
    
    for (const auto& current_file: files) {
        std::string command = ".L " + current_file + "+"; 
        gROOT->ProcessLine(command.c_str());
    }

    gROOT->ProcessLine(".X simulation.C");
    gROOT->ProcessLine(".X reconstruction.C");
    gROOT->ProcessLine(".X graphs.C");
    clock_t end=clock();
    	cout<<"Execution time: "<<((double)(end-start)/CLOCKS_PER_SEC)<<"sec"<<endl;
}
