#include <TROOT.h>
#include <Riostream.h>
#include <vector>


void loader() {   
    std::vector<std::string> files {
        "Cilindro.cxx", 
        "Utility.cxx",
        "hit.cxx", 
        "Particella.cxx"
    };
    
    for (const auto& current_file: files) {
        std::string command = ".L " + current_file + "+"; 
        gROOT->ProcessLine(command.c_str());
    }

    gROOT->ProcessLine(".X GeneraTree.C");
}
