#include "settings.h"

namespace Settings{
    void setNumThreads(int argc, char* argv[]){
        std::stringstream strStream {};
        std::string argOption {};
        int nThreads {0};
        switch (argc)
        {
        case 1:
            std::cout << "TO MULTI-THREAD (OpenMP) USE -nt OPTION FOLLOWED BY "
                << "THE NUMBER OF THREADS. Using a default thread count.\n";
            return;
        case 3:
            argOption = argv[1];
            if (argOption == "-nt"){
                strStream << argv[2];
                strStream >> nThreads;
                std::cout << "Number of threads set to: " << nThreads << '\n';
                omp_set_num_threads(nThreads);
            } else{
                std::cout << "Option " << argOption << " not recognized, "
                    << "USE -nt OPTION FOLLOWED BY THE NUMBER OF THREADS.\n";
            }
            return;
        default:
            // replace this with a throw 
            throw "OPTION NOT RECOGNIZED, USE -nt FOLLOWED BY THE NUMBER OF THREADS TO MULTI-THREAD (OpenMP).\n";
        }
    }
}