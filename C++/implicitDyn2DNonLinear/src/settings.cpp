#include "settings.h"

namespace Settings{
    void setNumThreadsFromUserInput(int argc, char* argv[]){
        std::string option {};
        switch (argc)
        {
        case 1:
            std::cout << "TO MULTI-THREAD (OpenMP) USE -nt OPTION FOLLOWED BY "
                << "THE NUMBER OF THREADS. Using a SINGLE THREAD.\n";
            omp_set_num_threads(1);
            nThreads = 1;
            return;
        case 3:
            option = *(argv+1);
            if (option == "-nt"){
                std::stringstream strStream;
                strStream << *(argv+2);
                int numThreads {};
                strStream >> numThreads;
                omp_set_num_threads(numThreads);
                nThreads = numThreads;
                #pragma omp parallel master
                {
                std::cout << "Number of threads set to: " 
                    << omp_get_num_threads() << " as requested.\n";
                }
            }else throw "Option not recognized, available ones: '-nt X'\n";
            return;
        default:
            throw "OPTION NOT RECOGNIZED, USE -nt FOLLOWED BY THE NUMBER OF THREADS TO MULTI-THREAD (OpenMP).\n";
        }
    }
}