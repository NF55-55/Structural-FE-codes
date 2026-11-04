#include "estimateMemoryUsage.h"

double estimateMemUsageMB(const int nDof, const int nElements){
    int meshSize = nElements*Settings::nNodesEl*sizeof(int) 
        + nDof*sizeof(double);
    int forceAndDispSize = 3*nDof*sizeof(Settings::precision); // considering 
    // temp for splitting solution
    int mapDofSize = 2*nDof*sizeof(int); // included the permutation matrix
    int tripletListSize = 3*(2*nDof*Settings::nDofEl)*sizeof(Settings::
        precision);
    int kGSize = 2*((2*nDof*Settings::nDofEl)*(sizeof(Settings::precision)
        +sizeof(int)) + nDof*3*sizeof(int)); // considering temp copy for 
        // splitting solution
    return ((meshSize + forceAndDispSize + mapDofSize + tripletListSize 
        + kGSize)/1.0e6*Settings::multiplierMemory4Solving);
}