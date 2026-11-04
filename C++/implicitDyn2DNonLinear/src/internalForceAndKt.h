#ifndef INTERNAL_FORCE_AND_KT_H
#define INTERNAL_FORCE_AND_KT_H

#include <array>
#include <omp.h>
#include <algorithm>
#include "Eigen/Core"
#include "settings.h"
#include "material.h"
#include "Mesh.h"
#include "massMatrix.h"

namespace FintKtInt{
    using SizT = std::size_t;
    using nDofVectf = Eigen::Vector<Settings::precision,Eigen::Dynamic>;
    using nDofElVectf = Eigen::Vector<Settings::precision,Settings::nDofEl>;
    using nNodesRowVecti = Eigen::RowVector<int,Settings::nNodesEl>;
    using sparseMat = Eigen::SparseMatrix<Settings::precision>;
    using elMat = Eigen::Matrix<Settings::precision,Settings::nDofEl
        ,Settings::nDofEl>;
    using dNEigen = Eigen::Matrix<Settings::precision,Settings::nNodesEl
        ,Settings::nDofEl>;
    using invJEigen = Eigen::Matrix<Settings::precision,Settings::nNodesEl
        ,Settings::nNodesEl>;
    using B0Eigen = Eigen::Matrix<Settings::precision,Settings::nNodesEl
        ,Settings::nDofEl>;
    using voigtTensorFull = Eigen::Vector<Settings::precision
        ,Settings::dofNode*Settings::dofNode>;
    using voigtTensor = Eigen::Vector<Settings::precision,3>;
    using tensor = Eigen::Matrix<Settings::precision,Settings::dofNode
        ,Settings::dofNode>;
    using matMultipliedByVel = Eigen::Matrix<Settings::precision,3
        ,Settings::nDofEl>;
    using matMultipliedByVelFull = Eigen::Matrix<Settings::precision,4
        ,Settings::nDofEl>;

    class DNprecomputed{
    private:
        std::array<dNEigen,4uz> m_dNs;
    public:
        DNprecomputed()
        {
            SizT indDN {0uz};
            for (const auto xi : mMat::gaussianPts){
                for (const auto eta : mMat::gaussianPts){
                    m_dNs[indDN] << (-1.0+eta),0.0,(1.0-eta),0.0,(1.0+eta),0.0
                        ,(-1.0-eta),0.0,(-1.0+xi),0.0,(-1.0-xi),0.0,(1.0+xi)
                        ,0.0,(1.0-xi),0.0,0.0,(-1.0+eta),0.0,(1.0-eta),0.0
                        ,(1.0+eta),0.0,(-1.0-eta),0.0,(-1.0+xi),0.0,(-1.0-xi)
                        ,0.0,(1.0+xi),0.0,(1.0-xi);
                    m_dNs[indDN++] *= 0.25;
                }
            }
        }
        const dNEigen& getdNi(SizT indx) const{ return m_dNs[indx]; }
    };

    class KtIntCache{
    private:
        int m_nElements {};
        int m_rows {};
        int m_cols {};
        int m_nnz  {};
        sparseMat::StorageIndex* m_innerIndices {};
        sparseMat::StorageIndex* m_outerIndices {};
        Settings::precision* m_values {};
        bool m_isLoopInnerIndxCached {false};
        std::vector<sparseMat::StorageIndex> m_cachedLoopInnerIndx;
    public:
        KtIntCache() = default;

        bool isCached() const{
            return (m_innerIndices == nullptr) ? false : true;
        }
        bool isInnerIndxCached() const{ return m_isLoopInnerIndxCached; }
        void setInnerIndxToCached(){ m_isLoopInnerIndxCached = true; }
        void setupCache(sparseMat& ktInt){
            m_rows = ktInt.rows();
            m_cols = ktInt.cols();
            m_nnz  = ktInt.nonZeros();
            m_innerIndices = ktInt.innerIndexPtr();
            m_outerIndices = ktInt.outerIndexPtr();
            m_values = ktInt.valuePtr();
        }
        void checkAndResetValues(const sparseMat& ktInt, const Mesh& mesh);

        void addContribution2El(const sparseMat::StorageIndex row
            , const sparseMat::StorageIndex col, const SizT elIndx
            , const sparseMat::StorageIndex localContigIndx
            , const Settings::precision val);
    };

    void compute(const Mesh& mesh, const nDofVectf& dispNew, nDofVectf& fInt
        , sparseMat& ktInt, KtIntCache& ktIntCache);
}


#endif