# include<iostream>
# include<random>

# include "FaceAddressedMesh2D.hpp"
# include "MeshFileReader2D.hpp"
# include "LinearVector.hpp"

# include "DenseLinearSystem.hpp"
# include "DenseMatrix.hpp"

# include "SparseAddress.hpp"
# include "SparseLinearSystem.hpp"
# include "SparseMatrix.hpp"


# include "Diagnostics.hpp"
# include "PathConfig.hpp"

using namespace std ;

/* *******************************************************************
    Main driver for running all diagnostics for the dense and 
    sparse matrix performance study.
******************************************************************* */

int main(){

    // Stores the input and output paths used by the diagnostics routines.
    PathConfig paths ;

    // Creates the diagnostics object that runs all algebraic operation benchmarks.
    Diagnostics diagnostics ;

    // Runs vector-operation diagnostics, including addition, subtraction, scaling, and norm evaluation.
    diagnostics.vectorAdditionDiagnostics() ;
    diagnostics.vectorSubstractionDiagnostics() ;
    diagnostics.vectorScalarMultiplicationDiagnostics() ;
    diagnostics.vectorOneNormDiagnostics() ;
    diagnostics.vectorTwoNormDiagnostics() ;
    diagnostics.vectorInfinityNormDiagnostics() ;

    // Runs sparse-matrix diagnostics, including addition, subtraction, and scalar multiplication.
    diagnostics.sparseMatrixAdditionDiagnostics() ;
    diagnostics.sparseMatrixSubstractionDiagnostics() ;
    diagnostics.sparseMatrixScalarMultiplicationDiagnostics() ;

    // Runs dense-matrix diagnostics, including addition, subtraction, and scalar multiplication.
    diagnostics.denseMatrixAdditionDiagnostics() ;
    diagnostics.denseMatrixSubstractionDiagnostics() ;
    diagnostics.denseMatrixScalarMultiplicationDiagnostics() ;

    // Runs sparse linear-system diagnostics, including algebraic operations, matrix-vector multiplication, and residual evaluation.
    diagnostics.slsAddition() ;
    diagnostics.slsSubstraction() ;
    diagnostics.slsScalarMultiplication() ;
    diagnostics.slsNominalMultiplication() ;
    diagnostics.slsOptimisedMultiplication() ;
    diagnostics.slsResidual() ;
    diagnostics.slsOptimisedResidual() ;

    // Runs dense linear-system diagnostics, including algebraic operations, nominal and BLAS multiplication, and residual evaluation.
    diagnostics.dlsAddition() ;
    diagnostics.dlsSubstraction() ;
    diagnostics.dlsScalarMultiplication() ;
    diagnostics.dlsNominalMultiplication() ;
    diagnostics.dlsBLASMultiplication() ;
    diagnostics.dlsResidual() ;
    diagnostics.dlsOptimisedResidual() ;

    // Writes representative sample matrix-vector output for visualisation.
    diagnostics.SampleMatrixVectorOutput() ;
}

