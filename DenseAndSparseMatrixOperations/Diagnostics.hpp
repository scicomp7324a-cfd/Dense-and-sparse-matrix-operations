# include <iostream>
# include <ostream>
# include <istream>
# include <chrono>
# include <cstdint>
# include <random>

# include "FaceAddressedMesh2D.hpp"
# include "MeshFileReader2D.hpp"
# include "LinearVector.hpp"

# include "DenseLinearSystem.hpp"
# include "DenseMatrix.hpp"

# include "SparseAddress.hpp"
# include "SparseLinearSystem.hpp"
# include "SparseMatrix.hpp"
# include "PathConfig.hpp"

using namespace std ;

class Diagnostics{

    private:

    const vector<int> n_ = {10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 15000, 20000, 30000, 40000, 50000} ;
    const vector<int> ndense_ = {10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 15000, 20000, 30000, 40000} ;
    PathConfig paths ;
    
    public:

    // vector diagnostics

    /*
        Vector addition diagnostics
    */

    void vectorAdditionDiagnostics(){
        
        ofstream out(paths.outputDir() / "outputVectorAddition.csv") ;
        out << "Vector Size,Time Taken" << endl ;
        for(int dimension: n_){

            vector<mt19937> rngs = {make_rng(dimension + 46), 
                                    make_rng(dimension + 100), 
                                    make_rng(dimension + 135),
                                    make_rng(dimension + 56),
                                    make_rng(dimension + 87)} ;
            vector<double> timeForEachCombination ;

            for(size_t i = 0; i < rngs.size(); ++i){
                for(size_t j = i+1; j < rngs.size(); ++j){

                    LinearVector X1 = setupVector(dimension,rngs[i]);
                    LinearVector X2 = setupVector(dimension,rngs[j]) ;

                    size_t repeats = 2 ;
                    double time_sum = 0.0 ;
                    for(size_t i=0; i < repeats; ++i){

                        using clock = std::chrono::steady_clock;
                        auto t0 = clock::now() ;
                        X1+X2 ;
                        auto t1 = clock::now() ;

                        std::chrono::duration<double> time_taken = t1-t0 ;
                        time_sum += time_taken.count() ;
                    }

                    timeForEachCombination.push_back((time_sum/repeats)) ;
                }
            }

            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    
    }

    /*
        Vector substraction diagnostics
    */

    void vectorSubstractionDiagnostics(){
        
        ofstream out(paths.outputDir() / "outputVectorSubstraction.csv") ;
        out << "Vector Size,Time Taken" << endl ;

        for(int dimension: n_){

            vector<mt19937> rngs = {make_rng(dimension + 46), 
                                    make_rng(dimension + 100), 
                                    make_rng(dimension + 135),
                                    make_rng(dimension + 56),
                                    make_rng(dimension + 87)} ;
            vector<double> timeForEachCombination ;

            for(size_t i = 0; i < rngs.size(); ++i){
                for(size_t j = i+1; j < rngs.size(); ++j){

                    LinearVector X1 = setupVector(dimension,rngs[i]);
                    LinearVector X2 = setupVector(dimension,rngs[j]) ;

                    size_t repeats = 2 ;
                    double time_sum = 0.0 ;
                    for(size_t i=0; i < repeats; ++i){

                        using clock = std::chrono::steady_clock;
                        auto t0 = clock::now() ;
                        X1-X2 ;
                        auto t1 = clock::now() ;

                        std::chrono::duration<double> time_taken = t1-t0 ;
                        time_sum += time_taken.count() ;
                    }

                    timeForEachCombination.push_back((time_sum/repeats)) ;
                }
            }

            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    
    }


    /*
        Vector multiplication by scalar
    */

    void vectorScalarMultiplicationDiagnostics(){

        ofstream out(paths.outputDir() / "outputVectorScalarMultiplication.csv") ;
        out << "Vector Size,Time Taken" << endl ;

        for(int dimension: n_){

            vector<mt19937> rngs = {make_rng(dimension + 46), 
                                    make_rng(dimension + 100), 
                                    make_rng(dimension + 135),
                                    make_rng(dimension + 56),
                                    make_rng(dimension + 87)} ;
            vector<double> timeForEachCombination ;

            for(size_t i = 0; i < rngs.size(); ++i){
                for(size_t j = 0; j < rngs.size(); ++j){

                    LinearVector X1 = setupVector(dimension,rngs[i]);
                    double s = rand_uniform(rngs[j], -999.999, 999.999) ;

                    size_t repeats = 2 ;
                    double time_sum = 0.0 ;
                    for(size_t i=0; i < repeats; ++i){

                        using clock = std::chrono::steady_clock;
                        auto t0 = clock::now() ;
                        X1*s ;
                        auto t1 = clock::now() ;

                        std::chrono::duration<double> time_taken = t1-t0 ;
                        time_sum += time_taken.count() ;
                    }

                    timeForEachCombination.push_back((time_sum/repeats)) ;
                }
            }

            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    
    }

    /*
        Vector one norm
    */

    void vectorOneNormDiagnostics(){

        ofstream out(paths.outputDir() / "outputVectorOneNorm.csv") ;
        out << "Vector Size,Time Taken" << endl ;

        for(int dimension: n_){

            vector<mt19937> rngs = {make_rng(dimension + 46), 
                                    make_rng(dimension + 100), 
                                    make_rng(dimension + 135),
                                    make_rng(dimension + 56),
                                    make_rng(dimension + 87)} ;
            vector<double> timeForEachCombination ;

            for(size_t i = 0; i < rngs.size(); ++i){

                LinearVector X1 = setupVector(dimension,rngs[i]);
                size_t repeats = 2 ;
                double time_sum = 0.0 ;
                for(size_t i=0; i < repeats; ++i){
                    using clock = std::chrono::steady_clock;
                    auto t0 = clock::now() ;
                    X1.oneNorm() ;
                    auto t1 = clock::now() ;

                    std::chrono::duration<double> time_taken = t1-t0 ;
                    time_sum += time_taken.count() ;
                }
                timeForEachCombination.push_back((time_sum/repeats)) ;

            }
            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    
    }

    /*
        Vector two norm
    */

    void vectorTwoNormDiagnostics(){

        ofstream out(paths.outputDir() / "outputVectorTwoNorm.csv") ;
        out << "Vector Size,Time Taken" << endl ;

        for(int dimension: n_){

            vector<mt19937> rngs = {make_rng(dimension + 46), 
                                    make_rng(dimension + 100), 
                                    make_rng(dimension + 135),
                                    make_rng(dimension + 56),
                                    make_rng(dimension + 87)} ;
            vector<double> timeForEachCombination ;

            for(size_t i = 0; i < rngs.size(); ++i){

                LinearVector X1 = setupVector(dimension,rngs[i]);
                size_t repeats = 2 ;
                double time_sum = 0.0 ;
                for(size_t i=0; i < repeats; ++i){
                    using clock = std::chrono::steady_clock;
                    auto t0 = clock::now() ;
                    X1.twoNorm() ;
                    auto t1 = clock::now() ;

                    std::chrono::duration<double> time_taken = t1-t0 ;
                    time_sum += time_taken.count() ;
                }
                timeForEachCombination.push_back((time_sum/repeats)) ;

            }
            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    
    }

    /*
        Vector infinity norm
    */

    void vectorInfinityNormDiagnostics(){

        ofstream out(paths.outputDir() / "outputVectorInfinityNorm.csv") ;
        out << "Vector Size,Time Taken" << endl ;

        for(int dimension: n_){

            vector<mt19937> rngs = {make_rng(dimension + 46), 
                                    make_rng(dimension + 100), 
                                    make_rng(dimension + 135),
                                    make_rng(dimension + 56),
                                    make_rng(dimension + 87)} ;
            vector<double> timeForEachCombination ;

            for(size_t i = 0; i < rngs.size(); ++i){

                LinearVector X1 = setupVector(dimension,rngs[i]);
                size_t repeats = 2 ;
                double time_sum = 0.0 ;
                for(size_t i=0; i < repeats; ++i){
                    using clock = std::chrono::steady_clock;
                    auto t0 = clock::now() ;
                    X1.infinityNorm() ;
                    auto t1 = clock::now() ;

                    std::chrono::duration<double> time_taken = t1-t0 ;
                    time_sum += time_taken.count() ;
                }
                timeForEachCombination.push_back((time_sum/repeats)) ;

            }
            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    
    }

    // Sparse Matrix diagnostics

    /*
        Sparse matrix addition
    */

    void sparseMatrixAdditionDiagnostics(){

        ofstream out(paths.outputDir() / "outputSparseMatrixAddition.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: n_){

            SparseAddress addr_(paths.getMeshFilenames(dimension)) ;
            SparseMatrix A1(addr_) ;
            SparseMatrix A2(addr_) ;
            vector<double> timeForEachCombination ;
            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;

            for(size_t i = 0; i < mtxFilenames.size(); ++i){
                for(size_t j = i+1; j < mtxFilenames.size(); ++j){

                    ifstream in1(mtxFilenames[i]) ;
                    ifstream in2(mtxFilenames[j]) ;
                    A1.readMatrix(in1) ;
                    A2.readMatrix(in2) ;

                    size_t repeats = 2 ;
                    double time_sum = 0.0 ;
                    for(size_t i=0; i < repeats; ++i){

                        auto rng1 = make_rng(dimension + 46)  ;
                        auto rng2 = make_rng(dimension + 100)  ;
                        double s1 = rand_uniform(rng1, -999.999, 999.999) ;
                        double s2 = rand_uniform(rng2, -999.999, 999.999) ;
                        SparseMatrix M1 = A1*s1 ;
                        SparseMatrix M2 = A2*s2 ;

                        using clock = std::chrono::steady_clock;
                        auto t0 = clock::now() ;
                        M1 + M2 ;
                        auto t1 = clock::now() ;

                        std::chrono::duration<double> time_taken = t1-t0 ;
                        time_sum += time_taken.count() ;

                    }
                    timeForEachCombination.push_back((time_sum/repeats)) ;
                }
            }

            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    }

    /*
        Sparse matrix substraction
    */

    void sparseMatrixSubstractionDiagnostics(){

        ofstream out(paths.outputDir() / "outputSparseMatrixSubstraction.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: n_){

            SparseAddress addr_(paths.getMeshFilenames(dimension)) ;
            SparseMatrix A1(addr_) ;
            SparseMatrix A2(addr_) ;
            vector<double> timeForEachCombination ;
            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;

            for(size_t i = 0; i < mtxFilenames.size(); ++i){
                for(size_t j = i+1; j < mtxFilenames.size(); ++j){

                    ifstream in1(mtxFilenames[i]) ;
                    ifstream in2(mtxFilenames[j]) ;
                    A1.readMatrix(in1) ;
                    A2.readMatrix(in2) ;

                    size_t repeats = 2 ;
                    double time_sum = 0.0 ;
                    for(size_t i=0; i < repeats; ++i){

                        auto rng1 = make_rng(dimension + 46)  ;
                        auto rng2 = make_rng(dimension + 100)  ;
                        double s1 = rand_uniform(rng1, -999.999, 999.999) ;
                        double s2 = rand_uniform(rng2, -999.999, 999.999) ;
                        SparseMatrix M1 = A1*s1 ;
                        SparseMatrix M2 = A2*s2 ;

                        using clock = std::chrono::steady_clock;
                        auto t0 = clock::now() ;
                        M1 - M2 ;
                        auto t1 = clock::now() ;

                        std::chrono::duration<double> time_taken = t1-t0 ;
                        time_sum += time_taken.count() ;

                    }
                    timeForEachCombination.push_back((time_sum/repeats)) ;
                }
            }

            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    }

    /*
        Sparse matrix scalar multiplication
    */

    void sparseMatrixScalarMultiplicationDiagnostics(){

        ofstream out(paths.outputDir() / "outputSparseMatrixScalarMultiplication.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: n_){

            SparseAddress addr_(paths.getMeshFilenames(dimension)) ;
            SparseMatrix A1(addr_) ;
            vector<double> timeForEachCombination ;
            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;

            for(size_t i = 0; i < mtxFilenames.size(); ++i){

                ifstream in1(mtxFilenames[i]) ;
                A1.readMatrix(in1) ;
                
                size_t repeats = 2 ;
                double time_sum = 0.0 ;
                for(size_t i=0; i < repeats; ++i){

                    auto rng1 = make_rng(dimension + 46)  ;
                    double s1 = rand_uniform(rng1, -999.999, 999.999) ;
                    
                    using clock = std::chrono::steady_clock;
                    auto t0 = clock::now() ;
                    A1*s1 ;                
                    auto t1 = clock::now() ;

                    std::chrono::duration<double> time_taken = t1-t0 ;
                    time_sum += time_taken.count() ;

                }
                timeForEachCombination.push_back((time_sum/repeats)) ;
                
            }
            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    }


    // Dense Matrix Diagnostics

    /*
        Dense Matrix Addition
    */

    void denseMatrixAdditionDiagnostics(){

        ofstream out(paths.outputDir() / "outputDenseMatrixAddition.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: ndense_ ){

            DenseMatrix A1(dimension) ;
            DenseMatrix A2(dimension) ;

            vector<double> timeForEachCombination ;
            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;

            for(size_t i = 0; i < mtxFilenames.size(); ++i){
                for(size_t j = i+1; j < mtxFilenames.size(); ++j){

                    ifstream in1(mtxFilenames[i]) ;
                    ifstream in2(mtxFilenames[j]) ;
                    A1.readMatrix(in1) ;
                    A2.readMatrix(in2) ;

                    size_t repeats = 2 ;
                    double time_sum = 0.0 ;
                    for(size_t i=0; i < repeats; ++i){

                        auto rng1 = make_rng(dimension + 46)  ;
                        auto rng2 = make_rng(dimension + 100)  ;
                        double s1 = rand_uniform(rng1, -999.999, 999.999) ;
                        double s2 = rand_uniform(rng2, -999.999, 999.999) ;
                        DenseMatrix M1 = A1*s1 ;
                        DenseMatrix M2 = A2*s2 ;

                        using clock = std::chrono::steady_clock;
                        auto t0 = clock::now() ;
                        M1 + M2 ;
                        auto t1 = clock::now() ;

                        std::chrono::duration<double> time_taken = t1-t0 ;
                        time_sum += time_taken.count() ;

                    }
                    timeForEachCombination.push_back((time_sum/repeats)) ;
                }
            }

            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    }

    /*
        Dense Matrix Substraction
    */

    void denseMatrixSubstractionDiagnostics(){

        ofstream out(paths.outputDir() / "outputDenseMatrixSubstraction.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: ndense_ ){

            DenseMatrix A1(dimension) ;
            DenseMatrix A2(dimension) ;

            vector<double> timeForEachCombination ;
            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;

            for(size_t i = 0; i < mtxFilenames.size(); ++i){
                for(size_t j = i+1; j < mtxFilenames.size(); ++j){

                    ifstream in1(mtxFilenames[i]) ;
                    ifstream in2(mtxFilenames[j]) ;
                    A1.readMatrix(in1) ;
                    A2.readMatrix(in2) ;

                    size_t repeats = 2 ;
                    double time_sum = 0.0 ;
                    for(size_t i=0; i < repeats; ++i){

                        auto rng1 = make_rng(dimension + 46)  ;
                        auto rng2 = make_rng(dimension + 100)  ;
                        double s1 = rand_uniform(rng1, -999.999, 999.999) ;
                        double s2 = rand_uniform(rng2, -999.999, 999.999) ;
                        DenseMatrix M1 = A1*s1 ;
                        DenseMatrix M2 = A2*s2 ;

                        using clock = std::chrono::steady_clock;
                        auto t0 = clock::now() ;
                        M1 - M2 ;
                        auto t1 = clock::now() ;

                        std::chrono::duration<double> time_taken = t1-t0 ;
                        time_sum += time_taken.count() ;

                    }
                    timeForEachCombination.push_back((time_sum/repeats)) ;
                }
            }

            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    }

    /*
        Dense Matrix Scalar Multiplication
    */

    void denseMatrixScalarMultiplicationDiagnostics(){

        ofstream out(paths.outputDir() / "outputDenseMatrixScalarMultiplication.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: ndense_ ){

            DenseMatrix A1(dimension) ;
            vector<double> timeForEachCombination ;
            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;

            for(size_t i = 0; i < mtxFilenames.size(); ++i){

                ifstream in1(mtxFilenames[i]) ;
                A1.readMatrix(in1) ;
                
                size_t repeats = 2 ;
                double time_sum = 0.0 ;
                for(size_t i=0; i < repeats; ++i){

                    auto rng1 = make_rng(dimension + 46)  ;
                    double s1 = rand_uniform(rng1, -999.999, 999.999) ;
                    
                    using clock = std::chrono::steady_clock;
                    auto t0 = clock::now() ;
                    A1*s1 ;                
                    auto t1 = clock::now() ;

                    std::chrono::duration<double> time_taken = t1-t0 ;
                    time_sum += time_taken.count() ;

                }
                timeForEachCombination.push_back((time_sum/repeats)) ;
                
            }
            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    }

    // Sparse linear system diagnostics

    /*
        Sparse linear system addition
    */

    void slsAddition(){

        ofstream out(paths.outputDir() / "outputSLSAddition.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: n_){

            SparseAddress addr_(paths.getMeshFilenames(dimension)) ;
            SparseLinearSystem sys1(addr_) ;
            SparseLinearSystem sys2(addr_) ;

            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;
            vector<double> timeForEachCombination ;
           
            for(size_t i = 0; i < mtxFilenames.size(); ++i){
                for(size_t j = i+1; j < mtxFilenames.size(); ++j){

                    ifstream in1(mtxFilenames[i]) ;
                    ifstream in2(mtxFilenames[j]) ;
                    sys1.A().readMatrix(in1) ;
                    sys2.A().readMatrix(in2) ;

                    size_t repeats = 2 ;
                    double time_sum = 0.0 ;
                    for(size_t i=0; i < repeats; ++i){

                        auto rng1 = make_rng(dimension + 46) ;
                        auto rng2 = make_rng(dimension + 112) ;
                        auto rng3 = make_rng(dimension + 54) ;
                        auto rng4 = make_rng(dimension + 77) ;

                        sys1.X() = setupVector(dimension, rng1) ;
                        sys1.b() = setupVector(dimension, rng2) ;
                        sys2.X() = setupVector(dimension, rng3) ;
                        sys2.b() = setupVector(dimension, rng4) ;
                        
                        using clock = std::chrono::steady_clock;
                        auto t0 = clock::now() ;
                        sys1 + sys2 ;
                        auto t1 = clock::now() ;

                        std::chrono::duration<double> time_taken = t1-t0 ;
                        time_sum += time_taken.count() ;

                    }
                    timeForEachCombination.push_back((time_sum/repeats)) ;
                }
            }

            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    
    }

    /*
        Sparse linear system substraction
    */

    void slsSubstraction(){
        ofstream out(paths.outputDir() / "outputSLSSubstraction.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: n_){

            SparseAddress addr_(paths.getMeshFilenames(dimension)) ;
            SparseLinearSystem sys1(addr_) ;
            SparseLinearSystem sys2(addr_) ;

            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;
            vector<double> timeForEachCombination ;
           
            for(size_t i = 0; i < mtxFilenames.size(); ++i){
                for(size_t j = i+1; j < mtxFilenames.size(); ++j){

                    ifstream in1(mtxFilenames[i]) ;
                    ifstream in2(mtxFilenames[j]) ;
                    sys1.A().readMatrix(in1) ;
                    sys2.A().readMatrix(in2) ;

                    size_t repeats = 2 ;
                    double time_sum = 0.0 ;
                    for(size_t i=0; i < repeats; ++i){

                        auto rng1 = make_rng(dimension + 46) ;
                        auto rng2 = make_rng(dimension + 112) ;
                        auto rng3 = make_rng(dimension + 54) ;
                        auto rng4 = make_rng(dimension + 77) ;

                        sys1.X() = setupVector(dimension, rng1) ;
                        sys1.b() = setupVector(dimension, rng2) ;
                        sys2.X() = setupVector(dimension, rng3) ;
                        sys2.b() = setupVector(dimension, rng4) ;
                        
                        using clock = std::chrono::steady_clock;
                        auto t0 = clock::now() ;
                        sys1 - sys2 ;
                        auto t1 = clock::now() ;

                        std::chrono::duration<double> time_taken = t1-t0 ;
                        time_sum += time_taken.count() ;

                    }
                    timeForEachCombination.push_back((time_sum/repeats)) ;
                }
            }

            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    
    }

    /*
        Sparse linear system scalar multiplication
    */

    void slsScalarMultiplication(){

        ofstream out(paths.outputDir() / "outputSLSScalarMultiplication.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: n_){

            SparseAddress addr_(paths.getMeshFilenames(dimension)) ;
            SparseLinearSystem sys1(addr_) ;

            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;
            vector<double> timeForEachCombination ;
           
            for(size_t i = 0; i < mtxFilenames.size(); ++i){

                ifstream in1(mtxFilenames[i]) ;
                sys1.A().readMatrix(in1) ;

                size_t repeats = 2 ;
                double time_sum = 0.0 ;
                for(size_t i=0; i < repeats; ++i){

                    auto rng1 = make_rng(dimension + 46) ;
                    auto rng2 = make_rng(dimension + 112) ;
                    auto rng3 = make_rng(dimension + 54) ;

                    sys1.X() = setupVector(dimension, rng1) ;
                    sys1.b() = setupVector(dimension, rng2) ;
                    double s1 = rand_uniform(rng3, -999.999, 999.999) ;
                    
                    using clock = std::chrono::steady_clock;
                    auto t0 = clock::now() ;
                    sys1*s1 ;
                    auto t1 = clock::now() ;

                    std::chrono::duration<double> time_taken = t1-t0 ;
                    time_sum += time_taken.count() ;

                }
                timeForEachCombination.push_back((time_sum/repeats)) ;
            }
            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    }

    /*
        Sparse linear system multiplication (nominal)
    */

    void slsNominalMultiplication(){

        ofstream out(paths.outputDir() / "outputSLSNominalMultiplication.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: n_){

            SparseAddress addr_(paths.getMeshFilenames(dimension)) ;
            SparseLinearSystem sys1(addr_) ;

            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;
            vector<double> timeForEachCombination ;

            for(size_t i=0; i < mtxFilenames.size(); ++i){

                ifstream in1(mtxFilenames[i]) ;
                sys1.A().readMatrix(in1) ;

                size_t repeats = 2 ;
                double time_sum = 0.0 ;
                for(size_t i=0; i < repeats; ++i){

                    auto rng1 = make_rng(dimension + 46) ;
                    auto rng2 = make_rng(dimension + 112) ;

                    sys1.X() = setupVector(dimension, rng1) ;
                    sys1.b() = setupVector(dimension, rng2) ;
                    
                    using clock = std::chrono::steady_clock;
                    auto t0 = clock::now() ;
                    sys1.A().multiplyNominal(sys1.X()) ;
                    auto t1 = clock::now() ;

                    std::chrono::duration<double> time_taken = t1-t0 ;
                    time_sum += time_taken.count() ;

                }
                timeForEachCombination.push_back((time_sum/repeats)) ;
            }
            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    }

    /*
        Sparse linear system multiplication (optimised)
    */

    void slsOptimisedMultiplication(){

        ofstream out(paths.outputDir() / "outputSLSOptimisedMultiplication.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: n_){

            SparseAddress addr_(paths.getMeshFilenames(dimension)) ;
            SparseLinearSystem sys1(addr_) ;

            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;
            vector<double> timeForEachCombination ;

            for(size_t i=0; i < mtxFilenames.size(); ++i){

                ifstream in1(mtxFilenames[i]) ;
                sys1.A().readMatrix(in1) ;

                size_t repeats = 2 ;
                double time_sum = 0.0 ;
                for(size_t i=0; i < repeats; ++i){

                    auto rng1 = make_rng(dimension + 46) ;
                    auto rng2 = make_rng(dimension + 112) ;

                    sys1.X() = setupVector(dimension, rng1) ;
                    sys1.b() = setupVector(dimension, rng2) ;
                    
                    using clock = std::chrono::steady_clock;
                    auto t0 = clock::now() ;
                    sys1.A().multiplyOptimised(sys1.X()) ;
                    auto t1 = clock::now() ;

                    std::chrono::duration<double> time_taken = t1-t0 ;
                    time_sum += time_taken.count() ;

                }
                timeForEachCombination.push_back((time_sum/repeats)) ;
            }
            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    }

    /*
        Sparse matrix residual
    */

    void slsResidual(){

        ofstream out(paths.outputDir() / "outputSLSResidual.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: n_){

            SparseAddress addr_(paths.getMeshFilenames(dimension)) ;
            SparseLinearSystem sys1(addr_) ;

            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;
            vector<double> timeForEachCombination ;

            for(size_t i=0; i < mtxFilenames.size(); ++i){

                ifstream in1(mtxFilenames[i]) ;
                sys1.A().readMatrix(in1) ;

                size_t repeats = 2 ;
                double time_sum = 0.0 ;
                for(size_t i=0; i < repeats; ++i){

                    auto rng1 = make_rng(dimension + 46) ;
                    auto rng2 = make_rng(dimension + 112) ;

                    sys1.X() = setupVector(dimension, rng1) ;
                    sys1.b() = setupVector(dimension, rng2) ;
                    
                    using clock = std::chrono::steady_clock;
                    auto t0 = clock::now() ;
                    sys1.calculateResidual() ;
                    auto t1 = clock::now() ;

                    std::chrono::duration<double> time_taken = t1-t0 ;
                    time_sum += time_taken.count() ;

                }
                timeForEachCombination.push_back((time_sum/repeats)) ;
            }
            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }

    }

    /*
        Sparse matrix optimised residual 
    */

    void slsOptimisedResidual(){

        ofstream out(paths.outputDir() / "outputSLSOptimisedResidual.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: n_){

            SparseAddress addr_(paths.getMeshFilenames(dimension)) ;
            SparseLinearSystem sys1(addr_) ;

            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;
            vector<double> timeForEachCombination ;

            for(size_t i=0; i < mtxFilenames.size(); ++i){

                ifstream in1(mtxFilenames[i]) ;
                sys1.A().readMatrix(in1) ;

                size_t repeats = 2 ;
                double time_sum = 0.0 ;
                for(size_t i=0; i < repeats; ++i){

                    auto rng1 = make_rng(dimension + 46) ;
                    auto rng2 = make_rng(dimension + 112) ;

                    sys1.X() = setupVector(dimension, rng1) ;
                    sys1.b() = setupVector(dimension, rng2) ;
                    
                    using clock = std::chrono::steady_clock;
                    auto t0 = clock::now() ;
                    sys1.calculateResidualOptimised() ;
                    auto t1 = clock::now() ;

                    std::chrono::duration<double> time_taken = t1-t0 ;
                    time_sum += time_taken.count() ;

                }
                timeForEachCombination.push_back((time_sum/repeats)) ;
            }
            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
        
    }

    // Dense linear system diagnostics

    /*
        Dense linear system addition
    */

    void dlsAddition(){

        ofstream out(paths.outputDir() / "outputDLSAddition.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: ndense_ ){

            DenseLinearSystem sys1(dimension) ;
            DenseLinearSystem sys2(dimension) ;

            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;
            vector<double> timeForEachCombination ;
           
            for(size_t i = 0; i < mtxFilenames.size(); ++i){
                for(size_t j = i+1; j < mtxFilenames.size(); ++j){

                    ifstream in1(mtxFilenames[i]) ;
                    ifstream in2(mtxFilenames[j]) ;
                    sys1.A().readMatrix(in1) ;
                    sys2.A().readMatrix(in2) ;

                    size_t repeats = 2 ;
                    double time_sum = 0.0 ;
                    for(size_t i=0; i < repeats; ++i){

                        auto rng1 = make_rng(dimension + 46) ;
                        auto rng2 = make_rng(dimension + 112) ;
                        auto rng3 = make_rng(dimension + 54) ;
                        auto rng4 = make_rng(dimension + 77) ;

                        sys1.X() = setupVector(dimension, rng1) ;
                        sys1.b() = setupVector(dimension, rng2) ;
                        sys2.X() = setupVector(dimension, rng3) ;
                        sys2.b() = setupVector(dimension, rng4) ;
                        
                        using clock = std::chrono::steady_clock;
                        auto t0 = clock::now() ;
                        sys1 + sys2 ;
                        auto t1 = clock::now() ;

                        std::chrono::duration<double> time_taken = t1-t0 ;
                        time_sum += time_taken.count() ;

                    }
                    timeForEachCombination.push_back((time_sum/repeats)) ;
                }
            }
            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    
    }

    /*
        Dense linear system substraction
    */

    void dlsSubstraction(){

        ofstream out(paths.outputDir() / "outputDLSSubstraction.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: ndense_ ){

            DenseLinearSystem sys1(dimension) ;
            DenseLinearSystem sys2(dimension) ;

            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;
            vector<double> timeForEachCombination ;
           
            for(size_t i = 0; i < mtxFilenames.size(); ++i){
                for(size_t j = i+1; j < mtxFilenames.size(); ++j){

                    ifstream in1(mtxFilenames[i]) ;
                    ifstream in2(mtxFilenames[j]) ;
                    sys1.A().readMatrix(in1) ;
                    sys2.A().readMatrix(in2) ;

                    size_t repeats = 2 ;
                    double time_sum = 0.0 ;
                    for(size_t i=0; i < repeats; ++i){

                        auto rng1 = make_rng(dimension + 46) ;
                        auto rng2 = make_rng(dimension + 112) ;
                        auto rng3 = make_rng(dimension + 54) ;
                        auto rng4 = make_rng(dimension + 77) ;

                        sys1.X() = setupVector(dimension, rng1) ;
                        sys1.b() = setupVector(dimension, rng2) ;
                        sys2.X() = setupVector(dimension, rng3) ;
                        sys2.b() = setupVector(dimension, rng4) ;
                        
                        using clock = std::chrono::steady_clock;
                        auto t0 = clock::now() ;
                        sys1 - sys2 ;
                        auto t1 = clock::now() ;

                        std::chrono::duration<double> time_taken = t1-t0 ;
                        time_sum += time_taken.count() ;

                    }
                    timeForEachCombination.push_back((time_sum/repeats)) ;
                }
            }

            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    
    }


    /*
        Dense linear system scalar multiplication 
    */

    void dlsScalarMultiplication(){

        ofstream out(paths.outputDir() / "outputDLSScalarMultiplication.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: ndense_ ){

            DenseLinearSystem sys1(dimension) ;
            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;
            vector<double> timeForEachCombination ;

            for(size_t i = 0; i < mtxFilenames.size(); ++i){

                ifstream in1(mtxFilenames[i]) ;
                sys1.A().readMatrix(in1) ;

                size_t repeats = 2 ;
                double time_sum = 0.0 ;
                for(size_t i=0; i < repeats; ++i){

                    auto rng1 = make_rng(dimension + 46) ;
                    auto rng2 = make_rng(dimension + 112) ;
                    auto rng3 = make_rng(dimension + 54) ;

                    sys1.X() = setupVector(dimension, rng1) ;
                    sys1.b() = setupVector(dimension, rng2) ;
                    double s1 = rand_uniform(rng3, -999.999, 999.999) ;
                    
                    using clock = std::chrono::steady_clock;
                    auto t0 = clock::now() ;
                    sys1*s1 ;
                    auto t1 = clock::now() ;

                    std::chrono::duration<double> time_taken = t1-t0 ;
                    time_sum += time_taken.count() ;
                }
                timeForEachCombination.push_back((time_sum/repeats)) ;
            }
            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    }

    /*
        Dense linear system nominal multiplication
    */

    void dlsNominalMultiplication(){

        ofstream out(paths.outputDir() / "outputDLSNominalMultiplication.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: ndense_ ){

            DenseLinearSystem sys1(dimension) ;

            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;
            vector<double> timeForEachCombination ;

            for(size_t i=0; i < mtxFilenames.size(); ++i){

                ifstream in1(mtxFilenames[i]) ;
                sys1.A().readMatrix(in1) ;

                size_t repeats = 2 ;
                double time_sum = 0.0 ;
                for(size_t i=0; i < repeats; ++i){

                    auto rng1 = make_rng(dimension + 46) ;
                    auto rng2 = make_rng(dimension + 112) ;

                    sys1.X() = setupVector(dimension, rng1) ;
                    sys1.b() = setupVector(dimension, rng2) ;
                    
                    using clock = std::chrono::steady_clock;
                    auto t0 = clock::now() ;
                    sys1.A().multiplyNominal(sys1.X()) ;
                    auto t1 = clock::now() ;

                    std::chrono::duration<double> time_taken = t1-t0 ;
                    time_sum += time_taken.count() ;

                }
                timeForEachCombination.push_back((time_sum/repeats)) ;
            }
            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    
    }

    /*
        Dense linear system optimised multiplication
    */

    void dlsBLASMultiplication(){

        ofstream out(paths.outputDir() / "outputDLSBLASMultiplication.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: ndense_ ){

            DenseLinearSystem sys1(dimension) ;

            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;
            vector<double> timeForEachCombination ;

            for(size_t i=0; i < mtxFilenames.size(); ++i){

                ifstream in1(mtxFilenames[i]) ;
                sys1.A().readMatrix(in1) ;

                size_t repeats = 2 ;
                double time_sum = 0.0 ;
                for(size_t i=0; i < repeats; ++i){

                    auto rng1 = make_rng(dimension + 46) ;
                    auto rng2 = make_rng(dimension + 112) ;

                    sys1.X() = setupVector(dimension, rng1) ;
                    sys1.b() = setupVector(dimension, rng2) ;
                    
                    using clock = std::chrono::steady_clock;
                    auto t0 = clock::now() ;
                    sys1.A().multiplyBLAS(sys1.X()) ;
                    auto t1 = clock::now() ;

                    std::chrono::duration<double> time_taken = t1-t0 ;
                    time_sum += time_taken.count() ;

                }
                timeForEachCombination.push_back((time_sum/repeats)) ;
            }
            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
    
    }

    /*
        Dense linear system residual
    */

    void dlsResidual(){

        ofstream out(paths.outputDir() / "outputDLSResidual.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: ndense_ ){

            DenseLinearSystem sys1(dimension) ;

            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;
            vector<double> timeForEachCombination ;

            for(size_t i=0; i < mtxFilenames.size(); ++i){

                ifstream in1(mtxFilenames[i]) ;
                sys1.A().readMatrix(in1) ;

                size_t repeats = 2 ;
                double time_sum = 0.0 ;
                for(size_t i=0; i < repeats; ++i){

                    auto rng1 = make_rng(dimension + 46) ;
                    auto rng2 = make_rng(dimension + 112) ;

                    sys1.X() = setupVector(dimension, rng1) ;
                    sys1.b() = setupVector(dimension, rng2) ;
                    
                    using clock = std::chrono::steady_clock;
                    auto t0 = clock::now() ;
                    sys1.calculateResidual() ;
                    auto t1 = clock::now() ;

                    std::chrono::duration<double> time_taken = t1-t0 ;
                    time_sum += time_taken.count() ;

                }
                timeForEachCombination.push_back((time_sum/repeats)) ;
            }
            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }

    }

    /*
        Dense linear system optimised residual
    */

    void dlsOptimisedResidual(){

        ofstream out(paths.outputDir() / "outputDLSOptimisedResidual.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: ndense_ ){

            DenseLinearSystem sys1(dimension) ;

            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;
            vector<double> timeForEachCombination ;

            for(size_t i=0; i < mtxFilenames.size(); ++i){

                ifstream in1(mtxFilenames[i]) ;
                sys1.A().readMatrix(in1) ;

                size_t repeats = 2 ;
                double time_sum = 0.0 ;
                for(size_t i=0; i < repeats; ++i){

                    auto rng1 = make_rng(dimension + 46) ;
                    auto rng2 = make_rng(dimension + 112) ;

                    sys1.X() = setupVector(dimension, rng1) ;
                    sys1.b() = setupVector(dimension, rng2) ;
                    
                    using clock = std::chrono::steady_clock;
                    auto t0 = clock::now() ;
                    sys1.calculateResidualOptimised() ;
                    auto t1 = clock::now() ;

                    std::chrono::duration<double> time_taken = t1-t0 ;
                    time_sum += time_taken.count() ;

                }
                timeForEachCombination.push_back((time_sum/repeats)) ;
            }
            double time = average(timeForEachCombination) ;
            out << dimension << "," << time << endl ;
        }
        
    }


    /*
        Sample matrix and vector visualisation script.
    */

    void SampleMatrixVectorOutput(){
        PathConfig paths ;
        vector<int> ndense_ = {10, 20, 50, 100, 200} ;
        ofstream out(paths.outputDir() / "outputDLSBLASMultiplication.csv") ;
        out << "Matrix Size,Time Taken" << endl ;

        for(int dimension: ndense_ ){

            DenseLinearSystem sys1(dimension) ;

            vector<string> mtxFilenames = paths.getMTXFilenames(dimension) ;
            vector<double> timeForEachCombination ;
            ifstream in1(mtxFilenames[0]) ;
            sys1.A().readMatrix(in1) ;

            auto rng1 = make_rng(dimension + 46) ;
            auto rng2 = make_rng(dimension + 112) ;

            sys1.X() = setupVector(dimension, rng1) ;
            sys1.b() = setupVector(dimension, rng2) ;

            fs::path newpathmatrix = paths.outputDir() / ("_matrix_" + to_string(dimension) + ".csv") ;
            fs::path newpathvector = paths.outputDir() / ("_vector_" + to_string(dimension) + ".csv") ;

            sys1.writeMatrixCSV(sys1.A(), newpathmatrix) ;
            sys1.writeVectorCSV(sys1.X(), newpathvector) ;
        }
    }

   
    LinearVector setupVector(int size, std::mt19937& rng){
        LinearVector x(size) ;
        for(int i=0; i<size; ++i){
            x(i) = rand_uniform(rng, -999999, 999999) ;
        }
        return x ;
    }

    mt19937 make_rng(unsigned seed = 12345u) {
        return mt19937(seed);
    }

    double rand_uniform(mt19937& rng, double lo = -1.0, double hi = 1.0) {
        uniform_real_distribution<double> dist(lo, hi);
        return dist(rng);
    }

    double average(vector<double> time){
        if(time.empty()) return 0.0;
        double sum = 0.0 ;
        for(size_t i=0; i<time.size(); ++i){
            sum += time[i] ;
        }
        return sum/time.size() ;
    }

} ;
