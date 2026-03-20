# pragma once

# include <iostream>
# include <vector>
# include <filesystem>
# include <stdexcept>

# include "SparseMatrix.hpp"
# include "LinearVector.hpp"

using namespace std ;
namespace fs = filesystem ;

/* *******************************************************************
    SparseLinearSystem.hpp - This class contains all the algebraic 
    operations for sparse linear system and residual calculation  
    functions.
******************************************************************* */

class SparseLinearSystem{

    private:

        SparseMatrix A_ ; 
        int n_ ;

        LinearVector X_ ;
        LinearVector b_ ;
        LinearVector r_ ;

    public:

        // Constructor builds the sparse linear system using sparse addressing
        SparseLinearSystem(const SparseAddress& addr): A_(addr),
                                                 n_(A_.dim()), 
                                                 X_(n_),
                                                 b_(n_),
                                                 r_(n_) {}

        // Returns the dimension of the sparse matrix
        int size() const{
            return n_ ;
        }

        // Non-constant reference to the sparse matrix
        SparseMatrix& A(){ return A_ ;}

        // Constant reference to the sparse matrix
        const SparseMatrix& A() const{ return A_ ;}

        // Non-constant access to the solution vector
        LinearVector& X(){ return X_ ;}

        // Constant access to the solution vector
        const LinearVector& X() const{ return X_ ;}

        // Non-constant access to the source vector
        LinearVector& b(){ return b_ ;}

        // Constant access to the source vector
        const LinearVector& b() const{ return b_ ;}

        // Non-constant access to the residual
        LinearVector& r(){ return r_ ;}

        // Constant access to the residual
        const LinearVector& r() const{ return r_ ;}

        // Operator overload for addition of two sparse linear system
        SparseLinearSystem operator+(const SparseLinearSystem& L) const{

            if (n_ != L.size()) throw std::runtime_error("System size mismatch");

            SparseLinearSystem Lr((*this).A().addressing()) ;
            SparseMatrix Ar = (*this).A() + L.A() ;
            for(int i = 0; i < n_; ++i){
                Lr.A()(i,i) = Ar(i,i) ;
                for(int j = (*this).A().addressing().rowArray()[i]; j < (*this).A().addressing().rowArray()[i+1]; ++j){
                    Lr.A()(i, (*this).A().addressing().columnIndexArray()[j]) = Ar(i, (*this).A().addressing().columnIndexArray()[j]) ;
                }
            }

            Lr.b() = (*this).b() + L.b() ;

            return Lr ;
        }

        // Operator overload for substraction of two sparse linear system
        SparseLinearSystem operator-(const SparseLinearSystem& L) const{

            if (n_ != L.size()) throw std::runtime_error("System size mismatch");

            SparseLinearSystem Lr((*this).A().addressing()) ;
            SparseMatrix Ar = (*this).A() - L.A() ;
            for(int i = 0; i < n_; ++i){
                Lr.A()(i,i) = Ar(i,i) ;
                for(int j = (*this).A().addressing().rowArray()[i]; j < (*this).A().addressing().rowArray()[i+1]; ++j){
                    Lr.A()(i, (*this).A().addressing().columnIndexArray()[j]) = Ar(i, (*this).A().addressing().columnIndexArray()[j]) ;
                }
            }

            Lr.b() = (*this).b() - L.b() ;

            return Lr ;
            
        }

        // Operator overload for multiplication of sparse linear system with scalar
        SparseLinearSystem operator*(const double s) const{

            SparseLinearSystem Lr((*this).A().addressing()) ;
            SparseMatrix Ar = (*this).A()*s ;
            for(int i = 0; i < n_; ++i){
                Lr.A()(i,i) = Ar(i,i) ;
                for(int j = (*this).A().addressing().rowArray()[i]; j < (*this).A().addressing().rowArray()[i+1]; ++j){
                    Lr.A()(i, (*this).A().addressing().columnIndexArray()[j]) = Ar(i, (*this).A().addressing().columnIndexArray()[j]) ;
                }
            }

            Lr.b() = (*this).b()*s ;
            Lr.X() = (*this).X() ;

            return Lr ;
        }

        // Operator overload for multiplication of sparse matrix with linear vector in the linear system
        LinearVector operator*(const LinearVector& x) const{
            if(x.size() != n_) throw std::runtime_error("System size mismatch");
            return (*this).A()*x ;
        }

        // Function for multiplication using nominal algorithm
        LinearVector multiplyNominal(const LinearVector& x) const{
            if(x.size() != n_) throw std::runtime_error("System size mismatch");
            return (*this).A().multiplyNominal(x);
        }

        // Function for multiplication using optimised algorithm
        LinearVector multiplyOptimised(const LinearVector& x) const{
            if(x.size() != n_) throw std::runtime_error("System size mismatch");
            return (*this).A().multiplyOptimised(x);
        }

        // Function to calculate residual 
        void calculateResidual(){
            int sizeOfArray = A_.dim() ;      
            LinearVector b_star_ = (*this).A()*((*this).X_) ;
            for(int i=0; i < sizeOfArray; ++i){
                r_(i) = b_(i) - b_star_(i) ;
            }
        }

        // Function to calculate residual using optimised algorithm
        void calculateResidualOptimised(){
            int sizeOfArray = A_.dim() ;      
            LinearVector b_star_ = (*this).A().multiplyOptimised((*this).X_) ;
            for(int i=0; i < sizeOfArray; ++i){
                r_(i) = b_(i) - b_star_(i) ;
            }
        }

        // Function to calculate the average of all the elements in x
        double average(const LinearVector x) const {
            int size = x.size() ;
            double sum = 0 ;
            for(int i=0; i<size; ++i){
                sum += x(i) ;
            }
            return sum/size ;
        }

        // Function to calculate normalisation factor
        int calculateNormalizationFactor(){
            int nf_ = 0 ;
            int sizeOfArray = X_.size() ;
            double avg = average(X_) ;

            LinearVector Xref_(sizeOfArray);
            for(int i=0; i<sizeOfArray; ++i){
                Xref_(i) = avg ;
            }

            LinearVector Ax_ = A_*X_ ;
            LinearVector Axref_ = A_*Xref_ ;

            for(int i = 0; i < sizeOfArray; ++i){
                nf_ += (abs(Ax_(i)-Axref_(i))+abs(b_(i)-Axref_(i))) ;
            }

            return nf_ ;
        }

        // Function to calculate normalised one norm
        double normalisedOneNorm(){
            return (r_.oneNorm()/calculateNormalizationFactor()) ;
        }

        // Function to calculate normalised two norm
        double normalisedTwoNorm(){
            return (r_.twoNorm()) ;
        }

        // Function to calculate normalised infinity norm
        double normalisedInfinityNorm(){
            return (r_.infinityNorm()/calculateNormalizationFactor()) ;
        }

        // Function to write vector values as csv file for visualisation
        void writeVectorCSV(const LinearVector& X, fs::path filename){
            ofstream out(filename) ;
            out << "Index,Value,Z" << endl ;
            for(int i=0; i<X.size(); ++i){
                out << i << "," << X(i) << "," << 0 << endl ;
            }
        }

        // Function to write matrix values as csv file for visualisation
        void writeMatrixCSV(const SparseMatrix& D, fs::path filename){
            ofstream out(filename) ;
            out << "Row,Column,Value,Z" << endl ;
            for(int i=0; i<D.dim(); ++i){
                for(int j=0; j<D.dim(); ++j){
                    out << D.dim() - 1 - i << "," << j << "," << D(i,j) << "," << 0 << endl ;
                }
            }
        }

} ;

// Operator to read matrix from mtx file
inline std::istream& operator>>(std::istream& in, SparseLinearSystem& sys) {
    if (!in) {
        throw std::runtime_error("Failed to open input stream");
    }
    sys.A().readMatrix(in);
    return in;
}

// Operator to write matrix to mtx file
inline std::ostream& operator<<(std::ostream& out, const SparseLinearSystem& sys) {
    if (!out) {
        throw std::runtime_error("Failed to open output stream");
    }
    sys.A().writeMatrix(out);
    return out;
}