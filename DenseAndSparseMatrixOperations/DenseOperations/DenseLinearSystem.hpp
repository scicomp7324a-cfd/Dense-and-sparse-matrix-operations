# pragma once

# include <iostream>
# include <vector>
# include <filesystem>
# include <stdexcept>

# include "DenseMatrix.hpp"
# include "LinearVector.hpp"

using namespace std ;
namespace fs = filesystem ;


/* *******************************************************************
    DenseLinearSystem.hpp - This class contains all the algebraic 
    operations for dense linear system and residual calculation  
    functions.
******************************************************************* */

class DenseLinearSystem{

    private:

        int n_ ;
        DenseMatrix A_ ;
        LinearVector X_ ;
        LinearVector b_ ;
        LinearVector r_ ;


    public:

        DenseLinearSystem(int n): n_(n),
                                A_(n),
                                X_(n),
                                b_(n),
                                r_(n){}

        int size() const{
            return n_ ;
        }

        // Non-constant reference to the dense matrix
        DenseMatrix& A(){ return A_ ; }

        // Constant reference to the dense matrix
        const DenseMatrix& A() const{ return A_ ; }

        // Non-constant access to the solution vector
        LinearVector& X(){ return X_ ;}

        // Constant access to the solution vector
        const LinearVector& X() const { return X_ ;}
        
        // Non-constant access to the source vector
        LinearVector& b(){ return b_; }

        // Constant access to the source vector
        const LinearVector& b() const { return b_ ;}

        // Non-constant access to the residual
        LinearVector& r(){ return r_; }

        // Constant access to the residual
        const LinearVector& r() const { return r_ ;}

        // Operator overload for addition of two dense linear system
        DenseLinearSystem operator+(const DenseLinearSystem& L) const {

            if (n_ != L.size()) throw std::runtime_error("System size mismatch");

            DenseLinearSystem Lr(n_) ;

            Lr.A() = (*this).A() + L.A() ;
            Lr.b() = (*this).b() + L.b() ;

            return Lr ;
        }

        // Operator overload for substraction of two dense linear system
        DenseLinearSystem operator-(const DenseLinearSystem& L) const{

            if (n_ != L.size()) throw std::runtime_error("System size mismatch");

            DenseLinearSystem Lr(n_) ;

            Lr.A() = (*this).A() - L.A() ;
            Lr.b() = (*this).b() - L.b() ;

            return Lr ;
            
        }

        // Operator overload for multiplication of dense linear system with scalar
        DenseLinearSystem operator*(const double s) const{

            DenseLinearSystem Lr(n_) ;

            Lr.A() = (*this).A()*s ;
            Lr.b() = (*this).b()*s ;
            Lr.X() = (*this).X() ;

            return Lr ;
        }

        // Operator overload for multiplication of dense matrix with linear vector in the linear system
        LinearVector operator*(const LinearVector& x) const{
            if(x.size() != n_) throw std::runtime_error("System size mismatch");
            return (*this).A()*x ;
        }

        // Function for multiplication using nominal algorithm
        LinearVector multiplyNominal(const LinearVector& x) const{
            if(x.size() != n_) throw std::runtime_error("System size mismatch");
            return (*this).A().multiplyNominal(x) ;
        }

        // Function for multiplication using optimised algorithm
        LinearVector multiplyOptimised(const LinearVector& x) const{
            if(x.size() != n_) throw std::runtime_error("System size mismatch");
            return (*this).A().multiplyBLAS(x) ;
        }

        // Function to calculate residual 
        void calculateResidual(){
            LinearVector Ax = (*this)*(*this).X() ;
            for(int i=0; i < n_; ++i){
                (*this).r()(i) = (*this).b()(i) - Ax(i) ;
            }
        }

        // Function to calculate residual using optimised algorithm
        void calculateResidualOptimised(){
            LinearVector Ax = (*this).A().multiplyBLAS((*this).X()) ;
            for(int i=0; i < n_; ++i){
                (*this).r()(i) = (*this).b()(i) - Ax(i) ;
            }
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
        void writeMatrixCSV(const DenseMatrix& D, fs::path filename){
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
inline std::istream& operator>>(std::istream& in, DenseLinearSystem& sys) {
    if (!in) {
        throw std::runtime_error("Failed to open input stream");
    }
    sys.A().readMatrix(in);
    return in;
}

// Operator to write matrix to mtx files
inline std::ostream& operator<<(std::ostream& out, const DenseLinearSystem& sys) {
    if (!out) {
        throw std::runtime_error("Failed to open output stream");
    }
    sys.A().writeMatrix(out);
    return out;
}