# pragma once

# include<iostream>
# include<vector>
# include<istream>
# include<ostream>
# include<string>
# include<sstream>
# include "LinearVector.hpp"
# include <cblas.h>

using namespace std ;

/* *******************************************************************
    DenseMatrix.hpp - This class contains all the algebraic 
    operations for dense matrix and calculation of other matrix 
    related metrics functions.
******************************************************************* */

class DenseMatrix{

    private:

        int n_ ;
        vector<double> A_ ;       

    public: 

        // Constructs the dense matrix using the matrix dimension
        DenseMatrix(int n):n_(n){
            A_.resize(n_*n_) ;
        }

        // Constant access to matrix data
        double* data() {
            return A_.data();
        }

        // Non-constant access to matrix data
        const double* data() const {
            return A_.data();
        }

        // Non-constant operator to access elements of dense matrix
        double& operator()(int i, int j) {
            return A_[i*n_ + j] ;
        }

        // Constant operator to access elements of dense matrix
        const double& operator()(int i, int j) const{
            return A_[i*n_ + j] ;
        }

        // Function to return dimension of the dense matrix
        const int& dim() const{
            return n_ ;
        }

        // Operator overload for adding dense matrices
        DenseMatrix operator+(const DenseMatrix& B) const{

            DenseMatrix R(n_) ;
            for(int i=0 ; i < n_ ; ++i){
                for(int j=0 ; j < n_; ++j){
                    R(i,j) = (*this)(i,j) + B(i,j) ;
                }
            }

            return R ; 
        }

        // Operator overload for substracting dense matrices
        DenseMatrix operator-(const DenseMatrix& B) const{

            DenseMatrix R(n_) ;
            for(int i=0 ; i < n_ ; ++i){
                for(int j=0 ; j < n_; ++j){
                    R(i,j) = (*this)(i,j) - B(i,j) ;
                }
            }

            return R ; 
        }

        // Operator overload for multiplying dense matrix with scalar
        DenseMatrix operator*(const double s) const{
            
            DenseMatrix R(n_) ;

            for(int i=0 ; i < n_ ; ++i){
                for(int j=0 ; j < n_; ++j){
                    R(i,j) = (*this)(i,j)*s;
                }
            }

            return R ; 
        }

        // Function to check if the dense matrix is symmetric
        bool symmetric() const{
            for(int i = 0; i < n_; ++i){
                for(int j = i+1 ; j < n_; ++j){
                    if((*this)(i,j) != (*this)(j,i)){
                        return false ;
                    }
                }
            }
            return true ;
        }

        // Function to check diagonal dominance of a dense matrix
        bool diagonalDominance() const{

            double sum = 0;
            for(int i=0; i < n_; ++i){
                for(int j= 0; j < n_; ++j){
                    sum += abs((*this)(i,j)) ;
                }
                sum -= abs((*this)(i,i)) ;
                if(sum > abs((*this)(i,i))){
                    return false ;
                }
                sum = 0 ;
            }
            return true ;
        }

        // Function to check strict diagonal dominance of a dense matrix
        bool strictDiagonalDominance() const{

            double sum = 0;
            for(int i=0; i < n_; ++i){
                for(int j= 0; j < n_; ++j){
                    sum += abs((*this)(i,j)) ;
                }
                sum -= abs((*this)(i,i)) ;
                if(sum >= abs((*this)(i,i))){
                    return false ;
                }
                sum = 0 ;
            }
            return true ;
        }

        // Function to check if all the diagonal elements of the matrix are positive
        bool allPositiveDiagonal() const{
            for(int i=0; i<n_; ++i){
                if((*this)(i,i) <= 0){
                    return false ;
                }
            }
            return true ;
        }

        // Operator overload to multiply dense matrix with a linear vector
        LinearVector operator*(const LinearVector& X) const{

            LinearVector res(n_) ;

            for(int i = 0; i < n_; ++i){
                for(int j = 0; j < n_; ++j){
                    res(i) += (*this)(i,j)*X(j) ;
                }
            }

            return res ;
        }

        // Function for vector-matrix multiplication using a nominal algorithm
        LinearVector multiplyNominal(const LinearVector& X) const{

            LinearVector res(n_) ;

            for(int i = 0; i < n_; ++i){
                for(int j = 0; j < n_; ++j){
                    res(i) += (*this)(i,j)*X(j) ;
                }
            }

            return res ;
        }

        // Function for vector-matrix multiplication using a optimised algorithm
        LinearVector multiplyBLAS(const LinearVector& X) const{

            LinearVector res(n_) ;

            cblas_dgemv(
                CblasRowMajor,   // your matrix uses A_[i*n_ + j]
                CblasNoTrans,
                n_,             // rows
                n_,             // cols
                1.0,            // alpha
                A_.data(),      // matrix data
                n_,             // leading dimension for row-major
                X.data(),       // input vector
                1,              // stride of X
                0.0,            // beta
                res.data(),     // output vector
                1               // stride of res
            );

            return res;
        }

        // Utility function to split a string at white space
        vector<string> splitWhiteSpace(string s) {
            istringstream iss(s);
            vector<string> out;
            for (string tok; iss >> tok; ) out.push_back(tok);
            return out;
        }

        // Utility function to strip white saces at the stat and end of the string
        string stripWhiteSpace(string s){
            auto first = s.find_first_not_of(" \t\r");
            if (first == string::npos) return ""; 
            auto last  = s.find_last_not_of(" \t\r");
            return s.substr(first, last - first + 1);
        }

        // Function to read matrix coefficients from mtx file
        void readMatrix(istream& in){

            string line ;
            vector<string> strVec ;

            while(getline(in,line)){
                
                if(line.empty() || line[0] == '%'){
                    continue ;
                }
                else{
                    line = stripWhiteSpace(line) ;
                    strVec = splitWhiteSpace(line) ;

                    if(strVec[0]!=(strVec[1]) || stoi(strVec[0])!=n_) throw std::runtime_error("inconsistent matrix dimensions") ;
                    
                    int nonZeros  = stoi(strVec[2]) ;
                    for(int i=0; i<nonZeros; ++i){

                        getline(in,line) ;
                        line = stripWhiteSpace(line) ;
                        strVec = splitWhiteSpace(line) ;

                        (*this)(stoi(strVec[0])-1, stoi(strVec[1])-1) = stod(strVec[2]) ;
                    }
                }
            }
        }

        // Function to write matrix coefficients to mtx file
        void writeMatrix(ostream& out) const{

            int nonZeros = findNumNonZeros() ;
            out << "%%MatrixMarket matrix coordinate real general" << endl ;
            out << "%Comments" << endl ;
            out << "%There are " << nonZeros << " non zero entries in the following matrix" << endl ;
            out << n_ << " " << n_ << " " << nonZeros << endl ;

            for(int i = 0; i < n_; ++i){
                for(int j = 0; j < n_; ++j){
                    if((*this)(i,j) != 0){
                        out << i+1 << " " << j+1 << " " << (*this)(i,j) << endl ;
                    }
                }
            }

            out << endl ;
        }

        // Function to find the number of non-zeros in the matrix. 
        int findNumNonZeros() const{
            int count = 0 ;
            for(int i=0; i<n_; ++i){
                for(int j=0; j<n_; ++j){
                    if((*this)(i,j) != 0){
                        ++count ;
                    }
                }
            }
            return count ;
        }

} ;