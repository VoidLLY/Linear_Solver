#pragma once
#include "Matrix.h"

template <class T>
class CSRMatrix : public Matrix<T>
{
public:

    // constructor where we want to preallocate ourselves
    CSRMatrix(int rows, int cols, int nnzs, bool preallocate);
    // constructor where we already have allocated memory outside
    CSRMatrix(int rows, int cols, int nnzs, T* values_ptr, int* row_position, int* col_index);

    //copy constructor 
    CSRMatrix(const CSRMatrix<T>& mat);
    // destructor
    ~CSRMatrix();

    // Print out the values in our matrix
    virtual void printMatrix();

    // Perform some operations with our matrix
    void deepcopy(const CSRMatrix<T>& mat);
    
    // Perform some operations with our matrix
    void matVecMult(double* input, double* output);

    //Judge if the sparse matrix is SPD

    //bool judge_SPD(CSRMatrix Umat);
    bool judge_SPD(CSRMatrix<T>* U);
    //Transpose the sparse matrix
    void CSRtranspose(CSRMatrix<T>& transposed);

    //Cholesky Decompostion
    void Chol_sparse(double* b, double* output);

    //LU factorisation related function 
    virtual void LU_sparse(double* b, double* output);
    virtual void swapRows(int j, int k);
    virtual int find_max_index(int k);
    virtual void back_subtitution(T* Y, T* output);
    virtual void forward_subtitution(CSRMatrix<T>& L, T* Y, T* b);
    
    //Several iterative solvers related function
    virtual T MAX(T* a);
    virtual void Jacobi(T* mat_b, T* output, T* initial);
    virtual void Gauss_Seidel(T* mat_b, T* output, T* initial);

    // Explicitly using the C++11 nullptr here
    int* row_position = nullptr;
    int* col_index = nullptr;

    // How many non-zero entries we have in the matrix
    int nnzs = -1;

    // Private variables - there is no need for other classes 
    // to know about these variables
private:

};