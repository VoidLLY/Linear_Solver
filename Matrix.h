#pragma once
#include<vector>

template <class T>
class Matrix
{
public:

    // constructor where we want to preallocate ourselves
    Matrix(int rows, int cols, bool preallocate);
    // constructor where we already have allocated memory outside
    Matrix(int rows, int cols, T* values_ptr);
    //copy constructor
    Matrix(const Matrix& source);

    // destructor
    virtual ~Matrix();

    // Print out the values in our matrix
    void printValues();
    virtual void printMatrix();

    // Perform some operations with our matrix
    void deepcopy(const Matrix<T>& source);
    void matMatMult(Matrix<T>& mat_right, Matrix<T>& output);
    //void MtoArr(vector<vector<T>> &submat);
    bool judge_symmetry();
    T MAX(Matrix<T>& a);
    T NORM(Matrix<T>& a, Matrix<T>& b);



    void swapRows(int j, int k);


    //void swapRowsPartialPivoting(int j, int k, T** P);
    int find_max_index(int k);
    void back_subtitution(Matrix<T>& Y, Matrix<T>& output);
    void forward_subtitution(Matrix<T>& L, Matrix<T>& Y, Matrix<T>& mat_b);
    //T** find_transpose(T** P);
    void get_transpose();

    //Perform linear solver
    void LU(Matrix<T>& mat_b, Matrix<T>& output);
    void Chol(Matrix<T>& mat_b, Matrix<T>& output);
    void Jacobi(Matrix<T>& mat_b, Matrix<T>& output, Matrix<T>& initial);
    void Gauss_Seidel(Matrix<T>& mat_b, Matrix<T>& output, Matrix<T>& initial);
    void Gauss_Seidel_SOR(Matrix<T>& mat_b, Matrix<T>& output, Matrix<T>& initial, double w);

    // Explicitly using the C++11 nullptr here
    T* values = nullptr;
    int rows = -1;
    int cols = -1;

    // We want our subclass to know about this
protected:
    bool preallocated = false;

    // Private variables - there is no need for other classes 
    // to know about these variables
private:

    int size_of_values = -1;
};