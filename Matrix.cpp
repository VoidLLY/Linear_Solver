#include <iostream>
#include "Matrix.h"
#include <cstdlib>
#include <vector>

using namespace std;

// Constructor - using an initialisation list here
template <class T>
Matrix<T>::Matrix(int rows, int cols, bool preallocate) : rows(rows), cols(cols), size_of_values(rows* cols), preallocated(preallocate)
{
    // If we want to handle memory ourselves
    if (this->preallocated)
    {
        // Must remember to delete this in the destructor
        this->values = new T[size_of_values];
    }
}

// Constructor - now just setting the value of our T pointer
template <class T>
Matrix<T>::Matrix(int rows, int cols, T* values_ptr) : rows(rows), cols(cols), size_of_values(rows* cols), values(values_ptr)
{}


//copy constructor
template <class T>
Matrix<T>::Matrix(const Matrix<T>& source)
{
    deepcopy(source);
}


// destructor
template <class T>
Matrix<T>::~Matrix()
{
    // Delete the values array
    if (this->preallocated) {
        delete[] this->values;
    }
}

// Just print out the values in our values array
template <class T>
void Matrix<T>::printValues()
{
    std::cout << "Printing values" << std::endl;
    for (int i = 0; i < this->size_of_values; i++)
    {
        std::cout << this->values[i] << " ";
    }
    std::cout << std::endl;
}

// Explicitly print out the values in values array as if they are a matrix
template <class T>
void Matrix<T>::printMatrix()
{
    std::cout << "Printing matrix" << std::endl;
    for (int j = 0; j < this->rows; j++)
    {
        std::cout << std::endl;
        for (int i = 0; i < this->cols; i++)
        {
            // We have explicitly used a row-major ordering here
            std::cout << this->values[i + j * this->cols] << " ";
        }
    }
    std::cout << std::endl;
}

template <class T>
void Matrix<T>::deepcopy(const Matrix<T>& source)
{
    // deallocate the value the pointer points to
    if (this->preallocated) {
        delete[] this->values;
    }
    //shallow copy ordinary paramaters
    rows = source.rows;
    cols = source.cols;
    size_of_values = source.size_of_values;
    if (source.values)
    {
        //allocate memory for the copy
        values = new T[size_of_values];
        // do the copy
        for (int i = 0; i < size_of_values; i++)
            values[i] = source.values[i];
    }
    else
        values = nullptr;
}

// Do matrix matrix multiplication
// output = this * mat_right
template <class T>
void Matrix<T>::matMatMult(Matrix<T>& mat_right, Matrix<T>& output)
{

    // Check our dimensions match
    if (this->cols != mat_right.rows)
    {
        std::cerr << "Input dimensions for matrices don't match" << std::endl;
        return;
    }

    // Check if our output matrix has had space allocated to it
    if (output.values != nullptr)
    {
        // Check our dimensions match
        if (this->rows != output.rows || mat_right.cols != output.cols)
        {
            std::cerr << "Input dimensions for matrices don't match" << std::endl;
            return;
        }
    }
    // The output hasn't been preallocated, so we are going to do that
    else
    {
        output.values = new T[this->rows * mat_right.cols];
    }

    // Set values to zero before hand
    for (int i = 0; i < output.size_of_values; i++)
    {
        output.values[i] = 0;
    }

    // Now we can do our matrix-matrix multiplication
    // CHANGE THIS FOR LOOP ORDERING AROUND
    // AND CHECK THE TIME SPENT
    // Does the ordering matter for performance. Why??
    for (int i = 0; i < this->rows; i++)
    {
        for (int k = 0; k < this->cols; k++)
        {
            for (int j = 0; j < mat_right.cols; j++)
            {
                output.values[i * output.cols + j] += this->values[i * this->cols + k] * mat_right.values[k * mat_right.cols + j];
            }
        }
    }
}

//template <class T>
//void Matrix<T>::MtoArr(vector<vector<T>> &submat)
//{
//    for (int i = 0;i < this->rows;i++)
//    {
//        for (int j = 0;j < this->rows;j++)
//        {
//            submat[i][j] = this->values[i * rows + j];
//        }
//    }
//}

//template <class T>
//void Matrix<T>::computedet()
//{
//    vector<vector<double>> swp(this->rows);
//    for (int i = 0;i < swp.size();i++)
//    {
//        swp[i].resize(this->cols);
//    }
//    for (int i = 0;i < swp.size();i++)
//    {
//        for (int j = 0;j < swp[0].size();j++)
//        {
//            swp[i][j] = this->values[i * this->cols + j];
//        }
//    }
//    double d = 0;
//    double det(int n,vector<vector<double>>mat)
//    {
//        int c, subi, i, j, subj;
//        vector<vector<double>> submat(n,vector<double>(n));
//        if (n == 1)
//        {
//            return mat[0][0];
//        }
//        if (n == 2)
//        {
//            return((mat[0][0] * mat[1][1]) - (mat[1][0] * mat[0][1]));
//        }
//        else
//        {
//            for (c = 0; c < n; c++)
//            {
//                subi = 0;
//                for (i = 1; i < n; i++)
//                {
//                    subj = 0;
//                    for (j = 0; j < n; j++)
//                    {
//                        if (j == c)
//                        {
//                            continue;
//                        }
//                        submat[subi][subj] = mat[i][j];
//                        subj++;
//                    }
//                    subi++;
//                }
//                d = d + (pow(-1, c) * mat[0][c] * det(n - 1, submat));
//            }return d;
//        }
//        
//    }
//    cout << "The determinant of the matirx is:"<<det(this->rows,use_mat)<<endl;
//
//  
//}
template <class T>
void Matrix<T>::LU(Matrix<T>& mat_b, Matrix<T>& output)
{
    // Set values to zero before hand
    for (int i = 0; i < output.size_of_values; i++)
    {
        output.values[i] = 0;
    }

    // initiate partial pivoting matrix
    auto* P = new Matrix<T>(this->rows, this->cols, true);

    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < this->cols; j++)
        {
            P->values[i * this->cols + j] = 0;
            if (i == j)
            {
                P->values[i * this->cols + j] = 1;
            }
        }
    }

    // initialize L matrix to zero
    auto* L = new Matrix<T>(this->rows, this->cols, true);

    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < this->cols; j++)
        {
            L->values[i * this->cols + j] = 0;
        }
    }

    for (int k = 0; k < this->rows - 1; k++)
    {
        // find the index of maximum magnitude
        int j = find_max_index(k);

        swapRows(j, k); // do the swap of A's rows

        P->swapRows(j, k);  // record the swap in the matrix P
        L->swapRows(j, k);
        for (int i = k + 1; i < this->rows; i++)
        {
            T s = this->values[i * this->cols + k] / this->values[k * this->cols + k];
            for (int m = k; m < this->cols; m++)
            {
                this->values[i * this->cols + m] -= this->values[k * this->cols + m] * s;
            }
            L->values[i * this->cols + k] = s;
        }
    }

    for (int col_n = 0; col_n < this->rows; col_n++)
    {
        L->values[col_n * this->cols + col_n] += 1;
    }


    // update mat_b with P*mat_b

    auto* matpb = new Matrix<T>(mat_b.rows, mat_b.cols, true);
    P->matMatMult(mat_b, *matpb);

    auto* Y = new Matrix<T>(mat_b.rows, mat_b.cols, true);
    forward_subtitution(*L, *Y, *matpb);

    back_subtitution(*Y, output);

    // print output x
    cout << "This is LU decomposition method:" << endl;
    cout << " X = :" << endl;
    output.printMatrix();
    delete P;
    P = nullptr;
    delete L;
    L = nullptr;

    delete matpb;
    matpb = nullptr;
}

// method to swap rows of A
template <class T>
void Matrix<T>::swapRows(int j, int k)
{
    for (int col = 0; col < this->cols; col++)
    {
        auto temp = this->values[j * this->cols + col];
        this->values[j * this->cols + col] = this->values[k * this->cols + col];
        this->values[k * this->cols + col] = temp;
        //swap(this->values[j * this->cols + col], this->values[k * this->cols + col]);
    }
}

template <class T>

int Matrix<T>::find_max_index(int k)
{
    int start = k;
    int max_index = start;
    for (int r = start; r < this->rows; r++)
    {
        if (abs(this->values[r * this->cols + k]) > abs(this->values[max_index * this->cols + k]))
        {
            max_index = r;
        }
    }

    return max_index;
}

template <class T>
// matrix A is on upper triangular form
void Matrix<T>::back_subtitution(Matrix<T>& Y, Matrix<T>& output)
{
    output.values[this->rows - 1] = Y.values[this->rows - 1] / this->values[this->size_of_values - 1];
    for (int k = this->rows - 2; k > -1; k--)
    {
        T s = 0;
        for (int j = k + 1; j < this->rows; j++)
        {
            s = s + this->values[k * this->cols + j] * output.values[j];
        }
        output.values[k] = (Y.values[k] - s) / this->values[k * this->cols + k];
    }

}


// for lower triangular form matrix
template <class T>
void Matrix<T>::forward_subtitution(Matrix<T>& L, Matrix<T>& Y, Matrix<T>& mat_b)
{
    Y.values[0] = mat_b.values[0] / L.values[0];
    for (int i = 1; i < L.rows; i++)
    {
        T s = 0;
        for (int j = 0; j < i; j++)
        {
            s = s + L.values[i * L.cols + j] * Y.values[j];
        }
        Y.values[i] = (mat_b.values[i] - s) / L.values[i * L.cols + i];
    }
}

template <class T>
void Matrix<T>::get_transpose()
{
    T temp;
    for (int i = 0; i < this->rows; i++)
    {
        for (int j = i; j < this->cols; j++)
        {
            temp = this->values[i * this->cols + j];
            this->values[i * this->cols + j] = this->values[j * this->cols + i];
            this->values[j * this->cols + i] = temp;
        }
    }
}

//function to implement Cholesky Decompostion method
template <class T>
void Matrix<T>::Chol(Matrix<T>& mat_b, Matrix<T>& output)
{
    cout << "The is Cholesky Decomposition method:" << endl;
    // Set values to zero before hand
    for (int i = 0; i < output.size_of_values; i++)
    {
        output.values[i] = 0;
    }

    auto* U = new Matrix<T>(this->rows, this->cols, true);
    //initiate upper triangle matrix U and ouput matrix
    //for transpose(U)*U = A
    for (int i = 0; i < U->size_of_values; i++)
    {
        U->values[i] = 0;
    }

    if (this->judge_symmetry() == 0)
    {
        cout << "We can't use Cholesky Decomposition method to solve this problem!" << endl;
        return;
    }


    if (this->values[0] <= 0)
    {
        cerr << "This matrix A isn't positive definete" << endl;
        return;
    }
    U->values[0] = sqrt(this->values[0]);//U[0,0]=sqrtA[0,0]
    for (int i = 0; i < this->cols; i++)
    {
        U->values[i] = this->values[i] / U->values[0];
        //U[0,i]=A[0,i]/U[0,0]
        //the first row of U
    }

    for (int i = 1; i < this->rows; i++)
    {
        // starts from the second row

        //compute U[i,i]
        T subsum = U->values[0] - U->values[0];
        int k = 0;
        do {
            subsum += pow(U->values[k * cols + i], 2);
            k++;
        } while (k == i - 1);


        T sum = this->values[i * this->cols + i] - subsum;//A[i,i]-sum(U[k,i]^2),(k=1 to i-1)
        if (sum <= 0)
        {

            cerr << "This matrix A isn't positive definete" << endl;
            return;
        }
        else
        {
            U->values[i * this->cols + i] = sqrt(sum);
            for (int j = i + 1; j < this->cols; j++)
            {
                T ssubsum = U->values[0] - U->values[0];
                int kk = 0;
                do {
                    ssubsum += U->values[kk * cols + i] * U->values[kk * cols + j];
                    kk++;
                } while (kk == i - 1);

                T ssum = this->values[i * this->cols + j] - ssubsum;
                U->values[i * this->cols + j] = ssum / U->values[i * this->cols + i];
            }
        }
    } 
    auto* L = new Matrix<T>(this->rows, this->cols, true);
    *L = Matrix<T>(*U);//deep copy 
    L->get_transpose();//transpose U to get L

    auto* Y = new Matrix<T>(mat_b.rows, mat_b.cols, true);
    forward_subtitution(*L, *Y, mat_b);

    U->back_subtitution(*Y, output);

    // print output x

    cout << " X = :" << endl;
    output.printMatrix();
    delete L;
    L = nullptr;
    delete U;
    U = nullptr;

}

template<class T>
bool Matrix<T>::judge_symmetry()
{
    bool symmetry = true;
    for (int i = 0; i < this->rows; i++)
    {
        for (int j = 0; j < this->cols; j++)
        {

            if (this->values[i * this->cols + j] != this->values[j * this->cols + i])
            {
                cerr << "!Warning! The matrix A isn't symmetrical!" << endl;
                symmetry = false;
            }
        }
    }
    return symmetry;
}




template<class T>
T Matrix<T>::MAX(Matrix<T>& a)
{
    T max = 0;
    for (int i = 0; i < a.rows * a.cols; i++)
    {
        if (fabs(a.values[i]) > max)
            max = fabs(a.values[i]);
    }
    return max;
}

template<class T>
T Matrix<T>::NORM(Matrix<T>& a, Matrix<T>& b)
{
    T norm = 0;
    T sum = 0;
    for (int i = 0; i < a.rows * a.cols; i++)
    {
        sum += pow((a.values[i] - b.values[i]), 2);
    }
    norm = sqrt(sum);
    return norm;
}

template<class T>
void Matrix<T>::Jacobi(Matrix<T>& mat_b, Matrix<T>& output, Matrix<T>& initial)
{
    int count = 0;

    T deviation;
    bool flag;

    // Set values to zero before hand
    for (int i = 0; i < output.size_of_values; i++)
    {
        output.values[i] = 0;
    }

    do {
        // to store resursive times
        count++;
        for (int j = 0; j < this->cols; j++)
            initial.values[j] = output.values[j];

        for (int i = 0; i < this->cols; i++)
        {
            T sum = 0;
            for (int j = 0; j < this->cols; j++)
            {
                if (j == i)
                    continue;
                else
                    sum += this->values[i * this->cols + j] * initial.values[j];
            }
            output.values[i] = (mat_b.values[i] - sum) / this->values[i * this->cols + i];
        }

        deviation = fabs(MAX(initial) - MAX(output));

        //set accurary = 0.000001, where the iteration should stop
        if (deviation < 0.000001)
        {
            flag = 0;
        }
        else
        {
            flag = 1;
        }

        if (count > 5000)
        {
            flag = 0;
            cout << "The solution doesn't converge. " << endl;
        }
    } while (flag);

    cout << "this is Jacobi method  " << endl;
    for (int i = 0; i < this->cols; i++)
        cout << "x[" << i + 1 << "]= " << output.values[i] << endl;
    cout << "recursive times:  " << count << endl;
}

template<class T>
void Matrix<T>::Gauss_Seidel(Matrix<T>& mat_b, Matrix<T>& output, Matrix<T>& initial)
{
    int count = 0;

    T deviation;
    bool flag;

    // Set values to zero before hand
    for (int i = 0; i < output.size_of_values; i++)
    {
        output.values[i] = 0;
    }

    do {
        // to store resursive times
        count++;

        for (int j = 0; j < this->cols; j++)
            initial.values[j] = output.values[j];

        for (int i = 0; i < this->cols; i++)
        {
            T sum1 = 0;
            T sum2 = 0;
            for (int j = 0; j < i; j++)
            {
                sum1 += this->values[i * this->cols + j] * output.values[j];
            }
            for (int j = i + 1; j < this->cols; j++)
            {
                sum2 += this->values[i * this->cols + j] * initial.values[j];
            }
            output.values[i] = (mat_b.values[i] - sum1 - sum2) / this->values[i * this->cols + i];
        }

        deviation = fabs(MAX(initial) - MAX(output));

        //set accurary = 0.000001, where the iteration should stop
        if (deviation < 0.000001)
        {
            flag = 0;
        }
        else
        {
            flag = 1;
        }

        if (count > 100)
        {
            flag = 0;
            cout << "The solution doesn't converge. " << endl;
        }
    } while (flag);

    cout << "this is Gauss-Seidel method  " << endl;
    for (int i = 0; i < this->cols; i++)
        cout << "x[" << i + 1 << "]= " << output.values[i] << endl;
    cout << "recursive times:  " << count << endl;
}

template<class T>
void Matrix<T>::Gauss_Seidel_SOR(Matrix<T>& mat_b, Matrix<T>& output, Matrix<T>& initial, double w)
{
    int count = 0;

    T deviation;
    T dnorm;
    T dnorm_pre;
    bool flag;

    // Set values to zero before hand
    for (int i = 0; i < output.size_of_values; i++)
    {
        output.values[i] = 0;
        initial.values[i] = 0;
    }

    do {
        // to store resursive times
        count++;

        dnorm_pre = fabs(NORM(initial, output));

        for (int j = 0; j < this->cols; j++)
            initial.values[j] = output.values[j];

        for (int i = 0; i < this->cols; i++)
        {
            T sum1 = 0;
            T sum2 = 0;
            for (int j = 0; j < i; j++)
            {
                sum1 += this->values[i * this->cols + j] * output.values[j];
            }
            for (int j = i + 1; j < this->cols; j++)
            {
                sum2 += this->values[i * this->cols + j] * initial.values[j];
            }
            output.values[i] = w * ((mat_b.values[i] - sum1 - sum2) / this->values[i * this->cols + i]) + (1 - w) * initial.values[i];
        }

        deviation = fabs(MAX(initial) - MAX(output));
        dnorm = fabs(NORM(initial, output));

        //set accurary = 0.000001, where the iteration should stop
        if (deviation > 0.000001)
        {
            flag = 1;
            if ((count > 5) && (count < 100))
            {
                w = 2 / (1 + sqrt(1 - dnorm / dnorm_pre));
            }
            else if (count >= 100)
            {
                flag = 0;
                cout << "The solution doesn't converge. " << endl;
            }
        }
        else
        {
            flag = 0;
        }
    } while (flag);

    cout << "this is Gauss-Seidel-SOR method  " << endl;
    for (int i = 0; i < this->cols; i++)
        cout << "x[" << i + 1 << "]= " << output.values[i] << endl;
    cout << "recursive times:  " << count << endl;
}