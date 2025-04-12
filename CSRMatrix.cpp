#include <iostream>
#include "CSRMatrix.h"
#include <cstdlib>

using namespace std;
// Constructor - using an initialisation list here
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, bool preallocate) : Matrix<T>(rows, cols, false), nnzs(nnzs)
{
    // If we don't pass false in the initialisation list base constructor, it would allocate values to be of size
    // rows * cols in our base matrix class
    // So then we need to set it to the real value we had passed in
    this->preallocated = preallocate;

    // If we want to handle memory ourselves
    if (this->preallocated)
    {
        // Must remember to delete this in the destructor
        this->values = new T[this->nnzs];
        this->row_position = new int[this->rows + 1];
        this->col_index = new int[this->nnzs];
    }
}

// Constructor - now just setting the value of our T pointer
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, T* values_ptr, int* row_position, int* col_index) : Matrix<T>(rows, cols, values_ptr), nnzs(nnzs), row_position(row_position), col_index(col_index)
{}

////copy constructor
//template <class T>
//CSRMatrix<T>::CSRMatrix(const CSRMatrix& mat)
//{
//    this->rows = mat.rows;
//    this->cols = mat.cols;
//    this->nnzs = mat.nnzs;
//    this->preallocated = mat.preallocated;
//    this->values = new T[this->nnzs];
//    this->row_position = new int[this->rows + 1];
//    this->col_index = new int[this->nnzs];
//
//    for (int i = 0; i < nnzs; i++)
//    {
//        this->values[i] = mat->values[i];
//        this->col_index[i] = mat->col_index[i];
//        this->row_position[i] = mat->row_position[i];
//
//    }
//    this->row_position[nnzs] = nnzs;
//
//}

template <class T>
CSRMatrix<T>::CSRMatrix(const CSRMatrix<T>& mat) :Matrix<T>(mat)
{
    deepcopy(mat);

}


// destructor
template <class T>
CSRMatrix<T>::~CSRMatrix()
{
    // Delete the values array
    if (this->preallocated) {
        delete[] this->row_position;
        delete[] this->col_index;
        //delete[] this->values;
    }
    // The super destructor is called after we finish here
    // This will delete this->values if preallocated is true
}

// Explicitly print out the values in values array as if they are a matrix
template <class T>
void CSRMatrix<T>::printMatrix()
{
    std::cout << "Printing matrix" << std::endl;
    std::cout << "Values: ";
    for (int j = 0; j < this->nnzs; j++)
    {
        std::cout << this->values[j] << " ";
    }
    std::cout << std::endl;
    std::cout << "row_position: ";
    for (int j = 0; j < this->rows + 1; j++)
    {
        std::cout << this->row_position[j] << " ";
    }
    std::cout << std::endl;
    std::cout << "col_index: ";
    for (int j = 0; j < this->nnzs; j++)
    {
        std::cout << this->col_index[j] << " ";
    }
    std::cout << std::endl;
}

//copy function
template <class T>
void CSRMatrix<T>::deepcopy(const CSRMatrix<T>& mat)
{
    // deallocate the value the pointers point to
    if (this->preallocated)
    {
        delete[] this->values;
        delete[] this->row_position;
        delete[] this->col_index;
    }
    this->rows = mat.rows;
    this->cols = mat.cols;
    this->nnzs = mat.nnzs;
    this->preallocated = mat.preallocated;
    this->values = new T[this->nnzs];


    if (mat.values)
    {
        //allocate memory for the copy
        this->values = new T[mat.nnzs];
        this->row_position = new int[mat.rows + 1];
        this->col_index = new int[mat.nnzs];
        // do the copy
        for (int i = 0;i < mat.nnzs;i++)
        {
            this->values[i] = mat.values[i];
            this->col_index[i] = mat.col_index[i];

        }
        for (int j = 0;j < mat.rows + 1;j++)
        {
            this->row_position[j] = mat.row_position[j];
        }

    }

    else
    {
        this->values = nullptr;
        this->col_index = nullptr;
        this->row_position = nullptr;
    }
}

// Do a matrix-vector product
// output = this * input
template<class T>
void CSRMatrix<T>::matVecMult(double* input, double* output)
{
    if (input == nullptr || output == nullptr)
    {
        std::cerr << "Input or output haven't been created" << std::endl;
        return;
    }

    // Set the output to zero
    for (int i = 0; i < this->rows; i++)
    {
        output[i] = 0.0;
    }

    int val_counter = 0;
    // Loop over each row
    for (int i = 0; i < this->rows; i++)
    {
        // Loop over all the entries in this col
        for (int val_index = this->row_position[i]; val_index < this->row_position[i + 1]; val_index++)
        {
            // This is an example of indirect addressing
            // Can make it harder for the compiler to vectorise!
            output[i] += this->values[val_index] * input[this->col_index[val_index]];

        }
    }
}


template <class T>
bool CSRMatrix<T>::judge_SPD(CSRMatrix<T>* U)
{

    bool judge = false;
    //bool judge_symmetry == true;//assume the matrix is symmetrical
    if (this->rows != this->cols)
    {
        cerr << "!Warning! The matrix A isn't a square matrix!" << endl;
        return judge;
    }
    else
    {
        int nii = 0;//used to check if all elements on A's diagonal are nonzero
        int diapositive = 0;//used to check if all elements on A's diagonal are positive 
        bool symmetry = false;//used to check if A is symmetrical
        for (int i = 0;i < this->rows;i++) {
            for (int j = this->row_position[i];j < this->row_position[i + 1];j++)
            {

                int srow = this->col_index[j];//srow represents the column index for the element in the ith row
                //i.e. A[i,srow]=this->value[j];
                //srow is corresponding to the inverted row index

                if (srow == i)
                {
                    nii++;//if A[i,i] is nonzero, nii++
                    if (this->values[j] > 0)
                        diapositive++;//if A[i,i]is positive          
                }

                //loop over the nonzero elements in the [srow]th row
                for (int scol = 0;scol < (this->row_position[srow + 1] - this->row_position[srow]);scol++)
                {
                    //to judge if A[i,srow]==A[srow,i])
                    //pick the location of A[srow,i]
                    if (this->col_index[this->row_position[srow] + scol] == i)
                    {
                        if (this->values[j] == this->values[this->row_position[srow] + scol])//symmetrical 
                            symmetry = true;
                    }

                }


            }
        }
        if ((nii != this->rows) || (diapositive != this->rows) || (!symmetry))
        {
            cerr << "!Warning! The matrix A isn't a positive definete matrix!" << endl;
            return judge;
        }
        else
        {
            //use cholesky decomposition to judge if matrix A is positive definete
            //the upper triangle matrix U
            // assume all the elements in the upper area are nonzero
            //i.e. nnzs=rows*(rows+1)/2
            
            // initialize U matrix to zero

            for (int i = 0; i < this->rows + 1; i++)
            {
                U->row_position[i] = (this->rows * 2 - i + 1) * i / 2;//only consider the upper triangle value
            }
            for (int i = 0; i < this->rows; i++)
            {
                int cnt = i;
                for (int val_index = U->row_position[i]; val_index < U->row_position[i + 1]; val_index++)
                {
                    U->values[val_index] = 0;
                    U->col_index[val_index] = cnt;
                    cnt++;

                }
            }

            U->values[0] = sqrt(this->values[0]);//U[0,0]=sqrtA[0,0]

            for (int i = 0;i < this->row_position[1]; i++)
            {

                U->values[this->col_index[i]] = this->values[i] / U->values[0];
                //U[0,i]=A[0,i]/U[0,0]
                //the first row of U
            }
            //U->printMatrix();
            for (int i = 1; i < this->rows; i++)
            {
                // starts from the second row

                //compute U[i,i]
                T subsum = 0;
                int k = 0;
                do {
                    subsum += pow(U->values[k * (this->rows * 2 - k + 1) / 2 + i-k], 2);//sum(U[k,i]^2),(k=1 to i-1)
                    k++;
                } while (k == i - 1);

                // get the diagnal element Aii=A[i, i]
                T Aii = 0;
                for (int val_index = this->row_position[i]; val_index < this->row_position[i + 1]; val_index++)
                {
                    if (this->col_index[val_index] == i)
                    {
                        Aii = this->values[val_index];
                    }
                }

                T sum = Aii - subsum;//A[i,i]-sum(U[k,i]^2),(k=1 to i-1)
                if (sum <= 0)
                {
                    cerr << "This matrix A isn't positive definete" << endl;
                    return judge;
                }
                else
                {
                    judge = true;
                    U->values[i * (this->rows * 2 - i + 1) / 2] = sqrt(sum);//U[i,i]

                    //to compute U[i,j] j>i
                    for (int j = i + 1; j < this->rows; j++)
                    {
                        T ssubsum = 0;
                        int kk = 0;
                        do {
                            //sum(U[kk,i]*U[kk,j]),(kk=1 to i-1)

                            ssubsum += U->values[kk * (this->rows * 2 - kk + 1) / 2 + i -kk] * U->values[kk * (this->rows * 2 - kk + 1) / 2 + j -kk];
                            kk++;
                        } while (kk == i - 1);

                        // get the element Aij=A[i, j]
                        T Aij = 0;
                        for (int val_index = this->row_position[i]; val_index < this->row_position[i + 1]; val_index++)
                        {
                            if (this->col_index[val_index] == j)
                            {
                                Aij = this->values[val_index];
                            }
                        }
                        T ssum = Aij - ssubsum;
                        U->values[i * (this->rows * 2 - i + 1) / 2 + j - i] = ssum / sqrt(sum);//U[i,j]


                    }//derive the Matrix U 



                }
            }
            cout << "This matrix A is positive definete and could be Cholesky decomposition" << endl;

        }

    }
    return judge;
}

template <class T>
void CSRMatrix<T>::CSRtranspose(CSRMatrix<T>& transposed)
{
    if ((this->rows != transposed.cols) || (this->cols != transposed.rows) || (this->nnzs != transposed.nnzs))
    {
        cerr << "The transposed matrix needs to be resized" << endl;
        return;
    }
    int* p = new int[this->cols];//an array used to store the times of appearance
    // of transposed.row_position during the loop
    for (int N = 0;N < this->cols;N++)
    {
        p[N] = 0;
    }//initialize 

    for (int i = 0;i < this->cols + 1;i++)
    {
        transposed.row_position[i] = 0;


    }

    for (int i = 0;i < this->rows;i++)
    {
        for (int j = this->row_position[i];j < this->row_position[i + 1];j++)
        {
            int srow = this->col_index[j];
            //value[j]=A[i,srow];
            //the [srow]th column
            //count the elements at the [srow]th row of the inverted matrix
            (transposed.row_position[srow + 1])++;
        }

    }

    for (int j = 1;j < this->cols + 1;j++)
    {
        transposed.row_position[j] += transposed.row_position[j - 1];

    }//get new transposed.row_position

    for (int i = 0;i < this->rows;i++)
    {
        for (int j = this->row_position[i];j < this->row_position[i + 1];j++)
        {
            int srow = this->col_index[j];//
            //this->value[j]=A[i,srow];
            //the [srow]th column
            if (transposed.row_position[srow + 1] - transposed.row_position[srow])
            {

                //transposed[srow,i]=A[i,srow]=this->values[j]

                transposed.values[transposed.row_position[srow] + p[srow]] = this->values[j];
                //transposed.values[j] = this->values[this->row_position[srow] + i];
                transposed.col_index[transposed.row_position[srow] + p[srow]] = i;
                p[srow]++;
            }
          
            
            //transposed[srow,i]=A[i,srow]=this->values[j]

            
        }

    }
    delete[]p;
}

template <class T>

void CSRMatrix<T>::Chol_sparse(double* b, double* output)
{
    cout << "The is Cholesky Decomposition method:" << endl;

    // Set values to zero before hand
    for (int i = 0; i < this->rows; i++)
    {
        output[i] = 0;
    }

    //to create a matrix to store the upper triangle matrix
    auto* U = new CSRMatrix<double>(this->rows, this->rows, this->rows * (this->rows + 1) / 2, true);

    // initialize U matrix to zero

    for (int i = 0; i < this->rows + 1; i++)
    {
        U->row_position[i] = (this->rows * 2 - i + 1) * i / 2;//only consider the upper triangle value
    }
    for (int i = 0; i < this->rows; i++)
    {
        int cnt = i;
        for (int val_index = U->row_position[i]; val_index < U->row_position[i + 1]; val_index++)
        {
            U->values[val_index] = 0;
            U->col_index[val_index] = cnt;
            cnt++;

        }
    }

    this->judge_SPD(U);//derive the upper matrix
    
    
    CSRMatrix<double>L(*U);

    U->CSRtranspose(L);//transpose U to get L
    double* Y = new double[this->rows];
    
    forward_subtitution(L, Y, b);
    U->back_subtitution(Y, output);
    
   
    for (int i = 0;i < this->rows;i++)
    {
        cout << "X" << "[" << i + 1 << "]" <<"="<<output[i]<< endl;
    }
    
    delete U;
    U = nullptr;
    delete[]Y;
}



template <class T>

void CSRMatrix<T>::LU_sparse(double* b, double* output)
{
    // Set values to zero before hand
    for (int i = 0; i < this->rows; i++)
    {
        output[i] = 0;
    }

    // initiate partial pivoting matrix (identity matrix)
    auto* P = new CSRMatrix<T>(this->rows, this->cols, this->rows, true);

    for (int i = 0; i < this->rows; i++)
    {
        P->values[i] = 1;
        P->col_index[i] = i;
        P->row_position[i] = i;
    }
    P->row_position[this->rows] = this->rows;

    auto* L = new CSRMatrix<T>(this->rows, this->cols, this->rows * this->cols, true);

    // initialize L matrix to zero
    for (int i = 0; i < this->rows + 1; i++)
    {
        L->row_position[i] = this->rows * i;
    }
    for (int i = 0; i < this->rows; i++)
    {
        int cnt = 0;
        for (int val_index = L->row_position[i]; val_index < L->row_position[i + 1]; val_index++)
        {
            L->values[val_index] = 0;
            L->col_index[val_index] = cnt;
            cnt++;
        }
    }

    // cout << "check L:" << endl;
     //L->printMatrix();
     //cout << "check P:" << endl;
     //P->printMatrix();

    for (int k = 0; k < this->rows - 1; k++)
    {
        //find the index of maximum magnitude
        int j = find_max_index(k);
        swapRows(j, k); // do the swap of A's rows
        P->swapRows(j, k); // record the swap in the matrix P
        L->swapRows(j, k);
       

        // get the diagnal element A[k, k]
        T Akk = 0;
        for (int val_index = this->row_position[k]; val_index < this->row_position[k + 1]; val_index++)
        {
            if (this->col_index[val_index] == k)
            {
                Akk = this->values[val_index];
            }
        }


        for (int i = k + 1; i < this->rows; i++)
        {
            // get A[i,k]
            T Aik = 0;
            for (int val_index = this->row_position[i]; val_index < this->row_position[i + 1]; val_index++)
            {
                if (this->col_index[val_index] == k)
                {
                    Aik = this->values[val_index];
                }
            }
            T s = Aik / Akk;

            // vectorises the row update from before for j in range(k,n)
            for (int val_index = this->row_position[i]; val_index < this->row_position[i + 1]; val_index++)
            {
                int colIndex = this->col_index[val_index];
                // for A[i,colIndex]
                if (colIndex >= k)
                {
                    // get A[k, colIndex]
                    T AcolIndex = 0;
                    for (int val_ii = this->row_position[k]; val_ii < this->row_position[k + 1]; val_ii++)
                    {
                        if (this->col_index[val_ii] == colIndex)
                        {
                            AcolIndex = this->values[val_ii];
                        }
                    }
                    this->values[val_index] -= AcolIndex * s;

                }

            }
            L->values[i * this->cols + k] = s;
        }

    }

    // add in the diagonal of ones into L
    for (int i = 0; i < this->rows; i++)
    {
        L->values[i * (this->rows + 1)] += 1;
    }

   

    // update mat_b with P*mat_b
    double* matpb = new double[this->rows];
    P->matVecMult(b, matpb);

    double* Y = new double[this->rows];
    forward_subtitution(*L, Y, matpb);
    back_subtitution(Y, output);

    for (int i = 0;i < this->rows;i++)
    {
        cout << "X" << "[" << i + 1 << "]" << "=" << output[i] << endl;
    }

}

// although CSR matrix reduces the memory, the number of computations increase for swap 
template <class T>

void CSRMatrix<T>::swapRows(int j, int k)
{
    // make j < k
    if (j >= k)
    {
        int temp = j;
        j = k;
        k = temp;
    }
    // get the updated rowIndex array
    int* temp_row_index = new int[this->rows + 1];
    for (int i = 0; i < this->rows + 1; i++)
    {
        temp_row_index[i] = this->row_position[i];
    }

    // update rowIndex array
    this->row_position[j + 1] = this->row_position[j] + (temp_row_index[k + 1] - temp_row_index[k]);

    for (int i = j + 2; i < k + 1; i++)
    {
        this->row_position[i] = this->row_position[i - 1] + (temp_row_index[i] - temp_row_index[i - 1]);
    }
    // this->row_position[j + 2] = this->row_position[j + 1] + temp_row_index[j + 2] - temp_row_index[j + 1];

    this->row_position[k + 1] = this->row_position[k] + temp_row_index[j + 1] - temp_row_index[j];

    // copy the values between row j and row k
    T* temp_values = new T[this->nnzs];

    // copy the col index between row j and row k
    T* temp_col = new T[this->nnzs];

    for (int i = 0; i < this->nnzs; i++)
    {
        temp_values[i] = this->values[i];
        temp_col[i] = this->col_index[i];
    }



    int cnt2 = 0;

    for (int val_index = this->row_position[j]; val_index < this->row_position[j + 1]; val_index++)
    {
        this->values[val_index] = temp_values[temp_row_index[k] + cnt2];
        this->col_index[val_index] = temp_col[temp_row_index[k] + cnt2];
        cnt2++;
    }
    cnt2 = 0;
    for (int val_index = this->row_position[k]; val_index < this->row_position[k + 1]; val_index++)
    {
        this->values[val_index] = temp_values[temp_row_index[j] + cnt2];
        this->col_index[val_index] = temp_col[temp_row_index[j] + cnt2];
        cnt2++;
    }

    for (int i = j + 1; i < k; i++)
    {
        int cnt3 = 0;
        for (int val_index = this->row_position[i]; val_index < this->row_position[i + 1]; val_index++)
        {
            // repopulate the value array between row j and row k (j, k not included)

            this->values[val_index] = temp_values[temp_row_index[i] + cnt3];
            // repopulate the col index array between row j and row k (j, k not included)
            this->col_index[val_index] = temp_col[temp_row_index[i] + cnt3];
            cnt3++;
        }
    }
}

template <class T>
int CSRMatrix<T>::find_max_index(int k)
{
    int maxIndex = k;
    T maxValue = 0;
    for (int i = k; i < this->rows; i++)
    {
        for (int val_index = this->row_position[i]; val_index < this->row_position[i + 1]; val_index++)
        {
            if (this->col_index[val_index] == k)
            {
                if (abs(this->values[val_index]) > abs(maxValue))
                {
                    maxIndex = i;
                    maxValue = this->values[val_index];
                }
            }
        }
    }
    return maxIndex;
}

template <class T>
void CSRMatrix<T>::back_subtitution(T* Y, T* output)
{   // store diagonal element in an array
    T* diag = new T[this->rows];

    for (int i = 0; i < this->rows; i++)
    {
        for (int val_index = this->row_position[i]; val_index < this->row_position[i + 1]; val_index++)
        {
            if (this->col_index[val_index] == i)
            {
                diag[i] = this->values[val_index];
            }
        }
    }

    output[this->rows - 1] = Y[this->rows - 1] / this->values[this->nnzs - 1];
    for (int k = this->rows - 2; k > -1; k--)
    {
        T s = 0;
        for (int val_index = this->row_position[k]; val_index < this->row_position[k + 1]; val_index++)
        {
            if (this->col_index[val_index] > k)
            {
                s = s + this->values[val_index] * output[this->col_index[val_index]];
            }
        }
        output[k] = (Y[k] - s) / diag[k];
    }
}

template <class T>

void CSRMatrix<T>::forward_subtitution(CSRMatrix<T>& L, T* Y, T* b)
{
  


    T* diag = new T[this->rows];

    for (int i = 0; i < L.rows; i++)
    {
        for (int val_index = L.row_position[i]; val_index < L.row_position[i + 1]; val_index++)
        {
            if (L.col_index[val_index] == i)
            {
                diag[i] = L.values[val_index];
               
            }
        }
    }

    Y[0] = b[0] / diag[0];

    for (int i = 1; i < L.rows; i++)
    {
        T s = 0;
        for (int val_index = L.row_position[i]; val_index < L.row_position[i + 1]; val_index++)
        {
            if (L.col_index[val_index] < i)
            {
                s = s + L.values[val_index] * Y[L.col_index[val_index]];
            }
        }
        Y[i] = (b[i] - s) / diag[i];
    }

    
}


//Find the maximum value of the array for each iteration
template<class T>
T CSRMatrix<T>::MAX(T* a)
{
    T max = 0;
    for (int i = 0; i < this->rows; i++)
    {
        if (fabs(a[i]) > max)
            max = fabs(a[i]);
    }
    return max;
}


template<class T>
void CSRMatrix<T>::Jacobi(T* mat_b, T* output, T* initial)
{
    // to store recursive times
    int count = 0;
    // to store diagonal elements
    T* diag = new T[this->rows];
    // to store the deviation between iterations
    T deviation;
    // check when the iteration stops
    bool flag;

    // Set values to zero before hand
    for (int i = 0; i < this->rows; i++)
    {
        output[i] = 0;
        diag[i] = 0;
    }

    do {
        count++;
        // store the results of the last iteration
        for (int j = 0; j < this->rows; j++)
            initial[j] = output[j];

        // loop over each row
        for (int i = 0; i < this->rows; i++)
        {
            T sum = 0;
            for (int val_index = this->row_position[i]; val_index < this->row_position[i + 1]; val_index++)
            {
                // when (row = col), store the diagonal elements 
                if (i == this->col_index[val_index])
                {
                    diag[i] = this->values[val_index];
                    continue;
                }
                // when (row != col), do some multiplication and add them altogether
                else
                {
                    sum += this->values[val_index] * initial[this->col_index[val_index]];
                }
            }
            output[i] = (mat_b[i] - sum) / diag[i];
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

    cout << "This is Jacobi method  " << endl;
    for (int i = 0; i < this->cols; i++)
        cout << "x[" << i + 1 << "]= " << output[i] << endl;
    cout << "recursive times:  " << count << endl;

    delete[] diag;
}

template<class T>
void CSRMatrix<T>::Gauss_Seidel(T* mat_b, T* output, T* initial)
{
    int count = 0;

    T* diag = new T[this->rows];
    T deviation;
    bool flag;

    // Set values to zero before hand
    for (int i = 0; i < this->rows; i++)
    {
        output[i] = 0;
        diag[i] = 0;
    }

    do {
        // to store resursive times
        count++;
        for (int j = 0; j < this->rows; j++)
            initial[j] = output[j];

        for (int i = 0; i < this->rows; i++)
        {
            T sum1 = 0;
            T sum2 = 0;
            for (int val_index = this->row_position[i]; val_index < this->row_position[i + 1]; val_index++)
            {
                if (i > this->col_index[val_index])
                {
                    sum1 += this->values[val_index] * output[this->col_index[val_index]];
                }
                else if (i == this->col_index[val_index])
                {
                    diag[i] = this->values[val_index];
                    continue;
                }
                else
                {
                    sum2 += this->values[val_index] * initial[this->col_index[val_index]];
                }
            }
            output[i] = (mat_b[i] - sum1 - sum2) / diag[i];
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

    cout << "this is Gauss Seidel method  " << endl;
    for (int i = 0; i < this->cols; i++)
        cout << "x[" << i + 1 << "]= " << output[i] << endl;
    cout << "recursive times:  " << count << endl;

    delete[] diag;
}