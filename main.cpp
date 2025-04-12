#include <iostream>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <stdlib.h>
#include "Matrix.h"
#include "CSRMatrix.h"
#include "Matrix.cpp"
#include "CSRMatrix.cpp"
#include <fstream>
#include <string>
#include<iostream>
#include<fstream>
#include<stdlib.h>
#include<string.h>
#include <typeinfo>
#include <omp.h>


using namespace std;
void ReadLineData(const char* fileName, int lineNum, char* data);
int split_double(double* arr, char* str);
int split_int(int* arr, char* str);
int CountLines(char* filename);
int user_dense();
int user_sparse();
int test_dense();

int main()
{
	user_dense();
	//test_dense();
	user_sparse();
	system("pause");
	return 0;
}


void ReadLineData(const char* fileName, int lineNum, char* data)
{
	ifstream in;
	in.open(fileName);

	int line = 0;
	while (in.getline(data, 1024))
	{
		if (lineNum == line)
		{
			break;
		}
		line++;
	}

	in.close();
}



int split_double(double* arr, char* str) {
	char* buf;
	char* s = strtok_s(str, " ", &buf);
	//cout<<"S1";
	int c = 0;
	while (s != NULL) {
		arr[c] = atof(s);
		s = strtok_s(NULL, " ", &buf);
		c++;
	}

	return c;
}

int split_int(int* arr, char* str) {
	char* buf;
	char* s = strtok_s(str, " ", &buf);
	//cout<<"S1";
	int c = 0;
	while (s != NULL) {
		arr[c] = atoi(s);
		s = strtok_s(NULL, " ", &buf);
		c++;
	}

	return c;
}


int CountLines(char* filename)
{
	ifstream ReadFile;
	int count = 0;
	string tmp;
	ReadFile.open(filename, ios::in);
	if (ReadFile.fail())
	{
		return 0;
	}
	else
	{
		while (getline(ReadFile, tmp, '\n'))
		{
			count++;
		}
		ReadFile.close();
		return count;
	}
}

int user_dense()
{
	cout << "Solving dense Matrix..." << endl;
	//int line;
	char filename[] = "data.txt";
	ifstream infile;
	infile.open("data.txt");

	int lines = CountLines(filename);
	int rows = lines - 1;
	int cols = rows;

	double* matrix_values = new double[rows * cols];
	double* b = new double[rows];
	while (!infile.eof())
	{
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				infile >> matrix_values[i * rows + j];
			}
		}

		for (int i = 0; i < rows; i++)
		{
			infile >> b[i];
		}
	}
	infile.close();

	//auto* user_dense_mat = new Matrix<double>(rows, cols, matrix_values);
	unique_ptr<Matrix<double>> user_dense_mat(new Matrix<double>(rows, cols, matrix_values));
	unique_ptr<Matrix<double>> b_mat(new Matrix<double>(rows, 1, b));
	//auto* b_mat = new Matrix<double>(rows, 1, b);
	unique_ptr<Matrix<double>> output(new Matrix<double>(rows, 1, true));
	unique_ptr<Matrix<double>> initial(new Matrix<double>(rows, 1, true));
	// call Jacobi
	double dtime = omp_get_wtime();
	user_dense_mat->Jacobi(*b_mat, *output, *initial);
	dtime = omp_get_wtime() - dtime;
	printf("Jacobi method takes time %f ms\n", dtime * 1000);

	// call Gauss Seidel
	dtime = omp_get_wtime();
	user_dense_mat->Gauss_Seidel(*b_mat, *output, *initial);
	dtime = omp_get_wtime() - dtime;
	printf("Gauss Seidel method takes time %f ms\n", dtime * 1000);

	// call Gauss Seidel SOR
	dtime = omp_get_wtime();
	user_dense_mat->Gauss_Seidel_SOR(*b_mat, *output, *initial, 1);
	dtime = omp_get_wtime() - dtime;
	printf("Gauss Seidel SOR method takes time %f ms\n", dtime * 1000);

	// call Cholesky factorisation
	dtime = omp_get_wtime();
	user_dense_mat->Chol(*b_mat, *output);
	dtime = omp_get_wtime() - dtime;
	printf("Cholesky factorisation  method takes time %f ms\n", dtime * 1000);
	// call LU
	dtime = omp_get_wtime();
	user_dense_mat->LU(*b_mat, *output);
	dtime = omp_get_wtime() - dtime;
	cout << "X =" << endl;

	output->printMatrix(); printf("LU method takes time %f ms\n", dtime * 1000);
	cout << endl;

	delete[] b;
	delete[] matrix_values;

	return 0;
}


int test_dense()
{
	cout << "This test uses random matrix and b" << endl;
	cout << "-------------------------------------------------" << endl;
	int rows = 10;
	int cols = 10;
	//auto* dense_mat = new Matrix<double>(rows, cols, true);
	unique_ptr<Matrix<double>> dense_mat(new Matrix<double>(rows, cols, true));
	srand((unsigned int)time(NULL));
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < i; j++)
		{
			//set random values for elements below the diagonal
			dense_mat->values[i * cols + j] = rand() % 10 + 1;
			//set corresponding values for elements beyond the diagonal, resulting in a symmetric matrix
			dense_mat->values[j * rows + i] = dense_mat->values[i * cols + j];
		}
	}
	//diagonally dominant
	for (int i = 0; i < rows; i++)
	{
		int sum = 0;
		dense_mat->values[i * cols + i] = rand() % 10;
		for (int j = 0; j < cols; j++)
		{
			sum += dense_mat->values[i * cols + j];
		}
		dense_mat->values[i * cols + i] = sum;
	}
	dense_mat->printMatrix();
	unique_ptr<Matrix<double>> b_mat(new Matrix<double>(rows, 1, true));
	srand((unsigned int)time(NULL));
	for (int i = 0; i < rows; i++)
	{
		b_mat->values[i] = rand() % 20 + 3;
	}
	//auto* b_mat = new Matrix<double>(rows, 1, b);
	unique_ptr<Matrix<double>> output(new Matrix<double>(rows, 1, true));
	unique_ptr<Matrix<double>> initial(new Matrix<double>(rows, 1, true));
	// call Jacobi
	double dtime = omp_get_wtime();
	dense_mat->Jacobi(*b_mat, *output, *initial);
	dtime = omp_get_wtime() - dtime;
	printf("Jacobi method takes time %f ms\n", dtime * 1000);
	// call Gauss Seidel
	dtime = omp_get_wtime();
	dense_mat->Gauss_Seidel(*b_mat, *output, *initial);
	dtime = omp_get_wtime() - dtime;
	printf("Gauss Seidel method takes time %f ms\n", dtime * 1000);
	// call Gauss Seidel SOR
	dtime = omp_get_wtime();
	dense_mat->Gauss_Seidel_SOR(*b_mat, *output, *initial, 1);
	dtime = omp_get_wtime() - dtime;
	printf("Gauss Seidel SOR method takes time %f ms\n", dtime * 1000);
	// call Cholesky factorisation
	dtime = omp_get_wtime();
	dense_mat->Chol(*b_mat, *output);
	dtime = omp_get_wtime() - dtime;
	printf("Cholesky factorisation  method takes time %f ms\n", dtime * 1000);
	// call LU
	dtime = omp_get_wtime();
	dense_mat->LU(*b_mat, *output);
	dtime = omp_get_wtime() - dtime;
	cout << "X =" << endl;
	output->printMatrix(); printf("LU method takes time %f ms\n", dtime * 1000);
	cout << endl;
	return 0;
}

int user_sparse()
{
	cout << "Solving sparse Matrix..." << endl;
	//int line;
	char filename[] = "sparse_data.txt";


	char nn_zs[10] = { 0 };
	char n_rows[10] = { 0 };
	char val[1024] = { 0 };
	char row[1024] = { 0 };
	char col[1024] = { 0 };
	char mat_b[1024] = { 0 };

	ReadLineData(filename, 0, nn_zs);
	ReadLineData(filename, 1, n_rows);
	ReadLineData(filename, 2, val);
	ReadLineData(filename, 3, col);
	ReadLineData(filename, 4, row);
	ReadLineData(filename, 5, mat_b);

	//n[0] = number of non-zeros
	int* n = new int[10];
	int s1 = split_int(n, nn_zs);

	//n_r[0] = number of rows
	int* n_r = new int[10];
	int s2 = split_int(n_r, n_rows);

	double* v = new double[n[0]];
	int ss1 = split_double(v, val);
	int* c = new int[n[0]];
	int ss2 = split_int(c, col);
	int* r = new int[n_r[0] + 1];
	int ss3 = split_int(r, row);
	double* b = new double[n_r[0]];
	int ss4 = split_double(b, mat_b);

	int rows = ss3 - 1;//the number of rows
	int nnzs = ss1;//the number of nonzero elements

	//auto* user_dense_mat = new Matrix<double>(rows, cols, matrix_values);
	unique_ptr<CSRMatrix<double>> user_CSR_mat(new CSRMatrix<double>(rows, rows, nnzs, true));
	//unique_ptr <double> b_mat(new double[rows]);
	auto* b_mat = new double[rows];
	//unique_ptr <double> output(new double[rows]);
	auto* output = new double[rows];
	//unique_ptr <double> initial(new double[rows]);
	auto* initial = new double[rows];


	for (int i = 0; i < ss1; i++)
	{
		user_CSR_mat->values[i] = v[i];
	}

	for (int i = 0; i < ss2; i++)
	{
		user_CSR_mat->col_index[i] = c[i];
	}

	for (int i = 0; i < ss3; i++)
	{
		user_CSR_mat->row_position[i] = r[i];
	}

	for (int i = 0; i < rows; i++)
	{
		b_mat[i] = b[i];
	}

	// call Jacobi
	double dtime = omp_get_wtime();
	user_CSR_mat->Jacobi(b_mat, output, initial);
	dtime = omp_get_wtime() - dtime;
	printf("Jacobi method takes time %f ms\n", dtime * 1000);

	// call Gauss Seidel
	dtime = omp_get_wtime();
	user_CSR_mat->Gauss_Seidel(b_mat, output, initial);
	dtime = omp_get_wtime() - dtime;
	printf("Gauss Seidel method takes time %f ms\n", dtime * 1000);



	//// call Cholesky factorisation
	//dtime = omp_get_wtime();
	//user_CSR_mat->Chol_sparse(b_mat, output);
	//dtime = omp_get_wtime() - dtime;
	//printf("Cholesky factorisation  method takes time %f ms\n", dtime * 1000);
	// call LU
	dtime = omp_get_wtime();
	user_CSR_mat->LU_sparse(b_mat, output);
	dtime = omp_get_wtime() - dtime;


	printf("LU method takes time %f ms\n", dtime * 1000);
	cout << endl;

	delete[]b_mat;
	delete[]output;
	delete[]initial;
	delete[]v;
	delete[]c;
	delete[]r;
	delete[]b;
	return 0;
}