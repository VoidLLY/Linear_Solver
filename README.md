This software implements algorithms to solve the linear system Ax=b, where A is a positive definite matrix, and x and b. Methods include Jacobi, Guass Seidel, SOR, LU decomposition and Cholesky factorisation.

### Installation Guide

To install the tool, run the following code on a machine that has Git installed on it:

```
git clone https://github.com/VoidLLY/Linear_Solver.git
```
### User instructions

If user wants to solve a dense matrix:

write the values in data.txt file and make sure the file is in the same folder of the source code. Then run the program.
An example file is as follows:
```
3 1 1 
1 4 1 
1 1 5
2 3 6
```
The above file solves a 3x3 matrix with b = 2, 3, 6.

If user wants to solve a sparse matrix:

write the values in sparse_data.txt file and make sure the file is in the same folder of the source code. Then run the program.
The input of the file should follows the order: number of non-zeros, size, values, column index of values, index of the first non-zero of each row, and values of b matrix.
An example file is as follows:
```
10
10
1 2 3 4 5 6 7 8 9 10
0 1 2 3 4 5 6 7 8 9
0 1 2 3 4 5 6 7 8 9 10
6 7 10 11 14 18 -1 9 12 7 
```
The above file solves a 10x10 diagonal matrix of which values are 1,2,3,4,5,6,7,8,9,10 and b values are 6,7,10,11,14,18,-1,9,12,7.
### Testing

The program uses a random generator which ensures that the matrix is positive definite, symmetric and diagonal dominant, to feed values into the solvers. The results with computation time of each solvers are printed to the screen so that we could compare the solutions and performance. Testing function is included in main.cpp file with name test_dense(). If you want to run the test_dense() function, uncomment the code test_dense() in main() function.
