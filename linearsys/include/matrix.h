/*
 * Copyright (c) 2021-2022 Ismael Mosquera Rivera
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#ifndef ___MATRIX_H___
#define ___MATRIX_H___

#ifdef __cplusplus
extern "C" {
	#endif

#include "vector.h"

/*
* Matrix type.
*/
typedef struct
{
int _rows;
int _columns;
double* _data;
}Matrix;

/*
* Creates a nnew Matrix with all his components initialized to zero.
* param: int r => number of rows.
* param: int c => number of columns.
*
* returns:
* A pointer to a new created Matrix.
*/
Matrix* create_matrix(int r, int c);

/*
* Destroy a Matrix.
* Release the memory previously allocated.
* param: Matrix* m => Matrix to destroy.
*/
void destroy_matrix(Matrix* m);

/*
* Makes a clone of the Matrix passed as parameter.
* param: const Matrix* m => pointer to a Matrix to clone.
*
* returns:
* A pointer to a clone of the Matrix.
*/
Matrix* clone_matrix(const Matrix* m);

/*
* Loads a Matrix from a file.
* The format of the stored Matrix is well defined.
*
* #r #c
* a00 a01 a02 ... a0n
* a10 a11 a12 ... a1n
* ...
* am0 am1 am2 ... amn
*
* Where:
* #r => number of roes.
* #c => number of columns.
* a00..a0n
* a10..a1n
* ...
* am0..amn
* are the components of the Matrix.
* m => r-1
* n => c-1
*
* param: const char* filename => file where the Matrix is stored.
*
* returns:
* A pointer to the loaded Matrix.
*/
Matrix* load_matrix(const char* filename);

/*
* Stores a Matrix in a file.
* The format of the Matrix is the same as explained in the above function.
*
* param: const Matrix m => Pointer to a Matrix to store.
* param: const char* filename => name of the file to store the Matrix.
*/
void store_matrix(const Matrix* m, const char* filename);

/*
* Prints a Matrix to the console.
* The format to print the Matrix is as follows:
* [a00, a01, a02, ... a0n]
* [a10, a11, a12, ... a1n]
* ...
* [am0, am1, am2, ... amn]
*
* m => r-1
* n => c-1
*
* param: const Matrix* m => pointer to a Matrix to print to the console.
*/
void print_matrix(const Matrix* m);

/*
* Gets a value from a Matrix.
* param: const Matrix* m => a pointer to a Matrix
* param: int i => index i ( row )
* param: int j => index j ( column )
*
* returns:
* value in matrix(i, j)
*/
double get_matrix(const Matrix* m, int i, int j);

/*
* Sets a value to a Matrix (i, j)
* param: Matrix* m => a pointer to a Matrix
* param: double value => value to set
* param: int i => index i ( row )
* param: int j => index j ( column )
*/
void set_matrix(Matrix* m, double value, int i, int j);

/*
* Gets a chunk from a matrix.
* param: const Matrix* m => Matrix to get a chunk.
* param: int from_row => first row.
* param: int to_row => last row.
* param: int from_column => first column.
* param: int to_column => last column.
* return a matrix containing the requested chunk as content.
* If the operation cannot be performed this function returns NULL.
*/
Matrix* get_matrix_chunk(const Matrix* m, int from_row, int to_row, int from_column, int to_column);

/*
* Resizes the matrix to the number of rows and columns passed as parameter.
* param: const Matrix* m => Matrix to be resized.
* param: int nrows => number of rows.
* param: int ncolumns => number of columns.
*
* returns:
* a resized Matrix.
* If the operation cannot be performed this function returns NULL.
*/
Matrix* resize_matrix(const Matrix* m, int nrows, int ncolumns);

/*
* Resizes only the rows of a matrix.
* param: const Matrix*m => a Matrix to resize its rows.
* param: int nrows => number of rows.
*
* returns:
* a resized Matrix.
* If the operation cannot be performed this function returns NULL.
*/
Matrix* resize_matrix_rows(const Matrix* m, int nrows);

/*
* Resizes only the columns of a matrix.
* param: const Matrix* m => a Matrix to resize.
* param: int ncolumns => number of columns.
*
* returns:
* a resized Matrix.
* If the operation cannot be performed this function returns NULL.
*/
Matrix* resize_matrix_columns(const Matrix* m, int ncolumns);

/*
* Matrix addition.
* param: const Matrix* m1 => a pointer to a matrix.
* param: const Matrix* m2 => a pointer to a matrix.
* The number of rows and columns for m1 must be the same as for m2.
*
* returns:
* A pointer to a matrix = m1 + m2.
*/
Matrix* add_matrix(const Matrix* m1, const Matrix* m2);

/*
* Matrix substraction.
* param: const Matrix* m1 => a pointer to a matrix.
* param: const Matrix* m2 => a pointer to a matrix.
* The number of rows and columns for m1 must be the same as for m2.
*
* returns:
* A pointer to a matrix = m1 - m2.
*/
Matrix* sub_matrix(const Matrix* m1, const Matrix* m2);

/*
* Matrix product.
* param: const Matrix* m1 => a pointer to a matrix.
* param: const Matrix* m2 => a pointer to a matrix.
* columns m1 must be equal to rows m2.
*
* returns:
* A pointer to a matrix = m1 * m2.
* The returned matrix will have m1 rows and m2 columns.
*/
Matrix* mul_matrix(const Matrix* m1, const Matrix* m2);

/*
* Scales a matrix by a scalar passed as second parameter.
* param: const Matrix* m A matrix to be scaled.
* param: double value A scalar to scale the matrix.
*
* returns: A scaled matrix.
*/
Matrix* scale_matrix(const Matrix* m, double value);

/*
* Build an identity matrix.
* The identity matrix has 1s in the diagonal and zeros in the rest.
* param: int n => order of the matrix, the dimension wil be n*n ( n rows and n columns ).
*
* returns:
* A pointer to an identity matrix of order n.
*/
Matrix* identity_matrix(int n);

/*
* Gets the transpose of a matrix.
* param: const Matrix* m => a pointer to a matrix to transpose.
*
* returns:
* A pointer to a transposed matrix.
* The returned matrix will have:
* number of rows = m columns.
* number of columns = m rows.
* where m is the matrix passed as parameter to be transposed.
*/
Matrix* transpose_matrix(const Matrix* m);

/*
* Gets the trace of a square matrix.
* The matrix passed as parameter must be square.
* The trace of a matrix is the product of its diagonal entries.
* param: const Matrix* m A matrix to trace.
*
* returns: the trace of the matrix passed as parameter.
*
*/
double trace_matrix(const Matrix* m);

/*
* Builds a matrix composed of the elements of the diagonal of the matrix passed as parameter.
* The matrix passed as parameter must ve square.
* param: const Matrix*m A matrix from to gets the elements of its diagonal.
*
* returns: A matrix with its diagonal composed from the matrix passed as parameter.
*
*/
Matrix* diag_matrix(const Matrix* m);

/*
* Gets the inverse of the matrix passed as parameter.
* param: const Matrix* m => a pointer to a square matrix.
*
* returns:
* inverse matrix.
*/
Matrix* inv_matrix(const Matrix* m);

/*
* Gets the determinant of a matrix.
* param: const Matrix* m => a pointer to a matrix.
* The matrix passed as parameter must have the same number of rows than columns.
* That is, it must be a square matrix.
*
* returns:
* The determinant of the matrix passed as parameter.
*/
double det_matrix(const Matrix* m);

/*
* Eulers's rotation matrix.
* param x the x component.
* param y the y component.
* param z the z component.
* param w angle.
* return a rotation matrix.
* x, y and z are the components of a vector and w is the rotation angle.
*/
Matrix* rotation_matrix(double x,double y,double z,double w);

/*
* Gets the ith row from the matrix passed as parameter as a matrix.
* param: m Matrix from get a row.
* param: i Index of the row.
*
* returns: ith row of the matrix passed as parameter or NULL if the operation canot be done.
*/
Matrix* get_row_matrix(const Matrix* m, int i);

/*
* Gets the jth column from the matrix passed as parameter as a matrix.
* param: m Matrix from get a column.
* param: j Index of the column.
*
* returns: jth column of the matrix passed as parameter or NULL if the operation canot be done.
*/
Matrix* get_column_matrix(const Matrix* m, int j);

/*
* Gets the ith row from the matrix passed as parameter as a vector.
* param: m Matrix from get a row.
* param: i Index of the row.
*
* returns: ith row of the matrix passed as parameter or NULL if the operation canot be done.
*/
Vector* get_row_vector(const Matrix* m, int i);

/*
* Gets the jth column from the matrix passed as parameter as a vector.
* param: m Matrix from get a column.
* param: j Index of the column.
*
* returns: jth column of the matrix passed as parameter or NULL if the operation canot be done.
*/
Vector* get_column_vector(const Matrix* m, int j);

/*
* Checks matrix orthogonality.
* param: m Matrix to check for orthogonality.
*
 returns: 1 if the matrix passed as parameter is orthogonal or 0 otherwise.
 */
 int is_orthogonal_matrix(const Matrix* m);

/*
* Checks matrix orthonormality.
* param: m Matrix to check for orthonormality.
*
 returns: 1 if the matrix passed as parameter is orthonormal or 0 otherwise.
 */
 int is_orthonormal_matrix(const Matrix* m);

/*
* Macros to get the number of rows and columns from a Matrix.
*/
#define rows_matrix(matrix) ((matrix)->_rows)
#define columns_matrix(matrix) ((matrix)->_columns)

/*
* Macro to get the dimension of a matrix.
*/
#define dimension_matrix(matrix) (rows_matrix(matrix) * columns_matrix(matrix))

#ifdef __cplusplus
}
#endif

#endif
