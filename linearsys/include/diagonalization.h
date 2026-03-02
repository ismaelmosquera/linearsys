/*
 * Copyright (c) 2026 Ismael Mosquera Rivera
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

#ifndef ___MATRIX_DIAGONALIZATION___
#define ___MATRIX_DIAGONALIZATION___

#ifdef __cplusplus
extern "C" {
	#endif

#include "matrix.h"

/*
* This header has functions suitable to perform matrix diagonalization.
* Matrix diagonalization can be done by finding a P orthogonal matrix and a D diagonal matrix so that:
* A = PDP^-1
* where A is the matrix to diagonalize.
* Take in account that maybe it could not be possible to find such a P and D matrices.
* In that case, we say that the matrix A is defective.
*/

/*
* Diagonalization type definition.
*/
typedef struct
{
	Matrix* _p; /* Orthogonal matrix P */
	Matrix* _d; /* Diagonal matrix D */
	Matrix* _pt; /* Inverse of P */
}Diagonalization;


/*
* Destroys a Diagonalization previously allocated structure.
* param: diag
* A pointer to a previously allocated Diagonalization structure.
*/
void destroy_diagonalization(Diagonalization* diag);

/*
* Makes and returns a clone of the Diagonalization passed as parameter.
* param: diag
* A Diagonalization to be cloned.
*
* returns: A clone of the Diagonalization passed as parameter or NULL if the operation cannot be done.
*
*/
Diagonalization* clone_diagonalization(const Diagonalization* diag);

/*
* Performs the diagonalization of the matrix passed as parameter.
* param: m
* A square matrix to diagonalize.
*
* returns: Diagonalization of the matrix passed as parameter or NULL if the operation cannot be done.
*
*/
Diagonalization* diagonalize_matrix(const Matrix* m);

/*
* Computes the determinat of a previously diagonalized square matrix.
* param: diag A matrix diagonalization.
*
* returns: Determinant of the diagonalized matrix or NaN if the operation cannot be done.
*
*/
double det_diagonalized_matrix(const Diagonalization* diag);

/*
* Prints to the console the Diagonalization passed as parameter.
* param diag
* A pointer to a Diagonalization structure.
*
*/
void print_diagonalization(const Diagonalization* diag);


/*
* Macros to access the matrices from a Diagonalization structure.
*/
#define diagonalization_p(diag) ((diag)->_p)
#define diagonalization_d(diag) ((diag)->_d)
#define diagonalization_pt(diag) ((diag)->_pt)

#ifdef __cplusplus
}
#endif

#endif
