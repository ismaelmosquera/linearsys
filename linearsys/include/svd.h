/*
 * Copyright (c) 2025 Ismael Mosquera Rivera
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

#ifndef ___SVD_H___
#define ___SVD_H___

#ifdef __cplusplus
extern "C" {
	#endif

#include "matrix.h"

/*
* SVD type definition.
*
* A SVD ( singular Value Decomposition ) is a factorization for a MxN matrix so that:
* M = UDV^t
* were
* U is an orthogonal matrix.
* D ( or sigma ) is a diagonal matrix.
* V is an orthogonal matrix.
*
*/
typedef struct
{
	Matrix* _u; /* U matrix */
	Matrix* _sigma; /* D ( or sigma ) matrix */
	Matrix* _v; /* V matrix */
}SVD;

/*
* Destroys a SVD.
* param: svd SVD to be destroyed.
*
*/
void destroy_svd(SVD* svd);

/*
* Prints a SVD to the console.
* param: svd SVD to print.
*
*/
void print_svd(const SVD* svd);

/*
* Computes a SVD factorization for a MxN matrix.
* param: m Matrix to be factorized.
*
* The result of this factorization is as follows:
* U => an orthogonal matrix having the left singular vectors of m ( matrix passed as parameter ).
* Sigma => a diagonal matrix having the singular values of m ( matrix passed as parameter ).
* V an orthogonal matrix having the right singular vectors of m ( matrix passed as parameter ).
*
* so that:
* M = USigmaV^t
*
* in case that the above expression is not true,
* it means that the factorization for the matrix passed as parameter cannot be done properly.
*
* returns: a SVD factorization for the matrix passed as parameter.
*/
SVD* svd_factorization(const Matrix* m);

/*
* Macros to access SVD data members.
*/
#define svd_u(svd) ((svd)->_u)
#define svd_sigma(svd) ((svd)->_sigma)
#define svd_v(svd) ((svd)->_v)

#ifdef __cplusplus
}
#endif

#endif
