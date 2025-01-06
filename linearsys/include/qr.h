/*
 * Copyright (c) 2024 Ismael Mosquera Rivera
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

#ifndef ___QR_H___
#define ___QR_H___

#ifdef __cplusplus
extern "C" {
	#endif

#include "matrix.h"

/*
* QR type definition.
*/
typedef struct
{
Matrix* _q;
Matrix* _r;
}QR;

/*
* Releases the memory previously allocated for the QR type passed as parameter
* param: QR* qr => a pointer to a QR type
*/
void destroy_qr(QR* qr);

/*
* Computes a QR factorization for a matrix passed as parameter,
* using Givens Rotations method.
* The result must be an orthogonal matrix => Q and a upper matrix => R,
* where m = QR.
* param: Matrix* m -> a matrix to be factorized.
*
* returns: a QR factorization for m or NULL if the operation cannot be done.
*/
QR* qr_factorization(const Matrix* m);

/*
* Prints a QR to the console.
* param: const QR* qr => a pointer to a QR
*
* The output format is as follows:
*
* q:
*[ ... ]
* [ ... ]
* ...
* [ ... ]
*
* r:
*[ ... ]
* [ ... ]
* ...
* [ ... ]
*
*/
void print_qr(const QR* qr);

/*
* Computes the determinant of a matrix using QR factorization.
* param: QR* qr => a pointer to an already computed QR.
*
* returns: determinant.
*/
double qr_det(const QR* qr);

/*
* Macros to get Q and R matrices.
*/
#define qr_q(qr) ((qr)->_q)
#define qr_r(qr) ((qr)->_r)

#ifdef __cplusplus
}
#endif

#endif
