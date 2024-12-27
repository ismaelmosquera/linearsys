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

#ifndef ___LINEARSYS_H___
#define ___LINEARSYS_H___

#ifdef __cplusplus
extern "C" {
	#endif

#include "qr.h"
#include "lu.h"
#include "matrix.h"
#include "vector.h"

/*
	* Function to solve a linear system with a mxn matrix passed as parameter where n = m+1.
	* This function uses the Cramer's rule to solve the system.
	* param m a mxn system matrix.
	* return a vector with the solution of the system.
	* The matrix passed as parameter is the nxn coeficient matrix and the last column has the coeficient vector.
	* So it is a mxn matrix where n = m+1.
	*
	* If the system cannot be solved, this function returns NULL.
	*/
Vector* smcramer_system_solver(const Matrix* m);

/*
* Function to solve a linear system passing a nxn coeficient matrix and a n coeficient vector as parameter.
* This function uses the Cramer's rule to solve the system.
* param m a nxn coeficient matrix.
* param v a n coeficient vector.
* return a vector with the solution of the system.
* If the system cannot be solved, this function returns NULL.
*/
Vector* mvcramer_system_solver(const Matrix* m, const Vector* v);

/*
	* Function to solve a linear system with a mxn matrix passed as parameter where n = m+1.
	* This function uses Gaussian elimination to solve the system.
	* param m a mxn system matrix.
	* return a vector with the solution of the system.
	* The matrix passed as parameter is the nxn coeficient matrix and the last column has the coeficient vector.
	* So it is a mxn matrix where n = m+1.
	*
	* If the system cannot be solved, this function returns NULL.
	*/
Vector* smgauss_system_solver(const Matrix* m);

/*
* Function to solve a linear system passing a nxn coeficient matrix and a n coeficient vector as parameter.
* This function uses Gaussian elimination to solve the system.
* param m a nxn coeficient matrix.
* param v a n coeficient vector.
* return a vector with the solution of the system.
* If the system cannot be solved, this function returns NULL.
*/
Vector* mvgauss_system_solver(const Matrix* m, const Vector* v);

/*
* Function to solve a linear system using LU decomposition.
* param lu, a LU decomposition type.
* param v vector of n coeficients.
* return a vector with the solution of the system.
* If the system cannot be solved, this function returns null.
*/
Vector* lu_system_solver(const LU* lu, const Vector* v);

/*
* Function to solve a linear system using QR factorization.
* param qr, a qr factorization type.
* param v vector of n coeficients.
* return a vector with the solution of the system.
*/
Vector* qr_system_solver(const QR* qr, const Vector* v);

#ifdef __cplusplus
}
#endif

#endif
