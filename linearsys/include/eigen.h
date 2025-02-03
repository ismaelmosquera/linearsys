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

#ifndef ___EIGEN_H___
#define ___EIGEN_H___

#ifdef __cplusplus
extern "C" {
	#endif

#include "matrix.h"
#include "vector.h"

/*
* Eigen type definition.
* An Eigen is a pair ( eigenvalue, eigenvector )
*/
typedef struct
{
double _value;
Vector* _vector;
}Eigen;

/*
* EigenSystem type definition.
* An EigenSystem is the set of eigens of a matrix.
*/
typedef struct
{
int _size; /* number of Eigens in the EigenSystem */
Eigen** _eigen; /* array of Eigens */
}EigenSystem;

/*
* Creates an Eigen.
* param: double value eigenvalue.
* param: Vector* vector Asociated eigenvector.
*
* returns: A pointer to the newly created Eigen.
*
*/
Eigen* create_eigen(double value, const Vector* vector);

/*
* Destroys an Eigen.
* param: Eigen* eigen An Eigen to destroy.
*
*/
void destroy_eigen(Eigen* eigen);

/*
* Makes a clone of the Eigen passed as parameter.
* param: const Eigen* eigen An Eigen to clone.
*
* returns: cloned Eigen.
*
*/
Eigen* clone_eigen(const Eigen* eigen);

/*
* Prints an Eigen to the console.
* param: const Eigen* eigen An Eigen to print.
*
*/
void print_eigen(const Eigen* eigen);

/*
* Creates an EigenSystem.
* Initially, all the eigens are set to NULL.
* param: size Number of Eigens in the system.
*
* returns: An EigenSystem with all its Eigens set to NULL.
*
*/
EigenSystem* create_eigensystem(int size);

/*
* Destroys an EigenSystem.
* param: eigsys An EigenSystem to destroy.
*
*/
void destroy_eigensystem(EigenSystem* eigsys);

/*
* Clones an EigenSystem.
* param: eigsys An EigenSystem to clone.
*
* returns: A clone of the EigenSystem passed as parameter.
*
*/
EigenSystem* clone_eigensystem(const EigenSystem* eigsys);

/*
* Prints an EigenSystem to the console.
* param: eigsys An EigenSystem to print.
*
*/
void print_eigensystem(const EigenSystem* eigsys);

/*
* Computes the maximum Eigen from a square matrix using the Power Method.
* param: const Matrix* m A square matrix to compute its maximum eigen.
*
* returns: Maximum eigen or NULL if the operation cannot be done.
*
*/
Eigen* max_eigen_power_method(const Matrix* m);

/*
* computes the EigenSystem for a square matrix using QR algorithm.
* param: m a square matrix.
*
* returns: EigenSystem for the matrix passed as parameter.
*
*/
EigenSystem* eigen_system(const Matrix* m);


  /*
  * Macros to access an Eigen structure.
  */
  #define eigen_value(eigen) ((eigen)->_value)
  #define eigen_vector(eigen) ((eigen)->_vector)

/*
* Macro to get the size of an EigenSystem.
*/
#define size_eigensystem(eigsys) ((eigsys)->_size)

#ifdef __cplusplus
}
#endif

#endif
