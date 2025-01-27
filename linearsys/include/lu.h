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

#ifndef ___LU_H___
#define ___LU_H___

#ifdef __cplusplus
extern "C" {
	#endif

#include "matrix.h"

/*
* LU type definition.
*/
typedef struct
{
int* _permutation;
Matrix* _lower;
Matrix* _upper;
}LU;

/*
* Releases the memory previously allocated for the LU type passed as parameter
* param: LU* lu => a pointer to a LU type
*/
void destroy_lu(LU* lu);

/*
* Performs a LU decomposition of the matrix passed as parameter
* param: const Matrix* m => a pointer to a square Matrix.
* returns:
* A pointer to a LU of the Matrix passed as parameter
*/
LU* lu_decomposition(const Matrix* m);

/*
* Prints a LU to the console.
* param: const LU* lu => a pointer to a LU
*
* The output format is as follows:
*
* lower:
*[ ... ]
* [ ... ]
* ...
* [ ... ]
*
* upper:
*[ ... ]
* [ ... ]
* ...
* [ ... ]
*
* permutation:
* [ ... ]
*
*/
void print_lu(const LU* lu);

/*
* Macro to get the permutation vector.
*/
#define lu_permutation(lu) ((lu)->_permutation)

/*
* Macros to get lower and upper matrices.
*/
#define lu_lower(lu) ((lu)->_lower)
#define lu_upper(lu) ((lu)->_upper)

#ifdef __cplusplus
}
#endif

#endif
