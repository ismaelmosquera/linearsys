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

#include <stdio.h>
#include "linearsys.h"


int main()
{
	LU* lu = NULL;
	Vector* v = NULL;
	Vector* x = NULL;
Matrix* m = NULL;
Matrix* inv_m = NULL;
Matrix* r1 = NULL;
Matrix*r2 = NULL;

printf("Linear System Solver using Cramer\'s rule.\n");
printf("loading system\'s matrix...\n");
m = load_matrix("system.dat");
printf("Extended matrix of the system:\n");
print_matrix(m);
printf("\n");
x = smcramer_system_solver(m);
printf("Solution:\n");
print_vector(x);
printf("\n");
destroy_matrix(m);
destroy_vector(x);

printf("loading coefficient\'s matrix...\n");
m = load_matrix("mcoef.dat");
printf("Matrix of coefficients:\n");
print_matrix(m);
printf("\n");

printf("loading coefficient\'s vector...\n");
v = load_vector("vcoef.dat");
printf("Vector of coefficients:\n");
print_vector(v);
printf("\n");
x = mvcramer_system_solver(m, v);
printf("solution:\n");
print_vector(x);

destroy_matrix(m);
destroy_vector(v);
destroy_vector(x);

printf("\nLinear System Solver using Gaussian elimination.\n");
printf("loading system\'s matrix...\n");
m = load_matrix("system.dat");
printf("Extended matrix of the system:\n");
print_matrix(m);
printf("\n");
x = smgauss_system_solver(m);
printf("Solution:\n");
print_vector(x);
printf("\n");

destroy_matrix(m);
destroy_vector(x);

printf("loading coefficient\'s matrix...\n");
m = load_matrix("mcoef.dat");
printf("Matrix of coefficients:\n");
print_matrix(m);
printf("\n");

printf("loading coefficient\'s vector...\n");
v = load_vector("vcoef.dat");
printf("Vector of coefficients:\n");
print_vector(v);
printf("\n");
x = mvgauss_system_solver(m, v);
printf("solution:\n");
print_vector(x);

destroy_matrix(m);
destroy_vector(v);
destroy_vector(x);

printf("\nLinear system solver using LU decomposition.\n");
v = load_vector("vcoef.dat");
lu = lu_decomposition(load_matrix("mcoef.dat"));
printf("\n");
print_lu(lu);
printf("\ncoeficient vector:\n");
print_vector(v);
x = lu_system_solver(lu, v);
printf("\nsolution:\n");
print_vector(x);

destroy_lu(lu);
destroy_vector(v);
destroy_vector(x);

	printf("\nbye.\n");

return 0;
}

/* END */
