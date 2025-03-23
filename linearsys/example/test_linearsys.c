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
#include "eigen.h"
#include "svd.h"

#define PI 3.14159265359

int main()
{
	SVD* svd = NULL;
EigenSystem* eigsys = NULL;
Eigen* maxeigen = NULL;
	QR* qr = NULL;
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

printf("\nInverse matrix computation.\n\n");
m = load_matrix("mcoef.dat");
printf("M:\n");
print_matrix(m);
inv_m = inverse_matrix(m);
printf("\nM^-1:\n");
print_matrix(inv_m);
printf("\nM * M^-1 = I:\n");
print_matrix(mul_matrix(m, inv_m));
printf("\n");

destroy_matrix(m);
destroy_matrix(inv_m);

printf("\nBuilding a rotation matrix r1 ...\n");
r1 = rotation_matrix(1.0, 1.0, 0.0, PI / 4.0);
printf("Rotation matrix r1:\n");
print_matrix(r1);
printf("\n");
printf("Compute the transpose matrix of r1...\n");
r2 = transpose_matrix(r1);
printf("Print the transpose matrix r2:\n");
print_matrix(r2);
printf("\n");
printf("The transposed of a rotation matrix is its inverse.\n");
printf("So, r1*r2=I\n");
print_matrix(mul_matrix(r1, r2));
printf("\n");

destroy_matrix(r1);
destroy_matrix(r2);

printf("Compute QR factorization:\n");
m = load_matrix("mcoef.dat");
printf("original matrix (A):\n");
print_matrix(m);
qr = qr_factorization(m);
print_qr(qr);
printf("\n");
printf("A = QR:\n");
print_matrix(mul_matrix(qr_q(qr), qr_r(qr)));
printf("\n");
printf("q is orthogonal = %d\n", is_orthogonal_matrix(qr_q(qr)));
printf("Q is orthonormal = %d\n", is_orthonormal_matrix(qr_q(qr)));
printf("det(A) = %.2lf\n", qr_det(qr));
printf("\n");

printf("Solve linear system using QR factorization:\n");
printf("Coefficient's vector:\n");
print_vector(load_vector("vcoef.dat"));
printf("\n");
printf("System solution:\n");
print_vector(qr_system_solver(qr, load_vector("vcoef.dat")));
printf("\n");

printf("Eigen:\n");
printf("Matrix:\n");
print_matrix(m);
printf("\n");
printf("Find the max Eigen using the Power Method:\n");
maxeigen = max_eigen_power_method(m);
print_eigen(maxeigen);
printf("\n");

/* compute eigensystem */
printf("Find the eigensystem using the QR Algorithm:\n");
eigsys = eigen_system(m);
printf("EigenSystem:\n");
print_eigensystem(eigsys);

printf("singular Value Decomposition\n");
m = load_matrix("m.dat");
printf("M:\n");
print_matrix(m);
printf("\n");

printf("SVD:\n");
svd = svd_factorization(m);
print_svd(svd);

printf("U is orthogonal = %d\n", is_orthogonal_matrix(svd_u(svd)));
printf("V is orthogonal = %d\n", is_orthogonal_matrix(svd_v(svd)));
printf("U is orthonormal = %d\n", is_orthonormal_matrix(svd_u(svd)));
printf("V is orthonormal = %d\n", is_orthonormal_matrix(svd_v(svd)));
printf("\n");

printf("det(U) = %.2lf\n", det_matrix(svd_u(svd)));
printf("det(V) = %.2lf\n", det_matrix(svd_v(svd)));
printf("\n");

printf("M = UDV^t:\n");
print_matrix(mul_matrix(mul_matrix(svd_u(svd), svd_sigma(svd)), transpose_matrix(svd_v(svd))));
printf("\n");

destroy_eigen(maxeigen);
destroy_eigensystem(eigsys);
destroy_matrix(m);
destroy_qr(qr);
destroy_svd(svd);

	printf("bye.\n");

return 0;
}

/* END */
