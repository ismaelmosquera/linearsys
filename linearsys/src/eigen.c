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

#include <stdio.h>
#include <stdlib.h>
#include "qr.h"
#include "eigen.h"


#define THRESHOLD 1E-3
#define MAX_ITERATIONS 50000

/*
* helper functions.
*/
static double abs_value(double x)
{
return (x < 0.0) ? -x : x;
}

static int done(const Matrix* m)
{
	int i;
int ret = 1;
for(i = 0; i < rows_matrix(m); i++)
{
	if(abs_value(get_matrix(m, i, i)) != 1.0) return 0;
}
return ret;
}

static Matrix* upper_triangular(const Matrix* m)
{
int n, row;
int i, j, k;
double pivot, tmp, remove, t;
Matrix* upper = NULL;
if(rows_matrix(m) != columns_matrix(m)) return NULL; /* m must be square */
upper = clone_matrix(m);
n = rows_matrix(upper);
i = 0;
while(i < n)
{
pivot = abs_value(upper->_data[i*n+i]);
row = i;
for(k = i+1; k < n; k++)
{
if(abs_value(upper->_data[k*n]) > pivot)
{
	pivot = abs_value(upper->_data[k*n]);
	row = k;
}
}
if(i != row)
{
for(j = 0; j < n; j++)
{
	tmp = upper->_data[i*n+j];
	upper->_data[i*n+j] = upper->_data[row*n+j];
	upper->_data[row*n+j] = tmp;
}
}
for(j = i+1; j < n; j++)
{
remove = -(upper->_data[j*n+i]) / upper->_data[i*n+i];
for(k = 0; k < n; k++)
{
t = upper->_data[j*n+k]+remove*upper->_data[i*n+k];
upper->_data[j*n+k] = t;
}
}
i++;
}
set_matrix(upper, 1.0, n-1, n-1);
return upper;
}

/*
* Function to solve upper system.
*/
static Vector* upper_system_solver(const Matrix* m, const Vector* v)
{
int i,j;
Vector* x = create_vector(v->_size);
int n = rows_matrix(m);
x->_data[n-1] = v->_data[n-1] / m->_data[(n-1)*m->_columns+(n-1)];
for(i = n-2; i >= 0; i--)
{
x->_data[i] = v->_data[i];
for(j = i+1; j < n; j++)
{
x->_data[i] -= m->_data[i*m->_columns+j] * x->_data[j];
}
x->_data[i] /= m->_data[i*m->_columns+i];
}
return x;
}

/* end helper functions */

/* implementation */

Eigen* create_eigen(double value, const Vector* vector)
{
Eigen* eigen = (Eigen*)malloc(sizeof(Eigen));
eigen_value(eigen) = value;
eigen_vector(eigen) = clone_vector(vector);
return eigen;
}

void destroy_eigen(Eigen* eigen)
{
if(eigen == NULL) return;
eigen_value(eigen) = 0.0;
if(eigen_vector(eigen) != NULL) free(eigen_vector(eigen));
if(eigen != NULL) free(eigen);
}

Eigen* clone_eigen(const Eigen* eigen)
{
if(eigen == NULL) return NULL;
return create_eigen(eigen_value(eigen), eigen_vector(eigen));
}

void print_eigen(const Eigen* eigen)
{
if(eigen == NULL)
{
printf("\n[]\n");
return;
}
printf("eigenvalue = %.2lf\n", eigen_value(eigen));
printf("eigenvector:\n");
print_vector(eigen_vector(eigen));
}

EigenSystem* create_eigensystem(int size)
{
	int i;
EigenSystem* eigsys = (EigenSystem*)malloc(sizeof(EigenSystem));
eigsys->_size = size;
eigsys->_eigen = (Eigen**)malloc(size*sizeof(Eigen*));
for(i = 0; i < size; i++) eigsys->_eigen[i] = NULL;
return eigsys;
}

void destroy_eigensystem(EigenSystem* eigsys)
{
if(eigsys == NULL) return;
if(eigsys->_eigen != NULL) free(eigsys->_eigen);
free(eigsys);
eigsys = NULL;
}

EigenSystem* clone_eigensystem(const EigenSystem* eigsys)
{
	int i;
	EigenSystem* eigensystem = NULL;
if(eigsys == NULL) return NULL;
eigensystem = create_eigensystem(size_eigensystem(eigsys));
for(i = 0; i < size_eigensystem(eigsys); i++) eigensystem->_eigen[i] = clone_eigen(eigsys->_eigen[i]);
return eigensystem;
}

void print_eigensystem(const EigenSystem* eigsys)
{
	int i;
if(eigsys == NULL)
{
printf("\n[]\n");
return;
}
for(i = 0; i < size_eigensystem(eigsys); i++)
{
	print_eigen(eigsys->_eigen[i]);
	printf("\n");
}
}

Eigen* max_eigen_power_method(const Matrix* m)
{
	Eigen* eigen = NULL;
	Matrix* a = NULL;
	Matrix* b = NULL;
Matrix* aux = NULL;
Vector* eigenvector = NULL;
double value = 0.0;
double ant, lim;
int i, j, n, major;
if(rows_matrix(m) != columns_matrix(m)) return NULL; /* m must be square */
a = clone_matrix(m);
n = columns_matrix(a);
aux = create_matrix(n, 1);
for(i = 0; i < n; i++) aux->_data[i] = 1.0;
i = 1;
do
{
ant = value;
b = mul_matrix(a, aux);
major = 0.0;
for(j = 0; j < n; j++)
{
if(abs_value(b->_data[major]) < abs_value(b->_data[j])) major = j;
}
value = b->_data[major];
for(j = 0; j < n; j++)
{
	b->_data[j] /= value;
}
lim = abs_value(value - ant);
if(lim < THRESHOLD) break; /* done */
for(j = 0; j < n; j++)
{
aux->_data[j] = b->_data[j];
}
i++;
}while(i < MAX_ITERATIONS);
if(i >= MAX_ITERATIONS)
{
destroy_matrix(a);
destroy_matrix(b);
destroy_matrix(aux);
return NULL;
}
eigenvector = create_vector(n);
for(i = 0; i < n; i++) set_vector(eigenvector, get_matrix(b, i, 0), i);
eigen = create_eigen(value, eigenvector);

destroy_matrix(a);
destroy_matrix(b);
destroy_matrix(aux);
destroy_vector(eigenvector);

return eigen;
}

EigenSystem* eigen_system(const Matrix* m)
{
	double param, trace;
int i, j, k, n;
Matrix* upper = NULL;
Matrix* aux = NULL;
Matrix* lambda = NULL;
Vector* v = NULL;
EigenSystem* eigensys = NULL;
QR* qr = NULL;
if(rows_matrix(m) != columns_matrix(m)) return NULL; /* m must be square */
n = rows_matrix(m);
eigensys = create_eigensystem(n);
qr = qr_factorization(m); /* get QR factorization of m */
/* Compute eigenvalues */
k = 0;
while(1)
{
	qr = qr_factorization(mul_matrix(qr_r(qr), qr_q(qr)));
		if(done(qr_q(qr)) || k > MAX_ITERATIONS) break;
		k++;
}
/* eigenvalues are listed in the diagonal of lambda */
lambda = mul_matrix(qr_q(qr), qr_r(qr));
aux = create_matrix(n, n);
	v = create_vector(n);
/* compute eigenvectors and build eigensystem */
for(i = 0; i < n; i++)
{
	for(j = 0; j < n; j++)
	{
	set_matrix(aux, get_matrix(lambda, i, i), j, j);
}
upper = upper_triangular(sub_matrix(m, aux));
trace = trace_matrix(upper)-1.0;
param = (trace == 0.0) ? 1.0 : (abs_value(trace)>=1.0) ? 1.0/trace : trace;
set_vector(v, param, n-1);
eigensys->_eigen[i] = create_eigen(get_matrix(lambda, i, i),upper_system_solver(upper, v));
}
/* release previously allocated memory */
destroy_matrix(upper);
destroy_matrix(lambda);
destroy_matrix(aux);
destroy_vector(v);
destroy_qr(qr);

return eigensys;
}


/* END */
