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

#include <stddef.h>
#include "linearsys.h"

/*
* Declaration of determination system threshold.
*/
static const double __threshold__ = 1E-6;

/*
* Helper function to compute the absolute value.
*/
static double __abs_(double x)
{
return (x < 0) ? -x : x;
}

/*
* Helper function to use for solve linear systems using Cramer's rule.
*/
static double __k_det_(const Matrix* m, const Vector* v, int k)
{
	int i;
Matrix* a = clone_matrix(m);
for(i = 0; i < a->_rows; i++)
{
	*(a->_data + i*a->_columns + k) = v->_data[i];
}
return det_matrix(a);
}

/*
* Helper function to perform Gaussian elimination.
*/
static Matrix* __gaussian_elimination_(const Matrix* m)
{
int row, i, j, k;
double pivot, elim, tmp;
Matrix* u = clone_matrix(m);
i = 0;
 while(i < u->_rows)
 {
	 /* find pivot */
  pivot = __abs_(*(u->_data + i*u->_columns + i));
  row = i;
  for(k = i+1; k < u->_rows; k++)
   {
   if(__abs_(*(u->_data + k*m->_columns)) > pivot)
   {
    pivot = __abs_(*(u->_data + k*m->_columns));
    row = k;
   }
}
   /* check for determination */
   if(pivot < __threshold__) return NULL;
   if(i != row) /* swap rows if necessary */
  {
    for(j = 0; j < u->_columns; j++)
    {
     tmp = *(u->_data + i*u->_columns + j);
     *(u->_data + i*u->_columns + j) = *(u->_data + row*u->_columns + j);
     *(u->_data + row*u->_columns + j) = tmp;
}
}
  for(j = i+1; j < u->_rows; j++)
  {
	  /* find zeros */
   elim = -(*(u->_data + j*u->_columns + i)) / *(u->_data + i*u->_columns + i);
   for(k = 0; k < u->_columns; k++)
   {
    *(u->_data + j*u->_columns + k) += elim * *(u->_data + i*u->_columns + k);
}
   }
 i++;
}
return u;
}

/*
* Helper function to solve Gaussian triangular system.
*/
static Vector* __triangular_system_solver_(const Matrix* m)
{
int i, j;
int dim = dimension_matrix(m);
Vector* x = create_vector(m->_rows);
/*
* compute xn.
*/
x->_data[m->_rows-1] = *(m->_data + dim-1) / *(m->_data + dim-2);
/*
* Compute xn-1 xn-2 xn-3 ... x0
*/
for(i = m->_rows-2; i >= 0; i--)
{
	x->_data[i] = *(m->_data + i*m->_columns + m->_rows);
	for(j = i+1; j < m->_rows; j++) x->_data[i] -= *(m->_data + i*m->_columns + j) * x->_data[j];
x->_data[i] /= *(m->_data + i*m->_columns + i);
}
return x;
}

/* Helper functions to solve LU systems */
/*
* Function to handle permutation.
*/
static Vector* handle_permutation(int* p, Vector* v)
{
	int i;
	Vector* x = clone_vector(v);
	for(i = 0; i < x->_size; i++)
	{
		x->_data[i] = v->_data[p[i]];
	}
	return x;
}

/*
* Function to solve lower system.
*/
static Vector* lower_system_solver(const Matrix* m, const Vector* v)
{
int i, j;
Vector* x = create_vector(v->_size);
int n = rows_matrix(m);
x->_data[0] = v->_data[0] / m->_data[0];
for(i = 1; i < n; i++)
{
x->_data[i] = v->_data[i];
for(j = 0; j < i; j++)
{
x->_data[i] -= m->_data[i*m->_columns+j] * x->_data[j];
}
x->_data[i] /= m->_data[i*m->_columns+i];
}
return x;
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

/* End helper functions */

/* implementation */

Vector* smcramer_system_solver(const Matrix* m)
{
	double d;
int i, j;
Matrix* a = NULL;
Vector* x = NULL;
Vector* s = NULL;
if(m->_rows+1 != m->_columns) return NULL;
a = create_matrix(m->_rows, m->_rows);
for(i = 0; i < a->_rows; i++)
{
for(j = 0; j < a->_columns; j++)
{
	*(a->_data + i*a->_columns + j) = *(m->_data + i*m->_columns + j);
}
}
d = det_matrix(a);
if(!d) return NULL;
x = create_vector(m->_rows);
for(i = 0; i < x->_size; i++)
{
x->_data[i] = *(m->_data + i*m->_columns + m->_columns-1);
}

s = create_vector(x->_size);
for(i = 0; i < s->_size; i++)
{
s->_data[i] = __k_det_(a, x, i) / d;
}
destroy_matrix(a);
destroy_vector(x);
return s;
}

Vector* mvcramer_system_solver(const Matrix* m, const Vector* v)
{
	int i;
	double d;
Vector* x = NULL;
if(m->_rows != m->_columns || v->_size != m->_rows) return NULL;
d = det_matrix(m);
if(!d) return NULL;
x = create_vector(v->_size);
for(i = 0; i < x->_size; i++)
{
x->_data[i] = __k_det_(m, v, i) / d;
}
return x;
}

Vector* smgauss_system_solver(const Matrix* m)
{
if(m->_rows+1 != m->_columns) return NULL;
return __triangular_system_solver_(__gaussian_elimination_(m));
}

Vector* mvgauss_system_solver(const Matrix* m, const Vector* v)
{
	int i, j;
Matrix* s = NULL;
if(m->_rows != m->_columns || v->_size != m->_rows) return NULL;
s = create_matrix(m->_rows, m->_rows+1);
/*
* Build extended system's matrix.
*/
for(i = 0; i < m->_rows; i++)
{
for(j = 0; j < m->_columns; j++)
{
	*(s->_data + i*s->_columns + j) = *(m->_data + i*m->_columns + j);
}
}
j = s->_columns - 1;
for(i = 0; i < s->_rows; i++)
{
	*(s->_data + i*s->_columns + j) = v->_data[i];
}
return __triangular_system_solver_(__gaussian_elimination_(s));
}

Vector* lu_system_solver(const LU* lu, const Vector* v)
{
int i, pivoting = 0;
Vector* x = NULL;
if(rows_matrix(lu->_lower) != v->_size) return NULL;
x = clone_vector(v);
for(i = 0; i < x->_size; i++)
{
if(x->_data[i] != lu->_permutation[i])
{
	pivoting = 1;
	break;
}
}
if(pivoting) x = handle_permutation(lu->_permutation, x);
return upper_system_solver(lu->_upper, lower_system_solver(lu->_lower, x));
}


/* END */
