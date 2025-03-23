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
#include <stdlib.h>
#include <math.h>
#include "numio.h"
#include "matrix.h"
#include "linearsys.h"
#include "svd.h"

#define THRESHOLD 1E-6

/* Declare a Not a Number constant */
static const double __NaN__ = 1E-9;


/* Helper functions */

/*
* Static function to compute the determinant of a square Matrix.
*/
static double __det_(Matrix* m, int n)
{
int i, j, k;
int ai, aj;
double sign = 1.0;
double d = 0.0;
Matrix* adj = NULL;
if(n == 1)
{
return *(m->_data);
}
else
{
	adj = create_matrix(n-1, n-1);
	for(k = 0; k < n; k++)
	{
ai=0;
for(i=0;i<n;i++)
{
aj=0;
for(j=1;j<n;j++)
if(i != k)
{
	*(adj->_data + ai*adj->_columns + aj) = *(m->_data + i*m->_columns + j);
aj++;
}
if(i != k) ai++;
}
		d+=sign * *(m->_data + k*m->_columns) * __det_(adj, n-1);
		sign *= -1.0;
	}
}
return d;
}

/*
* Helper function to get a minimum value from the parameters.
*/
static double __get_min(double a, double b)
{
return (a <= b) ? a : b;
}

/*
* Helper function to set values in src to dest.
*/
static void __fill_matrix(const Matrix* src, Matrix* dest, int rows, int columns)
{
int i, j;
for(i = 0; i < rows; i++)
{
	for(j = 0; j < columns; j++)
	{
		set_matrix(dest, get_matrix(src, i, j), i, j);
	}
}
}

/*
* Helper function to get the square of a value.
*/
static double _square_(double x)
{
return x*x;
}

/* end helper functions */

/* implementation */
Matrix* create_matrix(int r, int c)
{
int i, j, dim;
Matrix* m = (Matrix*)malloc(sizeof(Matrix));
m->_rows = r;
m->_columns = c;
dim = m->_rows * m->_columns;
m->_data = (double*)malloc(dim*sizeof(double));
for(i = 0; i < m->_rows; i++)
{
for(j = 0; j < m->_columns; j++)
{
	*(m->_data + i*m->_columns + j) = 0.0;
}
}
return m;
}

void destroy_matrix(Matrix* m)
{
if(!m) return;
if(m->_data) free(m->_data);
m->_data = NULL;
free(m);
m = NULL;
}

Matrix* clone_matrix(const Matrix* m)
{
	int i, j;
	Matrix* result = NULL;
if(!m || m->_rows < 1 || m->_columns <1) return NULL;
result = create_matrix(m->_rows, m->_columns);
for(i = 0; i < result->_rows; i++)
{
for(j = 0; j < result->_columns; j++)
{
	*(result->_data + i*result->_columns + j) = *(m->_data + i*m->_columns + j);
}
}
return result;
}

Matrix* load_matrix(const char* filename)
{
	double d;
	int i, j, r, c;
	Matrix* result = NULL;
FILE* file = fopen_read(filename);
fread_int(file, &r);
fread_int(file, &c);
if(r < 1 || c < 1) return NULL;
result = create_matrix(r, c);
for(i = 0; i < result->_rows; i++)
{
for(j = 0; j < result->_columns; j++)
{
	fread_double(file, &d);
	*(result->_data + i*result->_columns + j) = d;
}
}
fclose(file);
return result;
}

void store_matrix(const Matrix* m, const char* filename)
{
int i, j;
FILE* file = fopen_write(filename);
fwrite_int(file, m->_rows);
fprintf(file, " ");
fwrite_int(file, m->_columns);
fprintf(file, "\n");
for(i = 0; i < m->_rows; i++)
{
	for(j = 0; j < m->_columns; j++)
	{
		if(j > 0) fprintf(file, " ");
		fwrite_double(file, *(m->_data + i*m->_columns +j));
	}
	fprintf(file, "\n");
}
fclose(file);
}

void print_matrix(const Matrix* m)
{
int i, j;
for(i = 0; i < m->_rows; i++)
{
	printf("[");
	for(j = 0; j < m->_columns; j++)
	{
		if(j > 0) printf(", ");
		printf("%.2lf", *(m->_data + i*m->_columns +j));
	}
	printf("]\n");
}
}

double get_matrix(const Matrix* m, int i, int j)
{
if(i<0 || i>m->_rows-1 || j<0 || j>m->_columns-1) return __NaN__;
return *(m->_data + i*m->_columns + j);
}

void set_matrix(Matrix* m, double value, int i, int j)
{
if(i<0 || i>m->_rows-1 || j<0 || j>m->_columns-1) return;
*(m->_data + i*m->_columns + j) = value;
}

Matrix* get_matrix_chunk(const Matrix* m, int from_row, int to_row, int from_column, int to_column)
{
	int i, j;
	int s, t;
	int new_rows;
	int new_columns;
	Matrix* r = NULL;
	if(from_row < 0 || to_row > rows_matrix(m)-1 || from_row > to_row ||
	from_column < 0 || to_column > columns_matrix(m)-1 || from_column > to_column) return NULL;
new_rows = (to_row - from_row) + 1;
new_columns = (to_column - from_column) + 1;
r = create_matrix(new_rows, new_columns);
i = 0;
for(s = from_row; s <= to_row; s++)
{
	j = 0;
	for(t = from_column; t <= to_column; t++)
	{
	set_matrix(r, get_matrix(m, s, t), i, j);
	j++;
	}
	i++;
}
return r;
}

Matrix* resize_matrix(const Matrix* m, int nrows, int ncolumns)
{
	int rows, columns;
Matrix* r = NULL;
if(nrows < 1 || ncolumns < 1) return NULL;
r = create_matrix(nrows, ncolumns);
rows = __get_min(rows_matrix(m), rows_matrix(r));
columns = __get_min(columns_matrix(m), columns_matrix(r));
__fill_matrix(m, r, rows, columns);
return r;
}

Matrix* resize_matrix_rows(const Matrix* m, int nrows)
{
int rows;
Matrix* r = NULL;
if(nrows < 1) return NULL;
r = create_matrix(nrows, columns_matrix(m));
rows = __get_min(rows_matrix(m), rows_matrix(r));
__fill_matrix(m, r, rows, columns_matrix(r));
return r;
}

Matrix* resize_matrix_columns(const Matrix* m, int ncolumns)
{
int columns;
Matrix* r = NULL;
if(ncolumns < 1) return NULL;
r = create_matrix(rows_matrix(m), ncolumns);
columns = __get_min(columns_matrix(m), columns_matrix(r));
__fill_matrix(m, r, rows_matrix(r), columns);
return r;
}

Matrix* add_matrix(const Matrix* m1, const Matrix* m2)
{
int i, j;
Matrix* m = NULL;
if(m1->_rows!=m2->_rows || m1->_columns!=m2->_columns) return NULL;
m = create_matrix(m1->_rows, m1->_columns);
for(i = 0; i < m->_rows; i++)
{
for(j = 0; j < m->_columns; j++)
{
	*(m->_data + i*m->_columns + j) = *(m1->_data + i*m1->_columns + j) + *(m2->_data + i*m2->_columns + j);
}
}
return m;
}

Matrix* sub_matrix(const Matrix* m1, const Matrix* m2)
{
int i, j;
Matrix* m = NULL;
if(m1->_rows!=m2->_rows || m1->_columns!=m2->_columns) return NULL;
m = create_matrix(m1->_rows, m1->_columns);
for(i = 0; i < m->_rows; i++)
{
for(j = 0; j < m->_columns; j++)
{
	*(m->_data + i*m->_columns + j) = *(m1->_data + i*m1->_columns + j) - *(m2->_data + i*m2->_columns + j);
}
}
return m;
}

Matrix* mul_matrix(const Matrix* m1, const Matrix* m2)
{
	double accum = 0.0;
int i, j, k;
Matrix* m = NULL;
if(m1->_columns != m2->_rows) return NULL;
m = create_matrix(m1->_rows, m2->_columns);
for(i = 0; i < m->_rows; i++)
{
for(j = 0; j < m->_columns; j++)
{
	accum = 0.0;
	for(k = 0; k < m1->_columns; k++)
	{
		accum += *(m1->_data + i*m1->_columns + k) * *(m2->_data + k*m2->_columns + j);
	}
	*(m->_data + i*m->_columns + j) = accum;
}
}
return m;
}

Matrix* scale_matrix(const Matrix* m, double value)
{
	int i, j;
Matrix* out = NULL;
if(m == NULL) return NULL;
out = clone_matrix(m);
for(i = 0; i < rows_matrix(out); i++)
{
	for(j = 0; j < columns_matrix(out); j++)
	{
		out->_data[i*columns_matrix(out)+j] *= value;
	}
}
return out;
}

Matrix* identity_matrix(int n)
{
int i;
Matrix* m = create_matrix(n, n);
for(i = 0; i < m->_rows; i++)
{
	*(m->_data + i*m->_columns +i) = 1.0;
}
return m;
}

Matrix* transpose_matrix(const Matrix* m)
{
int i, j;
Matrix* t = create_matrix(m->_columns, m->_rows);
for(i = 0; i < t->_rows; i++)
{
for(j = 0; j < t->_columns; j++)
{
	*(t->_data + i*t->_columns +j) = *(m->_data + j*m->_columns +i);
}
}
return t;
}

double trace_matrix(const Matrix* m)
{
	int i;
double t = 0.0;
if(m == NULL) return 0.0;
if(rows_matrix(m) != columns_matrix(m)) return 0.0; /* m must be square */
for(i = 0; i < columns_matrix(m); i++)
{
	t += get_matrix(m, i, i);
}
return t;
}

Matrix* diag_matrix(const Matrix* m)
{
	int i;
	Matrix* out = NULL;
if(m == NULL) return NULL;
if(rows_matrix(m) != columns_matrix(m)) return NULL; /* m must be square */
	out = create_matrix(rows_matrix(m), columns_matrix(m));
	for(i = 0; i < columns_matrix(m); i++)
	{
		set_matrix(out, get_matrix(m, i, i), i, i);
	}
	return out;
}

Matrix* inverse_matrix(const Matrix* m)
{
int i, j;
Matrix* out = NULL;
Vector* v = NULL;
Vector* x = NULL;
LU* lu = NULL;
if(rows_matrix(m) != columns_matrix(m)) return NULL;
lu = lu_decomposition(m);
v = create_vector(rows_matrix(m));
x = create_vector(rows_matrix(m));
out = create_matrix(rows_matrix(m), rows_matrix(m));

for(i = 0; i < rows_matrix(m); i++)
{
	if(i > 0) v->_data[i-1] = 0.0;
v->_data[i] = 1.0;
x = lu_system_solver(lu, v);
for(j = 0; j < out->_columns; j++)
{
	out->_data[j*out->_columns+i] = x->_data[j];
}
}

destroy_lu(lu);
destroy_vector(v);
destroy_vector(x);

return out;
}

double det_matrix(const Matrix* m)
{
Matrix* t = NULL;
if(m->_rows != m->_columns) return 0.0;
t = clone_matrix(m);
return __det_(t, t->_rows);
}

Matrix* rotation_matrix(double x,double y,double z,double w)
{
Matrix* r = create_matrix(3, 3);
double d = sqrt(x*x+y*y+z*z);
if(d == 0) d = 1.0;

double u1=x/d;
double u2=y/d;
double u3=z/d;

double pos00 = u1*u1+(cos(w)*(1.0-u1*u1));
	 double pos01 = (u1*u2*(1.0-cos(w)))-(u3*sin(w));
	 double pos02 = (u3*u1*(1.0-cos(w)))+(u2*sin(w));
	 double pos10 = (u1*u2*(1.0-cos(w)))+(u3*sin(w));
	 double pos11 = u2*u2+(cos(w)*(1.0-u2*u2));
	 double pos12 = (u2*u3*(1.0-cos(w)))-(u1*sin(w));
	 double pos20 = (u3*u1*(1.0-cos(w)))-(u2*sin(w));
	 double pos21 = (u2*u3*(1.0-cos(w)))+(u1*sin(w));
	 double pos22 = u3*u3+(cos(w)*(1.0-u3*u3));

set_matrix(r, pos00, 0, 0);
set_matrix(r, pos01, 0, 1);
set_matrix(r, pos02, 0, 2);
set_matrix(r, pos10, 1, 0);
set_matrix(r, pos11, 1, 1);
set_matrix(r, pos12, 1, 2);
set_matrix(r, pos20, 2, 0);
set_matrix(r, pos21, 2, 1);
set_matrix(r, pos22, 2, 2);

return r;
}

Matrix* get_row_matrix(const Matrix* m, int i)
{
	Matrix* out = NULL;
int j;
int size = columns_matrix(m);
if(i<0 || i > rows_matrix(m)-1) return NULL;
out = create_matrix(1, size);
for(j = 0; j < size; j++)
{
set_matrix(out, get_matrix(m, i, j), 0, j);
}
return out;
}

Matrix* get_column_matrix(const Matrix* m, int j)
{
Matrix* out = NULL;
int i;
int size = rows_matrix(m);
if(j<0 || j>columns_matrix(m)-1) return NULL;
out = create_matrix(size, 1);
for(i = 0; i < size; i++)
{
set_matrix(out, get_matrix(m, i, j), i, 0);
}
return out;
}

Vector* get_row_vector(const Matrix* m, int i)
{
	Vector* v = NULL;
int j;
int size = columns_matrix(m);
if(i<0 || i > rows_matrix(m)-1) return NULL;
v = create_vector(size);
for(j = 0; j < size; j++)
{
set_vector(v, get_matrix(m, i, j), j);
}
return v;
}

Vector* get_column_vector(const Matrix* m, int j)
{
Vector* v = NULL;
int i;
int size = rows_matrix(m);
if(j<0 || j>columns_matrix(m)-1) return NULL;
v = create_vector(size);
for(i = 0; i < size; i++)
{
set_vector(v, get_matrix(m, i, j), i);
}
return v;
}

int is_orthogonal_matrix(const Matrix* m)
{
	int i, j;
if(rows_matrix(m) != columns_matrix(m)) return 0;
for(i= 0; i < columns_matrix(m)-1; i++)
{
for(j = i+1; j < columns_matrix(m); j++)
{
if(dot_product_vector(get_column_vector(m, i), get_column_vector(m, j)) > THRESHOLD) return 0;
}
}
return 1;
}

int is_orthonormal_matrix(const Matrix* m)
{
	int j;
	if(m == NULL) return 0;
	if(!is_orthogonal_matrix(m)) return 0;
	for(j = 0; j < columns_matrix(m); j++)
	{
		if(round(module_vector(get_column_vector(m, j))) != 1) return 0;
		}
	return 1;
}

Matrix* set_row_matrix(const Matrix* m, const Vector* v, int row)
{
	int j;
	Matrix* out = NULL;
if(m == NULL) return NULL;
if(row < 0 || row > rows_matrix(m)-1) return clone_matrix(m);
if(size_vector(v) != columns_matrix(m)) return clone_matrix(m);
out = clone_matrix(m);
for(j = 0; j < columns_matrix(m); j++)
{
	set_matrix(out, get_vector(v, j), row, j);
}
return out;
}

Matrix* set_column_matrix(const Matrix* m, const Vector* v, int column)
{
int i;
	Matrix* out = NULL;
if(m == NULL) return NULL;
if(column < 0 || column > columns_matrix(m)-1) return clone_matrix(m);
if(size_vector(v) != rows_matrix(m)) return clone_matrix(m);
out = clone_matrix(m);
for(i = 0; i < rows_matrix(m); i++)
{
	set_matrix(out, get_vector(v, i), i, column);
}
return out;
}

Matrix* pseudoinverse_matrix(const Matrix* m)
{
int i, k;
SVD* svd = NULL;
Matrix* sigma_t = NULL;
Matrix* out = NULL;
if(m == NULL) return NULL;
svd = svd_factorization(m);
k = (rows_matrix(m) <= columns_matrix(m)) ? rows_matrix(m) : columns_matrix(m);
sigma_t = transpose_matrix(svd_sigma(svd));
for(i = 0; i < k; i++)
{
set_matrix(sigma_t, 1.0 / get_matrix(sigma_t, i, i), i, i);
}
out = mul_matrix(mul_matrix(svd_v(svd), sigma_t), transpose_matrix(svd_u(svd)));

/* release previously alocated memory */
destroy_svd(svd);
destroy_matrix(sigma_t);

/* return pseudoinverse matrix */
return out;
}

Matrix* nearest_orthogonal_matrix(const Matrix* m)
{
Matrix* out = NULL;
SVD* svd = NULL;
if(rows_matrix(m) != columns_matrix(m) || m == NULL) return NULL;
svd = svd_factorization(m);
out = mul_matrix(svd_u(svd), transpose_matrix(svd_v(svd)));

/* release previously allocated memory */
destroy_svd(svd);

/* return nearest orthogonal matrix of 'm' */
return out;
}

double forbinius_norm_matrix(const Matrix* m)
{
	int i, j;
double sum = 0.0;
if(m == NULL) return -1.0;
for(i = 0; i < rows_matrix(m); i++)
{
for(j = 0; j < columns_matrix(m); j++)
{
	sum += _square_(get_matrix(m, i, j));
}
}
return sqrt(sum);
}


/* END */
