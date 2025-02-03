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
#include "vector.h"

/* Declare a Not a Number constant */
static const double __NaN__ = 1E-9;

/*
* Helper function to compute the square of a value.
*/
static double __square_(double x)
{
return x*x;
}

/* Implementation */
Vector* create_vector(int size)
{
	int i;
Vector* v = (Vector*)malloc(sizeof(Vector));
v->_size = size;
v->_data = (double*)malloc(v->_size * sizeof(double));
for(i = 0; i < v->_size; i++) v->_data[i] = 0.0;
return v;
}

void destroy_vector(Vector* v)
{
if(!v) return;
if(v->_data) free(v->_data);
v->_data = NULL;
free(v);
v = NULL;
}

Vector* clone_vector(const Vector* v)
{
int i;
Vector* result = create_vector(size_vector(v));
for(i = 0; i < result->_size; i++)
{
	result->_data[i] = v->_data[i];
}
return result;
}

Vector* resize_vector(const Vector* v, int new_size)
{
	int i;
Vector* t = NULL;
if(size_vector(v) <1 || new_size < 1) return NULL;
if(new_size == size_vector(v)) return clone_vector(v);
if(new_size < size_vector(v)) return get_vector_chunk(v, 0, new_size-1);
t = create_vector(new_size);
for(i = 0; i < size_vector(v); i++) t->_data[i] = v->_data[i];
return t;
}

Vector* get_vector_chunk(const Vector* v, int from_index, int to_index)
{
int i;
int k = 0;
int new_size;
Vector* t = NULL;
if(from_index < 0 || to_index > size_vector(v)-1 || from_index > to_index) return NULL;
new_size = (to_index-from_index) + 1;
t = create_vector(new_size);
for(i = from_index; i <= to_index; i++) t->_data[k++] = v->_data[i];
return t;
}

Vector* load_vector(const char* filename)
{
int i, size;
double d;
Vector* v = NULL;
FILE* file = fopen_read(filename);
fread_int(file, &size);
v = create_vector(size);
for(i = 0; i < v->_size; i++)
{
fread_double(file, &d);
v->_data[i] = d;
}
fclose(file);
return v;
}

void store_vector(const Vector* v, const char* filename)
{
int i;
FILE* file = fopen_write(filename);
fwrite_int(file, v->_size);
fprintf(file, "\n");
for(i = 0; i < v->_size; i++)
{
if(i > 0) fprintf(file, " ");
fwrite_double(file, v->_data[i]);
}
fprintf(file, "\n");
fclose(file);
}

void print_vector(const Vector* v)
{
int i;
printf("[");
for(i = 0; i < v->_size; i++)
{
	if(i > 0) printf(", ");
	printf("%.2lf", v->_data[i]);
}
printf("]\n");
}

double get_vector(const Vector* v, int pos)
{
if(pos < 0 || pos > size_vector(v)-1) return __NaN__;
return v->_data[pos];
}

void set_vector(Vector* v, double value, int pos)
{
if(pos < 0 || pos > size_vector(v)-1) return;
v->_data[pos] = value;
}

Vector* add_vector(const Vector* v1, const Vector* v2)
{
int i;
Vector* result = NULL;
if(size_vector(v1) != size_vector(v2)) return NULL;
result = create_vector(size_vector(v1));
for(i = 0; i < result->_size; i++)
{
result->_data[i] = v1->_data[i] + v2->_data[i];
}
return result;
}

Vector* sub_vector(const Vector* v1, const Vector* v2)
{
int i;
Vector* result = NULL;
if(size_vector(v1) != size_vector(v2)) return NULL;
result = create_vector(size_vector(v1));
for(i = 0; i < result->_size; i++)
{
result->_data[i] = v1->_data[i] - v2->_data[i];
}
return result;
}

Vector* scale_vector(const Vector* v, double value)
{
int i;
Vector* result = create_vector(size_vector(v));
for(i = 0; i < result->_size; i++)
{
result->_data[i] = v->_data[i] * value;
}
return result;
}

double module_vector(const Vector* v)
{
int i;
double sum = 0.0;
for(i = 0; i < size_vector(v); i++)
{
	sum += __square_(v->_data[i]);
}
return sqrt(sum);
}

Vector* normalize_vector(const Vector* v)
{
int i;
Vector* out = clone_vector(v);
double d = module_vector(out);
if(!d) return NULL;
for(i = 0; i < size_vector(v); i++)
{
out->_data[i] /= d;
}
return out;
}

double dot_product_vector(const Vector* v1, const Vector* v2)
{
int i;
double dp = 0.0;
if(size_vector(v1) != size_vector(v2)) return __NaN__;
for(i = 0; i < size_vector(v1); i++)
{
dp += v1->_data[i] * v2->_data[i];
}
return dp;
}

Vector* cross_product_vector(const Vector* v1, const Vector* v2)
{
Vector* result = NULL;
if(size_vector(v1) != 3 || size_vector(v1) != size_vector(v2)) return NULL;
result = create_vector(size_vector(v1));
result->_data[0] = v1->_data[1]*v2->_data[2] - v1->_data[2]*v2->_data[1];
result->_data[1] = v1->_data[2]*v2->_data[0] - v1->_data[0]*v2->_data[2];
result->_data[2] = v1->_data[0]*v2->_data[1] - v1->_data[1]*v2->_data[0];
return result;
}

/* END */
