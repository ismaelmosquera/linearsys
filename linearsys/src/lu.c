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
#include "lu.h"

/* declare a threshold constant to check for determination */
static const double __threshold_ = 1E-6;

/*
* static function to get the absolute value.
*/
static double __abs_(double x)
{
return (x < 0.0) ? -x : x;
}

/* Implementation */

void destroy_lu(LU* lu)
{
if(lu == NULL) return;
if(lu_permutation(lu) != NULL)
{
free(lu->_permutation);
lu->_permutation = NULL;
}
if(lu->_lower != NULL) destroy_matrix(lu_lower(lu));
if(lu->_upper != NULL) destroy_matrix(lu_upper(lu));
free(lu);
lu = NULL;
}

LU* lu_decomposition(const Matrix* m)
{
int n, row;
int i, j, k;
double pivot, tmp, remove, t;
LU* lu = NULL;
if(rows_matrix(m) != columns_matrix(m)) return NULL;
lu = (LU*)malloc(sizeof(LU));
lu->_upper = clone_matrix(m);
n = rows_matrix(lu->_upper);
lu->_lower = identity_matrix(n);
lu->_permutation = (int*)malloc(n*sizeof(int));
for(i = 0; i < n; i++) lu->_permutation[i] = i;
i = 0;
while(i < n)
{
pivot = __abs_(lu->_upper->_data[i*n+i]);
row = i;
for(k = i+1; k < n; k++)
{
if(__abs_(lu->_upper->_data[k*n]) > pivot)
{
	pivot = __abs_(lu->_upper->_data[k*n]);
	row = k;
}
}
if(pivot < __threshold_) return NULL;
if(i != row)
{
for(j = 0; j < n; j++)
{
	tmp = lu->_upper->_data[i*n+j];
	lu->_upper->_data[i*n+j] = lu->_upper->_data[row*n+j];
	lu->_upper->_data[row*n+j] = tmp;
	lu->_permutation[i] = row;
	lu->_permutation[row] = i;
}
}
for(j = i+1; j < n; j++)
{
remove = -(lu->_upper->_data[j*n+i]) / lu->_upper->_data[i*n+i];
lu->_lower->_data[j*n+i] = -remove;
for(k = 0; k < n; k++)
{
t = lu->_upper->_data[j*n+k]+remove*lu->_upper->_data[i*n+k];
lu->_upper->_data[j*n+k] = t;
}
}
i++;
}
return lu;
}

void print_lu(const LU* lu)
{
	int i;
if(lu == NULL) return;
printf("lower:\n");
print_matrix(lu_lower(lu));
printf("\nupper:\n");
print_matrix(lu_upper(lu));
	printf("\npermutation:\n");
	printf("[");
	for(i = 0; i < rows_matrix(lu->_lower); i++)
	{
		if(i > 0) printf(", ");
		printf("%d", lu->_permutation[i]);
	}
	printf("]\n");
}

/* END */
