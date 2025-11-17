/*
 * Copyright (c) 2025 Ismael Mosquera Rivera
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
#include "eigen.h"
#include "svd.h"

#define LEFT_SIDE 0 /* left side */
#define RIGHT_SIDE 1 /* right side */

/*
* Helper functions.
*/

/* get min ( LEFT or RIGHT ). */
static int get_min(int m, int n)
{
return (m <= n) ? LEFT_SIDE : RIGHT_SIDE;
}

/* end helper functions */


/* implementation */

void destroy_svd(SVD* svd)
{
if(svd == NULL) return;
if(svd_u(svd) != NULL) destroy_matrix(svd_u(svd));
if(svd_sigma(svd) != NULL) destroy_matrix(svd_sigma(svd));
if(svd_v(svd) != NULL) destroy_matrix(svd_v(svd));
free(svd);
svd = NULL;
}

void print_svd(const SVD* svd)
{
if(svd == NULL)
{
printf("[]");
return;
}
printf("U:\n");
print_matrix(svd_u(svd));
printf("\n");
printf("Sigma:\n");
print_matrix(svd_sigma(svd));
printf("\n");
printf("V:\n");
print_matrix(svd_v(svd));
printf("\n");
}

SVD* svd_factorization(const Matrix* m)
{
int i, n, min;
SVD* svd = NULL;
EigenSystem* left = NULL;
EigenSystem* right = NULL;
Matrix* u = NULL;
Matrix* sigma = NULL;
Matrix* v = NULL;
if(m == NULL) return NULL;
u = create_matrix(rows_matrix(m), rows_matrix(m));
sigma = create_matrix(rows_matrix(m), columns_matrix(m));
v = create_matrix(columns_matrix(m), columns_matrix(m));
/* compute left and right eigen */
left = eigen_system(mul_matrix(m, transpose_matrix(m)));
right = eigen_system(mul_matrix(transpose_matrix(m), m));

/* compute UDV */
n = size_eigensystem(left);
for(i = 0; i < n; i++)
{
u = set_column_matrix(u, normalize_vector(eigen_vector(left->_eigen[i])), i);
}
n = size_eigensystem(right);
for(i = 0; i < n; i++)
{
v = set_column_matrix(v, normalize_vector(eigen_vector(right->_eigen[i])), i);
}
min = get_min(rows_matrix(m), columns_matrix(m));
if(min == LEFT_SIDE)
{
	n = size_eigensystem(left);
	for(i = 0; i < n; i++)
	{
set_matrix(sigma, sqrt(eigen_value(left->_eigen[i])), i, i);
}
}
else
{
	n = size_eigensystem(right);
		for(i = 0; i < n; i++)
		{
	set_matrix(sigma, sqrt(eigen_value(right->_eigen[i])), i, i);
	}
}
/* build SVD */
svd = (SVD*)malloc(sizeof(SVD));
svd_u(svd) = clone_matrix(u);
svd_sigma(svd) = clone_matrix(sigma);
svd_v(svd) = clone_matrix(v);

/* release previously allocated memory */
destroy_eigensystem(left);
destroy_eigensystem(right);
destroy_matrix(u);
destroy_matrix(sigma);
destroy_matrix(v);

return svd;
}

/* END */

