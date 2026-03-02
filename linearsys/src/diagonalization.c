/*
 * Copyright (c) 2026 Ismael Mosquera Rivera
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
#include "eigen.h"
#include "diagonalization.h"

/*
* Helper private function to check for eigen value multiplicity.
*/
static int multiplicity(const double* v, int n)
{
int i, j, retval;
double current;
retval = 0;
for(i = 0; i < n-1; i++)
{
current = v[i];
for(j = i+1; j < n; j++)
{
	if(v[j] == current)
	{
		retval = 1;
		break;
	}
}
}
return retval;
}

/* End helper functions */

/* Implementation */

void destroy_diagonalization(Diagonalization* diag)
{
if(diag == NULL) return;
if(diagonalization_p(diag) != NULL) destroy_matrix(diagonalization_p(diag));
if(diagonalization_d(diag) != NULL) destroy_matrix(diagonalization_d(diag));
if(diagonalization_pt(diag) != NULL) destroy_matrix(diagonalization_pt(diag));
free(diag);
diag = NULL;
}

Diagonalization* clone_diagonalization(const Diagonalization* diag)
{
Diagonalization* out = (Diagonalization*)malloc(sizeof(Diagonalization));
diagonalization_p(out) = clone_matrix(diagonalization_p(diag));
diagonalization_d(out) = clone_matrix(diagonalization_d(diag));
diagonalization_pt(out) = clone_matrix(diagonalization_pt(diag));
return out;
}

Diagonalization* diagonalize_matrix(const Matrix* m)
{
int i, n;
double* x = NULL;
Matrix* p = NULL;
Matrix* d = NULL;
EigenSystem* eigsys = NULL;
Diagonalization* diag = NULL;
if(m == NULL) return NULL;
n = rows_matrix(m);
if(n != columns_matrix(m)) return NULL; /* m must be square */
eigsys = eigen_system(m);
if(eigsys == NULL) return NULL; /* fat chance */
n = size_eigensystem(eigsys);
x = (double*)malloc(n * sizeof(double));
for(i = 0; i < n; i++) x[i] = eigen_value(eigen_eigensystem(eigsys)[i]);
if(multiplicity(x, n)) return NULL;
p = create_matrix(n, n);
d = create_matrix(n, n);
for(i = 0; i < n; i++)
{
set_matrix(d, x[i], i, i);
p = set_column_matrix(p, normalize_vector(eigen_vector(eigen_eigensystem(eigsys)[i])), i);
}
/* create and build diagonalization */
diag = (Diagonalization*)malloc(sizeof(Diagonalization));
diagonalization_p(diag) = clone_matrix(p);
diagonalization_d(diag) = clone_matrix(d);
diagonalization_pt(diag) = transpose_matrix(p);

/* release previously allocated memory */
free(x);
destroy_matrix(p);
destroy_matrix(d);
destroy_eigensystem(eigsys);

return diag; /* done */
}

double det_diagonalized_matrix(const Diagonalization* diag)
{
	double d;
	int i, n;
if(diag == NULL) return NaN;
n = rows_matrix(diagonalization_d(diag));
d = 1.0;
for(i = 0; i < n; i++)
{
	d *= get_matrix(diagonalization_d(diag), i, i);
}
return d;
}

void print_diagonalization(const Diagonalization* diag)
{
if(diag == NULL)
{
printf("[]\n");
return;
}
printf("Diagonalization:\n");
printf("\n");
printf("P:\n");
print_matrix(diagonalization_p(diag));
printf("\n");
printf("D:\n");
print_matrix(diagonalization_d(diag));
printf("\n");
printf("P^-1:\n");
print_matrix(diagonalization_pt(diag));
printf("\n");
}


/* END */
