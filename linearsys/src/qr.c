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
#include <math.h>
#include "qr.h"

/* implementation */

void destroy_qr(QR* qr)
{
if(qr == NULL) return;
if(qr_q(qr) != NULL) free(qr_q(qr));
if(qr_r(qr) != NULL) free(qr_r(qr));
if(qr != NULL) free(qr);
}

QR* qr_factorization(const Matrix* m)
{
QR* qr = NULL;
Matrix* q = NULL;
Matrix* q1 = NULL;
Matrix* r = NULL;
double w;
int i, j, k, n;

if(rows_matrix(m) < 2) return NULL; /* order must be at least 2 */
if(rows_matrix(m) != columns_matrix(m)) return NULL; /* m must be square */
r = clone_matrix(m);
n = columns_matrix(m);
k = 0; /* sentinel */
for(i = 0; i < n-1; i++)
{
for(j = n-1; j > i; j--)
{
q1 = identity_matrix(n);
w = atan2(-get_matrix(r, j, i), get_matrix(r, i, i));
set_matrix(q1, cos(w), i, i);
set_matrix(q1, -sin(w), i, j);
set_matrix(q1, sin(w), j, i);
set_matrix(q1, cos(w), j, j);
r = mul_matrix(q1, clone_matrix(r));
q = (k == 0) ? clone_matrix(q1) : mul_matrix(q1, clone_matrix(q));
k++;
}
}
qr = (QR*)malloc(sizeof(QR));
qr_q(qr) = transpose_matrix(q);
qr_r(qr) = clone_matrix(r);

destroy_matrix(q);
destroy_matrix(q1);
destroy_matrix(r);

return qr;
}

void print_qr(const QR* qr)
{
if(qr == NULL)
{
printf("\n[]\n");
return;
}
printf("\nq:\n");
print_matrix(qr_q(qr));
printf("\nr:\n");
print_matrix(qr_r(qr));
}

double qr_det(const QR* qr)
{
Matrix* r = NULL;
int i;
double d = 1.0;

if(qr == NULL) return 0.0;

r = clone_matrix(qr_r(qr));
for(i = 0; i < columns_matrix(r); i++)
{
d *= get_matrix(r, i, i);
}
destroy_matrix(r);
return d;
}


/* END */

