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

#include "numio.h"

FILE* fopen_read(const char* filename)
{
return fopen(filename, "r");
}

FILE* fopen_write(const char* filename)
{
return fopen(filename, "w");
}

void fread_int(FILE* pf, int* value)
{
int n = 0;
fscanf(pf, "%d", &n);
*value = n;
}

void fwrite_int(FILE* pf, int value)
{
fprintf(pf, "%d", value);
}

void fread_double(FILE* pf, double* value)
{
double d = 0.0;
fscanf(pf, "%lf", &d);
*value = d;
}

void fwrite_double(FILE* pf, double value)
{
fprintf(pf, "%lf", value);
}

/* END */

