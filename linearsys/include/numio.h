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

#ifndef ___NUMIO_H___
#define ___NUMIO_H___

#ifdef __cplusplus
extern "C" {
	#endif

#include <stdio.h>

/*
* Basic helper functions to read and write int and double values from and to a file.
*/

/*
* Get a pointer to an open file to read.
* param: file name.
* returns:
* A pointer to read from a file.
*/
FILE* fopen_read(const char* filename);

/*
* Get a pointer to an open file to write.
* param: file name.
* returns:
* A pointer to write to a file.
*/
FILE* fopen_write(const char* filename);

/*
* Read the next int value from a file.
* param: FILE* pf => pointer to a file.
* param: int* value => value to read.
*/
void fread_int(FILE* pf, int* value);

/*
* Write an int value to a file.
* param: FILE* pf => pointer to a file.
* param: int value => value to write.
*/
void fwrite_int(FILE* pf, int value);

/*
* Read the next double value from a file.
* param: FILE* pf => pointer to a file.
* param: double* value => value to read.
*/
void fread_double(FILE* pf, double* value);

/*
* Write a double value to a file.
* param: FILE* pf => pointer to a file.
* param: double value => value to write.
*/
void fwrite_double(FILE* pf, double value);

#ifdef __cplusplus
}
#endif

#endif
