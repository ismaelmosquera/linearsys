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

#ifndef ___VECTOR_H___
#define ___VECTOR_H___

#ifdef __cplusplus
extern "C" {
	#endif

/*
* Vector type definition.
*/
typedef struct
{
int _size;
double* _data;
}Vector;

/*
* Creates a Vector with all the components initialized to 0.
* param: int size => size of the Vector ( number of elements ).
* returns:
* A pointer to a new created Vector.
*/
Vector* create_vector(int size);

/*
* Destroys the Vector passed as parameter.
* param: Vector*v => a Vpointer to a Vector to destroy.
*/
void destroy_vector(Vector* v);

/*
* Gets a clone of the Vector passed as parameter.
* param: const Vector* v => a pointer to a Vector to clone.
* returns:
* A pointer to a clone of the Vector passed as parameter.
*/
Vector* clone_vector(const Vector* v);

/*
* Resizes a Vector.
* param: const Vector* v => a Vector to resize.
* param: int new_size => size for the resized Vector.
*
* returns:
* NULL if the new size is less than 1.
* a Vector with positions 0..bew_size-1 if the new size is less than the size of the Vector passed as parameter.
* a clone of the vector passed as parameter if new size equals size(v).
* a Vector with positions 0..size(v)-1 and the new positions set to 0 if new_size greater than size(v).
*/
Vector* resize_vector(const Vector* v, int new_size);

/*
* Gets a chunk from the Vector passed as parameter.
* param: const Vector* v => a Vector from to get a chunck.
* param: int from_index => first index of the chunk.
* param: int to_index => last index of the chunk.
*
* returns:
* NULL if some index is out of range.
* a Vector with range from_index .. to_index from the Vector passed as parameter.
*/
Vector* get_vector_chunk(const Vector* v, int from_index, int to_index);

/*
* Loads a Vector from a file.
* param: const char* filename => name of the file.
* returns:
* A pointer to a Vector loaded from the file passed as parameter.
*
* Format of the file:
*
* n
* v0 v1 v2 ... vn
*
* where n is the size ( number of elements in the vector )
* and v0 v1 v2 ... are the components of the vector.
*/
Vector* load_vector(const char* filename);

/*
* Writes a Vector to a file.
* param: const Vector* v => vector to store.
* param: const char* filename => the name of the file to store the vector.
*
* The format of the file is the same as above in the load_vector function.
*/
void store_vector(const Vector* v, const char* filename);

/*
* Prints a Vector to the console.
* param: const Vector* v => vector to be printed.
*
* The format to print must be as follows:
* [v0, v1, v2, ... vn]
*
* Each double value has to have a decimal precision of 2 digits.
* modifier to print just 2 decimal digits: "%.2lf"
*
*/
void print_vector(const Vector* v);

/*
* Gets a value from a Vector.
* param: const Vector* v => pointer to a Vector.
* param: int pos => position from to get the value.
*
* returns:
* double value in position 'pos' from Vector 'v'.
*/
double get_vector(const Vector* v, int pos);

/*
* Sets a value to a Vector.
* param: Vector*v => a pointer to a Vector.
* param: double value => value to set in the vector.
* param: int pos => position to set in the vector.
*/
void set_vector(Vector* v, double value, int pos);

/*
* Vector addition.
* param: const Vector* v1 => a pointer to a vector.
* param: const Vector* v2 => a pointer to a vector.
*
* returns:
* A pointer to a Vector result of the addition.
* size v1 must be equal to size v2
*/
Vector* add_vector(const Vector* v1, const Vector* v2);

/*
* Vector substraction.
* param: const Vector* v1 => a pointer to a vector.
* param: const Vector* v2 => a pointer to a vector.
*
* returns:
* A pointer to a Vector result of the substraction.
* size v1 must be equal to size v2
*/
Vector* sub_vector(const Vector* v1, const Vector* v2);

/*
* Multiply all components of a Vector by a scalar.
* param: const Vector* v => a pointer to a vector.
* param: double value => a scalar.
*
* returns:
* A pointer to a Vector result of the product.
*/
Vector* mul_vector(const Vector* v, double value);

/*
* Divide all components of a Vector by a scalar.
* param: const Vector* v => a pointer to a vector.
* param: double value => a scalar.
*
* returns:
* A pointer to a Vector result of the division.
*/
Vector* div_vector(const Vector* v, double value);

/*
* Gets a module ( magnitude ) of a Vector.
* param: const Vector* v => a pointer to a Vector.
*
* returns:
* The magnitude of the vector.
*/
double module_vector(const Vector* v);

/*
* Normalize the Vector passed as parameter.
* param: Vector* v => a pointer to a Vector.
* The parameter is passed by reference, so it is modified to get it normalized.
*/
void normalize_vector(Vector* v);

/*
* Gets the dot product of 2 vectors.
* param: const Vector* v1 => a pointer to a Vector.
* param: const Vector* v2 => a pointer to a Vector.
*
* returns:
* The dot product of v1 and v2.
* size v1 must be equal to size v2
*/
double dot_product_vector(const Vector* v1, const Vector* v2);

/*
* Gets the cross product of 2 vectors.
* param: const Vector* v1 => a pointer to a Vector of size 3
* param: const Vector* v2 => a pointer to a Vector of size 3
*
* returns:
* A pointer to a Vector result of the cross product.
* size v1 = size v2 = 3
*/
Vector* cross_product_vector(const Vector* v1, const Vector* v2);

/*
* Macro to get the number of elements ( size ) of a vector.
*/
#define size_vector(vector) ((vector)->_size)

#ifdef __cplusplus
}
#endif

#endif
