/*
Stochastically ranked evolutionary strategy sampler for zebrafish segmentation
Copyright (C) 2013 Ahmet Ay, Jack Holland, Adriana Sperlea, Sebastian Sangervasi

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
memory.h contains function declarations for memory.cpp
*/

#ifndef MEMORY_H
#define MEMORY_H

#include <stdlib.h> // Needed for size_t

// These memory functions have been added to the libSRES code, which is compiled as C code, so 'extern "C"' must be added to prevent C++'s signature mangling from hiding the names from C
#ifdef __cplusplus
extern "C" {
#endif
void* mallocate(size_t);
void* callocate(size_t, size_t);
void* reallocate(void*, size_t);
void mfree(void*);
#ifdef __cplusplus
}
#endif
#if defined(MEMTRACK)
	void print_heap_usage();
#endif

#endif

