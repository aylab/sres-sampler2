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

#include "memory.h"
#include "macros.h"
#include "structs.h"

extern terminal* term;

#if defined(MEMTRACK)
	size_t heap_current = 0;
	size_t heap_total = 0;
#endif

void* mallocate (size_t size) {
	void* block;
	#if defined(MEMTRACK)
		block = malloc(sizeof(size_t) + size);
	#else
		block = malloc(size);
	#endif
	if (block == NULL) {
		term->no_memory();
		exit(EXIT_MEMORY_ERROR);
	}
	#if defined(MEMTRACK)
		heap_current += size;
		heap_total += size;
		size_t* sizeblock = (size_t*)block;
		*sizeblock = size;
		return (void*)(sizeblock + 1);
	#else
		return block;
	#endif
}

void* callocate (size_t elements, size_t size) {
	void* mem = mallocate(elements * size);
	memset(mem, 0, elements * size);
	return mem;
}

void* reallocate (void* mem, size_t size) {
	#if defined(MEMTRACK)
		if (mem == NULL) {
			return mallocate(size);
		} else {
			size_t* sizeblock = (size_t*)mem - 1;
			if (*sizeblock < size) {
				void* newmem = mallocate(size);
				memcpy(newmem, mem, *sizeblock);
				mfree(sizeblock);
				return newmem;
			} else {
				return mem;
			}
		}
	#else
		return realloc(mem, size);
	#endif
}

void mfree (void* mem) {
	#if defined(MEMTRACK)
		if (mem != NULL) {
			size_t* memblock = (size_t*)mem - 1;
			heap_current -= *memblock;
			free(memblock);
		}
	#else
		free(mem);
	#endif
}

void* operator new (size_t size) {
	return mallocate(size);
}

void* operator new[] (size_t size) {
	return mallocate(size);
}

void operator delete (void* mem) {
	mfree(mem);
}

void operator delete[] (void* mem) {
	mfree(mem);
}

#if defined(MEMTRACK)
static void print_mem_amount (size_t mem) {
	static size_t kB = 1024;
	static size_t MB = square(1024);
	static size_t GB = cube(1024);
	
	if (mem > GB) {
		cout << mem / GB << " GB";
	} else if (mem > MB) {
		cout << mem / MB << " MB";
	} else if (mem > kB) {
		cout << mem / kB << " kB";
	} else {
		cout << mem << " B";
	}
	cout << endl;
}

void print_heap_usage () {
	cout << term->blue << "Current heap usage:\t" << term->reset;
	print_mem_amount(heap_current);
	cout << term->blue << "Total heap usage:\t" << term->reset;
	print_mem_amount(heap_total);
}
#endif

