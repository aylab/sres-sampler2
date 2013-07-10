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

#ifndef STRUCTS_H
#define STRUCTS_H

#include <cstring>
#include <iostream>
#include <fstream>

#include "ESES.h"

#include "memory.h"

using namespace std;

// Stores values and shortcut functions for coloring terminal output and printing common messages (set -c or --no-color to disable colors)
struct terminal {
	// Escape codes
	const char* code_blue;
	const char* code_red;
	const char* code_reset;
	
	// Colors
	char* blue;
	char* red;
	char* reset;
	
	// Verbose stream
	streambuf* verbose_streambuf;
	ostream* verbose_stream;
	
	terminal () {
		this->code_blue = "\x1b[34m";
		this->code_red = "\x1b[31m";
		this->code_reset = "\x1b[0m";
		this->blue = (char*)mallocate(strlen(this->code_blue) + 1);
		this->red = (char*)mallocate(strlen(this->code_red) + 1);
		this->reset = (char*)mallocate(strlen(this->code_reset) + 1);
		this->verbose_stream = new ostream(cout.rdbuf());
	}
	
	~terminal () {
		free(this->blue);
		free(this->red);
		free(this->reset);
		delete verbose_stream;
	}
	
	// Used to indicate a process is done
	void done (ostream& stream) {
		stream << this->blue << "Done" << this->reset << endl;
	}
	
	void done () {
		done(cout);
	}
	
	// Used to indicate the program is out of memory
	void no_memory () {
		cout << this->red << "Not enough memory!" << this->reset << endl;
	}
	
	// Used to indicate the program couldn't read from a pipe
	void failed_pipe_read () {
		cout << this->red << "Couldn't read from the pipe!" << this->reset << endl;
	}
	
	// Used to indicate the program couldn't read from a pipe
	void failed_pipe_write () {
		cout << this->red << "Couldn't write to the pipe!" << this->reset << endl;
	}
	
	ostream& verbose () {
		return *(this->verbose_stream);
	}
	
	void set_verbose_streambuf (streambuf* sb) {
		this->verbose_stream->rdbuf(sb);
		this->verbose_streambuf = this->verbose_stream->rdbuf();
	}
};

// Stores all of the input parameters passed via command-line
struct input_params {
	int pop_parents; // The population of parent simulations to use each generation, default=30
	int pop_children; // The population of child simulations to use each generation, default=200
	int generations; // The number of generations to run before returning results, default=1
	int seed; // The seed used in the evolutionary strategy, default=current UNIX time
	
	char** sim_args; // Arguments to be passed to the simulation
	
	bool verbose; // Whether or not the program is verbose, i.e. prints many messages about program and simulation state
	bool quiet; // Whether or not the program is quiet, i.e. redirects cout to /dev/null
	streambuf* cout_orig; // cout's original buffer to be restored at program completion
	ofstream* null_stream; // A stream to /dev/null that cout is redirected to if quiet mode is set
	
	input_params () {
		this->pop_parents = 30;
		this->pop_children = 200;
		this->generations = 1;
		this->seed = time(0);
		this->sim_args = NULL;
		this->null_stream = new ofstream("/dev/null");
	}
	
	~input_params () {
		if (sim_args != NULL) {
			for (int i = 0; sim_args[i] != NULL; i++) {
				free(sim_args[i]);
			}
			free(sim_args);
		}
		delete this->null_stream;
	}
};

// Stores SRES structs and variables
struct sres_params {
	ESParameter* param;
	ESPopulation* population;
	ESStatistics* stats;
	double pf;
	
	sres_params () {
		this->param = NULL;
		this->population = NULL;
		this->stats = NULL;
		this->pf = 0;
	}
};

#endif

