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

#include <cmath>

#include "init.h"
#include "macros.h"
#include "main.h"

using namespace std; 

terminal* term; 

char* copy_str (const char* str) { 
	if (str == NULL) {
		return NULL;
	} else {
		char* newstr = (char*)mallocate(sizeof(char) * (strlen(str) + 1));
		return strcpy(newstr, str);
	}
}

void init_terminal () {
	term = new terminal();
	if (term->blue == NULL || term->red == NULL || term->reset == NULL) {
		term->no_memory();
		exit(EXIT_MEMORY_ERROR);
	}
	strcpy(term->blue, term->code_blue);
	strcpy(term->red, term->code_red);
	strcpy(term->reset, term->code_reset);
}

void free_terminal () {
	cout << term->reset;
	delete term;
}

void accept_input_params (int num_args, char** args, input_params& ip) {
	ip.argc = num_args;
	ip.argv = args;
	
	if (num_args > 1) { // if arguments were given and each argument option is followed by a value
		for (int i = 1; i < num_args; i += 2) { // iterate through each argument pair
			char* option = args[i];
			char* value;
			if (i < num_args - 1) {
				value = args[i + 1];
			} else {
				value = NULL;
			}
			
			/*
			 Check for each possible argument option and overrides the default value for each specified option. If the option isn't recognized or the value given for an option doesn't appear valid then the usage information for the program is printed with an error message and no simulations are run. The code should be fairly self-explanatory with a few exceptions:
			 1) atoi converts a string to an integer, atof converts a string to a floating point number (i.e. rational)
			 2) strings should always be compared using strcmp, not ==, and strcmp returns 0 if the two strings match
			 3) usage(true) prints the usage information with an error message while usage(false) prints it without one
			*/
			
			if (option_set(option, "-d", "--dimensions")) {
				ensure_nonempty(option, value);
				ip.num_dims = atoi(value);
				if (ip.num_dims < 1) {
					usage("The simulation must include a positive number of dimensions. Set -d or --dimensions to at least 1.");
				}
			} else if (option_set(option, "-P", "--parent-population")) {
				ensure_nonempty(option, value);
				ip.pop_parents = atoi(value);
				if (ip.pop_parents < 1) {
					usage("The parent population must be at least one simulation. Set -P or --parent-population to at least 1.");
				}
			} else if (option_set(option, "-p", "--child-population")) {
				ensure_nonempty(option, value);
				ip.pop_children = atoi(value);
				if (ip.pop_children < 1) {
					usage("The child population must be at least one simulation. Set -p or --child-population to at least 1.");
				}
			} else if (option_set(option, "-g", "--generations")) {
				ensure_nonempty(option, value);
				ip.generations = atoi(value);
				if (ip.generations < 1) {
					usage("The population must exist for at least one generation. Set -g or --generations to at least 1.");
				}
			} else if (option_set(option, "-s", "--seed")) {
				ensure_nonempty(option, value);
				ip.seed = atof(value);
				if (ip.seed <= 0) {
					usage("The seed to generate random numbers must be a positive integer. Set -s or --seed to at least 1.");
				}
			} else if (option_set(option, "-f", "--simulation")) {
				ensure_nonempty(option, value);
				mfree(ip.sim_path);
				ip.sim_path = copy_str(value);
			} else if (option_set(option, "-a", "--arguments")) {
				ensure_nonempty(option, value);
				++i;
				ip.num_sim_args = num_args - i + 8;
				ip.sim_args = (char**)mallocate(sizeof(char*) * (ip.num_sim_args));
				ip.sim_args[0] = copy_str("deterministic");
				for (int j = 1; j < ip.num_sim_args - 7; j++) {
					char* arg = args[i + j - 1];
					ip.sim_args[j] = (char*)mallocate(sizeof(char) * (strlen(arg) + 1));
					sprintf(ip.sim_args[j], "%s", arg);
				}
				ip.sim_args[ip.num_sim_args - 7] = copy_str("--pipe-in");
				ip.sim_args[ip.num_sim_args - 6] = copy_str("0");
				ip.sim_args[ip.num_sim_args - 5] = copy_str("--pipe-out");
				ip.sim_args[ip.num_sim_args - 4] = copy_str("0");
				ip.sim_args[ip.num_sim_args - 3] = copy_str("--gradients-file");
				ip.sim_args[ip.num_sim_args - 2] = NULL;
				ip.sim_args[ip.num_sim_args - 1] = NULL;
				i = num_args;
			} else if (option_set(option, "-c", "--no-color")) {
				term->blue = copy_str("");
				term->red = copy_str("");
				term->reset = copy_str("");
				i--;
			} else if (option_set(option, "-v", "--verbose")) {
				if (!ip.verbose) {
					ip.verbose = true;
				}
				i--;
			} else if (option_set(option, "-q", "--quiet")) {
				if (!ip.quiet) {
					ip.quiet = true;
					ip.cout_orig = cout.rdbuf();
					cout.rdbuf(ip.null_stream->rdbuf());
					term->set_verbose_streambuf(ip.null_stream->rdbuf());
				}
				i--;
			} else if (option_set(option, "-h", "--help")) {
				usage("");
				i--;
			} else if (option_set(option, "-l", "--licensing")) { 
				licensing();
				i--;
			} else {
				usage("One of the given command-line arguments is not a valid option. Please check that every argument matches one available in the following usage information.");
			}
		}
	}
}

inline bool option_set (const char* option, const char* short_name, const char* long_name) {
	return strcmp(option, short_name) == 0 || strcmp(option, long_name) == 0;
}

void ensure_nonempty (const char* flag, const char* arg) {
	if (arg == NULL) {
		char* message = (char*)mallocate(strlen("Missing the argument for the '' flag.") + strlen(flag) + 1);
		sprintf(message, "Missing the argument for the '%s' flag.", flag);
		usage(message);
	}
}

void init_sim_args (input_params& ip) {
	if (ip.num_sim_args == 0) {
		ip.num_sim_args = 8;
		ip.sim_args = (char**)mallocate(sizeof(char*) * ip.num_sim_args);
		ip.sim_args[0] = copy_str("deterministic");
		for (int i = 1; i < ip.num_sim_args; i++) {
			ip.sim_args[i] = NULL;
		}
	}
}

char** copy_args (char** args, int num_args) {
	char** new_args = (char**)mallocate(sizeof(char*) * num_args);
	for (int i = 0; i < num_args; i++) {
		new_args[i] = copy_str(args[i]);
	}
	return new_args;
}

void store_pipe (char** args, int index, int pipe) {
	mfree(args[index]);
	int int_size = log10(pipe > 0 ? pipe : 1) + 1;
	args[index] = (char*)mallocate(sizeof(char) * (int_size + 1));
	sprintf(args[index], "%d", pipe);
}

void reset_cout (input_params& ip) {
	if (ip.quiet) {
		cout.rdbuf(ip.cout_orig);
	}
}

