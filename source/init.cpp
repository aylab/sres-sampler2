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
	char* newstr = (char*)mallocate(sizeof(char) * strlen(str) + 1);
	return strcpy(newstr, str);
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
			
			if (strcmp(option, "-P") == 0 || strcmp(option, "--parent-population") == 0) {
				ensure_nonempty(option, value);
				ip.pop_parents = atoi(value);
				if (ip.pop_parents < 1) {
					usage("The parent population must be at least one simulation. Set -P or --parent-population to at least 1.");
				}
			} else if (strcmp(option, "-p") == 0 || strcmp(option, "--child-population") == 0) {
				ensure_nonempty(option, value);
				ip.pop_children = atoi(value);
				if (ip.pop_children < 1) {
					usage("The child population must be at least one simulation. Set -p or --child-population to at least 1.");
				}
			} else if (strcmp(option, "-g") == 0 || strcmp(option, "--generations") == 0) {
				ensure_nonempty(option, value);
				ip.generations = atoi(value);
				if (ip.generations < 1) {
					usage("The population must exist for at least one generation. Set -g or --generations to at least 1.");
				}
			} else if (strcmp(option, "-s") == 0 || strcmp(option, "--seed") == 0) {
				ensure_nonempty(option, value);
				ip.seed = atof(value);
				if (ip.seed <= 0) {
					usage("The seed to generate random numbers must be a positive integer. Set -s or --seed to at least 1.");
				}
			} else if (strcmp(option, "-a") == 0 || strcmp(option, "--arguments") == 0) {
				ensure_nonempty(option, value);
				++i;
				ip.num_sim_args = num_args - i + 6;
				ip.sim_args = (char**)mallocate(sizeof(char*) * (ip.num_sim_args));
				ip.sim_args[0] = copy_str("deterministic");
				for (int j = 1; j < ip.num_sim_args - 5; j++) {
					ip.sim_args[j] = (char*)mallocate(sizeof(char) * (strlen(args[i]) + 1));
					sprintf(ip.sim_args[j], "%s", args[i + j - 1]);
				}
				ip.sim_args[ip.num_sim_args - 5] = copy_str("--pipe-in");
				ip.sim_args[ip.num_sim_args - 4] = copy_str("0");
				ip.sim_args[ip.num_sim_args - 3] = copy_str("--pipe-out");
				ip.sim_args[ip.num_sim_args - 2] = copy_str("0");
				ip.sim_args[ip.num_sim_args - 1] = NULL;
				i = num_args;
			} else if (strcmp(option, "-c") == 0 || strcmp(option, "--no-color") == 0) {
				strcpy(term->blue, "");
				strcpy(term->red, "");
				strcpy(term->reset, "");
				i--;
			} else if (strcmp(option, "-v") == 0 || strcmp(option, "--verbose") == 0) {
				if (!ip.verbose) {
					ip.verbose = true;
				}
				i--;
			} else if (strcmp(option, "-q") == 0 || strcmp(option, "--quiet") == 0) {
				if (!ip.quiet) {
					ip.quiet = true;
					ip.cout_orig = cout.rdbuf();
					cout.rdbuf(ip.null_stream->rdbuf());
					term->set_verbose_streambuf(ip.null_stream->rdbuf());
				}
				i--;
			} else if (strcmp(option, "-h") == 0 || strcmp(option, "--help") == 0) {
				usage("");
				i--;
			} else if (strcmp(option, "-l") == 0 || strcmp(option, "--licensing") == 0) { 
				licensing();
				i--;
			} else {
				usage("One of the given command-line arguments is not a valid option. Please check that every argument matches one available in the following usage information.");
			}
		}
	}
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
		ip.num_sim_args = 6;
		ip.sim_args = (char**)mallocate(sizeof(char*) * 6);
		ip.sim_args[0] = copy_str("deterministic");
		for (int i = 1; i < 6; i++) {
			ip.sim_args[i] = NULL;
		}
	}
}

void store_pipe (input_params& ip, int index, int pipe) {
	free(ip.sim_args[index]);
	int int_size = log10(pipe > 0 ? pipe : 1) + 1;
	ip.sim_args[index] = (char*)mallocate(sizeof(char) * (int_size + 1));
	sprintf(ip.sim_args[index], "%d", pipe);
}

