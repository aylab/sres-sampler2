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

#include <cstdio>
#include <cstdlib>

#include "main.h"
#include "macros.h"
#include "structs.h"
#include "init.h"
#include "sres.h"

using namespace std;

extern terminal* term;

input_params ip;

int main (int argc, char** argv) {
	for (int i = 0; i < argc; i++) {
		cerr << argv[i] << " ";
	}
	cerr << endl;
	exit(0);
	init_terminal();
	accept_input_params(argc, argv, ip);
	init_sim_args(ip);
	sres_params sp;
	init_sres(ip, sp);
	run_sres(sp);
	free_sres(sp);
	free_terminal();
	return 0;
}

void usage (const char* message) {
	cout << endl;
	bool error = strlen(message) != 0;
	if (error) {
		cout << term->red << message << term->reset << endl << endl;
	}
	cout << "Usage: [-option [value]]. . . [--option [value]]. . ." << endl;
	cout << "-P, --parent-population [int]        : the population of parent simulations to use each generation, min=1, default=30" << endl;
	cout << "-p, --child-population  [int]        : the population of child simulations to use each generation, min=1, default=200" << endl;
	cout << "-g, --generations       [int]        : the number of generations to run before returning results, min=1, default=1000" << endl;
	cout << "-s, --seed              [int]        : the seed used in the evolutionary strategy (not simulations), min=1, default=time" << endl;
	cout << "-a, --arguments         [N/A]        : every argument following this will be sent to the deterministic simulation" << endl;
	cout << "-c, --no-color          [N/A]        : disable coloring the terminal output, default=unused" << endl;
	cout << "-v, --verbose           [N/A]        : print detailed messages about the program state" << endl;
	cout << "-q, --quiet             [N/A]        : hide the terminal output, default=unused" << endl;
	cout << "-l, --licensing         [N/A]        : view licensing information (no simulations will be run)" << endl;
	cout << "-h, --help              [N/A]        : view usage information (i.e. this)" << endl;
	cout << endl << term->blue << "Example: ./sres-sampler " << term->reset << endl << endl;
	if (error) {
		exit(EXIT_INPUT_ERROR);
	} else {
		exit(EXIT_SUCCESS);
	}
}

void licensing () {
	cout << endl;
	cout << "Stochastically ranked evolutionary strategy sampler for zebrafish segmentation" << endl;
	cout << "Copyright (C) 2013 Ahmet Ay (aay@colgate.edu), Jack Holland (jholland@colgate.edu), Adriana Sperlea (asperlea@colgate.edu), Sebastian Sangervasi (ssangervasi@colgate.edu)" << endl;
	cout << "This program comes with ABSOLUTELY NO WARRANTY" << endl;
	cout << "This is free software, and you are welcome to redistribute it under certain conditions;" << endl;
	cout << "You can use this code and modify it as you wish under the condition that you refer to the article: ???" << endl;
	cout << endl;
	exit(EXIT_SUCCESS);
}

