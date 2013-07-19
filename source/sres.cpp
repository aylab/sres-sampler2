/*
Stochastically ranked evolutionary strategy sampler for zebrafish segmentation
Copyright (C) 2013 Ahmet Ay, Jack Holland, Adriana Sperlea, Sebastian Sangervasi

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.	If not, see <http://www.gnu.org/licenses/>.
*/

#include <ctime>
#include <stdio.h>

#if defined(MPI)
	#include "../libsres-mpi/sharefunc.h"
	#include "../libsres-mpi/ESSRSort.h"
	#include "../libsres-mpi/ESES.h"
#else
	#include "../libsres/sharefunc.h"
	#include "../libsres/ESSRSort.h"
	#include "../libsres/ESES.h"
#endif

#include "sres.h"
#include "macros.h"
#include "io.h"

extern terminal* term;

void init_sres (input_params& ip, sres_params& sp) {
	int es = esDefESSlash;
	int constraint = 0;
	int dim = ip.num_dims;
	int miu = ip.pop_parents;
	int lambda = ip.pop_children;
	int gen = ip.generations;
	double gamma = esDefGamma;
	double alpha = esDefAlpha;
	double varphi = esDefVarphi;
	int retry = 0;
	sp.pf = essrDefPf;
	
	sp.trsfm = (ESfcnTrsfm*)mallocate(sizeof(ESfcnTrsfm) * dim);
	
	for (int i = 0; i < dim; i++) {
		sp.trsfm[i] = transform;
	}
	
	sp.lb = (double*)mallocate(sizeof(double) * dim);
	sp.ub = (double*)mallocate(sizeof(double) * dim);
	double* lb = sp.lb;
	double* ub = sp.ub;
	for (int i = 0; i < dim; i++) {
		lb[i] = 0;
		ub[i] = 0;
	}
	
	if (dim == 27) {
		lb[0] = 30,		ub[0] = 65;
		lb[1] = 30,		ub[1] = 65;
		lb[2] = 30,		ub[2] = 65;
		lb[3] = 30,		ub[3] = 65;
		lb[4] = 0.1,	ub[4] = 0.46;
		lb[5] = 0.1,	ub[5] = 0.46;
		lb[6] = 0.1,	ub[6] = 0.46;
		lb[7] = 0.1,	ub[7] = 0.46;
		lb[8] = 10,		ub[8] = 60;
		lb[9] = 10,		ub[9] = 60;
		lb[10] = 10,	ub[10] = 60;
		lb[11] = 20,	ub[11] = 60;
		lb[12] = 0.11,	ub[12] = 0.35;
		lb[13] = 0.11,	ub[13] = 0.35;
		lb[14] = 0.11,	ub[14] = 0.35;
		lb[15] = 0.07,	ub[15] = 0.35;
		lb[16] = 8,		ub[16] = 12;
		lb[17] = 8,		ub[17] = 12;
		lb[18] = 0,		ub[18] = 0;
		lb[19] = 6,		ub[19] = 12;
		lb[20] = 0.4,	ub[20] = 2;
		lb[21] = 0.4,	ub[21] = 2;
		lb[22] = 0.4,	ub[22] = 2;
		lb[23] = 10,	ub[23] = 25;
		lb[24] = 150,	ub[24] = 900;
		lb[25] = 150,	ub[25] = 900;
		lb[26] = 200,	ub[26] = 800;
	} else if (dim == 45) {
		lb[0] = 30,		ub[0] = 65;
		lb[1] = 30,		ub[1] = 65;
		lb[2] = 30,		ub[2] = 65;
		lb[3] = 30,		ub[3] = 65;
		lb[4] = 0.1,	ub[4] = 0.46;
		lb[5] = 0.1,	ub[5] = 0.46;
		lb[6] = 0.1,	ub[6] = 0.46;
		lb[7] = 0.1,	ub[7] = 0.46;
		lb[8] = 10,		ub[8] = 60;
		lb[9] = 10,		ub[9] = 60;
		lb[10] = 10,	ub[10] = 60;
		lb[11] = 20,	ub[11] = 60;
		lb[12] = 0.11,	ub[12] = 0.35;
		lb[13] = 0.11,	ub[13] = 0.35;
		lb[14] = 0.11,	ub[14] = 0.35;
		lb[15] = 0.07,	ub[15] = 0.35;
		lb[16] = 0.0003,ub[16] = 0.03;
		lb[17] = 0.0003,ub[17] = 0.03;
		lb[18] = 0.0003,ub[18] = 0.03;
		lb[19] = 0.0003,ub[19] = 0.03;
		lb[20] = 0.0003,ub[20] = 0.03;
		lb[21] = 0.0003,ub[21] = 0.03;
		lb[22] = 0.003,	ub[22] = 0.3;
		lb[23] = 0.003,	ub[23] = 0.3;
		lb[24] = 0.003,	ub[24] = 0.3;
		lb[25] = 0.003,	ub[25] = 0.3;
		lb[26] = 0.003,	ub[26] = 0.3;
		lb[27] = 0.003,	ub[27] = 0.3;
		lb[28] = 0.11,	ub[28] = 0.35;
		lb[29] = 0.11,	ub[29] = 0.35;
		lb[30] = 0.11,	ub[30] = 0.35;
		lb[31] = 0.11,	ub[31] = 0.35;
		lb[32] = 0.11,	ub[32] = 0.35;
		lb[33] = 0.11,	ub[33] = 0.35;
		lb[34] = 8,		ub[34] = 12;
		lb[35] = 8,		ub[35] = 12;
		lb[36] = 0,		ub[36] = 0;
		lb[37] = 6,		ub[37] = 12;
		lb[38] = 0.4,	ub[38] = 2;
		lb[39] = 0.4,	ub[39] = 2;
		lb[40] = 0.4,	ub[40] = 2;
		lb[41] = 10,	ub[41] = 25;
		lb[42] = 150,	ub[42] = 900;
		lb[43] = 150,	ub[43] = 900;
		lb[44] = 200,	ub[44] = 800;
	} else if (dim == 65) {
		lb[0] = 30, 	ub[0] = 65;
		lb[1] = 30, 	ub[1] = 65;
		lb[2] = 30, 	ub[2] = 65;
		lb[3] = 30, 	ub[3] = 65;
		lb[4] = 30, 	ub[4] = 65;
		lb[5] = 0.1, 	ub[5] = 0.46;
		lb[6] = 0.1, 	ub[6] = 0.46;
		lb[7] = 0.1, 	ub[7] = 0.46;
		lb[8] = 0.1, 	ub[8] = 0.46;
		lb[9] = 0.1, 	ub[9] = 0.46;
		lb[10] = 10, 	ub[10] = 60;
		lb[11] = 10, 	ub[11] = 60;
		lb[12] = 10, 	ub[12] = 60;
		lb[13] = 10, 	ub[13] = 60;
		lb[14] = 20, 	ub[14] = 60;
		lb[15] = 0.11, 	ub[15] = 0.35;
		lb[16] = 0.11, 	ub[16] = 0.35;
		lb[17] = 0.11, 	ub[17] = 0.35;
		lb[18] = 0.11, 	ub[18] = 0.35;
		lb[19] = 0.07, 	ub[19] = 0.35;
		lb[20] = 0.0003,ub[20] = 0.03;
		lb[21] = 0.0003,ub[21] = 0.03;
		lb[22] = 0.0003,ub[22] = 0.03;
		lb[23] = 0.0003,ub[23] = 0.03;
		lb[24] = 0.0003,ub[24] = 0.03;
		lb[25] = 0.0003,ub[25] = 0.03;
		lb[26] = 0.0003,ub[26] = 0.03;
		lb[27] = 0.0003,ub[27] = 0.03;
		lb[28] = 0.0003,ub[28] = 0.03;
		lb[29] = 0.0003,ub[29] = 0.03;
		lb[30] = 0.003,	ub[30] = 0.3;
		lb[31] = 0.003,	ub[31] = 0.3;
		lb[32] = 0.003,	ub[32] = 0.3;
		lb[33] = 0.003,	ub[33] = 0.3;
		lb[34] = 0.003,	ub[34] = 0.3;
		lb[35] = 0.003,	ub[35] = 0.3;
		lb[36] = 0.003,	ub[36] = 0.3;
		lb[37] = 0.003,	ub[37] = 0.3;
		lb[38] = 0.003,	ub[38] = 0.3;
		lb[39] = 0.003,	ub[39] = 0.3;
		lb[40] = 0.11,	ub[40] = 0.35;
		lb[41] = 0.11,	ub[41] = 0.35;
		lb[42] = 0.11,	ub[42] = 0.35;
		lb[43] = 0.11,	ub[43] = 0.35;
		lb[44] = 0.11,	ub[44] = 0.35;
		lb[45] = 0.11,	ub[45] = 0.35;
		lb[46] = 0.11,	ub[46] = 0.35;
		lb[47] = 0.11,	ub[47] = 0.35;
		lb[48] = 0.11,	ub[48] = 0.35;
		lb[49] = 0.11,	ub[49] = 0.35;
		lb[50] = 8,		ub[50] = 12;
		lb[51] = 8,		ub[51] = 12;
		lb[52] = 8,		ub[52] = 12;
		lb[53] = 0,		ub[53] = 0;
		lb[54] = 6,		ub[54] = 12;
		lb[55] = 0.4,	ub[55] = 2;
		lb[56] = 0.4,	ub[56] = 2;
		lb[57] = 0.4,	ub[57] = 2;
		lb[58] = 0.4,	ub[58] = 2;
		lb[59] = 10,	ub[59] = 25;
		lb[60] = 150,	ub[60] = 900;
		lb[61] = 150,	ub[61] = 900;
		lb[62] = 1000,	ub[62] = 9000;
		lb[63] = 1000,	ub[63] = 9000;
		lb[64] = 200,	ub[64] = 800;
	} else {
		cout << term->red << "The given number of dimensions does not have ranges programmed in! Please check that the given number (" << dim << ") is correct or add ranges to sres.cpp." << term->reset << endl;
	}
	
	ESInitial( // (AAy: ??? Initialize the ES method)
	#if defined(MPI)
		&(ip.argc), &(ip.argv),
	#endif
		ip.seed, &(sp.param), sp.trsfm, fitness, es, constraint, dim, ub, lb, miu, lambda, gen, gamma, alpha, varphi, retry, &(sp.population), &(sp.stats));
}

void run_sres (sres_params& sp) { // (AAy: ??? Run the SRES method)
	while (sp.stats->curgen < sp.param->gen) {
		ESStep(sp.population, sp.param, sp.stats, sp.pf);
	}
}

void free_sres (sres_params& sp) { // (AAy: ??? Free SRES arguments)
	mfree(sp.trsfm);
	mfree(sp.lb);
	mfree(sp.ub);
	ESDeInitial(sp.param, sp.population, sp.stats);
}

void fitness (double* parameters, double* score, double* constraints) { // (AAy: ??? Fitness function)
	*score = simulate_set(parameters);
}

double transform (double x) { // (AAy: ??? What is this code doing?)
	return x;
}

