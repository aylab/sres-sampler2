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
	
	if (dim == 3) {
		lb[0] = 0,	ub[0] = 45;
		lb[1] = 0,	ub[1] = 3;
		lb[2] = 5,	ub[2] = 300;
	} else {
		cout << term->red << "The given number of dimensions does not have ranges programmed in! Please check that the given number (" << dim << ") is correct or add ranges to sres.cpp." << term->reset << endl;
	}
	
	ESInitial(
	#if defined(MPI)
		&(ip.argc), &(ip.argv),
	#endif
		ip.seed, &(sp.param), sp.trsfm, fitness, es, constraint, dim, ub, lb, miu, lambda, gen, gamma, alpha, varphi, retry, &(sp.population), &(sp.stats));
}

void run_sres (sres_params& sp) {
	while (sp.stats->curgen < sp.param->gen) {
		ESStep(sp.population, sp.param, sp.stats, sp.pf);
	}
}

void free_sres (sres_params& sp) {
	mfree(sp.trsfm);
	mfree(sp.lb);
	mfree(sp.ub);
	ESDeInitial(sp.param, sp.population, sp.stats);
}

void fitness (double* parameters, double* score, double* constraints) {
	parameters[0] = (int)parameters[0];
	if (parameters[1] < 1) {
		parameters[1] = 23;
	} else if (parameters[1] < 2) {
		parameters[1] = 31;
	} else {
		parameters[1] = 39;
	}
	if (parameters[2] <= 5) {
		parameters[2] = 5;
	} else if (parameters[2] <= 25) {
		parameters[2] = 25;
	} else if (parameters[2] <= 50) {
		parameters[2] = 50;
	} else if (parameters[2] <= 150) {
		parameters[2] = 150;
	} else if (parameters[2] <= 200) {
		parameters[2] = 200;
	} else {
		parameters[2] = 300;
	}
	*score = simulate_set(parameters);
}

double transform (double x) {
	return x;
}

