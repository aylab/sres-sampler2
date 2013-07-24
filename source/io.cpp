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
#include <fstream>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include "io.h"
#include "macros.h"
#include "structs.h"
#include "init.h"

extern terminal* term;
extern input_params ip;

double simulate_set (double parameters[]) {
	int pipes[2];
	if (pipe(pipes) == -1) {
		term->failed_pipe_create();
		exit(EXIT_PIPE_CREATE_ERROR);
	}
	
	pid_t pid = fork();
	if (pid == -1) {
		term->failed_fork();
		exit(EXIT_FORK_ERROR);
	}
	
	if (pid == 0) {
		pid = getpid();
		int pid_strlen = log10(pid > 0 ? pid : 1) + 1;
		char* grad_fname = (char*)mallocate(sizeof(char) * (strlen("input-.gradients") + pid_strlen + 1));
		sprintf(grad_fname, "input-%d.gradients", pid);
		ofstream grad_file(grad_fname);
		grad_file << (int)parameters[0] << " (7 1) (" << (int)parameters[1] << " " << (int)parameters[2] << ")" << endl;
		grad_file.close();
		
		char** sim_args = copy_args(ip.sim_args, ip.num_sim_args);
		store_pipe(sim_args, ip.num_sim_args - 6, pipes[0]);
		store_pipe(sim_args, ip.num_sim_args - 4, pipes[1]);
		sim_args[ip.num_sim_args - 2] = grad_fname;
		
		if (execv(ip.sim_path, sim_args) == -1) {
			term->failed_exec();
			exit(EXIT_EXEC_ERROR);
		}
	} else {
		double par_set[45] = {56.796071,32.909131,56.996694,62.820661,0.352096,0.330384,0.101593,0.144759,50.563904,28.176744,47.758870,53.972291,0.330960,0.158198,0.320945,0.204394,0.027315,0.027341,0.026838,0.021798,0.019523,0.020656,0.289651,0.292611,0.126843,0.006280,0.072729,0.150319,0.272693,0.167980,0.122839,0.150364,0.342542,0.120830,11.600145,9.575987,0.000000,8.224160,0.719038,1.282242,0.407751,10.482257,703.313280,153.899249,365.043112};
		write_pipe(pipes[1], par_set);
	}
	
	int status = 0; 
	waitpid(pid, &status, WUNTRACED);
	if (WIFEXITED(status) == 0) {
		term->failed_child();
		exit(EXIT_CHILD_ERROR);
	}
	if (close(pipes[1]) == -1) {
		term->failed_pipe_write();
		exit(EXIT_PIPE_WRITE_ERROR);
	}
	
	int max_score;
	int score;
	read_pipe(pipes[0], &max_score, &score);
	if (close(pipes[0]) == -1) {
		term->failed_pipe_read();
		exit(EXIT_PIPE_WRITE_ERROR);
	}
	
	return 1 - ((double)score / max_score);
}

void write_pipe (int fd, double parameters[]) {
	write_pipe_int(fd, 45);
	write_pipe_int(fd, 1);
	if (write(fd, parameters, sizeof(double) * ip.num_dims) == -1) {
		term->failed_pipe_write();
		exit(EXIT_PIPE_WRITE_ERROR);
	}
}

void write_pipe_int (int fd, int value) {
	if (write(fd, &value, sizeof(int)) == -1) {
		term->failed_pipe_write();
		exit(EXIT_PIPE_WRITE_ERROR);
	}
}

void read_pipe (int fd, int* max_score, int* score) {
	read_pipe_int(fd, max_score);
	read_pipe_int(fd, score);
}

void read_pipe_int (int fd, int* address) {
	if (read(fd, address, sizeof(int)) == -1) {
		term->failed_pipe_read();
		exit(EXIT_PIPE_READ_ERROR);
	}
}

