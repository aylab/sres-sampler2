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
	
	char** sim_args = copy_args(ip.sim_args, ip.num_sim_args);
	store_pipe(sim_args, ip.num_sim_args - 4, pipes[0]);
	store_pipe(sim_args, ip.num_sim_args - 2, pipes[1]);
	
	pid_t pid = fork();
	if (pid == -1) {
		term->failed_fork();
		exit(EXIT_FORK_ERROR);
	}
	if (pid == 0) {
		if (execv(ip.sim_path, sim_args) == -1) {
			term->failed_exec();
			exit(EXIT_EXEC_ERROR);
		}
	} else {
		write_pipe(pipes[1], parameters);
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
	
	for (int i = 0; sim_args[i] != NULL; i++) {
		free(sim_args[i]);
	}
	free(sim_args);
	
	return 1 - ((double)score / max_score);
}

double abs (double val) {
	return val >= 0 ? val : -val;
}

void write_pipe (int fd, double parameters[]) {
	write_pipe_int(fd, ip.num_dims);
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

