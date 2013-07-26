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
io.cpp contains functions for input and output of files and pipes. All I/O related functions should be placed in this file.
*/

#include <sys/wait.h> // Needed for waitpid
#include <unistd.h> // Needed for pipe, read, write, close, fork, execv

#include "io.h"
#include "macros.h"
#include "structs.h"
#include "init.h"

extern terminal* term; // Declared in init.cpp
extern input_params ip; // Declared in main.cpp

/* simulate_set performs the required piping to setup and run a simulation with the given parameters
	parameters:
		parameters: the parameters to pass as a parameter set to the simulation
	returns: the score the simulation received
	notes:
	todo:
*/
double simulate_set (double parameters[]) {
	// Create a pipe
	int pipes[2];
	if (pipe(pipes) == -1) {
		term->failed_pipe_create();
		exit(EXIT_PIPE_CREATE_ERROR);
	}
	
	// Copy the user-specified simulation arguments and fill the copy with the pipe's file descriptors
	char** sim_args = copy_args(ip.sim_args, ip.num_sim_args);
	store_pipe(sim_args, ip.num_sim_args - 4, pipes[0]);
	store_pipe(sim_args, ip.num_sim_args - 2, pipes[1]);
	
	// Fork the process so the child can run the simulation
	pid_t pid = fork();
	if (pid == -1) {
		term->failed_fork();
		exit(EXIT_FORK_ERROR);
	}
	if (pid == 0) { // The child runs the simulation
		if (execv(ip.sim_path, sim_args) == -1) {
			term->failed_exec();
			exit(EXIT_EXEC_ERROR);
		}
	} else { // The parent pipes in the parameter set to run
		write_pipe(pipes[1], parameters);
	}
	
	// Wait for the child to finish simulating
	int status = 0; 
	waitpid(pid, &status, WUNTRACED);
	if (WIFEXITED(status) == 0) {
		term->failed_child();
		exit(EXIT_CHILD_ERROR);
	}
	
	// Close the writing end of the pipe
	if (close(pipes[1]) == -1) {
		term->failed_pipe_write();
		exit(EXIT_PIPE_WRITE_ERROR);
	}
	
	// Pipe in the simulation's score
	int max_score;
	int score;
	read_pipe(pipes[0], &max_score, &score);
	
	// Close the reading end of the pipe
	if (close(pipes[0]) == -1) {
		term->failed_pipe_read();
		exit(EXIT_PIPE_WRITE_ERROR);
	}
	
	// Free the simulation arguments
	for (int i = 0; sim_args[i] != NULL; i++) {
		mfree(sim_args[i]);
	}
	mfree(sim_args);
	
	// libSRES requires scores from 0 to 1 with 0 being a perfect score so convert the simulation's score format into libSRES's
	return 1 - ((double)score / max_score);
}

/* write_pipe writes the given parameter set to the given pipe
	parameters:
		fd: the file descriptor of the pipe to write to
		parameters: the parameter set to pipe
	returns: nothing
	notes:
	todo:
*/
void write_pipe (int fd, double parameters[]) {
	write_pipe_int(fd, ip.num_dims); // Write the number of dimensions, i.e. parameters per set, being sent
	write_pipe_int(fd, 1); // Write that one parameter set is being sent
	if (write(fd, parameters, sizeof(double) * ip.num_dims) == -1) {
		term->failed_pipe_write();
		exit(EXIT_PIPE_WRITE_ERROR);
	}
}

/* write_pipe_int writes the given integer to the given pipe
	parameters:
		fd: the file descriptor of the pipe to write to
		value: the integer to pipe
	returns: nothing
	notes:
	todo:
*/
void write_pipe_int (int fd, int value) {
	if (write(fd, &value, sizeof(int)) == -1) {
		term->failed_pipe_write();
		exit(EXIT_PIPE_WRITE_ERROR);
	}
}

/* read_pipe reads the maximum score and the received score from the given pipe
	parameters:
		fd: the file descriptor of the pipe to write to
		max_score: a pointer to store the maximum score the simulation could have received
		score: a pointer to store the score the simulation actually received
	returns: nothing
	notes:
	todo:
*/
void read_pipe (int fd, int* max_score, int* score) {
	read_pipe_int(fd, max_score);
	read_pipe_int(fd, score);
}

/* read_pipe_int writes an integer from the given pipe
	parameters:
		fd: the file descriptor of the pipe to write to
		address: a pointer to store the received integer
	returns: nothing
	notes:
	todo:
*/
void read_pipe_int (int fd, int* address) {
	if (read(fd, address, sizeof(int)) == -1) {
		term->failed_pipe_read();
		exit(EXIT_PIPE_READ_ERROR);
	}
}

