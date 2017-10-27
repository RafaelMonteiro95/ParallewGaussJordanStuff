/* Trabalho 2 - Gauss Jordan Matrix MPI
	Lucas Alexandre Soares 		9293265
	Giovanna Oliveira Guimarães 9293692
	Rafael Augusto Monteiro		9293095
	Choyoung Francisco Lim 		6436060
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>

#include "matrix.h"

typedef enum { false, true } bool;

void kill(int mpi_error){
	MPI_Abort(MPI_COMM_WORLD, mpi_error);
	MPI_Finalize();
	exit(1);
}

void usage(){

	printf("Usage: ./")	
}

int main(int argc, char *argv[]){

	int nproc, rank;
	Matrix *matrix = NULL;
	
	int *recv = NULL;

	int *sendVec;
	int *displacement;

	if(argc != 3){
		usage();
		return 0;
	}

	matrix = CreateMatrix(atoi(argv[1]), atoi(argv[2]));

	/* Initialization */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank == 0){
		i = 0;

		/* Get input */

	}


	if(rank == 0){
		printf("Found number %d %d times\n", number, result);
		free(vector);
	}

	MPI_Finalize();
	free(recv);
	free(sendVec);
	free(displacement);

	return 0;
}