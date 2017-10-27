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

typedef enum { false, true } bool;

void kill(int mpi_error){
	MPI_Abort(MPI_COMM_WORLD, mpi_error);
	MPI_Finalize();
	exit(1);
}

void printvec(int *v, int n){
	
	int i;
	
	printf("v = {");
	for(i = 0; i < n; i++)
		printf("%d, ", v[i]);
	printf("\b\b}\n");
}

int main(int argc, char *argv[]){

	int nproc, rank;
	int i = 0;
	int number, repetitions = 0;
	int *vector = (int *) malloc(sizeof(int)*VEC_SIZE);
	int *recv = NULL;

	int *sendVec;
	int *displacement;

	/* Initialization */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int count = VEC_SIZE/nproc;
	int rest = VEC_SIZE%nproc;
	sendVec = (int *) malloc(sizeof(int)*nproc);
	displacement = (int *) malloc(sizeof(int)*nproc);
	
	if(rank == 0){
		i = 0;

		/* Get input */
	}


	recv = (int *) malloc(sizeof(int)*VEC_SIZE/nproc);

	// Send vector - scatterv support non integer workload division
	MPI_Scatterv(vector, sendVec, displacement, MPI_INTEGER, 
				 recv, sendVec[rank], MPI_INTEGER, 0, MPI_COMM_WORLD);
	// Send vector - this scatter must have a divisible vector
	/*MPI_Scatter(vector, count, MPI_INTEGER, recv, count, 
				MPI_INTEGER, 0, MPI_COMM_WORLD);*/
	
	// No need for barrier

	// Start searching
	#pragma omp parallel for reduction(+:repetitions)
	for(i = 0; i < sendVec[rank]; i++){
		if(recv[i] == number){ // Found number, increment repetition counter
			repetitions++;
#ifdef debug
			fprintf(stderr, "[debug] Found number %d %d times on process %d\n", 
					number, repetitions, rank);
#endif
		}
	}

	// Reduce repetitions
	int result;
	MPI_Reduce(&repetitions, &result, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD);

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