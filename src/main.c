/* Trabalho 2 - Gauss Jordan Matrix MPI
	Lucas Alexandre Soares 		9293265
	Giovanna Oliveira Guimarães 9293692
	Rafael Augusto Monteiro		9293095
	Choyoung Francisco Lim 		6436060

	NOTE:
	Modelo do código: https://docs.google.com/document/d/129-iAEgjICFppIovFr41jjV9S2ORI6DCErZhCD0yePE/edit?usp=sharing
	Implementar seguindo o passo a passo modelado
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

void usage(const char *progName){

	printf("Usage: %s [nrows] [ncols]", progName);
}

int main(int argc, char *argv[]){

	int nproc, rank;
	Matrix *matrix = NULL;
	
	int i;
	int *recv = NULL;

	int *sendVec;
	int *displacement;

	if(argc != 3){
		usage(argv[1]);
		return 0;
	}

	// Get nrows and ncols from command line
	matrix = CreateMatrix(atoi(argv[1]), atoi(argv[2]));

	/* Initialization */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* Read matrix */
	if(rank == 0){

	}

	// For each column
	for(i = 0; i < matrix->cols; i++){
		
		/* Master only */
		if(rank == 0){

			/* Find pivot - Use OpenMP here */
			// FindPivot(matrix);

			/* Position pivot */
			// SwapLines(matrix, line1, line2);

			/* Reduce pivot line (divide line by pivot) - Use OpenMP here */
			// MultiplyLineByScalar(matrix, line, value);

			/* Send pivot and their line (indexed by rank) to each slaves */
			// NOTE: dá pra usar aqueles tipos de dados q o psergio usou no
			// exemplo da multiplicação de matriz (vetor) pra separar cada linha
			// da matriz e mandar tudo de uma vez

		// Slaves
		} else {
			/* Receive message (pivot and rank lines) */

		}

		/* Sum each line (indexed by rank) with the pivot line multiplied by a 
		scalar. The scalar shall be the oposite of the element in the same 
		column as the pivot	of that line, i.e, if the pivot line is [0, 1, 2, 3] 
		and the  line [0, 2, 3, 4] should be solved by multiplying the pivot 
		line by -2 and then adding these two lines: 
			 [0, 1, 2, 3]*(-2)
			+[0, 2, 3, 4]
			---------------
			 [0, 0, -1, -2]

			Use OpenMP here.
		*/
		// MultiplyLineByScalar(matrix, line, value);
		// AddLines(matrix, line1, line2);

		/* Send result back to master */

		/* Master can debug print */
		if(rank == 0){

		}
	}

	MPI_Finalize();
	free(recv);
	free(sendVec);
	free(displacement);

	return 0;
}