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

inline void kill(int mpi_error){
	MPI_Abort(MPI_COMM_WORLD, mpi_error);
	MPI_Finalize();
	exit(1);
}

int main(int argc, char *argv[]){

	int i, j;
	int r, c;
	int nproc, rank, pline;
	
	int *recv = NULL;
	int *sendVec;
	int *pivotline;
	
	Matrix *matrix = NULL;

	/* Initialization */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* Read matrix */
	if(rank == 0){
		
		printf("Rows and cols: ");
		scanf("%d%d", &r, &c);
		matrix = CreateMatrix(r, c);
		
		for(i = 0; i < matrix->rows; i++){
			for(j = 0; j < matrix->cols; j++){

				scanf("%lf", &(matrix->values[i][j]));
			}
		}
	}

	recv = (int *) malloc(matrix->cols * sizeof(int));
	pivot = (int *) malloc(matrix->cols * sizeof(int));
	
	// For each row
	for(i = 0; i < matrix->rows; i++){
		
		/* Master only */
		if(rank == 0){

			sendVec = ToArray(matrix, matrix->cols*matrix->rows);

			#ifdef DEBUG
				printf("Searching pivot in col: %d\n", i);
				PrintMatrix(matrix);
			#endif

			/* Find pivot - Use OpenMP here */
			pline = FindPivot(matrix, i);

			#ifdef DEBUG
				if(pline == -1) printf("Pivot not found\n");
				else printf("Pivot line: %d\n", pline);
			#endif
			
			// Pivot not found (all values are 0 or matrix is reduced), 
			// skip to next column
			if(pline == -1) continue;

			/* Reduce pivot line (divide line by pivot) - Use OpenMP here */
			// Its easies to first divide pivot line and then swap it
			MultiplyLineByScalar(matrix, pline, 1.0/matrix->values[pline][i]);
				
			#ifdef DEBUG
				printf("Multplying matrix by 1/%d\n", matrix->values[pline][i]);
				PrintMatrix(matrix);
			#endif

			/* Position pivot */
			SwapLines(matrix, pline, i);
			PrintMatrix(matrix);
			
			#ifdef DEBUG
				printf("Swapping lines %d and %d\n", pline, i);	
				PrintMatrix(matrix);
			#endif

			/* Send pivot and their line (indexed by rank) to each slaves */
			MPI_Bcast(&matrix->values[pline], 1, MPI_INT, rank, MPI_COMM_WORLD);
      		MPI_Scatter(&sendVec[1], matrix->cols, MPI_DOUBLE, recv, matrix->cols,
      											MPI_DOUBLE, rank, MPI_COMM_WORLD);

      		free(sendVec);

		// Slaves
		} else {
			/* Receive message (pivot and rank lines) */
			MPI_Bcast(pivotline, 1, MPI_INT, 0, MPI_COMM_WORLD);
      		MPI_Scatter(&sendVec[1], matrix->cols, MPI_DOUBLE, recv, matrix->cols,
      											   MPI_DOUBLE, 0, MPI_COMM_WORLD);

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
	free(pivotline);

	return 0;
}