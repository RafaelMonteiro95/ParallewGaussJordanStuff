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
#include <unistd.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>

#include "matrix.h"

inline void ApocalipseTrombose(const char *msg, int mpi_error) {
	fprintf(stderr, "%s\n", msg);
	MPI_Abort(MPI_COMM_WORLD, mpi_error);
	MPI_Finalize();
	exit(1);
}

int main(int argc, char *argv[]) {

	int i, j;
	int r, c;
	int nproc, rank;
	int cols, pline, prevPivotLine;
	
	double *recv = NULL;
	double *sendVec;
	double *pivotline;
	
	Matrix *matrix = NULL;

	/* Initialization */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#ifdef DEBUG

	char *tmp = (char *) malloc(sizeof(char)*100000);
	strcpy(tmp, "logs/log-p");
	char nprocAscii[10];
	sprintf(nprocAscii, "%d", rank);

	FILE **logs = (FILE **) malloc(sizeof(FILE *)*nproc);
	logs[rank] = fopen(strcat(strcat(tmp, nprocAscii), "-log"), "w");
#endif

	/* Read matrix */
	if(rank == 0) {

		if(isatty(STDIN_FILENO)) printf("Rows and cols: ");
		scanf("%d%d", &r, &c);

		if(nproc != r) {
			ApocalipseTrombose("Number of processes must be the same number of rows in matrix", MPI_ERR_ASSERT);
		}

		matrix = CreateMatrix(r, c);
		
		for(i = 0; i < matrix->rows; i++) {
			for(j = 0; j < matrix->cols; j++) {

				scanf("%lf", &(matrix->values[i][j]));
			}
		}

		MPI_Bcast(&matrix->cols, 1, MPI_INT, rank, MPI_COMM_WORLD);
		cols = matrix->cols;
	} else {
		MPI_Bcast(&cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
		fprintf(logs[rank], "[debug #%d]: cols: %d\n", rank, cols);
	}

	recv = (double *) malloc(cols * sizeof(double));
	pivotline = (double *) malloc(cols * sizeof(double));
	prevPivotLine = 0;

	// For each col in matrix
	for(i = 0; i < cols; i++) {
		
		/* Master only */
		if(rank == 0) {

			int size = cols*matrix->rows;

			#ifdef DEBUG
				fprintf(logs[rank], "[debug #0]: Searching pivot in col: %d\n", i);
				PrintMatrix(matrix);
			#endif

			/* Find pivot - Use OpenMP here */
			pline = FindPivot(matrix, i, prevPivotLine);

			#ifdef DEBUG
				if(pline == -1) fprintf(logs[rank], "[debug #0]: Pivot not found\n");
				else fprintf(logs[rank], "[debug #0]: Pivot line: %d\n", pline);
			#endif
			
			// Pivot not found (all values are 0 or matrix is reduced), 
			// skip to next column
			if(pline == -1) continue;
			prevPivotLine++;

			/* Reduce pivot line (divide line by pivot) - Use OpenMP here */
			// Its easies to first divide pivot line and then swap it
			MultiplyLineByScalar(matrix, pline, 1.0/matrix->values[pline][i]);
				
			#ifdef DEBUG
				fprintf(logs[rank], "[debug #0]: Multplying matrix[%d] by 1/%lf\n", pline, matrix->values[pline][i]);
				PrintMatrix(matrix);
			#endif

			/* Send pivot and their line (indexed by rank) to each slaves */
			MPI_Bcast(matrix->values[pline], cols, MPI_DOUBLE, rank, MPI_COMM_WORLD);
			#ifdef DEBUG
				fprintf(logs[rank], "[debug #0]: Broadcasting: {");
				for(int k = 0; k < cols; k++) {
					fprintf(logs[rank], "%lf, ", matrix->values[pline][k]);
				} 
				fprintf(logs[rank], "\b\b}\n");
			#endif
			
			/* Position pivot */
			SwapLines(matrix, pline, i);
			sendVec = ToArray(matrix, &size); // Create array after matrix operations
			
			#ifdef DEBUG
				fprintf(logs[rank], "[debug #0]: Swapping lines %d and %d\n", pline, i);	
				PrintMatrix(matrix);
			#endif
			
      		MPI_Scatter(&sendVec[cols], cols, MPI_DOUBLE, recv, cols,
      											MPI_DOUBLE, rank, MPI_COMM_WORLD);
			#ifdef DEBUG
				fprintf(logs[rank], "\n[debug #0]: Scattering: {");
				for(int k = cols; k < cols*matrix->rows; k++) {
					fprintf(logs[rank], "%lf, ", sendVec[k]);
				} 
				fprintf(logs[rank], "\b\b} with %d elements per chunk\n", cols);
			#endif

			free(sendVec);

		// Slaves
		} else {

			// Each iteration discards one process
			if(rank < nproc - prevPivotLine) 
				break;

			/* Receive message (pivot and rank lines) */
			MPI_Bcast(pivotline, cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
			#ifdef DEBUG
				fprintf(logs[rank], "[debug #%d]: Received broadcast: {", rank);
				for(int k = 0; k < cols; k++) {
					fprintf(logs[rank], "%lf, ", pivotline[k]);
				} 
				fprintf(logs[rank], "\b\b}\n");
			#endif
			
      		MPI_Scatter(NULL, cols, MPI_DOUBLE, recv, cols,
      								MPI_DOUBLE, 0, MPI_COMM_WORLD);
			#ifdef DEBUG
				fprintf(logs[rank], "\n[debug #%d]: Received scatter: {", rank);
				for(int k = 0; k < cols; k++) {
					fprintf(logs[rank], "%lf, ", recv[k]);
				} 
				fprintf(logs[rank], "\b\b} with %d elements per chunk\n", cols);
			#endif

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
		if(rank == 0) {

		}
	}

#ifdef DEBUG
	free(tmp);
	fclose(logs[rank]);
#endif

	MPI_Finalize();
	free(recv);
	free(pivotline);

	return 0;
}