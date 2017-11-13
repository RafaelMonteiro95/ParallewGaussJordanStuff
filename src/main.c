/* Trabalho 2 - Gauss Jordan Matrix MPI
	Lucas Alexandre Soares 		9293265
	Giovanna Oliveira Guimarães 9293692
	Rafael Augusto Monteiro		9293095
	Choyoung Francisco Lim 		6436060
	NOTE:
	Modelo do código: https://docs.google.com/document/d/129-iAEgjICFppIovFr41jjV9S2ORI6DCErZhCD0yePE/edit?usp=sharing
	Implementar seguindo o passo a passo modelado
	https://matrix.reshish.com/ todos os casos
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>

#include "matrix.h"

inline void kill(const char *msg, int mpi_error) {
	fprintf(stderr, "%s\n", msg);
	MPI_Abort(MPI_COMM_WORLD, mpi_error);
	MPI_Finalize();
	exit(1);
}

int main(int argc, char *argv[]) {

	int i, j;
	int r, c;
	int nproc, rank;
	int pcol, prow;
	int chunkSize, lastChunkSize;
	int cols, prevPivotLine;

	time_t time = clock();
	
	double scalar;
	double *chunk;
	double *pivotRow = NULL;
	double *backupRow = NULL;
	double *currentRow = NULL; 

	FILE* matrixFp;
	FILE* vectorFp;

	Matrix *matrix = NULL;

	/* Initialization */
	omp_set_num_threads(atoi(argv[3]));
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int displs[nproc], sendcounts[nproc];

#ifdef DEBUG

	char *tmp = (char *) malloc(sizeof(char)*1000);
	strcpy(tmp, "logs/log-p");
	char nprocAscii[10];
	sprintf(nprocAscii, "%d", rank);

	FILE **logs = (FILE **) malloc(sizeof(FILE *)*nproc);
	logs[rank] = fopen(strcat(strcat(tmp, nprocAscii), "-log"), "w");
	free(tmp);
#endif

	/* Read matrix */
	if(rank == 0){

		// opening files with equation system
		matrixFp = fopen(argv[1], "r");
		vectorFp = fopen(argv[2], "r");

		if(!matrixFp || !vectorFp) {
			char buf[255];
			sprintf(buf, "Failed to open %s.\n", matrixFp ? "matriz.txt" : "vetor.txt");
			kill(buf, 1);
		}

		// getting size of matrix by counting the number of '\n' in the file.
		char cbuff;
		c = 0, r = 0;
		while(!feof(matrixFp)){
			cbuff = fgetc(matrixFp);
			if(cbuff == '\n'){
				c++;
				r++;
			}
		}
		
		// I'm using c+1 because there's the results column which isn't counted
		rewind(matrixFp);

		// Creating the matrix to be reduced.
		matrix = CreateMatrix(r, c+1);
		
		for(i = 0; i < matrix->rows; i++){
			for(j = 0; j < matrix->cols - 1; j++){
				fscanf(matrixFp, "%lf", &(matrix->values[mat2vec(matrix->cols, i, j)]));
			}
		}

		// VectorFp contains our matrix last column
		for(i = 0; i < matrix->rows; i++){
			fscanf(vectorFp, "%lf", &(matrix->values[mat2vec(matrix->cols, i, c)]));
		}
	
		// Sending matrix size
		c++;
		MPI_Bcast(&r, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// Calculates chunk size
		backupRow = (double *) malloc(sizeof(double) * c);
		chunkSize = r/nproc; // Integer division
		
		// If not divisible, last chunk is smaller
		lastChunkSize = r - chunkSize*(nproc-1);
		chunk = (double *) malloc(sizeof(double)*chunkSize);

		#ifdef DEBUG
			PrintMatrix(matrix);
		#endif
			
		fclose(matrixFp);
		fclose(vectorFp);

		// Scatter chunks
		int sum = 0;
		for(int i = 0; i < nproc-1; i++){
			
			if(chunkSize == 0 && i < r)
				sendcounts[i] = c;
			else 
				sendcounts[i] = chunkSize*c;
			
			displs[i] = sum;
			sum += sendcounts[i];
		}
		if(chunkSize != 0)
			sendcounts[nproc-1] = lastChunkSize*c;
		displs[nproc-1] = sum;

	} else {

		MPI_Bcast(&r, 1, MPI_INT, 0, MPI_COMM_WORLD);
		c = r+1; // Number of cols is always rows+1
	
		// Calculates chunk size
		backupRow = (double *) malloc(sizeof(double) * c);
		chunkSize = r/nproc; // Integer division
		lastChunkSize = r - chunkSize*(nproc-1); // If not divisible, last chunk is smaller
				
		// Scatter chunks
		int sum = 0;
		for(int i = 0; i < nproc-1; i++){
			
			if(chunkSize == 0 && i < r)
				sendcounts[i] = c;
			else 
				sendcounts[i] = chunkSize*c;
			
			displs[i] = sum;
			sum += sendcounts[i];
		}
		if(chunkSize != 0)
			sendcounts[nproc-1] = lastChunkSize*c;
		displs[nproc-1] = sum;

		chunk = (double *) malloc(sizeof(double)*sendcounts[rank]);			
	}
	
	// For each row in the matrix
	for(i = 0; i < r; i++){
		
		/* Master only */
		if(rank == 0) {

			#ifdef DEBUG
				printf("Searching pivot in col: %d\n", i);
				PrintMatrix(matrix);
			#endif

			/* Find pivot - Use OpenMP here */
			pcol = i;
			prow = FindPivot(matrix, pcol);

			#ifdef DEBUG
				if(prow == -1) printf("Pivot not found\n");
				else printf("Pivot line: %d (value: %lf)\n", prow, matrix->values[mat2vec(c, prow, pcol)]);
			#endif
			
			// Pivot not found (all values are 0 or matrix is reduced), 
			// skip to next column
			if(prow == -1) continue;

			#ifdef DEBUG
				printf("Multplying matrix by 1/%lf\n", matrix->values[mat2vec(c, prow, i)]);
				// PrintMatrix(matrix);
			#endif
			
			/* Reduce pivot line (divide line by pivot) - Use OpenMP here */
			// Its easiest to first divide pivot line and then swap it
			MultiplyLineByScalar(matrix, prow, 1.0/matrix->values[mat2vec(c, prow, i)]);

			/* Send pivot and their line (indexed by rank) to each slaves */
			#ifdef DEBUG
				fprintf(logs[rank], "[debug #0]: Broadcasting: {");
				for(int k = 0; k < cols; k++)
					fprintf(logs[rank], "%lf, ", matrix->values[mat2vec(c, prow, k)]);
				fprintf(logs[rank], "\b\b}\n\n");
			#endif

			/* Swap pivot line */
			SwapLines(matrix, prow, i);
			prow = i;

			pivotRow = matrix->values + (prow * c);
			chunk = matrix->values;

			#ifdef DEBUG
				FPrintMatrix(matrix, logs[rank]);
			#endif
		
		} /* End master */

		else { /* Slaves */
			if(!pivotRow) pivotRow = (double *) malloc(sizeof(double)*c);
		} /* End slaves */

		
		/* Everyone */
		
		// Broadcast pivot line
		MPI_Bcast(pivotRow, c, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		// Scatter each process chunk
		double *tmp = NULL;
		if(rank == 0) tmp = matrix->values;
		
		MPI_Scatterv(tmp, sendcounts, displs,
            MPI_DOUBLE, chunk, sendcounts[rank],
            MPI_DOUBLE, 0, MPI_COMM_WORLD);

		#ifdef DEBUG
			fprintf(logs[rank], "[debug#%d] chunk: {", rank);
			for(int i = 0; i < sendcounts[rank];){
				fprintf(logs[rank], "%lf, ", chunk[i]);
				if(++i%c == 0) fprintf(logs[rank], "\n");
			} fprintf(logs[rank], "\b\b}\n");
		#endif

		#ifdef DEBUG
			fprintf(logs[rank], "[debug#%d] (i: %d) pivot row = {", rank, i);
			for (int k = 0; k < c; k++){
				fprintf(logs[rank], "%lf, ", pivotRow[k]);
			} fprintf(logs[rank], "\b\b}\n");
		#endif

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

		// Searching for pcol
		int k;
		for(k = 0; k < c; k++){
			if(pivotRow[k] == 1.0){
				pcol = k;
				prow = k;
				break;
			}
		}

		// For each row in the chunk
		for(j = 0; j < sendcounts[rank]/c; j++){

			#ifdef DEBUG
				// printf("Swapping lines %d and %d\n", prow, i);	
				// PrintMatrix(matrix);
			#endif

			// I don't need to process the pivot row
			if(j+displs[rank]/c == prow) continue;

			/* we have to multiply the pivot row by a value 
			that would zero the value in this current row that is
			below or above the current pivot.*/

			// calculating the scalar value of line product
			scalar = -chunk[(j*c)+pcol];
			
			// Creating a auxiliar vector to store the prod. value
			memcpy(backupRow, pivotRow, sizeof(double) * c);

			#ifdef DEBUG
				fprintf(logs[rank], "[debug#%d] scalar: %lf\n", rank, scalar);
				for (int k = 0; k < c; k++){
					fprintf(logs[rank], "[debug#%d] pivotrow: %lf\n", rank, pivotRow[k]);
				}for (int k = 0; k < c; k++){
					fprintf(logs[rank], "[debug#%d] backuprow: %lf\n", rank, backupRow[k]);
				}
			#endif

			// Multipying row by scalar
			_MultiplyLineByScalar(backupRow, c, scalar);

			// Sum of currently selected row and multiplied row
			currentRow = chunk + (c*j);
			_AddLines(currentRow, backupRow, c);
		}

		/* Send result back to master */
		tmp = NULL;
		if(rank == 0) tmp = matrix->values;

		#ifdef DEBUG
			fprintf(logs[rank], "[debug#%d]: Gathering chunk: {", rank);
			for(int k = 0; k < sendcounts[rank]; k++)
				fprintf(logs[rank], "%lf, ", chunk[k]);
			fprintf(logs[rank], "\b\b}\n\n");
		#endif

		MPI_Gatherv(chunk, sendcounts[rank], MPI_DOUBLE,
            tmp, sendcounts, displs, MPI_DOUBLE,
            0, MPI_COMM_WORLD);

	}
	
	/* Master can print */
	if(rank == 0) {
		FILE *result_fp = fopen("resultado.txt", "w");
		if(!result_fp){
			printf("Failed to open result file, printing to stdout\n");
			result_fp = stdout;
		}
		FPrintMatrix(matrix, result_fp);
		if(result_fp != stdout) fclose(result_fp);

		printf("time: %lf", (double) (clock()-time)/CLOCKS_PER_SEC/atoi(argv[3]));
	}

	#ifdef DEBUG
		fclose(logs[rank]);
		free(logs);
	#endif

	free(backupRow);
	free(chunk);
	MPI_Finalize();

	return 0;
}