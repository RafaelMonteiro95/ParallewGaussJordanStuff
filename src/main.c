/* Trabalho 2 - Gauss Jordan Matrix MPI
	Lucas Alexandre Soares 		9293265
	Giovanna Oliveira Guimarães 9293692
	Rafael Augusto Monteiro		9293095
	Choyoung Francisco Lim 		6436060
	NOTE:
	Modelo do código: https://docs.google.com/document/d/129-iAEgjICFppIovFr41jjV9S2ORI6DCErZhCD0yePE/edit?usp=sharing
	Implementar seguindo o passo a passo modelado
	http://www.resolvermatrices.com/ para testar casos de teste até 16x16
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
	
	double scalar;
	double *chunk;
	double *pivotRow = NULL;
	double *currentRow = NULL;
	double *backupRow = NULL;

	FILE* matrixFp;
	FILE* vectorFp;

	Matrix *matrix = NULL;

	/* Initialization */
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
#endif

	/* Read matrix */
	if(rank == 0){

		// opening files with equation system
		matrixFp = fopen("matriz.txt", "r");
		vectorFp = fopen("vetor.txt", "r");

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
		lastChunkSize = r - chunkSize*(nproc-1); // If not divisible, last chunk is smaller
		chunk = (double *) malloc(sizeof(double)*chunkSize);

		#ifdef DEBUG
			fprintf(stderr, "[debug] chunkSize: %d\n", chunkSize);
		#endif

		PrintMatrix(matrix);

		fclose(matrixFp);
		fclose(vectorFp);

		// Scatter chunks
		int sum = 0;
		for(int i = 0; i < nproc-1; i++){
			sendcounts[i] = chunkSize*c;
			displs[i] = sum;
			sum += sendcounts[i];
			fprintf(stderr, "[debug#0] displs[%d]: %d\n", i, displs[i]);
		}
		sendcounts[nproc-1] = lastChunkSize*c;
		displs[nproc-1] = sum;
		fprintf(stderr, "[debug#%d] displs[%d]: %d\n", rank, nproc-1, displs[nproc-1]);

		#ifdef DEBUG
			for (int i = 0; i < nproc; i++){
				fprintf(stderr, "[debug#0] sendcounts[%d]: %d\n", i, sendcounts[i]);
			}
		#endif

		MPI_Scatterv(matrix->values, sendcounts, displs,
            MPI_DOUBLE, chunk, chunkSize,
            MPI_DOUBLE, 0, MPI_COMM_WORLD);

	} else {

		MPI_Bcast(&r, 1, MPI_INT, 0, MPI_COMM_WORLD);
		c = r+1; // Number of cols is always rows+1
	
		// Calculates chunk size
		backupRow = (double *) malloc(sizeof(double) * c);
		chunkSize = r/nproc; // Integer division
		lastChunkSize = r - chunkSize*(nproc-1); // If not divisible, last chunk is smaller
		chunk = (double *) malloc(sizeof(double)*chunkSize);
				
		// Scatter chunks
		int sum = 0;
		for(int i = 0; i < nproc-1; i++){
			sendcounts[i] = chunkSize*c;
			displs[i] = sum;
			sum += sendcounts[i];
			fprintf(logs[rank], "[debug#%d] displs[%d]: %d\n", rank, i, displs[i]);
		}
		sendcounts[nproc-1] = lastChunkSize*c;
		displs[nproc-1] = sum;
		fprintf(logs[rank], "[debug#%d] displs[%d]: %d\n", rank, nproc-1, displs[nproc-1]);

		fprintf(logs[rank], "[debug#%d] allocating %d elements for chunk\n", rank, sendcounts[rank]);
		chunk = (double *) malloc(sizeof(double)*sendcounts[rank]);

		MPI_Scatterv(NULL, sendcounts, displs,
            MPI_DOUBLE, chunk, sendcounts[rank],
            MPI_DOUBLE, 0, MPI_COMM_WORLD);

		fprintf(logs[rank], "[debug#%d] sendcounts[%d]: %d\n", rank, rank, sendcounts[rank]);
		fprintf(logs[rank], "[debug#%d] receiving chunk.\n", rank);
		fprintf(logs[rank], "[debug#%d] chunk: {", rank);
		for(int i = 0; i < sendcounts[rank];){
			fprintf(logs[rank], "%lf, ", chunk[i]);
			if(++i%c == 0) fprintf(logs[rank], "\n");
		} fprintf(logs[rank], "\b\b}\n");
	}
	
	free(backupRow);
	MPI_Finalize();

	return 0;

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
				PrintMatrix(matrix);
			#endif
			
			/* Reduce pivot line (divide line by pivot) - Use OpenMP here */
			// Its easiest to first divide pivot line and then swap it
			MultiplyLineByScalar(matrix, prow, 1.0/matrix->values[mat2vec(c, prow, i)]);

			/* Send pivot and their line (indexed by rank) to each slaves */
			#ifdef DEBUG
				fprintf(logs[rank], "[debug #0]: Broadcasting: {");
				for(int k = 0; k < cols; k++)
					fprintf(logs[rank], "%lf, ", matrix->values[mat2vec(c, prow, k)]);
				fprintf(logs[rank], "\b\b}\n");
			#endif

			// MPI_Bcast(matrix->values[prow], cols, MPI_DOUBLE, rank, MPI_COMM_WORLD); // FIXME

			/* Swap pivot line */
			SwapLines(matrix, prow, i);
			prow = i;

			// /* Sends pivot and chunks to each slaves */
			// MPI_Bcast(pivotRow, c, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			// MPI_Scatter(matrix->values, chunkSize * c, MPI_DOUBLE, currentRow, chunkSize * c, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		// Slaves
		} else {

			/* Receive message (pivot and rank lines) */
			// MPI_Bcast(pivotRow, c, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			// MPI_Scatter(matrix->values, chunkSize * c, MPI_DOUBLE, currentRow, chunkSize * c, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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

		/* SEQUENTIAL VERSION OF CODE*/

		// Alocating memory for processing each row
		// Selecting the pivot row
		// pivotRow = matrix->values[prow]; // FIXME

		// For each row in the matrix:
		// #pragma omp parallel for
		for(j = 0; j < r; j++){

			#ifdef DEBUG
				// printf("Swapping lines %d and %d\n", prow, i);	
				// PrintMatrix(matrix);
			#endif

			// I don't need to process the pivot row
			if (j == prow) continue;

			/* we have to multiply the pivot row by a value 
			that would zero the value in this current row that is
			below or above the current pivot.*/

			// calculating the scalar value of line product
			scalar = -matrix->values[mat2vec(c, j, pcol)];
			pivotRow = matrix->values + (prow * c);

			// Creating a auxiliar vector to store the prod. value
			memcpy(backupRow, pivotRow, sizeof(double) * c);

			// Multipying row by scalar
			_MultiplyLineByScalar(backupRow, c, scalar);

			// Sum of currently selected row and multiplied row
			currentRow = matrix->values + (c*j);

			_AddLines(currentRow, backupRow, c);
		}
		
		/* END SEQUENTIAL CODE */

		/* PARALLEL CODE*/


		/* END PARALLEL CODE */

		/* Send result back to master */

		/* Master can debug print */
		if(rank == 0) {
			PrintMatrix(matrix);
		}
	}

#ifdef DEBUG
	free(tmp);
	fclose(logs[rank]);
#endif

	MPI_Finalize();
	free(backupRow);

	return 0;
}