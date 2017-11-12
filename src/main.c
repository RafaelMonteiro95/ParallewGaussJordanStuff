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
#include <string.h>
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
	int nproc, rank;
	int chunkSize;
	
	int *recv = NULL;
	int *sendVec;

	FILE* matrixFp;
	FILE* vectorFp;
	
	int prow, pcol;
	double scalar;

	double* pivotRow = NULL;
	double* currentRow = NULL;
	
	Matrix *matrix = NULL;

	/* Initialization */
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* Read matrix */
	if(rank == 0){

		// opening files with equation system
		matrixFp = fopen("matriz.txt","r");
		vectorFp = fopen("vetor.txt","r");
		if(!matrixFp || !vectorFp) {
			printf("Could not open %s.\n",matrixFp? "matriz.txt":"vetor.txt");
			kill(1);
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
		rewind(matrixFp);

		// Creating the matrix to be reduced.
		// I'm using c+1 because there's the results column which isn't counted when counting to c
		matrix = CreateMatrix(r, c + 1);
		
		for(i = 0; i < matrix->rows; i++){
			for(j = 0; j < matrix->cols - 1; j++){
				fscanf(matrixFp,"%lf", &(matrix->values[i][j]) );
			}
		}

		// VectorFp contains our matrix last column
		for(i = 0; i < matrix->rows; i++){
			fscanf(vectorFp,"%lf", &(matrix->values[i][c]) );
		}

		//sending matrix size
		MPI_Bcast(&r, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// Calculates chunk size
		chunkSize = (r / nproc) + 1; 

		fclose(matrixFp);
		fclose(vectorFp);
	} else {

		//receiving matrix size
		MPI_Bcast(&r, 1, MPI_INT, 0, MPI_COMM_WORLD);
		c = r+1;
	}

	// For each row in the matrix
	for(i = 0; i < matrix->rows; i++){
		
		/* Master only */
		if(rank == 0){

			#ifdef DEBUG
				printf("Searching pivot in col: %d\n", i);
				PrintMatrix(matrix);
			#endif

			/* Find pivot - Use OpenMP here */
			pcol = i;
			prow = FindPivot(matrix, pcol);

			#ifdef DEBUG
				if(prow == -1) printf("Pivot not found\n");
				else printf("Pivot line: %d\n", prow);
			#endif
			
			// Pivot not found (all values are 0 or matrix is reduced), 
			// skip to next column
			if(prow == -1) continue;

			/* Reduce pivot line (divide line by pivot) - Use OpenMP here */
			// Its easies to first divide pivot line and then swap it
			MultiplyLineByScalar(matrix, prow, 1.0/matrix->values[prow][i]);
				
			#ifdef DEBUG
				printf("Multplying matrix by 1/%d\n", matrix->values[prow][i]);
				PrintMatrix(matrix);
			#endif

			/* Position pivot */
			SwapLines(matrix, prow, i);
			prow = i;

			#ifdef DEBUG
				printf("Swapping lines %d and %d\n", prow, i);	
				PrintMatrix(matrix);
			#endif


			// /* Sends pivot and chunks to each slaves */
			MPI_Bcast(pivotRow, c, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Scatter(matrix->values, chunkSize * c, MPI_Double, currentRow, chunksize * c, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		// Slaves
		} else {
			/* Receive message (pivot and rank lines) */
			MPI_Bcast(pivotRow, c, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Scatter(matrix->values, chunkSize * c, MPI_Double, currentRow, chunksize * c, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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

		// Selecting the pivot row
		pivotRow = matrix->values[prow];
		// Alocating memory for processing each row
		double* backupRow = malloc(sizeof(double)*matrix->cols);

		// For each row in the matrix:
		int j;
		for(j = 0; j < matrix->rows; j++){

			// I don't need to process the pivot row
			if (j == prow) continue;

			/* we have to multiply the pivot row by a value 
			that would zero the value in this current row that is
			below or above the current pivot.*/

			// calculating the scalar that multiplies the selected row
			scalar = -matrix->values[j][pcol];

			// Creating a auxiliar vector to store the prod. value
			memcpy(backupRow, pivotRow, sizeof(double) * matrix->cols);
			// Multipying row by scalar
			_MultiplyLineByScalar(backupRow, matrix->cols, scalar);

			// Sum of currently selected row and multiplied row
			currentRow = matrix->values[j];

			_AddLines(currentRow,backupRow,matrix->cols);
		}
		
		/* END SEQUENTIAL CODE */

		/* PARALLEL CODE*/


		/* END PARALLEL CODE */

		/* Send result back to master */

		/* Master can debug print */
		if(rank == 0){

		}
	}


	if(backupRow) free(backupRow);
	MPI_Finalize();
	free(recv);
	free(sendVec);

	return 0;
}