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
	
	int *recv = NULL;
	int *sendVec;

	FILE* matrixFp;
	FILE* vectorFp;
	
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

		fclose(matrixFp);
		fclose(vectorFp);
	}

	int prow, pcol;
	// For each row
	for(i = 0; i < matrix->rows; i++){
		
		/* Master only */
		if(rank == 0){

			// #ifdef DEBUG
			// 	printf("Searching pivot in col: %d\n", i);
			// 	PrintMatrix(matrix);
			// #endif

			/* Find pivot - Use OpenMP here */
			pcol = i;
			prow = FindPivot(matrix, pcol);
			// printf("prow:%d\n",prow);

			// #ifdef DEBUG
			// 	if(prow == -1) printf("Pivot not found\n");
			// 	else printf("Pivot line: %d\n", prow);
			// #endif
			
			// Pivot not found (all values are 0 or matrix is reduced), 
			// skip to next column
			if(prow == -1) continue;

			/* Reduce pivot line (divide line by pivot) - Use OpenMP here */
			// Its easies to first divide pivot line and then swap it
			MultiplyLineByScalar(matrix, prow, 1.0/matrix->values[prow][i]);
				
			// #ifdef DEBUG
			// 	printf("Multplying matrix by 1/%d\n", matrix->values[prow][i]);
			// 	PrintMatrix(matrix);
			// #endif

			/* Position pivot */
			SwapLines(matrix, prow, i);

			// #ifdef DEBUG
			// 	printf("Swapping lines %d and %d\n", prow, i);	
			// 	PrintMatrix(matrix);
			// #endif

			/* Send pivot and their line (indexed by rank) to each slaves */

		// Slaves
		} else {
			/* Receive message (pivot and rank lines) */

		}
		

		/* SEQUENTIAL VERSION OF CODE*/

		double* pivotRow = matrix->values[prow];
		double* backupRow = malloc(sizeof(double)*matrix->cols);
		double* currentRow = NULL;

		int j;
		for(j = 0; j < matrix->rows; j++){

			if (j == prow) continue;
			// calculating the scalar value of line product
			double value = -matrix->values[j][pcol];

			// printf("matrix[%d][%d] = %lf ,value: %lf\n",j,pcol,matrix->values[j][pcol],value);

			// Creating a auxiliar vector to store the prod. value
			memcpy(backupRow, pivotRow, sizeof(double) * matrix->cols);
			// Multipying row by scalar
			_MultiplyLineByScalar(backupRow, matrix->cols, value);

			// Sum of currently selected row and multiplied row
			currentRow = matrix->values[j];
			_AddLines(currentRow,backupRow,matrix->cols);
		}
		printf("Final Result:\n");
		PrintMatrix(matrix);
		/* END SEQUENTIAL CODE */

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


		/* PARALLEL CODE*/


		/* END PARALLEL CODE */


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

	return 0;
}