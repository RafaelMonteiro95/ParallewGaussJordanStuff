#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

#define PARALLEL_THRESHOLD 40

Matrix *CreateMatrix(int r, int c){

	Matrix *m = (Matrix *) malloc(sizeof(Matrix));
	
	if(!m) return NULL;

	m->values = (int **) malloc(sizeof(int *)*r);
	if(!m->values) return NULL;

	m->rows = r;
	m->cols = c;

	for(int i = 0; i < r; i++){
		m->values[i] = (int *) malloc(sizeof(int)*c);
		if(!m->values[i]) return NULL;
	}

	return m;
}

void DestroyMatrix(Matrix **m){

	if(!m) return;

	for(int i = 0; i < (*m)->rows; i++)
		free((*m)->values[i]);
	free((*m)->values);
	free(*m);
}

void PrintMatrix(Matrix *m){

	if(!m) return (void) printf("(nil)");
	
	for(int i = 0; i < m->rows; i++){
		for(int j = 0; j < m->cols; j++){
			printf("%d ", m->values[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

int FindPivot(Matrix *matrix, int col){
	
	int pivotLine = -1;

	#pragma omp parallel for if(matrix->rows - col > PARALLEL_THRESHOLD)
	for(int i = col; i < matrix->rows; i++){
		if(matrix->values[i][col] != 0){
			pivotLine = i;
		}
	}

	return pivotLine;
}

void SwapLines(Matrix *matrix, int line1, int line2){
	int *aux = matrix->values[line1];
	matrix->values[line1] = matrix->values[line2];
	matrix->values[line2] = aux;
}

void MultiplyLineByScalar(Matrix *matrix, int line, float value){

	#pragma omp parallel for
	for(int i = 0; i < matrix->cols; i++)
		matrix->values[line][i] *= value;
}

void AddLines(Matrix *matrix, int line1, int line2){
	
	#pragma omp parallel for
	for(int i = 0; i < matrix->cols; i++)
		matrix->values[line1][i] += matrix->values[line2][i];
}