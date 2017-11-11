#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "matrix.h"

#define PARALLEL_THRESHOLD 40

Matrix *CreateMatrix(int r, int c){

	Matrix *m = (Matrix *) malloc(sizeof(Matrix));
	
	if(!m) return NULL;

	m->values = (double **) malloc(sizeof(double *)*r);
	if(!m->values) return NULL;

	m->rows = r;
	m->cols = c;

	for(int i = 0; i < r; i++){
		m->values[i] = (double *) malloc(sizeof(double)*c);
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

	FPrintMatrix(m,stdout);
}

void FPrintMatrix(Matrix *m, FILE* fp){

	if(!m) return (void) fprintf(fp,"(nil)");
	
	for(int i = 0; i < m->rows; i++){
		for(int j = 0; j < m->cols; j++){
			// Sometimes -0.0 is stored instead of 0.0.
			if(fabs(m->values[i][j]) < 0.000001) m->values[i][j] = 0.0;
			
			fprintf(fp,"%-3.2lf ", m->values[i][j]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
}

int FindPivot(Matrix *matrix, int col){
	
	int pivotLine = -1;

	// #pragma omp parallel for if(matrix->rows - col > PARALLEL_THRESHOLD)
	for(int i = col; i < matrix->rows; i++){
		if(matrix->values[i][col] != 0){
			pivotLine = i;
			break;
		}
	}

	return pivotLine;
}

void SwapLines(Matrix *matrix, int line1, int line2){
	double *aux = matrix->values[line1];
	matrix->values[line1] = matrix->values[line2];
	matrix->values[line2] = aux;
}

void MultiplyLineByScalar(Matrix *matrix, int line, double value){

	#pragma omp parallel for
	for(int i = 0; i < matrix->cols; i++)
		matrix->values[line][i] *= value;
}

void _MultiplyLineByScalar(double* line, int size, double value){

	#pragma omp parallel for
	for(int i = 0; i < size; i++)
		line[i] *= value;
}

void AddLines(Matrix *matrix, int line1, int line2){
	
	#pragma omp parallel for
	for(int i = 0; i < matrix->cols; i++)
		matrix->values[line1][i] += matrix->values[line2][i];
}

void _AddLines(double* destline, double* oline, int size){

	#pragma omp parallel for
	for(int i = 0; i < size; i++)
		destline[i] += oline[i];
}

double *ToArray(Matrix *matrix, int *length){

	int _length = matrix->rows*matrix->cols; // Resulting array length
	size_t lineSize = sizeof(double)*matrix->cols; // Matrix's lines size in bytes
	
	double *array = (double *) malloc(sizeof(double)*_length);

	for(int i = 0; i < matrix->rows; i++)
		memcpy(&array[i*matrix->cols], matrix->values[i], lineSize);

	if(length != NULL) *length = _length;
	return array;
}
