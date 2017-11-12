#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "matrix.h"

#define PARALLEL_THRESHOLD 40

Matrix *CreateMatrix(int r, int c){

	Matrix *matrix = (Matrix *) malloc(sizeof(Matrix));
	
	if(!matrix) return NULL;

	matrix->values = (double *) malloc(sizeof(double)*(c*r));
	if(!matrix->values) return NULL;

	matrix->rows = r;
	matrix->cols = c;

	return matrix;
}

void DestroyMatrix(Matrix **matrix){

	if(!matrix) return;
	free((*matrix)->values);
	free(*matrix);
}

void PrintMatrix(Matrix *matrix){ FPrintMatrix(matrix, stdout); }
void FPrintMatrix(Matrix *matrix, FILE* fp){

	if(!matrix) return (void) fprintf(fp, "(nil)");
	
	for(int i = 0; i < matrix->rows; i++){
		for(int j = 0; j < matrix->cols; j++){
			// Sometimes -0.0 is stored instead of 0.0.
			if(fabs(matrix->values[mat2vec(matrix->cols, i, j)]) < 0.000001) matrix->values[mat2vec(matrix->cols, i, j)] = 0.0;

			fprintf(fp, "%-5.2lf ", matrix->values[mat2vec(matrix->cols, i, j)]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");
}

int FindPivot(Matrix *matrix, int col){
	
	int pivotLine = -1;

	#pragma omp parallel for if(matrix->rows - col > PARALLEL_THRESHOLD)
	for(int i = col; i < matrix->rows; i++){
		if(matrix->values[mat2vec(matrix->cols, i, col)] != 0){
			pivotLine = i;
		}
	}

	return pivotLine;
}

void SwapLines(Matrix *matrix, int line1, int line2){
	
	double *aux = (double *) malloc(sizeof(double)*matrix->cols); 
	double *line1p = &(matrix->values[matrix->cols*line1]);
	double *line2p = &(matrix->values[matrix->cols*line2]);
	
	memcpy(aux   , line1p, matrix->cols*sizeof(double));
	memcpy(line1p, line2p, matrix->cols*sizeof(double));
	memcpy(line2p, aux   , matrix->cols*sizeof(double));
	
	free(aux);
}

void MultiplyLineByScalar(Matrix *matrix, int line, double value){

	#pragma omp parallel for
	for(int i = 0; i < matrix->cols; i++)
		matrix->values[mat2vec(matrix->cols, line, i)] *= value;
}

void _MultiplyLineByScalar(double* line, int size, double value){

	#pragma omp parallel for
	for(int i = 0; i < size; i++)
		line[i] *= value;
}

void AddLines(Matrix *matrix, int line1, int line2){
	
	#pragma omp parallel for
	for(int i = 0; i < matrix->cols; i++)
		matrix->values[mat2vec(matrix->cols, line1, i)] += 
		matrix->values[mat2vec(matrix->cols, line2, i)];
}

void _AddLines(double* destline, double* oline, int size){

	#pragma omp parallel for
	for(int i = 0; i < size; i++)
		destline[i] += oline[i];
}
