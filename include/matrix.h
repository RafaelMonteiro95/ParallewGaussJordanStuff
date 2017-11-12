#ifndef _MATRIX_H_
#define _MATRIX_H_

typedef struct matrix {

	int rows, cols;
	double **values;

} Matrix;

Matrix *CreateMatrix(int r, int c);
void DestroyMatrix(Matrix **m);
void PrintMatrix(Matrix *m);
void FPrintMatrix(Matrix *m, FILE* fp);
void SwapLines(Matrix *matrix, int line1, int line2);

// Returns the pivot's line (-1 if not found or out of bounds)
int FindPivot(Matrix *matrix, int col, int prevPivotLine);

// This function overwrite the matrix's destLine with the new values.
void AddLines(Matrix *matrix, int destLine, int line2);
void _AddLines(double* destline, double* oline, int size);

// This function overwrite the matrix's destLine with the new values.
void MultiplyLineByScalar(Matrix *matrix, int line, double value);
void _MultiplyLineByScalar(double *line, int size, double value);

// The user must free the returned buffer
double *ToArray(Matrix *matrix, int *length);

#endif