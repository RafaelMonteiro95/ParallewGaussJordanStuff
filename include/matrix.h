#ifndef _MATRIX_H_
#define _MATRIX_H_

typedef struct matrix {

	int rows, cols;
	int **values;

} Matrix;

Matrix *CreateMatrix(int r, int c);
void DestroyMatrix(Matrix **m);
void PrintMatrix(Matrix *m);
void SwapLines(Matrix *matrix, int line1, int line2);

// Returns the pivot's line (-1 if not found or out of bounds)
int FindPivot(Matrix *matrix, int col);

// This function overwrite the matrix's destLine with the new values.
void AddLines(Matrix *matrix, int destLine, int line2);

// This function overwrite the matrix's destLine with the new values.
void MultiplyLineByScalar(Matrix *matrix, int line, float value);

#endif