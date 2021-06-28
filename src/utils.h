#ifndef UTILS_H
#define UTILS_H

#include "definitions.h"

void _fillMatrixDefaultValue(double _matrix[SM][SM], double _initValue);
void _fillMatrixNull(double _matrix[SM][SM]);
void _fillDynIntMatrix(int** _matrix, int _nR, int _nC, int _initValue);
void _fillDynDoubleMatrix(double** _matrix, int _nR, int _nC, double _initValue);

#endif // UTILS_H