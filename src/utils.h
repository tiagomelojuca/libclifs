#ifndef UTILS_H
#define UTILS_H

#include "definitions.h"

void _fillMatrixDefaultValue(double _matrix[SM][SM], double _initValue);
void _fillMatrixNull(double _matrix[SM][SM]);
void _fillDynIntMatrix(int** _matrix, int _nR, int _nC, int _initValue);
void _fillDynDoubleMatrix(double** _matrix, int _nR, int _nC, double _initValue);
void _copyDynDoubleMatrix(double** _matrix, double** _other, int _nR, int _nC);
void _multiplyDynMatrix(double** _matrix, double** _other1, double** _other2,
                        int _rOther1, int _cOther1, int _rOther2, int _cOther2);

#endif // UTILS_H