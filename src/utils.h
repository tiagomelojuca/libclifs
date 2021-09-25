// --------------------------------------------------------------------------------
// LIBCLIFS - TIAGO MELO JUCA (2021)                                              |
// --------------------------------------------------------------------------------
// This library is for studying purposes only - please don't take it seriously :)
// You'll see a lot of comments here. Yep, I know Clean Code, but, as I said, it's
// just for learning; the purpose is a reference for myself, nothing more
// --------------------------------------------------------------------------------

#ifndef UTILS_H
#define UTILS_H

#include <stdbool.h>

// --------------------------------------------------------------------------------
#define MTX_SIZE 12
// --------------------------------------------------------------------------------

void _fillMatrixDefaultValue(double _matrix[MTX_SIZE][MTX_SIZE], double _initVal);
void _fillMatrixNull(double _matrix[MTX_SIZE][MTX_SIZE]);
void _fillDynIntMatrix(int** _matrix, int _nR, int _nC, int _initVal);
void _fillDynDoubleMatrix(double** _matrix, int _nR, int _nC, double _initVal);
void _copyDynDoubleMatrix(double** _matrix, double** _other, int _nR, int _nC);
void _multiplyDynMatrix(double** _matrix, double** _other1, double** _other2,
                        int _rOther1, int _cOther1, int _rOther2, int _cOther2);

double _getAbsValue(double _val);
void _lupDecompose(double **_matrix, int _matrixDimension, int *_vecPermut);
void _lupInvert(double **_matrix, double **_other,
                int *_vecPermut, int _matrixDimension);
double _calcDeterminant(double **_matrix, int *_vecPermut, int _matrixDimension);

bool _isMatrixSingular(double** _matrix, int _dimension);

#endif // UTILS_H
