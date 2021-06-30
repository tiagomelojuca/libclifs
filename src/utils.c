#include "utils.h"

#include <stdlib.h>

// --------------------------------------------------------------------------------

void _fillMatrixDefaultValue(double _matrix[SM][SM], double _initValue)
{
    for(int i = 0; i < SM; i++) {
        for(int j = 0; j < SM; j++) {
            _matrix[i][j] = _initValue;
        }
    }
}

// --------------------------------------------------------------------------------

void _fillMatrixNull(double _matrix[SM][SM])
{
    _fillMatrixDefaultValue(_matrix, 0.0);
}

// --------------------------------------------------------------------------------

void _fillDynIntMatrix(int** _matrix, int _nR, int _nC, int _initValue)
{
    for(int i = 0; i < _nR; i++) {
        for(int j = 0; j < _nC; j++) {
            _matrix[i][j] = _initValue;
        }
    }
}

// --------------------------------------------------------------------------------

void _fillDynDoubleMatrix(double** _matrix, int _nR, int _nC, double _initValue)
{
    for(int i = 0; i < _nR; i++) {
        for(int j = 0; j < _nC; j++) {
            _matrix[i][j] = _initValue;
        }
    }
}

// --------------------------------------------------------------------------------

void _copyDynDoubleMatrix(double** _matrix, double** _other, int _nR, int _nC)
{
    for(int i = 0; i < _nR; i++) {
        for(int j = 0; j < _nC; j++) {
            _matrix[i][j] = _other[i][j];
        }
    }
}

// --------------------------------------------------------------------------------

void _multiplyDynMatrix(double** _matrix, double** _other1, double** _other2,
                        int _rOther1, int _cOther1, int _rOther2, int _cOther2)
{   
    for(int i = 0; i < _rOther1; i++) {
        for(int j = 0; j < _cOther2; j++) {
            _matrix[i][j] = 0;
        }
    }
    
    for(int i = 0; i < _rOther1; i++) {
        for(int j = 0; j < _cOther2; j++) {
            for(int k = 0; k < _cOther1; k++) {
                _matrix[i][j] += _other1[i][k] * _other2[k][j];
            }
        }
    }
}

// --------------------------------------------------------------------------------

double _getAbsValue(double _val) {
    if(_val < 0) {
        return -_val;
    } else {
        return _val;
    }
}

// --------------------------------------------------------------------------------

void _lupDecompose(double **_matrix, int _matrixDimension, int *_vecPermut)
{
    const double tolerance = 0.00001;
    
    int maxIndex;
    double maxValueMatrix, absValueMatrix;
    double* pMatrix;

    for (int i = 0; i <= _matrixDimension; i++) {
        _vecPermut[i] = i;
    }

    for (int i = 0; i < _matrixDimension; i++) {
        maxValueMatrix = 0.0;
        maxIndex = i;

        for (int k = i; k < _matrixDimension; k++) {
            absValueMatrix = _getAbsValue(_matrix[k][i]);

            if (absValueMatrix > maxValueMatrix) { 
                maxValueMatrix = absValueMatrix;
                maxIndex = k;
            }
        }

        if (maxValueMatrix > tolerance) {
            if (maxIndex != i) {
                int swap = _vecPermut[i];
                _vecPermut[i] = _vecPermut[maxIndex];
                _vecPermut[maxIndex] = swap;

                pMatrix = _matrix[i];
                _matrix[i] = _matrix[maxIndex];
                _matrix[maxIndex] = pMatrix;

                _vecPermut[_matrixDimension]++;
            }

            for (int j = i + 1; j < _matrixDimension; j++) {
                _matrix[j][i] /= _matrix[i][i];

                for (int k = i + 1; k < _matrixDimension; k++) {
                    _matrix[j][k] -= _matrix[j][i] * _matrix[i][k];
                }
            }
        }
    }
}

// --------------------------------------------------------------------------------

void _lupInvert(double **_matrix, double **_other,
                int *_vecPermut, int _matrixDimension)
{
    for (int i = 0; i < _matrixDimension; i++) {
        for (int j = 0; j < _matrixDimension; j++) {

            if (_vecPermut[j] == i) {
                _matrix[j][i] = 1;
            } else {
                _matrix[j][i] = 0;
            }

            for (int k = 0; k < j; k++) {
                _matrix[j][i] -= _other[j][k] * _matrix[k][i];
            }
        }

        for (int j = _matrixDimension - 1; j >= 0; j--) {
            for (int k = j + 1; k < _matrixDimension; k++) {
                _matrix[j][i] -= _other[j][k] * _matrix[k][i];
            }

            _matrix[j][i] = _matrix[j][i] / _other[j][j];
        }
    }
}

// --------------------------------------------------------------------------------

double _calcDeterminant(double **_matrix, int *_vecPermut, int _matrixDimension)
{
    double det = _matrix[0][0];

    for (int i = 1; i < _matrixDimension; i++) {
        det *= _matrix[i][i];
    }

    if ((_vecPermut[_matrixDimension] - _matrixDimension) % 2 == 0) {
        return det;
    } else {
        return -det;
    }
}

// --------------------------------------------------------------------------------

bool _isMatrixSingular(double** _matrix, int _dimension)
{
    const int sizeVecPermut = _dimension + 1;
    int* vecPermut = malloc(sizeVecPermut * sizeof(int));

    double** mtxCopy = malloc(_dimension * sizeof(double*));
    for(int i = 0; i < _dimension; i++) {
        mtxCopy[i] = malloc(_dimension * sizeof(double));
    } _copyDynDoubleMatrix(mtxCopy, _matrix, _dimension, _dimension);

    _lupDecompose(mtxCopy, _dimension, vecPermut);
    double det = _calcDeterminant(mtxCopy, vecPermut, _dimension);

    return det == 0;
}

// --------------------------------------------------------------------------------
