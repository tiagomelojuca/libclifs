#include "utils.h"

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
