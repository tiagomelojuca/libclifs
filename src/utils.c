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
