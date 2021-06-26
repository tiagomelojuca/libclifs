#include "libclifs.h"

#include <math.h>

void fillLocalStiffnessMatrix(double _matrix[12][12], Bar* _bar)
{
    double L  = _bar->L;
    double E  = _bar->E;
    double G  = _bar->G;
    double A  = _bar->A;
    double Iy = _bar->Iy;
    double Iz = _bar->Iz;
    double J  = _bar->J;

    for(int i = 0; i < 12; i++)
    {
        for(int j = 0; j < 12; j++)
        {
            _matrix[i][j] = 0;
        }
    }
        
    _matrix[0][0]   = (E*A)/L;
    _matrix[6][0]   = (-1)*((E*A)/L);
    _matrix[1][1]   = (12*E*Iz)/pow(L,3);
    _matrix[5][1]   = (6*E*Iz)/pow(L,2);
    _matrix[7][1]   = (-1)*((12*E*Iz)/pow(L,3));
    _matrix[11][1]  = (6*E*Iz)/pow(L,2);
    _matrix[2][2]   = (12*E*Iy)/pow(L,3);
    _matrix[4][2]   = (-1)*((6*E*Iy)/pow(L,2));
    _matrix[8][2]   = (-1)*((12*E*Iy)/pow(L,3));
    _matrix[10][2]  = (-1)*((6*E*Iy)/pow(L,2));
    _matrix[3][3]   = (G*J)/L;
    _matrix[9][3]   = (-1)*((G*J)/L);
    _matrix[2][4]   = (-1)*((6*E*Iy)/pow(L,2));
    _matrix[4][4]   = (4*E*Iy)/L;
    _matrix[8][4]   = (6*E*Iy)/pow(L,2);
    _matrix[10][4]  = (2*E*Iy)/L;
    _matrix[1][5]   = (6*E*Iz)/pow(L,2);
    _matrix[5][5]   = (4*E*Iz)/L;
    _matrix[7][5]   = (-1)*((6*E*Iz)/pow(L,2));
    _matrix[11][5]  = (2*E*Iz)/L;
    _matrix[0][6]   = (-1)*((E*A)/L);
    _matrix[6][6]   = (E*A)/L;
    _matrix[1][7]   = (-1)*((12*E*Iz)/pow(L,3));
    _matrix[5][7]   = (-1)*((6*E*Iz)/pow(L,2));
    _matrix[7][7]   = (12*E*Iz)/pow(L,3);
    _matrix[11][7]  = (-1)*((6*E*Iz)/pow(L,2));
    _matrix[2][8]   = (-1)*((12*E*Iy)/pow(L,3));
    _matrix[4][8]   = (6*E*Iy)/pow(L,2);
    _matrix[8][8]   = (12*E*Iy)/pow(L,3);
    _matrix[10][8]  = (6*E*Iy)/pow(L,2);
    _matrix[3][9]   = (-1)*((G*J)/L);
    _matrix[9][9]   = (G*J)/L;
    _matrix[2][10]  = (-1)*((6*E*Iy)/pow(L,2));
    _matrix[4][10]  = (2*E*Iy)/L;
    _matrix[8][10]  = (6*E*Iy)/pow(L,2);
    _matrix[10][10] = (4*E*Iy)/L;
    _matrix[1][11]  = (6*E*Iz)/pow(L,2);
    _matrix[5][11]  = (2*E*Iz)/L;
    _matrix[7][11]  = (-1)*((6*E*Iz)/pow(L,2));
    _matrix[11][11] = (4*E*Iz)/L;
}