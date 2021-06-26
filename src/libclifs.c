#include "libclifs.h"

#include <math.h>

// --------------------------------------------------------------------------------

void _filllocalStiffnessMatrix(double _matrix[12][12], Bar* _bar)
{
    double l  = _bar->l;
    double e  = _bar->e;
    double g  = _bar->g;
    double a  = _bar->a;
    double iy = _bar->iy;
    double iz = _bar->iz;
    double j  = _bar->j;

    for(int i = 0; i < 12; i++)
    {
        for(int j = 0; j < 12; j++)
        {
            _matrix[i][j] = 0;
        }
    }
        
    _matrix[0][0]   = (e*a)/l;
    _matrix[6][0]   = (-1)*((e*a)/l);
    _matrix[1][1]   = (12*e*iz)/pow(l,3);
    _matrix[5][1]   = (6*e*iz)/pow(l,2);
    _matrix[7][1]   = (-1)*((12*e*iz)/pow(l,3));
    _matrix[11][1]  = (6*e*iz)/pow(l,2);
    _matrix[2][2]   = (12*e*iy)/pow(l,3);
    _matrix[4][2]   = (-1)*((6*e*iy)/pow(l,2));
    _matrix[8][2]   = (-1)*((12*e*iy)/pow(l,3));
    _matrix[10][2]  = (-1)*((6*e*iy)/pow(l,2));
    _matrix[3][3]   = (g*j)/l;
    _matrix[9][3]   = (-1)*((g*j)/l);
    _matrix[2][4]   = (-1)*((6*e*iy)/pow(l,2));
    _matrix[4][4]   = (4*e*iy)/l;
    _matrix[8][4]   = (6*e*iy)/pow(l,2);
    _matrix[10][4]  = (2*e*iy)/l;
    _matrix[1][5]   = (6*e*iz)/pow(l,2);
    _matrix[5][5]   = (4*e*iz)/l;
    _matrix[7][5]   = (-1)*((6*e*iz)/pow(l,2));
    _matrix[11][5]  = (2*e*iz)/l;
    _matrix[0][6]   = (-1)*((e*a)/l);
    _matrix[6][6]   = (e*a)/l;
    _matrix[1][7]   = (-1)*((12*e*iz)/pow(l,3));
    _matrix[5][7]   = (-1)*((6*e*iz)/pow(l,2));
    _matrix[7][7]   = (12*e*iz)/pow(l,3);
    _matrix[11][7]  = (-1)*((6*e*iz)/pow(l,2));
    _matrix[2][8]   = (-1)*((12*e*iy)/pow(l,3));
    _matrix[4][8]   = (6*e*iy)/pow(l,2);
    _matrix[8][8]   = (12*e*iy)/pow(l,3);
    _matrix[10][8]  = (6*e*iy)/pow(l,2);
    _matrix[3][9]   = (-1)*((g*j)/l);
    _matrix[9][9]   = (g*j)/l;
    _matrix[2][10]  = (-1)*((6*e*iy)/pow(l,2));
    _matrix[4][10]  = (2*e*iy)/l;
    _matrix[8][10]  = (6*e*iy)/pow(l,2);
    _matrix[10][10] = (4*e*iy)/l;
    _matrix[1][11]  = (6*e*iz)/pow(l,2);
    _matrix[5][11]  = (2*e*iz)/l;
    _matrix[7][11]  = (-1)*((6*e*iz)/pow(l,2));
    _matrix[11][11] = (4*e*iz)/l;
}

// --------------------------------------------------------------------------------

double _calcDistBetweenPoints(Point p1, Point p2)
{
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;
    
    return sqrt(
        pow(dx, 2) + pow(dy, 2) + pow(dz, 2)
    );
}

// --------------------------------------------------------------------------------
