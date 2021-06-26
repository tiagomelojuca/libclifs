#include "libclifs.h"

#include <math.h>

 //--------------------------------------------------------------------------------

void _setBarPropsLowLevel(Bar* _bar, Node _n1, Node _n2, Point _auxvec,
                          double _e, double _g,
                          double _a, double _iy, double _iz, double _j)
{
    _bar->node1 = _n1;
    _bar->node2 = _n2;
    _bar->auxvec = _auxvec;
    _bar->e = _e;
    _bar->g = _g;
    _bar->a = _a;
    _bar->iy = _iy;
    _bar->iz = _iz;
    _bar->j = _j;
}

 //--------------------------------------------------------------------------------

void _setBarPropsHighLevel(Bar* _bar, Node _n1, Node _n2, Point _auxvec,
                           Material* _material, Section* _section)
{
    _setBarPropsLowLevel(_bar, _n1, _n2, _auxvec,
                         _material->e, _material->g,
                         _section->a, _section->iy, _section->iz, _section->j);
}

 //--------------------------------------------------------------------------------

void _filllocalStiffnessMatrix(double _matrix[12][12], Bar* _bar)
{
    double L  = _bar->l;
    double E  = _bar->e;
    double G  = _bar->g;
    double A  = _bar->a;
    double Iy = _bar->iy;
    double Iz = _bar->iz;
    double J  = _bar->j;
        
    _matrix[0][0]   = (E * A) / L;
    _matrix[6][0]   = (-1) * ((E * A) / L);
    _matrix[1][1]   = (12 * E * Iz) / pow(L, 3);
    _matrix[5][1]   = (6 * E * Iz) / pow(L, 2);
    _matrix[7][1]   = (-1) * ((12 * E * Iz) / pow(L, 3));
    _matrix[11][1]  = (6 * E * Iz) / pow(L, 2);
    _matrix[2][2]   = (12 * E * Iy) / pow(L, 3);
    _matrix[4][2]   = (-1) * ((6 * E * Iy) / pow(L, 2));
    _matrix[8][2]   = (-1) * ((12 * E * Iy) / pow(L, 3));
    _matrix[10][2]  = (-1) * ((6 * E * Iy) / pow(L, 2));
    _matrix[3][3]   = (G * J) / L;
    _matrix[9][3]   = (-1) * ((G * J) / L);
    _matrix[2][4]   = (-1) * ((6 * E * Iy) / pow(L, 2));
    _matrix[4][4]   = (4 * E * Iy) / L;
    _matrix[8][4]   = (6 * E * Iy) / pow(L, 2);
    _matrix[10][4]  = (2 * E * Iy) / L;
    _matrix[1][5]   = (6 * E * Iz) / pow(L, 2);
    _matrix[5][5]   = (4 * E * Iz) / L;
    _matrix[7][5]   = (-1) * ((6 * E * Iz) / pow(L, 2));
    _matrix[11][5]  = (2 * E * Iz) / L;
    _matrix[0][6]   = (-1) * ((E * A) / L);
    _matrix[6][6]   = (E * A) / L;
    _matrix[1][7]   = (-1) * ((12 * E * Iz) / pow(L, 3));
    _matrix[5][7]   = (-1) * ((6 * E * Iz) / pow(L, 2));
    _matrix[7][7]   = (12 * E * Iz) / pow(L, 3);
    _matrix[11][7]  = (-1) * ((6 * E * Iz) / pow(L, 2));
    _matrix[2][8]   = (-1) * ((12 * E * Iy) / pow(L, 3));
    _matrix[4][8]   = (6 * E * Iy) / pow(L, 2);
    _matrix[8][8]   = (12 * E * Iy) / pow(L, 3);
    _matrix[10][8]  = (6 * E * Iy) / pow(L, 2);
    _matrix[3][9]   = (-1) * ((G * J) / L);
    _matrix[9][9]   = (G * J) / L;
    _matrix[2][10]  = (-1) * ((6 * E * Iy) / pow(L, 2));
    _matrix[4][10]  = (2 * E * Iy) / L;
    _matrix[8][10]  = (6 * E * Iy) / pow(L, 2);
    _matrix[10][10] = (4 * E * Iy) / L;
    _matrix[1][11]  = (6 * E * Iz) / pow(L, 2);
    _matrix[5][11]  = (2 * E * Iz) / L;
    _matrix[7][11]  = (-1) * ((6 * E * Iz) / pow(L, 2));
    _matrix[11][11] = (4 * E * Iz) / L;
}

 //--------------------------------------------------------------------------------

double _calcDistBetweenPoints(Point p1, Point p2)
{
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;
    
    return sqrt(
        pow(dx, 2) + pow(dy, 2) + pow(dz, 2)
    );
}

 //--------------------------------------------------------------------------------
