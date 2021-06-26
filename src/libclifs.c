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

void _fillNullMatrix(double _matrix[12][12])
{
    for(int i = 0; i < 12; i++) {
        for(int j = 0; j < 12; j++) {
            _matrix[i][j] = 0;
        }
    }
}

//--------------------------------------------------------------------------------

void _fillLocalStiffnessMatrix(double _matrix[12][12], Bar* _bar)
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

void _fillReducedRotationMatrix(double _matrix[3][3], Bar* _bar)
{
    // This is all civil engineering business logic, doesn't bother at all
    _matrix[0][0] = (_bar->node2.position.x - _bar->node1.position.x) / _bar->l;
    _matrix[0][1] = (_bar->node2.position.y - _bar->node1.position.y) / _bar->l;
    _matrix[0][2] = (_bar->node2.position.z - _bar->node1.position.z) / _bar->l;

    double dx = _bar->auxvec.x - _bar->node1.position.x;
    double dy = _bar->auxvec.y - _bar->node1.position.y;
    double dz = _bar->auxvec.z - _bar->node1.position.z;
    double auxvecLength = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));

    double cosAlfa = (_bar->auxvec.x - _bar->node1.position.x) / auxvecLength;
    double cosBeta = (_bar->auxvec.y - _bar->node1.position.y) / auxvecLength;
    double cosGama = (_bar->auxvec.z - _bar->node1.position.z) / auxvecLength;

    double c = sqrt(pow((_matrix[0][1] * cosGama - _matrix[0][2] * cosBeta), 2) +
                    pow((_matrix[0][2] * cosAlfa - _matrix[0][0] * cosGama), 2) +
                    pow((_matrix[0][0] * cosBeta - _matrix[0][1] * cosAlfa), 2));

    _matrix[2][0] = (_matrix[0][1] * cosGama - _matrix[0][2] * cosBeta) / c;
    _matrix[2][1] = (_matrix[0][2] * cosAlfa - _matrix[0][0] * cosGama) / c;
    _matrix[2][2] = (_matrix[0][0] * cosBeta - _matrix[0][1] * cosAlfa) / c;

    _matrix[1][0] = _matrix[0][2] * _matrix[2][1] - _matrix[0][1] * _matrix[2][2];
    _matrix[1][1] = _matrix[0][0] * _matrix[2][2] - _matrix[0][2] * _matrix[2][0];
    _matrix[1][2] = _matrix[0][1] * _matrix[2][0] - _matrix[0][0] * _matrix[2][1];
}

//--------------------------------------------------------------------------------

double _calcDistBetweenPoints(Point _p1, Point _p2)
{
    double dx = _p1.x - _p2.x;
    double dy = _p1.y - _p2.y;
    double dz = _p1.z - _p2.z;
    
    return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
}

//--------------------------------------------------------------------------------
