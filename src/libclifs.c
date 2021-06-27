#include "libclifs.h"

#include <math.h>

// --------------------------------------------------------------------------------

void setPointCoords(Point* _point, double _x, double _y, double _z)
{
    _point->x = _x;
    _point->y = _y;
    _point->z = _z;
}

// --------------------------------------------------------------------------------

Point createPoint(double _x, double _y, double _z)
{
    Point p;
    setPointCoords(&p, _x, _y, _z);

    return p;
}

// --------------------------------------------------------------------------------

double _calcDistBetweenPoints(Point _p1, Point _p2)
{
    double dx = _p1.x - _p2.x;
    double dy = _p1.y - _p2.y;
    double dz = _p1.z - _p2.z;
    
    return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
}

// --------------------------------------------------------------------------------

void setDegreesOfFreedomProps(DegreesOfFreedom* _dof, bool _x, bool _y, bool _z)
{
    _dof->x = _x;
    _dof->y = _y;
    _dof->z = _z;
}

// --------------------------------------------------------------------------------

DegreesOfFreedom createDegreesOfFreedom(bool _x, bool _y, bool _z)
{
    DegreesOfFreedom dof;
    setDegreesOfFreedomProps(&dof, _x, _y, _z);

    return dof;
}

// --------------------------------------------------------------------------------

void setNodalLoadValues(NodalLoad* _load, double _fx, double _fy, double _fz,
                        double _mx, double _my, double _mz)
{
    _load->fx = _fx;
    _load->fy = _fy;
    _load->fz = _fz;
    _load->mx = _mx;
    _load->my = _my;
    _load->mz = _mz;
}

// --------------------------------------------------------------------------------

NodalLoad createNodalLoad(double _fx, double _fy, double _fz,
                          double _mx, double _my, double _mz)
{
    NodalLoad l;
    setNodalLoadValues(&l, _fx, _fy, _fz, _mx, _my, _mz);

    return l;
}

// --------------------------------------------------------------------------------

void setNodeProps(Node* _node, Point _position,
                  DegreesOfFreedom _translation, DegreesOfFreedom _rotation,
                  NodalLoad _load)
{
    _node->position = _position;
    _node->translation = _translation;
    _node->rotation = _rotation;
    _node->load = _load;
}

// --------------------------------------------------------------------------------

Node createNode(Point _position,
                DegreesOfFreedom _translation, DegreesOfFreedom _rotation,
                NodalLoad _load)
{
    Node n;
    setNodeProps(&n, _position, _translation, _rotation, _load);

    return n;
}

// --------------------------------------------------------------------------------

void _setBarPropsLowLevel(Bar* _bar, Node _n1, Node _n2, Point _auxvec,
                          double _e, double _g,
                          double _a, double _iy, double _iz, double _j)
{
    _bar->node1 = _n1;
    _bar->node2 = _n2;
    _bar->auxvec = _auxvec;
    _bar->l = _calcDistBetweenPoints(_n1.position, _n2.position);
    _bar->e = _e;
    _bar->g = _g;
    _bar->a = _a;
    _bar->iy = _iy;
    _bar->iz = _iz;
    _bar->j = _j;
}

// --------------------------------------------------------------------------------

void setBarProps(Bar* _bar, Node _n1, Node _n2, Point _auxvec,
                 Material* _material, Section* _section)
{
    _setBarPropsLowLevel(_bar, _n1, _n2, _auxvec,
                         _material->e, _material->g,
                         _section->a, _section->iy, _section->iz, _section->j);
}

// --------------------------------------------------------------------------------

void setMaterialProps(Material* _material, double _e, double _g)
{
    _material->e = _e;
    _material->g = _g;
}

// --------------------------------------------------------------------------------

Material createMaterial(double _e, double _g)
{
    Material m;
    setMaterialProps(&m, _e, _g);

    return m;
}

// --------------------------------------------------------------------------------

void setSectionProps(Section* _section, double _a,
                     double _iy, double _iz, double _j)
{
    _section->a = _a;
    _section->iy = _iy;
    _section->iz = _iz;
    _section->j = _j;
}

// --------------------------------------------------------------------------------

Section createSection(double _a, double _iy, double _iz, double _j)
{
    Section s;
    setSectionProps(&s, _a, _iy, _iz, _j);

    return s;
}

// --------------------------------------------------------------------------------

void _fillMatrixDefaultValue(double _matrix[12][12], double _defaultValue)
{
    for(int i = 0; i < 12; i++) {
        for(int j = 0; j < 12; j++) {
            _matrix[i][j] = _defaultValue;
        }
    }
}

// --------------------------------------------------------------------------------

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

// --------------------------------------------------------------------------------

void _fillReducedRotationMatrix(double _matrix[3][3], Bar* _bar)
{
    // This is all civil engineering business logic, don't bother at all
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

// --------------------------------------------------------------------------------

void _fillRotationMatrix(double _matrix[12][12], double _reducedMatrix[3][3])
{
    // Civil engineering business logic, don't bother at all
    _matrix[0][0]   = _reducedMatrix[0][0];
    _matrix[3][3]   = _reducedMatrix[0][0];
    _matrix[6][6]   = _reducedMatrix[0][0];
    _matrix[9][9]   = _reducedMatrix[0][0];

    _matrix[0][1]   = _reducedMatrix[0][1];
    _matrix[3][4]   = _reducedMatrix[0][1];
    _matrix[6][7]   = _reducedMatrix[0][1];
    _matrix[9][10]  = _reducedMatrix[0][1];

    _matrix[0][2]   = _reducedMatrix[0][2];
    _matrix[3][5]   = _reducedMatrix[0][2];
    _matrix[6][8]   = _reducedMatrix[0][2];
    _matrix[9][11]  = _reducedMatrix[0][2];

    _matrix[1][0]   = _reducedMatrix[1][0];
    _matrix[4][3]   = _reducedMatrix[1][0];
    _matrix[7][6]   = _reducedMatrix[1][0];
    _matrix[10][9]  = _reducedMatrix[1][0];

    _matrix[1][1]   = _reducedMatrix[1][1];
    _matrix[4][4]   = _reducedMatrix[1][1];
    _matrix[7][7]   = _reducedMatrix[1][1];
    _matrix[10][10] = _reducedMatrix[1][1];

    _matrix[1][2]   = _reducedMatrix[1][2];
    _matrix[4][5]   = _reducedMatrix[1][2];
    _matrix[7][8]   = _reducedMatrix[1][2];
    _matrix[10][11] = _reducedMatrix[1][2];

    _matrix[2][0]   = _reducedMatrix[2][0];
    _matrix[5][3]   = _reducedMatrix[2][0];
    _matrix[8][6]   = _reducedMatrix[2][0];
    _matrix[11][9]  = _reducedMatrix[2][0];

    _matrix[2][1]   = _reducedMatrix[2][1];
    _matrix[5][4]   = _reducedMatrix[2][1];
    _matrix[8][7]   = _reducedMatrix[2][1];
    _matrix[11][10] = _reducedMatrix[2][1];

    _matrix[2][2]   = _reducedMatrix[2][2];
    _matrix[5][5]   = _reducedMatrix[2][2];
    _matrix[8][8]   = _reducedMatrix[2][2];
    _matrix[11][11] = _reducedMatrix[2][2];
}

// --------------------------------------------------------------------------------

void _fillTransposeRotationMatrix(double _matrix[12][12], double _other[12][12])
{
    for(int i = 0; i < 12; i++) {
        for(int j = 0; j < 12; j++) {
            _matrix[i][j] = _other[j][i];
        }
    }
}

// --------------------------------------------------------------------------------

void _fillGlobalStiffnessMatrix(double _matrix[12][12],
                                double _tRotation[12][12], double _local[12][12])
{
    for(int i = 0; i < 12; i++) {
        for(int j = 0; j < 12; j++) {
            for(int k = 0; k < 12; k++) {
                _matrix[i][j] += _tRotation[i][k] * _local[k][j];
            }
        }
    }
}

// --------------------------------------------------------------------------------

void setStiffnessMatrix(StiffnessMatrix* _sMatrix, Bar* _associatedBar)
{
    _fillMatrixDefaultValue(_sMatrix->local, 1.0);
    _fillMatrixDefaultValue(_sMatrix->rotation, 3.0);
    _fillMatrixDefaultValue(_sMatrix->transposeRotation, 4.0);
    _fillMatrixDefaultValue(_sMatrix->global, 5.0);

    _fillLocalStiffnessMatrix(_sMatrix->local, _associatedBar);
    _fillReducedRotationMatrix(_sMatrix->reducedRotation, _associatedBar);
    _fillRotationMatrix(_sMatrix->rotation, _sMatrix->reducedRotation);
    _fillTransposeRotationMatrix(_sMatrix->transposeRotation, _sMatrix->rotation);
    _fillGlobalStiffnessMatrix(_sMatrix->global,
                               _sMatrix->transposeRotation,
                               _sMatrix->local);
}

// --------------------------------------------------------------------------------

void setFrameBarProps(FrameBar* _frameBar, Node _n1, Node _n2, Point _auxvec,
                      Material* _material, Section* _section)
{
    Bar* pBar = &(_frameBar->bar);
    setBarProps(pBar, _n1, _n2, _auxvec, _material, _section);

    StiffnessMatrix* pStiffnessMatrix = &(_frameBar->stiffnessMatrix);
    setStiffnessMatrix(pStiffnessMatrix, pBar);
}

// --------------------------------------------------------------------------------

FrameBar createFrameBar(Node _n1, Node _n2, Point _auxvec,
                        Material* _material, Section* _section)
{
    FrameBar fb;
    setFrameBarProps(&fb, _n1, _n2, _auxvec, _material, _section);

    return fb;
}

// --------------------------------------------------------------------------------
