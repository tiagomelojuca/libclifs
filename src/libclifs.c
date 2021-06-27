#include "libclifs.h"

#include <stdlib.h>
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

    _dof->constraints = 0;

    // Business Logic: note that if any _axis==true, it means freedom in that axis
    // Thus, we have a constraint if _axis==false (no freedom in axis)
    // So, we increment constraints only if !_axis
    if(!_x) { _dof->constraints++; }
    if(!_y) { _dof->constraints++; }
    if(!_z) { _dof->constraints++; }
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
    _node->constraints = _translation.constraints + _rotation.constraints;
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

void _setBarPropsLowLevel(Bar* _bar, Node* _n1, Node* _n2, Point _auxvec,
                          double _e, double _g,
                          double _a, double _iy, double _iz, double _j)
{
    _bar->node1 = _n1;
    _bar->node2 = _n2;
    _bar->auxvec = _auxvec;
    _bar->l = _calcDistBetweenPoints(_n1->position, _n2->position);
    _bar->e = _e;
    _bar->g = _g;
    _bar->a = _a;
    _bar->iy = _iy;
    _bar->iz = _iz;
    _bar->j = _j;
}

// --------------------------------------------------------------------------------

void setBarProps(Bar* _bar, Node* _n1, Node* _n2, Point _auxvec,
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

void _fillMatrixDefaultValue(double _matrix[SM][SM], double _defaultValue)
{
    for(int i = 0; i < SM; i++) {
        for(int j = 0; j < SM; j++) {
            _matrix[i][j] = _defaultValue;
        }
    }
}

// --------------------------------------------------------------------------------

void _fillMatrixNull(double _matrix[SM][SM])
{
    _fillMatrixDefaultValue(_matrix, 0.0);
}

// --------------------------------------------------------------------------------

void _fillLocalStiffnessMatrix(double _matrix[SM][SM], Bar* _bar)
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

void _fillReducedRotationMatrix(double _matrix[RM][RM], Bar* _bar)
{
    // This is all civil engineering business logic, don't bother at all
    _matrix[0][0] = (_bar->node2->position.x - _bar->node1->position.x) / _bar->l;
    _matrix[0][1] = (_bar->node2->position.y - _bar->node1->position.y) / _bar->l;
    _matrix[0][2] = (_bar->node2->position.z - _bar->node1->position.z) / _bar->l;

    double dx = _bar->auxvec.x - _bar->node1->position.x;
    double dy = _bar->auxvec.y - _bar->node1->position.y;
    double dz = _bar->auxvec.z - _bar->node1->position.z;
    double auxvecLength = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));

    double cosAlfa = (_bar->auxvec.x - _bar->node1->position.x) / auxvecLength;
    double cosBeta = (_bar->auxvec.y - _bar->node1->position.y) / auxvecLength;
    double cosGama = (_bar->auxvec.z - _bar->node1->position.z) / auxvecLength;

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

void _fillRotationMatrix(double _matrix[SM][SM], double _reducedMatrix[RM][RM])
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

void _fillTransposeRotationMatrix(double _matrix[SM][SM], double _other[SM][SM])
{
    for(int i = 0; i < SM; i++) {
        for(int j = 0; j < SM; j++) {
            _matrix[i][j] = _other[j][i];
        }
    }
}

// --------------------------------------------------------------------------------

void _fillGlobalStiffnessMatrix(double _matrix[SM][SM],
                                double _tRotation[SM][SM], double _local[SM][SM])
{
    for(int i = 0; i < SM; i++) {
        for(int j = 0; j < SM; j++) {
            for(int k = 0; k < SM; k++) {
                _matrix[i][j] += _tRotation[i][k] * _local[k][j];
            }
        }
    }
}

// --------------------------------------------------------------------------------

void setStiffnessMatrix(StiffnessMatrix* _sMatrix, Bar* _associatedBar)
{
    _fillMatrixNull(_sMatrix->local);
    _fillMatrixNull(_sMatrix->rotation);
    _fillMatrixNull(_sMatrix->transposeRotation);
    _fillMatrixNull(_sMatrix->global);

    _fillLocalStiffnessMatrix(_sMatrix->local, _associatedBar);
    _fillReducedRotationMatrix(_sMatrix->reducedRotation, _associatedBar);
    _fillRotationMatrix(_sMatrix->rotation, _sMatrix->reducedRotation);
    _fillTransposeRotationMatrix(_sMatrix->transposeRotation, _sMatrix->rotation);
    _fillGlobalStiffnessMatrix(_sMatrix->global,
                               _sMatrix->transposeRotation,
                               _sMatrix->local);
}

// --------------------------------------------------------------------------------

void setFrameBarProps(FrameBar* _frameBar, Node* _n1, Node* _n2, Point _auxvec,
                      Material* _material, Section* _section)
{
    _frameBar->material = *_material;
    _frameBar->section = *_section;

    Bar* pBar = &(_frameBar->bar);
    setBarProps(pBar, _n1, _n2, _auxvec, _material, _section);

    StiffnessMatrix* pStiffnessMatrix = &(_frameBar->stiffnessMatrix);
    setStiffnessMatrix(pStiffnessMatrix, pBar);
}

// --------------------------------------------------------------------------------

FrameBar createFrameBar(Node* _n1, Node* _n2, Point _auxvec,
                        Material* _material, Section* _section)
{
    FrameBar fb;
    setFrameBarProps(&fb, _n1, _n2, _auxvec, _material, _section);

    return fb;
}

// --------------------------------------------------------------------------------

void initNodeArray(NodeArray* _arr, size_t _initSize)
{
    _arr->nodes = malloc(_initSize * sizeof(Node));
    _arr->used = 0;
    _arr->size = _initSize;
}

// --------------------------------------------------------------------------------

void insertNodeArray(NodeArray* _arr, Node _node)
{
    if(_arr->used == _arr->size) {
        _arr->size *= 2;
        _arr->nodes = realloc(_arr->nodes, _arr->size * sizeof(Node));
    }
    _arr->nodes[_arr->used++] = _node;
}

// --------------------------------------------------------------------------------

void freeNodeArray(NodeArray *_arr)
{
    free(_arr->nodes);
    _arr->nodes = NULL;
    _arr->used = _arr->size = 0;
}

// --------------------------------------------------------------------------------

void initFrameBarArray(FrameBarArray* _arr, size_t _initSize)
{
    _arr->framebars = malloc(_initSize * sizeof(FrameBar));
    _arr->used = 0;
    _arr->size = _initSize;
}

// --------------------------------------------------------------------------------

void insertFrameBarArray(FrameBarArray* _arr, FrameBar _bar)
{
    if(_arr->used == _arr->size) {
        _arr->size *= 2;
        _arr->framebars = realloc(_arr->framebars, _arr->size * sizeof(FrameBar));
    }
    _arr->framebars[_arr->used++] = _bar;
}

// --------------------------------------------------------------------------------

void freeFrameBarArray(FrameBarArray *_arr)
{
    free(_arr->framebars);
    _arr->framebars = NULL;
    _arr->used = _arr->size = 0;
}

// --------------------------------------------------------------------------------

void initGlobalSystem(GlobalSystem* _gSys)
{
    NodeArray* pNodeArray = &(_gSys->nodeArray);
    FrameBarArray* pFrameBarArray = &(_gSys->framebarsArray);

    _gSys->numEquations = 0;
    _gSys->numEqFreedoms = 0;
    _gSys->numEqConstraint = 0;

    _gSys->mtxStiffness = NULL;
    _gSys->mtxNodalLoads = NULL;
    _gSys->mtxDisplacements = NULL;
    _gSys->mtxConstraints = NULL;
    _gSys->mtxFreedoms = NULL;
    _gSys->mtxSpreading = NULL;

    initNodeArray(pNodeArray, 1);
    initFrameBarArray(pFrameBarArray, 1);
}

// --------------------------------------------------------------------------------

void insertNodeGlobalSystem(GlobalSystem* _gSys, Node _node)
{
    NodeArray* pNodeArray = &(_gSys->nodeArray);
    insertNodeArray(pNodeArray, _node);

    _gSys->numEqConstraint += _node.constraints;
    _gSys->numEquations = _gSys->nodeArray.used * DOG;
    _gSys->numEqFreedoms = _gSys->numEquations - _gSys->numEqConstraint;
}

// --------------------------------------------------------------------------------

void insertFrameBarGlobalSystem(GlobalSystem* _gSys, FrameBar _bar)
{
    FrameBarArray* pFrameBarArray = &(_gSys->framebarsArray);
    insertFrameBarArray(pFrameBarArray, _bar);
}

// --------------------------------------------------------------------------------

void _fillDynamicMatrix(double** _matrix, int _nR, int _nC, double _initialValue)
{
    for(int i = 0; i < _nR; i++) {
        for(int j = 0; j < _nC; j++) {
            _matrix[i][j] = _initialValue;
        }
    }
}

// --------------------------------------------------------------------------------

void _initStiffnessMatrix(GlobalSystem* _gSys, double initialValue)
{
    const int nEquations = _gSys->numEquations;

    if(_gSys->mtxStiffness == NULL) {
        _gSys->mtxStiffness = malloc(nEquations * sizeof(double*));
        for(int i = 0; i < nEquations; i++) {
            _gSys->mtxStiffness[i] = malloc(nEquations * sizeof(double));
        }
    }
    
    _fillDynamicMatrix(_gSys->mtxStiffness, nEquations, nEquations, initialValue);
}

// --------------------------------------------------------------------------------

void _mountStiffnessMatrix(GlobalSystem* _gSys)
{
    for(int i = 0; i < _gSys->framebarsArray.used; i++) {
        FrameBar* fBars = _gSys->framebarsArray.framebars;
        StiffnessMatrix* fBarStiffMatrix = &(fBars[i].stiffnessMatrix);
        
        for(int j = 0; j < SM; j++) {
            int row = _gSys->mtxSpreading[i][j] - 1;
            for(int k = 0; k < SM; k++) {
                int col = _gSys->mtxSpreading[i][k] - 1;
                _gSys->mtxStiffness[row][col] += fBarStiffMatrix->global[row][col];
            }
        }
    }
}

// --------------------------------------------------------------------------------

void _freeStiffnessMatrix(GlobalSystem* _gSys)
{
    for(int i = 0; i < _gSys->numEquations; i++) {
        double* currentPtr = _gSys->mtxStiffness[i];
        free(currentPtr);
    } free(_gSys->mtxStiffness);

    _gSys->mtxStiffness = NULL;
}

// --------------------------------------------------------------------------------

void _initNodalLoadVector(GlobalSystem* _gSys, double initialValue)
{
    const int nEquations = _gSys->numEquations;

    if(_gSys->mtxNodalLoads == NULL) {
        _gSys->mtxNodalLoads = malloc(nEquations * sizeof(double*));
        for(int i = 0; i < nEquations; i++) {
            _gSys->mtxNodalLoads[i] = malloc(1 * sizeof(double));
        }
    }
    
    _fillDynamicMatrix(_gSys->mtxNodalLoads, nEquations, 1, initialValue);
}

// --------------------------------------------------------------------------------

void _mountNodalLoadVector(GlobalSystem* _gSys)
{
    for(int i = 0; i < _gSys->nodeArray.used; i++) {
        int adjustedIndex = 0;

        double loadVectorCurrentNode[6] = {_gSys->nodeArray.nodes[i].load.fx,
                                           _gSys->nodeArray.nodes[i].load.fy,
                                           _gSys->nodeArray.nodes[i].load.fz,
                                           _gSys->nodeArray.nodes[i].load.mx,
                                           _gSys->nodeArray.nodes[i].load.my,
                                           _gSys->nodeArray.nodes[i].load.mz};

        for(int j = 0; j < DOG; j++) {
            adjustedIndex = _gSys->mtxFreedoms[j][i] - 1;
            _gSys->mtxNodalLoads[adjustedIndex][0] += loadVectorCurrentNode[j];
        }
    }
}

// --------------------------------------------------------------------------------

void _freeNodalLoadVector(GlobalSystem* _gSys)
{
    for(int i = 0; i < _gSys->numEquations; i++) {
        double* currentPtr = _gSys->mtxNodalLoads[i];
        free(currentPtr);
    } free(_gSys->mtxNodalLoads);

    _gSys->mtxNodalLoads = NULL;
}

// --------------------------------------------------------------------------------

void _initDisplacementsMatrix(GlobalSystem* _gSys, double initialValue)
{
    const int nEquations = _gSys->numEquations;

    if(_gSys->mtxDisplacements == NULL) {
        _gSys->mtxDisplacements = malloc(nEquations * sizeof(double*));
        for(int i = 0; i < nEquations; i++) {
            _gSys->mtxDisplacements[i] = malloc(1 * sizeof(double));
        }
    }
    
    _fillDynamicMatrix(_gSys->mtxDisplacements, nEquations, 1, initialValue);
}

// --------------------------------------------------------------------------------

void _freeDisplacementsMatrix(GlobalSystem* _gSys)
{
    for(int i = 0; i < _gSys->numEquations; i++) {
        double* currentPtr = _gSys->mtxDisplacements[i];
        free(currentPtr);
    } free(_gSys->mtxDisplacements);

    _gSys->mtxDisplacements = NULL;
}

// --------------------------------------------------------------------------------

void _initConstraintsMatrix(GlobalSystem* _gSys, double initialValue)
{
    const int nNodes = _gSys->nodeArray.used;

    if(_gSys->mtxConstraints == NULL) {
        _gSys->mtxConstraints = malloc(DOG * sizeof(double*));
        for(int i = 0; i < DOG; i++) {
            _gSys->mtxConstraints[i] = malloc(nNodes * sizeof(double));
        }
    }
    
    _fillDynamicMatrix(_gSys->mtxConstraints, DOG, nNodes, initialValue);
}

// --------------------------------------------------------------------------------

void _mountConstraintsMatrix(GlobalSystem* _gSys)
{
    for(int i = 0; i < _gSys->nodeArray.used; i++) {
        const bool isTranslationFixedX = !_gSys->nodeArray.nodes[i].translation.x;
        const bool isTranslationFixedY = !_gSys->nodeArray.nodes[i].translation.y;
        const bool isTranslationFixedZ = !_gSys->nodeArray.nodes[i].translation.z;
        const bool isRotationFixedX = !_gSys->nodeArray.nodes[i].rotation.x;
        const bool isRotationFixedY = !_gSys->nodeArray.nodes[i].rotation.y;
        const bool isRotationFixedZ = !_gSys->nodeArray.nodes[i].rotation.z;

        if(isTranslationFixedX) _gSys->mtxConstraints[0][i] = 1;
        if(isTranslationFixedY) _gSys->mtxConstraints[1][i] = 1;
        if(isTranslationFixedZ) _gSys->mtxConstraints[2][i] = 1;
        if(isRotationFixedX) _gSys->mtxConstraints[3][i] = 1;
        if(isRotationFixedY) _gSys->mtxConstraints[4][i] = 1;
        if(isRotationFixedZ) _gSys->mtxConstraints[5][i] = 1;
    }
}

// --------------------------------------------------------------------------------

void _freeConstraintsMatrix(GlobalSystem* _gSys)
{
    for(int i = 0; i < DOG; i++) {
        double* currentPtr = _gSys->mtxConstraints[i];
        free(currentPtr);
    } free(_gSys->mtxConstraints);

    _gSys->mtxConstraints = NULL;
}

// --------------------------------------------------------------------------------

void _initFreedomsMatrix(GlobalSystem* _gSys, double initialValue)
{
    const int nNodes = _gSys->nodeArray.used;

    if(_gSys->mtxFreedoms == NULL) {
        _gSys->mtxFreedoms = malloc(DOG * sizeof(double*));
        for(int i = 0; i < DOG; i++) {
            _gSys->mtxFreedoms[i] = malloc(nNodes * sizeof(double));
        }
    }
    
    _fillDynamicMatrix(_gSys->mtxFreedoms, DOG, nNodes, initialValue);
}

// --------------------------------------------------------------------------------

void _mountFreedomsMatrix(GlobalSystem* _gSys)
{
    int countFreedoms = 0;
    int countConstraints = _gSys->numEqFreedoms;

    for(int i = 0; i < _gSys->nodeArray.used; i++) {
        for(int j = 0; j < DOG; j++) {
            if(_gSys->mtxConstraints[j][i] == 0) {
                countFreedoms++;
                _gSys->mtxFreedoms[j][i] = countFreedoms;
            } else {
                countConstraints++;
                _gSys->mtxFreedoms[j][i] = countConstraints;
            }
        }
    }
}

// --------------------------------------------------------------------------------

void _freeFreedomsMatrix(GlobalSystem* _gSys)
{
    for(int i = 0; i < DOG; i++) {
        double* currentPtr = _gSys->mtxFreedoms[i];
        free(currentPtr);
    } free(_gSys->mtxFreedoms);

    _gSys->mtxFreedoms = NULL;
}

// --------------------------------------------------------------------------------

void _initSpreadingMatrix(GlobalSystem* _gSys, double initialValue)
{
    const double nBars = _gSys->framebarsArray.used;

    if(_gSys->mtxSpreading == NULL) {
        _gSys->mtxSpreading = malloc(nBars * sizeof(double*));
        for(int i = 0; i < nBars; i++) {
            _gSys->mtxSpreading[i] = malloc(SM * sizeof(double));
        }
    }
    
    _fillDynamicMatrix(_gSys->mtxSpreading, nBars, SM, initialValue);
}

// --------------------------------------------------------------------------------

void _mountSpreadingMatrix(GlobalSystem* _gSys)
{
    for(int i = 0; i < _gSys->framebarsArray.used; i++) {
        for(int j = 0; j < DOG; j++) {
            Node* pNode = _gSys->framebarsArray.framebars[i].bar.node1;
            int nIndex = pNode - _gSys->nodeArray.nodes;

            _gSys->mtxSpreading[i][j] = _gSys->mtxFreedoms[j][nIndex];
        }

        for(int j = DOG; j < SM; j++) {
            Node* pNode = _gSys->framebarsArray.framebars[i].bar.node2;
            int nIndex = pNode - _gSys->nodeArray.nodes;
            
            int adjustedIndex = j - DOG;
            _gSys->mtxSpreading[i][j] = _gSys->mtxFreedoms[adjustedIndex][nIndex];
        }
    }
}

// --------------------------------------------------------------------------------

void _freeSpreadingMatrix(GlobalSystem* _gSys)
{
    for(int i = 0; i < _gSys->framebarsArray.used; i++) {
        double* currentPtr = _gSys->mtxSpreading[i];
        free(currentPtr);
    } free(_gSys->mtxSpreading);

    _gSys->mtxSpreading = NULL;
}

// --------------------------------------------------------------------------------

void _initAllGlobalMatrix(GlobalSystem* _gSys)
{
    const double defaultValue = 0.0;

    _initStiffnessMatrix(_gSys, defaultValue);
    _initNodalLoadVector(_gSys, defaultValue);
    _initDisplacementsMatrix(_gSys, defaultValue);
    _initConstraintsMatrix(_gSys, defaultValue);
    _initFreedomsMatrix(_gSys, defaultValue);
    _initSpreadingMatrix(_gSys, defaultValue);
}

// --------------------------------------------------------------------------------

void _mountAllGlobalMatrix(GlobalSystem* _gSys)
{
    _mountConstraintsMatrix(_gSys);
    _mountFreedomsMatrix(_gSys);
    _mountSpreadingMatrix(_gSys);
    _mountStiffnessMatrix(_gSys);
    _mountNodalLoadVector(_gSys);
}

// --------------------------------------------------------------------------------

void _freeAllGlobalMatrix(GlobalSystem* _gSys)
{
    _freeStiffnessMatrix(_gSys);
    _freeNodalLoadVector(_gSys);
    _freeDisplacementsMatrix(_gSys);
    _freeConstraintsMatrix(_gSys);
    _freeFreedomsMatrix(_gSys);
    _freeSpreadingMatrix(_gSys);
}

// --------------------------------------------------------------------------------

void mountGlobalSystem(GlobalSystem* _gSys)
{
    _initAllGlobalMatrix(_gSys);
    _mountAllGlobalMatrix(_gSys);
}

// --------------------------------------------------------------------------------

void freeGlobalSystem(GlobalSystem* _gSys)
{
    NodeArray* pNodeArray = &(_gSys->nodeArray);
    FrameBarArray* pFrameBarArray = &(_gSys->framebarsArray);

    _freeAllGlobalMatrix(_gSys);
    freeFrameBarArray(pFrameBarArray);
    freeNodeArray(pNodeArray);

    _gSys->numEqConstraint = 0;
    _gSys->numEqFreedoms = 0;
    _gSys->numEquations = 0;
}

// --------------------------------------------------------------------------------
