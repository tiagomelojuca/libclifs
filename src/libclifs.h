// --------------------------------------------------------------------------------
// LIBCLIFS - TIAGO MELO JUCA (2021)                                              |
// --------------------------------------------------------------------------------
// This library is for studying purposes only - please don't take it seriously :)
// You'll see a lot of comments here. Yep, I know Clean Code, but, as I said, it's
// just for learning; the purpose is a reference for myself, nothing more
// --------------------------------------------------------------------------------

// HEADER GUARD (BEGIN) -----------------------------------------------------------
#ifndef LIBCLIFS_H
#define LIBCLIFS_H

// MUST BE DEFINED AT COMPILE TIME, BUT ONLY WHEN BUILDING THE DLL ITSELF ---------
#ifdef _WIN32
    #ifdef BUILD_DLL
        #define CLIFS_API __declspec(dllexport)
    #else
        #define CLIFS_API __declspec(dllimport)
    #endif
#else
    #define CLIFS_API
#endif

// MAKE THIS HEADER FILE COMPATIBLE WITH C++ CODE TOO (BEGIN) ---------------------
#ifdef __cplusplus
extern "C" {
#endif

// HEADER CONTENT -----------------------------------------------------------------

#include <stddef.h>
#include <stdbool.h>

#include "definitions.h"

// TYPES DECLARATION
// REMINDER TO MYSELF (SOURCE-> www.tutorialspoint.com/cprogramming/c_typedef.htm)
// By convention, uppercase letters are used for these definitions to
// remind the user that the type name is really a symbolic abbreviation
typedef unsigned char BYTE;
typedef enum { ALL_FREE = 0, ALL_FIX } Constraint;

// You can use typedef to give a name to your user defined data types as well
// For example, you can use typedef with structure to define a new data type
// and then use that data type to define structure variables directly
// We could declare our struct first, and later just do as following:
// typedef struct point Point;
// But let's just do the short way...
typedef struct point {
    double x;
    double y;
    double z;
} Point;

typedef struct degreesOfFreedom {
    bool x;
    bool y;
    bool z;
    int constraints;
} DegreesOfFreedom;

typedef struct nodalLoad {
    double fx;
    double fy;
    double fz;
    double mx;
    double my;
    double mz;
} NodalLoad;

typedef struct node {
    Point position;
    DegreesOfFreedom translation;
    DegreesOfFreedom rotation;
    NodalLoad load;
    int constraints;
} Node;

typedef struct bar {
    Node* node1;
    Node* node2;
    Point auxvec;
    double l;
    double e;
    double g;
    double a;
    double iy;
    double iz;
    double j;
} Bar;

// Optionally, we can use high level abstractions
typedef struct material {
    double e;
    double g;
} Material;

typedef struct section {
    double a;
    double iy;
    double iz;
    double j;
} Section;

typedef struct stiffnessMatrix {
    double local[SM][SM];
    double reducedRotation[RM][RM];
    double rotation[SM][SM];
    double transposeRotation[SM][SM];
    double global[SM][SM];
} StiffnessMatrix;

typedef struct frameBar {
    Bar bar;
    Material material;
    Section section;
    StiffnessMatrix stiffnessMatrix;
} FrameBar;

typedef struct nodeArray {
    Node* nodes;
    size_t used;
    size_t size;
} NodeArray;

typedef struct frameBarArray {
    FrameBar* framebars;
    size_t used;
    size_t size;
} FrameBarArray;

typedef struct globalSystem {
    NodeArray nodeArray;
    FrameBarArray framebarsArray;
    int numEquations;
    int numEqFreedoms;
    int numEqConstraints;
    int** mtxConstraints;
    int** mtxFreedoms;
    int** mtxSpreading;
    double** mtxStiffness;
    double** mtxNodalLoads;
    double** mtxDOFFrees;
    double** mtxPartitionTop;
    double** mtxPartitionBot;
    double** mtxDOFConstrained;
    double** vecLoadsDOFFrees;
    double** vecLoadsDOFConstrained;
    double** vecDisplacementsConstrained;
    double** vecDisplacementsFree;
} GlobalSystem;

// --------------------------------------------------------------------------------

void setPointCoords(Point* _point, double _x, double _y, double _z);
Point createPoint(double _x, double _y, double _z);
double _calcDistBetweenPoints(Point _p1, Point _p2);

void setDegreesOfFreedomProps(DegreesOfFreedom* _dof, bool _x, bool _y, bool _z);
DegreesOfFreedom createDegreesOfFreedom(bool _x, bool _y, bool _z);

void setNodalLoadValues(NodalLoad* _load, double _fx, double _fy, double _fz,
                        double _mx, double _my, double _mz);
NodalLoad createNodalLoad(double _fx, double _fy, double _fz,
                          double _mx, double _my, double _mz);

void setNodeProps(Node* _node, Point _position,
                  DegreesOfFreedom _translation, DegreesOfFreedom _rotation,
                  NodalLoad _load);
Node createNode(Point _position,
                DegreesOfFreedom _translation, DegreesOfFreedom _rotation,
                NodalLoad _load);

void _setBarPropsLowLevel(Bar* _bar, Node* _n1, Node* _n2, Point _auxvec,
                          double _e, double _g,
                          double _a, double _iy, double _iz, double _j);
void setBarProps(Bar* _bar, Node* _n1, Node* _n2, Point _auxvec,
                 Material* _material, Section* _section);

void setMaterialProps(Material* _material, double _e, double _g);
Material createMaterial(double _e, double _g);

void setSectionProps(Section* _section, double _a,
                     double _iy, double _iz, double _j);
Section createSection(double _a, double _iy, double _iz, double _j);

void _fillLocalStiffnessMatrix(double _matrix[SM][SM], Bar* _bar);
void _fillReducedRotationMatrix(double _matrix[RM][RM], Bar* _bar);
void _fillRotationMatrix(double _matrix[SM][SM], double _reducedMatrix[RM][RM]);
void _fillTransposeRotationMatrix(double _matrix[SM][SM], double _other[SM][SM]);
void _fillGlobalStiffnessMatrix(double _matrix[SM][SM],
                                double _tRotation[SM][SM], double _local[SM][SM]);
void setStiffnessMatrix(StiffnessMatrix* _sMatrix, Bar* _associatedBar);

void setFrameBarProps(FrameBar* _frameBar, Node* _n1, Node* _n2, Point _auxvec,
                      Material* _material, Section* _section);
FrameBar createFrameBar(Node* _n1, Node* _n2, Point _auxvec,
                        Material* _material, Section* _section);

void initNodeArray(NodeArray* _arr, size_t _initSize);
void insertNodeArray(NodeArray* _arr, Node _node);
void freeNodeArray(NodeArray *_arr);

void initFrameBarArray(FrameBarArray* _arr, size_t _initSize);
void insertFrameBarArray(FrameBarArray* _arr, FrameBar _bar);
void freeFrameBarArray(FrameBarArray *_arr);

void initGlobalSystem(GlobalSystem* _gSys);
void insertNodeGlobalSystem(GlobalSystem* _gSys, Node _node);
void insertFrameBarGlobalSystem(GlobalSystem* _gSys, FrameBar _bar);

void _initConstraintsMatrix(GlobalSystem* _gSys, int _initValue);
void _mountConstraintsMatrix(GlobalSystem* _gSys);
void _freeConstraintsMatrix(GlobalSystem* _gSys);

void _initFreedomsMatrix(GlobalSystem* _gSys, int _initValue);
void _mountFreedomsMatrix(GlobalSystem* _gSys);
void _freeFreedomsMatrix(GlobalSystem* _gSys);

void _initSpreadingMatrix(GlobalSystem* _gSys, int _initValue);
void _mountSpreadingMatrix(GlobalSystem* _gSys);
void _freeSpreadingMatrix(GlobalSystem* _gSys);

void _initStiffnessMatrix(GlobalSystem* _gSys, double _initValue);
void _mountStiffnessMatrix(GlobalSystem* _gSys);
void _freeStiffnessMatrix(GlobalSystem* _gSys);

void _initNodalLoadVector(GlobalSystem* _gSys, double _initValue);
void _mountNodalLoadVector(GlobalSystem* _gSys);
void _freeNodalLoadVector(GlobalSystem* _gSys);

void _initDOFFreesMatrix(GlobalSystem* _gSys, double _initValue);
void _mountDOFFreesMatrix(GlobalSystem* _gSys);
void _freeDOFFreesMatrix(GlobalSystem* _gSys);

void _initPartitionTopMatrix(GlobalSystem* _gSys, double _initValue);
void _mountPartitionTopMatrix(GlobalSystem* _gSys);
void _freePartitionTopMatrix(GlobalSystem* _gSys);

void _initPartitionBotMatrix(GlobalSystem* _gSys, double _initValue);
void _mountPartitionBotMatrix(GlobalSystem* _gSys);
void _freePartitionBotMatrix(GlobalSystem* _gSys);

void _initDOFConstrainedMatrix(GlobalSystem* _gSys, double _initValue);
void _mountDOFConstrainedMatrix(GlobalSystem* _gSys);
void _freeDOFConstrainedMatrix(GlobalSystem* _gSys);

void _initVecLoadsDOFFrees(GlobalSystem* _gSys, double _initValue);
void _mountVecLoadsDOFFrees(GlobalSystem* _gSys);
void _freeVecLoadsDOFFrees(GlobalSystem* _gSys);

void _initVecLoadsDOFConstrained(GlobalSystem* _gSys, double _initValue);
void _mountVecLoadsDOFConstrained(GlobalSystem* _gSys);
void _freeVecLoadsDOFConstrained(GlobalSystem* _gSys);

void _initVecDisplacementsConstrained(GlobalSystem* _gSys, double _initValue);
void _mountVecDisplacementsConstrained(GlobalSystem* _gSys);
void _freeVecDisplacementsConstrained(GlobalSystem* _gSys);

void _initVecDisplacementsFree(GlobalSystem* _gSys, double _initValue);
void _mountVecDisplacementsFree(GlobalSystem* _gSys);
void _freeVecDisplacementsFree(GlobalSystem* _gSys);

void _initAllGlobalMatrix(GlobalSystem* _gSys);
void _mountAllGlobalMatrix(GlobalSystem* _gSys);
void _freeAllGlobalMatrix(GlobalSystem* _gSys);

void mountGlobalSystem(GlobalSystem* _gSys);
void freeGlobalSystem(GlobalSystem* _gSys);

// MAKE THIS HEADER FILE COMPATIBLE WITH C++ CODE TOO (END) -----------------------
#ifdef __cplusplus
}
#endif

// HEADER GUARD (END) -------------------------------------------------------------
#endif // LIBCLIFS_H

// --------------------------------------------------------------------------------
