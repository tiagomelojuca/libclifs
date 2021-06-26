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

// TYPES DECLARATION
// REMINDER TO MYSELF (SOURCE-> www.tutorialspoint.com/cprogramming/c_typedef.htm)
// By convention, uppercase letters are used for these definitions to
// remind the user that the type name is really a symbolic abbreviation
typedef unsigned char BYTE;

#ifndef __cplusplus
    typedef enum { false = 0, true = !false } bool;
#endif

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
} Node;

typedef struct bar {
    Node node1;
    Node node2;
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
    double local[12][12];
    double reducedRotation[3][3];
    double rotation[12][12];
    double transposeRotation[12][12];
    double global[12][12];
} StiffnessMatrix;

// --------------------------------------------------------------------------------

double _calcDistBetweenPoints(Point _p1, Point _p2);

void _setBarPropsLowLevel(Bar* _bar, Node _n1, Node _n2, Point _auxvec,
                          double _e, double _g,
                          double _a, double _iy, double _iz, double _j);

void setBarProps(Bar* _bar, Node _n1, Node _n2, Point _auxvec,
                 Material* _material, Section* _section);

void _fillMatrixDefaultValue(double _matrix[12][12], double _defaultValue);
void _fillLocalStiffnessMatrix(double _matrix[12][12], Bar* _bar);
void _fillReducedRotationMatrix(double _matrix[3][3], Bar* _bar);
void _fillRotationMatrix(double _matrix[12][12], double _reducedMatrix[3][3]);
void _fillTransposeRotationMatrix(double _matrix[12][12], double _other[12][12]);
void _fillGlobalStiffnessMatrix(double _matrix[12][12],
                                double _tRotation[12][12], double _local[12][12]);

void setStiffnessMatrix(StiffnessMatrix* _sMatrix, Bar* _associatedBar);

// MAKE THIS HEADER FILE COMPATIBLE WITH C++ CODE TOO (END) -----------------------
#ifdef __cplusplus
}
#endif

// HEADER GUARD (END) -------------------------------------------------------------
#endif // LIBCLIFS_H

// --------------------------------------------------------------------------------
