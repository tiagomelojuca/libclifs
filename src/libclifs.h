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

// REMINDER TO MYSELF (SOURCE-> www.tutorialspoint.com/cprogramming/c_typedef.htm)
// By convention, uppercase letters are used for these definitions to
// remind the user that the type name is really a symbolic abbreviation
typedef unsigned char BYTE;

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

typedef struct bar {
    double L;
    double E;
    double G;
    double A;
    double Iy;
    double Iz;
    double J;
} Bar;

void fillLocalStiffnessMatrix(double _matrix[12][12], Bar* _bar);

// MAKE THIS HEADER FILE COMPATIBLE WITH C++ CODE TOO (END) -----------------------
#ifdef __cplusplus
}
#endif

// HEADER GUARD (END) -------------------------------------------------------------
#endif // LIBCLIFS_H