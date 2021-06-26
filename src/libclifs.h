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

// libclifs declarations goes here

// MAKE THIS HEADER FILE COMPATIBLE WITH C++ CODE TOO (END) -----------------------

#ifdef __cplusplus
}
#endif

// HEADER GUARD (END) -------------------------------------------------------------

#endif // LIBCLIFS_H