CMAKE_MINIMUM_REQUIRED(VERSION 3.12)

PROJECT(libclifs)
SET(SRC_DIR "src")

ADD_DEFINITIONS(-DBUILD_DLL)

ADD_LIBRARY(
    clifs
    SHARED
    "${SRC_DIR}/libclifs.c"
    "${SRC_DIR}/utils.c"
)

set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/exe_obj")
install(TARGETS clifs DESTINATION "lib")
install(FILES "${SRC_DIR}/libclifs.h" DESTINATION "inc")
