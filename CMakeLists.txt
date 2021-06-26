CMAKE_MINIMUM_REQUIRED(VERSION 3.12)

PROJECT(libclifs)
SET(SRC_DIR "src")

ADD_LIBRARY(
    clifs
    SHARED
    "${SRC_DIR}/libclifs.c"
)