@echo off
echo.

echo -- Testing elements.c --
cmd /c "gcc src/libclifs.c tests/elements.c -o build/a.exe && .\build\a"
echo.
echo.

echo -- Testing datastructs.c --
cmd /c "gcc src/libclifs.c tests/datastructs.c -o build/a.exe && .\build\a"
echo.
echo.

echo -- Testing example1.c --
cmd /c "gcc src/libclifs.c tests/example1.c -o build/a.exe && .\build\a"
echo.
echo.