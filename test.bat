@echo off
echo.

echo -- Testing example1.c --
cmd /c "gcc src/libclifs.c tests/example1.c -o build/a.exe && .\build\a"
echo.
echo.