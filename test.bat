@echo off
echo.

echo -- Testing example1.c --
cmd /c "gcc -g src/libclifs.c src/utils.c tests/example1.c -o build/a.exe && .\build\a"
echo.
echo.
