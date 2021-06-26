@echo off
cmd /c "gcc tests/test.c src/libclifs.c -o build/a.exe && .\build\a"