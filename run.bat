@echo off
REM One-shot setup + launch for Windows (CMD).
setlocal
set ROOT=%~dp0
set ROOT=%ROOT:~0,-1%

if not exist "%ROOT%\.venv" (
    echo Creating virtual environment...
    python -m venv "%ROOT%\.venv"
    if errorlevel 1 goto :error
)
set PY=%ROOT%\.venv\Scripts\python.exe

echo Installing C++ core (omp_core) + backend deps...
"%PY%" -m pip install --quiet --upgrade pip pybind11
if errorlevel 1 goto :error
"%PY%" -m pip install --quiet -e "%ROOT%\core"
if errorlevel 1 goto :error
"%PY%" -m pip install --quiet -r "%ROOT%\backend\requirements.txt"
if errorlevel 1 goto :error

echo.
echo Starting server at http://127.0.0.1:5000 ...
"%PY%" "%ROOT%\backend\app.py"
goto :eof

:error
echo.
echo Setup failed. Ensure Python 3.9+ and a C++17 compiler (MSVC) are installed.
exit /b 1
